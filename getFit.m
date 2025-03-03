function [mdl,beta]=getFit(X,y,fx,beta0)
    %% Scale response and predictor variables
    % xmean=mean(X,1,'omitmissing');

    ymean=mean(y,'omitmissing');
    y=y./ymean;


    %% Normalize function
    %Turn function into string
    fxstr=func2str(fx);


    %Get function variable names, erase function handle definition
    vars=regexp(fxstr,'^@\((\w),(\w)\)','tokens');
    vars=vars{1};

    fxstr=erase(fxstr,"@("+wildcardPattern+")");


    %Replace instances of first variable with beta0 and second one with xmean
    fxstr=regexprep(fxstr,[vars{1},'\((\d+)\)'],'beta0($1)');
    fxstr=regexprep(fxstr,[vars{2},'\(:,(\d+)\)'],'xmean(:,$1)');


    % %Replace additive factors to variables with changed variables
    % sums=regexp(fxstr,'\((\d+)([-+])(xmean\(:,\d+\))\)','tokens');
    % if ~isempty(sums)
    %     sums=vertcat(sums{:});
    %     sums(:,1)=cellfun(@str2double,sums(:,1),'UniformOutput',false);
    % 
    %     idx=cellfun(@(x) textscan(x,'xmean(:,%f)'),sums(:,3));
    %     idx=vertcat(idx{:});
    %     for i=1:length(idx)
    %         switch sums{i,2}
    %             case '+'
    %                 X(:,idx(i))=sums{i,1}+X(:,idx(i));
    %             case '-'
    %                 X(:,idx(i))=sums{i,1}-X(:,idx(i));
    %         end
    %     end
    %     fxstr=regexprep(fxstr,'\(\d+[-+](xmean\(:,\d+\))\)','$1');
    % end


    %% Scale beta0
    %Basic function string properties
    l=length(fxstr);    %Length
    x=1:l;              %Position vector


    %Detect signs, operators, and levels delineated by braces
    isSign=getSigns(fxstr);
    isOp=getOps(fxstr);
    levels=getLevels(fxstr);

    
    %Detect starts of multiplication and power betas
    betaMultStart=regexp(fxstr,'(?<!(\^|\^-))beta0');
    betaPowerStart=regexp(fxstr,'(?<=(\^|\^-))beta0');

    betaMultEnd=betaMultStart;
    betaPowerEnd=betaPowerStart;

    betaMultLevel=levels(betaMultStart);
    betaPowerLevel=levels(betaPowerStart);


    %Shift start and end indices for power betas
    for i=1:length(betaPowerStart)
        betaPowerStart(i)=find(...
            (levels==betaPowerLevel(i) & ...
            x<betaPowerStart(i) & ...
            isOp) | ...
            (levels<betaPowerLevel(i) & ...
            x<betaPowerStart(i)),...
            1,'last')+1;
        % betaPowerStart(i)=find(...
        %     (levels==betaPowerLevel(i) & ...
        %     x<betaPowerStart(i) & ...
        %     isOp) | ...
        %     (levels<betaPowerLevel(i) & ...
        %     x<betaPowerStart(i)),...
        %     1,'last')+1;

        betaPowerLength=regexp(fxstr(x>=betaPowerEnd(i)),...
            '^(beta0\(\d+\))','tokenExtents');
        
        betaPowerEnd(i)=betaPowerEnd(i)+betaPowerLength{1}(2)-1;
    end


    %Get start and end indices for multiplication betas
    betaMultStart=[betaMultStart,l];
    
    for i=1:length(betaMultStart)-1
        %Level ends when brace level becomes smaller than beta's level
        levelEnd=find(levels<betaMultLevel(i) & ...
            x>betaMultStart(i) & ...
            x<betaMultStart(i+1),...
            1);


        %Level must end before next multiplication beta
        if isempty(levelEnd)
            levelEnd=betaMultStart(i+1);
        end

        
        %End level sooner if there is a sign on the same level
        sign=find(levels==betaMultLevel(i) & ...
            isSign & ...
            x>betaMultStart(i) & ...
            x<levelEnd,...
            1);

        if isempty(sign)
            betaMultEnd(i)=find(levels==betaMultLevel(i) & ...
                x<levelEnd,...
                1,'last');
        else
            betaMultEnd(i)=sign-1;
        end


        %Check if multiplication beta includes power beta on the same level
        hasPowerBeta=betaMultStart(i)<betaPowerStart & ...
            betaPowerEnd<=betaMultEnd(i) & ...
            betaMultLevel(i)==betaPowerLevel;

        if any(hasPowerBeta)
            %Split substring at power betas and include operators (-2) 
            bias=betaMultStart(i)-1;
            startIdx=betaPowerStart(hasPowerBeta)-2-bias;
            endIdx=betaPowerEnd(hasPowerBeta)-bias;

            newStr=splitBetween(fxstr(betaMultStart(i):betaMultEnd(i)),...
                startIdx,endIdx);


            %Reorder power betas to the end, except those ending with ^
            newStr=newStr([1:2:end,2:2:end]);

            powerEnd=endsWith(newStr,'^');
            newStr=[newStr(~powerEnd);newStr(powerEnd)];


            %Update function string
            fxstr(betaMultStart(i):betaMultEnd(i))=[newStr{:}];


            %Update indices, exclude operators for power betas (+2)
            newIdx=betaMultStart(i)+cumsum(cellfun(@length,newStr));
            n=nnz(hasPowerBeta);
            nP=nnz(powerEnd);

            betaPowerStart(hasPowerBeta)=newIdx(end-n-nP:end-nP-1)+2;
            betaPowerEnd(hasPowerBeta)=newIdx(end-n-nP+1:end-nP)-1;

            betaMultEnd(i)=newIdx(end-n-nP)-1+2;


            %Redo sign, operator and level detection
            isSign=getSigns(fxstr);
            isOp=getOps(fxstr);
            levels=getLevels(fxstr);
        end


        %Power beta on a different level needs to be separated
        hasPowerBeta=betaMultStart(i)<betaPowerStart & ...
            betaPowerEnd<=betaMultEnd(i) & ...
            betaMultLevel(i)~=betaPowerLevel;

        if any(hasPowerBeta)
            betaMultEnd(i)=find(levels==betaMultLevel(i) & ...
                x>betaMultStart(i) & ...
                x<min(betaPowerStart(hasPowerBeta)),...
                1,'last');
        end




        %Ensure that beta is not followed by a division operator
        op=find(isOp & x>betaMultStart(i),1);
        if matches(fxstr(op),'/')
            %Find end of division
            divEnd=find(levels==betaMultLevel(i) & ...
                x>op & ...
                isOp,1)-2;


            %Turn division into multiplication, update function string
            newStr=['(',fxstr(op+1:divEnd),').^-1'];
            fxstr=[fxstr(1:op-1),'*',newStr,fxstr(divEnd+1:end)];


            %Update indices
            betaMultStart(i+1:end)=betaMultStart(i+1:end)+6;
            betaMultEnd(i:end)=betaMultEnd(i:end)+6;

            isLater=betaPowerStart>op;
            betaPowerStart(isLater)=betaPowerStart(isLater)+6;
            betaPowerEnd(isLater)=betaPowerEnd(isLater)+6;


            %Redo sign, operator and level detection
            x=1:length(fxstr);
            isSign=getSigns(fxstr);
            isOp=getOps(fxstr);
            levels=getLevels(fxstr);
        end
    end

    betaMultStart=betaMultStart(1:end-1);

    
    %Split function string at betas
    betaStart=[betaMultStart,betaPowerStart];
    betaEnd=[betaMultEnd,betaPowerEnd];

    [betaStart,idx]=sort(betaStart);
    betaEnd=betaEnd(idx);

    fxsplit=splitBetween(fxstr,betaStart,betaEnd);
    fxsplit(cellfun(@isempty,fxsplit))=[];


    %Replace scaled variables with new variables (X2) and scale beta0
    X2idx=size(X,2)+1;
    X=[X,NaN(size(X,1),numel(betaStart))];
    beta0norm=beta0;
    xmean=ones(1,size(X,2));
    scaleStr=cell(size(fxsplit));
    hasMultBeta=false(size(fxsplit));
    for i=1:length(fxsplit)
        if startsWith(fxsplit{i},'beta0')
            %Is multiplication beta

            %Split expression
            betaSplit=regexp(fxsplit{i},...
                '^(beta0\(\d+\))([.*/\^]+)(.*?)([.*/]*\w*)$','tokens');

            
            if ~isempty(betaSplit) && ...
                    ~isempty(betaSplit{1}{3}) && ...
                    ~contains(betaSplit{1}{2},'^')

                betaSplit=betaSplit{1};


                %Add constant powers
                if ~contains(betaSplit{4},lettersPattern) && ...
                        endsWith(betaSplit{3},'^')
                    betaSplit{3}=[betaSplit{3:4}];
                    betaSplit{4}='';
                end


                %Calculate new variable and its mean
                X2str=strrep(betaSplit{3},'xmean','X');
                X(:,X2idx)=eval(X2str);
    
                X2mean=mean(X(:,X2idx),1,'omitmissing');
                X(:,X2idx)=X(:,X2idx)./X2mean;
                xmean(X2idx)=X2mean;
    
    
                %Scale beta0. It was ensured before that beta0 is followed
                %by a multiplication operator
                idx=betaIdx(betaSplit{1});
                beta0norm(idx)=beta0(idx).*X2mean;
    
    
                %Update function string
                betaSplit{3}=['xmean(:,',num2str(X2idx),')'];
                fxsplit{i}=[betaSplit{:}];
                X2idx=X2idx+1;
            end
    

        elseif endsWith(fxsplit{i},"beta0("+digitsPattern+")") && ... 
            ((i>1 && startsWith(fxsplit{i-1},'beta0')) || ...
            (i>2 && endsWith(fxsplit{i-2},"beta0("+digitsPattern+")") && ...
            matches(fxsplit{i-1},{'.*','./'}) && ...
            hasMultBeta(i-2)))
            %Is power beta and either:
            % fxsplit before is multiplication beta, or 
            % fxsplit two positions before is power beta
            %and those splits are on the same level

            %Split expression
            betaSplit=regexp(fxsplit{i},...
                '(.*?)\.\^(beta0\(\d+\))$','tokens');
            betaSplit=betaSplit{1};
            

            %Calculate new variable and its mean
            X2str=strrep(betaSplit{1},'xmean','X');
            X(:,X2idx)=eval(X2str);

            X2mean=mean(X(:,X2idx),1);
            X(:,X2idx)=X(:,X2idx)./X2mean;
            xmean(X2idx)=X2mean;


            %Find and scale last multiplication beta0
            multBeta=find(startsWith(fxsplit(1:i),'beta0'),1,'last');
            idx=betaIdx(fxsplit{multBeta});
            if endsWith(fxsplit{i-1},'*')
                beta0norm(idx)=beta0norm(idx).*X2mean.^eval(betaSplit{2});
            elseif endsWith(fxsplit{i-1},'/')
                beta0norm(idx)=beta0norm(idx)./X2mean.^eval(betaSplit{2});
            end


            %Update function string
            betaSplit{1}=['xmean(:,',num2str(X2idx),')'];
            betaSplit{2}=['.^',betaSplit{2}];
            fxsplit{i}=[betaSplit{:}];
            X2idx=X2idx+1;


            %Save scale formula for later rescaling
            if isempty(scaleStr{multBeta})
                scaleStr{multBeta}=fxsplit{multBeta};
            end

            if matches(fxsplit{i-1},{'.*','./'})
                scaleStr{multBeta}=[scaleStr{multBeta},fxsplit{i-1}];
            end

            scaleStr{multBeta}=[scaleStr{multBeta},fxsplit{i}];


            %Indicate that this power beta has a multiplication beta
            hasMultBeta(i)=true;
        end
    end


    %Remove unused scaled variables
    X(:,X2idx:end)=[];
    scaleStr(cellfun(@isempty,scaleStr))=[];


    %Add ymean scale factor to function
    fx=eval(['@(beta0,xmean) 1/ymean.*(',[fxsplit{:}],');']);






    % fig=figure(1);
    % clf(fig);
    % ax=gca();
    % colors=ax.ColorOrder;
    % hold(ax,'on');
    % 
    % plot(ax,levels,'Color',colors(1,:));
    % xline(ax,betaStart,'Color',colors(2,:));
    % xline(ax,betaEnd,'Color',colors(2,:),'LineStyle','--');










    % %Split function at "+" or "-" (but not "^-" or "^(-")
    % fxsplit=split(fxstr,regexpPattern('(?<!(\^|\^\())[-+]'));
    % fxsplit=cleanfx(fxsplit);
    % 
    % 
    % %Ensure that regressors without a beta remain unscaled
    % xmean=mean(X,1,'omitmissing');
    % 
    % hasBeta0=contains(fxsplit,'beta0');
    % hasXmean=contains(fxsplit,'xmean');
    % 
    % idxStr=arrayfun(@(i) regexp(fxsplit{i},'xmean\(:,(\d+)\)','tokens'),...
    %     find(hasXmean & ~hasBeta0)','UniformOutput',false);
    % if ~isempty(idxStr)
    %     idxStr=[idxStr{:}];
    %     xmean(str2double([idxStr{:}]))=1;
    % end
    % 
    % 
    % %Remove splits without a beta0
    % fxsplit=fxsplit(hasBeta0);
    % 
    % 
    % %Split again at "^(beta" or "^(-beta"
    % c=cellfun(@(x) count(x,'beta0'),fxsplit);
    % fxsplit2=cell(sum(c),1);
    % counter=1;
    % for i=1:length(fxsplit)
    %     idx=regexp(fxsplit{i},('(\^\(|\^\(-)beta0'));
    %     idx=[1,idx,length(fxsplit{i})]; %#ok<AGROW>
    %     betasplit=arrayfun(@(j) extractBetween(fxsplit{i},idx(j),idx(j+1)),...
    %         1:length(idx)-1)';
    % 
    %     fxsplit2(counter:counter+length(betasplit)-1)=betasplit;
    % 
    %     counter=counter+length(betasplit);
    % end
    % fxsplit2=fxsplit2(~cellfun(@isempty,fxsplit2));
    % fxsplit2=cleanfx(fxsplit2);
    % 
    % 
    % %Calculate conversion factors
    % beta0norm=NaN(size(beta0));
    % for i=1:length(fxsplit2)
    %     %Find beta0 index and calculate beta0norm
    %     idxStr=regexp(fxsplit2{i},'^beta0\((\d+)\)','tokens');
    %     beta0norm(str2double(idxStr{1}))=eval(fxsplit2{i});
    % end
    % 
    % 
    % %beta0norm=beta0 where there is no conversion
    % idx=isnan(beta0norm);
    % beta0norm(idx)=beta0(idx);
    % 
    % 
    % %Add ymean scale factor to function
    % % fx=@(b,x) 1/ymean.*fx(b,x);
    % % beta0norm(1)=beta0norm(1)/ymean;
    % fx=eval(['@(beta0,xmean) 1/ymean.*(',fxstr,');']);


    %% Calculate fit
    %Set options
    opts=statset('fitnlm');
    opts.Display='iter';
    opts.TolFun=1e-12;
    opts.TolX=1e-12;
    % opts.Robust='on';
    opts.MaxIter=10000;


    %Initial model
    % X=X./xmean;
    mdl=fitnlm(X,y,fx,beta0norm,'Options',opts);
    % [betanorm,R,J,CovB,MSE,ErrorModelInfo]=nlinfit(X,y,fx,beta0norm,opts);


    %Remove outliers, redo fit
    outliers=(mdl.Diagnostics.CooksDistance)>...
        3*mean(mdl.Diagnostics.CooksDistance,'omitmissing');
    y(outliers)=[];
    X(outliers,:)=[];

    lastwarn('','');    %Reset warning message stack for later analysis

    mdl=fitnlm(X,y,fx,beta0norm,'Options',opts);


    %% Rescale results
    %Normalized results
    betanorm=mdl.Coefficients.Estimate;


    %Rescale for multiplication betas
    beta=betanorm./(beta0norm./beta0);


    %Rescale for power betas: replace beta0 with normalized results
    for i=1:length(scaleStr)
        idx=betaIdx(scaleStr{i});

        scaleStr{i}=erase(scaleStr{i},regexpPattern('^beta0\(\d+\)[.*/]+'));
        scaleStr{i}=strrep(scaleStr{i},'beta0','betanorm');

        beta(idx)=betanorm(idx)./eval(scaleStr{i});
    end


end


%% Auxiliary functions
function splitStr=splitBetween(str,startPos,endPos)
    splitStr=cell(2*length(startPos)+1,1);

    splitStr{1}=str(1:startPos(1)-1);
    splitStr{2}=str(startPos(1):endPos(1));

    for i=2:length(startPos)
        splitStr{2*i-1}=str(endPos(i-1)+1:startPos(i)-1);
        splitStr{2*i}=str(startPos(i):endPos(i));
    end

    splitStr{end}=str(endPos(end)+1:end);
end


function isSign=getSigns(fxstr)
    %Detect signs that are not part of a power operator
    isSign=false(1,length(fxstr));
    isSign(regexp(fxstr,'(?<!(\^|\^\(|e))[-+]'))=true;
end


function isOp=getOps(fxstr)
    %Detect operators
    isOp=false(1,length(fxstr));
    isOp(regexp(fxstr,'[-+*/]'))=true;
end


function levels=getLevels(fxstr)
    %Get levels delineated by braces 
    levels=zeros(1,length(fxstr));
    levels(strfind(fxstr,'('))=1;
    levels(strfind(fxstr,')'))=-1;
    levels=cumsum(levels);
end


function idx=betaIdx(str)
    %Get beta index
    idxStr=regexp(str,'^beta0\((\d+)\)','tokens');
    idx=str2double(idxStr{1});
end


% function fxsplit=cleanfx(fxsplit)
%     for i=1:length(fxsplit)
%         %Remove leading operators
%         % idx=regexp(fxsplit{i},'(^[.*/^+-\(]+)','tokenExtents');
%         idx=regexp(fxsplit{i},'(^\W+)','tokenExtents');
%         if ~isempty(idx)
%             idx=idx{1};
%             fxsplit{i}(idx(1):idx(2))=[];
%         end
% 
% 
%         %Find braces
%         brace=zeros(size(fxsplit{i}));
%         brace(strfind(fxsplit{i},'('))=1;
%         brace(strfind(fxsplit{i},')'))=-1;
% 
% 
%         %Last brace cannot be an open one
%         lastOpen=find(brace==1,1,'last');
%         lastClose=find(brace==-1,1,'last');
%         if lastClose<lastOpen | (isempty(lastClose) && ~isempty(lastOpen))
%             fxsplit{i}(lastOpen:end)=[];
%             brace(lastOpen:end)=[];
%         end
% 
% 
%         %Ensure consistent bracing
%         brace=cumsum(brace);
%         % lastEqual=find(brace==0,1,'last');
%         % fxsplit{i}=fxsplit{i}(1:lastEqual);
%         firstNeg=find(brace<0,1);
%         if ~isempty(firstNeg)
%             fxsplit{i}=fxsplit{i}(1:firstNeg-1);
%         end
% 
% 
%         %Remove trailing operators
%         idx=regexp(fxsplit{i},'[.*/^+-]+$');
%         fxsplit{i}(idx:end)=[];
%     end
% end


% function s=splitDelim(str,delimiter)
%     %Split, but keep delimiters
%     s=split(str,delimiter);
%     s(cellfun(@isempty,s))=[];
%     s=cellfun(@(x) [delimiter,x],s,'UniformOutput',false);
% end