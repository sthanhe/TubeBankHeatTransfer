%% Conduct regression analysis
% GNU General Public License v3.0
% By Stefan Thanheiser: https://orcid.org/0000-0003-2765-1156
%
% Part of the paper:
%
% Thanheiser, S.
% Molerus and Wirth's Heat Transfer Model for Bubbling Fluidized Beds: 
% Proposal for an Extended Model Including Immersed Tube Banks and Particle 
% Cross-Flow
%
% All data, along with methodology reports and supplementary documentation, 
% is published in the data repository:
% https://doi.org/10.5281/zenodo.15576311
%
% All required files for this function can be found in the software
% repository: 
% https://doi.org/10.5281/zenodo.15576950
%
%
%
% This function analyses the model function of a regression analysis, 
% scales the starting values of the regression coefficients, conducts the
% regression analysis, removes outliers, and re-conducts the regression
% analysis with the outliers removed. This improves stability of the
% regression analysis.
%
%
%Requires all auxiliary classes and functions on the MATLAB path
%
%Required products, version 24.1:
%   - MATLAB
%   - Statistics and Machine Learning Toolbox
%Necessary classes, functions, files, and scripts:
%   - None


function [mdl,beta]=getFit(X,y,fx,beta0)
    % Inputs:
    % X: regressor matrix, double [m,n]
    % y: response variable, double [m,1]
    % fx: model function in the form fx(b,X), where b are the regression
    %       coefficients, function handle
    % beta0: starting values of the regression coefficients, double
    % 
    % 
    % Outputs:
    % mdl: results of the regression with scaled variables, NonLinearModel
    % beta: resulting regression coefficients (scaled back), double


    %% Scale response variable
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


    %% Identify multiplication and power betas
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

    
    %% Split function string at betas
    betaStart=[betaMultStart,betaPowerStart];
    betaEnd=[betaMultEnd,betaPowerEnd];

    [betaStart,idx]=sort(betaStart);
    betaEnd=betaEnd(idx);

    fxsplit=splitBetween(fxstr,betaStart,betaEnd);
    fxsplit(cellfun(@isempty,fxsplit))=[];


    %% Restructure function and scale beta0
    %Check if betas are scaled more than once
    idx=cellfun(@(x) betaIdx(x),fxsplit);
    multiscale=arrayfun(@(i) nnz(idx==i)>1,1:length(beta0));


    %Replace scaled variables with new variables (X2)
    X2idx=size(X,2)+1;
    X=[X,NaN(size(X,1),numel(betaStart))];
    beta0norm=beta0;
    xmean=ones(1,size(X,2));
    scaleStr=cell(size(fxsplit));
    hasMultBeta=false(size(fxsplit));
    for i=1:length(fxsplit)
        %Skip if beta is scaled multiple times
        idx=betaIdx(fxsplit{i});
        if ~isnan(idx) && multiscale(idx)
            continue;
        end


        %Identify type of beta (multiplication or power)
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
    
                infinite=isinf(X(:,X2idx));
                X2mean=mean(X(~infinite,X2idx),1,'omitmissing');
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


    %% Calculate fit
    %Set options
    opts=statset('fitnlm');
    opts.Display='iter';
    opts.TolFun=1e-12;
    opts.TolX=1e-12;
    opts.MaxIter=10000;


    %Initial regression: adapt derivative step if it is not converging
    isconv=false;
    while ~isconv && opts.DerivStep<10
        try
            mdl=fitnlm(X,y,fx,beta0norm,'Options',opts);

            isconv=true;
        catch ME
            if strcmp(ME.identifier,'stats:nlinfit:NonFiniteFunOutput')
                opts.DerivStep=opts.DerivStep.^(1/3);
            else
                rethrow(ME);
            end
        end
    end


    %Reset warning message stack for later analysis
    lastwarn('','');


    %Remove outliers, redo fit
    if isconv
        outliers=(mdl.Diagnostics.CooksDistance)>...
            3*mean(mdl.Diagnostics.CooksDistance,'omitmissing');
        y(outliers)=[];
        X(outliers,:)=[];
    
        mdl=fitnlm(X,y,fx,beta0norm,'Options',opts);
    else
        throw(ME);
    end


    %% Rescale results
    %Normalized results
    betanorm=mdl.Coefficients.Estimate;


    %Rescale for multiplication betas
    beta=betanorm./(beta0norm./beta0);


    %Rescale for power betas: replace beta0 with normalized results
    for i=1:length(scaleStr)
        idx=betaIdx(scaleStr{i});
        idx=idx(1);

        scaleStr{i}=erase(scaleStr{i},regexpPattern('^beta0\(\d+\)[.*/]+'));
        scaleStr{i}=strrep(scaleStr{i},'beta0','betanorm');

        beta(idx)=betanorm(idx)./eval(scaleStr{i});
    end


end


%% Auxiliary functions
function splitStr=splitBetween(str,startPos,endPos)
    %Split a string between the indicated numerical indices
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
    idxStr=regexp(str,'beta0\((\d+)\)','tokens');
    idx=str2double([idxStr{:}]);
end




