%% Compare results to other secondary data
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
% This function compares the results from particle-convective regressions 
% in "calcPC" to collected secondary data other than the data from Grewal
% and Saxena.
%
%
%Requires all auxiliary classes and functions on the MATLAB path
%
%Required products, version 24.1:
%   - MATLAB
%   - Statistics and Machine Learning Toolbox
%Necessary classes, functions, files, and scripts:
%   - Rsq.m
%   - secData.mat --> created by the script "prepSec"


function R2=compSec(fx,beta,pis,sec,dirFigs,figidx,small)
    % Inputs:
    % fx: model function in the form fx(b,X), where b are the regression
    %       coefficients, function handle
    % beta: regression coefficients, double
    % dirFigs: path to directore where figures should be stored, char
    % figidx: index of figure window, double
    % small: indicator whether to create the small figure version, logical
    % 
    % Outputs:
    % R2: table of charasteristic values with respect to the goodness of fit
    %       (largely the coefficient of determination R²)
    
    
    %% Remove training data (Grewal)
    idx=strcmp(sec.Author,'Grewal');
    sec(idx,:)=[];
    pis(idx,:)=[];
    
    authors=unique(sec.Author);
    
    
    %% Response variables
    %General
    y=pis.pi1;                              %Measured values
    yMW=sec.Nu_mixMW;                       %Molerus and Wirth (original)
    yExt=fx(beta,pis{:,:})+sec.Nu_gcMW;     %Extended model
    
    
    %Olsson / Wiman: good fit by extended model
    idx=contains(sec.Author,{'Olsson','Wiman'});
    xFit=y(idx);
    yFit=yExt(idx);
    
    
    %Kim / Eder: apparent bias in the estimates
    idx=contains(sec.Author,{'Kim','Eder'});
    xBias=y(idx);
    yBias=yExt(idx);
    

    %Biases
    biasMW=fitnlm(y,yMW,@(b,x) b+x,0);          %Molerus and Wirth
    biasExt=fitnlm(xBias,yBias,@(b,x) b+x,0);   %Extended model    
    
    
    %% Coefficients of determination
    names={'extFit','extBias','MWfit','MWbias'};
    nExt=size(pis,2);   %Number of regressors in the Extended Model
    nMW=7;              %Number of regressors in Moleruns and Wirth's model
    
    R2=table('Size',[length(names),3],...
        'VariableTypes',{'string','double','double'},...
        'VariableNames',{'Model','Rsquared','RMSE'});
    
    R2.Model=names';
    
    R2.Rsquared(strcmp(R2.Model,'extFit'))=Rsq(xFit,yFit,nExt);
    R2.Rsquared(strcmp(R2.Model,'extBias'))=Rsq(xBias,yBias-biasExt.Coefficients.Estimate,nExt);
    R2.Rsquared(strcmp(R2.Model,'MWfit'))=Rsq(y,yMW,nMW);
    R2.Rsquared(strcmp(R2.Model,'MWbias'))=Rsq(y,yMW-biasMW.Coefficients.Estimate,nMW);
    
    R2.RMSE(strcmp(R2.Model,'extFit'))=rmse(yFit,xFit,'omitmissing');
    R2.RMSE(strcmp(R2.Model,'extBias'))=biasExt.RMSE;
    R2.RMSE(strcmp(R2.Model,'MWfit'))=rmse(yMW,y,'omitmissing');
    R2.RMSE(strcmp(R2.Model,'MWbias'))=biasMW.RMSE;
    
    
    %% Set up figure
    if figidx==0
        fig=figure(100);
    else
        fig=figure(figidx);
    end
    clf(fig);


    %Set up axes and tiles
    if small
        ax=cell(1,1);
        t=tiledlayout(fig,1,1,'Padding','tight');
    else
        ax=cell(1,2);
        t=tiledlayout(fig,1,2);
        t.TileIndexing='columnmajor';
        t.TileSpacing='none';
    end
    
    
    %% Figure, left side: extended model
    ax{1}=nexttile(t,1);
    colors=ax{1}.ColorOrder;
    hold(ax{1},'on');
    
    
    %Scatter plot for each author
    mkr={'o','+','*','x','square','diamond','^','v','>','<'};
    legItems=cell(1,numel(authors));
    
    for i=1:length(authors)
        idx=strcmp(sec.Author,authors(i));
        if contains(authors(i),{'Eder','Kim'})
            color=colors(2,:);
        else
            color=colors(1,:);
        end
    
        legItems{i}=scatter(ax{1},y(idx),yExt(idx),18,mkr{i},...
            'MarkerEdgeColor',color);
    end
    

    %Plot equivalence lines
    eq=linspace(0,max([y,yExt,yMW],[],'all'),100);
    
    plot(ax{1},eq,eq,'Color',colors(1,:));
    plot(ax{1},eq,predict(biasExt,eq'),'Color',colors(2,:),'LineStyle','--');
    

    %Add arrow indicating bias
    x0=-predict(biasExt,0);
    quiver(ax{1},x0,0,0,x0,'off','Color','k','MaxHeadSize',0.5);
    quiver(ax{1},x0,x0,0,-x0,'off','Color','k','MaxHeadSize',0.5);
    
    hold(ax{1},'off');
    
    
    %Add arrow annotation
    text(ax{1},1.2*x0,0.5*x0,compose('bias=%.3f',x0),...
        'BackgroundColor','w',...
        'FontSize',7,...
        'Margin',1e-3);
    
    
    %Add axes labels, legend, and title
    if ~small
        legItems=[legItems{:}];
        legend(ax{1},legItems,authors,...
            'Location','northwest',...
            'FontSize',7);
        
        
        ylabel(ax{1},'Estimated Nu (-)');
        
        title(ax{1},'Extended Model');
    end
    
    
        %% Figure, right side: Molerus / Wirth
    if ~small
        ax{2}=nexttile(t,2);
        hold(ax{2},'on');
        

        %Plot data and lines
        legItems=cell(1,2);
        
        legItems{1}=scatter(ax{2},y,yMW,18,'MarkerEdgeColor','k');
        plot(ax{2},eq,eq,'Color','k');
        plot(ax{2},eq,predict(biasMW,eq'),'Color','k','LineStyle','--')
        
        idx=sec.Ar>1e4;
        legItems{2}=scatter(ax{2},y(idx),yMW(idx),18,'x',...
            'MarkerEdgeColor',colors(1,:));
        

        %Add arrow indicating bias
        x0=-predict(biasMW,0);
        quiver(ax{2},x0,0,0,x0,'off','Color','k','MaxHeadSize',1);
        quiver(ax{2},x0,x0,0,-x0,'off','Color','k','MaxHeadSize',1);
        
        hold(ax{2},'off');


        %Add arrow annotation
        text(ax{2},2.3*x0,1*x0,compose('bias=%.3f',x0),...
            'BackgroundColor','w',...
            'FontSize',7);
        
        
        %Add legend
        legItems=[legItems{:}];
        legend(ax{2},legItems,{'All authors','Ar > 1e4'},...
            'Location','northwest',...
            'FontSize',7);
        

        %Remove tick labels
        ax{2}.YTick=[];
        ax{2}.YTickLabel=[];
        

        %Add title
        title(ax{2},'Molerus / Wirth');
    end
    
    
    %% Axes configuration
    ax=[ax{:}];
    

    %Axes limits
    lim=[min(eq),max(eq)];
    set(ax,'XLim',lim);
    set(ax,'YLim',lim);
    

    %Axes appearance
    if small
        %Turn off axes visibility
        ax.Visible='off';
    else
        %Set font size, link axes, set labels
        set(ax,'FontSize',7);
        linkaxes(ax,'y');
        xlabel(ax,'Observed Nu (-)');
    end
    
    
    %% Print figure
    if small
        %Export figure for repository
        fname=[dirFigs,filesep,'sec_s',num2str(figidx)];
        
        t.Units='centimeters';
        t.InnerPosition=[0.5,0.5,4,4];
        
        fig.Units=t.Units;
        fig.Position(3:4)=t.OuterPosition(3:4)+0.5;
        
        exportgraphics(fig,[fname,'.tiff'],'Resolution',600);


    else
        %Export figure for manuscript
        fname=[dirFigs,filesep,'Figure',num2str(figidx)];
        
        t.Units='centimeters';
        t.InnerPosition=[1.5,1,14,7];

        ax(2).XTickLabel{1}=[];
        
        fig.Units=t.Units;
        fig.Position(3:4)=t.OuterPosition(3:4)+0.5;
        
        exportgraphics(fig,[fname,'.tiff'],'Resolution',600);
        
        
        %Export figure for Elsevier
        t.InnerPosition=[1.5,1,14,7];
        fig.Position(3:4)=t.OuterPosition(3:4)+0.5;
        
        exportgraphics(fig,[fname,'.eps']);
        savefig(fig,fname);
    end


end




