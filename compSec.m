function R2=compSec(fx,beta,dirFigs,figidx,small)
    
    
    fname='secData.mat';
    
    
    %% Load data
    load(fname,'sec','pis');
    
    
    %Remove training data (Grewal)
    idx=strcmp(sec.Author,'Grewal');
    sec(idx,:)=[];
    pis(idx,:)=[];
    
    authors=unique(sec.Author);
    
    
    %% Response variables
    y=pis.pi1;
    yMol=sec.Nu_mixMW;
    yEst=fx(beta,pis{:,:})+sec.Nu_gcMW;
    
    
    
    %Olsson / Wiman: good fit by extended model
    idx=contains(sec.Author,{'Olsson','Wiman'});
    xFit=y(idx);
    yFit=yEst(idx);
    
    % extFit=fitnlm(xFit,yFit,@(b,x) b+x,0);
    
    
    %Kim / Eder: apparent bias in the estimates
    idx=contains(sec.Author,{'Kim','Eder'});
    xBias=y(idx);
    yBias=yEst(idx);
    
    extBias=fitnlm(xBias,yBias,@(b,x) b+x,0);
    
    
    %Molerus/Wirth
    molBias=fitnlm(y,yMol,@(b,x) b+x,0);
    
    
    %% Coefficients of determination
    names={'extFit','extBias','molFit','molBias'};
    
    R2=table('Size',[length(names),3],...
        'VariableTypes',{'string','double','double'},...
        'VariableNames',{'Model','Rsquared','RMSE'});
    
    R2.Model=names';
    
    R2.Rsquared(strcmp(R2.Model,'extFit'))=Rsq(xFit,yFit,10);
    R2.Rsquared(strcmp(R2.Model,'extBias'))=Rsq(xBias,yBias-extBias.Coefficients.Estimate,10);
    R2.Rsquared(strcmp(R2.Model,'molFit'))=Rsq(y,yMol,7);
    R2.Rsquared(strcmp(R2.Model,'molBias'))=Rsq(y,yMol-molBias.Coefficients.Estimate,7);
    
    R2.RMSE(strcmp(R2.Model,'extFit'))=rmse(yFit,xFit,'omitmissing');
    R2.RMSE(strcmp(R2.Model,'extBias'))=extBias.RMSE;
    R2.RMSE(strcmp(R2.Model,'molFit'))=rmse(yMol,y,'omitmissing');
    R2.RMSE(strcmp(R2.Model,'molBias'))=molBias.RMSE;
    
    
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
    
    
    mkr={'o','+','*','x','square','diamond','^','v','>','<','pentagram','hexagram'};
    legItems=cell(1,numel(authors));
    
    for i=1:length(authors)
        idx=strcmp(sec.Author,authors(i));
        if contains(authors(i),{'Eder','Kim'})
            color=colors(2,:);
        else
            color=colors(1,:);
        end
    
        legItems{i}=scatter(ax{1},y(idx),yEst(idx),18,mkr{i},...
            'MarkerEdgeColor',color);
    end
    
    eq=linspace(0,max([y,yEst,yMol],[],'all'),100);
    
    plot(ax{1},eq,eq,'Color',colors(1,:));
    plot(ax{1},eq,predict(extBias,eq'),'Color',colors(2,:),'LineStyle','--');
    
    x0=-predict(extBias,0);
    quiver(ax{1},x0,0,0,x0,'off','Color','k','MaxHeadSize',0.5);
    quiver(ax{1},x0,x0,0,-x0,'off','Color','k','MaxHeadSize',0.5);
    
    hold(ax{1},'off');
    
    
    text(ax{1},1.5*x0,0.5*x0,compose('bias=%.3f',x0),...
        'BackgroundColor','w',...
        'FontSize',7);
    
    
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
        
        legItems=cell(1,2);
        
        legItems{1}=scatter(ax{2},y,yMol,18,'MarkerEdgeColor','k');
        plot(ax{2},eq,eq,'Color','k');
        plot(ax{2},eq,predict(molBias,eq'),'Color','k','LineStyle','--')
        
        idx=sec.Ar>1e4;
        legItems{2}=scatter(ax{2},y(idx),yMol(idx),18,'x','MarkerEdgeColor',colors(1,:));
        
        x0=-predict(molBias,0);
        quiver(ax{2},x0,0,0,x0,'off','Color','k','MaxHeadSize',1);
        quiver(ax{2},x0,x0,0,-x0,'off','Color','k','MaxHeadSize',1);
        
        hold(ax{2},'off');
        
        
        legItems=[legItems{:}];
        legend(ax{2},legItems,{'All authors','Ar > 1e4'},...
            'Location','northwest',...
            'FontSize',7);
        
        text(ax{2},1.7*x0,0.8*x0,compose('bias=%.3f',x0),...
            'BackgroundColor','w',...
            'FontSize',7);
        
        ax{2}.YTick=[];
        ax{2}.YTickLabel=[];
        
        title(ax{2},'Molerus / Wirth');
    end
    
    
    %% Axes configuration
    ax=[ax{:}];
    
    lim=[min(eq),max(eq)];
    set(ax,'XLim',lim);
    set(ax,'YLim',lim);
    

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
        
        fig.Units=t.Units;
        fig.Position(3:4)=t.OuterPosition(3:4)+0.5;
        
        exportgraphics(fig,[fname,'.tiff'],'Resolution',600);
        
        
        %Export figure for Elsevier
        t.InnerPosition=[1.5,1,14,7];
        fig.Position(3:4)=t.OuterPosition(3:4)+0.5;
        
        exportgraphics(fig,[fname,'.eps']);
    end


end




