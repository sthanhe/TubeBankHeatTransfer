function checkFit(X,y,fx,beta,name,figidx,small,dirFigs)
    yPred=fx(beta,X);
    % isNeg=yPred<-1e-3;
    % yPred(isNeg)=[];
    % y(isNeg)=[];


    if small
        mkrsz=18;
        inPos=[1.5,1,4.7,4.7];
        fname=[dirFigs,filesep,'fit_',name];
    else
        mkrsz=36;
        inPos=[2,1.5,12,12];
        fname=[dirFigs,filesep,'Model_',name];


        [~,id]=lastwarn();
        txt=['Model ',name];
        switch id
            case 'stats:nlinfit:ModelConstantWRTParam'
                txt=[txt,', Insensitive Parameters'];
            case 'stats:nlinfit:IllConditionedJacobian'
                txt=[txt,', Ill-Conditioned Jacobian'];
            case ''
        end
    
        isNeg=find(beta<0);
        negTitle=compose(', \\beta_%d=%.3f',isNeg,beta(isNeg));
        txt=[txt,negTitle{:}];
    end


    if figidx==0
        fig=figure(100);
    else
        fig=figure(figidx);
    end
    clf(fig);
    t=tiledlayout(fig,1,1,'Padding','tight');
    ax=nexttile(t);
    colors=ax.ColorOrder;
    hold(ax,'on');
    
    scatter(ax,y,yPred,mkrsz);
    
    lim=max([ax.XLim(2),ax.YLim(2)]);
    eq=linspace(0,lim,100);
    
    plot(ax,eq,eq,'Color',colors(2,:));
    plot(ax,eq,eq.*1.2,'Color',colors(2,:),'LineStyle','--');
    plot(ax,eq,eq./1.2,'Color',colors(2,:),'LineStyle','--');
    
    
    hold(ax,'off');
    
    
    % ax.XLim=[0,lim];
    % ax.YLim=[0,lim];

    
    
    if small
        ax.Visible='off';
    else
        xlabel(ax,'Measured \pi_1 (-)');
        ylabel(ax,'Predicted \pi_1 (-)');

        title(ax,txt);
    end


    %Export figure for repository
    t.Units='centimeters';
    t.InnerPosition=inPos;
    
    fig.Units=t.Units;
    fig.Position(3:4)=t.OuterPosition(3:4)+0.5;
    
    exportgraphics(fig,[fname,'.tiff'],'Resolution',600);
end




% function plotFit(X,y,fx,beta,name,figidx)
%     yPred=fx(beta,X);
%     % isNeg=yPred<-1e-3;
%     % yPred(isNeg)=[];
%     % y(isNeg)=[];
% 
% 
%     fig=figure(figidx);
%     clf(fig);
%     ax=gca;
%     colors=ax.ColorOrder;
%     hold(ax,'on');
% 
%     scatter(ax,y,yPred,18);
% 
%     lim=max([ax.XLim(2),ax.YLim(2)]);
%     eq=linspace(0,lim,100);
% 
%     plot(ax,eq,eq,'Color',colors(2,:));
%     plot(ax,eq,eq.*1.2,'Color',colors(2,:),'LineStyle','--');
%     plot(ax,eq,eq./1.2,'Color',colors(2,:),'LineStyle','--');
% 
% 
%     hold(ax,'off');
% 
% 
%     ax.XLim=[0,lim];
%     ax.YLim=[0,lim];
% 
%     ax.Visible='off';
% 
%     fig.Units='centimeters';
%     fig.Position=[10,5,6.23,5.85];
% 
%     exportgraphics(fig,[name,'.tiff'],'Resolution',600);
% 
% 
% 
%     % xlabel(ax,'Measured \pi_1 (-)');
%     % ylabel(ax,'Predicted \pi_1 (-)');
% 
% 
%     % [~,id]=lastwarn();
%     % t=['Model ',name];
%     % switch id
%     %     case 'stats:nlinfit:ModelConstantWRTParam'
%     %         t=[t,', Insensitive Parameters'];
%     %     case 'stats:nlinfit:IllConditionedJacobian'
%     %         t=[t,', Ill-Conditioned Jacobian'];
%     %     case ''
%     % end
%     % 
%     % isNeg=find(beta<0);
%     % negTitle=compose(', \\beta_%d=%.3f',isNeg,beta(isNeg));
%     % t=[t,negTitle{:}];
%     % 
%     % title(ax,t);
% end




