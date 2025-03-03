function checkFit(X,y,fx,beta,name,figidx)
    yPred=fx(beta,X);
    isNeg=yPred<-1e-3;
    yPred(isNeg)=[];
    y(isNeg)=[];


    fig=figure(figidx);
    clf(fig);
    ax=gca;
    colors=ax.ColorOrder;
    hold(ax,'on');
    
    scatter(ax,y,yPred);
    
    lim=max([ax.XLim(2),ax.YLim(2)]);
    eq=linspace(0,lim,100);
    
    plot(ax,eq,eq,'Color',colors(2,:));
    plot(ax,eq,eq.*1.2,'Color',colors(2,:),'LineStyle','--');
    plot(ax,eq,eq./1.2,'Color',colors(2,:),'LineStyle','--');
    
    
    hold(ax,'off');
    
    
    ax.XLim=[0,lim];
    ax.YLim=[0,lim];

    xlabel(ax,'Measured \pi_1 (-)');
    ylabel(ax,'Predicted \pi_1 (-)');


    [~,id]=lastwarn();
    t=['Model ',name];
    switch id
        case 'stats:nlinfit:ModelConstantWRTParam'
            t=[t,', Insensitive Parameters'];
        case 'stats:nlinfit:IllConditionedJacobian'
            t=[t,', Ill-Conditioned Jacobian'];
        case ''
    end

    isNeg=find(beta<0);
    negTitle=compose(', \\beta_%d=%.3f',isNeg,beta(isNeg));
    t=[t,negTitle{:}];

    title(ax,t);
end




