function para=hyptest(X,y,fx,beta,beta0,idx)
    %X: predictor variables (n,m)
    %y: Response variable (n,1)
    %beta: estimated parameters
    %beta0: estimated parameters of the null hypothesis
    %idx: parameter indices to be included


    %Set up table
    names={'Parameter','H0','Estimate','p'};
    para=table('Size',[length(idx)+1,length(names)],...
        'VariableTypes',[{'string'},repmat({'double'},1,length(names)-1)],...
        'VariableNames',names);
    
    
    %Fill basic values
    para.Parameter=[compose('P%d',idx),{'All'}]';   %Parameter names
    para.H0(1:end-1)=beta0(idx);                    %Null hypothesis
    para.Estimate(1:end-1)=beta;                    %Estimated parameters
    
    
    %p-values of individual parameters
    for i=1:height(para)-1
        beta_null=para.Estimate;
        beta_null(i)=para.H0(i);
    
        [~,p]=ttest(fx(beta_null,X),y);
        para.p(i)=p;
    end
    
    
    %p-value of all parameters combined
    beta_null=para.H0;
    [~,para.p(end)]=ttest(fx(beta_null,X),y);
end




