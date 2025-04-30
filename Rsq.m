function Rsq=Rsq(y,yHat,M)
    %y: observed data
    %yHat: estimated / predicted / modeled / fitted data
    %M: (optional) number of regressors excluding intercept
    

    %Normalize inputs
    n=numel(y);
    y=reshape(y,n,1);           
    yHat=reshape(yHat,n,1);     


    %Regular coefficient of determination
    epsHat=yHat-y;                  %Residuals
    yBar=mean(y,'omitmissing');     %Mean of observed data

    SSE=sum(epsHat.^2,'omitmissing');       %Sum of squared errors
    SST=sum((y-yBar).^2,'omitmissing');     %Total sum of squares

    Rsq=1-SSE./SST;     %Coefficient of determination


    %Adjusted coefficient of determination
    if nargin>2
        Rsq=1-(n-1)./(n-M).*(1-Rsq);
    end
end




