function Rsq=Rsq(y,yest,M)
    n=numel(y);
    y=reshape(y,n,1);           %Observed data
    yest=reshape(yest,n,1);     %Estimated / predicted / modeled / fitted data

    res=yest-y;                     %Residuals
    ymean=mean(y,'omitmissing');    %Mean of observed data

    SSE=sum(res.^2,'omitmissing');          %Sum of squared errors
    SST=sum((y-ymean).^2,'omitmissing');    %Total sum of squares

    Rsq=1-SSE./SST;     %Coefficient of determination


    if nargin>2
        %M=number of regressors excluding intercept
        n=nnz(y);
        Rsq=1-(n-1)./(n-M).*(1-Rsq);
    end
end