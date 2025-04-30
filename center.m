function [x,meancenters]=center(x,approx)
    cats=interp1(approx,approx,x,'nearest','extrap');

    meancenters=arrayfun(@(y) mean(x(cats==y)),approx);
    meancenters(isnan(meancenters))=[];

    if length(meancenters)>1
        x=interp1(meancenters,meancenters,x,'nearest','extrap');
    else
        x=repmat(meancenters,size(x));
    end
end