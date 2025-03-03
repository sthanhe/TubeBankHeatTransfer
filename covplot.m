function covplot(X,labels,figidx)
    %Create default labels if they are missing
    if isempty(labels)
        labels=compose('X%d',1:size(X,2));
    end


    %Get z-scores and calculate covariance
    X=normalize(X,1);
    isconstant=all(isnan(X),1);
    X(:,isconstant)=[];
    labels(isconstant)=[];
    
    C=cov(X,'omitrows');


    %Round covariance values close to 0 and remove constants (var=0)
    C(abs(C)<1e-6)=0;
    isconstant=arrayfun(@(i) C(i,i)==0,1:size(C,1));

    C(isconstant,:)=[];
    C(:,isconstant)=[];
    X(:,isconstant)=[];
    labels(isconstant)=[];


    %Final size properties
    sz=size(C);
    n=sz(1);

    
    %Set up figure
    fig=figure(figidx);
    clf(fig);

    t=tiledlayout(fig,n,n);
    t.TileIndexing='columnmajor';
    t.TileSpacing='none';


    %Plot individual tiles
    ax=cell(sz);
    for i=1:numel(ax)
        ax{i}=nexttile(t,i);


        [r,c]=ind2sub(sz,i);
        if r==c
            %Only plot labels in diagonal
            text(ax{i},0,0,labels{r},'HorizontalAlignment','center');
        else
            %Bin x-axis values, get x coordinates as means of bins
            [N,~,bin]=histcounts(X(:,c));
            x=arrayfun(@(i) mean(X(bin==i,c)),1:length(N));
        

            %Mean-y values of each bin + standard deviation
            ymean=arrayfun(@(i) mean(X(bin==i,r)),1:length(N));
            ystd=arrayfun(@(i) std(X(bin==i,r)),1:length(N));


            %Plot means + standard deviations as error bars
            errorbar(ax{i},x,ymean,ystd,'o');


            %Add cov-values as legend (best location)
            legend(ax{i},compose('cov=%.3f',C(i)),'Location','best');
            % leg.IconColumnWidth=0;    %Since 2024b


            %Format axis
            ax{i}.Box='off';
        end
    end


    %Make axes identical and set symmetric limits around 0
    linkaxes([ax{:}],'xy');
    ax{2}.XLim=(max(ax{2}.XLim)+0.2).*[-1,1];
    ax{2}.YLim=(max(ax{2}.YLim)+0.2).*[-1,1];


    %Remove ticks everywhere but the very left and bottom axes
    noX=[ax{1:end-1,:}];
    noY=[ax{:,2:end}];

    xticklabels(noX,{''});
    yticklabels(noY,{''});


    %Add title to figure
    title(t,'Covariance plot, z-scores');

    
end