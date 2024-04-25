function [fig,tiles,ax]=maineffects(Y,X,grps,secgrps,figidx,xnames,yname)
    n=length(Y);
    grpsRange=grps;
    

    %Center parameters around grouping values, readjust groups
    ngrps=length(grps);
    Xcentered=NaN(size(X));
    for i=1:ngrps
        [Xcentered(:,i),grps{i}]=center(X(:,i),grps{i});
    end
    

    %Get combinations of parameters
    combs=table2array(combinations(grps{:}));
    ncombs=size(combs,1);

    
    %Assign each observation to a set of parameters
    combbox=permute(combs,[3,2,1]);
    combbox=repmat(combbox,n,1,1);
    Xbox=repmat(Xcentered,1,1,ncombs);

    %assig(i,j)==true --> observation i was recorded at parameter set j
    assig=squeeze(all(Xbox==combbox,2));
    clear('combbox','Xbox');


    %Sort observations and calculate mean of each parameter combination
    nobs=arrayfun(@(x) nnz(assig(:,x)),1:ncombs);
    Ysorted=NaN(ncombs,max(nobs));
    for i=1:ncombs
        idx=find(assig(:,i));
        Ysorted(i,1:length(idx))=Y(idx);
    end
    Ysorted=mean(Ysorted,2,'omitnan');


    % %Throw warning when some of the parameter combinations are missing
    % if any(isnan(Ysorted))
    %     getstr=@(x) strjoin(compose('%s=%f',string(xnames)',combs(x,:)'),', ');
    %     str=arrayfun(getstr,find(isnan(Ysorted)),'UniformOutput',false);
    %     warning(['The following parameter sets are missing:\n%s\n',...
    %              'This will likely lead to misleading results'],strjoin(str,'\n'));
    % end


    
    %Set up figure
    fig=figure(figidx);
    clf(fig);
    colors=gca().ColorOrder;
    ax=cell(1,ngrps);

    
    %Create tiles
    tiles=tiledlayout(fig,'horizontal','TileSpacing','compact');
    for i=1:ngrps
        ax{i}=nexttile;
        hold(ax{i},'on');

        if ~isnan(secgrps(i))
            s=secgrps(i);
            n=length(grps{s});
            legItems=cell(n,1);
            for j=1:n
                %Plot mean of secondary group
                idxcell=arrayfun(@(prim) combs(:,i)==prim & combs(:,s)==grps{s}(j),grps{i},'UniformOutput',false);
                y=cellfun(@(x) mean(Ysorted(x),'omitnan'),idxcell);
                plot(ax{i},grps{i},y,'Color',colors(j,:));


                %Overlay actual values
                idx=any([idxcell{:}],2);
                idx=any(assig(:,idx),2);
                scatter(ax{i},X(idx,i),Y(idx),20,colors(j,:));


                %Create legend item (dummy plot)
                legItems{j}=plot(ax{i},NaN,NaN,'Color',colors(j,:),'Marker','o','MarkerSize',sqrt(20));
            end
            legItems=[legItems{:}];

        else
            %Scatter all values when there is no secondary group
            scatter(ax{i},X(:,i),Y,20,'k');
        end


        %Plot total mean
        y=arrayfun(@(x) mean(Ysorted(combs(:,i)==x),'omitnan'),grps{i});
        plot(ax{i},grps{i},y,'-.k');


        %Set legend
        if ~isnan(secgrps(i))
            name=strsplit(xnames{s},' ');
            if length(name)>1
                unit=string(extractBetween(name{2},'(',')'));
                if strcmp(unit,"-")
                    unit=string();
                end
            else
                unit=string();
            end
            legLabels=compose('%s=%.1f %s',repmat(string(name{1}),n,1),grps{s}',repmat(unit,n,1));
            legend(ax{i},legItems,legLabels,'Location','northoutside');
        end

        hold(ax{i},'off');


        %Set x label
        xlabel(ax{i},xnames{i});


        %Set x limit
        xlimit=[ax{i}.XLim;min(grpsRange{i}),max(grpsRange{i})];
        ax{i}.XLim=boundary(xlimit);
    end
    

    %Configure axes
    ax=[ax{:}];
    if length(ax)>1
        %Set common y limits
        ylimits=boundary(cell2mat(get(ax,'YLim')));
        set(ax,'YLim',ylimits);


        %Remove y tick labels except in first plot
        set(ax(2:end),'YTickLabelMode','manual');
        set(ax(2:end),'YTickLabel',[]);
    end


    %Set y axis label
    ylabel(ax(1),yname);
end


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


function x=boundary(x)
    x=[min(x(:,1)),max(x(:,2))];
end




