%% Covariance plot
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
% This function creates a plot to illustrate the covariance between
% variables.
%
%
%Requires all auxiliary classes and functions on the MATLAB path
%
%Required products, version 24.1:
%   - MATLAB
%Necessary classes, functions, files, and scripts:
%   - None


function [fig,t]=covplot(X,labels,figidx)
    % Inputs:
    % X: variable matrix, each column is a variable, double
    % labels: variable names, use empty array to use default values, cell
    %           array of character vectors
    % figidx: index of figure window, double
    % 
    % 
    % Outputs:
    % fig: figure handle
    % t: tile handle


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
        

            %Mean y-values of each bin + standard deviation
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




