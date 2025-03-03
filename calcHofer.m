%% Set storage folders
foldData=['..',filesep,'Data Repository',filesep,'Data',filesep,'Hofer'];
foldFigs='Figures';


%% Calculate means of stationary states
%Get file names
files=dir(foldData);
files=files(endsWith({files.name},'.csv'));
fnames={files.name};


%Initialize table for mean values
tab=readtable([foldData,filesep,fnames{1}]);
ncol=size(tab,2);
tab=table('Size',[length(fnames),ncol],...
    'VariableTypes',[{'datetime'},repmat({'double'},1,ncol-1)],...
    'VariableNames',tab.Properties.VariableNames);


%Read files, remove outliers, calculate means
for i=1:length(fnames)
    tabloc=readtable([foldData,filesep,fnames{i}]);

    outliers=isoutlier(tabloc{:,{'Tbed','eps','w'}});
    tabloc(any(outliers,2),:)=[];

    tab.Time(i)=tabloc.Time(1);
    tab{i,2:end}=mean(tabloc{:,2:end},1);
end


%Correct particle and tube diameter
tab.d_p=round(tab.d_p,9);
tab.d_t=round(tab.d_t,3);


%% Mean bed temperature of all experiments
Tbed=mean(tab.Tbed);


%% Bed porosities at minimum fluidization
%Combinations of particle and tube diameter
d_p=unique(tab.d_p);
combs=combinations(d_p,unique(tab.d_t));
combs.Properties.VariableNames={'d_p','d_t'};


%Set different markers and line styles for each particle diameter
combs.mkr=cell(height(combs),1);
combs.mkr(combs.d_p==d_p(1))=repmat({'o'},nnz(combs.d_p==d_p(1)),1);
combs.mkr(combs.d_p==d_p(2))=repmat({'x'},nnz(combs.d_p==d_p(2)),1);

combs.LineStyle=cell(height(combs),1);
combs.LineStyle(combs.d_p==d_p(1))=repmat({'-'},nnz(combs.d_p==d_p(1)),1);
combs.LineStyle(combs.d_p==d_p(2))=repmat({'--'},nnz(combs.d_p==d_p(2)),1);


%Storage for eps_mf
combs.eps_mf=NaN(height(combs),1);


%Set up global figure
fig=figure(1);
clf(fig);
ax=gca();
colors=ax.ColorOrder;
hold(ax,'on');


%Prepare global figure formatting
legItems=repmat(line(ax,'Visible','off'),height(combs),1);
x=linspace(0,max(tab.w),100);
mkrSize=18;


for i=1:height(combs)
    %Fit linear regression model to individual particle and tube diameter
    %data
    idx=find(tab.d_p==combs.d_p(i) & tab.d_t==combs.d_t(i));
    mdl=fitlm(tab.w(idx),tab.eps(idx));


    %Remove outliers based on Cook's distance
    outliers=mdl.Diagnostics.CooksDistance>3*mean(mdl.Diagnostics.CooksDistance);
    tab(idx(outliers),:)=[];


    %Refit model with clean data
    idx=tab.d_p==combs.d_p(i) & tab.d_t==combs.d_t(i);
    mdl=fitlm(tab.w(idx),tab.eps(idx));


    %Calculate eps_mf
    combs.eps_mf(i)=predict(mdl,0);


    %Plot global results
    scatter(ax,tab.w(idx),tab.eps(idx),mkrSize,colors(i,:),combs.mkr{i});
    plot(ax,x,predict(mdl,x'),...
        'Color',colors(i,:),...
        'LineStyle',combs.LineStyle{i});


    %Set up global legend
    legItems(i)=plot(ax,NaN,NaN,...
        'Color',colors(i,:),...
        'LineStyle',combs.LineStyle{i},...
        'Marker',combs.mkr{i},'MarkerSize',mkrSize);


    %Plot local results
    figLoc=figure(i+1);
    clf(figLoc);
    axLoc=gca();
    hold(axLoc,'on');

    scatter(axLoc,tab.w(idx),tab.eps(idx),mkrSize,colors(i,:),combs.mkr{i});
    plot(axLoc,x,predict(mdl,x'),...
        'Color',colors(i,:),...
        'LineStyle',combs.LineStyle{i});

    hold(axLoc,'off');


    %Format and save local figure     
    xlabel(axLoc,'w (m/s)');
    ylabel(axLoc,'\epsilon (-)');

    title(axLoc,sprintf('d_p=%.0f µm, d_T=%.0f mm',...
        combs.d_p(i).*10^6,combs.d_t(i).*10^3));
    
    figLoc.Units='centimeters';
    figLoc.Position=[10,5,17,8.5];

    exportgraphics(figLoc,[foldFigs,filesep,...
        sprintf('Hofer_dp%.0fum_dT%.0fmm',...
            combs.d_p(i).*10^6,combs.d_t(i).*10^3),...
        '.tiff']);
end


%Format and save global figure
hold(ax,'off');

legend(ax,legItems,compose('d_p=%.0f µm, d_T=%.0f mm',...
        combs.d_p.*10^6,combs.d_t.*10^3),...
    'Location','bestoutside');

xlabel(ax,'w (m/s)');
ylabel(ax,'\epsilon (-)');

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];

exportgraphics(fig,[foldFigs,filesep,'Hofer_epsOverW.tiff']);


%Calculate eps_mf for each particle diameter
eps_mf=arrayfun(@(i) mean(combs.eps_mf(combs.d_p==d_p(i))),1:numel(d_p));