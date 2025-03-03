%% Set storage folders
foldData=['..',filesep,'Data Repository',filesep,'Data',filesep,'Eder'];
foldTab=['..',filesep,'Data Repository',filesep,'Data',filesep,'Source',...
    filesep,'Eder',filesep,'temperatures.csv'];


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
    tab{i,2:end}=mean(tabloc{:,2:end},1,'omitmissing');
end


%Correct particle mass flows
tab.mDot=round(tab.mDot.*60^2)./60^2;


%% Mean bed temperatures
mDot=unique(tab.mDot);
TempsAll=table('Size',[numel(mDot),3],...
    'VariableTypes',{'double','double','double'},...
    'VariableNames',{'mDot','n','Tbed'});

TempsAll.mDot=mDot.*60^2;
TempsAll.n=arrayfun(@(i) nnz(tab.mDot==mDot(i)),1:numel(mDot))';
TempsAll.Tbed=arrayfun(@(i) mean(tab.Tbed(tab.mDot==mDot(i))),1:numel(mDot))';


%Approximate mean bed temperatures at reported circulation rates 
mDot=[0;88;206;323];
Temps=table('Size',[numel(mDot),2],...
    'VariableTypes',{'double','double'},...
    'VariableNames',{'mDot','Tbed'});

Temps.mDot=mDot;
Temps.Tbed(mDot==0)=TempsAll.Tbed(TempsAll.mDot==0);

idx=TempsAll.mDot==80 | TempsAll.mDot==100;
Temps.Tbed(mDot==88)=mean(TempsAll.Tbed(idx),"Weights",TempsAll.n(idx));

idx=TempsAll.mDot==175 | TempsAll.mDot==200 | TempsAll.mDot==250;
Temps.Tbed(mDot==206)=mean(TempsAll.Tbed(idx),"Weights",TempsAll.n(idx));

idx=TempsAll.mDot==300 | TempsAll.mDot==310 | TempsAll.mDot==400;
Temps.Tbed(mDot==323)=mean(TempsAll.Tbed(idx),"Weights",TempsAll.n(idx));


%Write table
writetable(Temps,foldTab);