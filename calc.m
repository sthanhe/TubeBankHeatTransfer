%% Set data directories
dirData='Data';             %Path to directory containing the data
dirFigures='Figures';               %Path to directory where figures should be stored

if ~isfolder(dirFigures)
    mkdir(dirFigures);
end


%% Prepare analysis
%Get constants
c=getConstants();


%Retrieve filenames
dirCont=dir(dirData);   %Content of directory containing the data
files={dirCont(~[dirCont.isdir]).name}';
files=files(contains(files,'heatTransfer_'));

%Retain natural order of runs
idx=str2double(extract(files,digitsPattern));
[~,idx]=sort(idx);
files=files(idx);


%Set up table for mean values
names={'Run','P','Tsurf','Tsand','TCO2','hVirt','FG','mDotSand','mode'};
htc=table('Size',[length(files),length(names)],'VariableTypes',[repmat({'double'},1,length(names)-1),'logical']);
htc.Properties.VariableNames=names;
clear('names');


%% Read individual files and do calculations
for i=1:length(files)
    %Get properties
    tab=readtable([dirData,filesep,files{i}]);
    TT=getProp(tab,c,htc.Properties.VariableNames(2:end));
    
    
    %Remove outliers
    outliers=isoutlier(TT.hVirt,1);
    TT{outliers,2:end-1}=NaN;


    %Record table for future analysis
    % writetable(TT,[dirStationary,filesep,'stat_Run',num2str(i),'.csv']);
    
    
    %Get means
    htc{i,2:length(TT.Properties.VariableNames)-1}=mean(TT{:,2:end-1},1,'omitnan');
    htc.Run(i)=i;
    htc.mode(i)=nnz(TT.mode)>height(htc)/2;
end


%Calculate dimensionless variables
% htc.pi1=(htc.wmf2.*c.d_p./htc.Gamma2).^-1;
% htc.pi2=htc.FG2-1;


%Record table for future analysis
htc(htc.hVirt<5,:)=[];
writetable(htc,[dirData,filesep,'htc_Sum.csv']);



%%
param=[htc.TCO2-273.15,htc.FG,htc.mDotSand,htc.mode];
grps={[70,185,300],[3,3.75,4.5],[0.3,0.75,1.5,3],[0,1]};
maineffects(htc.hVirt,param,grps,[3,3,1,3],1,{'TCO2','FG','mDotSand','mode'},'hVirt');











