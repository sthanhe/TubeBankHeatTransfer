%% Prepare Stationary Test Analysis
% GNU General Public License v3.0
% By Stefan Thanheiser: https://orcid.org/0000-0003-2765-1156
%
% Part of the paper:
%
% Thanheiser, S.; Haider, M.
% Molerus and Wirth's Heat Transfer Model for Bubbling Fluidized Beds: 
% Proposal for an Extended Model Including Immersed Tube Banks and Particle 
% Cross-Flow
% 
% Slight adaptation of the file of the same name in:
%
% S. Thanheiser, Particle Dispersion Model Software. (Feb. 07, 2025). 
% Zenodo. doi: 10.5281/zenodo.14833128.
%
%
%
%This script prepares the data of the stationary tests for further analysis
%by the "calcStatic" script.
% 
%
%Required products, version 24.1:
%   - MATLAB
%   - Simulink
%   - Requirements Toolbox
%   - Simulink Real-Time
%   - Stateflow
%Necessary files, classes, functions, and scripts:
%   - @DryAir
%   - @FluBed
%   - @implExp
%   - @Sinter
%   - baffleCalib.m
%   - getBCF.m
%   - getBIC.m
%   - mdlPostLoadFx.m
%   - loadGeometry.m
%   - getMdotSstatic.m
%   - dynamicModel.slx
%   - stat_SumPartDisp.csv


%% Prepare analysis
%Get constants
c=getConstants();


%Retrieve filenames
dirCont=dir(dirData);   %Content of directory containing the data
files={dirCont(~[dirCont.isdir]).name}';
files=files(startsWith(files,'heatTransfer_'));

%Retain natural order of runs
idx=str2double(extract(files,digitsPattern));
[~,idx]=sort(idx);
files=files(idx);


%Set up table for mean values
chambers=1:6;
nPi=5;
hNames=compose('h%d',chambers);
names=[{'Run','p0','mDotSand','mDotS',...
        'Tbed2','rho_g2','D2','Ar2',...
        'w_e2','wmf2','FG2','FG3',...
        'Tleft','Tcenter','Tright',...
        'epsLeft','epsCenter','epsRight'},...
        compose('AC%d',1:2),compose('AC%dset',1:2),...
        compose('air%d',1:4),...
        hNames,...
        compose('w_p%d',chambers),...
        compose('Phi%d',chambers),...
        compose('pi%d',1:nPi)];
flow=table('Size',[length(files),length(names)],...
            'VariableTypes',repmat({'double'},1,length(names)));
flow.Properties.VariableNames=names;
clear('names');


%% Read individual files and do calculations
for i=1:length(files)
    %Get properties
    tab=readtable([dirData,filesep,files{i}]);
    pDisp=getProp(tab,c,flow.Properties.VariableNames(2:end-nPi),chambers);
    
    
    %Get means
    flow{i,2:end-nPi}=mean(pDisp{:,2:end},1,'omitnan');
    flow.Run(i)=i;
end


%Add persistent bed levels
flow{:,hNames}=flow{:,hNames}+c.hBed;


%Add weir boundary condition and individual baffle correction factors
flow.Phigate=zeros(height(flow),1);

flow.baffleCorr1=ones(height(flow),1);
flow.baffleCorr2=ones(height(flow),1);
flow.baffleCorr3=ones(height(flow),1);


%% Do baffle calibration
baffleCalib;




