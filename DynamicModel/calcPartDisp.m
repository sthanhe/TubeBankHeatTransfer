%% Analyze Particle Dispersion Measurements
%GNU General Public License v3.0
%By Stefan Thanheiser: https://orcid.org/0000-0003-2765-1156
%
%Part of the paper:
%
%Thanheiser, S.; Haider, M.
%Dispersion Model for Level Control of Bubbling Fluidized Beds with 
%Particle Cross-Flow
%Chemical Engineering Research and Design 2025
%
%All data, along with methodology reports and supplementary documentation, 
%is published in the data repository:
%https://doi.org/10.5281/zenodo.7924693
%
%All required files for this script can be found in the software
%repository:
%https://doi.org/10.5281/zenodo.7948224
%
%
%
%This script analyzes the particle dispersion calibration data and creates 
%all published figures.
%
%
%Requires all particle dispersion data files ("partDisp_Run...") in the 
%folder configured below and all auxiliary classes and functions on the
%MATLAB path
%
%Required products, version 24.1:
%   - MATLAB
%   - Statistics and Machine Learning Toolbox
%Necessary files, classes, functions, and scripts:
%   - @DryAir
%   - @FluBed
%   - @implExp
%   - @Sinter
%   - @Orifice
%   - getConstants.m
%   - getProp.m
%   - maineffects.m
%   - partDisp_Run...


%% Set data directories
dirData='../Data/Own';             %Path to directory containing the data
dirStationary='DataStationary';     %Path to directory where stationary simulation data should be stored
dirFigures='../Figures/DynamicSims';               %Path to directory where figures should be stored

%Create storage folders if they do not exist
if ~isfolder(dirStationary)
    mkdir(dirStationary);
end

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
    
    
    % %Remove outliers
    % outliers=isoutlier(pDisp.D2,1);
    % pDisp{outliers,2:end}=NaN;


    %Record table for future analysis
    writetable(pDisp,[dirStationary,filesep,'stat_Run',num2str(i),'.csv']);
    
    
    %Get means
    flow{i,2:end-nPi}=mean(pDisp{:,2:end},1,'omitnan');
    flow.Run(i)=i;
end


%Add persistent bed levels
flow{:,hNames}=flow{:,hNames}+c.hBed;


%Calculate dimensionless variables
flow.pi1=flow.D2./sqrt(c.d_p.^3.*c.g);
flow.pi2=flow.w_e2./sqrt(c.d_p.*c.g);
flow.pi3=mean([flow.w_p4,flow.w_p5],2)./sqrt(c.d_p.*c.g);
flow.pi4=flow.rho_g2./(c.rho_p-flow.rho_g2);
flow.pi5=DryAir.eta(flow.Tbed2)./(sqrt(c.d_p.^3.*c.g).*(c.rho_p-flow.rho_g2));

flow.Ar2=flow.pi4./flow.pi5.^2;


%Record table for future analysis
writetable(flow,[dirStationary,filesep,'stat_SumPartDisp.csv']);
% flow{isnan(flow.D2),2:end}=NaN;


% %% Show relation between pi4, pi5, and Ar
% figidx=5;
% 
% %Main effects plot: group variables around mean test conditions
% pi4grp=[3.2,3.6,4.2].*1e-4;
% pi5invGrp=[787,868,966];
% Argrp=[170,275,420];
% 
% fig=maineffects(flow.pi1,[flow.pi4,flow.pi5.^-1,flow.Ar2],...
%             {pi4grp,pi5invGrp,Argrp},[NaN,NaN,NaN,NaN],figidx,...
%             {'\pi_4 (-)','1/\pi_5 (-)','Ar (-)'},'\pi_1 (-)');
% 
% 
% fig.Units='centimeters';
% fig.Position=[10,5,17,8.5];
% 
% name=[dirFigures,filesep,'Figure',num2str(figidx)];
% exportgraphics(fig,[name,'.eps']);
% exportgraphics(fig,[name,'.tiff'],'Resolution',600);
% savefig(fig,name);


% %% Fit model to data
% monofx=@(beta,X) beta(1).*prod(X.^(beta(2:end)),2);     %Monomial function, beta=[c,eps2,...]
% 
% X=[flow.pi2,1+flow.pi3,flow.Ar2];   %Dependent variables
% beta0=[1e4,1,-2,0.1];               %Initial guess for beta
% 
% mdl=fitnlm(X,flow.pi1,monofx,beta0);    %Fitted model


% %% Identify and remove outliers using Cook's distance
% %Plot Cook's distance
% fig=figure(907);
% clf(fig);
% 
% plotDiagnostics(mdl,'cookd');
% 
% fig.Units='centimeters';
% fig.Position=[10,5,17,8.5];
% 
% exportgraphics(fig,[dirFigures,filesep,'cookd.tiff'],'Resolution',600);
% 
% 
% %Identify and remove outliers = >3*mean of Cook's distance
% outliers=(mdl.Diagnostics.CooksDistance)>3*mean(mdl.Diagnostics.CooksDistance,'omitmissing');
% flow{outliers,2:end}=NaN;


% %% Refit model to clean data
% X1=[flow.pi2,1+flow.pi3,flow.Ar2];
% mdl1=fitnlm(X1,flow.pi1,monofx,beta0);
% 
% 
% %Relative standard errors
% SErel=mdl1.Coefficients.SE./mdl1.Coefficients.Estimate;




