%% Prepare Stationary Test Analysis
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
%This script prepares the data of the stationary tests for further analysis
%by the "calcStatic" script.
%
%
%Requires the file "stat_SumPartDisp.csv" that summarizes the particle 
%dispersion measurements, which gets created by the script "calcPartDisp" 
%and stored in the dirStationary folder ("../DataStationary" by default).
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


%% Set data directories
dirStationary='DataStationary';     %Path to directory where stationary simulation data should be stored
dirFigures='../Figures/DynamicSims';               %Path to directory where figures should be stored

%Create directory if it does not exist
if ~isfolder(dirFigures)
    mkdir(dirFigures);
end


%% Load dynamic model and activate fast restart
mdl='dynamicModel';
sys=load_system(mdl);
mdlPostLoadFx;

set_param(mdl,"FastRestart","on");
cleanup=onCleanup(@() set_param(mdl,"FastRestart","off"));


%% Load data
flow=readtable([dirStationary,filesep,'stat_SumPartDisp.csv']);

%Add variables: weir boundary condition and individual baffle correction factors
flow.Phigate=zeros(height(flow),1);

flow.baffleCorr1=ones(height(flow),1);
flow.baffleCorr2=ones(height(flow),1);
flow.baffleCorr3=ones(height(flow),1);


%% Do baffle calibration
baffleCalib;


%Record table for future analysis
writetable(flow,[dirStationary,filesep,'stat_SumPrep.csv']);




