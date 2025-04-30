%% Analyze Stationary Test Measurements
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
%This script analyzes the data of the stationary tests and creates all
%published figures.
%
%
%Requires the file "stat_SumPrep.csv" that summarizes the particle 
%dispersion measurements, which gets created by the script "prepStatic" 
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
%   - getBIC.m
%   - mdlPostLoadFx.m
%   - loadGeometry.m
%   - getMdotSstatic.m
%   - postStatic.m
%   - dynamicModel.slx
%   - stat_SumPrep.csv


%% Set data directories
dirStationary='DataStationary';     %Path to directory where stationary simulation data should be stored
dirFigures='../Figures/DynamicSims';               %Path to directory where figures should be stored

%Create directory if it does not exist
if ~isfolder(dirFigures)
    mkdir(dirFigures);
end


%% Load
%Dynamic model, activate fast restart
h2FG='dynamicModel';
sys=load_system(h2FG);
mdlPostLoadFx;

set_param(h2FG,"FastRestart","on");
cleanup=onCleanup(@() set_param(h2FG,"FastRestart","off"));


%Data
flow=readtable([dirStationary,filesep,'stat_SumPrep.csv']);
baffleMat=flow{:,compose('baffleCorr%d',1:nABs-1)};  %Baffle correction factor matrix

run=1:height(flow);   %Runs to analyze


%% Initial state for faster simulations
p0=flow.p0(1);              %Ambient pressure
baffleCorr=baffleMat(1,:);  %Baffle correction factors
Phigate=flow.Phigate(1);    %Weir boundary condition

[bc,Phi,mAC,HAC,mAB]=getBIC(flow(1,:),direction);     %Get other boundary and initial conditions


%Simulate
out=sim(sys,'LoadExternalInput','on','ExternalInput','bc',...
            'LoadInitialState','off');
xInit=out.xFinal;   %Initial state for other simulations = end state of this simulation


%% Run simulations
FG=NaN(height(flow),n);
h=NaN(height(flow),n);

for i=run
    p0=flow.p0(i);              %Ambient pressure
    Phigate=flow.Phigate(i);    %Weir boundary condition
    baffleCorr=baffleMat(i,:);  %Baffle correction factors

    bc=getBIC(flow(i,:),direction);   %Get other boundary conditions

    
    %Simulate
    out=sim(sys,'LoadExternalInput','on','ExternalInput','bc',...
                'LoadInitialState','on','InitialState','xInit');


    %Post processing
    postStatic(out,flow(i,:),x,xChambers,i,dirFigures);

    FGsim=timeseries2timetable(out.FG);
    hsim=timeseries2timetable(out.h);

    FG(i,:)=FGsim{end,:};
    h(i,:)=hsim{end,:};
end


%Deactivate fast restart
delete(cleanup);


%%
center=mean(posBedLevel(2:3));

DeltaFG=FG(:,posBedLevel(3))-FG(:,center);
Deltah=h(:,posBedLevel(3))-h(:,posBedLevel(2));

outliers=Deltah>0;
DeltaFG(outliers)=NaN;
Deltah(outliers)=NaN;

h2FG=fitlm(Deltah,DeltaFG);
save(['..',filesep,'h2FG.mat'],'h2FG');


fig=figure(914);
clf(fig);
t=tiledlayout(fig,1,1,'Padding','tight');
ax=nexttile(t);
colors=ax.ColorOrder;
hold(ax,'on');

mm=1000;
scatter(Deltah*mm,DeltaFG);
p=plot(Deltah*mm,predict(h2FG,Deltah));

hold(ax,'off');

xlabel(ax,'h_4 - h_5 (mm)')
ylabel(ax,'FG(h_4) - FG_2 (-)')

legend(ax,p,compose('R^2 = %.3f',h2FG.Rsquared.Ordinary));


%Export figure for repository
t.Units='centimeters';
t.OuterPosition=[0,0,17,8.5];

fig.Units=t.Units;
fig.Position(3:4)=t.OuterPosition(3:4)+0.5;

exportgraphics(fig,[dirFigures,filesep,'FGoverDeltaH.tiff'],...
    'Resolution',600);




