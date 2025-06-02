%% Analyze Stationary Test Measurements
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
%This script analyzes the data of the stationary tests and creates all
%published figures.
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
%   - getBIC.m
%   - mdlPostLoadFx.m
%   - loadGeometry.m
%   - getMdotSstatic.m
%   - postStatic.m
%   - dynamicModel.slx
%   - stat_SumPrep.csv


%% Runs to analyze
run=1:height(flow);


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


%% Dependence of fluidization on bed level gradients
%Get differences in second chamber
center=mean(posBedLevel(2:3));

DeltaFG=FG(:,posBedLevel(3))-FG(:,center);
Deltah=h(:,posBedLevel(3))-h(:,posBedLevel(2));


%Remove outliers: reversed bed level gradients
outliers=Deltah>0;
DeltaFG(outliers)=NaN;
Deltah(outliers)=NaN;


%Linear regression
h2FG=fitlm(Deltah,DeltaFG);
save([dirOrig,filesep,'h2FG.mat'],'h2FG');


%% Plot
%Set up figure
fig=figure(916);
clf(fig);
t=tiledlayout(fig,1,1,'Padding','tight');
ax=nexttile(t);
colors=ax.ColorOrder;
hold(ax,'on');


%Plot data and lines
mm=1000;
scatter(Deltah*mm,DeltaFG);
p=plot(Deltah*mm,predict(h2FG,Deltah));

hold(ax,'off');


%Set axes labels and legend
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




