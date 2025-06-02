%% Bed temperature sensor analysis
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
% All data, along with methodology reports and supplementary documentation, 
% is published in the data repository:
% https://doi.org/10.5281/zenodo.15576311
%
% All required files for this script can be found in the software
% repository: see the link to the supplemental release in the data 
% repository
%
%
%
% This script analyzes the measurements from the bed temperature sensors to
% decide which one accurately measured the bed temperature. It creates the
% corresponding figures üublished in the Methodology Report in the data
% repository.
%
%
%Requires all auxiliary classes and functions on the MATLAB path
%
%Required products, version 24.1:
%   - MATLAB
%Necessary classes, functions, files, and scripts:
%   - @figaux
%   - primData.mat --> created by the script "prepPrim"


%% Set data locations
dirFigs='Figures';      %Figure storage folder
fname='primData.mat';   %Primary data file


%% Make folders if they do not exist
if ~isfolder(dirFigs)
    mkdir(dirFigs);
end


%% Load data
load(fname);


%% Create table for comparison
%Set up table
names={'Strategy','mode','T1','T3','DeltaT','P_el','Tsurf'};
strat={'Pel=const.';'Tsurf-Tbed=const.'};
Tbed=table('Size',[length(strat),length(names)],...
    'VariableTypes',[{'string','logical'},repmat({'double'},1,length(names)-2)],...
    'VariableNames',names);


%Fill table
Tbed.Strategy=strat;
Tbed.mode(2)=true;

Tbed.T1=arrayfun(@(tf) mean(prim.T1(prim.mode==tf)),Tbed.mode);
Tbed.T3=arrayfun(@(tf) mean(prim.T3(prim.mode==tf)),Tbed.mode);
Tbed.Tsurf=arrayfun(@(tf) mean(prim.Tsurf(prim.mode==tf)),Tbed.mode);
Tbed.P_el=arrayfun(@(tf) mean(prim.P_el(prim.mode==tf)),Tbed.mode);
Tbed.DeltaT=Tbed.T3-Tbed.T1;


%% Plot
%Set up figure
fig=figure(912);
clf(fig);
t=tiledlayout(fig,1,1,'Padding','tight');
ax=nexttile(t);
colors=ax.ColorOrder;
hold(ax,'on');


%Plot values
plot(ax,Tbed.mode,Tbed.T3,'Color',colors(2,:));
plot(ax,Tbed.mode,Tbed.T1,'Color',colors(1,:));

hold(ax,'off');


%Format axes and add legend
ax.XTick=[0,1];
grid(ax,'on');

xlabel(ax,'Mode');
ylabel(ax,'Temperature (K)');

legend(ax,{'T_3','T_1'},'Location','north');


%Size figure for repository
t.Units='centimeters';
t.OuterPosition=[0,0,17,8.5];

fig.Units=t.Units;
fig.Position(3:4)=t.OuterPosition(3:4)+0.5;


%Add left arrow and label
x=[0,0];
y=[Tbed.T1(1),Tbed.T3(1)];
figaux.arrow(ax,x,y,'doublearrow');

text(ax,x(1),mean(y),compose(' \\DeltaT = %.2f K',Tbed.DeltaT(1)));


%Add right arrow and label
x=[1,1];
y=[Tbed.T1(2),Tbed.T3(2)];
figaux.arrow(ax,x,y,'doublearrow');

text(ax,x(1),mean(y),compose('\\DeltaT = %.2f K ',Tbed.DeltaT(2)),...
    'HorizontalAlignment','right');


%Export figure
exportgraphics(fig,[dirFigs,filesep,'T1T3.tiff'],'Resolution',600);




