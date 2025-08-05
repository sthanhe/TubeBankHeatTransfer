%% Analysis of particle size distribution
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
% All required files for this script can be found in the software
% repository: 
% https://doi.org/10.5281/zenodo.15576950
%
%
%
% This script calculates the mean particle diameter from the data in the 
% particle supplier's data sheet. It is largely a copy of a script with the
% same name from a previous study: https://doi.org/10.5281/zenodo.7948224
%
%
%Requires all auxiliary classes and functions on the MATLAB path
%
%Required products, version 24.1:
%   - MATLAB
%Necessary classes, functions, files, and scripts:
%   - None


%% Set data locations
dirFigs='Figures';  %Figure storage folder


%% Make folders if they do not exist
if ~isfolder(dirFigs)
    mkdir(dirFigs);
end


%% Analysis
mesh=[425,300,212,150,106,75,53,0]*10^-6;       %Mesh sizes
resid=[0,2.2,14.7,47.5,28.8,6.4,0.4,0]/100;     %Residues in each pan

meshmean=movmean(mesh,[1,0]);   %Mean particle size between each mesh, assuming linear distribution

d_p=sum(resid.*meshmean);    %Mean particle diameter


%Set up figure
fig=figure(908);
clf(fig);
t=tiledlayout(fig,1,1,'Padding','tight');
ax=nexttile(t);
box(ax,'on');


%Plot lines
plot(ax,[mesh;meshmean],resid);


%Set axis limits
ax.XLim=[0,max(mesh)];


%Set legend, axis labels, and title
legend(ax,{'Sieve','Linear'},'Location','best');

xlabel(ax,'Mesh size (m)');
ylabel(ax,'Retained mass fraction (-)');

title(ax,'Particle size distribution, GRANUSIL, Wedron IL #801, Grade 7020');


%Export figure for repository
t.Units='centimeters';
t.OuterPosition=[0,0,17,8.5];

fig.Units=t.Units;
fig.Position(3:4)=t.OuterPosition(3:4)+0.5;

exportgraphics(fig,[dirFigs,filesep,'particleSize.tiff'],...
    'Resolution',600);




