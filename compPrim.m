%% Compare results to primary data
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
% All required files for this function can be found in the software
% repository: 
% https://doi.org/10.5281/zenodo.15576950
%
%
%
% This script compares the results from the chosen particle cross-flow 
% model in "calcCF" to collected primary data. Caution: Make sure the
% coefficients (regression results) of the chosen model are implemented in
% the "extended" function of the "FluBed" class (by default, the
% coefficients from the main paper are implemented)
%
%
%Requires all auxiliary classes and functions on the MATLAB path
%
%Required products, version 24.1:
%   - MATLAB
%   - Curve Fitting Toolbox
%Necessary classes, functions, files, and scripts:
%   - @DryAir
%   - @FluBed
%   - @SiO2
%   - @figaux
%   - @implExp
%   - getConstants.m
%   - primData.mat --> created by the script "prepPrim"


%% Set data locations
dirFigs=['Figures',filesep,'CFmodels'];     %Figure storage folder

fname='primData.mat';   %Primary data file


%% Make folders if they do not exist
if ~isfolder(dirFigs)
    mkdir(dirFigs);
end


%% Load data
load(fname);


%Get constants
c=getConstants();


%% Response variables
%Particle convective and cross-flow Nu according to the extended model
[~,Nu_ext]=FluBed.extended(prim.w,prim.Tbed,prim.p,c.d_p,c.rho_p,c.phi_s,c.eps_mf,...
    @SiO2.c_p,c.d_t,c.p_hEff,prim.w_p);

yEst_ext=Nu_ext.cf+Nu_ext.pc;


%Particle convective Nu according to Molerus / Wirth
[~,Nu_MW]=FluBed.molWirth(prim.w,prim.Tbed,prim.p,c.d_p,c.rho_p,c.phi_s,c.eps_mf,...
    @SiO2.c_p);

yEst_MW=Nu_MW.pcMix;


%Measured effective Nu minus gas convection
y=prim.Nu_eff-Nu_ext.gc;


%Mean values of all estimates for each model
Nu_extMean=mean([Nu_ext.pc,Nu_ext.cf],1);   %External model
Nu_MWmean=[mean([Nu_MW.pcMix],1),0];        %Molerus / Wirht


%% Set up figure
figidx=7;
fig=figure(figidx);
clf(fig);

til=tiledlayout(fig,1,2,'Padding','tight');
til.TileIndexing='columnmajor';
til.TileSpacing='none';

ax=cell(1,2);


%% Left side: comparison between measurements and estimates
ax{1}=nexttile(til);
colors=ax{1}.ColorOrder;
hold(ax{1},'on');


%Plot data
legItems=cell(1,2);     %Legend item container

legItems{1}=scatter(ax{1},y,yEst_MW,18,colors(1,:),'o');
legItems{2}=scatter(ax{1},y,yEst_ext,18,colors(2,:),'+');


%Plot equivalence lines
lim=max([ax{1}.XLim(2),ax{1}.YLim(2)]);
eq=linspace(0,lim,100);

plot(ax{1},eq,eq,'Color','k');
plot(ax{1},eq,eq.*1.2,'Color','k','LineStyle','--');
plot(ax{1},eq,eq./1.2,'Color','k','LineStyle','--');

hold(ax{1},'off');


%Set legend
txt=figaux.subsz(...
    compose('%s (H_%d)',...
    ["Molerus / Wirth";"Extended Model"],...
    [0;1]),6);
lgd=legend(ax{1},[legItems{:}],...
    txt,...
    'Location','southeast',...
    'FontSize',7);


%Axes limits and labels
ax{1}.XLim=[0,lim];
ax{1}.YLim=[0,lim];

xlabel(ax{1},figaux.subsz('Measured Nu - Nu_{gc} (-)',6));
ylabel(ax{1},figaux.subsz('Estimated Nu - Nu_{gc} (-)',6));


%% Right side: composition of Nu numbers
ax{2}=nexttile(til);
colors=ax{2}.ColorOrder;


%Plot data
bar(ax{2},1:2,[Nu_MWmean;Nu_extMean],'stacked');


%Set legend
legend(ax{2},figaux.subsz({'Nu_{pc}','Nu_{cf}'},6),...
    'FontSize',7);


%Axes limits and labels
ax{2}.YLim=ax{1}.YLim;

txt=figaux.subsz(...
    compose('%s (H_%d)',...
    ["Molerus /\newlineWirth";"Extended\newlineModel"],...
    [0;1]),6);
ax{2}.XTickLabel=txt;


%Axes appearance
ax{2}.YTickLabel={};
ax{2}.YTick=[];

ax{2}.Box='off';


%Add descriptive texts
text(ax{2},1,sum(Nu_MWmean),...
    {'Particle Convection','Only'},...
    'HorizontalAlignment','center',...
    'VerticalAlignment','bottom',...
    'FontSize',7);

text(ax{2},2,sum(Nu_extMean),...
    {'Particle Cross-Flow','+ Particle Convection'},...
    'HorizontalAlignment','center',...
    'VerticalAlignment','bottom',...
    'FontSize',7);


%% Axes configuration
ax=[ax{:}];


%Axes limits
set(ax,'YLim',[0,lim]);
linkaxes(ax,'y');


%Text size
set(ax,'FontSize',7);


%Export figure for manuscript
fname=[dirFigs,filesep,'Figure',num2str(figidx)];

til.Units='centimeters';
til.InnerPosition=[1.5,1,14,7];

fig.Units=til.Units;
fig.Position(3:4)=til.OuterPosition(3:4)+0.5;

exportgraphics(fig,[fname,'.tiff'],'Resolution',600);


%Export figure for Elsevier
til.InnerPosition=[1.5,1,14,7];
fig.Position(3:4)=til.OuterPosition(3:4)+0.5;

exportgraphics(fig,[fname,'.eps']);
savefig(fig,fname);




