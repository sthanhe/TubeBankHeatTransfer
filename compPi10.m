%% Compare results to other measurements with particle cross-flow (pi10)
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
% This function compares the results from particle-convective regressions 
% in "calcPC" to other published models regarding the impact of probe size
% on the wall-to-bed HTC. 
%
%
%Requires all auxiliary classes and functions on the MATLAB path
%
%Required products, version 24.1:
%   - MATLAB
%Necessary classes, functions, files, and scripts:
%   - @figaux
%   - primData.mat --> created by the script "prepPrim"
%   - secData.mat --> created by the script "prepSec"


%% Set data locations
dirFigs=['Figures',filesep,'CFmodels'];     %Figure storage folder

fnamePrim='primData.mat';   %Primary data file
fnameSec='secData.mat';     %Secondary data file


%% Load data
%Primary data
s=load(fnamePrim);
pisOwn=s.pis;


%Secondary data: only data from Eder
s=load(fnameSec);
idx=strcmp(s.sec.Author,'Eder');
pisEder=s.pis(idx,:);


%% Total particle transport resistance
%Constants from the extended model
P4=6.45824254376922e-05;
P5=1.15226744565364;

C2=2.21729032264777;
C3=0.655360558124553;


%Individual functions
r=@(pis) (pis.pi6./pis.pi10).^(1/3)./pis.pi10;
t=@(pi5,pis) 1+0.28.*pis.pi7.^2.*sqrt(pis.pi4).*pi5.*pis.pi6;
s=@(pis) 1-exp(-P4.*pis.pi8);
d_cf=@(pi5,pis) 1+pi5.^C2.*pis.pi10.^C3.*(1-pis.pi9).^(C2*P5*3/4);

r_total=@(pi5,pis) r(pis).*t(pi5,pis).*s(pis).*d_cf(pi5,pis);


%Remove outliers
fxval=r_total(pisOwn.pi5,pisOwn);
outliers=(pisOwn.pi5<9 & fxval>19) | ...
    (pisOwn.pi5<7 & fxval>14);

pisOwn(outliers,:)=[];


%% Mean function values
%Mean pi-values for the respective source
pisMeanEder=pisEder(1,:);
pisMeanEder{:,:}=mean(pisEder{:,:},1);

pisMeanOwn=pisOwn(1,:);
pisMeanOwn{:,:}=mean(pisOwn{:,:},1);


%Set up table
authors={'Own','Eder'};
names={'Source','pi10','r','s','t','d_cf'};
tab=table('Size',[length(authors),length(names)],...
    'VariableTypes',[{'string'},repmat({'double'},1,length(names)-1)],...
    'VariableNames',names);


%Calculate mean function values
tab.Source=authors';
for i=authors
    switch i{1}
        case 'Own'
            pis=pisOwn;
        case 'Eder'
            pis=pisEder(pisEder.pi10>0,:);
    end

    idx=strcmp(tab.Source,i);

    tab.pi10(idx)=mean(pis.pi10);

    tab.r(idx)=mean(r(pis));
    tab.t(idx)=mean(t(pis.pi5,pis));
    tab.s(idx)=mean(s(pis));
    tab.d_cf(idx)=mean(d_cf(pis.pi5,pis));
end


%% Plot
%Set up figure
figidx=10;
fig=figure(figidx);
clf(fig);
til=tiledlayout(fig,1,1,'Padding','tight');
ax=nexttile(til);
colors=ax.ColorOrder;
hold(ax,'on');


%Plot data and lines
mkrSize=27;
pi5=linspace(0,30,1000)';

plot(ax,pi5,[r_total(pi5,pisMeanOwn),r_total(pi5,pisMeanEder)]);

scatter(ax,pisOwn.pi5,r_total(pisOwn.pi5,pisOwn),mkrSize,colors(1,:),'o');
scatter(ax,pisEder.pi5,r_total(pisEder.pi5,pisEder),mkrSize,colors(2,:),'x');


%Set up legend
legItems=repmat(line(ax,'Visible','off'),2,1);
legItems(1)=plot(ax,NaN,NaN,'Color',colors(2,:),...
    'Marker','x','MarkerSize',sqrt(mkrSize));
legItems(2)=plot(ax,NaN,NaN,'Color',colors(1,:),...
    'Marker','o','MarkerSize',sqrt(mkrSize));

hold(ax,'off');


%Configure axis scale
ax.YScale='log';


%Add legend and labels
lgd=legend(ax,legItems,{'Eder et al.','Test rig'},'Location','northwest');

xlabel(ax,figaux.subsz('\pi_5 (-)',6));
ylabel(ax,figaux.subsz('r_{total} (-)',6));


%Text size
ax.FontSize=7;
lgd.FontSize=7;


%Export figure for manuscript
fname=[dirFigs,filesep,'Figure',num2str(figidx)];

til.Units='centimeters';
til.OuterPosition=[0,0,9,9];

fig.Units=til.Units;
fig.Position(3:4)=til.OuterPosition(3:4)+0.5;

exportgraphics(fig,[fname,'.tiff'],'Resolution',600);


%Export figure for Elsevier
til.OuterPosition=[0,0,9,9];      %1 column
fig.Position(3:4)=til.OuterPosition(3:4)+0.5;

exportgraphics(fig,[fname,'.eps']);




