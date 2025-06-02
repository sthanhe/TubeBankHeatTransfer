%% Molerus and Wirth heat transfer demonstration
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
% This script demonstrates key features of the heat transfer correlation by
% Molerus and Wirth and creates the corresponding figure in the main paper. 
% It is largely a copy of a similar analysis from a previous study: 
% https://doi.org/10.5281/ZENODO.10207330
%
%
%Requires all auxiliary classes and functions on the MATLAB path
%
%Required products, version 24.1:
%   - MATLAB
%   - Curve Fitting Toolbox
%Necessary files, classes, functions, and scripts:
%   - @DryAir
%   - @FluBed
%   - @SiO2
%   - @implExp


%% Set data locations
dirFigs='Figures';  %Figure storage folder


%% Make folders if they do not exist
if ~isfolder(dirFigs)
    mkdir(dirFigs);
end


%% General parameters
n=1000;     %Number of cells for discretization

p=1e5;              %Bed pressure
T=20+273.15;        %Bed temperature
phi_s=0.8;          %Particle sphericity
eps_mf=0.45;        %Porosity at minimum fluidization
c_pfx=@SiO2.c_p;    %Specific heat capacity function: silica

Ar=logspace(2,6,3)';    %Archimedes numbers
wMax=2;                 %Maximum fluidization gas velocity


%% Heat transfer coefficients
%Derived particle and fluidization gas properties
rho_p=SiO2.rho(T);      %Particle density
rho_g=DryAir.rho(p,T);  %Fluidization gas density
c_p=SiO2.c_p(T);        %Particle specific isobaric heat capacity
k_g=DryAir.lambda(T);   %Fluidization gas thermal conductivity


%Particle diameter derived from Archimedes number
d_p=(rho_g.*(rho_p-rho_g).*FluBed.g./DryAir.eta(T).^2./Ar).^(-1/3);
d_p(1)=d_p(1)-1e-7;     %Fix rounding issue to ensure laminar regime


%Fluidization gas velocities
w_mf=FluBed.wmfErgun(d_p,rho_p,phi_s,eps_mf,p,T);   %Minimum fluidization
w_e=repmat(linspace(0,wMax,n),length(w_mf),1);      %Excess fluidization
pi5=(rho_p.*c_p./(k_g.*FluBed.g)).^(1/3).*w_e;      %Dimensionless excess fluidization


%Heat transfer coefficient
[~,Nu]=FluBed.molWirth(w_e+w_mf,T,p,d_p,rho_p,phi_s,eps_mf,c_pfx);


%% Plot
%Set up figure
figidx=2;
fig=figure(figidx);
clf(fig);
t=tiledlayout(fig,1,1,'Padding','tight');
ax=nexttile(t);
colors=ax.ColorOrder;
hold(ax,'on');


%Plot lines
plot(ax,pi5',Nu.total');
yline(ax,Nu.pcMax(1),'--','Color',colors(1,:));
yline(ax,Nu.gcMax(1),'--','Color',colors(3,:));

hold(ax,'off');


%Legend and axis labels
txt=compose('Ar=10^%.0f, ',log10(Ar));
txt=strcat(txt,{'laminar';'mixed';'turbulent'});
txt=figaux.subsz([txt;{'Nu_{max,lam}';'Nu_{max,turb}'}],6);
lgd=legend(ax,txt,'Location','east');

xlabel(ax,figaux.subsz('\pi_5 (-)',6));
ylabel(ax,figaux.subsz('Nu = \pi_1 (-)',6));


%Text size
ax.FontSize=7;
lgd.FontSize=7;


%Export figure for manuscript
fname=[dirFigs,filesep,'Figure',num2str(figidx)];

t.Units='centimeters';
t.OuterPosition=[0,0,9,9];

fig.Units=t.Units;
fig.Position(3:4)=t.OuterPosition(3:4)+0.5;

exportgraphics(fig,[fname,'.tiff'],'Resolution',600);


%Export figure for Elsevier
t.OuterPosition=[0,0,9,9];
fig.Position(3:4)=t.OuterPosition(3:4)+0.5;

exportgraphics(fig,[fname,'.eps']);




