%% Demonstration of packing density model
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
% This script demonstrates how the packing density impacts the vertical
% movement of particles and leads to an "equivalent" excess fluidization
% velocity.
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

d_t=25e-3;      %Tube diameter
p_h=1.5*d_t;    %Horizontal pitch

Ar=1e4;         %Archimedes numbers
wMax=0.35;      %Maximum fluidization gas velocity
pi5eq=[30,10];  %Equivalent dimensionless excess fluidization velocities (x-axis values)

scale=0.7;  %Scaling factor for Nusselt number of tube bank


%% Heat transfer coefficients
%Derived particle and fluidization gas properties
rho_p=SiO2.rho(T);      %Particle density
rho_g=DryAir.rho(p,T);  %Fluidization gas density
rho_e=rho_p-rho_g;      %Excess particle density
c_p=SiO2.c_p(T);        %Particle specific isobaric heat capacity
k_g=DryAir.lambda(T);   %Fluidization gas thermal conductivity
my_g=DryAir.eta(T);     %Fluidization gas dynamic viscosity


%Particle diameter derived from Archimedes number
d_p=(rho_g.*(rho_p-rho_g).*FluBed.g./DryAir.eta(T).^2./Ar).^(-1/3);


%Fluidization gas velocities
w_mf=FluBed.wmfErgun(d_p,rho_p,phi_s,eps_mf,p,T);   %Minimum fluidization
w_e=repmat(linspace(0,wMax,n),length(w_mf),1);      %Excess fluidization
pi5=(rho_p.*c_p./(k_g.*FluBed.g)).^(1/3).*w_e;      %Dimensionless excess fluidization


%Heat transfer coefficient
[~,Nu]=FluBed.molWirth(w_e+w_mf,T,p,d_p,rho_p,phi_s,eps_mf,c_pfx);


%Assumed impact of tube bank
Nu_bank=Nu.total'*scale;  


%Equivalent Nusselt numbers
[~,idx]=min(abs(pi5-pi5eq'),[],2);
pi1eq=Nu.total(idx);        %Point 1 on y-axis
pi1eqBank=Nu_bank(idx(2));  %Points 2 and 3 on y-axis


%% Plot
%Set up figure
figidx=3;
fig=figure(figidx);
clf(fig);
t=tiledlayout(fig,1,1,'Padding','tight');
ax=nexttile(t);
colors=ax.ColorOrder;
hold(ax,'on');


%Plot lines
legItems=cell(2,1);

legItems{1}=plot(ax,pi5',Nu.total');
legItems{2}=plot(ax,pi5',Nu.total'*0.7,'Color',colors(1,:),'LineStyle','--');


%Plot intersection points
xline(ax,pi5eq);
scatter(ax,pi5eq,pi1eq,36,colors(2,:),'x','LineWidth',1);
scatter(ax,pi5eq(2),pi1eqBank,36,colors(2,:),'x','LineWidth',1)

hold(ax,'off');


%Legend and axis labels
legItems=[legItems{:}];
lgd=legend(ax,legItems,...
    figaux.subsz({'H_0',['H_0 \times ',num2str(scale)]},6),...
    'Location','east');

xlabel(ax,figaux.subsz('\pi_5 (-)',6));
ylabel(ax,figaux.subsz('Nu = \pi_1 (-)',6));


%Text and figure size
fsz=7;
ax.FontSize=fsz;
lgd.FontSize=7;

t.Units='centimeters';
t.OuterPosition=[0,0,9,9];

fig.Units=t.Units;
fig.Position(3:4)=t.OuterPosition(3:4)+0.5;


%Arrow style and coordinates
style={'HeadLength',5,'HeadWidth',5};
deltaY=0.7e-3;

x=[pi5eq(1),pi5eq(2),pi5eq(2)];
y=[pi1eq(1),pi1eq(2),pi1eqBank];


%Arrow from point 1 to point 2 and labels
figaux.arrow(ax,x(1:2),[y(1),y(1)]+deltaY,'arrow',style{:});

text(ax,mean(x(1:2)),y(1)+deltaY,figaux.subsz('p (\pi_9)',6),...
    'VerticalAlignment','bottom',...
    'HorizontalAlignment','center',...
    'FontSize',fsz);

text(ax,x(1),y(1)+deltaY,figaux.subsz(' \pi_5',6),...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...
    'FontSize',fsz);

text(ax,x(2),y(1)+deltaY,figaux.subsz('\pi_{5,eq} ',6),...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','right',...
    'FontSize',fsz);

drawnow();


%Arrow from point 2 to point 3 and labels
figaux.arrow(ax,x(2:3),y(2:3),'arrow',style{:});

text(ax,x(2),mean(y(2:3)),figaux.subsz(' t (\pi_5)',6),...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...
    'FontSize',fsz);


text(ax,x(1),y(1),'  1',...
    'VerticalAlignment','top',...
    'HorizontalAlignment','left',...
    'FontSize',fsz);

text(ax,x(2),y(2),'  2',...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...
    'FontSize',fsz);

t1=text(ax,x(3),y(3),'  3',...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...
    'FontSize',fsz);


%Export figure for manuscript
fname=[dirFigs,filesep,'Figure',num2str(figidx)];

exportgraphics(fig,[fname,'.tiff'],'Resolution',600);


%Export figure for Elsevier
exportgraphics(fig,[fname,'.eps']);




