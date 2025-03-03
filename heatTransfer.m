%% Molerus heat transfer demonstration
%GNU General Public License v3.0
%By Stefan Thanheiser: https://orcid.org/0000-0003-2765-1156
%
%Part of the sandTES Engineering Manual
%
%All required files for this script can be found in the software
%repository:
%https://doi.org/10.5281/ZENODO.10207330
% 
%All parameters and results are in SI base units.
%
%
%
%This script demonstrates key features of the heat transfer correlation by
%Molerus and creates the Figures 13 and 14 in Section 3.4.4 of the 
%Engineering Manual.
%
%
%Requires all auxiliary classes and functions on the MATLAB path
%
%Required products:
%   - MATLAB, version 9.14
%   - Curve Fitting Toolbox, version 3.9
%Necessary files, classes, functions, and scripts:
%   - @DryAir
%   - @FluBed
%   - @SiO2
%   - @Sinter
%   - @implExp


%% General parameters
n=1000;     %Number of cells for discretization

p=1e5;          %Bed pressure
rho_p=2650;     %Particle density
phi_s=0.8;
eps_mf=0.45;    %Porosity at minimum fluidization
c_pfx=@SiO2.c_p;


%% Laminar, turbulent, and mixed heat transfer regimes
%Parameters
T=20+273.15;            %Bed temperature
Ar=logspace(2,6,3)';    %Archimedes numbers
wMax=2;                 %Maximum fluidization gas velocity


%Calculate particle diameters for the given Archimedes numbers
d_p=arrayfun(@(i) ...
        fzero(@(d_p) ...
            FluBed.Ar(d_p,rho_p,p,T)-Ar(i),...
        [50,5000].*10^-6),...
    1:length(Ar))';


%Fluidization gas velocities
w_mf=FluBed.wmfErgun(d_p,rho_p,phi_s,eps_mf,p,T);   %Minimum fluidization
w_e=repmat(linspace(0,wMax,n),length(w_mf),1);      %Excess fluidization

c_p=SiO2.c_p(T);
k_g=DryAir.lambda(T);
pi5=(rho_p.*c_p./(k_g.*FluBed.g)).^(1/3).*w_e;


%Heat transfer coefficient
[h,Nu]=FluBed.molerus(w_e+w_mf,T,p,d_p,rho_p,phi_s,eps_mf,c_pfx);


%Create figure
figidx=1;
fig=figure(figidx);
clf(fig);
ax=gca();
colors=ax.ColorOrder;
hold(ax,'on');

plot(ax,pi5',Nu.total');
yline(ax,Nu.pcMax(1),'--','Color',colors(1,:));
yline(ax,Nu.gcMax(1),'--','Color',colors(3,:));

hold(ax,'off');

legend(ax,[compose('Ar=10^%.0f',log10(Ar));{'Nu_{max,lam}';'Nu_{max,turb}'}],'Location','east');

xlabel(ax,'\pi_5 (-)');
ylabel(ax,'\pi_1 (-)');

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];

exportgraphics(fig,['Figures',filesep,'Figure1.tiff'],...
    'Resolution',600);


% %% Impact of bed temperature
% %Parameters
% T=[20,400]'+273.15;     %Bed temperatures
% d_p=400e-6;             %Particle diameter
% FG=linspace(1,10,n);    %Degrees of fluidization
% 
% 
% %Fluidization gas velocities
% w_mf=FluBed.wmf(d_p,rho_p,p,T);     %Minimum fluidization
% w=FG.*w_mf;                         %Actual fluidization velocity
% 
% 
% %Heat transfer coefficient
% h=FluBed.molerus(w,T,T,p,d_p,rho_p,eps_mf);
% Ar=FluBed.Ar(d_p,rho_p,p,T);                    %Archimedes number
% 
% 
% %Create figure
% fig=figure(14);
% clf(fig);
% ax=gca();
% colors=ax.ColorOrder;
% hold(ax,'on');
% 
% plot(ax,FG,h.total(2,:)','Color',colors(2,:));
% plot(ax,FG,h.total(1,:)','Color',colors(1,:));
% 
% hold(ax,'off');
% 
% 
% legend(ax,compose('T=%.0f°C',flipud(T-273.15)),'Location','best');
% 
% xlabel(ax,'Degree of fluidization FG (-)');
% ylabel(ax,'Heat transfer coefficient h (W/m²K)');
% 
% title(ax,sprintf('p=%.0f bar, d_p=%.0f µm, \\rho_p=%.0f kg/m³, \\epsilon_{mf}=%.2f',...
%                 p.*10^-5,d_p.*10^6,rho_p,eps_mf));
% 
% fig.Units='centimeters';
% fig.Position=[10,5,17,8.5];
% 
% exportgraphics(fig,['Figures',filesep,'HTCtemp.tiff']);




