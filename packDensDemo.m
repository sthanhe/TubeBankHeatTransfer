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
T=20+273.15;            %Bed temperature
rho_p=SiO2.rho(T);
rho_g=DryAir.rho(p,T);
rho_e=rho_p-rho_g;
my_g=DryAir.eta(T);
phi_s=0.8;
eps_mf=0.45;    %Porosity at minimum fluidization
c_pfx=@SiO2.c_p;

d_t=25e-3;
p_h=1.5*d_t;


%% Laminar, turbulent, and mixed heat transfer regimes
%Parameters

Ar=1e4;    %Archimedes numbers
wMax=0.35;                 %Maximum fluidization gas velocity
pi5eq=[30,10];


%Particle diameter for the given Archimedes number
d_p=(rho_g.*(rho_p-rho_g).*FluBed.g./DryAir.eta(T).^2./Ar).^(-1/3);


%Fluidization gas velocities
w_mf=FluBed.wmfErgun(d_p,rho_p,phi_s,eps_mf,p,T);   %Minimum fluidization
w_e=repmat(linspace(0,wMax,n),length(w_mf),1);      %Excess fluidization

c_p=SiO2.c_p(T);
k_g=DryAir.lambda(T);
pi5=(rho_p.*c_p./(k_g.*FluBed.g)).^(1/3).*w_e;


%Heat transfer coefficient
[~,Nu]=FluBed.molerus(w_e+w_mf,T,p,d_p,rho_p,phi_s,eps_mf,c_pfx);
Nu_bank=Nu.total'*0.7;

% [~,Nu]=FluBed.molExt(w_e+w_mf,T,p,d_p,rho_p,phi_s,eps_mf,c_pfx,...
%     d_t,Inf,0);

% [~,Nu_bank]=FluBed.molExt(w_e+w_mf,T,p,d_p,rho_p,phi_s,eps_mf,c_pfx,...
%     d_t,p_h,0);


[~,idx]=min(abs(pi5-pi5eq'),[],2);
pi1eq=Nu.total(idx);
pi1eqBank=Nu_bank(idx(2));



%%
% P1=0.0690771105844349;
% P2=18.9085028208424;
% P4=1.15226744565364;
% 
% 
% pi2=k_g./(2*c_p.*my_g);
% pi4=rho_g./rho_e;
% pi5=(rho_e.*c_p./(k_g.*FluBed.g)).^(1/3).*w_e;
% pi6=(rho_e.*c_p./(k_g.*FluBed.g)).^(1/3).*w_mf;
% pi7=1-eps_mf;
% pi9=d_t./p_h;
% 
% 
% 
% t=1+0.28.*pi7.^2.*sqrt(pi4).*pi5.*pi6;
% Nu_pcMax=P1.*pi7./(1+pi2.*t);
% 
% pfx=(1-pi9).^P4;
% d_pc=(1+P2.*(pi6./pi5).^(1/3)./pi5./pfx).^-1;
% Nu_pc=Nu_pcMax.*d_pc;


%% Create figure
figidx=3;
fig=figure(figidx);
clf(fig);
t=tiledlayout(fig,1,1,'Padding','tight');
ax=nexttile(t);
colors=ax.ColorOrder;
hold(ax,'on');

legItems=cell(2,1);

legItems{1}=plot(ax,pi5',Nu.total');
legItems{2}=plot(ax,pi5',Nu.total'*0.7,'Color',colors(1,:),'LineStyle','--');
% plot(ax,pi5',Nu_bank.total');
% plot(ax,pi5',Nu_pcMax');

xline(ax,pi5eq);
scatter(ax,pi5eq,pi1eq,36,colors(2,:),'x','LineWidth',1);
scatter(ax,pi5eq(2),pi1eqBank,36,colors(2,:),'x','LineWidth',1)


% quiver(ax,pi5eq(1),pi1eq(1),diff(pi5eq),0,...
%     'off','Color','k','MaxHeadSize',0.1/diff(pi5eq));
% ar=arrow(ax,pi5eq,[pi1eq(1),pi1eq(1)]);


hold(ax,'off');


%Legend and axis labels
legItems=[legItems{:}];
lgd=legend(ax,legItems,subsz({'H_0','H_0 \times 0.6'},6),...
    'Location','east');

xlabel(ax,subsz('\pi_5 (-)',6));
ylabel(ax,subsz('Nu = \pi_1 (-)',6));


%Text and figure size
fsz=7;
ax.FontSize=fsz;
lgd.FontSize=7;

t.Units='centimeters';
t.OuterPosition=[0,0,9,9];

fig.Units=t.Units;
fig.Position(3:4)=t.OuterPosition(3:4)+0.5;


%Add arrows
style={'HeadLength',5,'HeadWidth',5};
deltaY=0.7e-3;
% deltaY=0;

x=[pi5eq(1),pi5eq(2),pi5eq(2)];
y=[pi1eq(1),pi1eq(2),pi1eqBank];


arrow(ax,x(1:2),[y(1),y(1)]+deltaY,'arrow',style{:});

text(ax,mean(x(1:2)),y(1)+deltaY,subsz('p (\pi_9)',6),...
    'VerticalAlignment','bottom',...
    'HorizontalAlignment','center',...
    'FontSize',fsz);

text(ax,x(1),y(1)+deltaY,subsz(' \pi_5',6),...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','left',...
    'FontSize',fsz);

text(ax,x(2),y(1)+deltaY,subsz('\pi_{5,eq} ',6),...
    'VerticalAlignment','middle',...
    'HorizontalAlignment','right',...
    'FontSize',fsz);


drawnow();


arrow(ax,x(2:3),y(2:3),'arrow',style{:});

text(ax,x(2),mean(y(2:3)),subsz(' t (\pi_5)',6),...
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
fname=['Figures',filesep,'Figure',num2str(figidx)];



exportgraphics(fig,[fname,'.tiff'],'Resolution',600);


%Export figure for Elsevier
% t.OuterPosition=[0,0,9,9];
% fig.Position(3:4)=t.OuterPosition(3:4)+0.5;

exportgraphics(fig,[fname,'.eps']);


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




