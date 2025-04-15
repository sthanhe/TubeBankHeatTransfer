


%% Pi8 analysis
p=1e5;
T=293.15;

Ar=[1e3,1e4];
% Ar=1e4;
eps_mf=0.45;
phi_s=0.8;

p_h=Inf;
w_p=0;
c_pfx=@SiO2.c_p;

d_t=linspace(0,5e-3,100);
d_t=[d_t,linspace(d_t(end),60e-3,100),Inf];


%Gas and particle properties
rho_p=SiO2.rho(T);
rho_g=DryAir.rho(p,T);

d_p=(rho_g.*(rho_p-rho_g).*FluBed.g./DryAir.eta(T).^2./Ar).^(-1/3);

w_mf=FluBed.wmfErgun(d_p,rho_p,phi_s,eps_mf,p,T);

w=arrayfun(@(w_mf) ...
    linspace(w_mf,60*w_mf,10000)',...
    w_mf,'UniformOutput',false);
w=horzcat(w{:});

d_t=reshape(d_t,[1,1,numel(d_t)]);


[~,Nu]=FluBed.molExt(w,T,p,d_p,rho_p,phi_s,eps_mf,c_pfx,...
                d_t,p_h,w_p);

Nu_max=max(Nu.total,[],1,'omitmissing');
Nu_max=reshape(Nu_max,[length(Ar),length(d_t)]);

Nu_relExt=Nu_max./Nu_max(:,end);

d_t=reshape(d_t,[1,numel(d_t)]);


%Molerus / Wirth p. 17
C=0.85;
u_l=0.3e-2;
f_L=1;

Nu_relMW=C.*u_l./(f_L.*d_t)+1;


% Grewal / Saxena 1981, p. 113
d_tMaxGW=28.6e-3;
Nu_relGW=d_t.^-0.21./d_tMaxGW.^-0.21;


%Merzsch, p. 1043
d_tMaxMerzsch=33.7e-3; 
Nu_relMerzsch=d_t.^-0.3./d_tMaxMerzsch.^-0.3;



%Plot
figidx=8;
fig=figure(figidx);
clf(fig);
t=tiledlayout(fig,1,1,'Padding','tight');
ax=nexttile(t);
hold(ax,'on');


plot(ax,d_t.*10^3,Nu_relExt);
plot(ax,d_t.*10^3,[Nu_relMW;Nu_relGW;Nu_relMerzsch],'LineStyle','--');

hold(ax,'off');


ax.XLim=[0,40];
ax.YLim=[1,3.5];


lgd=legend(ax,subsz(...
    [compose('Extended Model, Ar=1e%d',log10(Ar'));...
    {'Molerus & Wirth';...
        ['Grewal & Saxena, d_{t,max} = ',num2str(d_tMaxGW*10^3), ' mm'];...
        ['Merzsch et al, d_{t,max} = ',num2str(d_tMaxMerzsch*10^3), ' mm']}],...
        6),...
    'Location','north');

xlabel(ax,subsz('Tube diameter d_t (mm)',6));
ylabel(ax,subsz('Nu_{max} / Nu_{max} (d_t \rightarrow \infty)',6));


%Text size
ax.FontSize=7;
lgd.FontSize=7;


%Export figure for manuscript
fname=['Figures',filesep,'Figure',num2str(figidx)];

t.Units='centimeters';
t.OuterPosition=[0,0,9,9];

fig.Units=t.Units;
fig.Position(3:4)=t.OuterPosition(3:4)+0.5;

exportgraphics(fig,[fname,'.tiff'],'Resolution',600);


%Export figure for Elsevier
t.OuterPosition=[0,0,9,9];
fig.Position(3:4)=t.OuterPosition(3:4)+0.5;

exportgraphics(fig,[fname,'.eps']);



