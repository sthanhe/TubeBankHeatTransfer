


%% Pi9 analysis
p=1e5;
T=293.15;

Ar=[1e3,1e4];
% Ar=1e4;
eps_mf=0.45;
phi_s=0.8;

d_t=25e-3;
w_p=0;
c_pfx=@SiO2.c_p;

s_h=linspace(1,0.95,100);
s_h=[s_h,linspace(0.95,sqrt(2)/2+1e-3,100)];
s_h=[s_h,linspace(sqrt(2)/2+1e-3,sqrt(2)/2-1e-2,1000)];
s_h=[s_h,linspace(sqrt(2)/2-1e-2,0,100)];
p_h=d_t./s_h;


%Gas and particle properties
rho_p=SiO2.rho(T);
rho_g=DryAir.rho(p,T);

d_p=(rho_g.*(rho_p-rho_g).*FluBed.g./DryAir.eta(T).^2./Ar).^(-1/3);

w_mf=FluBed.wmfErgun(d_p,rho_p,phi_s,eps_mf,p,T);

w=arrayfun(@(w_mf) ...
    linspace(w_mf,20*w_mf,1000)',...
    w_mf,'UniformOutput',false);
w=horzcat(w{:});

p_h=reshape(p_h,[1,1,numel(p_h)]);


[~,Nu]=FluBed.molExt(w,T,p,d_p,rho_p,phi_s,eps_mf,c_pfx,...
                d_t,p_h,w_p);

Nu_max=max(Nu.total,[],1,'omitmissing');
Nu_max=reshape(Nu_max,[length(Ar),length(p_h)]);

Nu_relExt=Nu_max./Nu_max(:,end);

p_h=reshape(p_h,[1,numel(p_h)]);

% plot(s_h,Nu_relExt);
% a=0;


% Grewal / Saxena 1983, p. 371
Nu_relGW=1-0.21.*(p_h./d_t).^-1.75;
Nu_relGW=Nu_relGW./Nu_relGW(:,end);


%Lechner, p. 17
% p_diag=sqrt((p_h./2).^2+p_h.^2);
% p_min=p_diag-d_t;
% Nu_relLechner=(1-d_t./p_h).^0.36.*(1-d_t./p_diag).^0.24.*(1-d_p'./p_min).^4;
% Nu_relLechner=Nu_relLechner./Nu_relLechner(:,end);


%Gelperin / Einstein
Nu_relGE=(1-(d_t./p_h).*(1+d_t./(d_t+p_h))).^0.25;
Nu_relGE(Nu_relGE~=real(Nu_relGE))=NaN;
Nu_relGE=Nu_relGE./Nu_relGE(:,end);

idx=isnan(Nu_relGE);
x_GE=s_h;
x_GE(idx)=NaN;

idx=find(idx,1,"last");
l=20;
x_GE(idx-l:idx)=linspace(sqrt(2)/2,x_GE(idx+1),l+1);
Nu_relGE(idx-l:idx)=linspace(0,Nu_relGE(idx+1),l+1);


%Natusch, see Hofer (2018), p. 161
Nu_relNatusch=(1-d_t./p_h).^0.25;
Nu_relNatusch=Nu_relNatusch./Nu_relNatusch(:,end);



%Plot
figidx=9;
fig=figure(figidx);
clf(fig);
t=tiledlayout(fig,1,1,'Padding','tight');
ax=nexttile(t);
hold(ax,'on');


plot(ax,s_h,Nu_relExt);
plot(ax,s_h,[Nu_relGW;Nu_relNatusch],'LineStyle','--');
plot(ax,x_GE,Nu_relGE,'LineStyle','--');

hold(ax,'off');


lgd=legend(ax,[compose('Extended Model, Ar = 1e%d',log10(Ar'));...
    {'Grewal & Saxena';...
        'Natusch et al.'};...
        'Gel''perin et al.'],...
    'Location','southwest');

xlabel(ax,subsz('\pi_9 = d_t / p_h (-)',6));
ylabel(ax,subsz('Nu_{max} / Nu_{max} (p_h \rightarrow \infty)',6));


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




