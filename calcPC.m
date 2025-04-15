fname='extLit.mat';


%% Load data
load(fname);


%Create response variable y and regressor matrix X: only data from Grewal
idx=strcmp(tab.Author,'Grewal');
y=pis.pi1(idx)-tab.Nu_gcMol(idx);
X=pis{idx,:};


%% Null model
fx0=@(b,x) ...
    b(1).*x(:,7)./...
    (1+x(:,2).*...
    (1-exp(-b(4).*x(:,8))).*...
    (1+b(2).*x(:,7).^2.*sqrt(x(:,4)).*x(:,5).*x(:,6))).*...
    (1+b(3).*(x(:,6)./x(:,5)).^(1/3)./x(:,5)./(1-x(:,9)).^b(5)).^-1;


beta_s0=[0.125,0.28,33.3,Inf,0]';

checkFit(X,y,fx0,beta_s0,'s0',600);
plotFit(X,y,fx0,beta_s0,['Figures',filesep,'Table6_s0'],601);


%Coefficient of determination
y_mean=mean(y);
y_est=fx0(beta_s0,X);
res=y-y_est;

SSE=sum(res.^2);
SST=sum((y-y_mean).^2);
R2=1-SSE./SST;

n=numel(y);
M=9;    %Number of regressors excluding intercept
R2_adj=1-(n-1)./(n-M).*(1-R2);




%% Model s1
fx1=@(b,x) ...
    b(1).*x(:,7)./...
    (1+x(:,2).*...
    (1-(1+b(3).*x(:,8).^b(5)).^-1).*...
    (1+0.28.*x(:,7).^2.*sqrt(x(:,4)).*x(:,5).*x(:,6))).*...
    (1+b(2).*(x(:,6)./x(:,5)).^(1/3)./x(:,5)./(1-x(:,9)).^b(4)).^-1;


beta0=[0.125,33.3,5e-8,1,2]';

[mdl_s1,beta_s1]=getFit(X,y,fx1,beta0);

checkFit(X,y,fx1,beta_s1,'s1',610);
plotFit(X,y,fx1,beta_s1,['Figures',filesep,'Table6_s1'],611);


%% Model s2
% factor 1/40 the expected value
fx2=@(b,x) ...
    b(1).*x(:,7)./...
    (1+x(:,2).*...
    (1-exp(-b(3).*x(:,8))).*...
    (1+0.28.*x(:,7).^2.*sqrt(x(:,4)).*x(:,5).*x(:,6))).*...
    (1+b(2).*(x(:,6)./x(:,5)).^(1/3)./x(:,5)./(1-x(:,9)).^b(4)).^-1;


beta0=[0.125,33.3,5e-5,1]';

[mdl_s2,beta_s2]=getFit(X,y,fx2,beta0);

checkFit(X,y,fx2,beta_s2,'s2',620);
plotFit(X,y,fx2,beta_s2,['Figures',filesep,'Table6_s2'],621);


%% Model s3
% factor still 1/10 the
%expected value. Slightly worse R²_adj, but fewer outliers. b(3)
%practically identical to Molerus (0.28). AIC=-148.9956
fx3=@(b,x) ...
    b(1).*x(:,7)./...
    (1+x(:,2).*...
    tanh(b(3).*x(:,8)).*...
    (1+0.28.*x(:,7).^2.*sqrt(x(:,4)).*x(:,5).*x(:,6))).*...
    (1+b(2).*(x(:,6)./x(:,5)).^(1/3)./x(:,5)./(1-x(:,9)).^b(4)).^-1;


beta0=[0.125,33.3,5e-5,1]';

[mdl_s3,beta_s3]=getFit(X,y,fx3,beta0);

checkFit(X,y,fx3,beta_s3,'s3',630);
plotFit(X,y,fx3,beta_s3,['Figures',filesep,'Table6_s3'],631);


%% s2 compared to Molerus
pidx=[1,3:5];
names={'Parameter','H0','Estimate','p'};
para=table('Size',[length(pidx)+1,length(names)],...
    'VariableTypes',[{'string'},repmat({'double'},1,length(names)-1)],...
    'VariableNames',names);


para.Parameter=[compose('P%d',pidx),{'All'}]';
para.H0(1:end-1)=beta_s0(pidx);
para.Estimate(1:end-1)=beta_s2;


for i=1:height(para)-1
    beta_null=para.Estimate;
    beta_null(i)=para.H0(i);

    [~,p]=ttest(fx2(beta_null,X),y);
    para.p(i)=p;
end


beta_null=para.H0;
[~,para.p(end)]=ttest(fx2(beta_null,X),y);


%% Plot graphic
figidx=5;
fig=figure(figidx);
clf(fig);
t=tiledlayout(fig,1,1,'Padding','tight');
ax=nexttile(t);
colors=ax.ColorOrder;
hold(ax,'on');


legItems=cell(1,2);

legItems{1}=scatter(ax,y,fx0(beta_s0,X),18,'o');
legItems{2}=scatter(ax,y,fx2(beta_s2,X),18,'+');

lim=max([ax.XLim(2),ax.YLim(2)]);
eq=linspace(0,lim,100);

plot(ax,eq,eq,'Color','k');
plot(ax,eq,eq.*1.2,'Color','k','LineStyle','--');
plot(ax,eq,eq./1.2,'Color','k','LineStyle','--');

hold(ax,'off');


legItems=[legItems{:}];
txt=subsz(compose('H_%d, R^2_{adj}=%.3f',...
    [0;1],[R2_adj;mdl_s2.Rsquared.Adjusted]),6);
lgd=legend(ax,legItems,txt,'Location','southeast');


ax.XLim=[0,lim];
ax.YLim=[0,lim];

xlabel(ax,subsz('Measured Nu_{pc} (-)',6));
ylabel(ax,subsz('Estimated Nu_{pc} (-)',6));


%Text size
ax.FontSize=7;
lgd.FontSize=7;


%Export figure for manuscript
fname=['Figures',filesep,'Figure',num2str(figidx)];

t.Units='centimeters';
t.InnerPosition=[1.5,1,7.4,7.4];

fig.Units=t.Units;
fig.Position(3:4)=t.OuterPosition(3:4)+0.5;

exportgraphics(fig,[fname,'.tiff'],'Resolution',600);


%Export figure for Elsevier
t.InnerPosition=[1.5,1,7.4,7.4];
fig.Position(3:4)=t.OuterPosition(3:4)+0.5;

exportgraphics(fig,[fname,'.eps']);


% %%
% idx=strcmp(tab.Author,'Grewal') & strcmp(tab.Material,'Silica');
% tab2=tab(idx,:);
% pis=pis(idx,:);
% 
% 
% p=mean(tab2.p);     %constant
% T=mean(tab2.T);     %constant
% 
% my_g=DryAir.eta(T);         %constant
% rho_g=DryAir.rho(p,T);      %constant
% rho_p=mean(tab2.rho_p);     %constant
% rho_e=rho_p-rho_g;          %constant
% c_p=mean(tab2.c_p);         %constant
% c_g=mean(tab2.c_g);         %constant
% l_l=(my_g./(rho_e.*sqrt(FluBed.g))).^(2/3);     %constant
% 
% pi2=mean(pis.pi2,'omitmissing');    %constant
% pi4=mean(pis.pi4,'omitmissing');    %constant
% pi5=mean(pis.pi5,'omitmissing');    %varying
% pi6=mean(pis.pi6,'omitmissing');    %varying
% pi7=mean(pis.pi7,'omitmissing');    %varying
% 
% d_t=[linspace(0,40e-3,100),Inf];
% pi8=d_t./l_l;   %varying with d_t only
% 
% d_p=mean(tab2.d_p);     %varying
% Ar=rho_g.*d_p.^3.*rho_e.*FluBed.g./my_g.^2;
% 
% 
% % my_g=1.96e-5;
% % k_g=0.149;
% % c_p=1000;
% % rho_g=0.1785;
% % rho_p=1000;
% % rho_e=rho_p-rho_g;
% % d_p=103e-6;
% % pi2=k_g./(2*c_p.*my_g);
% 
% 
% s1=1-(1+beta_s1(3).*pi8.^beta_s1(5)).^-1;
% s2=1-exp(-beta_s2(3).*pi8);
% s3=tanh(beta_s3(3).*pi8);
% 
% Nu_max1=beta_s1(1).*pi7./...
%     (1+pi2.*s1.*...
%     (1+0.28.*pi7.^2.*sqrt(pi4).*pi5.*pi6));
% 
% Nu_max2=beta_s2(1).*pi7./...
%     (1+pi2.*s2.*...
%     (1+0.28.*pi7.^2.*sqrt(pi4).*pi5.*pi6));
% 
% Nu_max3=beta_s3(1).*pi7./...
%     (1+pi2.*s3.*...
%     (1+0.28.*pi7.^2.*sqrt(pi4).*pi5.*pi6));
% 
% 
% Nu_rel1=Nu_max1./Nu_max1(end);
% Nu_rel2=Nu_max2./Nu_max2(end);
% Nu_rel3=Nu_max3./Nu_max3(end);
% 
% 
% %Molerus / Wirth p. 17
% C=0.85;
% u_l=0.3e-2;
% f_L=1;
% M=4.25e-3;
% 
% Nu_relMW=C.*u_l./(f_L.*d_t)+1;
% 
% 
% % Grewal / Saxena 1981, p. 113
% d_t127=12.7e-3;
% Nu_wpmax=0.9.*(Ar.*d_t127./d_t).^0.21.*(c_p./c_g).^(45.5.*Ar.^-0.7);
% Nu_relGW=Nu_wpmax./Nu_wpmax(end-1);
% 
% 
% 
% %Plot
% fig=figure(400);
% clf(fig);
% ax=gca();
% hold(ax,'on');
% 
% plot(ax,d_t.*10^3,[Nu_rel1;Nu_rel2;Nu_rel3]);
% plot(ax,d_t.*10^3,Nu_relMW,'Color','k','LineStyle','--');
% plot(ax,d_t.*10^3,Nu_relGW,'Color','k','LineStyle',':');
% 
% hold(ax,'off');
% 
% ax.YLim=[1,5];
% 
% legend(ax,[compose('s_%d',1:3),...
%     {'Molerus&Wirth','Grewal&Saxena'}],'Location','best');
% 
% xlabel(ax,'Tube diameter (mm)');
% ylabel(ax,'Nu_{max,pc} / Nu_{max,pc}(d_t\rightarrow\infty)');
% 
% fig.Units='centimeters';
% fig.Position=[10,5,17,8.5];
% 
% % exportgraphics(fig,['Figures',filesep,'Figure3.tiff'],'Resolution',600);

