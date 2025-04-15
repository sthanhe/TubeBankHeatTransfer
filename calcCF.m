%% Set data directories
dirData=['Data',filesep,'Own'];             %Path to directory containing the data
dirFigures='Figures';       %Path to directory where figures should be stored

if ~isfolder(dirFigures)
    mkdir(dirFigures);
end


%% Prepare analysis
%Get constants
c=getConstants();


%Retrieve filenames
files=dir(dirData);
files=files(startsWith({files.name},'heatTransfer_'));

fnames={files.name}';


%Retain natural order of runs
idx=str2double(extract(fnames,digitsPattern));
[~,idx]=sort(idx);
fnames=fnames(idx);


%Set up table for mean values
names={'Run','p','Tsurf','Tbed','w','w_p','hVirt','P','mode'};
htc=table('Size',[length(fnames),length(names)],'VariableTypes',[repmat({'double'},1,length(names)-1),'logical']);
htc.Properties.VariableNames=names;
clear('names');


%% Read individual files and do calculations
chambers=1:6;
for i=1:length(fnames)
    %Get properties
    tab=readtable([dirData,filesep,fnames{i}]);
    TT=getProp(tab,c,htc.Properties.VariableNames(2:end),chambers);
    
    
    %Remove outliers
    outliers=isoutlier(TT.hVirt,1);
    TT{outliers,2:end-1}=NaN;


    %Record table for future analysis
    % writetable(TT,[dirStationary,filesep,'stat_Run',num2str(i),'.csv']);
    
    
    %Get means
    htc{i,2:width(TT)-1}=mean(TT{:,2:end-1},1,'omitnan');
    htc.Run(i)=i;
    htc.mode(i)=nnz(TT.mode)>height(htc)/2;
end


%Record table for future analysis
htc(htc.hVirt<5,:)=[];
% writetable(htc,[dirData,filesep,'htc_Sum.csv']);


%% Effective HTC
D=c.d_t+2*c.h_f;                            %Outside tube diameter including fins
phi=(D./c.d_t-1).*(1+0.35*log(D./c.d_t));   %Geometry function


htc.h_eff=htc.hVirt;
htc.eta_f=ones(height(htc),1);
lambda=55;
for i=1:height(htc)
    err=1;
    counter=0;
    while err>1e-6 && counter<100
        X=phi.*c.d_t/2.*sqrt(2*htc.h_eff(i)./(lambda*c.s_f));
        htc.eta_f(i)=tanh(X)./X;
        A_eff=c.A_plain-c.A_bottom+htc.eta_f(i).*c.A_sides;

        deltaT=htc.Tsurf(i)-htc.Tbed(i);
        h_effNew=htc.P(i)./(A_eff.*deltaT);
        % h_effNew=htc.hVirt(i).*c.A_plain./A_eff;

        Tfin=htc.eta_f(i).*deltaT+htc.Tbed(i);
        lambda=DC04.lambda(Tfin);

        err=abs(htc.h_eff(i)-h_effNew);
        htc.h_eff(i)=h_effNew;

        counter=counter+1;
    end
end


%% Estimate bias due to fins
[h,Nu]=FluBed.molExt(htc.w,htc.Tbed,htc.p,c.d_p,c.rho_p,c.phi_s,c.eps_mf,...
    @SiO2.c_p,c.d_t,c.p_hEff,htc.w_p);

htc.h_mol=h.gc+h.pc;


X=htc.w_p;
y=htc.h_eff-htc.h_mol;

% fx=@(b,x) (b(1).*x(:,1)+b(2))./(x(:,1).^2+b(3).*x(:,1)+b(4));
% beta0=[-1.6902e+03,2.3766e+03,-4.3054e+03,-20.4956]'; %curve fitter results
% beta0=[-2e3,10e3,-4e3,-2e1]';
% [mdl,beta]=getFit(X,y,fx,beta0);
% htc.h_cf=htc.h_eff-htc.h_mol-fx(beta,0);


[fit,gof]=createFit(X,y);
htc.h_cf=htc.h_eff-htc.h_mol-fit(0);


fig=figure(100);
clf(fig);
ax=gca();
colors=ax.ColorOrder;
hold(ax,'on');

x2=linspace(0,10e-3,1000);
scatter(ax,X,htc.h_eff-htc.h_mol,18,colors(1,:));
% plot(ax,x2,fx(beta,x2'),'Color',colors(1,:));
plot(ax,x2,fit(x2),'Color',colors(1,:));

scatter(ax,X,htc.h_cf,18,colors(2,:));
% plot(ax,x2,fx(beta,x2')-fx(beta,0),'Color',colors(2,:));
plot(ax,x2,fit(x2)-fit(0),'Color',colors(2,:));

legItems=repmat(line(ax,'Visible','off'),2,1);
legItems(1)=plot(ax,NaN,NaN,'Color',colors(1,:),'Marker','o','MarkerSize',sqrt(18));
legItems(2)=plot(ax,NaN,NaN,'Color',colors(2,:),'Marker','o','MarkerSize',sqrt(18));

hold(ax,'off');

legend(ax,legItems,{'h_{eff}-h_{mol}','h_{cf}'},'Location','best');

xlabel(ax,'w_p (m/s)');
ylabel(ax,'HTC (W/m²K)');

fig.Units='centimeters';
fig.Position=[10,5,17,8.5];


%%
figure(101)

scatter(htc.h_cf,htc.h_cf./(htc.h_eff-fit(0)));


%% Dimensionless numbers
%Gas and particle properties
k_g=DryAir.lambda(htc.Tbed);
my_g=DryAir.eta(htc.Tbed);
c_p=SiO2.c_p(htc.Tbed);
rho_g=DryAir.rho(htc.p,htc.Tbed);

rho_e=c.rho_p-rho_g;
l_l=(my_g./(rho_e.*sqrt(FluBed.g))).^(2/3);
conv2cond=(rho_e.*c_p./(k_g.*FluBed.g)).^(1/3);


%Fluidization velocities
htc.w_mf=FluBed.wmfErgun(c.d_p,c.rho_p,c.phi_s,c.eps_mf,htc.p,htc.Tbed);
% w_mf=FluBed.wmf(c.d_p,c.rho_p,htc.p,htc.Tbed);
w_e=htc.w-htc.w_mf;
w_e(w_e<0)=NaN;


%Pi-factors
npi=10;
pis=table('Size',[height(htc),npi],...
    'VariableTypes',repmat({'double'},1,npi),...
    'VariableNames',compose('pi%d',1:npi));

pis.pi1=htc.h_cf.*l_l./k_g;
pis.pi2=k_g./(2*c_p.*my_g);
pis.pi3=DryAir.Pr(htc.Tbed);
pis.pi4=rho_g./rho_e;
pis.pi5=conv2cond.*w_e;
pis.pi6=conv2cond.*htc.w_mf;
pis.pi7=repmat(1-c.eps_mf,height(pis),1);
pis.pi8=c.d_t./l_l;
pis.pi9=repmat(c.d_t./c.p_hEff,height(pis),1);
pis.pi10=conv2cond.*htc.w_p;


%% Save for future analysis
save('intMeas','pis','htc');


%% Univariate plot
pisNorm=normalize(pis);

fig=figure(200);
clf(fig);
ax=gca();

boxplot(ax,pisNorm{:,:});

xlabel(ax,'\pi-index (-)');
ylabel(ax,'z-score (-)');


%% Bivariate plot
% covplot(pis{:,:},compose('\\pi_{%d}',1:size(pis,2)),300);


%% Model preparation
y=pis.pi1;
X=pis{:,:};


%% Model 1b1
% 
fx=@(b,x) ...
    b(1).*x(:,7)./...
    (1+b(2).*(x(:,6)./x(:,10)).^(1/3)./x(:,10).*...
    (1-exp(-6.45824254376922e-05.*x(:,8))).*...
    (1+0.28.*x(:,7).^2.*sqrt(x(:,4)).*x(:,5).*x(:,6))).*...
    ((1+(x(:,5)./x(:,10)).^b(3).*(1-x(:,9)).^b(4)).^-1);


beta0=[0.05,1,1,1]';

[mdl_1b1,beta_1b1]=getFit(X,y,fx,beta0);

checkFit(X,y,fx,beta_1b1,'1b1',1);


%% Model 1b2
% 
fx=@(b,x) ...
    b(1).*x(:,7)./...
    (1+b(2).*(x(:,6)./x(:,10)).^(1/3)./x(:,10).*...
    (1-exp(-6.45824254376922e-05.*x(:,8))).*...
    (1+0.28.*x(:,7).^2.*sqrt(x(:,4)).*x(:,5).*x(:,6))).*...
    ((1+(x(:,5)./x(:,10)).^b(3).*(1-x(:,9)).^b(3)).^-1);


beta0=[0.05,1,1]';

[mdl_1b2,beta_1b2]=getFit(X,y,fx,beta0);

checkFit(X,y,fx,beta_1b2,'1b2',2);


%% Model 1c1
% 
fx=@(b,x) ...
    b(1).*x(:,7)./...
    (1+b(2).*(x(:,6)./x(:,10)).^(1/3)./x(:,10).*...
    (1-exp(-6.45824254376922e-05.*x(:,8))).*...
    (1+0.28.*x(:,7).^2.*sqrt(x(:,4)).*x(:,5).*x(:,6))).*...
    ((1+(x(:,5)./x(:,6)).^b(3).*(x(:,5)./x(:,10)).^b(4).*(1-x(:,9)).^b(5)).^-1);


beta0=[0.05,1,1,1,1]';

[mdl_1c1,beta_1c1]=getFit(X,y,fx,beta0);

checkFit(X,y,fx,beta_1c1,'1c1',3);


%% Model 1c2
% 
fx=@(b,x) ...
    b(1).*x(:,7)./...
    (1+b(2).*(x(:,6)./x(:,10)).^(1/3)./x(:,10).*...
    (1-exp(-6.45824254376922e-05.*x(:,8))).*...
    (1+0.28.*x(:,7).^2.*sqrt(x(:,4)).*x(:,5).*x(:,6))).*...
    ((1+(x(:,5)./x(:,6)).^b(3).*(x(:,5)./x(:,10)).^b(4).*(1-x(:,9)).^(b(3)+b(4))).^-1);


beta0=[0.05,1,1,1]';

[mdl_1c2,beta_1c2]=getFit(X,y,fx,beta0);

checkFit(X,y,fx,beta_1c2,'1c2',4);


%% Model 2b1
% 
fx=@(b,x) ...
    b(1).*x(:,7)./...
    (1+b(2).*(x(:,6)./x(:,10)).^(1/3)./x(:,10).*...
    (1-exp(-6.45824254376922e-05.*x(:,8))).*...
    (1+0.28.*x(:,7).^2.*sqrt(x(:,4)).*x(:,5).*x(:,6))).*...
    exp(-(x(:,5)./x(:,10)).^b(3).*(1-x(:,9)).^b(4));


beta0=[0.05,1,1,1]';

[mdl_2b1,beta_2b1]=getFit(X,y,fx,beta0);

checkFit(X,y,fx,beta_2b1,'2b1',5);


%% Model 2b2
% 
fx=@(b,x) ...
    b(1).*x(:,7)./...
    (1+b(2).*(x(:,6)./x(:,10)).^(1/3)./x(:,10).*...
    (1-exp(-6.45824254376922e-05.*x(:,8))).*...
    (1+0.28.*x(:,7).^2.*sqrt(x(:,4)).*x(:,5).*x(:,6))).*...
    exp(-(x(:,5)./x(:,10)).^b(3).*(1-x(:,9)).^b(3));


beta0=[0.05,1,1]';

[mdl_2b2,beta_2b2]=getFit(X,y,fx,beta0);

checkFit(X,y,fx,beta_2b2,'2b2',6);


%% Model 2c1
% 
fx=@(b,x) ...
    b(1).*x(:,7)./...
    (1+b(2).*(x(:,6)./x(:,10)).^(1/3)./x(:,10).*...
    (1-exp(-6.45824254376922e-05.*x(:,8))).*...
    (1+0.28.*x(:,7).^2.*sqrt(x(:,4)).*x(:,5).*x(:,6))).*...
    exp(-(x(:,5)./x(:,6)).^b(3).*(x(:,5)./x(:,10)).^b(4).*(1-x(:,9)).^b(5));


beta0=[0.05,1,1,1,1]';

[mdl_2c1,beta_2c1]=getFit(X,y,fx,beta0);

checkFit(X,y,fx,beta_2c1,'2c1',7);


%% Model 2c2
% 
fx=@(b,x) ...
    b(1).*x(:,7)./...
    (1+b(2).*(x(:,6)./x(:,10)).^(1/3)./x(:,10).*...
    (1-exp(-6.45824254376922e-05.*x(:,8))).*...
    (1+0.28.*x(:,7).^2.*sqrt(x(:,4)).*x(:,5).*x(:,6))).*...
    exp(-(x(:,5)./x(:,6)).^b(3).*(x(:,5)./x(:,10)).^b(4).*(1-x(:,9)).^(b(3)+b(4)));


beta0=[0.05,1,1,1]';

[mdl_2c2,beta_2c2]=getFit(X,y,fx,beta0);

checkFit(X,y,fx,beta_2c2,'2c2',8);


%% Model 3b1
% 
fx=@(b,x) ...
    b(1).*x(:,7)./...
    (1+b(2).*(x(:,6)./x(:,10)).^(1/3)./x(:,10).*...
    (1-exp(-6.45824254376922e-05.*x(:,8))).*...
    (1+0.28.*x(:,7).^2.*sqrt(x(:,4)).*x(:,5).*x(:,6))).*...
    (1-tanh((x(:,5)./x(:,10)).^b(3).*(1-x(:,9)).^b(4)));


beta0=[0.05,1,1,1]';

[mdl_3b1,beta_3b1]=getFit(X,y,fx,beta0);

checkFit(X,y,fx,beta_3b1,'3b1',9);


%% Model 3b2
% 
fx=@(b,x) ...
    b(1).*x(:,7)./...
    (1+b(2).*(x(:,6)./x(:,10)).^(1/3)./x(:,10).*...
    (1-exp(-6.45824254376922e-05.*x(:,8))).*...
    (1+0.28.*x(:,7).^2.*sqrt(x(:,4)).*x(:,5).*x(:,6))).*...
    (1-tanh((x(:,5)./x(:,10)).^b(3).*(1-x(:,9)).^b(3)));


beta0=[0.05,1,1]';

[mdl_3b2,beta_3b2]=getFit(X,y,fx,beta0);

checkFit(X,y,fx,beta_3b2,'3b2',10);


%% Model 3c1
% 
fx=@(b,x) ...
    b(1).*x(:,7)./...
    (1+b(2).*(x(:,6)./x(:,10)).^(1/3)./x(:,10).*...
    (1-exp(-6.45824254376922e-05.*x(:,8))).*...
    (1+0.28.*x(:,7).^2.*sqrt(x(:,4)).*x(:,5).*x(:,6))).*...
    (1-tanh((x(:,5)./x(:,6)).^b(3).*(x(:,5)./x(:,10)).^b(4).*(1-x(:,9)).^b(5)));


beta0=[0.05,1,1,1,1]';

[mdl_3c1,beta_3c1]=getFit(X,y,fx,beta0);

checkFit(X,y,fx,beta_3c1,'3c1',11);


%% Model 3c2
% Of models 1-3, only c2 models have positive b(2). This one has slightly 
% worse R²_adj than the other c2 models, but its b(2) has the lowest 
% standard error / p-value
fx3c2=@(b,x) ...
    b(1).*x(:,7)./...
    (1+b(2).*(x(:,6)./x(:,10)).^(1/3)./x(:,10).*...
    (1-exp(-6.45824254376922e-05.*x(:,8))).*...
    (1+0.28.*x(:,7).^2.*sqrt(x(:,4)).*x(:,5).*x(:,6))).*...
    (1-tanh((x(:,5)./x(:,6)).^b(3).*(x(:,5)./x(:,10)).^b(4).*(1-x(:,9)).^(b(3)+b(4))));


beta0=[0.05,1,1,1]';

[mdl_3c2,beta_3c2]=getFit(X,y,fx3c2,beta0);

checkFit(X,y,fx3c2,beta_3c2,'3c2',12);


%% 3c2 compared to Molerus
yEst=fx3c2(beta_3c2,X);

R2=mdl_3c2.Rsquared.Adjusted;


cidx=1:4;
names={'Parameter','H0','Estimate','p'};
para=table('Size',[length(cidx)+1,length(names)],...
    'VariableTypes',[{'string'},repmat({'double'},1,length(names)-1)],...
    'VariableNames',names);


para.Parameter=[compose('C%d',cidx),{'C3+C4'}]';
para.Estimate(1:end-1)=beta_3c2;
para.Estimate(end)=sum(beta_3c2(3:4));


for i=1:height(para)-1
    beta_null=para.Estimate;
    beta_null(i)=para.H0(i);

    [~,p]=ttest(fx3c2(beta_null,X),y);
    para.p(i)=p;
end


beta_null=para.Estimate;
beta_null(3:4)=para.H0(3:4);
[~,para.p(end)]=ttest(fx3c2(beta_null,X),y);


%% Control strategy analysis
names={'Strategy','mode','Tsurf','Tbed','DeltaT','P','Nu_cf','ME'};
strat={'Pel=const.';'Tsurf-Tbed=const.'};
contr=table('Size',[length(strat),length(names)],...
    'VariableTypes',[{'string','logical'},repmat({'double'},1,length(names)-2)],...
    'VariableNames',names);


contr.Strategy=strat;
contr.mode(2)=true;

contr.Tsurf=arrayfun(@(tf) mean(htc.Tsurf(htc.mode==tf)),contr.mode);
contr.Tbed=arrayfun(@(tf) mean(htc.Tbed(htc.mode==tf)),contr.mode);
contr.P=arrayfun(@(tf) mean(htc.P(htc.mode==tf)),contr.mode);
contr.Nu_cf=arrayfun(@(tf) mean(y(htc.mode==tf)),contr.mode);

contr.ME=arrayfun(@(tf) mean(yEst(htc.mode==tf)-y(htc.mode==tf)),contr.mode);

contr.DeltaT=contr.Tsurf-contr.Tbed;


Nu_cfRel=contr.Nu_cf(2)./contr.Nu_cf(1);


%%
figidx=4;
fig=figure(figidx);
clf(fig);
t=tiledlayout(fig,1,1,'Padding','tight');
ax=nexttile(t);
colors=ax.ColorOrder;
hold(ax,'on');


legItems=cell(1,2);

legItems{1}=scatter(ax,y(~htc.mode),yEst(~htc.mode),18,colors(1,:),'o');
legItems{2}=scatter(ax,y(htc.mode),yEst(htc.mode),18,colors(2,:),'+');

lim=max([ax.XLim(2),ax.YLim(2)]);
eq=linspace(0,lim,100);

plot(ax,eq,eq,'Color','k');
plot(ax,eq,eq.*1.2,'Color','k','LineStyle','--');
plot(ax,eq,eq./1.2,'Color','k','LineStyle','--');

hold(ax,'off');


legItems=[legItems{:}];
txt=subsz({'P_{el} = const.','T_{surf} - T_{bed} = const.'},6);
lgd=legend(ax,legItems,txt,'Location','southeast');


ax.XLim=[0,lim];
ax.YLim=[0,lim];

xlabel(ax,subsz('Measured Nu_{cf} (-)',6));
ylabel(ax,subsz('Estimated Nu_{cf} (-)',6));


%Text size
ax.FontSize=7;
lgd.FontSize=7;


%Export figure for manuscript
fname=['Figures',filesep,'Figure',num2str(figidx)];

t.Units='centimeters';
t.InnerPosition=[1.5,1,8.2,8.2];

fig.Units=t.Units;
fig.Position(3:4)=t.OuterPosition(3:4)+1;

exportgraphics(fig,[fname,'.tiff'],'Resolution',600);


%Export figure for Elsevier
t.InnerPosition=[1.5,1,8.2,8.2];
fig.Position(3:4)=t.OuterPosition(3:4)+1;

exportgraphics(fig,[fname,'.eps']);


% %% Model 4b1
% % 
% fx=@(b,x) ...
%     b(1).*x(:,7)./...
%     (1+b(2).*(x(:,6)./x(:,10)).^(1/3)./x(:,10).*...
%     (1-exp(-6.45824254376922e-05.*x(:,8))).*...
%     (1+0.28.*x(:,7).^2.*sqrt(x(:,4)).*x(:,5).*x(:,6)).*...
%     (1+(x(:,5)./x(:,10)).^b(3).*(1-x(:,9)).^b(4)));
% 
% 
% beta0=[0.05,1,1,1]';
% 
% [mdl_4b1,beta_4b1]=getFit(X,y,fx,beta0);
% 
% checkFit(X,y,fx,beta_4b1,'4b1',13);
% 
% 
% %% Model 4b2
% % 
% fx=@(b,x) ...
%     b(1).*x(:,7)./...
%     (1+b(2).*(x(:,6)./x(:,10)).^(1/3)./x(:,10).*...
%     (1-exp(-6.45824254376922e-05.*x(:,8))).*...
%     (1+0.28.*x(:,7).^2.*sqrt(x(:,4)).*x(:,5).*x(:,6)).*...
%     (1+(x(:,5)./x(:,10)).^b(3).*(1-x(:,9)).^b(3)));
% 
% 
% beta0=[0.05,1,1]';
% 
% [mdl_4b2,beta_4b2]=getFit(X,y,fx,beta0);
% 
% checkFit(X,y,fx,beta_4b2,'4b2',14);
% 
% 
% %% Model 4c1
% % 
% fx=@(b,x) ...
%     b(1).*x(:,7)./...
%     (1+b(2).*(x(:,6)./x(:,10)).^(1/3)./x(:,10).*...
%     (1-exp(-6.45824254376922e-05.*x(:,8))).*...
%     (1+0.28.*x(:,7).^2.*sqrt(x(:,4)).*x(:,5).*x(:,6)).*...
%     (1+(x(:,5)./x(:,6)).^b(3).*(x(:,5)./x(:,10)).^b(4).*(1-x(:,9)).^b(5)));
% 
% 
% beta0=[0.05,1,1,1,1]';
% 
% [mdl_4c1,beta_4c1]=getFit(X,y,fx,beta0);
% 
% checkFit(X,y,fx,beta_4c1,'4c1',15);
% 
% 
% %% Model 4c2
% % All Model 4s: poor R²_adj compared to Model 3s
% fx=@(b,x) ...
%     b(1).*x(:,7)./...
%     (1+b(2).*(x(:,6)./x(:,10)).^(1/3)./x(:,10).*...
%     (1-exp(-6.45824254376922e-05.*x(:,8))).*...
%     (1+0.28.*x(:,7).^2.*sqrt(x(:,4)).*x(:,5).*x(:,6)).*...
%     (1+(x(:,5)./x(:,6)).^b(3).*(x(:,5)./x(:,10)).^b(4).*(1-x(:,9)).^(b(3)+b(4))));
% 
% 
% beta0=[0.05,1,1,1]';
% 
% [mdl_4c2,beta_4c2]=getFit(X,y,fx,beta0);
% 
% checkFit(X,y,fx,beta_4c2,'4c2',16);


% %% Model 5c2
% % 
% fx=@(b,x) ...
%     b(1).*x(:,7)./...
%     (1+b(2).*(x(:,6)./x(:,10)).^(1/3)./x(:,10).*...
%     (1-exp(-6.45824254376922e-05.*x(:,8))).*...
%     (1+0.28.*x(:,7).^2.*sqrt(x(:,4)).*x(:,5).*x(:,6)).*...
%     (x(:,5)./x(:,6)).^b(3).*(1-x(:,9)).^b(3));
% 
% 
% beta0=[0.05,1,1]';
% 
% [mdl_5c2,beta_5c2]=getFit(X,y,fx,beta0);
% 
% checkFit(X,y,fx,beta_5c2,'5c2',17);


% %% Model 6c2
% % p-value of b(2) and b(3) is fairly high compared to Model 5
% fx=@(b,x) ...
%     b(1).*x(:,7)./...
%     (1+b(2).*(x(:,6)./x(:,10)).^(1/3)./x(:,10).*...
%     (1-exp(-6.45824254376922e-05.*x(:,8))).*...
%     (1+0.28.*x(:,7).^2.*sqrt(x(:,4)).*x(:,5).*x(:,6)).*...
%     (x(:,5)./x(:,6)).^b(3).*(x(:,5)./x(:,10)).^b(4).*(1-x(:,9)).^(b(3)+b(4)));
% 
% 
% beta0=[0.05,1,1,1]';
% 
% [mdl_6c2,beta_6c2]=getFit(X,y,fx,beta0);
% 
% checkFit(X,y,fx,beta_6c2,'6c2',18);


% %% Model 7c2
% % 
% fx=@(b,x) ...
%     b(1).*x(:,7)./...
%     (1+b(2).*(x(:,6)./x(:,10)).^(1/3)./x(:,10).*...
%     tanh(b(4).*x(:,8)).*...
%     (1+b(5).*x(:,7).^2.*sqrt(x(:,4)).*x(:,5).*x(:,6)).*...
%     (x(:,5)./x(:,6)).^b(3).*(1-x(:,9)).^b(3));
% 
% 
% beta0=[0.05,1,1,4.5e-7,0.33]';
% 
% [mdl_7c2,beta_7c2]=getFit(X,y,fx,beta0);
% 
% checkFit(X,y,fx,beta_7c2,'7c2',19);



