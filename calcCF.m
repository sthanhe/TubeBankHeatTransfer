%% Set data directories
dirFigures=['Figures',filesep,'CFmodels'];   %Path to directory where figures should be stored
fname='primData.mat';


if ~isfolder(dirFigures)
    mkdir(dirFigures);
end


%% Load data
load(fname);



% %%
% figure(101)
% 
% scatter(prim.h_cf,prim.h_cf./(prim.h_eff-fit(0)));





%% Model preparation
idx=true(1,height(pis));
% idx=pis.pi1>0;
% idx=prim.mDot_p>1;
% idx=pis.pi1>0 & prim.mDot_p>0.5;

y=pis.pi1(idx);
X=pis{idx,:};


base=['@(b,x) ',...
    'b(1).*x(:,7)./',...
    '(1+b(2).*(x(:,6)./x(:,10)).^(1/3)./x(:,10).*',...
    '(1+0.28.*x(:,7).^2.*sqrt(x(:,4)).*x(:,5).*x(:,6)).*'];

s='(1-exp(-6.45824254376922e-05.*x(:,8)))).*';


fx_1start='((1+';
fx_1end=').^-1)';

fx_2start='exp(-';
fx_2end=')';

fx_3start='(1-tanh(';
fx_3end='))';


opt_a1='(x(:,5)./x(:,10)).^b(3).*(1-x(:,9)).^b(3)';
opt_a2='(x(:,5)./x(:,10)).^b(3).*(1-x(:,9)).^b(4)';

opt_b2='(x(:,5)./x(:,10)).^b(3).*(x(:,5)./x(:,6)).^b(4).*(1-x(:,9)).^(b(3)+b(4))';
opt_b3='(x(:,5)./x(:,10)).^b(3).*(x(:,5)./x(:,6)).^b(4).*(1-x(:,9)).^b(5)';


beta0_1=[0.05,1,1]';
beta0_2=[0.05,1,1,1]';
beta0_3=[0.05,1,1,1,1]';


base=[base,s];


%% Model 1a1
fx_1a1=eval([base,fx_1start,opt_a1,fx_1end]);


[mdl_1a1,beta_1a1]=getFit(X,y,fx_1a1,beta0_1);

checkFit(X,y,fx_1a1,beta_1a1,'1a1',111,false,dirFigures);


%% Model 1a2
fx_1a2=eval([base,fx_1start,opt_a2,fx_1end]);


[mdl_1a2,beta_1a2]=getFit(X,y,fx_1a2,beta0_2);

checkFit(X,y,fx_1a2,beta_1a2,'1a2',112,false,dirFigures);


%% Model 1b2
fx_1b2=eval([base,fx_1start,opt_b2,fx_1end]);


[mdl_1b2,beta_1b2]=getFit(X,y,fx_1b2,beta0_2);

checkFit(X,y,fx_1b2,beta_1b2,'1b2',122,false,dirFigures);


%% Model 1b3
fx_1b3=eval([base,fx_1start,opt_b3,fx_1end]);


[mdl_1b3,beta_1b3]=getFit(X,y,fx_1b3,beta0_3);

checkFit(X,y,fx_1b3,beta_1b3,'1b3',123,false,dirFigures);


%% Model 2a1
fx_2a1=eval([base,fx_2start,opt_a1,fx_2end]);


[mdl_2a1,beta_2a1]=getFit(X,y,fx_2a1,beta0_1);

checkFit(X,y,fx_2a1,beta_2a1,'2a1',211,false,dirFigures);


%% Model 2a2
fx_2a2=eval([base,fx_2start,opt_a2,fx_2end]);


[mdl_2a2,beta_2a2]=getFit(X,y,fx_2a2,beta0_2);

checkFit(X,y,fx_2a2,beta_2a2,'2a2',212,false,dirFigures);


%% Model 2b2
fx_2b2=eval([base,fx_2start,opt_b2,fx_2end]);


[mdl_2b2,beta_2b2]=getFit(X,y,fx_2b2,beta0_2);

checkFit(X,y,fx_2b2,beta_2b2,'2b2',222,false,dirFigures);


%% Model 2b3
fx_2b3=eval([base,fx_2start,opt_b3,fx_2end]);


[mdl_2b3,beta_2b3]=getFit(X,y,fx_2b3,beta0_3);

checkFit(X,y,fx_2b3,beta_2b3,'2b3',223,false,dirFigures);


%% Model 3a1
fx_3a1=eval([base,fx_3start,opt_a1,fx_3end]);


[mdl_3a1,beta_3a1]=getFit(X,y,fx_3a1,beta0_1);

checkFit(X,y,fx_3a1,beta_3a1,'3a1',311,false,dirFigures);


%% Model 3a2
fx_3a2=eval([base,fx_3start,opt_a2,fx_3end]);


[mdl_3a2,beta_3a2]=getFit(X,y,fx_3a2,beta0_2);

checkFit(X,y,fx_3a2,beta_3a2,'3a2',312,false,dirFigures);


%% Model 3b2
fx_3b2=eval([base,fx_3start,opt_b2,fx_3end]);


[mdl_3b2,beta_3b2]=getFit(X,y,fx_3b2,beta0_2);

checkFit(X,y,fx_3b2,beta_3b2,'3b2',322,false,dirFigures);


%% Model 3b3
fx_3b3=eval([base,fx_3start,opt_b3,fx_3end]);


[mdl_3b3,beta_3b3]=getFit(X,y,fx_3b3,beta0_3);

checkFit(X,y,fx_3b3,beta_3b3,'3b3',323,false,dirFigures);


%% Comparison to Molerus
% Choose best model
fx=fx2;
beta_s=beta_s2;
mdl=mdl_s2;


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


%%
figidx=4;
fig=figure(figidx);
clf(fig);
t=tiledlayout(fig,1,1,'Padding','tight');
ax=nexttile(t);
colors=ax.ColorOrder;
hold(ax,'on');


legItems=cell(1,2);

legItems{1}=scatter(ax,y(~prim.mode),yEst(~prim.mode),18,colors(1,:),'o');
legItems{2}=scatter(ax,y(prim.mode),yEst(prim.mode),18,colors(2,:),'+');

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


%% Control strategy analysis
names={'Strategy','mode','Tsurf','Tbed','DeltaT','P_el',...
    'Nu_cf','ME'};
strat={'Pel=const.';'Tsurf-Tbed=const.'};
contr=table('Size',[length(strat),length(names)],...
    'VariableTypes',[{'string','logical'},repmat({'double'},1,length(names)-2)],...
    'VariableNames',names);


contr.Strategy=strat;
contr.mode(2)=true;

contr.Tsurf=arrayfun(@(tf) mean(prim.Tsurf(prim.mode==tf)),contr.mode);
contr.Tbed=arrayfun(@(tf) mean(prim.Tbed(prim.mode==tf)),contr.mode);
contr.P_el=arrayfun(@(tf) mean(prim.P_el(prim.mode==tf)),contr.mode);
contr.Nu_cf=arrayfun(@(tf) mean(y(prim.mode==tf)),contr.mode);

contr.ME=arrayfun(@(tf) mean(yEst(prim.mode==tf)-y(prim.mode==tf)),contr.mode);

contr.DeltaT=contr.Tsurf-contr.Tbed;


Nu_cfRel=contr.Nu_cf(2)./contr.Nu_cf(1);







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



