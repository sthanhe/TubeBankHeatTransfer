dirFigures=['Figures',filesep,'PCmodels'];   %Path to directory where figures should be stored
fname='secData.mat';


%% Make figure and table folders if they do not exist
if ~isfolder(dirFigures)
    mkdir(dirFigures);
end


%% Load data
load(fname);


%Create response variable y and regressor matrix X: only data from Grewal
idx=strcmp(sec.Author,'Grewal');
y=pis.pi1(idx)-sec.Nu_gcMW(idx);
X=pis{idx,:};


%% Null model
fx0=@(b,x) ...
    b(1).*x(:,7)./...
    (1+x(:,2).*...
    (1-exp(-b(4).*x(:,8))).*...
    (1+b(2).*x(:,7).^2.*sqrt(x(:,4)).*x(:,5).*x(:,6))).*...
    (1+b(3).*(x(:,6)./x(:,5)).^(1/3)./x(:,5)./(1-x(:,9)).^b(5)).^-1;


beta_s0=[0.125,0.28,33.3,Inf,0]';

checkFit(X,y,fx0,beta_s0,'s0',601,false,dirFigures);
checkFit(X,y,fx0,beta_s0,'s0',600,true,dirFigures);

compSec(fx0,beta_s0,dirFigures,0,true);
compPi8(fx0,beta_s0,dirFigures,0,true);
compPi9(fx0,beta_s0,dirFigures,0,true);


%% Model s1
fx1=@(b,x) ...
    b(1).*x(:,7)./...
    (1+x(:,2).*...
    (1-(1+b(3).*x(:,8).^b(5)).^-1).*...
    (1+0.28.*x(:,7).^2.*sqrt(x(:,4)).*x(:,5).*x(:,6))).*...
    (1+b(2).*(x(:,6)./x(:,5)).^(1/3)./x(:,5)./(1-x(:,9)).^b(4)).^-1;


beta0=[0.125,33.3,5e-8,1,2]';

[mdl_s1,beta_s1]=getFit(X,y,fx1,beta0);

checkFit(X,y,fx1,beta_s1,'s1',611,false,dirFigures);
checkFit(X,y,fx1,beta_s1,'s1',610,true,dirFigures);

compSec(fx1,beta_s1,dirFigures,1,true);
compPi8(fx1,beta_s1,dirFigures,1,true);
compPi9(fx1,beta_s1,dirFigures,1,true);


%% Model s2
fx2=@(b,x) ...
    b(1).*x(:,7)./...
    (1+x(:,2).*...
    (1-exp(-b(3).*x(:,8))).*...
    (1+0.28.*x(:,7).^2.*sqrt(x(:,4)).*x(:,5).*x(:,6))).*...
    (1+b(2).*(x(:,6)./x(:,5)).^(1/3)./x(:,5)./(1-x(:,9)).^b(4)).^-1;


beta0=[0.125,33.3,5e-5,1]';

[mdl_s2,beta_s2]=getFit(X,y,fx2,beta0);

checkFit(X,y,fx2,beta_s2,'s2',621,false,dirFigures);
checkFit(X,y,fx2,beta_s2,'s2',620,true,dirFigures);

compSec(fx2,beta_s2,dirFigures,2,true);
compPi8(fx2,beta_s2,dirFigures,2,true);
compPi9(fx2,beta_s2,dirFigures,2,true);


%% Model s3
fx3=@(b,x) ...
    b(1).*x(:,7)./...
    (1+x(:,2).*...
    tanh(b(3).*x(:,8)).*...
    (1+0.28.*x(:,7).^2.*sqrt(x(:,4)).*x(:,5).*x(:,6))).*...
    (1+b(2).*(x(:,6)./x(:,5)).^(1/3)./x(:,5)./(1-x(:,9)).^b(4)).^-1;


beta0=[0.125,33.3,5e-5,1]';

[mdl_s3,beta_s3]=getFit(X,y,fx3,beta0);

checkFit(X,y,fx3,beta_s3,'s3',631,false,dirFigures);
checkFit(X,y,fx3,beta_s3,'s3',630,true,dirFigures);

compSec(fx3,beta_s3,dirFigures,3,true);
compPi8(fx3,beta_s3,dirFigures,3,true);
compPi9(fx3,beta_s3,dirFigures,3,true);


%% Evaluation
%Hypothesis tests (p-values of parameters)
para1=hyptest(X,y,fx1,beta_s1,[beta_s0;0],[1,3:6]);
para2=hyptest(X,y,fx2,beta_s2,beta_s0,[1,3:5]);
para3=hyptest(X,y,fx3,beta_s3,beta_s0,[1,3:5]);


%Goodness-of-fit parameters
mdls={mdl_s1;mdl_s2;mdl_s3};
names={'Model','R2_adj','AIC'};
gof=table('Size',[length(mdls)+1,length(names)],...
    'VariableTypes',[{'string'},repmat({'double'},1,length(names)-1)],...
    'VariableNames',names);

gof.Model=compose('s%d',0:length(mdls))';

gof.R2_adj(1)=Rsq(y,fx0(beta_s0,X),9);
gof.R2_adj(2:end)=cellfun(@(x) x.Rsquared.Adjusted,mdls);

gof.AIC(1)=NaN;
gof.AIC(2:end)=cellfun(@(x) x.ModelCriterion.AIC,mdls);


%% Plot graphics
%Choose best model
fx=fx2;
beta_s=beta_s2;
mdl=mdl_s2;


%Comparison to secondary data (for manuscript)
compSec(fx,beta_s,dirFigures,7,false);


%Set up figure
figidx=5;
fig=figure(figidx);
clf(fig);
t=tiledlayout(fig,1,1,'Padding','tight');
ax=nexttile(t);
colors=ax.ColorOrder;
hold(ax,'on');


legItems=cell(1,2);

legItems{1}=scatter(ax,y,fx0(beta_s0,X),18,'o');
legItems{2}=scatter(ax,y,fx(beta_s,X),18,'+');

lim=max([ax.XLim(2),ax.YLim(2)]);
eq=linspace(0,lim,100);

plot(ax,eq,eq,'Color','k');
plot(ax,eq,eq.*1.2,'Color','k','LineStyle','--');
plot(ax,eq,eq./1.2,'Color','k','LineStyle','--');

hold(ax,'off');


legItems=[legItems{:}];
txt=figaux.subsz(compose('H_%d, R^2_{adj}=%.3f',...
    [0;1],[gof.R2_adj(1);mdl.Rsquared.Adjusted]),6);
lgd=legend(ax,legItems,txt,'Location','southeast');


ax.XLim=[0,lim];
ax.YLim=[0,lim];

xlabel(ax,figaux.subsz('Measured Nu_{pc} (-)',6));
ylabel(ax,figaux.subsz('Estimated Nu_{pc} (-)',6));


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




