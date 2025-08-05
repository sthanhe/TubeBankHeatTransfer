%% Analyze particle-convection
% GNU General Public License v3.0
% By Stefan Thanheiser: https://orcid.org/0000-0003-2765-1156
%
% Part of the paper:
%
% Thanheiser, S.
% Molerus and Wirth's Heat Transfer Model for Bubbling Fluidized Beds: 
% Proposal for an Extended Model Including Immersed Tube Banks and Particle 
% Cross-Flow
%
% All data, along with methodology reports and supplementary documentation, 
% is published in the data repository:
% https://doi.org/10.5281/zenodo.15576311
%
% All required files for this script can be found in the software
% repository: 
% https://doi.org/10.5281/zenodo.15576950
%
%
%
% This script conducts the main analysis of the collected secondary data to
% find a suitable functional form of the size function s(pi8). It creates 
% all published figures and calculates the statistics mentioned in both the
% main paper and the Methodology Report. 
%
%
%Requires all auxiliary classes and functions on the MATLAB path
%
%Required products, version 24.1:
%   - MATLAB
%   - Statistics and Machine Learning Toolbox
%   - Curve Fitting Toolbox
%Necessary classes, functions, files, and scripts:
%   - @DryAir
%   - @FluBed
%   - @SiO2
%   - @figaux
%   - @implExp
%   - Rsq.m
%   - checkFit.m
%   - compPi8.m
%   - compPi9.m
%   - compSec.m
%   - getFit.m
%   - hypTest.m
%   - secData.mat --> created by the script "prepSec"


%% Set data locations
dirData=['Data',filesep,'Own'];             %Data storage folder
dirFigs=['Figures',filesep,'PCmodels'];     %Figure storage folder

fname='secData.mat';    %Secondary data file


%% Make folders if they do not exist
if ~isfolder(dirFigs)
    mkdir(dirFigs);
end


%% Load data
load(fname);


%Regressor matrix X and response variable y: only data from Grewal
idx=strcmp(sec.Author,'Grewal');
X=pis{idx,:};
y=pis.pi1(idx)-sec.Nu_gcMW(idx);


%% Null model
fx0=@(b,x) ...
    b(1).*x(:,7)./...
    (1+x(:,2).*...
    (1-exp(-b(4).*x(:,8))).*...
    (1+b(2).*x(:,7).^2.*sqrt(x(:,4)).*x(:,5).*x(:,6))).*...
    (1+b(3).*(x(:,6)./x(:,5)).^(1/3)./x(:,5)./(1-x(:,9)).^b(5)).^-1;


beta_s0=[0.125,0.28,33.3,Inf,0]';

checkFit(X,y,fx0,beta_s0,'s0',501,false,dirFigs);
checkFit(X,y,fx0,beta_s0,'s0',500,true,dirFigs);

compSec(fx0,beta_s0,pis,sec,dirFigs,0,true);
compPi8(fx0,beta_s0,dirFigs,0,true);
compPi9(fx0,beta_s0,dirFigs,0,true);


%% Model s1
fx1=@(b,x) ...
    b(1).*x(:,7)./...
    (1+x(:,2).*...
    (1-(1+b(3).*x(:,8).^b(5)).^-1).*...
    (1+0.28.*x(:,7).^2.*sqrt(x(:,4)).*x(:,5).*x(:,6))).*...
    (1+b(2).*(x(:,6)./x(:,5)).^(1/3)./x(:,5)./(1-x(:,9)).^b(4)).^-1;


beta0=[0.125,33.3,5e-8,1,2]';

[mdl_s1,beta_s1]=getFit(X,y,fx1,beta0);

checkFit(X,y,fx1,beta_s1,'s1',511,false,dirFigs);
checkFit(X,y,fx1,beta_s1,'s1',510,true,dirFigs);

compSec(fx1,beta_s1,pis,sec,dirFigs,1,true);
compPi8(fx1,beta_s1,dirFigs,1,true);
compPi9(fx1,beta_s1,dirFigs,1,true);


%% Model s2
fx2=@(b,x) ...
    b(1).*x(:,7)./...
    (1+x(:,2).*...
    (1-exp(-b(3).*x(:,8))).*...
    (1+0.28.*x(:,7).^2.*sqrt(x(:,4)).*x(:,5).*x(:,6))).*...
    (1+b(2).*(x(:,6)./x(:,5)).^(1/3)./x(:,5)./(1-x(:,9)).^b(4)).^-1;


beta0=[0.125,33.3,5e-5,1]';

[mdl_s2,beta_s2]=getFit(X,y,fx2,beta0);

checkFit(X,y,fx2,beta_s2,'s2',521,false,dirFigs);
checkFit(X,y,fx2,beta_s2,'s2',520,true,dirFigs);

compSec(fx2,beta_s2,pis,sec,dirFigs,2,true);
compPi8(fx2,beta_s2,dirFigs,2,true);
compPi9(fx2,beta_s2,dirFigs,2,true);


%% Model s3
fx3=@(b,x) ...
    b(1).*x(:,7)./...
    (1+x(:,2).*...
    tanh(b(3).*x(:,8)).*...
    (1+0.28.*x(:,7).^2.*sqrt(x(:,4)).*x(:,5).*x(:,6))).*...
    (1+b(2).*(x(:,6)./x(:,5)).^(1/3)./x(:,5)./(1-x(:,9)).^b(4)).^-1;


beta0=[0.125,33.3,5e-5,1]';

[mdl_s3,beta_s3]=getFit(X,y,fx3,beta0);

checkFit(X,y,fx3,beta_s3,'s3',531,false,dirFigs);
checkFit(X,y,fx3,beta_s3,'s3',530,true,dirFigs);

compSec(fx3,beta_s3,pis,sec,dirFigs,3,true);
compPi8(fx3,beta_s3,dirFigs,3,true);
compPi9(fx3,beta_s3,dirFigs,3,true);


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


%Comparison to secondary data and other models (for manuscript)
R2=compSec(fx,beta_s,pis,sec,dirFigs,8,false);
compPi8(fx,beta_s,dirFigs,9,false);
compPi9(fx,beta_s,dirFigs,10,false);


%Set up figure
figidx=5;
fig=figure(figidx);
clf(fig);
t=tiledlayout(fig,1,1,'Padding','tight');
ax=nexttile(t);
colors=ax.ColorOrder;
hold(ax,'on');


%Plot data
legItems=cell(1,2);     %Legend item container

legItems{1}=scatter(ax,y,fx0(beta_s0,X),18,'o');
legItems{2}=scatter(ax,y,fx(beta_s,X),18,'+');


%Plot equivalence lines
lim=max([ax.XLim(2),ax.YLim(2)]);
eq=linspace(0,lim,100);

plot(ax,eq,eq,'Color','k');
plot(ax,eq,eq.*1.2,'Color','k','LineStyle','--');
plot(ax,eq,eq./1.2,'Color','k','LineStyle','--');

hold(ax,'off');


%Set legend
txt=figaux.subsz(compose('H_%d, R^2_{adj}=%.3f',...
    [0;1],...
    [gof.R2_adj(1);mdl.Rsquared.Adjusted]),...
    6);
lgd=legend(ax,[legItems{:}],txt,'Location','southeast');


%Axes limits and labels
ax.XLim=[0,lim];
ax.YLim=[0,lim];

xlabel(ax,figaux.subsz('Measured Nu_{pc} (-)',6));
ylabel(ax,figaux.subsz('Estimated Nu_{pc} (-)',6));


%Text size
ax.FontSize=7;
lgd.FontSize=7;


%Export figure for manuscript
fname=[dirFigs,filesep,'Figure',num2str(figidx)];

t.Units='centimeters';
t.InnerPosition=[1.5,1,7.4,7.4];

fig.Units=t.Units;
fig.Position(3:4)=t.OuterPosition(3:4)+0.5;

exportgraphics(fig,[fname,'.tiff'],'Resolution',600);


%Export figure for Elsevier
t.InnerPosition=[1.5,1,7.4,7.4];
fig.Position(3:4)=t.OuterPosition(3:4)+0.5;

exportgraphics(fig,[fname,'.eps']);
savefig(fig,fname);




