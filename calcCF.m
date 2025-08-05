%% Analyze particle cross-flow
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
% This script conducts the main analysis of the collected primary data to
% find a suitable functional form of the cross-flow damping function 
% d_cf(pi5). It creates all published figures and calculates the statistics
% mentioned in both the main paper and the Methodology Report. 
%
%
%Requires all auxiliary classes and functions on the MATLAB path
%
%Required products, version 24.1:
%   - MATLAB
%   - Statistics and Machine Learning Toolbox
%Necessary classes, functions, files, and scripts:
%   - @figaux
%   - checkFit.m
%   - getFit.m
%   - primData.mat --> created by the script "prepPrim"


%% Set data locations
dirData=['Data',filesep,'Own'];             %Data storage folder
dirFigs=['Figures',filesep,'CFmodels'];     %Figure storage folder
dirTabs='Tables';                           %Table storage folder

fname='primData.mat';   %Primary data file


%% Make folders if they do not exist
if ~isfolder(dirFigs)
    mkdir(dirFigs);
end

if ~isfolder(dirTabs)
    mkdir(dirTabs);
end


%% Load data
load(fname);


%Regressor matrix X and response variable y
X=pis{:,:};
y=pis.pi1;


%% Function strings
%Basic
base='@(b,x) b().*x(:,7)./(1+';     %Start of function
b='b().*';                          %Coefficient


%Particle transport resistance r and turbulence function t
r='(x(:,6)./x(:,10)).^(1/3)./x(:,10).*';
t='(1+0.28.*x(:,7).^2.*sqrt(x(:,4)).*x(:,5).*x(:,6)).*';


%Size functions
s1='(1-(1+6.28756168714156e-11.*x(:,8).^2.58077440742616).^-1).*';
s2='(1-exp(-6.45824254376922e-05.*x(:,8))).*';
s3='tanh(5.37302742339688e-05.*x(:,8)).*';


%Particle-convection coefficients P5 for different size functions
P5s1=0.763002330862026;
P5s2=1.15226744565364;
P5s3=1.01505510898379;


%Choose size function
s=s2;
P5=P5s2;


%% Model 1
% All pi-factors, except pi3 (only gas-convection), pi7 (constant and 
% already included in nominator), pi8 (already included in size function),
% and (1-pi9) instead of pi9: insensitive parameters, R²=0.82, 
% AIC=-5.5, all p-values except those of pi5 and pi10 very high

fxstr=[base,r,t,s,'(1+',b,monostr([2,4,5,6,10]),...
    '.*(1-x(:,9)).^b()));'];
[fxstr,n]=bidx(fxstr);
fx1=eval(fxstr);

[mdl1,beta1]=getFit(X,y,fx1,[0.1,1e3,ones(1,n-2)]');

checkFit(X,y,fx1,beta1,'1',611,false,dirFigs);
checkFit(X,y,fx1,beta1,'1',610,true,dirFigs);


%% Model 2
% Like model 1, but removed b(2) and pi4: R²=0.82, AIC=-7.8, all p-values 
% except those of pi5 and pi10 very high

fxstr=[base,r,t,s,'(1+',monostr([2,5,6,10]),...
    '.*(1-x(:,9)).^b()));'];
[fxstr,n]=bidx(fxstr);
fx2=eval(fxstr);

[mdl2,beta2]=getFit(X,y,fx2,[0.1,ones(1,n-1)]');

checkFit(X,y,fx2,beta2,'2',621,false,dirFigs);
checkFit(X,y,fx2,beta2,'2',620,true,dirFigs);



%% Model 3
% Like model 2, but coupled exponent of pi9 to pi5 and scaled according to 
% particle convection results: R²=0.81, AIC=-7.1, all p-values except those
% of pi5 and pi10 very high

fxstr=[base,r,t,s,'(1+',monostr([2,5,6,10]),...
    '.*(1-x(:,9)).^(b()*',num2str(P5),'*3/4)));'];
fxstr=bidx(fxstr);
[fxstr,n]=couple(fxstr,[5,9]);
fx3=eval(fxstr);

[mdl3,beta3]=getFit(X,y,fx3,[0.1,ones(1,n-1)]');

checkFit(X,y,fx3,beta3,'3',631,false,dirFigs);
checkFit(X,y,fx3,beta3,'3',630,true,dirFigs);


%% Model 4
% Like model 2, but removed pi2: R²=0.81, AIC=-7.1, all p-values except 
% those of pi5 and pi10 very high

fxstr=[base,r,t,s,'(1+',monostr([5,6,10]),...
    '.*(1-x(:,9)).^b()));'];
[fxstr,n]=bidx(fxstr);
fx4=eval(fxstr);

[mdl4,beta4]=getFit(X,y,fx4,[0.1,ones(1,n-1)]');

checkFit(X,y,fx4,beta4,'4',641,false,dirFigs);
checkFit(X,y,fx4,beta4,'4',640,true,dirFigs);


%% Model 5
% Like model 2, but removed pi2 and coupled exponents of pi5 and pi9 like 
% in model 3: R²=0.81, AIC=-9.1, p-values of b(1) and pi6 about 2%

fxstr=[base,r,t,s,'(1+',monostr([5,6,10]),...
    '.*(1-x(:,9)).^(b()*',num2str(P5),'*3/4)));'];
fxstr=bidx(fxstr);
[fxstr,n]=couple(fxstr,[5,9]);
fx5=eval(fxstr);

[mdl5,beta5]=getFit(X,y,fx5,[0.1,ones(1,n-1)]');

checkFit(X,y,fx5,beta5,'5',651,false,dirFigs);
checkFit(X,y,fx5,beta5,'5',650,true,dirFigs);


%% Model 6
% Like model 5, but removed pi6: R²=0.83, AIC=-11.9, good p-values

fxstr=[base,r,t,s,'(1+',monostr([5,10]),...
    '.*(1-x(:,9)).^(b()*',num2str(P5),'*3/4)));'];
fxstr=bidx(fxstr);
[fxstr,n]=couple(fxstr,[5,9]);
fx6=eval(fxstr);

[mdl6,beta6]=getFit(X,y,fx6,[0.1,ones(1,n-1)]');

checkFit(X,y,fx6,beta6,'6',661,false,dirFigs,5e-3);
checkFit(X,y,fx6,beta6,'6',660,true,dirFigs,5e-3);


%% Comparison to Molerus
% Choose best model
fx=fx6;
beta=beta6;


%Estimated cross-flow Nusselt numbers
yEst=fx(beta,X);


%Set up figure
figidx=6;
fig=figure(figidx);
clf(fig);
til=tiledlayout(fig,1,1,'Padding','tight');
ax=nexttile(til);
colors=ax.ColorOrder;
hold(ax,'on');


%Plot data
legItems=cell(1,2);     %Legend item container

legItems{1}=scatter(ax,y(~prim.mode),yEst(~prim.mode),18,colors(1,:),'o');
legItems{2}=scatter(ax,y(prim.mode),yEst(prim.mode),18,colors(2,:),'+');


%Plot equivalence lines
lim=max([ax.XLim(2),ax.YLim(2)]);
eq=linspace(0,lim,100);

plot(ax,eq,eq,'Color','k');
plot(ax,eq,eq.*1.2,'Color','k','LineStyle','--');
plot(ax,eq,eq./1.2,'Color','k','LineStyle','--');

hold(ax,'off');


%Set legend
lgd=legend(ax,[legItems{:}],...
    figaux.subsz({'P_{el} = const.','T_{surf} - T_{bed} = const.'},6),...
    'Location','southeast');


%Axes limits and labels
ax.XLim=[0,lim];
ax.YLim=[0,lim];

xlabel(ax,figaux.subsz('Measured Nu_{cf} (-)',6));
ylabel(ax,figaux.subsz('Estimated Nu_{cf} (-)',6));


%Text size
ax.FontSize=7;
lgd.FontSize=7;


%Export figure for manuscript
fname=[dirFigs,filesep,'Figure',num2str(figidx)];

til.Units='centimeters';
til.InnerPosition=[1.5,1,8.2,8.2];

fig.Units=til.Units;
fig.Position(3:4)=til.OuterPosition(3:4)+1;

exportgraphics(fig,[fname,'.tiff'],'Resolution',600);


%Export figure for Elsevier
til.InnerPosition=[1.5,1,8.2,8.2];
fig.Position(3:4)=til.OuterPosition(3:4)+1;

exportgraphics(fig,[fname,'.eps']);
savefig(fig,fname);


%% Control strategy analysis
%Set up table
names={'Strategy','mode','Tsurf','Tbed','DeltaT','P_el',...
    'Nu_cf','ME'};
strat={'Pel=const.';'Tsurf-Tbed=const.'};
contr=table('Size',[length(strat),length(names)],...
    'VariableTypes',[{'string','logical'},repmat({'double'},1,length(names)-2)],...
    'VariableNames',names);


%Fill values
contr.Strategy=strat;
contr.mode(2)=true;

contr.Tsurf=arrayfun(@(tf) mean(prim.Tsurf(prim.mode==tf)),contr.mode);
contr.Tbed=arrayfun(@(tf) mean(prim.Tbed(prim.mode==tf)),contr.mode);
contr.P_el=arrayfun(@(tf) mean(prim.P_el(prim.mode==tf)),contr.mode);
contr.Nu_cf=arrayfun(@(tf) mean(y(prim.mode==tf)),contr.mode);

contr.ME=arrayfun(@(tf) mean(yEst(prim.mode==tf)-y(prim.mode==tf)),contr.mode);

contr.DeltaT=contr.Tsurf-contr.Tbed;


%Mean relative impact of Nusselt cross flow
Nu_cfRel=contr.Nu_cf(2)./contr.Nu_cf(1);


%% Auxiliary functions
function monostr=monostr(i)
    %Creates a monomial function string for pi-factor indices i

    monocell=compose('x(:,%d).^b()',i);
    monostr=strjoin(monocell,'.*');
end


function [str,n]=bidx(str)
    %Adds indices to the regression coefficient b in order of appearance
    %Output n: number of regression coefficients

    idx=strfind(str,'b()');
    n=numel(idx);
    for i=1:n
        str=[str(1:idx(i)+1),num2str(i),str(idx(i)+2:end)];
        idx=idx+1;
    end
end


function [str,n]=couple(str,i)
    %Couples the exponents of pi-factors i in the function string str
    %Output n: number of regression coefficients


    %Search string for regular expression: b-index of pi-factor x
    s=@(x) ['(?<=x\(:,',num2str(x),'\)\)?\.\^\(?b\()\d+'];


    %Search for first appearance of first b-index
    idx=regexp(str,s(i(1)));
    idx=str(idx);


    %Replace every b-index of coupled pi-factors with first b-index
    for k=i
        str=regexprep(str,s(k),idx);
    end


    %Number of regression coefficients remaining
    i=regexp(str,'b\((\d*)\)','tokens');
    i=[i{:}];
    n=max(str2double(i));
end




