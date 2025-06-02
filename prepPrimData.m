%% Prepare analysis of primary data
% GNU General Public License v3.0
% By Stefan Thanheiser: https://orcid.org/0000-0003-2765-1156
%
% Part of the paper:
%
% Thanheiser, S.; Haider, M.
% Molerus and Wirth's Heat Transfer Model for Bubbling Fluidized Beds: 
% Proposal for an Extended Model Including Immersed Tube Banks and Particle 
% Cross-Flow
%
% All data, along with methodology reports and supplementary documentation, 
% is published in the data repository:
% https://doi.org/10.5281/zenodo.15576311
%
% All required files for this script can be found in the software
% repository: see the link to the supplemental release in the data 
% repository
%
%
%
% This script conducts a basic analysis of the collected primary data and
% saves the results for the main analysis conducted by the script "calcCF".
% It also creates Figures published in the Methodology Report in the data
% repository.
%
%
%Requires all auxiliary classes and functions on the MATLAB path
%
%Required products, version 24.1:
%   - MATLAB
%   - Statistics and Machine Learning Toolbox
%   - Curve Fitting Toolbox
%Necessary classes, functions, files, and scripts:
%   - @DC04
%   - @DryAir
%   - @FluBed
%   - @Orifice
%   - @SiO2
%   - @figaux
%   - @implExp
%   - center.m
%   - covplot.m
%   - getConstants.m
%   - getProp.m
%   - h2FG.mat --> created by the script "prepFG" 


%% Set data locations
dirData=['Data',filesep,'Own'];             %Data storage folder
dirFigs=['Figures',filesep,'PrimData'];     %Figure storage folder
dirTabs='Tables';                           %Table storage folder


%% Make folders if they do not exist
if ~isfolder(dirFigs)
    mkdir(dirFigs);
end

if ~isfolder(dirTabs)
    mkdir(dirTabs);
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
names={'Run','T1','T3','Tsurf','Tbed',...
    'Tco2','p','w','FG','FG2',...
    'w_p','mDot_p','hVirt','P_el','mode'};
prim=table('Size',[length(fnames),length(names)],...
    'VariableTypes',[repmat({'double'},1,length(names)-1),'logical'],...
    'VariableNames',names);


%% Read individual files and do calculations
for i=1:length(fnames)
    %Get properties
    tab=readtable([dirData,filesep,fnames{i}]);
    TT=getProp(tab,c,prim.Properties.VariableNames(2:end),1:6);
    
    
    %Remove outliers
    outliers=isoutlier(TT.hVirt,1);
    TT{outliers,2:end-1}=NaN;
    
    
    %Get means
    prim{i,2:width(TT)-1}=mean(TT{:,2:end-1},1,'omitnan');
    prim.Run(i)=i;
    prim.mode(i)=nnz(TT.mode)>height(prim)/2;
end


%% Test matrix
%Set up and fill table
names={'Run','Tco2','FG2','mDot_p','mode','Tbed','FG','w_p'};
testmat=table('Size',[height(prim),length(names)],...
    'VariableTypes',repmat({'double'},1,length(names)),...
    'VariableNames',names);

testmat.Run=prim.Run;

testmat.Tco2=round(center(prim.Tco2,[350,460,570]),-1);
testmat.FG2=round(center(prim.FG2,[2.5,3.1,3.7]),1);
testmat.mDot_p=prim.mDot_p;
testmat.mode=prim.mode;

testmat.Tbed=round(prim.Tbed);
testmat.FG=round(prim.FG,1);
testmat.w_p=round(prim.w_p,4);


%Categorize parameters
TbedCat=center(prim.Tbed-273.15,[330,400,450]-273.15);
TbedCat=unique(TbedCat);
TbedCat=round(TbedCat/5)*5;

FGcat=center(prim.FG,[2.5,3.1,3.7]);
FGcat=round(unique(FGcat),1);

w_pCat=unique(testmat.w_p);


%Print test matrix
testmat.FG2=compose('%.1f',testmat.FG2);
testmat.mDot_p=compose('%.2f',testmat.mDot_p);
testmat.FG=compose('%.1f',testmat.FG);
testmat.w_p=compose('%.4f',testmat.w_p);
writetable(testmat,[dirTabs,filesep,'testmat.csv']);


%Remove unusable results
outliers=find(prim.hVirt<5);
prim(outliers,:)=[];


%% Effective HTC
prim.h_eff=prim.hVirt;
prim.eta_f=ones(height(prim),1);
k=55;   %Thermal conductivity of fin material
for i=1:height(prim)
    err=1;
    counter=0;
    deltaT=prim.Tsurf(i)-prim.Tbed(i);
    while err>1e-6 && counter<100
        X=c.phi.*c.d_t/2.*sqrt(2*prim.h_eff(i)./(k*c.s_f));
        prim.eta_f(i)=tanh(X)./X;
        A_eff=c.A_plain-c.A_bottom+prim.eta_f(i).*c.A_sides;
        
        h_effNew=prim.P_el(i)./(A_eff.*deltaT);

        T_f=prim.eta_f(i).*deltaT+prim.Tbed(i);
        k=DC04.lambda(T_f);

        err=abs(prim.h_eff(i)-h_effNew);
        prim.h_eff(i)=h_effNew;

        counter=counter+1;
    end
end


%% Dimensionless numbers
npi=10;
pis=table('Size',[height(prim),npi],...
    'VariableTypes',repmat({'double'},1,npi),...
    'VariableNames',compose('pi%d',1:npi));


[pis{:,:},k_g,l_lam]=FluBed.piFactors(prim.w,prim.Tbed,prim.p,c.d_p,c.rho_p,...
    c.phi_s,c.eps_mf,SiO2.c_p(prim.Tbed),c.d_t,c.p_hEff,prim.w_p);


%% Estimated bias due to fins
%Regressor: dimensionless mean horizontal particle velocity pi10
X=pis.pi10;


%Particle and gas convective HTC according to the extended model
[h,Nu]=FluBed.extended(prim.w,prim.Tbed,prim.p,c.d_p,c.rho_p,c.phi_s,c.eps_mf,...
    @SiO2.c_p,c.d_t,c.p_hEff,0);


%Effective Nusselt number
Nu_eff=prim.h_eff.*l_lam./k_g;


%Response variable: residual heat transfer (=cross-flow Nu)
y=Nu_eff-Nu.pc-Nu.gc;


%Bin and center regressor and response variable
[Xc,Xm,n]=center(X,[0.25,0.6,1.2,2.5]);
ym=arrayfun(@(i) mean(y(Xc==i)),Xm);


%Fit function to centered values
fx=@(b,x) b(1)./(x+b(2));   %rat01
beta0=[-0.1,5];             %Starting values (estimates)

mdl=fitnlm(Xm,ym,fx,beta0,'Weights',n);     %Fit model
beta=mdl.Coefficients.Estimate;             %Estimated coefficients


%Cross-flow Nusselt number
bias=predict(mdl,0);                    %Estimated bias
prim.Nu_cf=Nu_eff-Nu.pc-Nu.gc-bias;     %Estimated cross-flow Nu


%Set up figure
fig=figure(917);
clf(fig);
t=tiledlayout(fig,1,1,'Padding','tight');
ax=nexttile(t);
colors=ax.ColorOrder;
hold(ax,'on');


%Plot data and lines
x=linspace(0,max(X),1000);

scatter(ax,X,y,18,colors(1,:));
plot(ax,x,predict(mdl,x'),'Color',colors(1,:));

scatter(ax,X,prim.Nu_cf,18,colors(2,:));
plot(ax,x,predict(mdl,x')-bias,'Color',colors(2,:));

scatter(ax,Xm,ym,72,colors(4,:),'Marker','x','LineWidth',2);


%Set up legend
legItems=repmat(line(ax,'Visible','off'),2,1);
legItems(1)=plot(ax,NaN,NaN,'Color',colors(2,:),'Marker','o','MarkerSize',sqrt(18));
legItems(2)=plot(ax,NaN,NaN,'Color',colors(1,:),'Marker','o','MarkerSize',sqrt(18));

hold(ax,'off');


%Set legend, axes labels, and title
legend(ax,legItems,{'Nu_{cf}','Nu_{eff} - Nu_{pcExt} - Nu_{gc}'},...
    'Location','southeast');

xlabel(ax,'\pi_{10} (-)');
ylabel(ax,'Nu_{cf} (-)');

title(ax,compose('R² = %.4f',mdl.Rsquared.Ordinary));


%Size figure for repository
t.Units='centimeters';
t.OuterPosition=[0,0,17,8.5];

fig.Units=t.Units;
fig.Position(3:4)=t.OuterPosition(3:4)+0.5;


%Add arrow
xArr=[0,0];
yArr=[0,bias];
figaux.arrow(ax,xArr,yArr,'doublearrow');
text(ax,x(1),bias./2,compose(' bias = %.1e',bias));


%Export figure
exportgraphics(fig,[dirFigs,filesep,'finBias.tiff'],'Resolution',600);


%% Relative contribution of particle cross-flow
%Maximum contribution
maxContr=max(prim.Nu_cf./(Nu_eff-predict(mdl,0)));


%Set up figure
fig=figure(918);
clf(fig);
t=tiledlayout(fig,1,1,'Padding','tight');
ax=nexttile(t);
colors=ax.ColorOrder;
hold(ax,'on');


%Plot data and lines
pos=prim.Nu_cf>0;
scatter(ax,prim.Nu_cf(pos),prim.Nu_cf(pos)./(Nu_eff(pos)-predict(mdl,0)));
scatter(ax,prim.Nu_cf(~pos),prim.Nu_cf(~pos)./(Nu_eff(~pos)-predict(mdl,0)));

xline(ax,0);
yline(ax,0);

hold(ax,'off');


%Set axes labels
xlabel(ax,'Nu_{cf} (-)');
ylabel(ax,'Nu_{cf} / (Nu_{eff} - bias) (-)');


%Size figure for repository
t.Units='centimeters';
t.OuterPosition=[0,0,17,8.5];

fig.Units=t.Units;
fig.Position(3:4)=t.OuterPosition(3:4)+0.5;


%Export figure
exportgraphics(fig,[dirFigs,filesep,'crossFlowContr.tiff'],...
    'Resolution',600);


%% Remove outliers, save for future analysis
%Outliers=negative Nu_cf
neg=find(~pos);
prim(neg,:)=[];
pis(neg,:)=[];


%Correct outlier indices of test matrix by previously removed outliers
for i=flip(outliers)
    neg(neg>=i)=neg(neg>=i)+1;
end


%Add pi1=Nu_cf to pi-table
pis.pi1=prim.Nu_cf;


%Save data
save('primData','prim','pis');


%% Univariate plot
[fig,~,ax]=uniplot(pis{:,:},919);


%x-axis label
xlabel(ax,'\pi-index (-)');


%Export figure for repository
exportgraphics(fig,[dirFigs,filesep,'Univariate.tiff'],...
        'Resolution',600);


%% Bivariate plot
fig=covplot(pis{:,:},compose('\\pi_{%d}',1:npi),920);


%Export figure for repository (full screen)
fig.Units='normalized';
fig.Position=[0,0,1,1];

exportgraphics(fig,[dirFigs,filesep,'Bivariate.tiff'],...
    'Resolution',600);




