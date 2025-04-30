%% Set data directories
dirData=['Data',filesep,'Own'];             %Path to directory containing the data
dirFigures='Figures';       %Path to directory where figures should be stored
dirTabs='Tables';       %Path to directory where tables should be stored

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
names={'Run','T1','T3','Tsurf','Tbed',...
    'Tco2','p','w','FG','FG2',...
    'w_p','mDot_p','hVirt','P_el','mode'};
prim=table('Size',[length(fnames),length(names)],...
    'VariableTypes',[repmat({'double'},1,length(names)-1),'logical'],...
    'VariableNames',names);


%% Read individual files and do calculations
chambers=1:6;
for i=1:length(fnames)
    %Get properties
    tab=readtable([dirData,filesep,fnames{i}]);
    TT=getProp(tab,c,prim.Properties.VariableNames(2:end),chambers);
    
    
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
Tbedcat=center(prim.Tbed-273.15,[330,400,450]-273.15);
Tbedcat=unique(Tbedcat);
Tbedcat=round(Tbedcat/5)*5;

FGcat=center(prim.FG,[2.5,3.1,3.7]);
FGcat=round(unique(FGcat),1);

w_pcat=unique(testmat.w_p);


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


%% Estimated bias due to fins
%Particle and gas convective HTC according to the extended model
[h,Nu]=FluBed.molExt(prim.w,prim.Tbed,prim.p,c.d_p,c.rho_p,c.phi_s,c.eps_mf,...
    @SiO2.c_p,c.d_t,c.p_hEff,0);

prim.h_pcExt=h.pc;
prim.h_gcExt=h.gc;


%Fit relation between residual heat transfer (=cross-flow HTC) and particle
%velocity
X=prim.w_p;                 %Regressor: particle velocity
y=prim.h_eff-h.pc-h.gc;     %Response: residual HTC

fx=@(b,x) b(1)./(x+b(2));           %Model: rat01
beta0=[-2,0.02];                    %Starting values (estimates)
mdl=fitnlm(X,y,fx,beta0);           %Fit model
beta=mdl.Coefficients.Estimate;     %Estimated coefficients

bias=predict(mdl,0);                    %Estimated bias
prim.h_cf=prim.h_eff-h.pc-h.gc-bias;    %Estimated cross-flow HTC


%Set up figure
fig=figure(917);
clf(fig);
t=tiledlayout(fig,1,1,'Padding','tight');
ax=nexttile(t);
colors=ax.ColorOrder;
hold(ax,'on');


%Plot data and lines
x=linspace(0,max(prim.w_p),1000);
scatter(ax,X,prim.h_eff-h.pc-h.gc,18,colors(1,:));
plot(ax,x,predict(mdl,x'),'Color',colors(1,:));

scatter(ax,X,prim.h_cf,18,colors(2,:));
plot(ax,x,predict(mdl,x')-bias,'Color',colors(2,:));


%Set up legend
legItems=repmat(line(ax,'Visible','off'),2,1);
legItems(1)=plot(ax,NaN,NaN,'Color',colors(2,:),'Marker','o','MarkerSize',sqrt(18));
legItems(2)=plot(ax,NaN,NaN,'Color',colors(1,:),'Marker','o','MarkerSize',sqrt(18));

hold(ax,'off');


%Set legend, axes labels, and title
legend(ax,legItems,{'h_{cf}','h_{eff} - h_{pcExt} - h_{gcExt}'},...
    'Location','southeast');

xlabel(ax,'w_p (m/s)');
ylabel(ax,'HTC (W/m²K)');

title(ax,compose('R² = %.4f',mdl.Rsquared.Ordinary));


%Size figure for repository
t.Units='centimeters';
t.OuterPosition=[0,0,17,8.5];

fig.Units=t.Units;
fig.Position(3:4)=t.OuterPosition(3:4)+0.5;


%Add arrow
x=[0,0];
y=[0,bias];
figaux.arrow(ax,x,y,'doublearrow');
text(ax,x(1),mean(y),compose(' bias = %.0f W/m²K',bias));


%Export figure
exportgraphics(fig,[dirFigures,filesep,'finBias.tiff'],'Resolution',600);


%% Dimensionless numbers
npi=10;
pis=table('Size',[height(prim),npi],...
    'VariableTypes',repmat({'double'},1,npi),...
    'VariableNames',compose('pi%d',1:npi));


[pis{:,:},k_g,l_lam]=FluBed.piFactors(prim.w,prim.Tbed,prim.p,c.d_p,c.rho_p,...
    c.phi_s,c.eps_mf,SiO2.c_p(prim.Tbed),c.d_t,c.p_hEff,prim.w_p);


pis{:,1}=prim.h_cf.*l_lam./k_g;


%% Save for future analysis
save('primData','prim','pis');


%% Univariate plot
%Set directory and create it if it does not exist
dirPrim=[dirFigures,filesep,'PrimData'];

if ~isfolder(dirPrim)
    mkdir(dirPrim);
end


%Normalize (z-score)
pisNorm=normalize(pis);


%Set constants to 0
isconst=arrayfun(@(i) ...
    isscalar(unique(pis(:,i))),...
    1:size(pisNorm,2));
pisNorm{:,isconst}=0;


%Set up figure
fig=figure(918);
clf(fig);
t=tiledlayout(fig,1,1,'Padding','tight');
ax=nexttile(t);

boxchart(ax,pisNorm{:,:});


%Axes labels and title
xlabel(ax,'\pi-index (-)');
ylabel(ax,'z-score (-)');


%Export figure for repository
t.Units='centimeters';
t.OuterPosition=[0,0,17,8.5];

fig.Units=t.Units;
fig.Position(3:4)=t.OuterPosition(3:4)+0.5;

exportgraphics(fig,[dirPrim,filesep,'Univariate.tiff'],...
    'Resolution',600);


%% Bivariate plot
fig=covplot(pis{:,:},compose('\\pi_{%d}',1:npi),919);


%Export figure for repository (full screen)
fig.Units='normalized';
fig.Position=[0,0,1,1];

exportgraphics(fig,[dirPrim,filesep,'Bivariate.tiff'],...
    'Resolution',600);




