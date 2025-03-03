fname=['Data',filesep,'ExtLit.csv'];


%% Read table and calculate missing variables
tab=readtable(fname);


% Particle sphericity phi_s
idx=strcmp(tab.Author,'Kim (2013)') | strcmp(tab.Author,'Kim (2003)');
tab.phi_s(idx)=FluBed.eps2phi(tab.eps_mf(idx));


%Minimum fluidization velocity w_mf
idx=strcmp(tab.Author,'Olsson') | strcmp(tab.Author,'Hofer') | ...
    strcmp(tab.Author,'Kim (2013)');
tab.w_mf(idx)=FluBed.wmfErgun(tab.d_p(idx),tab.rho_p(idx),tab.phi_s(idx),...
                tab.eps_mf(idx),tab.p(idx),tab.T(idx));


%Fluidization velocity w
idx=strcmp(tab.Author,'Olsson') | strcmp(tab.Author,'Wiman');
tab.w(idx)=tab.w_e(idx)+tab.w_mf(idx);


idx=strcmp(tab.Author,'Kim (2013)');
tab.w(idx)=tab.FG(idx).*tab.w_mf(idx);


%% Remove values outside of scope
%HTCs in turbulent regime
Ar=FluBed.Ar(tab.d_p,tab.rho_p,tab.p,tab.T);
tab(1e5<=Ar,:)=[];


%HTCs with particle cross-flow
% tab(tab.w_p~=0,:)=[];


%No Hofer
tab(strcmp(tab.Author,'Hofer'),:)=[];


%% HTCs according to Molerus
tab.c_p=NaN(height(tab),1);
molerus=@(i,cpfx) FluBed.molerus(tab.w(i),tab.T(i),...
    tab.p(i),tab.d_p(i),tab.rho_p(i),tab.phi_s(i),...
    tab.eps_mf(i),cpfx);

molExt=@(i,cpfx) FluBed.molExt(tab.w(i),tab.T(i),...
    tab.p(i),tab.d_p(i),tab.rho_p(i),tab.phi_s(i),...
    tab.eps_mf(i),cpfx,tab.d_t(i),tab.p_h(i),tab.w_p(i));


mat={'Silica','Soda-lime glass','Alumina'};
cpfx={@SiO2.c_p,@SLglass.c_p,@Al2O3.c_p};
for i=1:length(mat)
    idx=strcmp(tab.Material,mat{i});
    tab.c_p(idx)=cpfx{i}(tab.T(idx));
    [h,Nu]=molerus(idx,cpfx{i});

    tab.h_gcMol(idx)=h.gc;
    tab.h_pcMol(idx)=h.pcMix;
    tab.h_mixMol(idx)=h.mix;

    tab.Nu_gcMol(idx)=Nu.gc;
    tab.Nu_pcMol(idx)=Nu.pcMix;
    tab.Nu_mixMol(idx)=Nu.mix;


    [h,Nu]=molExt(idx,cpfx{i});

    tab.hExt(idx)=h;
    tab.NuExt(idx)=Nu;
end


%% Dimensionless numbers
%Gas and particle properties
tab.c_g=DryAir.c_p(tab.T);
k_g=DryAir.lambda(tab.T);
my_g=DryAir.eta(tab.T);
rho_g=DryAir.rho(tab.p,tab.T);
rho_e=tab.rho_p-rho_g;
l_l=(my_g./(rho_e.*sqrt(FluBed.g))).^(2/3);


%Fluidization velocities
w_mf=FluBed.wmfErgun(tab.d_p,tab.rho_p,tab.phi_s,tab.eps_mf,tab.p,tab.T);
w_e=tab.w-w_mf;
w_e(w_e<0)=NaN;


%Pi-factors
npi=10;
pis=table('Size',[height(tab),npi],...
    'VariableTypes',repmat({'double'},1,npi),...
    'VariableNames',compose('pi%d',1:npi));

pis.pi1=tab.h.*l_l./k_g;
pis.pi2=k_g./(2*tab.c_p.*my_g);
pis.pi3=DryAir.Pr(tab.T);
pis.pi4=rho_g./rho_e;
pis.pi5=(rho_e.*tab.c_p./(k_g.*FluBed.g)).^(1/3).*w_e;
pis.pi6=(rho_e.*tab.c_p./(k_g.*FluBed.g)).^(1/3).*w_mf;
pis.pi7=1-tab.eps_mf;
pis.pi8=tab.d_t./l_l;
pis.pi9=tab.d_t./tab.p_h;
pis.pi10=(rho_e.*tab.c_p./(k_g.*FluBed.g)).^(1/3).*tab.w_p;


%% Save for future analysis
save('extLit','pis','tab');


%% Univariate plots
authors=unique(tab.Author);
figidx=100;
for i=authors'
    idx=find(strcmp(tab.Author,i));
    tabloc=tab(idx,:);
    pisNorm=normalize(pis(idx,:));


    %Univariate plot
    fig=figure(figidx);
    clf(fig);
    ax=gca();

    boxplot(ax,pisNorm{:,:});

    xlabel(ax,'\pi-index (-)');
    ylabel(ax,'z-score (-)');
    title(ax,i);


    %Increment figure index
    figidx=figidx+1;
end


%% Bivariate plots
% covplot(pis{:,:},compose('\\pi_{%d}',1:size(pis,2)),200);


%% Pi-factor variations without pi10
idx=pis.pi10==0;
pis2=pis(idx,1:9);
tab2=tab(idx,:);
npi=9;

authors=unique(tab2.Author);
paras=[{'Author','nObs'},compose('pi%d',1:npi)];
extLit=table('Size',[length(authors),length(paras)],...
    'VariableTypes',[{'string'},repmat({'double'},1,length(paras)-1)],...
    'VariableNames',paras);

extLit.Author=authors;
extLit.nObs=cellfun(@(x) nnz(strcmp(tab2.Author,x)),authors);


%Variation of pi-factors
for i=authors'
    for j=1:npi
        vars=unique(pis2{strcmp(tab2.Author,i),j});
        if length(vars)>1
            extLit{strcmp(extLit.Author,i),['pi',num2str(j)]}=1;
        end
    end
end


%% Weigh observations of each author equally
% authors=unique(tab.Author);
% n=cellfun(@(x) nnz(strcmp(tab.Author,x)),authors);
% nMax=lcm(n);
% 
% pisCell=arrayfun(@(i) ...
%     repmat(pis(strcmp(tab.Author,authors(i)),:),...
%         nMax./n(i),1),...
%     1:numel(n),...
%     'UniformOutput',false);
% 
% pis=vertcat(pisCell{:});
% pis=normalize(pis);
