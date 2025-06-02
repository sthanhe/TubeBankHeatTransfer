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
% This script conducts a basic analysis of the collected secondary data and
% saves the results for the main analysis conducted by the script "calcPC".
% It also creates Figures published in the Methodology Report in the data
% repository.
%
%
%Requires all auxiliary classes and functions on the MATLAB path
%
%Required products, version 24.1:
%   - MATLAB
%   - Curve Fitting Toolbox
%Necessary classes, functions, files, and scripts:
%   - @Al2O3
%   - @DryAir
%   - @FluBed
%   - @SLglass
%   - @SiO2
%   - @implExp
%   - covplot.m


%% Set data locations
dirData='Data';                     %Data storage folder
dirEder=[dirData,filesep,'Eder'];   %Storage folder of data from Eder
fname='secData.csv';                %Secondary data file

dirFigs=['Figures',filesep,'SecData'];     %Regular figure storage folder
dirUni=[dirFigs,filesep,'Univariate'];      %Univariate plots storage folder
dirBi=[dirFigs,filesep,'Bivariate'];        %Bivariate plots storage folder

dirTabs='Tables';   %Table storage folder


%% Make folders if they do not exist
if ~isfolder(dirFigs)
    mkdir(dirFigs);
end

if ~isfolder(dirUni)
    mkdir(dirUni);
end

if ~isfolder(dirBi)
    mkdir(dirBi);
end

if ~isfolder(dirTabs)
    mkdir(dirTabs);
end


%% Eder: calculate means of stationary states
%Get file names
foldEder=dirEder;

files=dir(foldEder);
files=files(endsWith({files.name},'.csv'));
fnames={files.name};


%Initialize table for mean values
eder=readtable([foldEder,filesep,fnames{1}]);
ncol=size(eder,2);
eder=table('Size',[length(fnames),ncol],...
    'VariableTypes',[{'datetime'},repmat({'double'},1,ncol-1)],...
    'VariableNames',eder.Properties.VariableNames);


%Read files, remove outliers, calculate means
for i=1:length(fnames)
    tabloc=readtable([foldEder,filesep,fnames{i}]);

    outliers=isoutlier(tabloc{:,{'T_bed','eps','w'}});
    tabloc{any(outliers,2),2:end}=NaN;

    eder.Time(i)=tabloc.Time(1);
    eder{i,2:end}=mean(tabloc{:,2:end},1,'omitmissing');
end


%Correct particle mass flows (rounding issue)
eder.mDot_p=round(eder.mDot_p.*60^2)./60^2;


%% Eder: mean bed temperatures
%Set up table for all recorded particle mass flows
mDot_p=unique(eder.mDot_p);
varnames={'mDot_p','n','T_bed'};
ederTempsAll=table('Size',[numel(mDot_p),length(varnames)],...
    'VariableTypes',repmat({'double'},1,length(varnames)),...
    'VariableNames',varnames);


%Calculate mean bed temperatures for all recorded particle mass flows (in
%kg/h)
ederTempsAll.mDot_p=mDot_p.*60^2;
ederTempsAll.n=arrayfun(@(i) nnz(eder.mDot_p==mDot_p(i)),1:numel(mDot_p))';
ederTempsAll.T_bed=arrayfun(@(i) mean(eder.T_bed(eder.mDot_p==mDot_p(i))),1:numel(mDot_p))';


%Approximate (weighted) mean bed temperatures at reported particle mass flows 
mDot_p=[0;88;206;323];  %reported particle mass flows  in kg/h
varnames={'mDot_p','T_bed'};
ederTemps=table('Size',[numel(mDot_p),length(varnames)],...
    'VariableTypes',repmat({'double'},1,length(varnames)),...
    'VariableNames',varnames);


ederTemps.mDot_p=mDot_p;
ederTemps.T_bed(mDot_p==0)=ederTempsAll.T_bed(ederTempsAll.mDot_p==0);

idx=ederTempsAll.mDot_p==80 | ederTempsAll.mDot_p==100;
ederTemps.T_bed(mDot_p==88)=mean(ederTempsAll.T_bed(idx),"Weights",ederTempsAll.n(idx));

idx=ederTempsAll.mDot_p==175 | ederTempsAll.mDot_p==200 | ederTempsAll.mDot_p==250;
ederTemps.T_bed(mDot_p==206)=mean(ederTempsAll.T_bed(idx),"Weights",ederTempsAll.n(idx));

idx=ederTempsAll.mDot_p==300 | ederTempsAll.mDot_p==310 | ederTempsAll.mDot_p==400;
ederTemps.T_bed(mDot_p==323)=mean(ederTempsAll.T_bed(idx),"Weights",ederTempsAll.n(idx));


%Save table for repository
ederTempsAll.T_bed=round(ederTempsAll.T_bed,2);
writetable(ederTempsAll,[dirTabs,filesep,'ederTempsAll.csv']);


%% Eder: add mean bed temperatures to secondary data
%Read secondary data
sec=readtable([dirData,filesep,fname]);


%Get particle mass flows in kg/h
mDot_p=sec.mDot_p;
isEder=strcmp(sec.Author,'Eder');
mDot_p(isEder)=round(mDot_p(isEder).*60^2);     %fix rounding issues


%Add mean bed temperatures
for i=1:height(ederTemps)
    idx=isEder & mDot_p==ederTemps.mDot_p(i);
    sec.T_bed(idx)=ederTemps.T_bed(i);
end


%% Calculate missing variables in secondary data
%Particle sphericity phi_s
idx=strcmp(sec.Author,'Kim (2003)') | strcmp(sec.Author,'Kim (2013)');
sec.phi_s(idx)=FluBed.eps2phi(sec.eps_mf(idx));


%Minimum fluidization velocity w_mf
idx=strcmp(sec.Author,'Olsson') | strcmp(sec.Author,'Kim (2013)');
sec.w_mf(idx)=FluBed.wmfErgun(sec.d_p(idx),sec.rho_p(idx),sec.phi_s(idx),...
                sec.eps_mf(idx),sec.p_bed(idx),sec.T_bed(idx));


%Fluidization velocity w
idx=strcmp(sec.Author,'Olsson') | strcmp(sec.Author,'Wiman');
sec.w(idx)=sec.w_e(idx)+sec.w_mf(idx);


idx=strcmp(sec.Author,'Kim (2013)');
sec.w(idx)=sec.FG(idx).*sec.w_mf(idx);


%% Remove values outside of mixed laminar/turbulent regime
%Archimedes number
sec.Ar=FluBed.Ar(sec.d_p,sec.rho_p,sec.p_bed,sec.T_bed);


%Remove values outside of scope
islam=sec.Ar<=1e2;      %Laminar heat transfer regime
isturb=1e5<=sec.Ar;     %Turbulent heat transfer regime

nlam=nnz(islam);    %Number of measurements in the laminar regime
nturb=nnz(isturb);  %Number of measurements in the turbulent regime
ntot=height(sec);   %Total number of collected measurements

sec(islam,:)=[];
sec(islam,:)=[];


%% Calculated Nusselt numbers and HTCs (original model)
%Auxiliary function handle
MW=@(i,cpfx) FluBed.molWirth(sec.w(i),sec.T_bed(i),...
    sec.p_bed(i),sec.d_p(i),sec.rho_p(i),sec.phi_s(i),...
    sec.eps_mf(i),cpfx);


%Calculate Nu and HTCs for each material
mat={'Silica','Soda-lime glass','Alumina'};
cpfx={@SiO2.c_p,@SLglass.c_p,@Al2O3.c_p};
c_p=NaN(height(sec),1);
for i=1:length(mat)
    idx=strcmp(sec.Material,mat{i});
    c_p(idx)=cpfx{i}(sec.T_bed(idx));
    [h,Nu]=MW(idx,cpfx{i});

    sec.h_gcMW(idx)=h.gc;
    sec.h_pcMW(idx)=h.pcMix;
    sec.h_mixMW(idx)=h.mix;

    sec.Nu_gcMW(idx)=Nu.gc;
    sec.Nu_pcMW(idx)=Nu.pcMix;
    sec.Nu_mixMW(idx)=Nu.mix;
end


%% Dimensionless numbers (pi-factors)
npi=10;
pis=table('Size',[height(sec),npi],...
    'VariableTypes',repmat({'double'},1,npi),...
    'VariableNames',compose('pi%d',1:npi));


[pis{:,:},k_g,l_lam]=FluBed.piFactors(sec.w,sec.T_bed,sec.p_bed,sec.d_p,...
    sec.rho_p,sec.phi_s,sec.eps_mf,c_p,sec.d_t,sec.p_h,sec.w_p);


pis{:,1}=sec.h.*l_lam./k_g;


%% Save for future analysis
save('secData','sec','pis');


%% Univariate plots
authors=[unique(sec.Author);'Eder9'];
figidx=801;
for i=authors'
    %Get figure
    switch i{1}
        case 'Eder9'
            idx=strcmp(sec.Author,'Eder') & pis.pi10==0;
        otherwise
            idx=strcmp(sec.Author,i);
    end
    [fig,~,ax]=uniplot(pis{idx,:},figidx);


    %Axes labels and title
    xlabel(ax,'\pi-index (-)');

    switch i{1}
        case 'Eder9'
            title(ax,'Eder, excluding \pi_{10}');
        otherwise
            title(ax,i);
    end


    %Export figure for repository
    exportgraphics(fig,[dirUni,filesep,i{1},'.tiff'],...
        'Resolution',600);


    %Increment figure index
    figidx=figidx+1;
end


%% Bivariate plots
pinames=compose('\\pi_{%d}',1:size(pis,2));
figidx=810;
for i=authors'
    %Get figure
    switch i{1}
        case 'Eder9'
            idx=strcmp(sec.Author,'Eder') & pis.pi10==0;
        otherwise
            idx=strcmp(sec.Author,i);
    end
    [fig,t]=covplot(pis{idx,:},pinames,figidx);


    %Add title to figure
    switch i{1}
        case 'Eder9'
            title(t,'Eder covariance plot excluding \pi_{10}, z-scores');
        otherwise
            title(t,[i{1},' covariance plot, z-scores']);
    end


    %Export figure for repository (full screen)
    fig.Units='normalized';
    fig.Position=[0,0,1,1];

    exportgraphics(fig,[dirBi,filesep,i{1},'.tiff'],...
        'Resolution',600);


    %Increment figure index
    figidx=figidx+1;
end


%% Pi-factor variations
% All 10 pi-factors
[tabMin10,tabMean10,tabMax10,tabVars10]=getVars(sec,pis,dirTabs);


% Without pi10, only mixed regime
idx=sec.Ar<1e5 & pis.pi10==0;
[tabMinMixed,~,tabMaxMixed,tabVarsMixed]=getVars(sec(idx,:),pis(idx,:),[]);


%Eder, without pi10
idx=strcmp(sec.Author,'Eder') & pis.pi10==0;
nEderCF=nnz(strcmp(sec.Author,'Eder') & pis.pi10>0);    %Number of Eder's measurements with cross-flow
[tabMin9,tabMean9,tabMax9,tabVars9]=getVars(sec(idx,:),pis(idx,1:9),dirTabs);


%% Auxiliary functions
%This function creates the tables that demonstrate pi-factor variations
function [tabMin,tabMean,tabMax,tabVars]=getVars(sec,pis,dirTabs)
    %Basic properties
    authors=unique(sec.Author);
    npi=size(pis,2);


    %Create minimum, mean, and maximum tables for Ar and pi-factors
    vars=[{'Author','Ar'},compose('pi%d',1:npi)];
    
    tabMin=table('Size',[length(authors),length(vars)],...
        'VariableTypes',[{'string'},repmat({'double'},1,length(vars)-1)],...
        'VariableNames',vars);
    
    tabMin.Author=authors;
    tabMean=tabMin;
    tabMax=tabMin;
    
    
    %Create table for number of observations and pi-factor variations
    vars=[{'Author','nObs'},compose('pi%d',1:npi)];
    tabVars=table('Size',[length(authors),length(vars)],...
        'VariableTypes',[{'string'},repmat({'logical'},1,length(vars)-1)],...
        'VariableNames',vars);
    
    tabVars.Author=authors;
    tabVars.nObs=cellfun(@(x) nnz(strcmp(sec.Author,x)),authors);
    
    
    %Add Archimedes numbers to min, mean, and max tables
    tabMin.Ar=cellfun(@(x) ...
        min(sec.Ar(strcmp(sec.Author,x))),...
        authors);

    tabMean.Ar=cellfun(@(x) ...
        mean(sec.Ar(strcmp(sec.Author,x)),'omitmissing'),...
        authors);
    
    tabMax.Ar=cellfun(@(x) ...
        max(sec.Ar(strcmp(sec.Author,x))),...
        authors);
    
    
    %Variation of pi-factors
    for i=authors'
        for j=1:npi
            pistr=['pi',num2str(j)];
    
            
            %Check if pi-factor varies or is constant
            vars=unique(pis{strcmp(sec.Author,i),j});
            if length(vars)>1
                tabVars{strcmp(tabVars.Author,i),pistr}=true;
            end
    
            
            %Calculate min, mean, and max of pi-factors
            authLoc=strcmp(tabVars.Author,i);
    
            tabMin{authLoc,pistr}=min(pis{strcmp(sec.Author,i),pistr});
            tabMean{authLoc,pistr}=mean(pis{strcmp(sec.Author,i),pistr},'omitmissing');
            tabMax{authLoc,pistr}=max(pis{strcmp(sec.Author,i),pistr});            
        end
    end


    %Convert numerical values to exponential notation for printing
    for i=2:width(tabMin)
        tabMin.(i)=compose('%.1e',tabMin{:,i});
        tabMean.(i)=compose('%.1e',tabMean{:,i});
        tabMax.(i)=compose('%.1e',tabMax{:,i});
    end


    %Write tables
    if ~isempty(dirTabs)
        writetable(tabMin,[dirTabs,filesep,'tabMin',num2str(npi),'.csv']);
        writetable(tabMean,[dirTabs,filesep,'tabMean',num2str(npi),'.csv']);
        writetable(tabMax,[dirTabs,filesep,'tabMax',num2str(npi),'.csv']);
    end
end




