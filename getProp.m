%% Calculate Basic Properties
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
% This function calculates the basic properties of the fluidized bed system
% needed for the subsequent analysis of measurements.
%
%
%Requires all auxiliary classes and functions on the MATLAB path
%
%Required products, version 24.1:
%   - MATLAB
%   - Statistics and Machine Learning Toolbox
%Necessary files, classes, functions, and scripts:
%   - @DryAir
%   - @FluBed
%   - @Orifice
%   - @implExp
%   - h2FG.mat --> created by the script "prepFG" 


function main=getProp(tab,c,htcNames,chambers)
    % Inputs:
    % tab: raw data measurements, table
    % c: constants of the test rig (see function "getConstants"), struct
    % htcNames: variable names of the output table
    % chambers: indices of chambers to include (between 1 and 6; indices 
    % are identical to the bed height indices), double


    %% Air flows
    nOrif=8;    %Number of orifice plates
    Onames=compose('O%d',1:nOrif);
    OscalNames=compose('O%dscal',1:nOrif);

    names=[{'Time','mDot','On'},...
        Onames,OscalNames,...
        compose('mDot%d',1:4)];
    air=table('Size',[height(tab),length(names)],...
            'VariableTypes',[{'datetime','double','logical'},repmat({'double'},1,length(names)-3)]);
    air.Properties.VariableNames=names;
    air.Time=tab.Time;
    clear('names');
    
    
    %Total air flow (anemometer). Ignore negative values
    %Used as reference for all other air flows
    air.mDot=max(zeros(height(tab),1),tab.w1*c.rho_N*c.Apipe);  
    air.On=air.mDot>c.OnLimit;  %Indicator whether air supply is running
    
    
    %Orifice plates: measured air mass flows, ignore negative values
    pos=tab.p10>=0;
    air.O1(pos)=Orifice.qm(tab.p20(pos),tab.p20(pos)-tab.p10(pos),tab.T8(pos),c.dOrif,c.DOrif,c.tap);
    
    pos=tab.p11>=0;
    air.O2(pos)=Orifice.qm(tab.p20(pos),tab.p20(pos)-tab.p11(pos),tab.T8(pos),c.dOrif,c.DOrif,c.tap);
    
    pos=tab.p12>=0;
    air.O3(pos)=Orifice.qm(tab.p20(pos),tab.p20(pos)-tab.p12(pos),tab.T8(pos),c.dOrif,c.DOrif,c.tap);
    
    pos=tab.p13>=0;
    air.O4(pos)=Orifice.qm(tab.p20(pos),tab.p20(pos)-tab.p13(pos),tab.T8(pos),c.dOrif,c.DOrif,c.tap);
    
    pos=tab.p14>=0;
    air.O5(pos)=Orifice.qm(tab.p20(pos),tab.p20(pos)-tab.p14(pos),tab.T8(pos),c.dOrif,c.DOrif,c.tap);
    
    pos=tab.p15>=0;
    air.O6(pos)=Orifice.qm(tab.p20(pos),tab.p20(pos)-tab.p15(pos),tab.T8(pos),c.dOrif,c.DOrif,c.tap);
    
    pos=tab.p16>=0;
    air.O7(pos)=Orifice.qm(tab.p20(pos),tab.p20(pos)-tab.p16(pos),tab.T8(pos),c.dOrif,c.DOrif,c.tap);
    
    pos=tab.p17>=0;
    air.O8(pos)=Orifice.qm(tab.p20(pos),tab.p20(pos)-tab.p17(pos),tab.T8(pos),c.dOrif,c.DOrif,c.tap);
    clear('pos');
    
    
    %Orifice plates: Mass flows scaled to total air mass flow (anemometer)
    S=air.mDot./sum(air{:,Onames},2);
    air{:,OscalNames}=air{:,Onames}.*S;
    clear('S');
    
    
    %Aggregated mass flows
    air.mDot1=air.O8scal;
    air.mDot2=air.O5scal+air.O6scal+air.O7scal;
    air.mDot3=air.O2scal+air.O3scal+air.O4scal;
    air.mDot4=air.O1scal;
    
    
    
    %% Bed properties
    hNames=compose('h%d',chambers);    
    names=[{'Time'},...
            compose('eps%d',1:3),...
            hNames];
    bed=table('Size',[height(tab),length(names)],...
            'VariableTypes',[{'datetime'},repmat({'double'},1,length(names)-1)]);
    bed.Properties.VariableNames=names;
    bed.Time=tab.Time;
    clear('names');
    
    
    %Bed porosities: only calculate values when air supply is on
    bed.eps1(air.On)=FluBed.eps(tab.p1(air.On),c.dh_eps1,c.rho_p);
    % bed.eps2(air.On)=FluBed.eps(tab.p2(air.On),c.dh_eps2,c.rho_p);
    bed.eps3(air.On)=FluBed.eps(tab.p3(air.On),c.dh_eps1,c.rho_p);
    
    %Assume minimum fluidization porosities when air supply is not on
    bed.eps1(~air.On)=c.eps_mf;
    % bed.eps2(~air.On)=c.eps_mf;
    bed.eps3(~air.On)=c.eps_mf;

    bed.eps1(:)=mean(bed.eps1);
    bed.eps3(:)=mean(bed.eps3);
    
    %Fix output of bed porosity 2: assume mean between porosities 1 and 3
    bed.eps2=mean([bed.eps1,bed.eps3],2);
    
    
    %Bed levels: ignore negative bed heights
    bed.h1=max(zeros(height(tab),1),FluBed.h(tab.p4,bed.eps1,c.rho_p));
    bed.h2=max(zeros(height(tab),1),FluBed.h(tab.p5,bed.eps1,c.rho_p));
    bed.h3=max(zeros(height(tab),1),FluBed.h(tab.p6,bed.eps2,c.rho_p));
    bed.h4=max(zeros(height(tab),1),FluBed.h(tab.p7,bed.eps2,c.rho_p));
    bed.h5=max(zeros(height(tab),1),FluBed.h(tab.p8,bed.eps3,c.rho_p));
    bed.h6=max(zeros(height(tab),1),FluBed.h(tab.p9,bed.eps3,c.rho_p));


    %% Fluidization
    names=[{'Time'},...
        compose('pBed%d',1:4),...
        compose('Tbed%d',1:4),...
        compose('rho_g%d',1:4),...
        compose('w%d',1:4),...
        compose('wmf%d',1:4),...
        compose('FG%d',1:4)];
    flu=table('Size',[height(tab),length(names)],'VariableTypes',[{'datetime'},repmat({'double'},1,length(names)-1)]);
    flu.Properties.VariableNames=names;
    flu.Time=tab.Time;
    clear('names');

    
    %Chamber 1: inlet chamber
    flu.pBed1=tab.p21+FluBed.deltaP((bed.h6+c.hBed)./2,bed.eps3,c.rho_p);
    flu.Tbed1=tab.T7;
    flu.rho_g1=DryAir.rho(flu.pBed1,flu.Tbed1);
    flu.w1=air.mDot1./(flu.rho_g1.*c.Afloor1);
    flu.wmf1=FluBed.wmf(c.d_p,c.rho_p,flu.pBed1,flu.Tbed1);
    flu.FG1=flu.w1./flu.wmf1;


    %Chamber 2
    flu.pBed2=tab.p21+tab.p19+FluBed.deltaP((mean([bed.h4,bed.h5],2)+c.hBed)./2,...
            mean([bed.eps2,bed.eps3],2),c.rho_p);
    flu.Tbed2=mean([tab.T7,tab.T1],2);
    flu.rho_g2=DryAir.rho(flu.pBed2,flu.Tbed2);
    flu.w2=air.mDot2./(flu.rho_g2.*c.Afloor2);
    flu.wmf2=FluBed.wmf(c.d_p,c.rho_p,flu.pBed2,flu.Tbed2);
    flu.FG2=flu.w2./flu.wmf2;


    %Chamber 3
    flu.pBed3=tab.p21+tab.p18+FluBed.deltaP((mean([bed.h2,bed.h3],2)+c.hBed)./2,...
            mean([bed.eps1,bed.eps2],2),c.rho_p);
    flu.Tbed3=mean([tab.T6,tab.T1],2);
    flu.rho_g3=DryAir.rho(flu.pBed3,flu.Tbed3);
    flu.w3=air.mDot3./(flu.rho_g3.*c.Afloor3);
    flu.wmf3=FluBed.wmf(c.d_p,c.rho_p,flu.pBed3,flu.Tbed3);
    flu.FG3=flu.w3./flu.wmf3;
    
    
    %Chamber 4: outlet chamber
    flu.pBed4=tab.p21+FluBed.deltaP((bed.h1+c.hBed)./2,bed.eps1,c.rho_p);
    flu.Tbed4=tab.T6;
    flu.rho_g4=DryAir.rho(flu.pBed4,flu.Tbed4);
    flu.w4=air.mDot4./(flu.rho_g4.*c.Afloor1);
    flu.wmf4=FluBed.wmf(c.d_p,c.rho_p,flu.pBed4,flu.Tbed4);
    flu.FG4=flu.w4./flu.wmf4;


    %% Particle flow
    names={'Time','mDot_p','w_p'};
    flow=table('Size',[height(tab),length(names)],...
            'VariableTypes',[{'datetime'},repmat({'double'},1,length(names)-1)]);
    flow.Properties.VariableNames=names;
    flow.Time=tab.Time;
    clear('names');


    %Particle mass flow; fix recording issue
    flow.mDot_p=tab.speed1*3/20;
    
    idx=flow.mDot_p<0.5;
    flow.mDot_p(idx)=round(flow.mDot_p(idx),1);
    flow.mDot_p(~idx)=round(flow.mDot_p(~idx)./0.25).*0.25;


    %Particle velocity = f(eps_mf)
    flow.w_p=flow.mDot_p./(c.rho_p.*(1-c.eps_mf).*c.l.*c.hFlow.*c.psi);


    %% Heat transfer
    names={'Time','Tsurf','Tbed','p','w','wmf','FG','P_el','hVirt','mode'};
    htc=table('Size',[height(tab),length(names)],...
            'VariableTypes',[{'datetime'},repmat({'double'},1,length(names)-2),'logical']);
    htc.Properties.VariableNames=names;
    htc.Time=tab.Time;
    clear('names');


    %Surface and bed temperatures
    htc.Tsurf=tab.T5;
    htc.Tbed=tab.T1;


    %Pressure reference point: at the test tube's height in the second
    %chamber. The test tube is slightly to the left (in the direction of
    %the second chamber) of the baffle separating the second and third
    %chambers (it is fully within the second chamber, right at the border)
    htc.p=tab.p21+tab.p19+FluBed.deltaP(bed.h4+c.hBed-c.d_t2f,...
        bed.eps2,c.rho_p);


    %Flow velocity reference point: at test tube location within chamber 2.
    %Correction of w2 based on simulation results and bed level differences
    s=load('h2FG.mat');
    h2FG=s.h2FG;

    DeltaFG=predict(h2FG,bed.h4-bed.h5);
    htc.FG=flu.FG2+DeltaFG;
    htc.wmf=FluBed.wmf(c.d_p,c.rho_p,htc.p,htc.Tbed);
    htc.w=htc.FG.*htc.wmf;


    %Electric power
    htc.P_el=tab.Pmean1;


    %Virtual HTC (based on plain tube surface)
    htc.hVirt=htc.P_el./(c.A_plain*(htc.Tsurf-htc.Tbed));


    %Control mode: true=Tsurf-Tbed~=100 K, false=constant P_el
    htc.mode=abs(htc.Tsurf-htc.Tbed-100)<2;
    

    %% Set up and fill main table
    names=[{'Time'},htcNames];
    main=table('Size',[height(tab),length(names)],...
            'VariableTypes',[{'datetime'},repmat({'double'},1,length(names)-2),'logical']);
    main.Properties.VariableNames=names;
    main.Time=tab.Time;

    
    main.Tsurf=htc.Tsurf;
    main.Tbed=htc.Tbed;
    main.T1=tab.T1;
    main.T3=tab.T3;
    main.Tco2=tab.Tco2hex1;

    main.p=htc.p;

    main.w=htc.w;
    main.FG=htc.FG;
    main.FG2=flu.FG2;

    main.w_p=flow.w_p;
    main.mDot_p=flow.mDot_p;

    main.P_el=htc.P_el;
    
    main.hVirt=htc.hVirt;

    main.mode=htc.mode;
end




