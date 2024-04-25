%% Calculate Basic Properties
%GNU General Public License v3.0
%By Stefan Thanheiser: https://orcid.org/0000-0003-2765-1156
%
%Part of the paper:
%
%Thanheiser, S.; Haider, M.
%Particle Mass Diffusion Model for Level Control of Bubbling Fluidized Beds
%with Horizontal Particle Flow
%Powder Technology 2023
%
%All data, along with methodology reports and supplementary documentation, 
%is published in the data repository:
%https://doi.org/10.5281/zenodo.7924694
%
%All required files for this function can be found in the software
%repository:
%https://doi.org/10.5281/zenodo.xxxxxxx
%
%
%
%This function calculates the basic properties of the fluidized bed system
%needed for the subsequent analysis of measurements
%
%
%Requires all auxiliary classes and functions on the MATLAB path
%
%Required products:
%   - MATLAB, version 9.14
%Necessary files, classes, functions, and scripts:
%   - @DryAir
%   - @FluBed
%   - @Orifice
%   - @implExp
%   - getProp.m


function main=getProp(tab,c,htcNames)
    %% Air flows
    names={'Time','mDot','On',...
        'O1','O2','O3','O4','O5','O6','O7','O8',...
        'O1scal','O2scal','O3scal','O4scal','O5scal','O6scal','O7scal','O8scal',...
        'mDot1','mDot2','mDot3','mDot4'};
    air=table('Size',[height(tab),length(names)],'VariableTypes',[{'datetime','double','logical'},repmat({'double'},1,length(names)-3)]);
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
    S=air.mDot./(air.O1+air.O2+air.O3+air.O4+air.O5+air.O6+air.O7+air.O8);
    air.O1scal=air.O1.*S;
    air.O2scal=air.O2.*S;
    air.O3scal=air.O3.*S;
    air.O4scal=air.O4.*S;
    air.O5scal=air.O5.*S;
    air.O6scal=air.O6.*S;
    air.O7scal=air.O7.*S;
    air.O8scal=air.O8.*S;
    clear('S');
    
    
    %Aggregated mass flows
    air.mDot1=air.O8scal;
    air.mDot2=air.O5scal+air.O6scal+air.O7scal;
    air.mDot3=air.O2scal+air.O3scal+air.O4scal;
    air.mDot4=air.O1scal;
    
    
    
    %% Bed properties
    names={'Time','eps1','eps2','eps3','h1','h2','h3','h4','h5','h6'};
    bed=table('Size',[height(tab),length(names)],'VariableTypes',[{'datetime'},repmat({'double'},1,length(names)-1)]);
    bed.Properties.VariableNames=names;
    bed.Time=tab.Time;
    clear('names');
    
    
    %Bed porosities: only calculate values when air supply is on
    bed.eps1(air.On)=FluBed.eps(tab.p1(air.On),c.dh_eps1,c.rho_p);
    % bed.eps2(air.On)=FluBed.eps(tab.p2(air.On),c.dh_eps2,c.rho_p);
    bed.eps3(air.On)=FluBed.eps(tab.p3(air.On),c.dh_eps1,c.rho_p);
    
    %Assume minimum fluidization porosities when air supply is not on
    bed.eps1(~air.On)=c.epsMf;
    % bed.eps2(~air.On)=c.epsMf;
    bed.eps3(~air.On)=c.epsMf;

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
    names={'Time','pBed1','pBed2','pBed3','pBed4',...
        'Tbed1','Tbed2','Tbed3','Tbed4',...
        'wmf1','wmf2','wmf3','wmf4',...
        'mDotMf1','mDotMf2','mDotMf3','mDotMf4',...
        'FG1','FG2','FG3','FG4'};
    flu=table('Size',[height(tab),length(names)],'VariableTypes',[{'datetime'},repmat({'double'},1,length(names)-1)]);
    flu.Properties.VariableNames=names;
    flu.Time=tab.Time;
    clear('names');

    Tcenter=mean([tab.T1,tab.T3,tab.T5],2);

    
    %Chamber 1: inlet chamber
    flu.pBed1=tab.p21+FluBed.deltaP((bed.h6+c.hBed)./2,bed.eps3,c.rho_p);
    flu.Tbed1=tab.T7;
    flu.wmf1=FluBed.wmf(c.d_p,c.rho_p,flu.pBed1,flu.Tbed1);
    flu.mDotMf1=flu.wmf1.*c.Afloor1.*DryAir.rho(flu.pBed1,flu.Tbed1);
    flu.FG1=air.mDot1./flu.mDotMf1;


    %Chamber 2
    flu.pBed2=tab.p21+tab.p19+FluBed.deltaP((mean([bed.h4,bed.h5],2)+c.hBed)./2,...
            mean([bed.eps2,bed.eps3],2),c.rho_p);
    flu.Tbed2=mean([tab.T7,Tcenter],2);
    flu.wmf2=FluBed.wmf(c.d_p,c.rho_p,flu.pBed2,flu.Tbed2);
    flu.mDotMf2=flu.wmf2.*c.Afloor2.*DryAir.rho(flu.pBed2,flu.Tbed2);
    flu.FG2=air.mDot2./flu.mDotMf2;


    %Chamber 3
    flu.pBed3=tab.p21+tab.p18+FluBed.deltaP((mean([bed.h2,bed.h3],2)+c.hBed)./2,...
            mean([bed.eps1,bed.eps2],2),c.rho_p);
    flu.Tbed3=mean([tab.T6,Tcenter],2);
    flu.wmf3=FluBed.wmf(c.d_p,c.rho_p,flu.pBed3,flu.Tbed3);
    flu.mDotMf3=flu.wmf3.*c.Afloor3.*DryAir.rho(flu.pBed3,flu.Tbed3);
    flu.FG3=air.mDot3./flu.mDotMf3;
    
    
    %Chamber 4: outlet chamber
    flu.pBed4=tab.p21+FluBed.deltaP((bed.h1+c.hBed)./2,bed.eps1,c.rho_p);
    flu.Tbed4=tab.T6;
    flu.wmf4=FluBed.wmf(c.d_p,c.rho_p,flu.pBed4,flu.Tbed4);
    flu.mDotMf4=flu.wmf4.*c.Afloor1.*DryAir.rho(flu.pBed4,flu.Tbed4);
    flu.FG4=air.mDot4./flu.mDotMf4;


    %% Main table
    names=[{'Time'},htcNames];
    main=table('Size',[height(tab),length(names)],'VariableTypes',[{'datetime'},repmat({'double'},1,length(names)-1)]);
    main.Properties.VariableNames=names;
    main.Time=tab.Time;


    %Fill table
    main.P=tab.Pmean1;
    main.Tsurf=tab.T5;
    main.Tsand=tab.T1;
    main.TCO2=tab.Tco2hex1;
    main.FG=flu.FG2;
    main.mDotSand=tab.speed1*3/20;
    
    main.hVirt=main.P./(c.AtestTube*(main.Tsurf-main.Tsand));
    main.mode=abs(main.Tsurf-main.Tsand-100)<2;


    % main.FG2=flu.FG2;
    % main.wmf2=flu.wmf2;
    % main.Tbed2=flu.Tbed2;
    % 
    % main.mDotSand=tab.speed1*3/20;
    % main.mDotSand=round(main.mDotSand./0.5).*0.5;     %Fix recording issue
    % main.mDotS=main.mDotSand./(c.hRef.*c.l);
    % 
    % main.p0=tab.p21;
    % 
    % main.Tleft=tab.T7;
    % main.Tcenter=Tcenter;
    % main.Tright=tab.T6;
    % 
    % main.epsLeft=bed.eps3;
    % main.epsCenter=bed.eps2;
    % main.epsRight=bed.eps1;
    % 
    % main.AC1=tab.AO1;
    % main.AC2=tab.AO2;
    % 
    % main.air1=air.mDot1;
    % main.air2=air.mDot2;
    % main.air3=air.mDot3;
    % main.air4=air.mDot4;
    % 
    % for j=chambers
    %     k=num2str(j);
    % 
    %     main{:,['h',k]}=bed{:,['h',k]};
    %     main{tab{:,['p',k+3]}<0,['h',k]}=NaN;
    % end


end




