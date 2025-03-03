%% Calculate Basic Properties
%GNU General Public License v3.0
%By Stefan Thanheiser: https://orcid.org/0000-0003-2765-1156
%
%Part of the paper:
%
%Thanheiser, S.; Haider, M.
%Dispersion Model for Level Control of Bubbling Fluidized Beds with 
%Particle Cross-Flow
%Chemical Engineering Science 2024
%
%All data, along with methodology reports and supplementary documentation, 
%is published in the data repository:
%https://doi.org/10.5281/zenodo.7924693
%
%All required files for this function can be found in the software
%repository:
%https://doi.org/10.5281/zenodo.7948224
%
%
%
%This function calculates the basic properties of the fluidized bed system
%needed for the subsequent analysis of measurements
%
%
%Requires all auxiliary classes and functions on the MATLAB path
%
%Required products, version 24.1:
%   - MATLAB
%Necessary files, classes, functions, and scripts:
%   - @DryAir
%   - @FluBed
%   - @Orifice
%   - @implExp


function main=getProp(tab,c,htcNames,chambers)
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


    %% Particle dispersion
    w_pNames=compose('w_p%d',chambers);
    PhiNames=compose('Phi%d',chambers);
    names=[{'Time',...
            'mDotSand','mDotS'},...
            w_pNames,...
            PhiNames,...
            compose('D%d',2:3)];
    disp=table('Size',[height(tab),length(names)],...
            'VariableTypes',[{'datetime'},repmat({'double'},1,length(names)-1)]);
    disp.Properties.VariableNames=names;
    disp.Time=tab.Time;
    clear('names');


    %Particle mass flows
    disp.mDotSand=tab.speed1*3/20;
    disp.mDotSand=round(disp.mDotSand./0.5).*0.5;     %Fix recording issue
    disp.mDotS=disp.mDotSand./(c.hRef.*c.l);


    %Fictional densities and particle velocities. Both=f(mean(eps))
    disp{:,PhiNames}=c.rho_p.*(bed{:,hNames}+c.hBed)/c.hRef;
    disp{:,w_pNames}=disp.mDotSand./(c.rho_p.*c.l.*(bed{:,hNames}+c.hBed));

    for i=chambers
        j=num2str(i);
        switch i
            case {1,2}
                disp{:,['Phi',j]}=disp{:,['Phi',j]}.*(1-mean(bed.eps1));
                disp{:,['w_p',j]}=disp{:,['w_p',j]}./(1-mean(bed.eps1));
            case {3,4}
                disp{:,['Phi',j]}=disp{:,['Phi',j]}.*(1-mean(bed.eps2));
                disp{:,['w_p',j]}=disp{:,['w_p',j]}./(1-mean(bed.eps2));
            case {5,6}
                disp{:,['Phi',j]}=disp{:,['Phi',j]}.*(1-mean(bed.eps3));
                disp{:,['w_p',j]}=disp{:,['w_p',j]}./(1-mean(bed.eps3));
        end
    end
    
    
    %Particle dispersion coefficients
    disp.D2=disp.mDotS./((disp.Phi5-disp.Phi4)./(c.x2-c.xPress));
    disp.D3=disp.mDotS./((disp.Phi3-disp.Phi2)./(c.x3-c.xPress));


    %% Fluidization
    names=[{'Time'},...
        compose('pBed%d',1:4),...
        compose('Tbed%d',1:4),...
        compose('rho_g%d',1:4),...
        compose('wmf%d',1:4),...
        compose('w%d',1:4),...
        compose('w_e%d',1:4),...
        compose('FG%d',1:4)];
    flu=table('Size',[height(tab),length(names)],'VariableTypes',[{'datetime'},repmat({'double'},1,length(names)-1)]);
    flu.Properties.VariableNames=names;
    flu.Time=tab.Time;
    clear('names');

    
    %Chamber 1: inlet chamber
    flu.pBed1=tab.p21+FluBed.deltaP((bed.h6+c.hBed)./2,bed.eps3,c.rho_p);
    flu.Tbed1=tab.T7;
    flu.wmf1=FluBed.wmf(c.d_p,c.rho_p,flu.pBed1,flu.Tbed1);
    flu.rho_g1=DryAir.rho(flu.pBed1,flu.Tbed1);
    flu.w1=air.mDot1./(flu.rho_g1.*c.Afloor1);
    flu.w_e1=flu.w1-flu.wmf1;
    flu.FG1=flu.w1./flu.wmf1;


    %Chamber 2
    flu.pBed2=tab.p21+tab.p19+FluBed.deltaP((mean([bed.h4,bed.h5],2)+c.hBed)./2,...
            mean([bed.eps2,bed.eps3],2),c.rho_p);
    flu.Tbed2=mean([tab.T7,tab.T1],2);
    flu.wmf2=FluBed.wmf(c.d_p,c.rho_p,flu.pBed2,flu.Tbed2);
    flu.rho_g2=DryAir.rho(flu.pBed2,flu.Tbed2);
    flu.w2=air.mDot2./(flu.rho_g2.*c.Afloor2);
    flu.w_e2=flu.w2-flu.wmf2;
    flu.FG2=flu.w2./flu.wmf2;


    %Chamber 3
    flu.pBed3=tab.p21+tab.p18+FluBed.deltaP((mean([bed.h2,bed.h3],2)+c.hBed)./2,...
            mean([bed.eps1,bed.eps2],2),c.rho_p);
    flu.Tbed3=mean([tab.T6,tab.T1],2);
    flu.wmf3=FluBed.wmf(c.d_p,c.rho_p,flu.pBed3,flu.Tbed3);
    flu.rho_g3=DryAir.rho(flu.pBed3,flu.Tbed3);
    flu.w3=air.mDot3./(flu.rho_g3.*c.Afloor3);
    flu.w_e3=flu.w3-flu.wmf3;
    flu.FG3=flu.w3./flu.wmf3;
    
    
    %Chamber 4: outlet chamber
    flu.pBed4=tab.p21+FluBed.deltaP((bed.h1+c.hBed)./2,bed.eps1,c.rho_p);
    flu.Tbed4=tab.T6;
    flu.wmf4=FluBed.wmf(c.d_p,c.rho_p,flu.pBed4,flu.Tbed4);
    flu.rho_g4=DryAir.rho(flu.pBed4,flu.Tbed4);
    flu.w4=air.mDot4./(flu.rho_g4.*c.Afloor1);
    flu.w_e4=flu.w4-flu.wmf4;
    flu.FG4=flu.w4./flu.wmf4;


    %% Set up and fill main table
    names=[{'Time'},htcNames];
    main=table('Size',[height(tab),length(names)],...
            'VariableTypes',[{'datetime'},repmat({'double'},1,length(names)-1)]);
    main.Properties.VariableNames=names;
    main.Time=tab.Time;



    
    main.Tsurf=tab.T5;
    main.Tbed=tab.T1;
    % main.T1=tab.T1;
    % main.T3=tab.T3;
    % main.T5=tab.T5;
    % main.T6=tab.T6;
    % main.T7=tab.T7;

    main.w_p=mean([disp.w_p3,disp.w_p4],2);


    % %Reference point: center of chamber 2
    % main.p=flu.pBed2;
    % main.w=flu.w2;


    %Reference point: border between chambers 1 and 2. Linear interpolation 
    %between chamber centers
    deltaX=(c.x2+c.x3)./2;
    deltaY=flu.w3-flu.w2;
    x=c.x2./2-c.tube2baffle;
    main.w=deltaY./deltaX.*x+flu.w2;

    % main.w2=flu.w2;
    % main.w3=flu.w3;

    



    %Pressure reference point: at the test tube's height in the second
    %chamber. The test tube is slightly to the left (in the direction of
    %the second chamber) of the baffle separating the second and third
    %chambers (it is fully within the second chamber, right at the border)
    % main.p=flu.pBed2;
    main.p=tab.p21+tab.p19+FluBed.deltaP(bed.h4+c.hBed-c.tube2floor,...
        bed.eps2,c.rho_p);


    
    
    main.hVirt=tab.Pmean1./(c.A_plain*(main.Tsurf-main.Tbed));
    main.P=tab.Pmean1;
    main.mode=abs(main.Tsurf-main.Tbed-100)<2;


    % main.P=tab.Pmean1;
    % main.TCO2=tab.Tco2hex1;
    % main.FG=flu.FG2;
    % main.mDotSand=tab.speed1*3/20;
    % main.hVirt=main.P./(c.A_plain*(main.Tsurf-main.Tbed));


    % main.Tbed=flu.Tbed2;
    % main.rho_g2=flu.rho_g2;
    
    % main.w_e2=flu.w_e2;
    % main.wmf2=flu.wmf2;
    % main.FG2=flu.FG2;
    
    % main.p0=tab.p21;
    
    % main.Tleft=tab.T7;
    % main.tab.T1=tab.T1;
    % main.Tright=tab.T6;
    
    % main.epsLeft=bed.eps3;
    % main.epsCenter=bed.eps2;
    % main.epsRight=bed.eps1;

    % main.mDotSand=disp.mDotSand;
    % main.mDotS=disp.mDotS;
    
    % main.AC1=tab.AO1;
    % main.AC2=tab.AO2;
    % main.AC1set=tab.ACset1+c.hBed;
    % main.AC2set=tab.ACset2+c.hBed;

    % main.AC1(main.AC1>0.99)=1;  %Fix recording issue
    % main.AC2(main.AC2>0.99)=1;  %Fix recording issue
    
    % main.air1=air.mDot1;
    % main.air2=air.mDot2;
    % main.air3=air.mDot3;
    % main.air4=air.mDot4;

    % main{:,hNames}=bed{:,hNames};
    % main{:,w_pNames}=disp{:,w_pNames};
    % main{:,PhiNames}=disp{:,PhiNames};
    % main.D2=disp.D2;


    %Set h, w_p, Phi and D2 values to NaN where h is negative
    % for i=chambers
    %     j=num2str(i);
    %     pneg=tab{:,['p',j+3]}<0;
    % 
    %     main{pneg,['h',j]}=NaN;
    %     main{pneg,['w_p',j]}=NaN;
    %     main{pneg,['Phi',j]}=NaN;
    % end
    % 
    % main.D2(isnan(main.Phi4) | isnan(main.Phi5))=NaN;


end




