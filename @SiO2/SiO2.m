%% Property functions of silicon dioxide (SiO2)
%GNU General Public License v3.0
%By Stefan Thanheiser: https://orcid.org/0000-0003-2765-1156
%
%
%This class describes the thermo-physical properties of silicon dioxide
%(SiO2) according to:
%
%Stefan Thanheiser, dissertation at TU Wien (2023)
%
%
%Requires all files packaged in the class folder and on the MATLAB path
%
%Required products:
%   - MATLAB, version 9.10
%   - Curve Fitting Toolbox, version 3.5.13
%Necessary files, classes and functions:
%   - createFits.m
%   - fits.mat


classdef SiO2
    %All parameters and results in SI base units
    
    %%
    properties(Constant)
        M=60.0843e-3;       %Molar mass

        TtransQuartz=847;   %Alpha-beta transition temperature of quartz
        TtransTri=390;      %Alpha-beta transition temperature of tridymite
        TtransCri=543;      %Alpha-beta transition temperature of cristobalite

        TmeltQuartz=1696;   %Melting temperature of quartz
        TmeltTri=1743;      %Melting temperature of tridymite
        TmeltCri=1996;      %Melting temperature of cristobalite
    end
    
    properties(Constant, Access=private)
        A=-6.076591;
        B=251.6755;
        C=-324.7964;
        D=168.5604;
        E=0.002548;
        F=-917.6893;

        A_beta=58.75340;
        B_beta=10.27925;
        C_beta=-0.131384;
        D_beta=0.025210;
        E_beta=0.025601;
        F_beta=-929.3292;

        TminQuartz=100;     %Minimum temperature of quartz functions
        TminTri=298;        %Minimum temperature of tridymite functions
        TminCri=0;          %Minimum temperature of cristobalite functions
    end
    
    
    %% Property Functions
    methods(Static)
        function c_p=c_p(T,phase)
            %Specific isobaric heat capacity
            persistent Atr Btr Ctr Dtr Etr Ftr Acr Bcr Ccr Dcr Ecr Fcr
            if isempty(Atr)
                Atr=3.27*4.184/SiO2.M;
                Btr=24.8e-3*4.184/SiO2.M;
                Ctr=74.904/SiO2.M;
                Dtr=3.0999e-3/SiO2.M;
                Etr=-2.3669e2/SiO2.M;
                Ftr=-1.174e6/SiO2.M;

                Acr=-0.0025;
                Bcr=3.2550;
                Ccr=-9.5255;
                Dcr=-8.5069e+06;
                Ecr=-1.6403;
                Fcr=1.2666e+03;
            end


            if nargin<2
                phase='quartz';
            end

            c_p=NaN(size(T));
            switch phase
                case 'quartz'
                    [alpha,beta]=SiO2.phasesQuartz(T);
                    T=T./1000;
                    
                    c_p(alpha)=c_pfx(T(alpha),SiO2.A,SiO2.B,SiO2.C,SiO2.D,SiO2.E);
                    c_p(beta)=c_pfx(T(beta),SiO2.A_beta,SiO2.B_beta,SiO2.C_beta,SiO2.D_beta,SiO2.E_beta);

                    % A=-122.698;
                    % B=592.567e-3;
                    % C=-29.755e5;
                    % D=-367.971e-6;
                    % E=-558.489e-9;
                    % F=590.415e-12;

                    % Ta=T(alpha);
                    % c_p(alpha)=A+B*Ta-C*Ta.^-2+D*Ta.^2+E*Ta.^3+F*Ta.^4;
                    % c_p=c_p./SiO2.M;
                
                case 'tridymite'
                    [alpha,beta]=SiO2.phasesTridymite(T);
                    Tbeta=T(beta);

                    c_p(alpha)=Atr+Btr*T(alpha);
                    c_p(beta)=Ctr+Dtr*Tbeta+Etr*Tbeta.^-0.5+Ftr*Tbeta.^-2;
                    
                case 'cristobalite'
                    [alpha,beta]=SiO2.phasesCristobalite(T);
                    Talpha=T(alpha);
                    
                    c_p(alpha)=Acr*Talpha.^2+Bcr*Talpha+Ccr;
                    c_p(beta)=Dcr*T(beta).^Ecr+Fcr;
            end
            
            
            function c_p=c_pfx(T,A,B,C,D,E)
                c_p=(A+B*T+C*T.^2+D*T.^3+E./T.^2)./SiO2.M;
            end
        end
        
        
        function h=h(T,phase)
            %Specific enthalpy
            %h(298.15)=0
            persistent H alphaTri betaTri alphaCri betaCri
            if isempty(H)
                H=-910.8568;

                vars=load('@SiO2\fits.mat','h_TalphaTri','h_TbetaTri','h_TalphaCri','h_TbetaCri');
                alphaTri=vars.h_TalphaTri;
                betaTri=vars.h_TbetaTri;
                alphaCri=vars.h_TalphaCri;
                betaCri=vars.h_TbetaCri;
            end


            if nargin<2
                phase='quartz';
            end

            h=NaN(size(T));
            switch phase
                case 'quartz'
                    [alpha,beta]=SiO2.phasesQuartz(T);
                    T=T./1000;
                    
                    h=NaN(size(T));
                    h(alpha)=hfx(T(alpha),SiO2.A,SiO2.B,SiO2.C,SiO2.D,SiO2.E,SiO2.F);
                    h(beta)=hfx(T(beta),SiO2.A_beta,SiO2.B_beta,SiO2.C_beta,SiO2.D_beta,SiO2.E_beta,SiO2.F_beta);
                
                case 'tridymite'
                    [alpha,beta]=SiO2.phasesTridymite(T);

                    h(alpha)=alphaTri(T(alpha));
                    h(beta)=betaTri(T(beta));
                    
                case 'cristobalite'
                    [alpha,beta]=SiO2.phasesCristobalite(T);
                    
                    h(alpha)=alphaCri(T(alpha));
                    h(beta)=betaCri(T(beta));
            end
            
            
            function h=hfx(T,A,B,C,D,E,F)
                h=(A*T+B*T.^2./2+C*T.^3./3+D*T.^4./4-E./T+F-H)./SiO2.M.*1000;
            end
        end
        
        
        function s=s(T)
            %Specific entropy
            persistent G G_beta
            if isempty(G)
                G=-27.96962;
                G_beta=105.8092;
            end
            
            [alpha,beta]=SiO2.phasesQuartz(T);
            T=T./1000;
            
            s=NaN(size(T));
            s(alpha)=sfx(T(alpha),SiO2.A,SiO2.B,SiO2.C,SiO2.D,SiO2.E,G);
            s(beta)=sfx(T(beta),SiO2.A_beta,SiO2.B_beta,SiO2.C_beta,SiO2.D_beta,SiO2.E_beta,G_beta);
            
            
            function s=sfx(T,A,B,C,D,E,G)
                s=(A*log(T)+B*T+C*T.^2./2+D*T.^3./3-E./(2*T.^2)+G)./SiO2.M;
            end
        end


        function rho=rho(T,phase)
            if nargin<2
                phase='quartz';
            end

            rho=NaN(size(T));
            switch phase
                case 'quartz'
                    [alpha,beta]=SiO2.phasesQuartz(T);
                    rho(alpha)=2648;
                    rho(beta)=2533;

                case 'tridymite'
                    [alpha,beta]=SiO2.phasesTridymite(T);
                    rho(alpha)=2265;
                    rho(beta)=2185;
                    
                case 'cristobalite'
                    [alpha,beta]=SiO2.phasesCristobalite(T);
                    rho(alpha)=2334;
                    rho(beta)=2448;
            end
        end
    end
    
    
    %% Backward Equations
    methods(Static)
        function T=T_h(h,phase)
            %Backwards-equation for temperature as function of specific
            %enthalpy
            persistent alphaQuartz betaQuartz alphaTri betaTri alphaCri betaCri hQuartzAlphaMin hQuartzAlphaMax hQuartzBetaMin hQuartzBetaMax hTriAlphaMin hTriAlphaMax hTriBetaMin hTriBetaMax hCriAlphaMin hCriAlphaMax hCriBetaMin hCriBetaMax
            if isempty(alphaQuartz)
                %Property functions
                vars=load('@SiO2\fits.mat','T_halphaQuartz','T_hbetaQuartz','T_halphaTri','T_hbetaTri','T_halphaCri','T_hbetaCri');
                alphaQuartz=vars.T_halphaQuartz;
                betaQuartz=vars.T_hbetaQuartz;
                alphaTri=vars.T_halphaTri;
                betaTri=vars.T_hbetaTri;
                alphaCri=vars.T_halphaCri;
                betaCri=vars.T_hbetaCri;


                %Boundaries
                hQuartzAlphaMin=SiO2.h(SiO2.TminQuartz,'quartz');
                hQuartzAlphaMax=SiO2.h(SiO2.TtransQuartz-eps(SiO2.TtransQuartz),'quartz');
                hQuartzBetaMin=SiO2.h(SiO2.TtransQuartz,'quartz');
                hQuartzBetaMax=SiO2.h(SiO2.TmeltQuartz,'quartz');

                hTriAlphaMin=SiO2.h(SiO2.TminTri,'tridymite');
                hTriAlphaMax=SiO2.h(SiO2.TtransTri-eps(SiO2.TtransTri),'tridymite');
                hTriBetaMin=SiO2.h(SiO2.TtransTri,'tridymite');
                hTriBetaMax=SiO2.h(SiO2.TmeltTri,'tridymite');

                hCriAlphaMin=SiO2.h(SiO2.TminCri,'cristobalite');
                hCriAlphaMax=SiO2.h(SiO2.TtransCri-eps(SiO2.TtransCri),'cristobalite');
                hCriBetaMin=SiO2.h(SiO2.TtransCri,'cristobalite');
                hCriBetaMax=SiO2.h(SiO2.TmeltCri,'cristobalite');
            end


            if nargin<2
                phase='quartz';
            end

            T=NaN(size(h));
            switch phase
                case 'quartz'
                    alpha=hQuartzAlphaMin<=h & h<=hQuartzAlphaMax;
                    trans=hQuartzAlphaMax<h & h<hQuartzBetaMin;
                    beta=hQuartzBetaMin<=h & h<=hQuartzBetaMax;

                    T(alpha)=alphaQuartz(h(alpha));
                    T(trans)=SiO2.TtransQuartz;
                    T(beta)=betaQuartz(h(beta));
                
                case 'tridymite'
                    alpha=hTriAlphaMin<=h & h<=hTriAlphaMax;
                    trans=hTriAlphaMax<h & h<hTriBetaMin;
                    beta=hTriBetaMin<=h & h<=hTriBetaMax;

                    T(alpha)=alphaTri(h(alpha));
                    T(trans)=SiO2.TtransTri;
                    T(beta)=betaTri(h(beta));
                    
                case 'cristobalite'
                    alpha=hCriAlphaMin<=h & h<=hCriAlphaMax;
                    trans=hCriAlphaMax<h & h<hCriBetaMin;
                    beta=hCriBetaMin<=h & h<=hCriBetaMax;

                    T(alpha)=alphaCri(h(alpha));
                    T(trans)=SiO2.TtransCri;
                    T(beta)=betaCri(h(beta));
            end
        end
        
        
        function createConstants()
            %Creates the curve fittings for lookup
            clear('SiO2');
            T0=298.15;  %Zero reference point for enthalpy
            
            TalphaQuartz=linspace(SiO2.TminQuartz,SiO2.TtransQuartz-eps(SiO2.TtransQuartz),1e4);
            TbetaQuartz=linspace(SiO2.TtransQuartz,SiO2.TmeltQuartz,1e4);
            TalphaTri=linspace(SiO2.TminTri,SiO2.TtransTri-eps(SiO2.TtransTri),1e4);
            TbetaTri=linspace(SiO2.TtransTri,SiO2.TmeltTri,1e4);
            TalphaCri=linspace(SiO2.TminCri,SiO2.TtransCri-eps(SiO2.TtransCri),1e4);
            TbetaCri=linspace(SiO2.TtransCri,SiO2.TmeltCri,1e4);

            halphaQuartz=SiO2.h(TalphaQuartz);
            hbetaQuartz=SiO2.h(TbetaQuartz);
            halphaTri=arrayfun(@(T) integral(@(T) SiO2.c_p(T,'tridymite'),TalphaTri(1),T),TalphaTri);
            hbetaTri=arrayfun(@(T) integral(@(T) SiO2.c_p(T,'tridymite'),TbetaTri(1),T),TbetaTri);
            halphaCri=arrayfun(@(T) integral(@(T) SiO2.c_p(T,'cristobalite'),TalphaCri(1),T),TalphaCri);
            hbetaCri=arrayfun(@(T) integral(@(T) SiO2.c_p(T,'cristobalite'),TbetaCri(1),T),TbetaCri);
            
            halphaTri0=interp1(TalphaTri,halphaTri,T0);
            halphaCri0=interp1(TalphaCri,halphaCri,T0);

            halphaTri=halphaTri-halphaTri0;
            halphaCri=halphaCri-halphaCri0;

            hbetaTri=hbetaTri+halphaTri(end)+2785;
            hbetaCri=hbetaCri+halphaCri(end)+22353;


            fitresult=SiO2.createFits(halphaQuartz,TalphaQuartz,...
                                hbetaQuartz,TbetaQuartz,...
                                halphaTri,TalphaTri,...
                                hbetaTri,TbetaTri,...
                                halphaCri,TalphaCri,...
                                hbetaCri,TbetaCri);

            
            T_halphaQuartz=fitresult{1};
            T_hbetaQuartz=fitresult{2};
            T_halphaTri=fitresult{3};
            T_hbetaTri=fitresult{4};
            T_halphaCri=fitresult{5};
            T_hbetaCri=fitresult{6};
            h_TalphaTri=fitresult{7};
            h_TbetaTri=fitresult{8};
            h_TalphaCri=fitresult{9};
            h_TbetaCri=fitresult{10};
            
            save('@SiO2\fits.mat','T_halphaQuartz','T_hbetaQuartz',...
                                    'T_halphaTri','T_hbetaTri',...
                                    'T_halphaCri','T_hbetaCri',...
                                    'h_TalphaTri','h_TbetaTri',...
                                    'h_TalphaCri','h_TbetaCri');
        end
    end
    
    
    %% Internal Functions
    methods(Static, Access=private)
        [fitresult,gof]=createFits(halphaQuartz,TalphaQuartz,...
                                    hbetaQuartz,TbetaQuartz,...
                                    halphaTri,TalphaTri,...
                                    hbetaTri,TbetaTri,...
                                    halphaCri,TalphaCri,...
                                    hbetaCri,TbetaCri)
        
        
        function [alpha,beta]=phasesQuartz(T)
            alpha=SiO2.TminQuartz<=T & T<SiO2.TtransQuartz;
            beta=SiO2.TtransQuartz<=T & T<=SiO2.TmeltQuartz;
        end


        function [alpha,beta]=phasesTridymite(T)
            alpha=SiO2.TminTri<=T & T<SiO2.TtransTri;
            beta=SiO2.TtransTri<=T & T<=SiO2.TmeltTri;
        end


        function [alpha,beta]=phasesCristobalite(T)
            alpha=SiO2.TminCri<=T & T<SiO2.TtransCri;
            beta=SiO2.TtransCri<=T & T<=SiO2.TmeltCri;
        end
    end
end




