%% Property functions of dry air at ambient pressure
%GNU General Public License v3.0
%By Stefan Thanheiser: https://orcid.org/0000-0003-2765-1156
%
%Modified from:
%Stefan Thanheiser, "sthanhe/HeatTransfer: Round 3 Release”. Zenodo, Jan. 
%22, 2022. doi: 10.5281/zenodo.5911319.
%
%
%Part of the paper:
%
%Thanheiser, S.; Haider, M.
%Particle Mass Diffusion Model for Level Control of Bubbling Fluidized Beds
%with Horizontal Particle Flow
%Powder Technology 2023
%
%All required files for this class can be found in the software
%repository:
%https://doi.org/10.5281/zenodo.xxxxxxx
%
%
%
%This class describes the thermo-physical properties of dry air as an ideal
%gas according to:
%
%Span, R. Properties of Dry Air. In VDI Heat Atlas, 2nd ed.; Stephan, P., 
%Kabelac, S., et al., Eds.; Springer: Berlin Heidelberg, Germany, 2010; 
%pp. 172–191. https://doi.org/10.1007/978-3-540-77877-6_11
%
%
%Requires all files packaged in the class folder and on the MATLAB path
%
%Required products:
%   - MATLAB, version 9.14
%Data files:
%   - dryAir.xls
%   - dryAirTable.mat


classdef DryAir 
    %All parameters and results in SI base units
    
    %% Constants
    properties(Constant)
        M=28.9583e-3;   %molar mass
        R=287.12;       %specific gas constant
    end
    
    
    %% prop(T) functions
    methods(Static)
        function rho=rho(p,T)
            %Density
            rho=p./(DryAir.R.*T);
        end
        
        
        function h=h(T)
            %Specific enthalpy
            persistent tab
            if isempty(tab)
                tabStruct=coder.load('@DryAir\dryAirTable.mat','tab');
                tab=tabStruct.tab;
            end
            
            h=interp1(tab.T,tab.h,T);
        end
        
        
        function s=s(T)
            %Specific entropy. Does not account for pressure variations!
            persistent tab
            if isempty(tab)
                tabStruct=coder.load('@DryAir\dryAirTable.mat','tab');
                tab=tabStruct.tab;
            end
            
            s=interp1(tab.T,tab.s,T);
        end
        
        
        function c_p=c_p(T)
            %Specific isobaric heat capacity
            persistent tab
            if isempty(tab)
                tabStruct=coder.load('@DryAir\dryAirTable.mat','tab');
                tab=tabStruct.tab;
            end
            
            c_p=interp1(tab.T,tab.c_p,T);
        end
        
        
        function c_v=c_v(T)
            %Specific isochoric heat capacity
            persistent tab
            if isempty(tab)
                tabStruct=coder.load('@DryAir\dryAirTable.mat','tab');
                tab=tabStruct.tab;
            end
            
            c_v=interp1(tab.T,tab.c_v,T);
        end
        
        
        function kappa=kappa(T)
            %Isentropic exponent
            kappa=DryAir.c_p(T)./DryAir.c_v(T);
        end
        
        
        function beta=beta(T)
            %Coefficient of thermal expansion
            persistent tab
            if isempty(tab)
                tabStruct=coder.load('@DryAir\dryAirTable.mat','tab');
                tab=tabStruct.tab;
            end
            
            beta=interp1(tab.T,tab.beta,T);
        end
        
        
        function w_s=w_s(T)
            %Speed of sound
            persistent tab
            if isempty(tab)
                tabStruct=coder.load('@DryAir\dryAirTable.mat','tab');
                tab=tabStruct.tab;
            end
            
            w_s=interp1(tab.T,tab.w_s,T);
        end
        
        
        function lambda=lambda(T)
            %Heat conductivity
            persistent tab
            if isempty(tab)
                tabStruct=coder.load('@DryAir\dryAirTable.mat','tab');
                tab=tabStruct.tab;
            end
            
            lambda=interp1(tab.T,tab.lambda,T);
        end
        
        
        function eta=eta(T) %#codegen
            %Dynamic viscosity
            persistent tab
            if isempty(tab)
                tabStruct=coder.load('@DryAir\dryAirTable.mat','tab');
                tab=tabStruct.tab;
            end
            
            eta=interp1(tab.T,tab.eta,T);
        end
        
        
        function ny=ny(p,T)
            %Kinematic viscosity
            persistent tab
            if isempty(tab)
                tabStruct=coder.load('@DryAir\dryAirTable.mat','tab');
                tab=tabStruct.tab;
            end
            
            ny=DryAir.eta(T)./DryAir.rho(p,T);
        end
        
        
        function a=a(p,T)
            %Thermal diffusivity            
            a=DryAir.ny(p,T)./DryAir.Pr(T);
        end
        
        
        function Pr=Pr(T)
            %Prandtl number
            persistent tab
            if isempty(tab)
                tabStruct=coder.load('@DryAir\dryAirTable.mat','tab');
                tab=tabStruct.tab;
            end
            
            Pr=interp1(tab.T,tab.Pr,T);
        end
    end
    
    
    %% Other property functions
    methods(Static)
        function T=T_h(h)
            %Backwards-equation for temperature as function of specific
            %enthalpy
            persistent tab
            if isempty(tab)
                tabStruct=coder.load('@DryAir\dryAirTable.mat','tab');
                tab=tabStruct.tab;
            end
            
            T=interp1(tab.h,tab.T,h);
        end
    end
    
    
    %% Auxilliary Functions
    methods(Static)
        function createConstants()
            %Creates the lookup-constants from the Excel-file
            clear('DryAir');
            
            tab=readtable('@DryAir\dryAir.xls','Range','A1:M53');

            tab{:,'T'}=tab{:,'T'}+273.15;
            tab{:,'h'}=tab{:,'h'}*1000;
            tab{:,'s'}=tab{:,'s'}*1000;
            tab{:,'c_p'}=tab{:,'c_p'}*1000;
            tab{:,'c_v'}=tab{:,'c_v'}*1000;
            tab{:,'beta'}=tab{:,'beta'}/1000;
            tab{:,'lambda'}=tab{:,'lambda'}/1000;
            tab{:,'eta'}=tab{:,'eta'}*10^-6;
            tab{:,'ny'}=tab{:,'ny'}*10^-7;
            tab{:,'a'}=tab{:,'a'}*10^-7;

            save('dryAirTable.mat','tab');
        end
    end
end




