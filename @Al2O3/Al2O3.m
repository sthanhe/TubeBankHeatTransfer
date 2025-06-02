%% Property functions of aluminum oxide (Al2O3)
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
% All required files for this class can be found in the software
% repository: see the link to the supplemental release in the data 
% repository here: https://doi.org/10.5281/zenodo.15576311
%
%
%
% This class describes the thermo-physical properties of aluminum oxide
% (corundum; only alpha-phase!) according to:
%
% NIST chemistry webbook
% https://webbook.nist.gov/cgi/cbook.cgi?ID=C1344281&Mask=2
%
%
%Requires all files packaged in the class folder and on the MATLAB path
%
%Required products, version 24.1:
%   - MATLAB


classdef Al2O3
    %Alpha-phase only! (corundum)
    %All parameters and results in SI base units

    
    %% Constants
    properties(Constant)
        M=101.9613e-3;     %Molar mass
        rho=3940;          %Density
    end
    
    properties(Constant, Access=private)
        A=102.4290;
        B=38.74980;
        C=-15.91090;
        D=2.628181;
        E=-3.007551;
        F=-1717.930;
        G=146.9970;
        H=-1675.690;
    end
    
    
    %% Property Functions
    methods(Static)
        function c_p=c_p(T)
            %Specific isobaric heat capacity
            T=T./1000;
            c_p=(Al2O3.A+Al2O3.B*T+Al2O3.C*T.^2+...
                Al2O3.D*T.^3+Al2O3.E./T.^2)...
                ./Al2O3.M;
        end
        
        
        function h=h(T)
            %Specific enthalpy
            %h(298.15)=0
            T=T./1000;
            h=(Al2O3.A*T+Al2O3.B*T.^2./2+Al2O3.C*T.^3./3+...
                Al2O3.D*T.^4./4-Al2O3.E./T+Al2O3.F-Al2O3.H)...
                ./Al2O3.M.*1000;
        end
        
        
        function s=s(T)
            %Specific entropy
            T=T./1000;
            s=(Al2O3.A*log(T)+Al2O3.B*T+Al2O3.C*T.^2./2+...
                Al2O3.D*T.^3./3-Al2O3.E./(2*T.^2)+Al2O3.G)...
                ./Al2O3.M;
        end


        function lambda=lambda(T)
            %Thermal conductivity (source: CRC-Handbook, 2012)
            p1=6.4919;
            p2=-14.2215;
            q1=-0.2313;
            q2=-1.9612;

            std=787.4;
            m=1173;

            lambda=NaN(size(T));

            idx=T<373.15;
            lambda(idx)=30;

            Tnorm=(T-m)./std;
            lambda(~idx)=(p1.*Tnorm(~idx)+p2)./(Tnorm(~idx).^2+q1.*Tnorm(~idx)+q2);
        end
    end
end




