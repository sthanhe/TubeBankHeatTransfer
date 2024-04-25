%% Orifice Plate Calculation 
%GNU General Public License v3.0
%By Stefan Thanheiser: https://orcid.org/0000-0003-2765-1156%
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
%This class calculates the mass flow of dry air through an orifice plate
%based on EN ISO 51677 (method qm)
%
%
%Requires all files packaged in the class folder and on the MATLAB path
%
%Required products:
%   - MATLAB, version 9.14
%Additional classes:
%   - DryAir


classdef Orifice
    %All parameters and results in SI base units

    methods(Static)
        function [qm,deltaOmega]=qm(p1,p2,T1,d,D,tap)  
            %qm: Mass flow through orifice plate
            %deltaOmega: Total pressure loss through orifice plate
            %p1: Absolute pressure upstream
            %p2: Absolute pressure downstream
            %T1: Temperature upstream
            %d: inner diameter orifice plate
            %D: inner diamater pipe
            %tap: type of pressure tap; 'corner', 'D-D/2' or 'flange'


            my1=DryAir.eta(T1);
            rho1=DryAir.rho(p1,T1);
            kappa=DryAir.kappa(T1);
            
            beta=d./D;
            deltaP=p1-p2;
            epsilon=epsOrifice;
            
            A1=epsilon.*d.^2.*sqrt(2.*deltaP.*rho1)./(my1.*D.*sqrt(1-beta.^4));
            C=0.6;
            err=1;
            while err>1e-3
                Re=C.*A1;
                C=C_orifice;
                err=abs((A1-Re./C)./A1);
            end
        
            qm=pi/4*my1.*D.*Re;
            
            
            x=sqrt(1-beta.^4*(1-C.^2));
            y=C.*beta.^2;
            deltaOmega=(x-y)./(x+y).*deltaP;
        
        
        
            function C=C_orifice()
                persistent a
                if isempty(a)        
                    a=NaN(1,16);
                    a(1)=0.5961;
                    a(2)=0.0261;
                    a(3)=0.216;
                    a(4)=0.000521;
                    a(5)=0.0188;
                    a(6)=0.0063;
                    a(7)=0.043;
                    a(8)=0.08;
                    a(9)=0.123;
                    a(10)=0.11;
                    a(11)=0.031;
                    a(12)=0.8;
                    a(13)=0.011;
                    a(14)=0.75;
                    a(15)=2.8;
                    a(16)=25.4;
                end
                
                
                switch tap
                    case 'corner'
                        L1=0;
                        L2=0;
                    case 'D-D/2'
                        L1=1;
                        L2=0.47;
                    case 'flange'
                        L1=25.4e-3./D;
                        L2=L1;
                end
                
                
                beta=d./D;
                M2=(2.*L2)./(1-beta);
                A=((19e3.*beta)./Re).^0.8;
                
                C=a(1)+a(2).*beta.^2-a(3).*beta.^8+a(4).*(10^6.*beta./Re).^0.7+ ...
                    (a(5)+a(6).*A).*beta.^3.5.*(10^6./Re).^0.3+ ...
                    (a(7)+a(8).*exp(-10.*L1)-a(9).*exp(-7.*L1)).*(1-a(10).*A).*beta.^4./(1-beta.^4)- ...
                    a(11).*(M2-a(12).*M2.^1.1).*beta.^1.3;
                
                if D<71.12e-3
                    C=C+a(13).*(a(14)-beta).*(a(15)-D./(1000*a(16)));
                end    
            end
            
            
            function epsilon=epsOrifice
                persistent a
                if isempty(a)
                    a=NaN(1,3);
                    a(1)=0.351;
                    a(2)=0.256;
                    a(3)=0.93;
                end
                
                epsilon=1-(a(1)+a(2).*beta.^4+a(3).*beta.^8)*(1-(p2./p1).^(kappa.^-1));
            end

        end

    end
end

