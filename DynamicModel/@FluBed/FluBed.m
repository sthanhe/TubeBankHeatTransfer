%% Functions for Fluidized Beds
%GNU General Public License v3.0
%By Stefan Thanheiser: https://orcid.org/0000-0003-2765-1156
%
%Part of the paper:
%
%Thanheiser, S.; Haider, M.
%Dispersion Model for Level Control of Bubbling Fluidized Beds with 
%Particle Cross-Flow
%Chemical Engineering Research and Design 2025
%
%All required files for this class can be found in the software
%repository:
%https://doi.org/10.5281/zenodo.7948224
%
%
%
%This class contains several functions for the calculation of properties in
%a fluidized bed. References are given in the respective functions.
%
%
%Requires all files packaged in the class folder and on the MATLAB path
%
%Required products, version 24.1:
%   - MATLAB
%Additional classes:
%   - DryAir
%   - implExp


classdef FluBed
    %All parameters and results in SI base units

    properties(Constant)
        g=9.81; %Gravitational acceleration, m/s²
    end
    
    
    methods(Static)
        function [wmf,Re]=wmf(d_p,rho_p,p,T)
            %Minimum fluidization velocity
            %Approximation to the Ergun equation by assuming / estimating
            %specific particle sphericity (phi_s) and porosity at minimum
            %fluidization (eps_mf)
            
            persistent C1 C2
            if isempty(C1)
                %C1, C2 according to Richardson (1971):
                %Kunii, D.; Levenspiel, O. Heat Transfer between Fluidized 
                %Beds and Surfaces. In Fluidization Engineering, 2nd ed.; 
                %Kunii, D.,Levenspiel, O., Eds.; Butterworth-Heinemann: 
                %Boston, MA, USA, 1991; p. 70.
                
                %Values correspond to phi_s=1, eps_mf=0.4
                %or roughly to phi_s=0.8, eps_mf=0.45
                C1=25.7;
                C2=0.0365;

                %Code to recalculate phi_s, eps_mf from other constants
%                 C1=25.7;
%                 C2=0.0365;
%                 K1=1/C2;
%                 K2=C1*2*K1;
%                 phi=@(eps_mf) 1.75/(eps_mf^3*K1);
%                 eps_mf=fzero(@(eps) K2-150*(1-eps)/(eps^3*phi(eps)^2),[0.1,0.9]);
%                 phi_s=phi(eps_mf);
            end

            Ar=FluBed.Ar(d_p,rho_p,p,T);
            Re=sqrt(C1.^2+C2.*Ar)-C1;

            wmf=Re.*DryAir.eta(T)./(d_p.*DryAir.rho(p,T));
        end
        
        
        function [wmf,Re_mf,Ar]=wmfErgun(d_p,rho_p,phi_s,eps_mf,p,T)
            %Minimum fluidization velocity according to the Ergun equation
            
            sz=implExp.size(d_p,rho_p,phi_s,eps_mf,p,T);
            [d_p,rho_p,phi_s,eps_mf,p,T]=implExp.normalize(sz,d_p,rho_p,phi_s,eps_mf,p,T);

            
            K1=1.75./(phi_s.*eps_mf.^3);
            K2=150.*(1-eps_mf)./(phi_s.^2.*eps_mf.^3);            

            Ar=FluBed.Ar(d_p,rho_p,p,T);

            C1=K2./(2*K1);
            C2=1./K1;
            Re_mf=-C1+sqrt(C1.^2+C2.*Ar);

            wmf=Re_mf.*DryAir.eta(T)./(d_p.*DryAir.rho(p,T));
            wmf=reshape(wmf,sz);
        end
        
        
        function eps=eps(deltaP,deltaH,rho_p)
            %Bed porosity when both pressure taps at a vertical distance of
            %deltaH are in the fluidized bed
            eps=1-deltaP./(rho_p.*FluBed.g.*deltaH);
        end
        
        
        function [eps,d_b]=porosity(w0,eps_mf,pitch,d_p,rho_p,p,T,d_H,z)
            %Floor area based on hydraulic diameter d_H=4*l*w/(l+w)
            d_H(d_H>1.2)=1.2;   %limited to 1.2 m according to Grace, p. 141
            A=d_H.^2*pi/4;  
            
            
            %Bubble diameter. Assumption: bubble diameter is equal to the 
            %maximum bubble diameter within the tube bank, but limited to
            %the horizontal pitch
            wmf=FluBed.wmf(d_p,rho_p,p,T);
            d_b0=3.685/FluBed.g*(w0-wmf).^2;  %Initial bubble diameter, K/L p. 131, based on Miwa (1972), but Choi (1998) gives a factor of 3.685 instead of 2.78 while referring to the same source
%             d_bm=0.65/100*(A*10^4.*(w0-wmf)*100).^0.4;   %Maximum bubble diameter, K/L p.146
            d_bm=2.59*(A.*(w0-wmf)./FluBed.g.^0.5).^0.4;    %Maximum bubble diameter, Grace p. 141
            d_b=d_bm-(d_bm-d_b0).*exp(-0.3*z./d_H);
            d_b=min([d_b;repmat(pitch,1,length(d_b))],[],1);
            
            
            %Bubble velocity, calculated from the rise velocity of a single
            %bubble and approximating bubble interaction and coalescence
            %with w0-wmf, according to Grace, p. 142
            w_br=0.711*sqrt(FluBed.g.*d_b);     %Single bubble rise velocity
            w_b=w0-wmf+w_br;
            
            
            %Fraction of bubble phase in the bed (Holdup), K/L p. 156-157
            if any(w_b<wmf./eps_mf,'all')
                warning(['w_b<wmf/eps_emf --> Bubbles may be slower than emulsion gas. \n' ...
                            'Fraction of the bed in bubbles (delta, "Holdup") may be overestimated. \n' ...
                            'See Kunii / Levenspiel, p. 156']);
            end
            c=(w_b.*eps_mf./wmf-1).*wmf/4;  %Smoothing factor
            c=min([c;repmat(2*wmf,1,length(c))]);
            delta=(w0-wmf)./(w_b+wmf-c);
            
            
            %Porosity
            eps_b=1;        %Porosity in bubble phase. Assumption that bubble phase is particle-free
            eps_e=eps_mf;   %Porosity in emulsion phase. Assumption that porosity in emulsion is about equal to porosity at minimum fluidization
            eps=delta.*eps_b+(1-delta).*eps_e;
        end
        
        
        function h=h(deltaP,eps,rho_p)
            %Bed level when the top pressure tap is not in the fluidized
            %bed
            h=deltaP./(rho_p.*FluBed.g.*(1-eps));
        end
        
        
        function deltaP=deltaP(deltaH,eps,rho_p)
            %Pressure drop across a fluidized bed of height deltaH
            deltaP=rho_p.*FluBed.g.*deltaH.*(1-eps);
        end


        function [phi_s,C1,C2]=phi_s(wmf,d_p,rho_p,eps_mf,p,T)
            %Calculation of effective sphericity phi_s based on the
            %measurements of minimum fluidization velocity wmf and
            %porosity eps_mf. All values at minimum fluidization conditions

            phi_s=fzero(@(phi) wmfErgun(d_p,rho_p,phi,eps_mf,p,T)-wmf,[0.1,0.9]);

            K1=1.75./(phi_s.*eps_mf.^3);
            K2=150.*(1-eps_mf)./(phi_s.^2.*eps_mf.^3);
            C2=1/K1;
            C1=K2/(2*K1);
        end


        function D=D(w,w_p,d_p,rho_p,p,T)
            %Particle dispersion coefficient
            persistent c eps2 eps3 epsAr
            if isempty(c)
                c=19317.8340011694;
                eps2=1.10169319091282;
                eps3=-2.04940397390448;
                epsAr=0.108571357573225;
            end

            den=sqrt(d_p.*FluBed.g);

            pi2=(w-FluBed.wmf(d_p,rho_p,p,T))./den;
            pi3=w_p./den;
            Ar=FluBed.Ar(d_p,rho_p,p,T);


            idx=pi2>0;
            D=zeros(size(pi2));
            if any(idx)
                D(idx)=c.*d_p.*den.*...
                        pi2(idx).^eps2.*...
                        (1+pi3(idx)).^eps3.*...
                        Ar(idx).^epsAr;
            end
        end


        function Ar=Ar(d_p,rho_p,p,T_A)
            %Archimedes number, assuming dry air as fluidizing gas
            rho_g=DryAir.rho(p,T_A);

            Ar=rho_g.*d_p.^3.*(rho_p-rho_g).*FluBed.g./DryAir.eta(T_A).^2;
        end
    end
    
    
    
    methods(Static, Access=protected)
        function Re=Re(d_p,w,p,T)
            %Reynolds number with respect to particle diameter
            Re=d_p.*w./(DryAir.ny(p,T));
        end
    end
end




