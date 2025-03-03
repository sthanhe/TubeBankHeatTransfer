%% Functions for Fluidized Beds
%GNU General Public License v3.0
%By Stefan Thanheiser: https://orcid.org/0000-0003-2765-1156
%
%Part of the paper:
%
%Thanheiser, S.; Haider, M.
%Dispersion Model for Level Control of Bubbling Fluidized Beds with 
%Particle Cross-Flow
%Applied Thermal Energy 2024
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


        function phi_s=eps2phi(eps_mf)
            %Calculates sphericity phi_s from the bed porosity at minimum
            %fluidization conditions eps_mf
            %Based on: 
            %Wen, C.Y. and Yu, Y.H. (1966), A generalized method for 
            %predicting the minimum fluidization velocity. 
            %AIChE J., 12: 610-612. https://doi.org/10.1002/aic.690120343

            phi_s=sqrt((1-eps_mf)./(11*eps_mf.^3));

        end


        function D=D(w,w_p,d_p,rho_p,p,T)
            %Particle dispersion coefficient
            persistent c eps2 eps3 epsAr
            if isempty(c)
                c=19030.1124724027;
                eps2=1.08175825313772;
                eps3=-3.62266055421733;
                epsAr=0.109738530498836;
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


        function [h,Nu]=molerus(w,T,p,d_p,rho_p,phi_s,eps_mf,c_pfx)
            persistent G1 G2 P1 P2 P1ast P2ast
            if isempty(G1)
                G1=0.05;
                G2=0.165;
                P1=25;
                P2=0.19;
                P1ast=33.3;
                P2ast=0.125;
                % P1ast=P1;
                % P2ast=P2;
            end


            %Implicit expansion
            sz=implExp.size(w,T,p,d_p,rho_p,phi_s,eps_mf);
            [w,T,p,d_p,rho_p,phi_s,eps_mf]=implExp.normalize(sz,w,T,p,d_p,rho_p,phi_s,eps_mf);

            
            %Gas and particle properties
            k_g=DryAir.lambda(T);     %Surface temperature?
            my=DryAir.eta(T);
            c_p=c_pfx(T);
            rho_g=DryAir.rho(p,T);
            rho_e=rho_p-rho_g;      %Excess density


            %Fluidization velocities
            % w_mf=FluBed.wmf(d_p,rho_p,p,T);     %Minimum
            w_mf=FluBed.wmfErgun(d_p,rho_p,phi_s,eps_mf,p,T);   %Minimum
            w_e=w-w_mf;                                         %Excess
            w_e(w_e<0)=NaN;                                     %Avoid complex results


            %Pi-factors, p. 70
            %pi1=Nu;
            pi2=k_g./(2*c_p.*my);
            pi3=DryAir.Pr(T);     %Surface temperature?
            pi4=rho_g./rho_e;
            pi5=(rho_p.*c_p./(k_g.*FluBed.g)).^(1/3).*w_e;
            pi6=w_e./w_mf;
            pi7=1-eps_mf;


            %Flow length scales
            l_t=(my./sqrt(FluBed.g*rho_e.*rho_g)).^(2/3);   %Eq. 4.25, p. 45
            turb2lam=pi4.^(1/3);                            %Eq. 7.16, p. 64
            l_l=l_t.*turb2lam;


            %Maximum particle convective HTC (laminar flow regime), Section 4.3
            Nu_pcMax=P2*pi7./(1+pi2);   %0.09=0.19*(1-eps_mf), Eq. 7.13, p. 62 
            h_pcMax=Nu_pcMax.*k_g./l_l;


            %Maximum gas convective HTC (turbulent flow regime), Section 4.4
            Nu_gcMax=G2*pi3.^(1/3).*turb2lam;
            h_gcMax=Nu_gcMax.*k_g./l_l;


            %Gas convection, Section 7.3
            gfx=(1+G1*pi6.^-1).^-1;     %g(w_e) damping function
            Nu_gc=Nu_gcMax.*gfx;        %Nusselt number for gas convection
            h_gc=Nu_gc.*k_g./l_l;       %HTC gas convection


            %Pure particle convection, Section 7.5
            pfx=@(P1) (1+P1*(pi6.^(1/3).*pi5).^-1).^-1;     %p(w_e) damping function
            Nu_pcPure=Nu_pcMax.*pfx(P1);                    %Nusselt number for pure particle convection
            h_pcPure=Nu_pcPure.*k_g./l_l;                   %HTC pure particle convection
            

            %Mixed particle convection, Section 7.6
            dfx=0.28*pi2.*pi7.^2.*sqrt(pi4).*pi5.^2.*pi6.^-1;   %d(w_e) damping function
            Nu_pcMix=P2ast*pi7./(1+pi2+dfx).*pfx(P1ast);        %Nusselt number for mixed particle convection
            h_pcMix=Nu_pcMix.*k_g./l_l;                         %HTC mixed particle convection


            %Mixed gas and particle convection, Section 7.6
            Nu_mix=Nu_pcMix+Nu_gc;
            h_mix=Nu_mix.*k_g./l_l;


            %Flow regimes, Section 4
            Ar=FluBed.Ar(d_p,rho_p,p,T);
            lam=Ar<=1e2;
            mix=1e2<Ar & Ar<1e5;
            turb=1e5<=Ar & Ar<=1e8;


            %Assign Nusselt numbers to their regimes
            Nu_total=NaN(size(Ar));
            Nu_total(lam)=Nu_pcPure(lam);
            Nu_total(mix)=Nu_mix(mix);
            Nu_total(turb)=Nu_gc(turb);
            h_total=Nu_total.*k_g./l_l;


            %Build output structures
            r=@(x) reshape(x,sz);
            Nu=struct('total',r(Nu_total),...
                    'gcMax',r(Nu_gcMax),'pcMax',r(Nu_pcMax),...
                    'gc',r(Nu_gc),'pcPure',r(Nu_pcPure),'pcMix',r(Nu_pcMix),...
                    'mix',r(Nu_mix));
            h=struct('total',r(h_total),...
                    'gcMax',r(h_gcMax),'pcMax',r(h_pcMax),...
                    'gc',r(h_gc),'pcPure',r(h_pcPure),'pcMix',r(h_pcMix),...
                    'mix',r(h_mix));
        end


        function [h,Nu]=molExt(w,T,p,d_p,rho_p,phi_s,eps_mf,c_pfx,...
                d_t,p_h,w_p)
            %Molerus' heat transfer correlation, extended (mixed regime
            %only)
            persistent G1 G2 P1 P2 P3 P4 %P5
            if isempty(G1)
                G1=0.165;
                G2=0.05;

                %s3 old
                % P1=0.0669009973829433;
                % P2=0.331251286752616;
                % P3=19.4345474590870;
                % 
                % P4=4.49720169564995e-05;
                % P5=-0.959145413661024;


                %s3
                % P1=0.0696085079075269;
                % P2=19.8830853275972;
                % 
                % P3=5.37302742339688e-05;
                % P4=1.01505510898379;

                %s2
                P1=0.0690771105844349;
                P2=18.9085028208424;

                P3=6.45824254376922e-05;
                P4=1.15226744565364;
            end
            

            %Gas and particle properties
            k_g=DryAir.lambda(T);
            my_g=DryAir.eta(T);
            rho_g=DryAir.rho(p,T);
            rho_e=rho_p-rho_g;
            l_l=(my_g./(rho_e.*sqrt(FluBed.g))).^(2/3);
            c_p=c_pfx(T);
            
            
            %Fluidization velocities
            w_mf=FluBed.wmfErgun(d_p,rho_p,phi_s,eps_mf,p,T);
            w_e=w-w_mf;
            w_e(w_e<0)=NaN;
            
            
            %Pi-factors            
            % pi1=h.*l_l./k_g;
            pi2=k_g./(2*c_p.*my_g);
            pi3=DryAir.Pr(T);
            pi4=rho_g./rho_e;
            pi5=(rho_e.*c_p./(k_g.*FluBed.g)).^(1/3).*w_e;
            pi6=(rho_e.*c_p./(k_g.*FluBed.g)).^(1/3).*w_mf;
            pi7=1-eps_mf;
            pi8=d_t./l_l;
            pi9=d_t./p_h;
            pi10=(rho_e.*c_p./(k_g.*FluBed.g)).^(1/3).*w_p;


            %Gas convection
            Nu_gcMax=G1*(pi3.*pi4).^(1/3);

            d_gc=(1+G2*pi6./pi5).^-1;
            Nu_gc=Nu_gcMax.*d_gc;


            %Particle convection
            %s3 old
            % t=1+P2.*pi7.^2.*sqrt(pi4).*pi5.*pi6;
            % s=tanh(P4.*pi8);
            % Nu_pcMax=P1.*pi7./(1+pi2.*t.*s);

            %s3
            % t=1+0.28.*pi7.^2.*sqrt(pi4).*pi5.*pi6;
            % s=tanh(P3.*pi8);
            % Nu_pcMax=P1.*pi7./(1+pi2.*t.*s);

            %s2
            t=1+0.28.*pi7.^2.*sqrt(pi4).*pi5.*pi6;
            s=1-exp(-P3.*pi8);
            Nu_pcMax=P1.*pi7./(1+pi2.*t.*s);


            %s3 old
            % pfx=(1-pi9).^P5;
            % d_pc=(1+P3.*(pi6./pi5).^(1/3)./pi5.*pfx).^-1;
            % Nu_pc=Nu_pcMax.*d_pc;
    
            %s3 and s2
            pfx=(1-pi9).^P4;
            d_pc=(1+P2.*(pi6./pi5).^(1/3)./pi5./pfx).^-1;
            Nu_pc=Nu_pcMax.*d_pc;


            %Heat transfer regimes: only use mixed
            Ar=FluBed.Ar(d_p,rho_p,p,T);
            mix=1e2<=Ar & Ar<=1e5;


            %Mixed gas and particle convection
            Nu=NaN(size(Ar));
            Nu(mix)=Nu_gc(mix)+Nu_pc(mix);
            h=Nu.*k_g./l_l;
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




