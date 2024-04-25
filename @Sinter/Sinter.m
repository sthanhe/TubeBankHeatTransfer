%% Sinter Plate Calculation 
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
%This class contains functions to calculate the flow of dry air through a 
%sinter plate
%
%
%Requires all files packaged in the class folder and on the MATLAB path
%
%Required products:
%   - MATLAB, version 9.14
%Data files:
%   - constants.xls
%   - constants.mat
%Additional classes:
%   - DryAir


classdef Sinter
    methods(Static)
        function deltaP=deltaP(w,p,T,s,name)
            persistent tab
            if isempty(tab)
                tabStruct=coder.load('@Sinter\constants.mat','tab');
                tab=tabStruct.tab;
            end
            
            idx=matches(tab.Name,name);
            alpha=tab.alpha(idx);
            beta=tab.beta(idx);

            eta=DryAir.eta(T);
            rho=DryAir.rho(p,T);

            deltaP=w.*s.*(eta./alpha+rho.*w./beta);
        end


        function [w2,w1,mDotS]=w(p1,p2,T,s,name)
            persistent tab
            if isempty(tab)
                tabStruct=coder.load('@Sinter\constants.mat','tab');
                tab=tabStruct.tab;
            end

            % sz=implExp.size(p1,p2,T,s,name);
            % [p1,p2,T,s,name]=implExp.normalize(sz,p1,p2,T,s,name);

            p1=reshape(p1,[1,numel(p1)]);
            p2=reshape(p2,[1,numel(p2)]);
            T=reshape(T,[1,numel(T)]);
            
            idx=matches(tab.Name,name);
            alpha=tab.alpha(idx);
            beta=tab.beta(idx);


            eta=DryAir.eta(T);
            rho=DryAir.rho(mean([p1;p2],1),T);

            p=eta.*beta./(rho.*alpha);
            q=-beta.*abs(p1-p2)./(rho.*s);
            w=-p./2+sqrt((p./2).^2-q);
            w=sign(p1-p2).*w;

            mDotS=w.*rho;
            w1=mDotS./DryAir.rho(p1,T);
            w2=mDotS./DryAir.rho(p2,T);

            % mDotS=reshape(mDotS,sz);    %in kg/(m²s)
            % w1=reshape(w1,sz);
            % w2=reshape(w2,sz);
        end


        function loadConstants()
            tab=readtable('@Sinter/constants.xls','Range','A:F');
            tab.alpha=tab.alpha*10^-12;
            tab.beta=tab.beta*10^-7;
            tab.x=tab.x*10^-6;
            tab.DeltaRho=tab.DeltaRho*10^2;
            tab.tau=tab.tau*10^-6;

            save('@Sinter/constants.mat','tab');
        end
    end
end




