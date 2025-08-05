%% Relative Nusselt number
% GNU General Public License v3.0
% By Stefan Thanheiser: https://orcid.org/0000-0003-2765-1156
%
% Part of the paper:
%
% Thanheiser, S.
% Molerus and Wirth's Heat Transfer Model for Bubbling Fluidized Beds: 
% Proposal for an Extended Model Including Immersed Tube Banks and Particle 
% Cross-Flow
%
% All data, along with methodology reports and supplementary documentation, 
% is published in the data repository:
% https://doi.org/10.5281/zenodo.15576311
%
% All required files for this function can be found in the software
% repository: 
% https://doi.org/10.5281/zenodo.15576950
%
%
%
% This function calculates the maximum Nusselt number depending on a certain
% parameter relative to the maximum Nusselt number at a specific parameter
% value. It is used by the scripts "calcPi8" and "calcPi9" to analyse the
% relative influence of tube diameter (pi8) and tube packing density (pi9)
% on the maximum Nusselt number.
%
%
%Requires all auxiliary classes and functions on the MATLAB path
%
%Required products, version 24.1:
%   - MATLAB
%   - Curve Fitting Toolbox
%Necessary classes, functions, files, and scripts:
%   - @DryAir
%   - @FluBed
%   - @SiO2
%   - @implExp


function Nu_relExt=Nu_rel(fx,beta,p,T,Ar,eps_mf,phi_s,d_t,w_p,c_p,p_h)
    % Inputs:
    % fx: model function in the form fx(b,X), where b are the regression
    %       coefficients, function handle
    % beta: regression coefficients, double
    % p: bed pressure, double
    % T: bed temperature, double
    % Ar: Archimedes number, double
    % eps_mf: bed voidage at minimum fluidization conditions, double
    % phi_s: particle sphericity, double
    % d_t: tube diameter, double
    % w_p: mean horizontal particle velocity, double
    % c_p: specific heat capacity of particles, double
    % p_h: horizontal pitch, double
    % 
    % 
    % Outputs:
    % Nu_relExt: maximum Nusselt number relative to the maximum Nusselt
    % number at a specific parameter


    %Particle diameter derived from Archimedes number
    rho_p=SiO2.rho(T);      %Particle density
    rho_g=DryAir.rho(p,T);  %Gas density
    
    d_p=(rho_g.*(rho_p-rho_g).*FluBed.g./DryAir.eta(T).^2./Ar).^(-1/3);


    %Fluidization velocities
    w_mf=FluBed.wmfErgun(d_p,rho_p,phi_s,eps_mf,p,T);   %Minimum    
    
    w=arrayfun(@(w_mf) ...
        linspace(w_mf,20*w_mf,1000)',...
        w_mf,'UniformOutput',false);
    w=horzcat(w{:});


    %Record size and create third dimension
    sz3=max(numel(d_t),numel(p_h));
    if isscalar(d_t)
        p_h=reshape(p_h,[1,1,numel(p_h)]);
    else
        d_t=reshape(d_t,[1,1,numel(d_t)]);
    end


    %Implicit expansion
    sz=implExp.size(w,T,p,d_p,rho_p,phi_s,eps_mf,...
        c_p,d_t,p_h,w_p);

    [w,T,p,d_p,rho_p,phi_s,eps_mf,...
        c_p,d_t,p_h,w_p]=implExp.normalize(sz,w,T,p,d_p,rho_p,phi_s,...
            eps_mf,c_p,d_t,p_h,w_p);


    %Gas-convective Nusselt number: from Molerus and Wirth
    [~,Nu]=FluBed.molWirth(w,T,p,d_p,rho_p,phi_s,eps_mf,@SiO2.c_p);
    Nu_gc=Nu.gc';


    %Particle-convective Nusselt number: from model
    pis=FluBed.piFactors(w,T,p,d_p,rho_p,phi_s,eps_mf,...
                c_p,d_t,p_h,w_p);

    Nu_pc=fx(beta,pis);

    
    %Total Nusselt number
    Nu=Nu_pc+Nu_gc;
    Nu=reshape(Nu,sz);
    

    %Maximum Nusselt number
    Nu_max=max(Nu,[],1,'omitmissing');
    Nu_max=reshape(Nu_max,[length(Ar),sz3]);
    

    %Relative maximum Nusselt number
    Nu_relExt=Nu_max./Nu_max(:,end);
end




