function Nu_relExt=Nu_rel(fx,beta,p,T,Ar,eps_mf,phi_s,d_t,w_p,c_p,p_h)
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
    [~,Nu]=FluBed.molerus(w,T,p,d_p,rho_p,phi_s,eps_mf,@SiO2.c_p);
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