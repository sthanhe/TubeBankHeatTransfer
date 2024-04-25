function h=molerus(u,T_b,T_s,p,d_p,phi_s,eps_mf)
    persistent g
    if isempty(g)
        g=9.81;
    end
    
    sz=implExp.size(u,T_b,T_s,p,d_p,phi_s,eps_mf);
    [u,T_b,T_s,p,d_p,phi_s,eps_mf]=implExp.normalize(sz,u,T_b,T_s,p,d_p,phi_s,eps_mf);
    
    
    %gas and particle properties
    k_g=DryAir.lambda(T_s);
    Pr=DryAir.Pr(T_s);
    my=DryAir.eta(T_b);
    rho_g=DryAir.rho(p,T_b);
    c_p=SiO2.c_p(T_b);
    rho_p=SiO2.rho;
    
    
    u_mf=wmfErgun(d_p,rho_p,phi_s,eps_mf,p,T_b);
    
    u_exc=u-u_mf;
    u_excND=(rho_p.*c_p./(k_g.*g)).^(1/3).*u_exc;
    u_excDL=u_exc./u_mf;
    
    densRel=rho_g./(rho_p-rho_g);
    cond2conv=k_g./(2*c_p.*my);
    partDens=1-eps_mf;
    l_l=(my./(sqrt(g)*(rho_p-rho_g))).^(2/3);
    
    
    %particle convection
    x=(1+33.3*(u_excDL.^(1/3).*u_excND).^-1).^-1;
    d=0.28*cond2conv.*partDens.^2.*sqrt(densRel).*u_excND.^2.*u_excDL.^-1;
    Nu_pc=0.125*partDens.*x./(1+cond2conv.*d);
    
    
    %gas convection
    Nu_gc=0.165*Pr.^(1/3).*densRel.^(1/3).*(1+0.05*u_excDL.^-1).^-1;
    
    
    h=(Nu_pc+Nu_gc).*k_g./l_l;
    h=reshape(h,sz);
    
end