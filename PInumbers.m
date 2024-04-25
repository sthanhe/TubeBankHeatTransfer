n=11;
dims=4;

A=table('Size',[n,dims],'VariableTypes',repmat({'double'},1,dims),...
        'VariableNames',{'kg','m','s','K'},...
        'RowNames',{'g','d_p','rho_p','c_p','my_g','k_g','ny','mDotS','h','d_T','eps'});

%Gravitational acceleration, m/s²
A.m('g')=1;       
A.s('g')=-2;

%Particle diameter, m
A.m('d_p')=1;

%Particle density, kg/m³
A.kg('rho_p')=1;
A.m('rho_p')=-3;

%Particle specific heat capcity, J/kgK=m²/s²K
A.m('c_p')=2;
A.s('c_p')=-2;
A.K('c_p')=-1;

%Dynamic viscosity gas, Pas=kg/ms
A.kg('my_g')=1;
A.m('my_g')=-1;
A.s('my_g')=-1; 

%Thermal conductivity gas, W/mK=kgm/s³K
A.kg('k_g')=1;
A.m('k_g')=1;
A.s('k_g')=-3;
A.K('k_g')=-1;     

%Specific mass flow, kg/m²s
A.kg('mDotS')=1;
A.m('mDotS')=-2;
A.s('mDotS')=-1;

%Density gradient, kg/m^4
% A.kg('densGrad')=1;
% A.m('densGrad')=-4;

%Mass diffusivity, m²/s
A.m('ny')=2;
A.s('ny')=-1;    

%Heat transfer coefficient, W/m²K=kg/s³K
A.kg('h')=1;
A.s('h')=-3;
A.K('h')=-1;

%Tube diameter, m
A.m('d_T')=1;    

%Porosity (1-eps_mf), -
% eps=[0,0,0,0];



K=table('Size',[n-dims,n],'VariableTypes',repmat({'double'},1,n),...
        'VariableNames',A.Properties.RowNames);

%d_p^3*g*rho_p^2/my_g^2
K.d_p(1)=3;
K.g(1)=1;
K.rho_p(1)=2;
K.my_g(1)=-2;

%h*d_T/k_g=Nu
K.h(2)=1;
K.d_T(2)=1;
K.k_g(2)=-1;

%mDotS/rho_p*tau/d_T
% k.ny(3)=1;
% k.densGrad(3)=1;
K.mDotS(3)=1;
K.c_p(3)=1;
K.d_p(3)=2;
K.k_g(3)=-1;
K.d_T(3)=-1;
K.eps(3)=-1;

%rho_p*c_p*ny/k_g
K.rho_p(4)=1;
K.c_p(4)=1;
K.ny(4)=1;
K.k_g(4)=-1;

%mDotS*c_p*d_T/k_g
% k.ny(5)=1;
% k.densGrad(5)=1;
K.mDotS(5)=1;
K.c_p(5)=1;
K.d_T(5)=1;
K.k_g(5)=-1;

%rho_p*ny/my_g
K.rho_p(6)=1;
K.ny(6)=1;
K.my_g(6)=-1;

%eps
K.eps(7)=1;


%Check
if all(K{:,:}*A{:,:}==0,'all')
    disp('Variable system is linearly independent');
else
    warning('Variable system is linearly dependent!')
end







