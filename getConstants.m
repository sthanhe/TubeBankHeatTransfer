%% Load Common Constants
%GNU General Public License v3.0
%By Stefan Thanheiser: https://orcid.org/0000-0003-2765-1156
%
%Part of the paper:
%
%Thanheiser, S.; Haider, M.
%Particle Mass Diffusion Model for Level Control of Bubbling Fluidized Beds
%with Horizontal Particle Flow
%Powder Technology 2023
%
%All data, along with methodology reports and supplementary documentation, 
%is published in the data repository:
%https://doi.org/10.5281/zenodo.7924694
%
%All required files for this function can be found in the software
%repository:
%https://doi.org/10.5281/zenodo.xxxxxxx
%
%
%
%This function loads the common constants needed for other calculations.
%
%
%Requires all auxiliary classes and functions on the MATLAB path
%
%Required products:
%   - MATLAB, version 9.14
%Necessary files, classes, functions, and scripts:
%   - @DryAir


function c=getConstants()
    c=struct();

    c.g=9.81;   %Gravitational acceleration
    
    c.d_p=174.494e-6;   %Particle diameter (GRANUSIL, Wedron IL #801, Grade 7020, to be confirmed)
    c.rho_p=2650;       %Particle density
    c.eps_mf=0.45;      %Bed porosity at minimum fluidization conditions
    c.phi_s=0.8;        %Particle sphericity
    
    c.dh_eps1=50e-3;    %Height diffc.xPresserence for porosity measurement, p1 and p3
    c.dh_eps2=100e-3;   %Height difference for porosity measurement, p2
    
    c.rho_N=DryAir.rho(1014e2,21+273.15);     %Reference density for flowmeter at 1014 hPa and 21°C
    c.Apipe=154.08e-3^2*pi/4;                 %Inner pipe area of main air line (6 inch, schedule 40 pipe)
    c.OnLimit=50e-3;                          %Minimum mass flow at which air supply is assumed to be on
    
    c.hBed=461.5e-3+14e-3;  %Persistent bed height, includes 14 mm from outlet weir to pressure probe that may not be filled with sand
    c.hRef=1;               %Reference bed height

    c.tube2floor=299.5e-3;  %Distance from the test tube's axis to the distributor floor
    c.tube2baffle=20e-3;    %Distance from the test tube's axis to the baffle between 2. and 3. chamber (in the direction of the 2. chamber)
    
    c.dOrif=14.75e-3;     %Orifice plate inner diameter
    c.DOrif=82.5e-3;      %Orifice plate outer diameter (=inner pipe diameter)
    c.tap='D-D/2';        %Orifice plate pressure tap type
    
    c.x1=202e-3;        %Length of inlet / outlet chamber (first and fourth chamber)
    c.x2=1.06;          %Length of second chamber
    c.x3=0.8;           %Length of third chamber
    c.xPress=37e-3;     %Distance between pressure taps for bed level measurement (next to baffles)
    c.l=0.5;            %Width of all chambers
    
    c.Afloor1=c.x1*c.l;	    %Distributor floor area, inlet / outlet chamber (first and fourth chamber)
    c.Afloor2=c.x2*c.l;	    %Distributor floor area, second chamber
    c.Afloor3=c.x3*c.l;     %Distributor floor area, third chamber

    c.p_h=40e-3;    %Tube bundle horizontal pitch

    c.d_t=20e-3;    %Outside tube diameter (without fins)
    c.h_f=10e-3;    %Fin height
    c.s_f=2e-3;     %Fin thickness
    c.pitch=9e-3;   %Fin pitch
    l_lead=5e-3;    %Lead distance from the end of the tube to the start of the fins

    t=(c.l-2*l_lead)./c.pitch;      %Number of turns
    k=c.pitch./(c.d_t*pi);          %Fin slope
    x=@(r) 2*r*pi*t.*sqrt(1+k.^2);  %Helix length depending on the reference radius r
    
    c.A_bottom=c.s_f.*x(c.d_t/2);           %Plain tube area where the fin is attached
    c.A_sides=2*c.h_f*x(c.d_t/2+c.h_f/2);   %Main fin area on the sides
    c.A_plain=c.d_t*pi.*c.l;                %Plain tube area


    %% Effective horizontal pitch
    %Projection into base plane as a function of axis variable z
    helfx=@(z) (c.d_t/2+c.h_f)*sin(2*pi*z./c.pitch);

    %Start of z where projection is not covered by plain tube anymore
    zStart=fzero(@(z) helfx(z)-c.d_t./2,[0,c.pitch./4]);

    %Projected fin and plain tube areas
    A1=integral(helfx,zStart,c.pitch./4);
    A2=c.s_f.*c.h_f;
    Atube=c.d_t.*c.l;

    Afins=2.*(2.*A1+A2).*t;     %Total projected fin area
    Aproj=Afins+Atube;          %Total projected area
    A0=c.p_h.*c.l;              %Free area

    %Effective pitch
    s_h=Aproj./A0;          %Relative pitch=relative blocked area
    c.p_hEff=c.d_t./s_h;    %Effective horizontal pitch


    % fig=figure(1);
    % clf(fig);
    % ax=gca();
    % colors=ax.ColorOrder;
    % hold(ax,'on');
    % 
    % legItems=cell(1,3);
    % 
    % n=1000;
    % z1=linspace(0,c.pitch,n);
    % legItems{1}=area(z1,helfx(z1),'FaceColor',colors(1,:));
    % area(z1+c.s_f,helfx(z1),'FaceColor',colors(1,:));
    % 
    % z2=linspace(c.pitch/4,c.pitch/4+c.s_f,n);
    % y2=repmat(c.d_t/2+c.h_f,1,n);
    % legItems{2}=area(z2,y2,'FaceColor',colors(2,:));
    % area(z2+c.pitch/2,-y2,'FaceColor',colors(2,:));
    % 
    % z3=linspace(0,c.pitch+c.s_f,n);
    % y3=repmat(c.d_t/2,1,n);
    % legItems{3}=area(z3,y3,'FaceColor',colors(3,:));
    % area(z3,-y3,'FaceColor',colors(3,:));
    % 
    % hold(ax,'off');
    % 
    % ax.XLim=[0,c.pitch+c.s_f];
    % 
    % legItems=vertcat(legItems{:});
    % legend(ax,legItems,{'A_1','A_2','A_{tube}'},'Location','best');
    % 
    % xlabel(ax,'z (m)');
    % ylabel(ax,'y (m)');
    % 
    % fig.Units='centimeters';
    % fig.Position=[10,5,17,8.5];


    
end




