%% Load common constants
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
% All data, along with methodology reports and supplementary documentation, 
% is published in the data repository:
% https://doi.org/10.5281/zenodo.15576311
%
% All required files for this script can be found in the software
% repository: see the link to the supplemental release in the data 
% repository
%
%
%
% This function loads the common constants needed for other calculations. 
% It also calculates the effective pitch and plots the corresponding graph
% published in the Methodology Report.
%
%
%Requires all auxiliary classes and functions on the MATLAB path
%
%Required products, version 24.1:
%   - MATLAB
%   - Curve Fitting Toolbox
%Necessary files, classes, functions, and scripts:
%   - @DryAir
%   - @SiO2


function c=getConstants(plotPitch)
    % Inputs:
    % plotPitch: whether or not to create the graph illustrating the
    % calculations of the effective pitch, logical. Optional, default=false
    % 
    % 
    % Outputs:
    % c: geometric constants of the test rig, struct

    
    %% Plot options
    %Logical controlling plot behavior: optional, default=false
    if nargin<1
        plotPitch=false;
    end


    %Set figure storage folder
    dirFigs='Figures';


    %Make folder if it does not exist
    if ~isfolder(dirFigs)
        mkdir(dirFigs);
    end


    %% General
    c=struct();

    c.g=9.81;   %Gravitational acceleration
    
    c.d_p=174.494e-6;           %Particle diameter (GRANUSIL, Wedron IL #801, Grade 7020, to be confirmed)
    c.rho_p=SiO2.rho(293.15);   %Particle density
    c.eps_mf=0.45;              %Bed porosity at minimum fluidization conditions
    c.phi_s=0.8;                %Particle sphericity
    
    c.dh_eps1=50e-3;    %Height difference for porosity measurement, p1 and p3
    c.dh_eps2=100e-3;   %Height difference for porosity measurement, p2
    
    c.rho_N=DryAir.rho(1014e2,21+273.15);     %Reference density for flowmeter at 1014 hPa and 21°C
    c.Apipe=154.08e-3^2*pi/4;                 %Inner pipe area of main air line (6 inch, schedule 40 pipe)
    c.OnLimit=50e-3;                          %Minimum mass flow at which air supply is assumed to be on
    
    c.hBed=461.5e-3+14e-3;  %Persistent bed height, includes 14 mm from outlet weir to pressure probe that may not be filled with sand
    c.hFlow=344e-3;         %Vertical distance between baffles separating the chambers

    c.d_t2f=299.5e-3;   %Distance from the test tube's axis to the distributor floor
    c.d_t2b=20e-3;      %Distance from the test tube's axis to the baffle between 2. and 3. chamber (in the direction of the 2. chamber)
    
    c.dOrif=14.75e-3;     %Orifice plate inner diameter
    c.DOrif=82.5e-3;      %Orifice plate outer diameter (=inner pipe diameter)
    c.tap='D-D/2';        %Orifice plate pressure tap type
    
    c.x1=202e-3;        %Length of inlet / outlet chamber (first and fourth chamber)
    c.x2=1.06;          %Length of second chamber
    c.x3=0.8;           %Length of third chamber
    c.l=0.5;            %Width of all chambers
    
    c.Afloor1=c.x1*c.l;	    %Distributor floor area, inlet / outlet chamber (first and fourth chamber)
    c.Afloor2=c.x2*c.l;	    %Distributor floor area, second chamber
    c.Afloor3=c.x3*c.l;     %Distributor floor area, third chamber


    %% Test tube heat transfer area
    c.d_t=20e-3;    %Outside tube diameter (without fins)
    c.h_f=10e-3;    %Fin height
    c.s_f=2e-3;     %Fin thickness
    c.p_f=9e-3;     %Fin pitch
    l_lead=5e-3;    %Lead distance from the end of the tube to the start of the fins

    D_t=c.d_t+2*c.h_f;                                  %Outside tube diameter including fins
    c.phi=(D_t./c.d_t-1).*(1+0.35*log(D_t./c.d_t));     %Geometry function

    t=(c.l-2*l_lead)./c.p_f;            %Number of turns
    k=@(r) c.p_f./(2*r*pi);             %Fin slope depending on the reference radius r
    x=@(r) 2*r*pi*t.*sqrt(1+k(r).^2);   %Helix length depending on the reference radius r
    
    c.A_bottom=c.s_f.*x(c.d_t/2);           %Plain tube area where the fin is attached
    c.A_sides=2*c.h_f*x(c.d_t/2+c.h_f/2);   %Main fin area on the sides
    c.A_plain=c.d_t*pi.*c.l;                %Plain tube area


    %% Tube bundle void fraction
    %Definition based on VDI Heat Atlas, G7, p. 726, with w=w_p
    c.p_h=40e-3;    %Tube bundle horizontal pitch (raw)
    c.p_v=50e-3;    %Tube bundle vertical pitch

    b=c.p_h./(2*c.d_t);
    idx=b<1;
    
    c.psi=NaN(size(b));
    c.psi(idx)=1-pi*c.d_t(idx).^2./(4*c.p_h(idx).*c.p_v(idx));
    c.psi(~idx)=1-pi*c.d_t(~idx)./(8*c.p_v(~idx));


    %% Effective horizontal pitch
    %Projection into base plane as a function of axis variable z
    helix=@(z) (c.d_t/2+c.h_f)*sin(2*pi*z./c.p_f);


    %Start of z where projection is not covered by plain tube anymore
    fin_start=fzero(@(z) helix(z)-c.d_t./2,[0,c.p_f./4]);
    % fin_start=c.p_f/12;


    %Projected fin and plain tube areas
    A1=integral(helix,fin_start,c.p_f./4)-(c.p_f./4-fin_start).*c.d_t/2;
    % A1=c.p_f*((c.d_t/2+c.h_f)*sqrt(3)/(4*pi)-c.d_t/12);
    A2=c.s_f.*c.h_f;
    Atube=c.d_t.*c.l;

    Afins=2.*(2.*A1+A2).*t;     %Total projected fin area
    Aproj=Afins+Atube;          %Total projected area
    A0=c.p_h.*c.l;              %Free flow area


    %Effective pitch
    pi9=Aproj./A0;          %Tube packing density=relative blocked area
    c.p_hEff=c.d_t./pi9;    %Effective horizontal pitch


    %Plot pitch if requested
    if plotPitch
        %Set up figure
        fig=figure(907);
        clf(fig);
        t=tiledlayout(fig,1,1,'Padding','tight');
        ax=nexttile(t);
        colors=ax.ColorOrder;
        hold(ax,'on');
    

        %Plot individual areas
        legItems=cell(1,3);
        n=1000;
        mm=1000;
        

        %Fin edges (A1)
        z1=linspace(0,c.p_f,n);
        legItems{1}=area(z1*mm,helix(z1)*mm,'FaceColor',colors(1,:));
        area((z1+c.s_f)*mm,helix(z1)*mm,'FaceColor',colors(1,:));
    

        %Area between edges (A2)
        z2=linspace(c.p_f/4,c.p_f/4+c.s_f,n);
        y2=repmat(c.d_t/2+c.h_f,1,n);
        legItems{2}=area(z2*mm,y2*mm,'FaceColor',colors(2,:));
        area((z2+c.p_f/2)*mm,-y2*mm,'FaceColor',colors(2,:));
    

        %Plain tube (A3)
        z3=linspace(0,c.p_f+c.s_f,n);
        y3=repmat(c.d_t/2,1,n);
        legItems{3}=area(z3*mm,y3*mm,'FaceColor',colors(3,:));
        area(z3*mm,-y3*mm,'FaceColor',colors(3,:));


        %Hidden fin edges
        plot(ax,z1*mm,helix(z1)*mm,'Color','k','LineStyle','--');
        plot(ax,(z3+c.s_f)*mm,helix(z3)*mm,'Color','k','LineStyle','--');


        %Start of fin affecting projected area
        xline(ax,fin_start*mm);
        text(ax,fin_start*mm,-15,' fin_{start}');


        %Peak of left fin edge 
        xline(ax,c.p_f/4*mm);
        text(ax,c.p_f/4*mm,-15,' p_f / 4');


        %Peak of right fin edge 
        xline(ax,(c.p_f/4+c.s_f)*mm);
        text(ax,(c.p_f/4+c.s_f)*mm,-15,' p_f / 4 + s_f');


        %Tube axis
        yline(ax,0,'LineStyle','-.');
    
        hold(ax,'off');
    

        %Axes limits
        ax.XLim=[0,(c.p_f+c.s_f)*mm];
    

        %Set legend
        legItems=vertcat(legItems{:});
        legend(ax,legItems,{'A_1','A_2','A_{tube}'},'Location','northeast');
    

        %Axes labels
        xlabel(ax,'z (mm)');
        ylabel(ax,'y (mm)');
    

        %Export figure for repository
        t.Units='centimeters';
        t.OuterPosition=[0,0,17,8.5];

        fig.Units=t.Units;
        fig.Position(3:4)=t.OuterPosition(3:4)+0.5;

        exportgraphics(fig,[dirFigs,filesep,'effectivePitch.tiff'],...
            'Resolution',600);
    end
end




