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
    
    c.d_p=174.494e-6;   %Particle diameter (GRANUSIL, Wedron IL #801, Grade 7020, to be confirmed)
    c.rho_p=2650;       %Particle density
    c.epsMf=0.45;       %Bed porosity at minimum fluidization conditions
    
    c.dh_eps1=50e-3;    %Height diffc.xPresserence for porosity measurement, p1 and p3
    c.dh_eps2=100e-3;   %Height difference for porosity measurement, p2
    
    c.rho_N=DryAir.rho(1014e2,21+273.15);     %Reference density for flowmeter at 1014 hPa and 21°C
    c.Apipe=154.08e-3^2*pi/4;                 %Inner pipe area of main air line (6 inch, schedule 40 pipe)
    c.OnLimit=50e-3;                          %Minimum mass flow at which air supply is assumed to be on
    
    c.hBed=461.5e-3+14e-3;  %Persistent bed height, includes 14 mm from outlet weir to pressure probe that may not be filled with sand
    c.hRef=1;               %Reference bed height
    c.hRef0=1;              %Reference bed height, normalized
    
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

    c.dTube=20e-3;                  %Tube diameter
    c.AtestTube=c.dTube*pi*c.l;     %Test tube plain outside surface area
end




