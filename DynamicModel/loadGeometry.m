%% Geometrical Parameters for the Dynamic Simulink Model
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
%All data, along with methodology reports and supplementary documentation, 
%is published in the data repository:
%https://doi.org/10.5281/zenodo.7924693
%
%All required files for this script can be found in the software
%repository:
%https://doi.org/10.5281/zenodo.7948224
%
%
%
%This script creates the geometrical constants in the base workspace, which
%are necessary to run the dynamic numerical model "dynamicModel.slx". 
%
%
%Requires all auxiliary classes and functions on the MATLAB path
%
%Required products, version 24.1:
%   - MATLAB
%Necessary files, classes, functions, and scripts:
%   - None


%% Geometrical parameters
x1=200e-3;  %Length of inlet / outlet chamber (first and fourth chamber). Rounded from 202 mm to 200 mm
x2=1.06;    %Length of second chamber
x3=0.8;     %Length of third chamber


xChambers=[x1,x2,x3,x1];    %Chamber lengths
xPress=37e-3;               %Distance between pressure taps for bed level measurement (next to baffles)
deltaX=20e-3;               %Element (discretization) length

l=0.5;              %Inner width of bed
hTotal=1.4715;      %Total height from floor to top of freeboard
hGate=461.5e-3;     %Height from floor to outlet weir

s=5e-3;                 %Sinter floor thickness
sinter="SIKA-R 15 AX";  %Sinter floor name (not functional)

Kvs=220;                    %Kvs value of air cushion valves in m³/h for water (1000 kg/m³) at a pressure difference of 1 bar
dACpipe=88.9e-3-2*3.2e-3;   %Inner diameter of air cushion valve connection pipe

d_p=174.494e-6;     %Particle diameter
rho_p=2650;         %Particle density
eps_mf=0.45;        %Porosity at minimum fluidization conditions

href=1;             %Reference height (mass diffusivity)


%Air box volumes, measured from CAD model
VABs=NaN(1,length(xChambers));
VABs(1)=109114.741e-6*160.5e-3+83595.112e-6*cos(0.15)*(200e-3-160.5e-3)/2+(109114.741e-6-83595.112e-6*cos(0.15))*(200e-3-160.5e-3);
VABs(2)=109114.741e-6*923.5e-3+88196.654e-6*cos(0.37)*(1023e-3-923.5e-3)/2+(109114.741e-6-88196.654e-6*cos(0.37))*(1023e-3-923.5e-3)+620e-3*32.5e-3*4e-3-4*33546.891e-6*4e-3;
VABs(3)=109114.741e-6*663e-3+88196.654e-6*cos(0.37)*(763e-3-663e-3)/2+(109114.741e-6-88196.654e-6*cos(0.37))*(763e-3-663e-3)+620e-3*32.5e-3*4e-3-3*33546.891e-6*4e-3;
VABs(4)=VABs(1);


%Derived values
n=round(sum(xChambers)./deltaX);                    %Number of elements
nChambers=round(xChambers./deltaX);                 %Number of elements per chamber
x=linspace(deltaX/2,sum(xChambers)-deltaX/2,n);     %x-coordinates of cell centers

%Airboxes
nABs=length(xChambers);         %Number of air boxes
posABcell=arrayfun(@(x) repmat(x,1,nChambers(x)),1:nABs,'UniformOutput',false);
posABs=cell2mat(posABcell);     %Position vector of airboxes. For every cell, the number indicates the index of the airbox assigned to the cell

%Air cushions
nACs=nABs-2;                            %Number of air cushions
posACs=posABs-1;
posACs(n-nChambers(end)+1:end)=0;       %Position vector of air cushion. For every cell, the number indicates the index of the air cushion assigned to the cell. 0=no air cushion

Kvs=Kvs.*(1/60^2*sqrt(1000/0.1e6));     %Transform from m³/h to m²
A=dACpipe.^2.*pi./4;                    %Cross section of air cushion pipe

%Bed level measurements
lCum=cumsum(xChambers(1:end-1));
posBedLevel=sort([lCum-xPress./2,lCum+xPress./2]);
posBedLevel=round(interp1(x,1:n,posBedLevel));      %Cell indices in which the bed level gets measured
posBedLevelForward=posBedLevel(3:2:end);            %Cell indices of bed level measurements when the direction is forward (left to right)
posBedLevelReverse=posBedLevel(2:2:end-1);          %Cell indices of bed level measurements when the direction is reverse (right to left)




