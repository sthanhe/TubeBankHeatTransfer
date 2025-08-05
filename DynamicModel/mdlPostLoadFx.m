%% Dynamic Model Post Load Function
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
% Slight adaptation of the file of the same name in:
%
% S. Thanheiser, Particle Dispersion Model Software. (Feb. 07, 2025). 
% Zenodo. doi: 10.5281/zenodo.14833128.
%
%
%
%This script creates the workspace variables and an initial set of 
%boundary and initial conditions necessary to run the dynamic numerical 
%model "dynamicModel.slx". 
%
%
%Requires all auxiliary classes and functions on the MATLAB path
%
%Required products, version 24.1:
%   - MATLAB
%Necessary files, classes, functions, and scripts:
%   - @DryAir
%   - @FluBed
%   - @implExp
%   - @Sinter
%   - getBIC.m
%   - loadGeometry.m
%   - getMdotSstatic.m


%% Basic geometry
loadGeometry;


%% Default values
direction=true;         %Flow direction
baffleCorr=ones(1,3);   %Baffle correction factors


%% Initial set of boundary and inital conditions
%Set up table
names={'p0','Tleft','Tcenter','Tright','epsLeft',...
    'epsCenter','epsRight','mDotS','air1','air2',...
    'air3','air4','AC1set','AC2set'};
init=table('Size',[1,length(names)],...
    'VariableTypes',[repmat({'double'},1,length(names))]);
init.Properties.VariableNames=names;
clear('names');


%Take values directly from a measurement point
init.p0=101322.321749582;
init.Tleft=321.622691993331;
init.Tcenter=319.331017705775;
init.Tright=330.764596384106;
init.epsLeft=0.467551270924290;
init.epsCenter=0.470101396127068;
init.epsRight=0.472651521329844;
init.mDotS=4;
init.air1=0.0141454434627567;
init.air2=0.0425907793606724;
init.air3=0.0370863059211219;
init.air4=0.0121968546724413;
init.AC1set=1;
init.AC2set=1;


%Get boundary and initial conditions
[bc,Phi,mAC,HAC,mAB]=getBIC(init(1,:),direction);


%Other boundary and initial conditions
p0=init.p0(1);  %Ambient pressure
Phigate=hGate./href.*rho_p.*(1-init.epsRight(1))+p0./(FluBed.g.*href);  %Weir boundary condition
Y0=0.6*ones(1,nACs);    %PID I-value


%% Simulation time
minSimTime=240;     %Minimum simulation time (seconds)
maxSimTime=600;     %Maximum simulation time (seconds)
statCond=1e-4;      %Condition for stationary status as a value of PhiDot


%% Static simulation parameters
mDotSstatic=getMdotSstatic(init.mDotS,direction);
isStatic=false;

YMan=[1,1];
SetMan=false;


%% Particle dispersion coefficients
c=19317.8340011694;
eps2=1.10169319091282;
eps3=-2.04940397390448;
epsAr=0.108571357573225;




