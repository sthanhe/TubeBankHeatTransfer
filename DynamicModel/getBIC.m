%% Retrieve Boundary and Initial Conditions for Dynamic Model
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
%All required files for this function can be found in the software
%repository:
%https://doi.org/10.5281/zenodo.7948224
%
%
%
%This function creates the boundary and initial conditions necessary to run
%the dynamic numerical model "dynamicModel.slx". 
%
%
%Requires all auxiliary classes and functions on the MATLAB path
%
%Required products, version 24.1:
%   - MATLAB
%Necessary files, classes, functions, and scripts:
%   - @DryAir
%   - loadGeometry.m


function [bc,Phi,mAC,HAC,mAB]=getBIC(flow,direction)
    %% Basic geometry
    loadGeometry;


    %% Boundary conditions
    p0=flow.p0;    %Ambient pressure
    
    Time=duration(0,0,0);   %Constant values
    
    Tbed=vals2cells(x,Time,flow.Tleft,flow.Tcenter,flow.Tright,nChambers);          %Temperature distribution
    poro=vals2cells(x,Time,flow.epsLeft,flow.epsCenter,flow.epsRight,nChambers);    %Porosity distribution
    
    
    mDotSand=zeros(1,n);
    if direction
        mDotSand(1:nChambers(1))=repmat(flow.mDotS./(deltaX*nChambers(1)),1,nChambers(1));
    else
        mDotSand(n-nChambers(end)+1:end)=repmat(flow.mDotS./(deltaX*nChambers(end)),1,nChambers(end)); %#ok<FNCOLND>
    end
    S=timetable(Time,mDotSand);     %Sand flow, distributed evenly across inlet chamber
    

    mDotABin=timetable(Time,flow{1,compose('air%d',1:nABs)});   %Air mass flows
    hSet=timetable(Time,flow{1,compose('AC%dset',1:nACs)});     %PID controller setpoint
    dirTT=timetable(Time,double(direction));                    %Sand flow direction


    %Set up dataset
    bc=Simulink.SimulationData.Dataset;
    bc=addElement(bc,Tbed,'Tbed');
    bc=addElement(bc,S,'S');
    bc=addElement(bc,poro,'poro');
    bc=addElement(bc,hSet,'hSet');
    bc=addElement(bc,mDotABin,'mDotABin');
    bc=addElement(bc,dirTT,'direction');
    

    %% Initial conditions
    h=repmat(500e-3,1,n);   %Sand levels
    Tbed1=Tbed{1,1};        %Bed temperature
    
    
    %Fluidized bed
    Phi=h./href.*rho_p.*(1-poro{1,1});     %Fictional density
    
    
    %Air cushions
    Vvec=(hTotal-h).*l.*deltaX;
    VAC=arrayfun(@(x) sum(Vvec(posACs==x)),1:nACs);
    TAC=arrayfun(@(x) mean(Tbed1(posACs==x)),1:nACs);
    
    mAC=VAC.*DryAir.rho(p0,TAC);    %Air mass in air cushions
    HAC=mAC.*DryAir.h(TAC);         %Air enthalpy in air cushions
    
    
    %Air boxes    
    TAB=arrayfun(@(x) mean(Tbed1(posABs==x)),1:nABs);
    mAB=p0.*VABs./(DryAir.R.*TAB);      %Air mass in air boxes
end


%% Auxiliary function: distribute measurements linearly across all cells
function TT=vals2cells(x,Time,left,center,right,nChambers)
    n=length(left);
    nCum=cumsum(nChambers);

    y=NaN(n,sum(nChambers));
    y(:,1:nCum(1))=repmat(left,1,nChambers(1));
    
    x1=[x(nCum(1)+1),x(nCum(2))];
    xq=x(nCum(1)+1:nCum(2));
    T1=arrayfun(@(i) interp1(x1,[left(i),center(i)],xq),1:n,'UniformOutput',false);
    y(:,nCum(1)+1:nCum(2))=cell2mat(T1');
    
    x1=[x(nCum(2)+1),x(nCum(3))];
    xq=x(nCum(2)+1:nCum(3));
    T2=arrayfun(@(i) interp1(x1,[center(i),right(i)],xq),1:n,'UniformOutput',false);
    y(:,nCum(2)+1:nCum(3))=cell2mat(T2');
    
    y(:,nCum(3)+1:end)=repmat(right,1,nChambers(4));
    
    TT=timetable(Time,y);
end




