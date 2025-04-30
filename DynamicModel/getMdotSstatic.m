%% Calculate Static Particle Mass Flow Boundary Condition
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
%This function creates the vector of specific particle mass flows for each 
%cell required to run static simulations 
%
%
%Requires all auxiliary classes and functions on the MATLAB path
%
%Required products, version 24.1:
%   - MATLAB
%Necessary files, classes, functions, and scripts:
%   - @DryAir
%   - loadGeometry.m


function mDotS=getMdotSstatic(mDotSin,direction)
    loadGeometry;   %Basic geometry


    %Find inlet chamber
    if direction
        idx=1;
    else
        idx=length(nChambers);
    end
    

    %Set up vectors in and out of cells based on particle flow direction
    mDotSR=[repmat(mDotSin/nChambers(idx),1,nChambers(idx)),zeros(1,n-nChambers(idx))];
    mDotSL=[0,mDotSR(1:end-1)];
    
    mDotSR=cumsum(mDotSR);
    mDotSL=cumsum(mDotSL);
    

    %Mean particle mass flow in each cell
    mDotS=mean([mDotSR;mDotSL],1);
    

    %Flip and set negative if flow direction is reversed
    if ~direction
        mDotS=-fliplr(mDotS);
    end

end