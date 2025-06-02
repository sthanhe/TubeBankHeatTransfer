%% Analyse degrees of fluidization (primary data)
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
% This script calls the dynamic numerical model contained in the folder
% "DynamicModel" to analyze the degrees of fluidization at the test tube.
% See the Methodology Report in the data repository for detailed
% explanations. It creates the file "h2FG.mat" required by the script
% "prepPrim". The dynamic numerical model is a slight adaptation of a model
% published in a previous study: https://doi.org/10.5281/zenodo.7948224
%
%
%Requires all auxiliary classes and functions on the MATLAB path
%
%Required products, version 24.1:
%   - MATLAB
%   - Simulink
%   - Requirements Toolbox
%   - Simulink Real-Time
%   - Stateflow
%Necessary classes, functions, files, and scripts:
%   - Everything contained in the folder "DynamicModel" in the software
%       repository


%% Set data directories
dirMdl='DynamicModel';              %Dynamic model storage folder
dirData='Data/Own';                 %Data storage folder
dirFigures='Figures/DynamicSims';   %Figure storage folder
dirTemp='Temp';                     %Temporary data storage folder


%% Create storage folders if they do not exist
if ~isfolder(dirFigures)
    mkdir(dirFigures);
end

if ~isfolder(dirTemp)
    mkdir(dirTemp);
end


%% Set preferences
%Simulink Data Inspector
storage=Simulink.sdi.getStorageLocation;
limit=Simulink.sdi.getArchiveRunLimit;

Simulink.sdi.setStorageLocation(what(dirTemp).path);
Simulink.sdi.setArchiveRunLimit(0);


%Get full paths of data directories
dirMdl=what(dirMdl).path;
dirData=what(dirData).path;
dirFigures=what(dirFigures).path;
if ~isempty(storage)
    storage=what(storage).path;
end


%Change current folder to dynamic model
dirOrig=cd(dirMdl);


%% Do analysis
%Prepare data
prepStatic;


%Do simulations
calcStatic;


%% Revert preferences
%Simulink Data Inspector
Simulink.sdi.setStorageLocation(storage);
Simulink.sdi.setArchiveRunLimit(limit);


%Change back to original current folder
cd(dirOrig);




