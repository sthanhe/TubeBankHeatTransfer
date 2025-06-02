%% Coefficient of determination R²
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
% This function calculates the coefficient of determination R².
%
%
%Requires all auxiliary classes and functions on the MATLAB path
%
%Required products, version 24.1:
%   - MATLAB
%Necessary classes, functions, files, and scripts:
%   - None


function Rsq=Rsq(y,yHat,M)
    % Inputs:
    % y: observed data
    % yHat: estimated / predicted / modeled / fitted data
    % M: (optional) number of regressors excluding intercept
    % 
    % 
    % Outputs:
    % Rsq: coefficient of determination. If the function is called with a
    %       number of regressors M (third parameter), Rsq is the adjusted
    %       coefficient of determination
    

    %Normalize inputs
    n=numel(y);
    y=reshape(y,n,1);           
    yHat=reshape(yHat,n,1);     


    %Regular coefficient of determination
    epsHat=yHat-y;                  %Residuals
    yBar=mean(y,'omitmissing');     %Mean of observed data

    SSE=sum(epsHat.^2,'omitmissing');       %Sum of squared errors
    SST=sum((y-yBar).^2,'omitmissing');     %Total sum of squares

    Rsq=1-SSE./SST;     %Coefficient of determination


    %Adjusted coefficient of determination
    if nargin>2
        Rsq=1-(n-1)./(n-M).*(1-Rsq);
    end
end




