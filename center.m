%% Bin and center values
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
% This function bins values around given center points.
%
%
%Requires all auxiliary classes and functions on the MATLAB path
%
%Required products, version 24.1:
%   - MATLAB
%Necessary classes, functions, files, and scripts:
%   - None


function [x,centers,n]=center(x,approx)
    % Inputs:
    % x         Values to bin and center
    % approx    Approximate centers of bins. The size of this array
    %           determines the number of bins
    % 
    % Outputs:
    % x         Mean (=center) of the bin to which each value was assigned
    % centers   Mean value of each bin
    % n         Number of values in each bin


    %Create bins by assigning each value to the nearest approximated center
    bins=interp1(approx,approx,x,'nearest','extrap');


    %Calculate centers = mean of each bin
    centers=arrayfun(@(y) mean(x(bins==y)),approx);
    centers(isnan(centers))=[];


    %Assign each value its center
    if length(centers)>1
        x=interp1(centers,centers,x,'nearest','extrap');
    else
        x=repmat(centers,size(x));
    end


    %Number of values in each bin
    n=arrayfun(@(i) nnz(x==i),centers);
end




