%% Univariate plot
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
% All data, along with methodology reports and supplementary documentation, 
% is published in the data repository:
% https://doi.org/10.5281/zenodo.15576311
%
% All required files for this function can be found in the software
% repository: 
% https://doi.org/10.5281/zenodo.15576950
%
%
%
% This function creates a plot to illustrate the variance within variables.
%
%
%Requires all auxiliary classes and functions on the MATLAB path
%
%Required products, version 24.1:
%   - MATLAB
%Necessary classes, functions, files, and scripts:
%   - None


function [fig,t,ax]=uniplot(X,figidx)
    % Inputs:
    % X: variable matrix, each column is a variable, double
    % figidx: index of figure window, double
    % 
    % 
    % Outputs:
    % fig: figure handle
    % t: tile handle
    % ax: axis handle


    %Normalize (z-score)
    Xnorm=normalize(X);
    
    
    %Set constants to 0
    isconst=arrayfun(@(i) ...
        isscalar(unique(X(:,i))),...
        1:size(Xnorm,2));
    Xnorm(:,isconst)=0;
    
    
    %Set up figure
    fig=figure(figidx);
    clf(fig);
    t=tiledlayout(fig,1,1,'Padding','tight');
    ax=nexttile(t);
    

    %Plot chart
    boxchart(ax,Xnorm);
    
    
    %y-axis label
    ylabel(ax,'z-score (-)');
    
    
    %Size figure
    t.Units='centimeters';
    t.OuterPosition=[0,0,17,8.5];
    
    fig.Units=t.Units;
    fig.Position(3:4)=t.OuterPosition(3:4)+0.5;
end




