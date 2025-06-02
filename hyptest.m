%% Hypothesis testing
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
% This function conducts the hypothesis tests as described in the 
% Methodology Report. 
%
%
%Requires all auxiliary classes and functions on the MATLAB path
%
%Required products, version 24.1:
%   - MATLAB
%   - Statistics and Machine Learning Toolbox
%Necessary classes, functions, files, and scripts:
%   - None


function para=hyptest(X,y,fx,beta,beta0,idx)
    % Inputs:
    % X: predictor variables (n,m)
    % y: Response variable (n,1)
    % beta: estimated parameters
    % beta0: estimated parameters of the null hypothesis
    % idx: parameter indices to be included
    % 
    % 
    % Outputs:
    % para: results of the hypothesis tests, table


    %Set up table
    names={'Parameter','H0','Estimate','p'};
    para=table('Size',[length(idx)+1,length(names)],...
        'VariableTypes',[{'string'},repmat({'double'},1,length(names)-1)],...
        'VariableNames',names);
    
    
    %Fill basic values
    para.Parameter=[compose('P%d',idx),{'All'}]';   %Parameter names
    para.H0(1:end-1)=beta0(idx);                    %Null hypothesis
    para.Estimate(1:end-1)=beta;                    %Estimated parameters
    
    
    %p-values of individual parameters
    for i=1:height(para)-1
        beta_null=para.Estimate;
        beta_null(i)=para.H0(i);
    
        [~,p]=ttest(fx(beta_null,X),y);
        para.p(i)=p;
    end
    
    
    %p-value of all parameters combined
    beta_null=para.H0;
    [~,para.p(end)]=ttest(fx(beta_null,X),y);
end




