%% Property functions of soda-lime glass
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
% All required files for this class can be found in the software
% repository: see the link to the supplemental release in the data 
% repository here: https://doi.org/10.5281/zenodo.15576311
% 
% 
%
% This class describes the thermo-physical properties of soda-lime glass 
% according to:
% 
% J. Huang, P.K. Gupta,
% Temperature dependence of the isostructural heat capacity of a soda lime 
% silicate glass, Journal of Non-Crystalline Solids, Volume 139, 1992, 
% Pages 239-247, https://doi.org/10.1016/S0022-3093(05)80831-2
%
% 
%Requires all files packaged in the class folder and on the MATLAB path
%
%Required products, version 24.1:
%   - MATLAB



classdef SLglass
    %All parameters and results in SI base units


    %% Property functions
    methods(Static)
        function c_p=c_p(T)
            %Specific isobaric heat capacity
            
            persistent A B C
            if isempty(A)
                A=0.828e3;
                B=4.418e-1;
                C=17185.9e3;
            end

            T(T<273.15 | T>650)=NaN;

            c_p=A+B.*T-C./T.^2;
        end
    end
end




