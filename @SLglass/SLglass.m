%% Property functions of soda-lime glass
%
%GNU General Public License v3.0
%By Stefan Thanheiser: https://orcid.org/0000-0003-2765-1156
%
%Based on:
%J. Huang, P.K. Gupta,
%Temperature dependence of the isostructural heat capacity of a soda lime silicate glass,
%Journal of Non-Crystalline Solids, Volume 139,1992,Pages 239-247,
%https://doi.org/10.1016/S0022-3093(05)80831-2
%
%All parameters and results in SI base units


classdef SLglass
    methods(Static)
        function c_p=c_p(T)
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