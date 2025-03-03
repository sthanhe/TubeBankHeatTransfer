classdef DC04
    methods(Static)
        function lambda=lambda(T)
            load(['@DC04',filesep,'fits.mat'],'T_lambda');

            valid=200<=T & T<=600;    %Estimated range of validity

            lambda=NaN(size(T));
            lambda(valid)=T_lambda(T(valid));
        end


        function createConstants()
            tab=readtable(['@DC04',filesep,'prop.xls']);
            tab.T=tab.T+273.15;

            T_lambda=DC04.createFit(tab.T,tab.lambda);

            save(['@DC04',filesep,'fits.mat'],'T_lambda');
        end
    end


    methods(Static, Access=private)
        [fitresult,gof]=createFit(T,lambda)
    end
end