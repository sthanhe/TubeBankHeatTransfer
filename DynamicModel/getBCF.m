%% Calculate Baffle Correction Factors
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
%All required files for this script can be found in the software
%repository:
%https://doi.org/10.5281/zenodo.7948224
%
%
%
%This script calculates the baffle correction factors for a specific
%baffle when the air cushion valve is inactive. It gets called by the
%"baffleCalib" script. 
%
%
%Required products, version 24.1:
%   - MATLAB
%   - Simulink
%   - Requirements Toolbox
%   - Simulink Real-Time
%   - Stateflow
%Necessary files, classes, functions, and scripts:
%   - @DryAir
%   - @FluBed
%   - @implExp
%   - @Sinter
%   - getBIC.m
%   - mdlPostLoadFx.m
%   - loadGeometry.m
%   - getMdotSstatic.m
%   - dynamicModel.slx


%% Baffle Correction Factor
for i=run
    if isnan(flow.h4(i))
        continue;
    end

    %Boundary conditions
    p0=flow.p0(i);
    baffleCorr=baffleMat(i,:);
    Phigate=flow.Phigate(i);
    bc=getBIC(flow(i,:),direction);


    %Simulate low baffleCorr estimate
    baffleCorr(baffleIdx)=tab.baffleLow(i);

    out=sim(sys,'LoadExternalInput','on','ExternalInput','bc',...
                'LoadInitialState','on','InitialState','xInit');

    hsim=timeseries2timetable(out.h);
    tab.hLow(i)=hsim.hOverX(end,bafflePos);     %Low h estimate


    %Continue if baffle has no effect
    if flow{i,hname}<tab.hLow(i)
        continue;
    end


    %Simulate high baffleCorr estimate
    baffleCorr(baffleIdx)=tab.baffleHigh(i);

    out=sim(sys,'LoadExternalInput','on','ExternalInput','bc',...
                'LoadInitialState','on','InitialState','xInit');

    hsim=timeseries2timetable(out.h);
    tab.hHigh(i)=hsim.hOverX(end,bafflePos);     %High h estimate


    %Get baffleCorr from estimates through linear interpolation
    baffleMat(i,baffleIdx)=interp1([tab.hLow(i),tab.hHigh(i)],...
                        [tab.baffleLow(i),tab.baffleHigh(i)],...
                        flow{i,hname},...
                        'linear','extrap'); %#ok<SAGROW>
end




