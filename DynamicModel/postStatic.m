%% Post Processing of Stationary Simulations
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
% Identical to the file of the same name in:
%
% S. Thanheiser, Particle Dispersion Model Software. (Feb. 07, 2025). 
% Zenodo. doi: 10.5281/zenodo.14833128.
%
%
%
%This function anaylzes the results of the stationary simulations conducted
%with the script "calcStationary" and creates all published figures.
%
%
%Requires all auxiliary classes and functions on the MATLAB path
%
%Required products, version 24.1:
%   - MATLAB
%Necessary files, classes, functions, and scripts:
%   - None


function postStatic(out,flow,x,xChambers,figidx,dirFigures)
    %% Common parameters
    xmeas=[181.5e-3,37e-3,1023e-3,37e-3,763e-3,37e-3];
    xmeas=cumsum(xmeas);    %x-coordinates of bed level measurements


    %Figure title:
    titleText=['Test ',num2str(figidx),...
                ', $FG_{2}$=',num2str(round(flow.FG2,1)),...
                ', $\dot{m}_p$=',num2str(flow.mDotSand),' kg/s',...
                ', $T_{bed2}$=',num2str(round(flow.Tbed2-273.15)),'$^{\circ}$C'];


    %% Bed level and degrees of fluidization
    hmeas=flow{:,compose('h%d',6:-1:1)};    %Measured bed levels
    hsim=timeseries2timetable(out.h);       %Simulated bed levels
    FG=timeseries2timetable(out.FG);        %Degree of fluidizations simulated


    %Set up figure
    fig915=figure(figidx);
    clf(fig915);
    t=tiledlayout(fig915,1,1,'Padding','tight');
    ax=nexttile(t);
    box(ax,'on');
    hold(ax,'on');
    colors=ax.ColorOrder;


    %Bed levels
    plot(ax,x,hsim{end,:},'Color',colors(1,:));
    plot(ax,xmeas,hmeas,'Color',colors(2,:));
    ylabel(ax,'Bed Level (m)');


    %Excess fluidization
    yyaxis(ax,'right');
    plot(ax,x,FG{end,:})
    ylabel(ax,'FG (-)');
    

    %Baffle positions
    xline(cumsum(xChambers(1:end-1)));
    hold(ax,'off');

    
    %Configure axes, legend, and title
    ax.XLim=[0,max(x)];
    legend(ax,{'Simulated','Measured'});

    title(ax,titleText,'Interpreter','latex');
    xlabel(ax,'Distance from Inlet (m)');
    

    %Export figure for repository
    t.Units='centimeters';
    t.OuterPosition=[0,0,17,8.5];

    fig915.Units=t.Units;
    fig915.Position(3:4)=t.OuterPosition(3:4)+0.5;

    exportgraphics(fig915,[dirFigures,filesep,'static',num2str(figidx),'.tiff'],...
        'Resolution',600);


    % %% Fluidization
    % 
    % D=timeseries2timetable(out.D);      %Mass diffusivity simulated
    % 
    % 
    % %Set up figure
    % fig916=figure(figidx+100);
    % clf(fig916);
    % ax=gca();
    % box(ax,'on');
    % 
    % 
    % %Excess fluidization
    % plot(ax,x,FG{end,:}-1)
    % ylabel(ax,'w_e/w_{mf} (-)');
    % 
    % 
    % %Particle dispersion
    % hold(ax,'on');
    % yyaxis(ax,'right');
    % plot(ax,x,squeeze(D{end,:}))
    % ylabel(ax,'D (m^2/s)');
    % 
    % 
    % %Baffle positions
    % xline(cumsum(xChambers(1:end-1)));
    % hold(ax,'off');
    % 
    % 
    % %Configure and save figure
    % title(ax,titleText,'Interpreter','latex');
    % xlabel(ax,'Distance from Inlet (m)');
    % 
    % fig916.Units='centimeters';
    % fig916.Position=[0.02,0.83,17,8.5];
    % ax.XLim=[0,max(x)];
    % 
    % exportgraphics(fig916,[dirFigures,filesep,'statFluidization',num2str(figidx),'.tiff'],...
    %     'Resolution',600);
end




