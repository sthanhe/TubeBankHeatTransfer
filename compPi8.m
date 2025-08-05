%% Compare results to other models for probe size (pi8)
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
% This function compares the results from particle-convective regressions 
% in "calcPC" to other published models regarding the impact of probe size
% on the wall-to-bed HTC. 
%
%
%Requires all auxiliary classes and functions on the MATLAB path
%
%Required products, version 24.1:
%   - MATLAB
%   - Curve Fitting Toolbox
%Necessary classes, functions, files, and scripts:
%   - @DryAir
%   - @FluBed
%   - @SiO2
%   - @figaux
%   - @implExp
%   - Nu_rel.m


function compPi8(fx,beta,dirFigs,figidx,small)
    % Inputs:
    % fx: model function in the form fx(b,X), where b are the regression
    %       coefficients, function handle
    % beta: regression coefficients, double
    % dirFigs: path to directore where figures should be stored, char
    % figidx: index of figure window, double
    % small: indicator whether to create the small figure version, logical


    %% Extended model
    %Basic parameters
    p=1e5;      %Pressure
    T=293.15;   %Temperature
    
    Ar=[1e3,1e4];   %Archimedes numbers
    eps_mf=0.45;    %Minimum fluidization voidage
    phi_s=0.8;      %Particle sphericity
    
    p_h=Inf;            %Horizontal tube spacing: single tube
    w_p=0;              %No particle cross-flow
    c_p=SiO2.c_p(T);    %Specific heat capacity: Silica particles
    

    %Tube diameter: higher density at low diameters
    d_t=linspace(0,5e-3,100);
    d_t=[d_t,linspace(d_t(end),60e-3,100),Inf];     
    
    
    %Relative maximum Nusselt number
    Nu_relExt=Nu_rel(fx,beta,p,T,Ar,eps_mf,phi_s,d_t,w_p,c_p,p_h);
    
    
    %% Alternative models: relative maximum Nusselt numbers
    %Molerus and Wirth
    C=0.85;
    u_l=0.3e-2;
    f_L=1;
    
    Nu_relMW=C.*u_l./(f_L.*d_t)+1;
    
    
    %Grewal and Saxena
    d_tMaxGS=28.6e-3;
    Nu_relGS=d_t.^-0.21./d_tMaxGS.^-0.21;
    
    
    %Merzsch et al.
    d_tMaxMerzsch=33.7e-3; 
    Nu_relMerzsch=d_t.^-0.3./d_tMaxMerzsch.^-0.3;
    
    
    %% Plot graphic
    %Set up figure
    if figidx==0
        fig=figure(100);
    else
        fig=figure(figidx);
    end
    clf(fig);
    t=tiledlayout(fig,1,1,'Padding','tight');
    ax=nexttile(t);
    hold(ax,'on');
    
    
    %Plot lines
    plot(ax,d_t.*10^3,Nu_relExt);
    plot(ax,d_t.*10^3,[Nu_relMW;Nu_relGS;Nu_relMerzsch],'LineStyle','--');
    
    hold(ax,'off');
    
    
    %Set axes limits
    ax.XLim=[0,40];
    ax.YLim=[1,3.5];
    
    
    %Configure and print figure: small for repository, regular for
    %manuscript
    if small
        %Turn off axes visibility
        ax.Visible='off';


        %Export figure for repository
        fname=[dirFigs,filesep,'pi8_s',num2str(figidx)];
        
        t.Units='centimeters';
        t.InnerPosition=[0.5,0.5,4,4];
        
        fig.Units=t.Units;
        fig.Position(3:4)=t.OuterPosition(3:4)+0.5;
        
        exportgraphics(fig,[fname,'.tiff'],'Resolution',600);


    else
        %Set legend
        lgd=legend(ax,figaux.subsz(...
            [compose('Extended Model, Ar=1e%d',log10(Ar'));...
            {'Molerus & Wirth';...
                ['Grewal & Saxena, d_{t,max} = ',num2str(d_tMaxGS*10^3), ' mm'];...
                ['Merzsch et al, d_{t,max} = ',num2str(d_tMaxMerzsch*10^3), ' mm']}],...
                6),...
            'Location','north');
        

        %Set axes labels
        xlabel(ax,figaux.subsz('Tube diameter d_t (mm)',6));
        ylabel(ax,figaux.subsz('Nu_{max} / Nu_{max} (d_t \rightarrow \infty)',6));
        
        
        %Text size
        ax.FontSize=7;
        lgd.FontSize=7;


        %Export figure for manuscript
        fname=[dirFigs,filesep,'Figure',num2str(figidx)];
        
        t.Units='centimeters';
        t.OuterPosition=[0,0,9,9];
        
        fig.Units=t.Units;
        fig.Position(3:4)=t.OuterPosition(3:4)+0.5;
        
        exportgraphics(fig,[fname,'.tiff'],'Resolution',600);
        
        
        %Export figure for Elsevier
        t.OuterPosition=[0,0,9,9];
        fig.Position(3:4)=t.OuterPosition(3:4)+0.5;
        
        exportgraphics(fig,[fname,'.eps']);
        savefig(fig,fname);
    end
    

end




