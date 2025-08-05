%% Compare results to other models for tube packing density (pi9)
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
% in "calcPC" to other published models regarding the impact of tube
% packing density on the wall-to-bed HTC. 
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


function compPi9(fx,beta,dirFigs,figidx,small)
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
    
    d_t=25e-3;          %Tube diameter
    w_p=0;              %No particle cross-flow
    c_p=SiO2.c_p(T);    %Specific heat capacity: Silica particles
    

    %Horizontal spacing: tighter around sqrt(2)/2
    s_h=linspace(1,0.95,100);
    s_h=[s_h,linspace(0.95,sqrt(2)/2+1e-3,100)];
    s_h=[s_h,linspace(sqrt(2)/2+1e-3,sqrt(2)/2-1e-2,1000)];
    s_h=[s_h,linspace(sqrt(2)/2-1e-2,0,100)];

    p_h=d_t./s_h;   %Horizontal pitch
    
    
    %Relative maximum Nusselt number
    Nu_relExt=Nu_rel(fx,beta,p,T,Ar,eps_mf,phi_s,d_t,w_p,c_p,p_h);
    
    
    %% Alternative models: relative maximum Nusselt numbers
    % Grewal and Saxena
    Nu_relGS=1-0.21.*(p_h./d_t).^-1.75;
    Nu_relGS=Nu_relGS./Nu_relGS(:,end);
    
    
    %Gelperin and Einstein
    Nu_relGE=(1-(d_t./p_h).*(1+d_t./(d_t+p_h))).^0.25;
    Nu_relGE(Nu_relGE~=real(Nu_relGE))=NaN;
    Nu_relGE=Nu_relGE./Nu_relGE(:,end);
    

    %Ensure Nu_relGE(s_h=sqrt(2)/2)=0
    idx=isnan(Nu_relGE);
    x_GE=s_h;
    x_GE(idx)=NaN;
    
    idx=find(idx,1,"last");
    l=20;
    x_GE(idx-l:idx)=linspace(sqrt(2)/2,x_GE(idx+1),l+1);
    Nu_relGE(idx-l:idx)=linspace(0,Nu_relGE(idx+1),l+1);
    
    
    %Natusch
    Nu_relNatusch=(1-d_t./p_h).^0.25;
    Nu_relNatusch=Nu_relNatusch./Nu_relNatusch(:,end);
    
    
    
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
    plot(ax,s_h,Nu_relExt);
    plot(ax,s_h,[Nu_relGS;Nu_relNatusch],'LineStyle','--');
    plot(ax,x_GE,Nu_relGE,'LineStyle','--');
    
    hold(ax,'off');
    

    %Configure and print figure: small for repository, regular for
    %manuscript
    if small
        %Turn off axes visibility
        ax.Visible='off';


        %Export figure for repository
        fname=[dirFigs,filesep,'pi9_s',num2str(figidx)];
        
        t.Units='centimeters';
        t.InnerPosition=[0.5,0.5,4,4];
        
        fig.Units=t.Units;
        fig.Position(3:4)=t.OuterPosition(3:4)+0.5;
        
        exportgraphics(fig,[fname,'.tiff'],'Resolution',600);


    else
        %Set legend
        lgd=legend(ax,[compose('Extended Model, Ar = 1e%d',log10(Ar'));...
            {'Grewal & Saxena';...
                'Natusch et al.'};...
                'Gel''perin et al.'],...
            'Location','southwest');
        

        %Set axes labels
        xlabel(ax,figaux.subsz('\pi_9 = d_t / p_h (-)',6));
        ylabel(ax,figaux.subsz('Nu_{max} / Nu_{max} (p_h \rightarrow \infty)',6));
        
        
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




