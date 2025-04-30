function compPi9(fx,beta,dirFigs,figidx,small)
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
    end


end




