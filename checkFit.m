%% Check results of regression graphically
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
% This function checks the results of any regression by plotting the 
% estimates against the measurements. It also provides some helpful 
% insights into the quality of the regression and possible issues.  
%
%
%Requires all auxiliary classes and functions on the MATLAB path
%
%Required products, version 24.1:
%   - MATLAB
%Necessary classes, functions, files, and scripts:
%   - None


function [fig,t,ax]=checkFit(X,y,fx,beta,name,figidx,small,dirFigs,lim)
    % Inputs:
    % X: regressor matrix, double [m,n]
    % y: response variable, double [m,1]
    % fx: model function in the form fx(b,X), where b are the regression
    %       coefficients, function handle
    % beta: regression coefficients, double
    % name: model name, char
    % figidx: index of figure window, double
    % small: indicator whether to create the small figure version, logical
    % dirFigs: path to directore where figures should be stored, char
    % lim: upper axes limit (optional)
    % 
    % Outputs:
    % fig: figure handle
    % t: tile handle
    % ax: axis handle
    % 
    % 
    % call this function right after the regression for best results


    %Model estimates: remove negative values
    yEst=fx(beta,X);

    isNeg=yEst<-1e-3;
    yEst(isNeg)=[];
    y(isNeg)=[];


    %Marker and figure size, figure name, figure title
    if small
        mkrsz=18;
        inPos=[1.5,1,4.7,4.7];
        fname=[dirFigs,filesep,'fit_',name];
    else
        mkrsz=36;
        inPos=[2,1.5,12,12];
        fname=[dirFigs,filesep,'Model_',name];


        %Display warning message from fit in title
        [~,id]=lastwarn();
        txt=['Model ',name];
        switch id
            case 'stats:nlinfit:ModelConstantWRTParam'
                txt=[txt,', Insensitive Parameters'];
            case 'stats:nlinfit:IllConditionedJacobian'
                txt=[txt,', Ill-Conditioned Jacobian'];
            case ''
        end
    

        %Display negative regression coefficients in title
        isNeg=find(beta<0);
        negTitle=compose(', \\beta_%d=%.3f',isNeg,beta(isNeg));
        txt=[txt,negTitle{:}];
    end


    %Set up figure
    if figidx==0
        fig=figure(100);
    else
        fig=figure(figidx);
    end
    clf(fig);
    t=tiledlayout(fig,1,1,'Padding','tight');
    ax=nexttile(t);
    colors=ax.ColorOrder;
    hold(ax,'on');
    

    %Plot data
    scatter(ax,y,yEst,mkrsz);
    

    %Plot equivalence lines
    if nargin<9
        lim=max([ax.XLim(2),ax.YLim(2)]);
    end
    eq=linspace(0,lim,100);
    
    plot(ax,eq,eq,'Color',colors(2,:));
    plot(ax,eq,eq.*1.2,'Color',colors(2,:),'LineStyle','--');
    plot(ax,eq,eq./1.2,'Color',colors(2,:),'LineStyle','--');
    
    hold(ax,'off');
    
    
    %Axes limits
    ax.XLim=[0,lim];
    ax.YLim=[0,lim];

    
    %Axes appearance
    if small
        ax.Visible='off';
    else
        xlabel(ax,'Measured \pi_1 (-)');
        ylabel(ax,'Predicted \pi_1 (-)');

        title(ax,txt);
    end


    %Export figure for repository
    t.Units='centimeters';
    t.InnerPosition=inPos;
    
    fig.Units=t.Units;
    fig.Position(3:4)=t.OuterPosition(3:4)+0.5;
    
    exportgraphics(fig,[fname,'.tiff'],'Resolution',600);
end




