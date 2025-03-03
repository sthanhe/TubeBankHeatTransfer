%% Dimensional Analysis
%GNU General Public License v3.0
%By Stefan Thanheiser: https://orcid.org/0000-0003-2765-1156
%
%Part of the paper:
%
%Thanheiser, S.; Haider, M.
%Dispersion Model for Level Control of Bubbling Fluidized Beds with 
%Particle Cross-Flow
%Chemical Engineering Science 2024
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
%This script sets up the dimensional set shown in the paper.
%
%
%Required products, version 24.1:
%   - MATLAB
%Necessary files, classes, functions, and scripts:
%   - None


%% Influencing factors
%Original influencing factors from Molerus
h=struct('name','h','unit','W/m²K');
g=struct('name','g','unit','m/s²');
rho_p_rho_g=struct('name','rho_p_rho_g','unit','kg/m³');
c_p=struct('name','c_p','unit','J/kgK');
my_g=struct('name','my_g','unit','Pas');
rho_g=struct('name','rho_g','unit','kg/m³');
c_g=struct('name','c_g','unit','J/kgK');
k_g=struct('name','k_g','unit','W/mK');
w_mf=struct('name','w_mf','unit','m/s');
eps_mf=struct('name','eps_mf','unit','');
w_e=struct('name','w_e','unit','m/s');


%Own (added) influencing factors
s_h=struct('name','s_h','unit','');
p_h=struct('name','p_h','unit','m');
d_t=struct('name','d_t','unit','m');
w_p=struct('name','w_p','unit','m/s');


%% Maximum heat transfer in laminar regime
% results in same dimensionless variables as in Molerus

[A,B,C,D]=dimMat(h,c_p,...
            k_g,rho_p_rho_g,g,my_g);

setLam=dimSet(A,B,C,D);


%% Maximum HTC in laminar regime including tube diameter
% results in same dimensionless variables as in Molerus, plus d_t/l_l

[A,B,C,D]=dimMat(h,c_p,d_t,...
            k_g,rho_p_rho_g,g,my_g);

setLamDt=dimSet(A,B,C,D);


%% Maximum HTC in turbulent regime
% Molerus suggests that pi_1 and pi_3 are coupled via: pi_1*pi_3^(-1/3)
%pi_1=h*l_l/k_g
%pi_3=rho_g/(rho_p-rho_g)
%pi_1*pi_3^(-1/3)=h*l_t/k_g


[A,B,C,D]=dimMat(h,c_g,rho_g,...
            k_g,rho_p_rho_g,g,my_g);

setTurb=dimSet(A,B,C,D);


%% Molerus dimensional analysis with all factors
% pi1...pi4 and pi_7 are identical to Molerus (pi_2 is reversed)
% Molerus' pi_6=replace pi_6 with pi_5/pi_6
% Molerus' pi_5=replace pi_5 with pi_5*pi_2^(1/3)
% or: small variations of the D-matrix (=linear combinations of pi-factors)


[A,B,C,D]=dimMat(h,c_p,c_g,rho_g,w_e,w_mf,eps_mf,...
            k_g,rho_p_rho_g,g,my_g);


D{'pi2','c_p'}=-1;
D{'pi5','c_p'}=1/3;
D{'pi6','w_e'}=1;
D{'pi6','w_mf'}=-1;


setMolerus=dimSet(A,B,C,D);


%% Mixed regime including particle cross-flow
%3 dimensionless velocities: w_e (pi_5), w_mf (pi_6) and w_p (pi_8)
%All include the same factor, relating particle convection to conduction
%Molerus' pi_6=pi_5/pi_6


[A,B,C,D]=dimMat(h,c_p,c_g,rho_g,w_e,w_mf,eps_mf,w_p,...
            k_g,rho_p_rho_g,g,my_g);


D{'pi2','c_p'}=-1;
D{'pi5','c_p'}=1/3;
D{'pi6','c_p'}=1/3;
% D{'pi6','w_e'}=1;
% D{'pi6','w_mf'}=-1;
D{'pi8','c_p'}=1/3;
% D{'pi8','w_mf'}=1/3;


setMolerusWp=dimSet(A,B,C,D);


%% New dimensional analysis with all new factors


[A,B,C,D]=dimMat(h,c_p,c_g,rho_g,w_e,w_mf,eps_mf,d_t,p_h,w_p,...
            k_g,rho_p_rho_g,g,my_g);


D{'pi2','c_p'}=-1;
D{'pi5','c_p'}=1/3;
D{'pi6','c_p'}=1/3;
D{'pi9','d_t'}=1;
D{'pi9','p_h'}=-1;
D{'pi10','c_p'}=1/3;


setNew=dimSet(A,B,C,D);


%% Auxiliary functions
%This function creates the matrices that make up the dimensional set
function [A,B,C,D]=dimMat(varargin)
    %Auxiliary functions to count dimensions from the created unit string
    count=@(str,dim) numel(regexp(str,['*?[^/]',dim]))-...
                    numel(regexp(str,['/',dim]));

    dimsum=@(nomden,dim) count(nomden{1},dim)-count(nomden{2},dim);


    %Get variables
    vars=[varargin{:}];
    n=numel(vars);

    
    %Set up dimensional matrix
    dimnames={'m','s','kg','K'};
    mat=table('Size',[length(dimnames),n],...
                'VariableTypes',repmat({'double'},1,n),...
                'VariableNames',{vars.name},...
                'RowNames',dimnames);


    %Extract dimensions from every unit
    for i=1:n
        %Normalize exponents
        vars(i).unit=strrep(vars(i).unit,'²','^2');
        vars(i).unit=strrep(vars(i).unit,'³','^3');


        %Replace derived SI units with basic SI units
        vars(i).unit=strrep(vars(i).unit,'Pa','N*m^-2');
        vars(i).unit=strrep(vars(i).unit,'W','J*s^-1');

        vars(i).unit=strrep(vars(i).unit,'J','N*m');

        vars(i).unit=strrep(vars(i).unit,'N','kg*m*s^-2');


        %Split into nominator and denominator
        nomden=strsplit(vars(i).unit,'/');
        if isscalar(nomden)
            nomden=[nomden,{''}]; %#ok<AGROW>
        end


        %Replace exponents with string of dimensions, normalize result
        for j=1:length(nomden)
            %Positive exponents
            nomden{j}=regexprep(nomden{j},...
                        '(\w+)\^(\d+)',...
                        '${repmat([$1,''*''],1,str2num($2))}');

            %Negative exponents
            nomden{j}=regexprep(nomden{j},...
                        '(\w+)\^-(\d+)',...
                        '${repmat([''/'',$1],1,str2num($2))}');

            %Cleanup and normalize
            nomden{j}=strrep(nomden{j},'*/','/');
            nomden{j}=['*',nomden{j}];
        end


        %Count units
        for j=1:length(dimnames)
            mat{dimnames{j},i}=dimsum(nomden,dimnames{j});
        end
    end
    

    %Remove empty dimensions
    mat(all(mat{:,:}==0,2),:)=[];


    %Constituting numbers
    nDims=height(mat);  %Number of dimensions
    nB=n-nDims;         %Number of variables in Matrix B
    % nA=nDims;         %Number of variables in Matrix A
    % nP=nB;            %Number of pi-factors
    
    
    %Matrices
    A=mat(:,nB+1:end);
    B=mat(:,1:nB);

    pinames=compose('pi%d',1:nB);
    C=table('Size',[nB,nDims],...
            'VariableTypes',repmat({'double'},1,nDims),...
            'VariableNames',A.Properties.VariableNames,...
            'RowNames',pinames);
    
    D=table('Size',[nB,nB],...
            'VariableTypes',repmat({'double'},1,nB),...
            'VariableNames',B.Properties.VariableNames,...
            'RowNames',pinames);
    
    D{:,:}=eye(nB);
end


%This function calculates the C matrix and creates the dimensional set
function set=dimSet(A,B,C,D)
    %Fundamental equation
    Cmat=-D{:,:}*(A{:,:}^-1*B{:,:})';
    
    
    %Fix rounding issues
    Cmat(abs(Cmat)<1e-6)=0;
    C{:,:}=Cmat;
    
    
    %Dimensional set
    set=[B,A;D,C];
end