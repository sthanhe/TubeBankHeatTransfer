%% Auxiliary functions for figures
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
% This class provides some auxiliary functions for plotting figures
%
%
%Requires all files packaged in the class folder and on the MATLAB path
%
%Required products, version 24.1:
%   - MATLAB


classdef figaux
    methods(Static)
        function ar=arrow(ax,x,y,type,varargin)
            %Adds an arrow to specified data coordinates
            %
            %ax: axis to add the arrow to
            %x: x-coordinates, [start,finish]
            %y: y-coordinates, [start,finish]
            %type: use types other than arrow (default), see "annotation"
            %documentation for possible types
            %varargin: property, value pairs to configure the arrow. See 
            %"annotation" documentation for possible values


            %Type=arrow if not specified otherwise
            if nargin<4
                type='arrow';
            end
        

            %Use dummy text for position indicator
            t=text(x,y,'',...
                'HorizontalAlignment','left',...
                'VerticalAlignment','middle');
            set(t,'Units','pixels');
        
        
            %Get figure object, count number of parents
            fig=ax.Parent;
            nParents=1;
            while ~isa(fig,'matlab.ui.Figure')
                fig=fig.Parent;
                nParents=nParents+1;
            end
        
            
            %Add object tree to cell array
            obj=cell(nParents+1,1);
            obj{1}=ax;
            for i=1:nParents
                obj{i+1}=obj{i}.Parent;
            end
        
        
            %Record object units, change to absolute values (pixels)
            Units=cellfun(@(x) x.Units,obj,'UniformOutput',false);
            for i=1:length(obj)
                obj{i}.Units='pixels';
            end
        
        
            %Get absolute position of axes inside figure
            pos=ax.Position;
            pos(1)=pos(1)-2;    %Bias: two pixels
            pos(2)=pos(2)-2;    %Bias: two pixels
        
            tPos={t.Position};
            tPos=vertcat(tPos{:});
        
            x=tPos(:,1)+pos(1);
            y=tPos(:,2)+pos(2);
        
        
            %x and y normalized
            x=x./fig.Position(3);
            y=y./fig.Position(4);
        
            
            %Reset units
            for i=1:length(obj)
                obj{i}.Units=Units{i};
            end
        
        
            %Draw arrow
            ar=annotation(fig,type,x,y,'Units',fig.Units,varargin{:});
        end


        function txt=subsz(txt,sz)
            %Set font size of all subscripts in txt to specified size sz
            % 
            %txt: char array or cell array of character vectors. Subscripts
            %are marked according to tex
            % 
            %sz: font size (points)

            if ischar(txt)
                txt={txt};
            end
        
            sz=num2str(sz);
            
            for i=1:numel(txt)
                %Wrap single super- and subscript in curly braces
                txt{i}=regexprep(txt{i},'(\^|_)(?!\{)(.)','$1{$2}');
            
                %Add fontsize
                txt{i}=regexprep(txt{i},'(\^|_)\{(.*?)\}',['$1{\\fontsize{',sz,'}$2}']);
            end
        end
    end
end




