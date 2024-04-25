%% Auxiliary Functions for Implicit Expansion
%GNU General Public License v3.0
%By Stefan Thanheiser: https://orcid.org/0000-0003-2765-1156
%
%Part of the paper:
%
%Thanheiser, S.; Haider, M.
%Particle Mass Diffusion Model for Level Control of Bubbling Fluidized Beds
%with Horizontal Particle Flow
%Powder Technology 2023
%
%All required files for this class can be found in the software
%repository:
%https://doi.org/10.5281/zenodo.xxxxxxx
%
%
%
%This class contains functions to assist in the use of implicit expansion.
%
%
%Requires all files packaged in the class folder and on the MATLAB path
%
%Required products:
%   - MATLAB, version 9.14


classdef implExp
    methods(Static)
        function sz=size(varargin)
            %Calculates the size of the resulting array when the inputs in
            %varargin are implicitly expanded

            dims=max(cellfun(@(x) ndims(x),varargin));
            szCell=cellfun(@(x) implExp.padsz(dims,x),varargin,'UniformOutput',false);
            szs=cell2mat(szCell');

            %all dimensions not equal to 1 must be the same
            not1=szs~=1;
            check=arrayfun(@(x) all(szs(not1(:,x),x)==szs(find(not1(:,x),1),x)),1:size(szs,2));

            if all(check)
                sz=max(szs,[],1);
                sz(any(szs==0,1))=0;
            else
                throw(MException('implExp:incompatibleArrays','Arrays have incompatible sizes'));
            end
        end


        function varargout=normalize(sz,varargin)
            %Converts all arrays in varargin to a line vectors when the
            %arrays are expected to have the size sz after implicit
            %expansion

            n=prod(sz);
            if n~=0
                pad=@(x) implExp.padsz(numel(sz),x);
                varargout=cellfun(@(x) reshape(repmat(x,sz./pad(x)),1,n),varargin,'UniformOutput',false);
            else
                varargout=repmat({NaN(sz)},1,length(varargin));
            end
        end
    end
    
    
    methods(Static, Access=private)
        function padsz=padsz(dims,x)
            %Makes size vectors x have the correct length based on the 
            %number of required dimensions dims
            padsz=[size(x),ones(1,dims-ndims(x))];
        end
    end
end




