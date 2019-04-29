classdef Gfm
    %Gfm Field Map Object for regularized least square estimation of the
    %field map
    %   Detailed explanation goes here
    
    properties
        mask = [];
        TEs = [];
        Nx;
        Ny;
        Nz;
        isTranspose = false;
    end
    
    methods
        function obj = Gfm(Nx,Ny,Nz,TEs,mask)
            %Gfm Construct an instance of this class
            %   Detailed explanation goes here
            obj.Nx = Nx;
            obj.Ny = Ny;
            obj.Nz = Nz;
            obj.TEs = TEs;
            obj.mask = logical(mask);
            
        end
        
        function y = mtimes(obj,x)
            %mtimes Matrix times operation for Gfm
            %   Detailed explanation goes here
            fmEst = reshape(x,[],length(obj.TEs));
            y = zeros(size(fmEst));
            if ~obj.isTranspose % Forward Operation
                
                for ii = 1:length(obj.TEs)
                    y(:,ii) = fmEst.*obj.TEs(ii);
                end
            else % Adjoint Operation
                for ii = 1:length(obj.TEs)
                    y(:,ii) = fmEst.*obj.TEs(ii);
                end
                y = sum(y,2);
            end
            
        end
        
        function obj = ctranspose(obj)
            obj.isTranspose = ~obj.isTranspose;
        end
    end
end

