classdef Wavelet3
    %WAVELET Simple 2D/3D Convenience class for implementing orthogonal
    %discrete wavelet transform. Inspired by Miki Lustig's @Wavelet class
    %from the sparseMRI package.
    %   Detailed explanation goes here
    
    properties
        isTrans = false;
        Nx;
        Ny;
        Nz;
        mask;
        S;
        Level;
    end
    
    methods
        function obj = Wavelet3(Nx,Ny,Nz,mask,Level)
            %WAVELET Construct an instance of this class
            %   Detailed explanation goes here
            
            obj.Nx = Nx;
            obj.Ny = Ny;
            obj.Nz = Nz;
            
            if nargin > 3
                mask = reshape(mask,obj.Nx,obj.Ny,obj.Nz);
                obj.mask = logical(mask);
            else
                obj.mask = true(obj.Nx,obj.Ny,obj.Nz);
            end
            if nargin > 4
                obj.Level = Level;
            else
                obj.Level = 4;
            end
            
            % Dummy wavelet transform to get S matrix for studpid matlab;
            [obj.S] = wavedec3(zeros(Nx,Ny,Nz),obj.Level,'haar');
            
        end
        
        function y = mtimes(obj,x)
            %mtimes Implements mtimes y = A*x operation
            %   Detailed explanation goes here
            if obj.isTrans % Adjoint Transformation
                img = waverec3(rebuildWavelet(obj,real(x))) + 1j*waverec3(rebuildWavelet(obj,imag(x)));
                y = col(img(obj.mask));
            else % Forward Transformation
                img = embed(x,logical(obj.mask));
                
                wd =  obj.collapseWavelet(wavedec3(real(img),obj.Level,'haar')) + 1j*obj.collapseWavelet(wavedec3(imag(img),obj.Level,'haar'));
                y = col(wd);
            end
            
        end
        
        function obj = ctranspose(obj)
            obj.isTrans = ~obj.isTrans;
        end
        
        function y = collapseWavelet(obj,x)
            %collapseWavelet Collapses cell arrays into a single vector for
            %output
            
            y = [];
            
            for ii = 1:length(obj.S.dec)
                y = vertcat(y, col(x.dec{ii}));
            end
        end
        
        function y = rebuildWavelet(obj,x)
            y = obj.S;
            ii = 0;
            for jj = 1:length(obj.S.dec)
                for ll = 1:length(col(obj.S.dec{jj}))
                    ii = ii + 1;
                    y.dec{jj}(ll) = x(ii);
                end
            end
                    
            
        end
            
        
        % DO we need a size function? - AMC 11/11/2017
    end
end

