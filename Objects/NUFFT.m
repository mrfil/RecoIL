classdef NUFFT
    %NUFFT  - Convenience class for creating an NUFFT transform object
    % Replaces Brad Sutton's Gnufft and Gnufft_v2 classes.
    % Syntax:  G = NUFFT(kx,ky,kz,Nx,Ny,Nz,varargin)
    %
    % Inputs:
    %    kx      - N-D Array containing kspace in unitless [-(Nx-1)/2,Nx/2]
    %              dimensions
    %    ky      - N-D Array containing kspace in unitless [-(Ny-1)/2,Ny/2]
    %              dimensions
    %    kz      - N-D Array containing kspace in unitless [-(Nz-1)/2,Nz/2]
    %              dimensions. May be array of zeros in 2D recons.
    %    Nx      - Scalar containing size of reconstructed image in x
    %    Ny      - Scalar containing size of reconstructed image in y
    %    Nz      - Scalar containing size of reconstructed image in z, for
    %              a 2D reconstruction or operation, this must be 1.
    %   
    % Optiona Name-Value Pair Inputs:
    %    'mask'  - Image domain support mask of logicals to restrict the
    %              transform only to be computed at the true points.
    %    'VoxelBasis' - 'delta'  (Default) Compute the transform assuming
    %                            each voxel is represented in image space 
    %                            by a delta function.
    %                   'boxcar' Compute the transform by assuming each
    %                            voxel is represented in image space by a
    %                            boxcar function, requiring a sinc
    %                            weighting to be applied in data/kspace
    %                            domain
    %    'InterpMethod' - 'sparse' (Default) Computing the NUFFT using a
    %                               sparse matrix-vector multiplication 
    %                               for the spatial interpolator step. This
    %                               sparse matrix can be quite large, 
    %                               requiring very high amounts of memory 
    %                               for large problems. This method is
    %                               multithreaded for the spatial
    %                               interpolation step, which is usually
    %                               the most computationally costly step in
    %                               the computation of the NUFFT. Use
    %                               'table' if you run out of memory when
    %                               initalizing the NUFFT object.
    %                   - 'table'   This method uses a compiled mex file to
    %                               perform the spatial interpolation step.
    %                               It requires very little memory even for
    %                               very large problems, but is single
    %                               threaded which can be slow for large
    %                               problems.
    %                        
    % Outputs:
    %    G     - NUFFT class object with overloaded mtimes operator for
    %            computing the forward and adjoint NUFFT operations through
    %            A*x and A'*x.
    %
    % Example:
    %   TBD - AMC
    %
    %
    % Other m-files required: IRT nufft implementation and private directories
    %   TimeSegmentation.m, NUFFT.m
    %   
    % Subfunctions: none
    % MAT-files required: none
    %
    % Author: Alex Cerjanic (1.0)
    % University of Illinois at Urbana-Champaign
    % email address:
    % Website:
    % 4-Sep-2017; Last revision: 4-Sep-2017
    
    properties
        Kx;
        Ky;
        Kz;
        Nx;
        Ny;
        Nz;
        st;
        dims;
        isTranspose = false;
        isMasked = false;
        mask = [];
        SincWeighting = false;
        VoxelBasisWeights = 1;
        InterpMethod = '';
    end
    
    methods
        function obj = NUFFT(kx,ky,kz,Nx,Ny,Nz,varargin)
            
            J = 5;
            Jz = 5;
            K = 2;
            %% Deal with Inputs
            obj.Kx = kx;
            obj.Ky = ky;
            obj.Kz = kz;
            
            obj.Nx = Nx;
            obj.Ny = Ny;
            obj.Nz = Nz;
            InterpMethodString = 'sparse';
            VoxelBasis = 'delta';
            mask = [];
            p = inputParser();
            
            addOptional(p,'InterpMethod',InterpMethodString);
            addOptional(p,'mask',mask);
            addOptional(p,'VoxelBasis',VoxelBasis);
            addOptional(p,'GridOversampling',K);
            parse(p,varargin{:});
            
            obj.mask = p.Results.mask;
            obj.InterpMethod = p.Results.InterpMethod;
            VoxelBasis = p.Results.VoxelBasis;
            K = p.Results.GridOversampling;
            if ~isempty(obj.mask)
                obj.isMasked = true;
            end
            
            %% Prepare the NUFFT - initialize the args cell for nufft_init_v2
            if(Nz > 1)
                fprintf('Performing 3D NUFFT...\n');
                if strcmp(obj.InterpMethod,'table')
                    args = {[2*pi/obj.Ny*obj.Ky(:),2*pi/obj.Nx*obj.Kx(:), 2*pi/obj.Nz*obj.Kz(:)],[obj.Ny, obj.Nx, obj.Nz],[J,J,Jz],[K*obj.Ny,K*obj.Nx,K*obj.Nz],[obj.Ny/2, obj.Nx/2, obj.Nz/2], 'table', 2^10, 'minmax:kb'};
                    fprintf('Using table lookup...\n');
                else
                    args = {[2*pi/obj.Ny*obj.Ky(:),2*pi/Nx*obj.Kx(:), 2*pi/obj.Nz*obj.Kz(:)],[obj.Ny, obj.Nx, obj.Nz],[J,J,Jz],[K*obj.Ny,K*obj.Nx,K*obj.Nz],[obj.Ny/2, obj.Nx/2, obj.Nz/2], 'minmax:kb', 2^10};
                    fprintf('Using sparse interpolator...\n');
                    % Estimate the size of sparse interpolator matrix
                    nnz = J*J*Jz*size(obj.Kx(:),1);
                    n = prod([K*obj.Ny,K*obj.Nx,K*obj.Nz]);
                    sparseMatrixSize = 2*(max(nnz,1)*(8+8) + (n+1) * 8);
                    sparseMatrixSizeGiB = sparseMatrixSize/(2^30);
                    fprintf('Minimum size of sparse matrix interpolator = %.1f GiB\n',sparseMatrixSizeGiB)
                    if sparseMatrixSizeGiB > get_free_mem/2
                        fprintf('Sparse Matrix is too large for system memory. Switching automatically to table method!');
                        args = {[2*pi/obj.Ny*obj.Ky(:),2*pi/obj.Nx*obj.Kx(:), 2*pi/obj.Nz*obj.Kz(:)],[obj.Ny, obj.Nx, obj.Nz],[J,J,Jz],[K*obj.Ny,K*obj.Nx,K*obj.Nz],[obj.Ny/2, obj.Nx/2, obj.Nz/2], 'table', 2^10, 'minmax:kb'};
                        fprintf('Using table lookup...\n');
                    end
                    
                end
                
            else
                if strcmp(obj.InterpMethod,'table')
                    args = {[2*pi/obj.Ny*obj.Ky(:),2*pi/obj.Nx*obj.Kx(:)],[obj.Ny,obj.Nx],[J,J],[K*obj.Ny,K*obj.Nx],[obj.Ny/2,obj.Nx/2]*(1), 'table', 2^10, 'minmax:kb'};
                    fprintf('Using table lookup...\n');
                else
                    args = {[2*pi/obj.Ny*obj.Ky(:),2*pi/obj.Nx*obj.Kx(:)],[obj.Ny,obj.Nx],[J,J],[K*obj.Ny,K*obj.Ny],[obj.Ny/2,obj.Nx/2]*(1), 'minmax:kb', 2^10};
                    fprintf('Using sparse interpolator...\n');
                    % Estimate the size of sparse interpolator matrix
                    nnz = J*J*size(obj.Kx(:),1);
                    n = prod([K*obj.Ny,K*obj.Nx]);
                    sparseMatrixSize = 2*(max(nnz,1)*(8+8) + (n+1) * 8);
                    sparseMatrixSizeGiB = sparseMatrixSize/(2^30);
                    fprintf('Minimum size of sparse matrix interpolator = %.1f GiB\n',sparseMatrixSizeGiB)
                    if sparseMatrixSizeGiB > get_free_mem/2
                        fprintf('Sparse Matrix is too large for system memory. Switching automatically to table method!');
                        args = {[2*pi/obj.Ny*obj.Ky(:),2*pi/obj.Nx*obj.Kx(:)],[obj.Ny,obj.Nx],[J,J],[K*obj.Ny,K*obj.Nx],[obj.Ny/2,obj.Nx/2]*(1), 'table', 2^10, 'minmax:kb'};
                        fprintf('Using table lookup...\n');
                    end
                end
            end
            
            %% Setup voxel basis function
            switch VoxelBasis
                case 'delta'
                    obj.SincWeighting = false;
                    obj.VoxelBasisWeights = ones(size(obj.Kx(:)));
                case 'boxcar'
                    obj.SincWeighting = true;
                    obj.VoxelBasisWeights = sinc(obj.Kx(:)/obj.Nx).*sinc(obj.Ky(:)/Ny).*sinc(obj.Kz(:)/Nz);
                otherwise
                    error('Unrecognized Voxel Basis function type!');
            end
            
            %% Initialize transform
            
            obj.st = nufft_init_v2(args{:});
            if obj.isMasked
                obj.dims = [size(obj.Kx(:),1), sum(obj.mask(:))];
            else
                obj.dims = [size(obj.Kx(:),1), obj.Nx*obj.Ny*obj.Nz];
            end
        end
        %% Overloaded Operators
        function y = mtimes(obj, x)
            %tic
            % Deal with forward transform first
            if ~obj.isTranspose
                % Apply mask to input data
                if obj.isMasked
                    x = embed(x, obj.mask);
                end
                % Reshape the data to be of the right size since it comes
                % in as a column vector
                x = reshape( x, obj.Ny,obj.Nx,obj.Nz);
                % Apply the transform using the interpolator matrix or
                % table operation
                y = nufft(double(x), obj.st);
                if (obj.SincWeighting == true)
                    y = y.*obj.VoxelBasisWeights;
                end
                
            else % Deal with the adjoint case
                % Apply the adjoint transform
                if (obj.SincWeighting == true)
                    x = x.*obj.VoxelBasisWeights;
                end
                y = nufft_adj(double(x), obj.st);
                if obj.isMasked
                    % Return only masked data
                    y = y(obj.mask);
                end
            end
            % Ensure we return a column vector
            y = y(:);
            %toc
        end
        
        function obj = ctranspose(obj)
            % Handle the transpose operator for adjoint operator
            obj.isTranspose = ~obj.isTranspose;
            obj.dims = fliplr(obj.dims);
        end
        
        function y = size(obj)
            % Return the size of the object
            y = obj.dims;
        end
        
        
        
    end
    
end

