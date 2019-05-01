classdef Gdft
   %Gdft Discrete Fourier Transform function for MRI
   %   MATLAB only function for Gdft
   
   properties
      transpose = false;
      kx;
      ky;
      kz;
      Nx;
      Ny;
      Nz;
      X;
      Y;
      Z;
      FM; % Field Map
      tt; % Timing Vector
      SincWeighting = false;
      VoxelBasisWeights = 1;
   end
   
   methods
      
      function obj = Gdft(kx,ky,kz,Nx,Ny,Nz,FM,tt,varargin)
         if nargin < 8
            obj.tt = 0;
         else
            %obj.tt = col(tt);
            nShots = length(kx(:))/length(tt(:));
            obj.tt = repmat(tt(:),nShots,1);
         end
         if nargin < 7
            obj.FM = 0;
         else
            obj.FM = col(rot90(FM));
         end
         
         VoxelBasis = 'delta';
         p = inputParser();
         
         addOptional(p,'VoxelBasis',VoxelBasis);
         parse(p,varargin{:});
         
         VoxelBasis = p.Results.VoxelBasis;
         
         % Generate coefficient locations
         [obj.X,obj.Y,obj.Z] = ndgrid(-(Nx)/2:(Nx-1)/2,-(Ny)/2:(Ny-1)/2,-(Nz)/2:(Nz-1)/2);
         obj.X = col(obj.X)/Nx;
         obj.Y = col(obj.Y)/Ny;
         obj.Z = col(obj.Z)/Nz;
         
         obj.kx = kx;
         obj.ky = ky;
         obj.kz = kz;
         obj.Nx = Nx;
         obj.Ny = Ny;
         obj.Nz = Nz;
         
         %% Setup voxel basis function
         switch VoxelBasis
            case 'delta'
               obj.SincWeighting = false;
               obj.VoxelBasisWeights = ones(size(obj.kx(:)));
            case 'boxcar'
               obj.SincWeighting = true;
               obj.VoxelBasisWeights = sinc(obj.kx(:)/obj.Nx).*sinc(obj.ky(:)/Ny).*sinc(obj.kz(:)/Nz);
            otherwise
               error('Unrecognized Voxel Basis function type!');
         end
      end
      
      function out = mtimes(obj,b)
         if ~obj.transpose % Forward transform
            out = zeros(size(obj.kx));
            for ii = 1:length(obj.kx)
               kxTP = obj.kx(ii)*2*pi;
               kyTP = obj.ky(ii)*2*pi;
               kzTP = obj.kz(ii)*2*pi;
               sumData = b.*exp(-1j*(kxTP*obj.X + ...
                              kyTP*obj.Y + ...
                              kzTP*obj.Z + ...
                              obj.FM.*obj.tt(ii)));
               out(ii) = sum(sumData(:));
            end
            if (obj.SincWeighting == true)
               out = out.*obj.VoxelBasisWeights;
            end
         else %Adjoint transform
            out = zeros(size(obj.X));
            if (obj.SincWeighting == true)
               b = b.*obj.VoxelBasisWeights;
            end
            for ii = 1:length(obj.X(:))
               ixTP = obj.X(ii)*2*pi;
               iyTP = obj.Y(ii)*2*pi;
               izTP = obj.Z(ii)*2*pi;
               sumData = b.*exp(1j*(ixTP*obj.kx + ...
                           iyTP*obj.ky + ... 
                           izTP*obj.kz + ...
                           obj.FM(ii).*obj.tt));
               out(ii) = sum(sumData(:));
            end
            out = 1/(obj.Nx*obj.Ny*obj.Nz)*out;
         end
      end
      
      
      function obj = ctranspose(obj)
         obj.transpose = ~obj.transpose;
      end
      
      function out = size(obj)
         out = zeros(1,2);
         out(1) = length(obj.kx(:));
         out(2) = obj.Nx*obj.Ny*obj.Nz;
      end
   end
   
end

