classdef Gdft_r2
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
      R2; % R2* Map
      tt; % Timing Vector
      Gx; % X Gradient Map
      Gy; % Y Gradient Map
      Gz; % Z Gradient Map
   end
   
   methods
      
      function obj = Gdft_r2(kx,ky,kz,Nx,Ny,Nz,FM,R2,tt,Gx,Gy,Gz)
         if nargin < 9
            obj.tt = 0;
            obj.Gx = 0;
            obj.Gy = 0;
            obj.Gz = 0;
            
         else
            %obj.tt = col(tt);
            nShots = length(kx(:))/length(tt(:));
            obj.tt = repmat(tt(:),nShots,1);
            obj.Gx = col(Gx);
            obj.Gy = col(Gy);
            obj.Gz = col(Gz);
         end
         if nargin < 7
            obj.FM = 0;
            obj.R2 = 0;
         else
            obj.FM = col(FM);
            obj.R2 = col(R2);
         end
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
         
         
         
      end
      
      function out = mtimes(obj,b)
         if ~obj.transpose % Forward transform
            out = zeros(size(obj.kx));
            for ii = 1:length(obj.kx)
               kxTP = obj.kx(ii);
               kyTP = obj.ky(ii);
               kzTP = obj.kz(ii);
               bfunc = sinc(kxTP./obj.Nx+obj.Gx.*obj.tt(ii)).* ...
                  sinc(kyTP./obj.Ny+obj.Gy.*obj.tt(ii)).* ...
                  sinc(kzTP./obj.Nz+obj.Gz.*obj.tt(ii));
               sumData = b.*bfunc.*exp(-obj.R2.*obj.tt(ii)).* ...
                  exp(-1j*(2*pi*(kxTP*obj.X+kyTP*obj.Y+kzTP*obj.Z ) + ...
                  obj.FM.*obj.tt(ii)));
               out(ii) = sum(sumData(:));
            end
         else %Adjoint transform
            out = zeros(size(obj.X));
            for ii = 1:length(obj.X(:))
               ixTP = obj.X(ii);
               iyTP = obj.Y(ii);
               izTP = obj.Z(ii);
               bfunc = sinc(obj.kx./obj.Nx+obj.Gx(ii)*obj.tt).* ...
                  sinc(obj.ky./obj.Ny+obj.Gy(ii)*obj.tt).* ...
                  sinc(obj.kz./obj.Nz+obj.Gz(ii)*obj.tt);
               sumData = b.*bfunc.*exp(-obj.R2(ii).*obj.tt).* ...
                  exp(1j*(2*pi*(ixTP*obj.kx+iyTP*obj.ky+izTP*obj.kz) + ...
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
