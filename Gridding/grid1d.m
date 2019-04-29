function dato = grid1d(kxin,datin,nx,W);
%function dato = grid1d(kxin,datin,nx,W);
% kxin 1D k-space trajectory
% datin is 1D data in
% nx is size of grid, ie. 64
% W is window width, ie. 2.5
%outputs image gridded onto evenly spaced points
%
%Brad Sutton
% Biomedical Imaging Center
% University of Illinois at Urbana-Champaign
% Sept. 28, 2005
  
 % besseli(0,val)  % bessel funciton of order 0 of first kind
kbeval = inline('1/Lx*besseli(0,beta*sqrt(1-4*(t.*t)./(Lx*Lx)))','t','beta','Lx');

if ~exist('W','var')
  W = 3;
  beta = 13.9;
else
  beta = pi*1.45*(W);
end


if ~exist('nx','var')
  nx = length(kxin);
end


ww = gradient(kxin);
datin = datin.*ww;

kxo = [-nx/2:0.5:nx/2-0.5];

tmp = zeros(size(kxo));

for ii = 1:length(kxin)
   distx = kxo-kxin(ii);
   ll = find(abs(distx)<(W/2));
   tmp(ll) = tmp(ll)+ datin(ii).*kbeval(abs(distx(ll)),beta,W);
end


x=linspace(-1,1,nx*2+1);
x=x(1:2*nx); %+1/(2*nx);

sinker = sqrt(pi^2*W^2.*(x.^2)-beta^2);
cx = sin(sinker)./sinker;

cx = cx./max(cx(:));
%keyboard
dato = fftshift(ifft(fftshift(tmp)))./cx;

dato = dato([nx/2+1:3/2*nx]);

