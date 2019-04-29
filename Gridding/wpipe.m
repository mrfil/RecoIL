% This file will determine the weighting coefficients
%   to be used with arbitrary kspace data, using the
% method of Pipe:99, for regridding & conj. phase

addpath('/afs/engin.umich.edu/u/b/p/bpsutton/Private/Data/');

load kspace;
kx = kspace(:,1);
ky = kspace(:,2);
ndat = length(kx);

% First, define a convolution filter.
% The convolution filter of choice is the Kaiser-Bessel
%  due to its FT going to zero quickly also

%Set width, D, of convlution filter:
D = max(max(kx)-min(kx), max(ky)-min(ky));
% set parameter of conv. filter
beta = 16;  % 16 is used in Pipe's paper

% Now the convolution filter will be determined by
% 1/D*(besseli(0,(beta*sqrt(1-((2/D)*u).^2))))
% where u is sqrt((kx(ii)-kx).^2+(ky(ii)-ky).^2)

% initialize variables
wc = zeros(ndat,1);  % this will be w convolved with c
epsilon = 1e-01;     % this will determine the stopping point for the iteration
iter = 0
w = ones(ndat,1);


tic
% iterate
while (abs(wc-ones(ndat,1)) > eps*ones(ndat,1))
     iter= iter+1
     conv = zeros(ndat,1);
     for ii = 1:ndat
         u = (((kx(ii)*ones(ndat,1)-kx).^2)+((ky(ii)*ones(ndat,1)-ky).^2)...
                              ).^(1/2);
         conv = 1/D*(besseli(0,(beta*sqrt(1-((2/D)*u).^2))));
         wc(ii) = (w.')*conv; 
     end;

     w = w./wc;  
end;

toc

save /afs/engin.umich.edu/u/b/p/bpsutton/Private/Data/wpipe w;

quit;

