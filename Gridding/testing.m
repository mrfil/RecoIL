path(path,'/net/fincher/home/bpsutton/Matlab/')
path(path,'/net/fincher/home/bpsutton/Matlab/FFT/')
path(path,'/net/fincher/home/bpsutton/Cfiles/')
%CV input area
FOV = 20;   %in cm
N = 64;     %reconstruction size
gamp = 2.2;
gslew = 180;
nl = 1;  %number of interleaves
TE = 20e-03;  %TE in seconds
tau = 2e-03;  %dTE between acquisitions for field map

if N == 128
     datseg = 14332; %points collected per spiral per coil
else 
     datseg = 4026;
end

%OBJECT FOR SIMULATION
% object to be used for simulation is close to phantom object - similar 
%    features
if 0   %Simple phantom-like object
  obj = masking(ones(N,N),23);
  obj(12:17,30:34) = 0;
  obj(47:52,30:34) = 0;
  obj(18:46,30:34) = 2;
  we = zeros(N*N,1);
end

if 1
  A = [1;-0;-0];
  a = [0.6;0.1;0.1];
  b = [0.6;0.4;0.4];
  x0 = [0;-0.27;0.27];
  y0 = [0;0;0];
  phi = [0;0;0];
obj = phantom([A a b x0 y0 phi],N);

A = [1;1];
a = [0.15;0.15];
b = [0.5;0.5];
x0 = [-0.27;0.27];
y0 = [0;0];
phi = [0;0];
we = 2*pi*0*phantom([A a b x0 y0 phi],N);
end


if 1
  Tmax = 5*16384*1e-6;   %5*16384*1e-6;
  %dts = 5e-6;    %4e-6
  [kx,ky] = genkspace(FOV,N,Tmax,nl*datseg,nl,gamp,gslew);
if nl== 2
  ndat = length(kx)/2;    %num data pts for one spiral
  kx1 = kx(1:ndat);       %kx for first spiral
  kx2 = kx((ndat+1):end); %kx for second interleave
  ky1 = ky(1:ndat);       %ky for first interleave
  ky2 = ky((ndat+1):end); %ky for second interleave
  datseg = ndat;          %set num data pts to one spiral
end
end


%TIMING
%********************************************************
% Set up timing of slice acquisition
tt = [(-datseg/2+1:datseg/2)*4e-06+TE]'; %Timing of non-delayed acquisition

tt_ext = [];
for jj = 1:nl
  tt_ext = [tt_ext;tt];
end


%SPACE COORDS
%********************************************************
% Set up coordinates of space domain
outsz = N;
[x,y] = meshgrid([-outsz/2:outsz/2-1]./outsz);
xvals = x(:);
yvals = y(:);


% FIELDMAP
%********************************************************
% Set up initial fieldmap
%we = zeros(size(x));

%load wefull
%we = wefull;

%MASK
%********************************************************
% Determine mask, to specify which values can be nonzero in image.
% We mask the fieldmap to make it zero outside the object
% Radius for mask determined by looking at field map

if 1
  mask = masking(ones(N,N),60);
elseif 0    %USING coil1 and coil2.mat
   load iterim1
   load iterim2
  DD = abs((iterim1.^2+iterim2.^2).^(1/2));
  SE = [1 1 1;
   1 1 1];
  mask = dilate(DD>2,SE,4);%erode(dilate(DD>1,SE,3),SE,3);
  ll = find(mask==1);
  mask = embed(ones(length(ll),1),mask);
  %  mask = (DD>2);
else
 mask = ones(outsz,outsz);
end


% Now mask out values of xval, yval, we
l = find(mask(:)>0);
xval = xvals(l);
yval = yvals(l);
we = we(:);
we = we(l);
npm = length(xval);   % gives number of points masked

%COIL SENSITIVITY
%******************************************************
% Get sensitivity data for coil
if 1
  numcoil = 1;
  s = ones(npm,numcoil);
end


if 1
   A2 = mri(kx, ky, tt_ext, we, xval, yval, s, FOV, 64);
   A1 = mri(kx, ky, tt_ext+tau, we, xval, yval, s, FOV, 64);
end
if 0
   A2 = mri_rect(kx, ky, tt_ext, we, xval, yval, s, FOV, 64);
   A1 = mri_rect(kx, ky, tt_ext+tau, we, xval, yval, s, FOV, 64);
end
if 1 
  objm = obj(l);
  sigma = 0;  %.0000000000000000001;
  dat2 = A2*objm+sigma/2*(randn(nl*datseg,1)+i*randn(nl*datseg,1));
  dat1 = A1*objm+sigma/2*(randn(nl*datseg,1)+i*randn(nl*datseg,1));
end

if 0
  dat1 = [s1t1r(1:datseg);s1t1r((2*datseg+1):(3*datseg))];
  dat2 = [s1t2r(1:datseg);s1t2r((2*datseg+1):(3*datseg))];
end

ww = weight_vor(kx,ky,nl);

result2 = k2image(kx,ky,dat1,ww,64);

imagesc(abs(result2))
