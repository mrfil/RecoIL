addpath('/home/bsutton/Matlab/EPIRaw/')

BW = 2894;
N = 64;
Echospus = 500;
fov = 220;

kvecr = grdcalcdti(BW,N,Echospus,fov).';

obj = zeros(64,1);
obj(20:30) = 1;

xx = linspace(-0.5,0.5,N+1);
xx = xx(1:end-1);
WW = exp(-i*2*pi*(kvecr*xx));

dat = WW*obj;


%W = 2;
%beta = 

% DONE IN grid1d ww = gradient(kvecr);
dato = grid1d(kvecr,dat,N);

% vs
kvecd = [-N/2:1:N/2-1].';

datr = interp1(kvecr,real(dat),kvecd,'cubic');
dati = interp1(kvecr,imag(dat),kvecd,'cubic');
datint = fftshift(fft(fftshift(datr+i*dati)));

