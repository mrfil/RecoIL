function [dcf] = calc3DDensityCompensation(kx,ky,kz,N,Nz)
%calc3DDensityCompensation Calculate 3D compensation function for an
%arbitrary trajectory. Uses approach sdc3_MAT.c and sdc3grid_kernel.c from Jim
%Pipe's group and Nick Zwart implemented with IRT functions
%   Inputs: 
%       kx - column vector of kx trajectory [-N/2:N/2-1]
%       ky - column vector of ky trajectory [-N/2:N/2-1]
%       kz - column vector of kz trajectory [-Nz/2:Nz/2-1]
%
%   Output:
%       dcf - density compensation function


% Settings from testmex.m
numIter = 25;
osf     = 3;

% Prepping inputs for input to sdc3_MAT. Important to avoid segfaulting

maxKx = max(kx);
minKx = min(kx);
maxKy = max(ky);
minKy = min(ky);
maxKz = max(kz);
minKz = min(kz);

maxes = [maxKx,maxKy,maxKz,abs(minKx),abs(minKy),abs(minKz)];

G = NUFFT(col(kx),col(ky),col(kz),N,N,Nz,'GridOverSampling',osf);

ww = ones(size(col(kx)));

for ii = 1:numIter
    C = G.st.p * (G.st.p' * ww);
    ww = ww ./ real(C);
end

dcf = ww;

end

