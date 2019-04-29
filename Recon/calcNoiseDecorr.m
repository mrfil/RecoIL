function noiseDecorrMatrix = calcNoiseDecorr(NoiseCorr)
% calcNoiseDecorr Calculate the Nosie decorrelation from the noise
% correlation matrix
% 
%   Inputs: 
%   NoiseCorr = Ncoils x Ncoils matrix of complex data 
%
%   Outputs: 
%   noiseDecorrMatrix = Ncoils X Ncoils matrix

sizes = size(NoiseCorr);

if length(sizes) < 2
    error('NoiseCorr must be a square 2D matrix!');
elseif sizes(1) ~= sizes(2)
    error('NoiseCorr must be a square 2D matrix!');
end

%Clean up any tiny imaginary numbers on the main diagonal
NoiseCorr(logical(eye(size(NoiseCorr)))) = abs(diag(NoiseCorr));

%Calculate the noise decorrelation matrix using the Choleski decomposition
noiseDecorrMatrix = inv(chol(NoiseCorr,'lower'));

end