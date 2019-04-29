function noiseCorrMatrix = calcNoiseCorr(NoiseScan)
% calcNoiseCorr Calculate the Nosie Correlation matrix from noise scan data
% 
%   Inputs: 
%   NoiseScan = Ncoils x Nadc x Nrep matrix of complex data 
%
%   Outputs: 
%   noiseCorrMatrix = Ncoils X Ncoils matrix of complex noise correlations

sizes = size(NoiseScan);

if length(sizes) < 2
    error('Noise Scan must be a 2D matrix and prefereably a 3D matrix!');
elseif length(sizes) < 3
    warning('Multiple Noise Scan ADCs are STRONGLY reccomended. One scan is insufficient to obtain a reliably noise correlation matrix!');
end

Nro_adc = sizes(2); %We assume the data is Ncoils x Nro_adc x Nreps
%PreAllocate data
NoiseCorr = zeros(sizes(1),sizes(1),sizes(3));

%Calculate the noise correlation matricices using by normalizing by the ADC
%length and then taking the outer product of the noise scan data. Resulting
%matrix will be Ncoils x Ncoils
for ii = 1:sizes(3)-1
    %NoiseCorr(:,:,ii) = 1./(Nro_adc-1)*(squeeze(NoiseScan(:,:,ii))*conj(squeeze(NoiseScan(:,:,ii)).'));
    NoiseCorr(:,:,ii) = cov(squeeze(NoiseScan(:,:,ii)).');
end

%Average the resulting Noise Correlation matrices
noiseCorrMatrix = mean(NoiseCorr,3);

end