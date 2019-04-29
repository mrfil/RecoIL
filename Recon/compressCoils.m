function [rInfo,compressedSen] = compressCoils(rInfo,sen,varargin)
%compressCoils Compress coils for use with recoInfo and sense object
%   
% Inputs:
%   rInfo       - 

%% Process inputs

if rInfo.isCoilCompressed
    error('recoInfo:alreadyCompressed','compressCoils: Cannot compress coils for already compressed recoInfo object.');
end

sizeSen = size(sen);

coilRank = -1;
coilRankValidationFcn = @(x) (x > 0) && isscalar(x) && (x < sizeSen(end));

energyLevel = 0.95;
energyLevelValidationFcn = @(x) isnumeric(x) && (x > 0) && (x <= 1) && isscalar(x);

p = inputParser();
p.addOptional('coilRank', coilRank, coilRankValidationFcn);
p.addOptional('energyLevel', energyLevel, energyLevelValidationFcn);

parse(p,varargin{:});

coilRank = p.Results.coilRank;
energyLevel = p.Results.energyLevel;


%% Calculate SVD of sense map

% Reshape sense map so that is is npts x ncoils, assuming coils is last
% dimension
sen = reshape(sen,[],sizeSen(end));

[~,sig,V] = svd(sen,0);

if (coilRank == -1)
    sig = diag(sig);
    coilRank = chooseCoilRank(sig,energyLevel);
    disp(['compressCoil: automatic coil rank selection returned ' num2str(coilRank) ' virtual coils based on ' num2str(energyLevel) ' energy level.']);
else
    disp(['compressCoil: automatic coil rank selection overridden to be ' num2str(coilRank) ' virtual coils.']);
end
rInfo.coilCompMat = V(:,1:coilRank);


compressedSen = sen*rInfo.coilCompMat;
compressedSen = reshape(compressedSen,[sizeSen(1:end-1) coilRank]);

% Write back parameters to the file
rInfo.nCoils = coilRank;
rInfo.isCoilCompressed = true;
end

