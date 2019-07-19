function [img, resid, roughness] = fieldCorrectedRecon(rInfo, sen, mask, FM, varargin)
%fieldCorrectedRecon - Basic Field Corrected Recon with SENSE
% Supports 2D, 3D, 3D SMS and theoretically 3D multslab reconstructions
% For phase corrected recons, use phaseCorrectedRecon instead.
%
% Syntax:  [img] = fieldCorrectedRecon(rInfo,sen,FM)
%
% Inputs:
%    rInfo   - intialized recoInfo object
%    sen     - SENSE Maps corresponding to images to be reconstructed
%              If the SENSE maps are not of the same size, they need to be
%              interpolated before calling this function.
%    FM      - Field Maps corresponding to images to be reconstructed
%              If the Field Maps are not of the same size, they need to be
%              interpolated before calling this function.
%
% Optional Name - Value Pair Inputs:
%    Rbeta   - Add beta for regularization function. Default Beta is zero
%              which is not recommended.
%
%    Penalty - String recognized by Robject (or Robj) to specify penalty
%              weighting. Default is 'quad'. Options, such as 'huber'
%              exist.
%
%    delta   - Some penalties require a delta option, such as 'huber'. This
%              has no effect if the penalty specified does not use the
%              delta value. Default is zero.
%
%    L       - Overrides the number of time segments used. Otherwise, grabs
%              the value from the recoInfo object rInfo.
%
%   Niter    - Sets the number of CG iterations used to solve for the image
%              by solve_pwls_pcg.m. Default is 10.
%
%   slicesToRecon - Array of indicies that correspond to images to be
%                   reconstructed. Indicies can be in any order or
%                   discontinuous. Default is 1:rInfo.nSlices.
%
%   averagesToRecon - Array of indicies that correspond to images to be
%                     reconstructed. Indicies can be in any order or
%                     discontinuous. Default is 1:rInfo.nAverages.
%
%   phasesToRecon - Array of indicies that correspond to images to be
%                   reconstructed. Indicies can be in any order or
%                   discontinuous. Default is 1:rInfo.nPhases.
%
%   echoesToRecon - Array of indicies that correspond to images to be
%                   reconstructed. Indicies can be in any order or
%                   discontinuous. Default is 1:rInfo.nEchoes.
%
%   repetitionsToRecon - Array of indicies corresponding to images to be
%                        reconstructed. Indicies can be in any order or
%                        discontinuous. Default is 1:rInfo.nRepetitions.
%
%   spatialInterp - Interpolation method recognized by NUFFT.m
%
%
% Outputs:
%    img     - Reconstructed Coil Images of "standard dimensions"
%
% Example:
%    [cImages, sen, FM, FMImages] = ReconSenFM();
%    rInfo  = rInfo(filename);
%    images = gridCoilImages(rInfo);
%    im(sum(abs(images).^2,10)); % Sum of Squares image for
%                                            % first echo time in field map
%
%
% Other m-files required: Too numerous to list
% Subfunctions: none
% MAT-files required: none
%
% Author: Alex Cerjanic
% University of Illinois at Urbana-Champaign
% email address:
% Website:
% 11-Apr-2017; Last revision: 4-Jan-2018

%% Deal with Optional Inputs

% Establish Defaults
L = rInfo.L;
LValidationFcn = @(x) isnumeric(x) && isscalar(x) && (x >= 0);

Rbeta = 0; %Use a zero beta, not a good idea in practice.
RBetaValidationFcn = @(x) isnumeric(x) && isscalar(x) && (x >= 0);

delta = 0; % Use a zero delta, not a good idea for huber.
DeltaValidationFcn = @(x) isnumeric(x) && isscalar(x) && (x >= 0);

penalty = 'quad'; % Specify defalt penalty weighting
PenaltyValidationFcn = @(x) ischar(x);

tempInterp = 'minmax'; % Specify default temporal interpolator
TempInterpValidationFcn = @(x) ischar(x);

Niter = 10;
NIterValidationFcn = @(x) isnumeric(x) && isscalar(x) && (x > 0);

dims2penalize = [1 1 1];
dims2penalizeValidationFcn = @(x) isnumeric(x);

slicesToRecon = 1:rInfo.nSlices;
slicesToReconValidationFcn = @(x) isnumeric(x);

averagesToRecon = 1:rInfo.nAverages;
averagesToReconValidationFcn = @(x) isnumeric(x);

phasesToRecon = 1:rInfo.nPhases;
phasesToReconValidationFcn = @(x) isnumeric(x);

echoesToRecon = 1:rInfo.nEchoes;
echoesToReconValidationFcn = @(x) isnumeric(x);

repetitionsToRecon = 1:rInfo.nRepetitions;
repetitionsToReconValidationFcn = @(x) isnumeric(x);

timingVec = rInfo.timingVec;
timingVecValidationFcn = @(x) isnumeric(x);

spatialInterp = 'sparse'; % Specify default temporal interpolator
spatialInterpValidationFcn = @(x) ischar(x);

p = inputParser();
addOptional(p,'Rbeta',Rbeta,RBetaValidationFcn);
addOptional(p,'delta',delta,DeltaValidationFcn);
addOptional(p,'penalty',penalty,PenaltyValidationFcn);
addOptional(p,'L',L,LValidationFcn);
addOptional(p,'Niter', Niter, NIterValidationFcn);
addOptional(p,'dims2penalize', dims2penalize, dims2penalizeValidationFcn);
addOptional(p,'slicesToRecon', slicesToRecon, slicesToReconValidationFcn);
addOptional(p,'averagesToRecon', averagesToRecon, averagesToReconValidationFcn);
addOptional(p,'phasesToRecon', phasesToRecon, phasesToReconValidationFcn);
addOptional(p,'echoesToRecon', echoesToRecon, echoesToReconValidationFcn);
addOptional(p,'repetitionsToRecon', repetitionsToRecon, repetitionsToReconValidationFcn);
addOptional(p,'temporalInterp', tempInterp, TempInterpValidationFcn);
addOptional(p,'timingVec', timingVec, timingVecValidationFcn);
addOptional(p,'spatialInterp', spatialInterp, spatialInterpValidationFcn);

parse(p,varargin{:});

% Override default parameters with their successfully parsed values.
Rbeta = p.Results.Rbeta;
delta = p.Results.delta;
penalty = p.Results.penalty;
L = p.Results.L;
Niter = p.Results.Niter;
dims2penalize = p.Results.dims2penalize;
slicesToRecon = p.Results.slicesToRecon;
averagesToRecon = p.Results.averagesToRecon;
phasesToRecon = p.Results.phasesToRecon;
echoesToRecon = p.Results.echoesToRecon;
repetitionsToRecon = p.Results.repetitionsToRecon;
tempInterp = p.Results.temporalInterp;
timingVec = p.Results.timingVec;
spatialInterp = p.Results.spatialInterp;

%Deal with the mask
mask = logical(mask);
if rInfo.multibandFactor > 1
    img = zeros(rInfo.N, rInfo.N, rInfo.multibandFactor, rInfo.nSlices, ...
                rInfo.nAverages, rInfo.nPhases, rInfo.nEchoes, rInfo.nRepetitions);
else
    img = zeros(rInfo.N, rInfo.N, rInfo.nPartitions, rInfo.nSlices, ...
                rInfo.nAverages, rInfo.nPhases, rInfo.nEchoes, rInfo.nRepetitions);
end

%% Execute Recon
if nargout > 2 % Preallocate roughness variable to return
    roughness = zeros(length(slicesToRecon), ...
                        length(averagesToRecon), ...
                        length(phasesToRecon), ...
                        length(echoesToRecon),...
                        length(repetitionsToRecon));
end
if rInfo.multibandFactor > 1 % 3D SMS case
    imginit = zeros(rInfo.N, rInfo.N, rInfo.multibandFactor);
elseif rInfo.nPartitions > 1 % General 3D case
    imginit = zeros(rInfo.N, rInfo.N, rInfo.nPartitions);
else  % General 2D case
    imginit = zeros(rInfo.N, rInfo.N, 1);
end

for kk = 1:length(slicesToRecon)
    % We only need to reestablish the FT objects and time segmentation
    % for new slices
    slc   = slicesToRecon(kk);
    
    for ll = 1:length(averagesToRecon)
        for mm = 1:length(phasesToRecon)
            for nn = 1:length(echoesToRecon)
                for oo = 1:length(repetitionsToRecon)
                    %% dealing with image reconstruction
                    
                    % Now we need to look up the image indexes
                    avg   = averagesToRecon(ll);
                    phs   = phasesToRecon(mm);
                    eco   = echoesToRecon(nn);
                    rep   = repetitionsToRecon(oo);
                    
                    data = rInfo.dataRead([],[],slc,avg,phs,eco,rep);
                    
                    tic
                    
                    % Setup Iterative Reconstruction
                    
                    if rInfo.multibandFactor > 1
                        G = NUFFT(col(rInfo.kRead(rInfo.dataMask,:,:,slc,avg,phs,eco,rep)), ...
                            col(rInfo.kPhase(rInfo.dataMask,:,:,slc,avg,phs,eco,rep)), ...
                            col(rInfo.kSlice(rInfo.dataMask,:,:,slc,avg,phs,eco,rep)), ...
                            rInfo.N, ...
                            rInfo.N, ...
                            rInfo.multibandFactor, ...
                            'mask',logical(mask(:,:,:,slc)), ...
                            'VoxelBasis','boxcar', ...
                            'InterpMethod',spatialInterp);
                        FMvol = FM(:,:,:,slc);
                    elseif rInfo.nPartitions > 1
                        G = NUFFT(col(rInfo.kRead(rInfo.dataMask,:,:,slc,avg,phs,eco,rep)), ...
                            col(rInfo.kPhase(rInfo.dataMask,:,:,slc,avg,phs,eco,rep)), ...
                            col(rInfo.kSlice(rInfo.dataMask,:,:,slc,avg,phs,eco,rep)), ...
                            rInfo.N, ...
                            rInfo.N, ...
                            rInfo.nPartitions, ...
                            'mask',logical(mask(:,:,:,slc)), ...
                            'VoxelBasis','boxcar', ...
                            'InterpMethod',spatialInterp);
                        FMvol = FM(:,:,:,slc);
                    else % 2D Case
                        G = NUFFT(col(rInfo.kRead(rInfo.dataMask,:,:,slc,avg,phs,eco,rep)), ...
                            col(rInfo.kPhase(rInfo.dataMask,:,:,slc,avg,phs,eco,rep)), ...
                            col(rInfo.kSlice(rInfo.dataMask,:,:,slc,avg,phs,eco,rep)), ...
                            rInfo.N,...
                            rInfo.N,...
                            1,...
                            'mask',logical(mask(:,:,:,slc)), ...
                            'InterpMethod',spatialInterp);
                        FMvol = FM(:,:,1,slc);
                    end
                    A = TimeSegmentation(G,...
                        col(timingVec(rInfo.dataMask,:,:,slc,avg,phs,eco,rep)), ...
                        col(FMvol(squeeze(mask(:,:,:,slc)))),...
                        L,'Interpolator',tempInterp);
                    R = Robj(logical(squeeze(mask(:,:,:,slc))),...
                        'beta',Rbeta,...
                        'potential',penalty,...
                        'delta',delta,...
                        'dims2penalize',dims2penalize);
                    if rInfo.multibandFactor > 1 % 3D SMS case
                        sen_tmp = reshape(sen(:,:,:,slc,:), rInfo.N*rInfo.N*rInfo.multibandFactor,[]);
                    elseif rInfo.nPartitions > 1 % General 3D case
                        sen_tmp = reshape(sen(:,:,:,slc,:), rInfo.N*rInfo.N*rInfo.nPartitions,[]);
                    else  % General 2D case
                        sen_tmp = reshape(sen(:,:,:,slc,:), rInfo.N*rInfo.N,[]);
                    end
                    
                    S = sense(A,sen_tmp(logical(col(mask(:,:,:,slc))),:));
                    data = col(data);
                    
                    xinit = col(imginit(squeeze(mask(:,:,:,slc))));
                    
                    [ colImage, ~, resid{kk,ll,mm,nn,oo}] = solve_pwls_pcg(xinit,...
                        S, ...
                        1, ...
                        data, ...
                        R, ...
                        'niter',Niter);
                    
                    if nargout > 2
                        roughness(kk,ll,mm,nn,oo) = R.penal(colImage);
                    end
                    
                    img(:,:,:,kk,ll,mm,nn,oo) = embed(colImage,mask(:,:,:,slc));
                    
                    %[test, ~, norms] = solve_pwls_pcg(col(imginit(logical(squeeze(mask(:,:,:,kk))))), S, 1, data, R, 'niter',Niter);
                    %img(:,:,:,kk,ll,mm,nn,oo) = embed(test,logical(mask(:,:,:,kk)));
                    toc
                end
            end
        end
    end
end

end


