function [img] = fieldCorrectedNavRecon(rInfo, sen, mask, FM, varargin)
%fieldCorrectedNavRecon - Basic Field Corrected Navigator Recon with SENSE
%
% Syntax:  [img] = fieldCorrectedNavRecon(rInfo,sen,FM)
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
% Optionl Name - Value Pair Inputs:
%    Rbeta   - Add beta for regularization function. Default Beta is zero
%              which is not recommended.
%
%    Penalty - String recognized by Robject (or Robj) to specify penalty
%              weighting. Default is 'quad'. Options, such as 'huber'
%              exist.
%
%    delta   - Some penalties require a delta option, such as 'huber'. This
%              has no effect if the penalty specified does not use the
%              delta value.
%
%    L       - Overrides the number of time segments used. Otherwise sets
%              the number of time segments to be one every 3 miliseconds of
%              readout length.
%
%   Niter    - Sets the number of CG iterations used to solve for the image
%              by solve_pwls_pcg.m. Default is 10.
%
%
% Outputs:
%    img     - Reconstructed Coil Images of "standard dimensions"
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
% 11-Apr-2017; Last revision: 10-Sep-2017

%% Deal with Optional Inputs

% Establish Defaults
L = ceil(abs(rInfo.timingVecNav(1) - rInfo.timingVecNav(end))/3E-3);
LValidationFcn = @(x) isnumeric(x) && isscalar(x) && (x >= 0);

Rbeta = 0; %Use a zero beta, not a good idea in practice.
RBetaValidationFcn = @(x) isnumeric(x) && isscalar(x) && (x >= 0);

delta = 0; % Use a zero delta, not a good idea for huber.
DeltaValidationFcn = @(x) isnumeric(x) && isscalar(x) && (x >= 0);

penalty = 'quad'; % Specify defalt penalty weighting
PenaltyValidationFcn = @(x) isstring(x);

Niter = 10;
NIterValidationFcn = @(x) isnumeric(x) && isscalar(x) && (x > 0);

dims2penalize = [1 1 1];
dims2penalizeValidationFcn = @(x) isnumeric(x) && ~isscalar(x) && (length(x) <= 3);

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

%% Handle reconstruction prereqs


% Test to see if navigator kspace is the same for all avg, phs, eco, and
% reps. Normally this is true. If so, we don't need to reinitialize the
% NUFFT, TimeSegementation for each set of images.

testKSpaceRead = col(rInfo.kReadNav(rInfo.dataMaskNav,1,1,1,1,1,1,1));
testKSpaceRead = repmat(testKSpaceRead,[1,rInfo.nShots,rInfo.nPartitions,rInfo.nSlices,rInfo.nAverages,rInfo.nPhases, rInfo.nEchoes, rInfo.nRepetitions]);
KSpaceReadSame = isequal(testKSpaceRead,rInfo.kReadNav(rInfo.dataMaskNav,:,:,:,:,:,:,:));

testKSpacePhase = col(rInfo.kPhaseNav(rInfo.dataMaskNav,1,1,1,1,1,1,1));
testKSpacePhase = repmat(testKSpacePhase,[1,rInfo.nShots,rInfo.nPartitions,rInfo.nSlices,rInfo.nAverages,rInfo.nPhases, rInfo.nEchoes, rInfo.nRepetitions]);
KSpacePhaseSame = isequal(testKSpacePhase,rInfo.kPhaseNav(rInfo.dataMaskNav,:,:,:,:,:,:,:));

testKSpaceSlice = col(rInfo.kSliceNav(rInfo.dataMaskNav,1,1,1,1,1,1,1));
testKSpaceSlice = repmat(testKSpaceSlice,[1,rInfo.nShots,rInfo.nPartitions,rInfo.nSlices,rInfo.nAverages,rInfo.nPhases, rInfo.nEchoes, rInfo.nRepetitions]);
KSpaceSliceSame = isequal(testKSpaceSlice,rInfo.kSliceNav(rInfo.dataMaskNav,:,:,:,:,:,:,:));

if (KSpaceReadSame && KSpacePhaseSame && KSpaceSliceSame)
    NavigatorTrajectoryDifferent = false;
    if rInfo.nPartitionsNav > 1
        maskUsed = true(rInfo.NNav, rInfo.NNav, rInfo.nPartitionsNav,rInfo.nSlices);
        G = NUFFT(col(rInfo.kReadNav(rInfo.dataMaskNav,1,1,1,1,1,1,1)), ...
            col(rInfo.kPhaseNav(rInfo.dataMaskNav,1,1,1,1,1,1,1)), ...
            col(rInfo.kSliceNav(rInfo.dataMaskNav,1,1,1,1,1,1,1)), ...
            rInfo.NNav, ...
            rInfo.NNav, ...
            rInfo.nPartitionsNav, ...
            'mask',maskUsed(:,:,:,1));
    else % 2D Case
        maskUsed = true(rInfo.NNav, rInfo.NNav,1,rInfo.nSlices);
        G = NUFFT(col(rInfo.kReadNav(rInfo.dataMaskNav,1,1,1,1,1,1,1)), ...
            col(rInfo.kPhaseNav(rInfo.dataMaskNav,1,1,1,1,1,1,1)), ...
            col(rInfo.kSliceNav(rInfo.dataMaskNav,1,1,1,1,1,1,1)), ...
            rInfo.NNav, ...
            rInfo.NNav, ...
            1, ...
            'mask',maskUsed(:,:,:,1));
    end
else
    NavigatorTrajectoryDifferent = true;
    maskUsed = logical(mask);
end

img = zeros(rInfo.nSlices,rInfo.NNav,rInfo.NNav,rInfo.nPartitionsNav,rInfo.nShots,rInfo.nPartitions,rInfo.nAverages,rInfo.nPhases,rInfo.nEchoes,rInfo.nRepetitions);

parfor kk = 1:length(slicesToRecon)
%for kk = 1:length(slicesToRecon)
    if(NavigatorTrajectoryDifferent == true) 
            img(kk,:,:,:,:,:,:,:,:,:) = loopBody(kk, slicesToRecon,averagesToRecon,phasesToRecon,echoesToRecon,repetitionsToRecon, ...
                    maskUsed,rInfo,penalty, Rbeta,delta, dims2penalize, NavigatorTrajectoryDifferent, FM, sen, Niter, L);
    else
           img(kk,:,:,:,:,:,:,:,:,:) = loopBody(kk, slicesToRecon,averagesToRecon,phasesToRecon,echoesToRecon,repetitionsToRecon, ...
                    maskUsed,rInfo,penalty, Rbeta,delta, dims2penalize, NavigatorTrajectoryDifferent, FM, sen, Niter, L, G); 
    end
end

img = permute(img,[2,3,4,5,6,1,7,8,9,10]);

end

function img = loopBody(slcIdx, slicesToRecon,averagesToRecon,phasesToRecon,echoesToRecon,repetitionsToRecon, ...
                    maskUsed,rInfo,penalty, Rbeta,delta, dims2penalize, NavigatorTrajectoryDifferent, FM, sen, Niter, L, G)

imginit = zeros(rInfo.NNav,rInfo.NNav,rInfo.nPartitionsNav);
img = zeros(rInfo.NNav,rInfo.NNav,rInfo.nPartitionsNav,rInfo.nShots,rInfo.nPartitions,rInfo.nAverages,rInfo.nPhases,rInfo.nEchoes,rInfo.nRepetitions);

slc   = slicesToRecon(slcIdx);
sen_tmp = reshape(sen(:,:,:,slc,:), ...
    rInfo.NNav*rInfo.NNav*rInfo.nPartitionsNav, ...
    []);
R = Robj(maskUsed(:,:,:,slc), ...
    'beta', Rbeta, ...
    'potential', penalty, ...
    'delta', delta, ...
    'dims2penalize', dims2penalize);

for ll = 1:length(averagesToRecon)
    for mm = 1:length(phasesToRecon)
        for nn = 1:length(echoesToRecon)
            for oo = 1:length(repetitionsToRecon) % CLJ: adding parfor here to speedup nav recons
                
                %% dealing with DataMask
                
                % Now we need to look up the image indexes
                
                avg   = averagesToRecon(ll);
                phs   = phasesToRecon(mm);
                eco   = echoesToRecon(nn);
                rep   = repetitionsToRecon(oo);
                
                if (NavigatorTrajectoryDifferent == true)
                    if rInfo.nPartitionsNav > 1
                        G = NUFFT(col(rInfo.kReadNav(rInfo.dataMaskNav,1,1,slc,avg,phs,eco,rep)), ...
                            col(rInfo.kPhaseNav(rInfo.dataMaskNav,1,1,slc,avg,phs,eco,rep)), ...
                            col(rInfo.kSliceNav(rInfo.dataMaskNav,1,1,slc,avg,phs,eco,rep)), ...
                            rInfo.NNav, ...
                            rInfo.NNav, ...
                            rInfo.nPartitionsNav, ...
                            'mask',maskUsed(:,:,:,slc));
                        A = TimeSegmentation(G, rInfo.timingVecNav(rInfo.dataMaskNav), ...
                            reshape(FM(maskUsed(:,:,:,slc)),rInfo.NNav,rInfo.NNav,rInfo.nPartitionsNav), ...
                            L);
                    else % 2D Case
                        G = NUFFT(col(rInfo.kReadNav(rInfo.dataMaskNav,1,1,slc,avg,phs,eco,rep)), ...
                            col(rInfo.kPhaseNav(rInfo.dataMaskNav,1,1,slc,avg,phs,eco,rep)), ...
                            col(rInfo.kSliceNav(rInfo.dataMaskNav,1,1,slc,avg,phs,eco,rep)), ...
                            rInfo.NNav, ...
                            rInfo.NNav, ...
                            1, ...
                            'mask',maskUsed(:,:,:,slc));
                        A = TimeSegmentation(G, rInfo.timingVecNav(rInfo.dataMaskNav), ...
                            reshape(FM(maskUsed(:,:,:,slc)),rInfo.NNav,rInfo.NNav,1,1), ...
                            L);
                    end
                else
                    if rInfo.nPartitionsNav > 1
                        maskTrue = true(rInfo.NNav, rInfo.NNav, rInfo.nPartitionsNav, rInfo.nSlices);
                        A = TimeSegmentation(G, rInfo.timingVecNav(rInfo.dataMaskNav), ...
                            reshape(FM(maskUsed(:,:,:,slc)),rInfo.NNav,rInfo.NNav,rInfo.nPartitionsNav), ...
                            L);
                    else
                        maskTrue = true(rInfo.NNav, rInfo.NNav, 1, rInfo.nSlices);
                        A = TimeSegmentation(G, rInfo.timingVecNav(rInfo.dataMaskNav), ...
                            reshape(FM(maskUsed(:,:,:,slc)),rInfo.NNav,rInfo.NNav,1,1), ...
                            L);
                    end
                    
                end
                
                data = rInfo.navRead([],[],slc,avg,phs,eco,rep);
                if (NavigatorTrajectoryDifferent == true)
                    S = sense(A,sen_tmp(col(maskUsed(:,:,:,slc)),:));
                else
                    S = sense(A,sen_tmp(col(maskUsed(:,:,:,slc)),:));
                end
                
                imgTemp = zeros(rInfo.NNav,rInfo.NNav,rInfo.nPartitionsNav,rInfo.nShots,rInfo.nPartitions);
                
                for ii = 1:rInfo.nShots
                    for jj = 1:rInfo.nPartitions
                        % Setup Iterative Reconstruction
                        dataShot = data(rInfo.dataMaskNav,ii,jj,:,:,:,:,:,:,:);
                        temp = embed(solve_pwls_pcg(col(imginit(maskUsed(:,:,:,slc))), S, 1, dataShot(:), R, 'niter',Niter),maskUsed(:,:,:,slc));
                        imgTemp(:,:,:,ii,jj) = temp;
                    end
                end
                
                img(:,:,:,:,:,ll,mm,nn,oo) = imgTemp;
            end
        end
    end
end
end
