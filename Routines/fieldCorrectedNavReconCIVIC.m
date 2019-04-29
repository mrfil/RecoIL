function [img] = fieldCorrectedNavReconCIVIC(rInfo, sen, mask, FM, dcf, varargin)
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
L = ceil(abs(rInfo.timingVecNav(1) - rInfo.timingVecNav(end))/1E-3);
LValidationFcn = @(x) isnumeric(x) && isscalar(x) && (x >= 0);

Rbeta = 0; %Use a zero beta, not a good idea in practice.
RBetaValidationFcn = @(x) isnumeric(x) && isscalar(x) && (x >= 0);

delta = 0; % Use a zero delta, not a good idea for huber.
DeltaValidationFcn = @(x) isnumeric(x) && isscalar(x) && (x >= 0);

penalty = 'quad'; % Specify defalt penalty weighting
PenaltyValidationFcn = @(x) ischar(x);

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

%% Handle reconstruction
img = zeros(rInfo.NNav,rInfo.NNav,rInfo.nPartitionsNav,rInfo.nShots,rInfo.nPartitions,rInfo.nSlices,rInfo.nAverages,rInfo.nPhases,rInfo.nEchoes,rInfo.nRepetitions);
%Temporary workaround for error in MDH indexes in sequence
navData = rInfo.DataObject.RTfeedback.unsorted();
navData = permute(navData,[1,3,2]);
% Note that phases (phase cycles) come before everything else for
% the CIVIC steam sequence. For PGSE, phases = 1, so it doesn't matter
% either way.
navData = reshape(navData,[],rInfo.nPhases,rInfo.nShots,rInfo.nPartitions,rInfo.nSlices,rInfo.nAverages,rInfo.nEchoes,rInfo.nRepetitions,rInfo.nCoils);
navData = navData(1+rInfo.ptsToDrop:rInfo.ptsToDrop+rInfo.shotLengthNav,:,:,:,:,:,:,:,:);
navData = permute(navData,[1,3,4,5,6,2,7,8,9,10]);
% Need to phase the data properly...
halfVoxelShift = 0;
% Add half voxel shift to the data if 3D encoded
if rInfo.nPartitionsNav > 1
    halfVoxelShift = -1./(2*rInfo.nPartitionsNav);
end

% Phase data for off center shifts
for ii = 1:rInfo.nShots
    for jj = 1:rInfo.nPartitions
        for kk = 1:rInfo.nSlices
            for ll = 1:rInfo.nAverages
                for mm = 1:rInfo.nPhases
                    for nn = 1:rInfo.nEchoes
                        for oo = 1:rInfo.nRepetitions
                            for qq = 1:rInfo.nCoils
                                curData = squeeze(navData(:,ii,jj,kk,ll,mm,nn,oo,qq));
                                navData(:,ii,jj,kk,ll,mm,nn,oo,qq) = curData.*exp(-1j*(rInfo.readShift(rInfo.ExcOrder(kk)).*rInfo.kxNav(:,ii,jj,kk,ll,mm,nn,oo)+rInfo.phaseShift(rInfo.ExcOrder(kk)).*rInfo.kyNav(:,ii,jj,kk,ll,mm,nn,oo)-(rInfo.sliceShift(rInfo.ExcOrder(kk))-halfVoxelShift).*rInfo.kzNav(:,ii,jj,kk,ll,mm,nn,oo))*2*pi);
                            end
                        end
                    end
                end
            end
        end
    end
end



dcf = repmat(dcf,[rInfo.nCoils,1]);

for kk = 1:length(slicesToRecon)
    slc   = slicesToRecon(kk);
    
    R = Robj(logical(mask(:,:,:,slc)), ...
        'edge_type','tight', ...
        'order', 2, ...
        'beta', Rbeta, ...
        'type_denom', 'matlab', ...
        'potential', penalty, ...
        'delta', delta, ...
        'dims2penalize', dims2penalize);
    
    for ll = 1:length(averagesToRecon)
        avg   = averagesToRecon(ll);
        for mm = 1:length(phasesToRecon)
            phs   = phasesToRecon(mm);
            for nn = 1:length(echoesToRecon)
                eco   = echoesToRecon(nn);
                for oo = 1:length(repetitionsToRecon)
                    rep   = repetitionsToRecon(oo);
                    
                    % Setup Iterative Reconstruction
                    
                    if rInfo.nPartitionsNav > 1
                        G = NUFFT(col(rInfo.kxNav(rInfo.dataMaskNav,ii,jj,slc,avg,phs,eco,rep)), ...
                            col(rInfo.kyNav(rInfo.dataMaskNav,ii,jj,slc,avg,phs,eco,rep)), ...
                            col(rInfo.kzNav(rInfo.dataMaskNav,ii,jj,slc,avg,phs,eco,rep)), ...
                            rInfo.NNav, ...
                            rInfo.NNav, ...
                            rInfo.nPartitionsNav, ...
                            'mask',logical(mask(:,:,:,slc)));
                        A = TimeSegmentation(G, rInfo.timingVecNav(rInfo.dataMaskNav), ...
                            FM(logical(mask(:,:,:,slc))), ...
                            L);
                        %                             [Gx, Gy, Gz] = gradient(FM/(2*pi));
                        %                              A = Gdft_r2(col(rInfo.kxNav(rInfo.dataMaskNav,ii,jj,slc,avg,phs,eco,rep)), ...
                        %                                      col(rInfo.kyNav(rInfo.dataMaskNav,ii,jj,slc,avg,phs,eco,rep)), ...
                        %                                      col(rInfo.kzNav(rInfo.dataMaskNav,ii,jj,slc,avg,phs,eco,rep)), ...
                        %                                      rInfo.NNav, ...
                        %                                      rInfo.NNav, ...
                        %                                      rInfo.nPartitionsNav, ...
                        %                                      FM, ...
                        %                                      zeros(size(FM)), ...
                        %                                      rInfo.timingVecNav(rInfo.dataMaskNav), ...
                        %                                      Gx,Gy,Gz);
                    else % 2D Case
                        G = NUFFT(col(rInfo.kxNav(rInfo.dataMaskNav,ii,jj,slc,avg,phs,eco,rep)), ...
                            col(rInfo.kyNav(rInfo.dataMaskNav,ii,jj,slc,avg,phs,eco,rep)), ...
                            col(rInfo.kzNav(rInfo.dataMaskNav,ii,jj,slc,avg,phs,eco,rep)), ...
                            rInfo.NNav, ...
                            rInfo.NNav, ...
                            1, ...
                            'mask',logical(mask(:,:,:,slc)));
                        A = TimeSegmentation(G, rInfo.timingVecNav(rInfo.dataMaskNav), ...
                            FM(logical(mask(:,:,:,slc))), ...
                            L);
                    end
                    
                    
                    sen_tmp = reshape(sen(:,:,:,slc,:), ...
                        rInfo.NNav*rInfo.NNav*rInfo.nPartitionsNav, ...
                        []);
                    %S = sense_svd(A,sen_tmp(logical(col(mask(:,:,:,kk))),:));
                    beta11 = 6e4; % regularization paramter for phase (rg2/rg4)
                    soft = @(t,a) (t - a * sign(t)) .* (abs(t) > a); % soft thresholding function
                    hard = @(t,a) (t.*(abs(t) > a));
                    curv = 40*40*20;
                    beta22 = 0.5*2^-4*curv;
                    U = Wavelet3(rInfo.NNav,rInfo.NNav,rInfo.nPartitionsNav,logical(col(mask(:,:,:,slc))));
                    % Need to apply the data mask
                    S = sense(A,sen_tmp(logical(col(mask(:,:,:,slc))),:));
                   
                    %data = data(rInfo.dataMaskNav,:,:,:,:,:,:,:,:,:);
                    %dataPrepped = col(prepData(S,col(data)));
                    imginit = zeros(rInfo.NNav,rInfo.NNav,rInfo.nPartitionsNav);
                    img_i = imginit(logical(col(mask(:,:,:,slc))));
                    C = C_sparse(logical(mask(:,:,:,slc)),dims2penalize);
                    for ii = 1:rInfo.nShots
                        for jj = 1:rInfo.nPartitions
                            tic
                            data = double(navData(rInfo.dataMaskNav,ii,jj,slc,avg,phs,eco,rep,:));
                            
                            for tt = 1:Niter
                                % Iterate once to obtain image estimate
                                tmp = img_i + 1/curv * (S' * (dcf .* (data(:) - S * img_i)));
                                %img_i = U' * soft(U * tmp, beta22 / curv);
                                img_i = U' * soft(U * tmp, 5E-12);
                                %img_i = U' * hard(U * tmp, 5E-9);
                            end
                            
                            mi = abs(img_i);
                            phi = angle(img_i);
                            phi = col(prelude(reshape(phi,40,40,20),ones(40,40,20)));
%                             
%                             del = 1E10; % parameter for edge-preserving
%                             beta11 = 6e10; % regularization paramter for phase (rg2/rg4)
% 
% 
%                             for tt = 1:Niter
%                                 % Update Phase
%                                 %phi = pcg_bls_exp_ep(S, C, data(:), mi, phi, beta11, del, ...
%                                 %2, 'mask', logical(col(mask(:,:,:,slc))));
%                                 phi = pcg_bls_exp(S,C,data(:),mi,phi,beta11,2,'mask',logical(col(mask(:,:,:,slc))));
%                                 % Update Magnitude 
%                                 tmp = U * mi + 1/curv * real(U * (S' * (dcf .* (data(:) - S * (mi.*exp(1j*phi))))));
%                                 %img_i = U' * soft(U * tmp, beta22 / curv);
%                                 mi = abs(U' * soft(tmp, 5E-10));
%                                 
%                                 %img_i = mi.*exp(1j*phi);
%                                 %mi = abs(img_i);
%                                 %phi = col(prelude(reshape(angle(img_i),40,40,20),reshape(mi,40,40,20)));
% 
%                             end
                            %[nav, ~, norms] = solve_pwls_pcg(col(imginit(logical(mask(:,:,:,slc)))), S, 1, data(:), R, 'niter',Niter);
                            %[test, ~, norms] = solve_pwls_pcg(col(imginit(logical(mask(:,:,:,kk)))), S, 1, data(:), R, 'niter',Niter),logical(mask(:,:,:,kk));
                            %nav = mi.*exp(1j*phi);
                            img(:,:,:,ii,jj,kk,ll,mm,nn,oo) = embed(phi,logical(mask(:,:,:,kk)));
                            toc
                        end
                    end
                end
            end
        end
    end
end

end


