function [img] = fieldCorrectedNavReconL1(rInfo, sen, mask, FM, varargin)
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
averagesToRecon = p.Results.echoesToRecon;
phasesToRecon = p.Results.phasesToRecon;
echoesToRecon = p.Results.echoesToRecon;
repetitionsToRecon = p.Results.repetitionsToRecon;

%% Handle reconstruction
img = zeros(rInfo.NNav,rInfo.NNav,rInfo.nPartitionsNav,rInfo.nShots,rInfo.nPartitions,rInfo.nSlices,rInfo.nAverages,rInfo.nPhases,rInfo.nEchoes,rInfo.nRepetitions);
% Temporary workaround for error in MDH indexes in sequence
% navData = rInfo.DataObject.RTfeedback.unsorted();
% navData = permute(navData,[1,3,2]);
% navData = reshape(navData,[],rInfo.nShots,rInfo.nPartitions,rInfo.nSlices,rInfo.nAverages,rInfo.nPhases,rInfo.nEchoes,rInfo.nRepetitions,rInfo.nCoils);
% navData = navData(1+rInfo.ptsToDrop:rInfo.ptsToDrop+rInfo.ShotLengthNav,:,:,:,:,:,:,:,:);

%dcf = repmat(dcf,[rInfo.nCoils,1]);
for ii = 1:rInfo.nShots
    for jj = 1:rInfo.nPartitions
        for kk = 1:length(slicesToRecon)
            for ll = 1:length(averagesToRecon)
                for mm = 1:length(phasesToRecon)
                    for nn = 1:length(echoesToRecon)
                        for oo = 1:length(repetitionsToRecon) % CLJ: adding parfor here to speedup nav recons
                            
                            %% dealing with DataMask
                            % Temporary workaround
                            
                            %data = navData(:,ii,jj,kk,ll,mm,nn,oo,:);
                            
                            % Now we need to look up the image indexes
                            slc   = slicesToRecon(kk);
                            avg   = averagesToRecon(ll);
                            phs   = phasesToRecon(mm);
                            eco   = echoesToRecon(nn);
                            rep   = repetitionsToRecon(oo);
                            
                            data = rInfo.navRead(ii,jj,slc,avg,phs,eco,rep);
                            tic
                            
                            % Setup Iterative Reconstruction
                            
                            if rInfo.nPartitionsNav > 1
                                G = NUFFT(col(rInfo.kxNav(rInfo.dataMaskNav,ii,jj,slc,avg,phs,eco,rep)), ...
                                    col(rInfo.kyNav(rInfo.dataMaskNav,ii,jj,slc,avg,phs,eco,rep)), ...
                                    col(rInfo.kzNav(rInfo.dataMaskNav,ii,jj,slc,avg,phs,eco,rep)), ...
                                    rInfo.NNav, ...
                                    rInfo.NNav, ...
                                    rInfo.nPartitionsNav, ...
                                    'mask',logical(mask(:,:,:,slc)));
                                A = TimeSegmentation(G, rInfo.timingVecNav, ...
                                    FM(logical(mask(:,:,:,slc))), ...
                                    L);
                            else % 2D Case
                                G = NUFFT(col(rInfo.kxNav(rInfo.dataMaskNav,ii,jj,slc,avg,phs,eco,rep)), ...
                                    col(rInfo.kyNav(rInfo.dataMaskNav,ii,jj,slc,avg,phs,eco,rep)), ...
                                    col(rInfo.kzNav(rInfo.dataMaskNav,ii,jj,slc,avg,phs,eco,rep)), ...
                                    rInfo.NNav, ...
                                    rInfo.NNav, ...
                                    1, ...
                                    'mask',logical(mask(:,:,:,slc)));
                                A = TimeSegmentation(G, rInfo.timingVecNav, ...
                                    FM(logical(mask(:,:,:,slc))), ...
                                    L);
                            end
                            R = Robj(logical(mask(:,:,:,slc)), ...
                                'edge_type','tight', ...
                                'order', 2, ...
                                'beta', Rbeta, ...
                                'type_denom', 'matlab', ...
                                'potential', penalty, ...
                                'delta', delta, ...
                                'dims2penalize', dims2penalize);
                            
                            sen_tmp = reshape(sen(:,:,:,slc,:), ...
                                rInfo.NNav*rInfo.NNav*rInfo.nPartitionsNav, ...
                                []);
                            %S = sense_svd(A,sen_tmp(logical(col(mask(:,:,:,kk))),:));
                            beta11 = 6e4; % regularization paramter for phase (rg2/rg4)
                            soft = @(t,a) (t - a * sign(t)) .* (abs(t) > a); % soft thresholding function
                            hard = @(t,a) (t.*(abs(t) > a));
                            curv = rInfo.NNav*rInfo.NNav*rInfo.nPartitionsNav;
                            beta22 = 0.5*2^-4*curv;
                            U = Wavelet3(rInfo.NNav,rInfo.NNav,rInfo.nPartitionsNav,logical(col(mask(:,:,:,slc))));
                            
                            % Need to apply the data mask
                            S = sense(A,sen_tmp(logical(col(mask(:,:,:,slc))),:));
                            data = data(rInfo.dataMaskNav,:,:,:,:,:,:,:,:,:);
                            %dataPrepped = col(prepData(S,col(data)));
                            imginit = zeros(rInfo.NNav,rInfo.NNav,rInfo.nPartitionsNav);
                            imginit = zeros(rInfo.NNav,rInfo.NNav,rInfo.nPartitionsNav);
                            img_i = imginit(logical(col(mask(:,:,:,slc))));
                            C = C_sparse(logical(mask(:,:,:,slc)),dims2penalize);
%                             for tt = 1:600
%                                 % Iterate once to obtain image estimate
%                                 tmp = img_i + 1/curv * (S' * (dcf .* (data(:) - S * img_i)));
%                                 %img_i = U' * soft(U * tmp, beta22 / curv);
%                                 %img_i = U' * soft(U * tmp, 5E-8);
%                                 img_i = U' * soft(U * tmp, 1E-8);
%                             end
                            img_i = solve_pwls_pcg(col(imginit(logical(mask(:,:,:,slc)))), S, 1, data(:), R, 'niter',Niter);
                            mi = abs(img_i);
                            phi = angle(img_i);       
                            phi = col(prelude(embed(phi,logical(mask(:,:,:,slc))),embed(ones(size(phi)),logical(mask(:,:,:,slc)))));
                            phi = phi(logical(mask(:,:,:,slc)));
%                             del = 1E5; % parameter for edge-preserving
%                             beta11 = 6e2; % regularization paramter for phase (rg2/rg4)
% 
%                             for tt = 1:Niter
%                                 % Update Phase
%                                 %phi = pcg_bls_exp_ep(S, C, data(:), mi, phi, beta11, del, ...
%                                 %2, 'mask', logical(col(mask(:,:,:,slc))));
%                                 phi = pcg_bls_exp(S,C,data(:),mi,phi,beta11,2,'mask',logical(col(mask(:,:,:,slc))));
%                                 % Update Magnitude 
%                                 tmp = U * mi + 1/curv * (U * (real((S' * (dcf .* ((data(:) - S * (mi.*exp(1j*phi)))))))));
%                                 %img_i = U' * soft(U * tmp, beta22 / curv);
%                                 mi = abs(U' * soft(tmp, 1E-8));
%                                 %mi = abs(U' * soft(tmp, 1E-9));
%                                 %img_i = mi.*exp(1j*phi);
%                                 %mi = abs(img_i);
%                                 %phi = col(prelude(reshape(angle(img_i),40,40,20),reshape(mi,40,40,20)));
%                             end
                            img_i = phi;
                            img(:,:,:,ii,jj,kk,ll,mm,nn,oo) = embed(img_i,logical(mask(:,:,:,slc)));

                            %img(:,:,:,ii,jj,kk,ll,mm,nn,oo) = embed(solve_pwls_pcg(col(imginit(logical(mask(:,:,:,slc)))), S, 1, data(:), R, 'niter',Niter),logical(mask(:,:,:,slc)));
                            %[test, ~, norms] = solve_pwls_pcg(col(imginit(logical(mask(:,:,:,kk)))), S, 1, data(:), R, 'niter',Niter),logical(mask(:,:,:,kk));
                            %img(:,:,:,ii,jj,kk,ll,mm,nn,oo) = embed(test,logical(mask(:,:,:,kk)));
                            toc
                        end
                    end
                end
            end
        end
    end
end

end


