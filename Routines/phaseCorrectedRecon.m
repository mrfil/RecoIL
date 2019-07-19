function [img, resid, roughness] = phaseCorrectedRecon(rInfo, sen, mask, FM, PMaps, varargin)
   %phaseCorrectedRecon - Basic Phase Corrected Recon with pcSENSE
   %
   % Syntax:  [img, resid, roughness] = phaseCorrectedRecon(rInfo, sen, ...
   %                                                          mask, FM, PMaps)
   %
   % Inputs:
   %    rInfo   - intialized recoInfo object
   %    sen     - SENSE Maps corresponding to images to be reconstructed
   %              If the SENSE maps are not of the same size, they need to be
   %              interpolated before calling this function.
   %    mask    - Support mask correspinding to images to be reconstructed
   %    FM      - Field Maps corresponding to images to be reconstructed
   %              If the Field Maps are not of the same size, they need to be
   %              interpolated before calling this function.
   %    PMaps   - Phase maps measuring motion induced phase error in same
   %              space as images to be reconstructed. Usually acquired via a
   %              navigator image.
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
   %   Dims2Penalize - Sets the dimensions where the penalty is applied. For
   %                   SMS techniques, the 3rd (Z) dimension often should not
   %                   be penalized across. Default is all 3 dimensions.
   %
   %
   %
   % Outputs:
   %    img     - Reconstructed Images of "standard dimensions"
   %    residuals - (Optional) Calculate residuals for each iteration
   %    roughness - (Optional) Calculated rougness penalty for final solution in
   %                    img.
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
   % 11-Apr-2017; Last revision: 10-Sep-2017
   
   %% Deal with Optional Inputs
   
   % Establish Defaults
   L = rInfo.L;
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
   averagesToRecon = p.Results.averagesToRecon; %AMC: Fixed error here
   phasesToRecon = p.Results.phasesToRecon;
   echoesToRecon = p.Results.echoesToRecon;
   repetitionsToRecon = p.Results.repetitionsToRecon;
   
   
   %% Deal SMS vs. 3D
   if rInfo.multibandFactor > 1 % We are in SMS
      Nz = rInfo.multibandFactor;
   else
      Nz = rInfo.nPartitions;
   end
   
   
   %% Deal with Recon.
   
   img = zeros(rInfo.N,rInfo.N,Nz,length(slicesToRecon), ...
            length(averagesToRecon),length(phasesToRecon), ...
            length(echoesToRecon),length(repetitionsToRecon));
   
   for kk = 1:length(slicesToRecon)
      for ll = 1:length(averagesToRecon)
         for mm = 1:length(phasesToRecon)
            for nn = 1:length(echoesToRecon)
               for oo = 1:length(repetitionsToRecon)
                  %% dealing with image reconstruction
                  
                  % Now we need to look up the image indexes
                  slc   = slicesToRecon(kk);
                  avg   = averagesToRecon(ll);
                  phs   = phasesToRecon(mm);
                  eco   = echoesToRecon(nn);
                  rep   = repetitionsToRecon(oo);
                  
                  FMtmp = FM(:,:,:,slc); % CLJ: dealing with the 4D mask indexing problem
                  
                  %% dealing with DataMask
                  data = rInfo.dataRead([],[],slc,avg,phs,eco,rep);
                  data = permute(data,[1,10,2:9]);
                  tic
                  if nargout > 2 % Preallocate roughness variable to return
                     roughness = zeros(length(slicesToRecon), ...
                              length(averagesToRecon), ...
                              length(phasesToRecon), ...
                              length(echoesToRecon),...
                              length(repetitionsToRecon));
                  end
                  % Setup Iterative Reconstruction
                  
                  % if rInfo.nPartitions > 1
                  if (rInfo.nPartitions > 1) || (rInfo.multibandFactor > 1) % CLJ: I think this is the correct config for MB imaging
                     G = NUFFT(reshape(rInfo.kRead(rInfo.dataMask,:,:,slc,avg,phs,eco,rep),rInfo.shotLength,[]), ...
                        reshape(rInfo.kPhase(rInfo.dataMask,:,:,slc,avg,phs,eco,rep),rInfo.shotLength,[]), ...
                        reshape(rInfo.kSlice(rInfo.dataMask,:,:,slc,avg,phs,eco,rep),rInfo.shotLength,[]), ...
                        rInfo.N, ...
                        rInfo.N, ...
                        Nz, ...
                        'mask',logical(mask(:,:,:,slc)));
                     A = TimeSegmentation(G, ...
                        col(rInfo.timingVec(rInfo.dataMask)), ...
                        col(FMtmp(logical(mask(:,:,:,slc)))), ... % CLJ: updated to FMtmp
                        L);
                  else % 2D Case
                      % CLJ note: does this work?
                     G = NUFFT(reshape(rInfo.kRead(rInfo.dataMask,:,:,slc,avg,phs,eco,rep),rInfo.shotLength,[]), ...
                        reshape(rInfo.kPhase(rInfo.dataMask,:,:,slc,avg,phs,eco,rep),rInfo.shotLength,[]), ...
                        reshape(rInfo.kSlice(rInfo.dataMask,:,:,slc,avg,phs,eco,rep),rInfo.shotLength,[]), ...
                        rInfo.N, ...
                        rInfo.N, ...
                        1, ...
                        'mask',logical(mask(:,:,:,slc)));
                     A = TimeSegmentation(G, ...
                        col(rInfo.timingVec(rInfo.dataMask)), ...
                        col(FMtmp(logical(mask(:,:,:,slc)))), ... % CLJ: updated to FMtmp
                        L);
                  end
                  
                  if rInfo.multibandFactor > 1
                     R = Robj(logical(mask(:,:,:,slc)),'beta',Rbeta,'potential',penalty,'delta',delta,'dims2penalize',dims2penalize);
                  else
                     R = Robject(logical(mask(:,:,:,slc)),'beta',Rbeta,'potential',penalty,'delta',delta);
                  end
                  
                  sen_tmp = reshape(sen(:,:,:,slc,:),rInfo.N*rInfo.N*Nz,[]);
                  %kRead = reshape(rInfo.kRead(rInfo.dataMask,:,:,slc,avg,phs,eco,rep),sum(rInfo.dataMask(:)),rInfo.nPartitions*rInfo.nShots);
                  %kPhase = reshape(rInfo.kPhase(rInfo.dataMask,:,:,slc,avg,phs,eco,rep),rInfo.nPartitions*rInfo.nShots);
                  %kSlice = reshape(rInfo.kSlice(rInfo.dataMask,:,:,slc,avg,phs,eco,rep),rInfo.nPartitions*rInfo.nShots);
                  %PMaps = zeros(size(PMaps)); % Hack for testing - AMC!!!!
                  P = reshape(PMaps(:,:,:,:,:,slc,avg,phs,eco,rep),rInfo.N,rInfo.N,Nz,rInfo.nShots*rInfo.nPartitions);
                  
                  % S = pcSENSE(A,P,sen_tmp(logical(col(mask(:,:,:,kk))),:));
                  S = pcSENSE(A,P,sen_tmp(logical(col(mask(:,:,:,slc))),:)); % CLJ: changing this to slc, instead of the index kk, for situations where slc ~= kk
                  data = col(data);
                  
                  % AMC: Dealing multiband and 3D case here
                  if (rInfo.multibandFactor > 1)
                    imginit = zeros(rInfo.N,rInfo.N,rInfo.multibandFactor);
                  else
                    imginit = zeros(rInfo.N,rInfo.N,rInfo.nPartitions);
                  end
                  
                  xinit = col(imginit(squeeze(logical(mask(:,:,:,slc)))));
                  
                  [colImage, ~, resid{kk,ll,mm,nn,oo}] = solve_pwls_pcg(xinit,...
                                                                     S, ...
                                                                     1, ...
                                                                     data, ...
                                                                     R, ...
                                                                     'niter',Niter);
                  if nargout > 2
                     roughness(kk,ll,mm,nn,oo) = R*colImage;
                  end
                  
                  img(:,:,:,kk,ll,mm,nn,oo) = embed(colImage,logical(mask(:,:,:,slc)));
                  toc
               end
            end
         end
      end
   end
   
end


