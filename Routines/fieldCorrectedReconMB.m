function [img, resid, roughness] = fieldCorrectedReconMB(rInfo, sen, mask, FM, varargin)
   %fieldCorrectedRecon - Basic Field Corrected Recon with SENSE
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
   % 11-Apr-2017; Last revision: 9-Sep-2017
   
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
   
   dims2penalize = [1 1 0];
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
   
   %Deal with the mask
   mask = logical(mask);
   if rInfo.multibandFactor > 1
    img = zeros(rInfo.N,rInfo.N,rInfo.multibandFactor,rInfo.nSlices, ...
                rInfo.nAverages,rInfo.nPhases,rInfo.nEchoes,rInfo.nRepetitions);
   else
    img = zeros(rInfo.N,rInfo.N,rInfo.nPartitions,rInfo.nSlices, ...
                rInfo.nAverages,rInfo.nPhases,rInfo.nEchoes,rInfo.nRepetitions);
   end

      %% Execute Recon

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
                  
                  data = rInfo.dataRead([],[],slc,avg,phs,eco,rep);
                  
                  tic
                  if nargout > 2 % Preallocate roughness variable to return
                     roughness = zeros(length(slicesToRecon), ...
                              length(averagesToRecon), ...
                              length(phasesToRecon), ...
                              length(echoesToRecon),...
                              length(repetitionsToRecon));
                  end
                  
                  % Setup Iterative Reconstruction
                  
                  if rInfo.multibandFactor > 1
                     G = NUFFT(col(rInfo.kx(rInfo.dataMask,:,:,slc,avg,phs,eco,rep)), ...
                               col(rInfo.ky(rInfo.dataMask,:,:,slc,avg,phs,eco,rep)), ...
                               col(rInfo.kz(rInfo.dataMask,:,:,slc,avg,phs,eco,rep)), ...
                               rInfo.N, ...
                               rInfo.N, ...
                               rInfo.multibandFactor, ...
                               'mask',mask(:,:,:,slc), ...
                               'VoxelBasis','boxcar');
                  elseif rInfo.nPartitions > 1
                     G = NUFFT(col(rInfo.kx(rInfo.dataMask,:,:,slc,avg,phs,eco,rep)), ...
                               col(rInfo.ky(rInfo.dataMask,:,:,slc,avg,phs,eco,rep)), ...
                               col(rInfo.kz(rInfo.dataMask,:,:,slc,avg,phs,eco,rep)), ...
                               rInfo.N, ...
                               rInfo.N, ...
                               rInfo.nPartitions, ...
                               'mask',mask(:,:,:,slc), ...
                               'VoxelBasis','boxcar');
                  else % 2D Case
                     G = NUFFT(col(rInfo.kx(rInfo.dataMask,:,:,slc,avg,phs,eco,rep)), ...
                               col(rInfo.ky(rInfo.dataMask,:,:,slc,avg,phs,eco,rep)), ...
                               col(rInfo.kz(rInfo.dataMask,:,:,slc,avg,phs,eco,rep)), ...
                               rInfo.N,...
                               rInfo.N,...
                               1,...
                               'mask',logical(mask(:,:,:,slc)));

                  end
                                       A = TimeSegmentation(G,...
                                          col(rInfo.timingVec(rInfo.dataMask)), ...
                                          col(FM(mask(:,:,:,slc))),...
                                          L);
                  R = Robj(logical(squeeze(mask(:,:,:,slc))),...
                           'edge_type','tight',...
                           'order',2,...
                           'beta',Rbeta,...
                           'type_denom','matlab',...
                           'potential',penalty,...
                           'delta',delta,...
                           'dims2penalize',dims2penalize);
                  
                  sen_tmp = reshape(sen(:,:,:,slc,:),rInfo.N*rInfo.N*rInfo.multibandFactor,[]);
                  
                  S = sense(A,sen_tmp(logical(col(mask(:,:,:,slc))),:));
                  data = col(data);
                  if rInfo.multibandFactor > 1 % 3D SMS case
                      imginit = zeros(rInfo.N,rInfo.N,rInfo.multibandFactor);
                  elseif rInfo.nPartitions > 1 % General 3D case
                      imginit = zeros(rInfo.N,rInfo.N,rInfo.nPartitions);
                  else  % General 2D case
                      imginit = zeros(rInfo.N,rInfo.N,1);
                  end
                  
                  xinit = col(imginit(squeeze(mask(:,:,:,slc))));

                  [ colImage, ~, resid{kk,ll,mm,nn,oo}] = solve_pwls_pcg(xinit,...
                                                                         S, ...
                                                                         1, ...
                                                                         data, ...
                                                                         R, ...
                                                                         'niter',Niter);
                  
                  if nargout > 2
                     roughness(kk,ll,mm,nn,oo) = R.penal(R,colImage);
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


