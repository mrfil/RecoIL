function [img] = fieldCorrectedReconOld(rInfo, sen, mask, FM, varargin)
%fieldCorrectedRecon - Basic Field Corrected Recon with SENSE
%
% Syntax:  [img] = fieldCorrectedRecon(rInfo,sen,mask,FM)
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

dims2penalize = [1 1 1];
dims2penalizeValidationFcn = @(x) isnumeric(x) && ~isscalar(x) && (length(x) <= 3);

p = inputParser();
addOptional(p,'Rbeta',Rbeta,RBetaValidationFcn);
addOptional(p,'delta',delta,DeltaValidationFcn);
addOptional(p,'penalty',penalty,PenaltyValidationFcn);
addOptional(p,'L',L,LValidationFcn);
addOptional(p,'Niter', Niter, NIterValidationFcn);
addOptional(p,'dims2penalize', dims2penalize, dims2penalizeValidationFcn);

parse(p,varargin{:});

% Override default parameters with their successfully parsed values.
Rbeta = p.Results.Rbeta;
delta = p.Results.delta;
penalty = p.Results.penalty;
L = p.Results.L;
Niter = p.Results.Niter;
dims2penalize = p.Results.dims2penalize;


%% Execute Recon
img = zeros(rInfo.N,rInfo.N,rInfo.nPartitions,rInfo.nSlices,rInfo.nAverages,rInfo.nPhases,rInfo.nEchoes,rInfo.nRepetitions);
for kk = 30
    for ll = 1:rInfo.nAverages
        for mm = 1:rInfo.nPhases
            for nn = 1:rInfo.nEchoes
                for oo = 1:rInfo.nRepetitions
                %for oo = 1:1
                    %% dealing with DataMask
                    data = rInfo.dataRead([],[],kk,ll,mm,nn,oo);
                    tic
                    
                    % Setup Iterative Reconstruction
                    
                    if rInfo.nPartitions > 1
                        G = NUFFT(col(rInfo.kx(rInfo.dataMask,:,:,kk,ll,mm,nn,oo)), col(rInfo.ky(rInfo.dataMask,:,:,kk,ll,mm,nn,oo)), col(rInfo.kz(rInfo.dataMask,:,:,kk,ll,mm,nn,oo)), rInfo.N,rInfo.N,rInfo.nPartitions,'mask',logical(mask(:,:,:,kk)));
                        A = TimeSegmentation(G,col(rInfo.timingVec(rInfo.dataMask)),col(FM(logical(mask(:,:,:,kk)))),L);
                    else % 2D Case
                        G = NUFFT(col(rInfo.kx(rInfo.dataMask,:,:,kk,ll,mm,nn,oo)), col(rInfo.ky(rInfo.dataMask,:,:,kk,ll,mm,nn,oo)), col(rInfo.kz(rInfo.dataMask,:,:,kk,ll,mm,nn,oo)), rInfo.N,rInfo.N,1,'mask',logical(mask(:,:,:,kk)));
                        %FMtemp = squeeze(FM(:,:,:,kk));
                        for ii = 1:rInfo.nSlices
                           FMtemp(:,:,1,ii) = squeeze(FM(:,:,1,ii));
                           %FMtemp(:,:,1,ii) = wiener2(squeeze(FM(:,:,1,ii)),[5,5]);
                           %FMtemp(:,:,1,ii) = squeeze(FM(:,:,1,ii));
                        end
                        [Gx,Gy,~,Gz] = gradient(FMtemp./(2*pi));
                        %A = TimeSegmentation(G,col(rInfo.timingVec(rInfo.dataMask)),col(FMtemp(:,:,1,kk)),L);
                        A = Gdft_r2(col(rInfo.ky(rInfo.dataMask,:,:,kk,ll,mm,nn,oo)), col(rInfo.kx(rInfo.dataMask,:,:,kk,ll,mm,nn,oo)), col(rInfo.kz(rInfo.dataMask,:,:,kk,ll,mm,nn,oo)),rInfo.N,rInfo.N,1,col(FMtemp(:,:,1,kk)),zeros(size(FMtemp(:,:,1,kk))),col(rInfo.timingVec(rInfo.dataMask)),col(Gy(:,:,1,kk)),col(Gx(:,:,1,kk)),col(Gz(:,:,1,kk)));
                    end
                    
                    R = Robj(logical(squeeze(mask(:,:,:,kk))),'edge_type','tight','order',2,'beta',Rbeta,'type_denom','matlab','potential',penalty,'delta',delta);
                    
                    sen_tmp = reshape(sen(:,:,:,kk,:),rInfo.N*rInfo.N*rInfo.nPartitions,[]);
                    
                    S = sense_svd(A,sen_tmp(logical(col(mask(:,:,:,kk))),:),'CoilRank',1);
                    %S = sense(A,sen_tmp(logical(col(mask(:,:,:,kk))),:));
                    data = col(data);
                    data = prepData(S,data);
                    imginit = zeros(rInfo.N,rInfo.N,rInfo.nPartitions);
                        
                    %img(:,:,:,kk,ll,mm,nn,oo) = embed(solve_pwls_pcg(col(imginit(logical(squeeze(mask(:,:,:,kk))))), S, 1, data, R, 'niter',Niter),logical(mask(:,:,:,kk)));
                    [test, ~, norms] = solve_pwls_pcg(col(imginit(logical(squeeze(mask(:,:,:,kk))))), S, 1, data, R, 'niter',Niter,'isave','all');
                    img(:,:,:,kk,ll,mm,nn,oo) = embed(test(:,end),logical(mask(:,:,:,kk)));
                    toc
                end
            end
        end
    end
end

end


