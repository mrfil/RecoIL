function img = gridCoilImages(rInfo, varargin)
%gridCoilImages - Grids coil images from data read in via mapVBVD and rInfo
%
% Syntax:  [img] = gridCoilImages(rInfo)
%
% Inputs:
%    rInfo - Object of class rInfo that has been initialized with a
%    Siemens VB/VD/VE *.dat file
%
% Optional Name - Value Pair Inputs:
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
%    img - Gridded Coil Images of size [N,N,Kz,NSlices (NSlabs),NCoils,
%                                         NEchos]
%
%
% Example:
%    rInfo  = rInfo(filename);
%    images = gridCoilImages(rInfo);
%    im(sum(abs(images(:,:,:,:,:,1)).^2,4)); % Sum of Squares image for
%                                            % first echo time in field map
%
%
% Other m-files required: rInfo.m, k2image.m, mapVBVD
% Subfunctions: none
% MAT-files required: none
%
% Author:
% University of Illinois at Urbana-Champaign
% email address:
% Website:
% July 2016; Last revision: 18-Mar-2018


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
addOptional(p,'slicesToRecon', slicesToRecon, slicesToReconValidationFcn);
addOptional(p,'averagesToRecon', averagesToRecon, averagesToReconValidationFcn);
addOptional(p,'phasesToRecon', phasesToRecon, phasesToReconValidationFcn);
addOptional(p,'echoesToRecon', echoesToRecon, echoesToReconValidationFcn);
addOptional(p,'repetitionsToRecon', repetitionsToRecon, repetitionsToReconValidationFcn);

parse(p,varargin{:});

slicesToRecon = p.Results.slicesToRecon;
averagesToRecon = p.Results.averagesToRecon;
phasesToRecon = p.Results.phasesToRecon;
echoesToRecon = p.Results.echoesToRecon;
repetitionsToRecon = p.Results.repetitionsToRecon;

if rInfo.multibandFactor > 1 % Multiband Case
    img = zeros(rInfo.N,rInfo.N,rInfo.multibandFactor,rInfo.nCoils,length(slicesToRecon),length(averagesToRecon),length(phasesToRecon),length(echoesToRecon),1,length(repetitionsToRecon));
else % Standard 2D or 3D case
	img = zeros(rInfo.N,rInfo.N,rInfo.nPartitions,rInfo.nCoils,length(slicesToRecon),length(averagesToRecon),length(phasesToRecon),length(echoesToRecon),1,length(repetitionsToRecon));
end

% Reconstruct the images
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
                    data = reshape(rInfo.dataRead([],[],slc,avg,phs,eco,rep,1),[],rInfo.nShots,rInfo.nPartitions,rInfo.nCoils); %read all data and phase for off center slices
                    if rInfo.multibandFactor > 1 % Multiband Case
                        G = NUFFT(col(rInfo.kRead(rInfo.dataMask,:,:,slc,avg,phs,eco,rep,1)),col(rInfo.kPhase(rInfo.dataMask,:,:,slc,avg,phs,eco,rep,1)),col(rInfo.kSlice(rInfo.dataMask,:,:,slc,avg,phs,eco,rep,1)),rInfo.N,rInfo.N,rInfo.multibandFactor);
                    elseif rInfo.nPartitions > 1 % 3D Case
                        %data = reshape(data,[],rInfo.nPartitions,rInfo.nSlices,rInfo.nCoils,rInfo.nEchoes); %Put the kz-encode dimension in the second dimension for tranforms
                        %data = fftshift(ifft(fftshift(data,3),[],3),3);
                        G = NUFFT(col(rInfo.kRead(rInfo.dataMask,:,:,slc,avg,phs,eco,rep,1)),col(rInfo.kPhase(rInfo.dataMask,:,:,slc,avg,phs,eco,rep,1)),col(rInfo.kSlice(rInfo.dataMask,:,:,slc,avg,phs,eco,rep,1)),rInfo.N,rInfo.N,rInfo.nPartitions);
                    else % 2D Case
                        %curData = col(data(:,:,jj,1,1,1,1,1,1,qq));
                        %img(:,:,jj,slc,ll,mm,nn,oo,pp,qq) = k2image(col(rInfo.kx(:,:,jj,slc,ll,mm,nn,oo,pp)),col(rInfo.ky(:,:,jj,slc,ll,mm,nn,oo,pp)),curData,col(rInfo.ww(:,:,jj,slc,ll,mm,nn,oo,pp)),rInfo.N,6);
                        G = NUFFT(col(rInfo.kRead(rInfo.dataMask,:,1,slc,avg,phs,eco,rep,1)),col(rInfo.kPhase(rInfo.dataMask,:,1,slc,avg,phs,eco,rep,1)),col(rInfo.kSlice(rInfo.dataMask,:,1,slc,avg,phs,eco,rep,1)),rInfo.N,rInfo.N,1);
                    end
                    for qq = 1:rInfo.nCoils
                        if rInfo.multibandFactor > 1 % Multiband Case
                            temp = G'*(col(rInfo.ww(rInfo.dataMask,:,:,slc,avg,phs,eco,rep,1)).*col(data(rInfo.dataMask,:,:,qq)));
                            img(:,:,:,qq,kk,ll,mm,nn,oo,1) = reshape(temp,rInfo.N,rInfo.N,rInfo.multibandFactor);
                        elseif rInfo.nPartitions > 1
                            temp = G'*(col(rInfo.ww(rInfo.dataMask,:,:,slc,avg,phs,eco,rep,1)).*col(data(rInfo.dataMask,:,:,qq)));
                            img(:,:,:,qq,kk,ll,mm,nn,oo,1) = reshape(temp,rInfo.N,rInfo.N,rInfo.nPartitions);
                        else
                            temp = G'*(col(rInfo.ww(rInfo.dataMask,:,:,slc,avg,phs,eco,rep,1)).*col(data(rInfo.dataMask,:,qq)));
                            img(:,:,:,qq,kk,ll,mm,nn,oo,1) = reshape(temp,rInfo.N,rInfo.N,1);
                        end
                    end
                end
            end
        end
    end
end

img = permute(img,[1,2,3,5,6,7,8,9,10,4]);

end

