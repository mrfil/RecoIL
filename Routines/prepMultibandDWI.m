function prepMultibandDWI(datadir,filename)
% function prepMultibandDWI(datadir,filename)
%
% Inputs (optional):
%
% - datadir:  string pointing to data directory
%             (if not entered, will prep current directory)
%
% - filename: string with filename for ISMRMRD data file
%             (if not entered, will use "data" as default)
%
% Wrapper function to prep multiband MRE data for reconstruction:
% 1) recons senfm and creates SENSE maps and field map
% 2) resamples and reshapes for MB
% 3) recons navigator images
% 4) calculates phase error maps
% 5) converts and saves as ismrmrd format for PowerGrid
%
% Authors:
% Alex Cerjanic - University of Illinois at Urbana-Champaign
% Curtis Johnson - University of Delaware
% Jun 2018
%
% Note: Built out of the former scratch script "reconMultibandMRE" created
% during MB sequence/recon development. Cleaned and consolidated by CLJ.







%% Proceed with preparing recon.
curdir = pwd;

if ~isempty(datadir)
    cd(datadir)
end

cd senfm
rInfoSen = recoInfo();
if (exist('FM.mat','file') ~= 2)
   reconSenFM;
end
load FM.mat
load FMImages.mat
load sen.mat
load mask.mat

if rInfoSen.nEchoes == 2 % currently a workaround to allow both Funai and PRELUDE versions
Fmask = FM~=0;
[~,FMx] = estimFieldMap_alt(rInfoSen,FMImages,Fmask);
%FMx = 2*pi*estimFieldMap(FMImages,Fmask,rInfoSen.TE/1E6);
FM = FMx;
end

cd ..

% parse data
rInfo = recoInfo();

if (exist('sen.mat','file') ~= 2)
    % Reshape the field map, mask, and sense maps
    senResampled    = complex(zeros(rInfo.N,rInfo.N,rInfo.nSlices*rInfo.multibandFactor,rInfo.nCoils));
    senResampledNav = complex(zeros(rInfo.NNav,rInfo.NNav,rInfo.nSlices*rInfo.multibandFactor,rInfo.nCoils));
    
    % resample the sense maps for image and nav data
    for ii = 1:rInfo.nCoils
        senResampled(:,:,:,ii)    = resampleMap(sen(:,:,:,ii),rInfo,rInfoSen);
        senResampledNav(:,:,:,ii) = resampleMapNav(sen(:,:,:,ii),rInfo,rInfoSen);
    end
    
    % resample the field maps for image and nav data
    FMResampled       = resampleMap(FM,rInfo,rInfoSen);
    FMResampledNav    = resampleMapNav(FM,rInfo,rInfoSen);
    FMImagesResampled = resampleMap(FMImages(:,:,:,2),rInfo,rInfoSen);
    
    % resample the masks for image and nav data
    maskResampled    = resampleMap(double(squeeze(mask)),rInfo,rInfoSen);
    maskResampled    = maskResampled > 0.5;
    maskResampledNav = resampleMapNav(double(squeeze(mask)),rInfo,rInfoSen);
    maskResampledNav = maskResampledNav > 0.5;
    
    % multiband pulse phase
    % from Wong ISMRM 2012, #2209
    lslice_phase = getMbRfPhases(rInfo.multibandFactor);

    
    FMMB         = zeros(rInfo.N,rInfo.N,rInfo.multibandFactor,rInfo.nSlices);
    %FMImagesMB   = zeros(rInfo.N,rInfo.N,rInfo.multibandFactor,rInfo.nSlices);
    
    FMMBNav      = zeros(rInfo.NNav,rInfo.NNav,rInfo.multibandFactor,rInfo.nSlices);
    
    senMB        = zeros(rInfo.N,rInfo.N,rInfo.multibandFactor,rInfo.nSlices,rInfo.nCoils);
    senMBNav     = zeros(rInfo.NNav,rInfo.NNav,rInfo.multibandFactor,rInfo.nSlices,rInfo.nCoils);
    
    maskMB       = ones(rInfo.N,rInfo.N,rInfo.multibandFactor,rInfo.nSlices);
    maskMBNav    = ones(rInfo.NNav,rInfo.NNav,rInfo.multibandFactor,rInfo.nSlices);
    
    % reorder everything for MB format and apply pulse phase
    sliceOrder   = 1:rInfo.nSlices*rInfo.multibandFactor;
    nn=0;
    for jj=1:rInfo.multibandFactor
        for ii=1:rInfo.nSlices
            nn=nn+1;
            FMMB(:,:,jj,ii) = FMResampled(:,:,sliceOrder(nn));
            FMMBNav(:,:,jj,ii) = FMResampledNav(:,:,sliceOrder(nn));
            FMImagesMB(:,:,jj,ii) = FMImagesResampled(:,:,sliceOrder(nn));
            
            maskMB(:,:,jj,ii) = maskResampled(:,:,sliceOrder(nn));
            maskMBNav(:,:,jj,ii) = maskResampledNav(:,:,sliceOrder(nn));
            for kk=1:rInfo.nCoils
                senMB(:,:,jj,ii,kk) = senResampled(:,:,sliceOrder(nn),kk).*exp(-1i*lslice_phase(jj));
                senMBNav(:,:,jj,ii,kk) = senResampledNav(:,:,sliceOrder(nn),kk).*exp(-1i*lslice_phase(jj));
            end
        end
    end
    
    % Flip Slices
     %senMB = flip(senMB,3);
      senMBNav = flip(senMBNav,3);
     %FMMB = flip(FMMB,3);
      FMMBNav = flip(FMMBNav,3);
     %FMImagesMB = flip(FMImagesMB,3);
      maskMBNav = flip(maskMBNav,3);
    
    % save sen, FM, and mask out formatted for MB
    save -v7.3 sen.mat  senMB  senMBNav
    save -v7.3 FM.mat   FMMB   FMMBNav  %FMImagesMB
    save -v7.3 mask.mat maskMB maskMBNav
else 
    load sen.mat;
    load FM.mat;
    load mask.mat;
end
rInfo.dataMaskNav(1:end/2) = false; % Drop spiral in part of nav
% reconstruct navigator images
if ~exist('imgNav.mat','file')
    [rInfoCC, senMBNavCC] = compressCoils(rInfo,senMBNav, 'energyLevel', 0.95);
    imgNav = fieldCorrectedNavRecon(rInfoCC,senMBNavCC,maskMBNav,FMMBNav,'Rbeta',200,'dims2penalize',[1,1,0],'Niter',10,'L',5); % CLJ: turned up the beta and iterations
    %imgNav = fieldCorrectedNavRecon(rInfo,senMBNav,maskMBNav,FMMBNav,'Rbeta',2000,'dims2penalize',[1,1,0],'Niter',20,'L',5); % CLJ: turned up the beta and iterations
    imgNav = flip(imgNav,3); % Deal with the fact that we flipped the sen map for the navigators only.
    save -v7.3 imgNav.mat imgNav
elseif ~exist('imgNav','var')
    load imgNav.mat
else
    error('imgNav.mat exists but no imgNav variable found!');
end

%% Now we need to reconstruct the first B_0 image to register that to the navigators
[rInfoCC, senMBCC] = compressCoils(rInfo, senMB, 'energyLevel',0.80); % Get these recons done quick since they are thrown away
imgB0 = fieldCorrectedRecon(rInfoCC, senMBCC, maskMB, FMMB,'Rbeta',200,'dims2penalize',[1,1,0],'Niter',20,'L',0,'averagesToRecon',1);

% Get interpolated navigators for registration
[~,imgNavFull] = calcDWINavPhaseErrors(rInfo, imgNav);

[optimizer,metric] = imregconfig('multimodal');

for ii = 1:rInfo.nSlices
    for jj = 1:rInfo.multibandFactor
        tform{jj,ii} = imregtform(abs(imgNavFull(:,:,jj,1,1,ii,1)),abs(imgB0(:,:,jj,ii,1)),'affine',optimizer,metric);
    end
end

% Now calculate the navigators using the transforms from the registration
[PMaps, imgNavFullReg, magErrors] = calcDWINavPhaseErrors(rInfo, imgNav, tform);
save -v7.3 imgNavFull.mat imgNavFull imgNavFullReg imgB0 
save -v7.3 PMaps.mat PMaps magErrors

% write data out to ISMRMRD format for PowerGrid
if isempty(filename)
    filename = 'data';
end

convertRecoInfoToIsmrmrd(sprintf('%s.h5',filename),rInfo,permute(senMB,[1 2 3 5 4]), FMMB, PMaps);

cd(curdir)
