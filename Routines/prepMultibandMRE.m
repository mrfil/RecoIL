function prepMultibandMRE(datadir,filename)
% function prepMultibandMRE(datadir,filename)
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
% Curtis Johnson - University of Delaware
% Alex Cerjanic - University of Illinois at Urbana-Champaign
% Nov 2017
%
% Note: Built out of the former scratch script "reconMultibandMRE" created
% during MB sequence/recon development. Cleaned and consolidated by CLJ.

curdir = pwd;

if ~isempty(datadir)
    cd(datadir)
end

cd senfm
rInfoSen = recoInfo;
if ~exist('FM.mat','file')
   reconSenFM;
end
load FM.mat
load FMImages.mat
load sen.mat
load mask.mat

if rInfoSen.nEchoes == 2 % currently a workaround to allow both Funai and PRELUDE versions
Fmask = FM~=0;
[~,FMx] = estimFieldMap_alt(rInfoSen,FMImages,Fmask);
FM = FMx;
end

cd ..

% parse data
rInfo = recoInfo;

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
FMImagesMB   = zeros(rInfo.N,rInfo.N,rInfo.multibandFactor,rInfo.nSlices);

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


% save sen, FM, and mask out formatted for MB
save -v7.3 sen.mat  senMB  senMBNav
save -v7.3 FM.mat   FMMB   FMMBNav
save -v7.3 mask.mat maskMB maskMBNav

% reconstruct navigator images
if ~exist('imgNav.mat','file')
    imgNav = fieldCorrectedNavRecon(rInfo,senMBNav,maskMBNav,FMMBNav,'Rbeta',100,'dims2penalize',[1,1,0],'Niter',10,'L',0); % CLJ: turned up the beta and iterations
    save -v7.3 imgNav.mat imgNav
elseif ~exist('imgNav','var')
    load imgNav.mat
else
    error('imgNav.mat exists but no imgNav variable found!');
end

% calculate phase error maps
if ~exist('PMaps.mat','file')
    PMaps = calcMRENavPhaseErrors(rInfo, imgNav);
    save -v7.3 PMaps.mat PMaps
elseif ~exist('PMaps','var')
    load PMaps.mat
else
    error('PMaps.mat exists but no PMaps variable found!');
end

% write data out to ISMRMRD format for PowerGrid
if isempty(filename)
    filename = 'data';
end

if exist(sprintf('%s.h5',filename),'file')
    eval(sprintf('!rm %s.h5',filename))
end

[rInfoCC,senCC] = compressCoils(rInfo,senMB,'energyLevel',0.95);

if rInfo.nCoils < 32
convertRecoInfoToIsmrmrd(sprintf('%s.h5',filename),rInfo,permute(senMB,[1 2 3 5 4]), FMMB, PMaps);
else
convertRecoInfoToIsmrmrd(sprintf('%s.h5',filename),rInfoCC,permute(senCC,[1 2 3 5 4]), FMMB, PMaps);
end

cd(curdir)
