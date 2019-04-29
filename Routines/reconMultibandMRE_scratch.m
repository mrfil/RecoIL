% wrapper function for recon multiband MRE data
% To Dos: Implement pcSENSE recon

cd senfm
rInfoSen = recoInfo();
if ~exist('FM.mat','file')
   reconSenFM;
end
load FM.mat
load FMImages.mat
load sen.mat
load mask.mat

cd ..

% parse data
rInfo = recoInfo();

% Reshape the field map, mask, and sense maps
senResampled    = complex(zeros(rInfo.N,rInfo.N,rInfo.nSlices*rInfo.multibandFactor,rInfo.nCoils));
senNavResampled = complex(zeros(rInfo.NNav,rInfo.NNav,rInfo.nSlices*rInfo.multibandFactor,rInfo.nCoils));

for ii = 1:rInfo.nCoils
   % senResampled(:,:,:,ii)    = resample_map_resolution(sen(:,:,:,ii),rInfo.N,rInfo.nSlices*rInfo.multibandFactor,rInfoSen.FOV*10,rInfo.FOV*10);
   % senResampledNav(:,:,:,ii) = resample_map_resolution(sen(:,:,:,ii),rInfo.NNav,rInfo.nSlices*rInfo.multibandFactor,rInfoSen.FOV*10,rInfo.FOV*10);
   senResampled(:,:,:,ii) =    resampleMap(sen(:,:,:,ii),rInfo,rInfoSen);
   senResampledNav(:,:,:,ii) = resampleMapNav(sen(:,:,:,ii),rInfo,rInfoSen);
end

%FMResampled      = resample_map_resolution(FM,rInfo.N,rInfo.nSlices*rInfo.multibandFactor,rInfoSen.FOV*10,rInfo.FOV*10);
%FMResampledNav   = resample_map_resolution(FM,rInfo.NNav,rInfo.nSlices*rInfo.multibandFactor,rInfoSen.FOV*10,rInfo.FOV*10);
FMResampled      = resampleMap(FM,rInfo,rInfoSen);
FMResampledNav   = resampleMapNav(FM,rInfo,rInfoSen);
FMImagesResampled      = resampleMap(FMImages(:,:,:,2),rInfo,rInfoSen);

maskResampled    = resampleMap(double(squeeze(mask)),rInfo,rInfoSen);
maskResampled = maskResampled > 0.5;
maskResampledNav = resampleMapNav(double(squeeze(mask)),rInfo,rInfoSen);
maskResampledNav = maskResampledNav > 0.5;
% this is for multiband!! 
lslice_phase = getMbRfPhases(rInfo.multibandFactor);

% reshape sen and FM and mask
FMMB         = zeros(rInfo.N,rInfo.N,rInfo.multibandFactor,rInfo.nSlices);
FMImagesMB   = zeros(rInfo.N,rInfo.N,rInfo.multibandFactor,rInfo.nSlices);

FMMBNav      = zeros(rInfo.NNav,rInfo.NNav,rInfo.multibandFactor,rInfo.nSlices);

senMB        = zeros(rInfo.N,rInfo.N,rInfo.multibandFactor,rInfo.nSlices,rInfo.nCoils);
senMBNav     = zeros(rInfo.NNav,rInfo.NNav,rInfo.multibandFactor,rInfo.nSlices,rInfo.nCoils);

maskMB       = ones(rInfo.N,rInfo.N,rInfo.multibandFactor,rInfo.nSlices);
maskMBNav    = ones(rInfo.NNav,rInfo.NNav,rInfo.multibandFactor,rInfo.nSlices);
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
            %senMB(:,:,jj,ii,kk) = senResampled(:,:,nn,kk);
        end
    end
end

%senMB       = flip(senMB,3);
%maskMB      = flip(maskMB,3);
%FMMB        = flip(FMMB,3);
%FMImagesMB  = flip(FMImagesMB,3);

%senNav      = flip(senMBNav,3);
%maskNav     = flip(maskMBNav,3);
%FMNav       = flip(FMMBNav,3);

save sen.mat  senMB  senMBNav
save FM.mat   FMMB   FMMBNav
save mask.mat maskMB maskMBNav

% rInfo.dataMask = true(rInfo.shotLength,1);
% rInfo.dataMaskNav = true(rInfo.shotLengthNav,1);

if ~exist('imgNav.mat','file')
    imgNav = fieldCorrectedNavRecon(rInfo,senMBNav,maskMBNav,FMMBNav,'Rbeta',1,'dims2penalize',[1,1,0],'Niter',4);
    save imgNav.mat imgNav
end
if ~exist('imgNav','var')
    load imgNav.mat
end

% PMaps = zeros(rInfo.N,rInfo.N,rInfo.multibandFactor,rInfo.nShots,rInfo.nPartitions,rInfo.nSlices,rInfo.nAverages,rInfo.nPhases,rInfo.nEchoes,rInfo.nRepetitions);
% for pp = 1:rInfo.nRepetitions
%     for oo = 1:rInfo.nEchoes
%         for nn = 1:rInfo.nPhases
%             for mm = 1:rInfo.nAverages
%                 for kk = 1:rInfo.nSlices
%                     for jj = 1:rInfo.nPartitions
%                         for ii = 1:rInfo.nShots
%                             for ll = 1:rInfo.multibandFactor
%                                 PMaps(:,:,ll,ii,jj,kk,mm,nn,oo,pp) = resample_map_resolution(imgNav(:,:,ll,ii,jj,kk,mm,nn,oo,pp),rInfo.N,1,1,1);
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end

PMaps = calcMRENavPhaseErrors(rInfo, imgNav);                
save -v7.3 PMaps.mat PMaps

convertRecoInfoToIsmrmrd('TBI_1707_AP.h5',rInfo,permute(senMB,[1 2 3 5 4]), FMMB, PMaps);


load sen.mat
load FM.mat
load mask.mat
load PMaps.mat
%img = fieldCorrectedReconMB(rInfo, senMB, maskMB, FMMB,'Rbeta',1,'dims2penalize',[1,1,0]);
img = phaseCorrectedRecon(rInfo, senMB, maskMB, FMMB, PMaps,'Niter',15,'L',6,'Rbeta',1,'dims2penalize',[1,1,0]);
save img.mat img

