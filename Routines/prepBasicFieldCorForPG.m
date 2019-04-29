function prepBasicFieldCorForPG(datadir,filename)
% function prepBasicFieldCorForPG(datadir,filename)
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
% April 2018
%

% Alex Cerjanic UIUC 2/2/2018

% Run from the directory with the data
if ~isempty(datadir)
    cd(datadir)
end

% Assume that our senfm data is in a subfolder of the current directory
% called ./senfm
cd senfm
rInfoSen = recoInfo();

%Skip rereconstructing the senfm if it exists
if ~exist('FM.mat','file')
   reconSenFM;
end

%Grab the data we need.
load FM.mat
load sen.mat
load mask.mat

%Get back to the data directory
cd ..

% parse data - we need a parser for the sequence for this to work
rInfo = recoInfo();

% Preallocate for the interpolated field map, mask, and sense maps 
if rInfoSen.nPartitions > 1
    senResampled    = complex(zeros(rInfo.N,rInfo.N,rInfo.nPartitions,rInfo.nSlices,rInfo.nCoils));
    for jj = 1:rInfoSen.nSlices
        for ii = 1:rInfo.nCoils
            senResampled(:,:,:,jj,ii) =    resampleMap(sen(:,:,:,jj,ii),rInfo,rInfoSen);
        end
    end
    FMResampled(:,:,jj,:)      = resampleMap(FM(:,:,jj,:),rInfo,rInfoSen);
%    maskResampled(:,:,jj,:) = resampleMap(double(squeeze(mask(:,:,jj,:))),rInfo,rInfoSen);

else
    senResampled    = complex(zeros(rInfo.N,rInfo.N,1,rInfo.nSlices,rInfo.nCoils));
    for ii = 1:rInfo.nCoils
        senResampled(:,:,1,:,ii) =    resampleMap(sen(:,:,:,ii),rInfo,rInfoSen);
    end
    FMResampled(:,:,1,:)  = resampleMap(FM,rInfo,rInfoSen);
%    maskResampled(:,:,1,:) = resampleMap(double(squeeze(mask)),rInfo,rInfoSen);

end

% Interpolate the sensemaps up to the full imaging resolution on a coil by
% coil basis


%maskResampled    =  maskResampled > 0.5;

convertRecoInfoToIsmrmrd(sprintf('%s.h5',filename),rInfo,permute(senResampled,[1 2 3 5 4]), FMResampled);

% 
% % If our data is 3D we need this annoying line. This is an open issue.
% if rInfo.nPartitions > 1
%     senResampled = flip(senResampled,3);
%     FMResampled = flip(FMResampled,3);
% elseif rInfo.multibandFactor > 1
%     FMResampled = flip(FMResampled,3);
%     senResampled = flip(senResampled,3);
% end

% This line seems silly at first glance, but it allows us to deal with
% 2D/3D/multislab/simultaneous multislice without needing a seperate
% version of the recon code for each case.

% Remember the convention is [Nx,Ny,N 3D encoded slices,N 2D Encoded
% Slices, N Coils]

% For the 2D case, N 3D encoded slices = 1
% For the 3D case, N 2D encoded slices = 1
% For the multislab case, neither N 2D slices nor N 3D slices = 1
% For the 3D SMS case, neither N 2D Slices nor N 3D slices = 1

%senResampled = permute(senResampled,[1,2,3,5,4]);

% Uncomment the following line if you want coil compression, and use the 
% compressed the rInfoCC and senResampledCC in the call to the recon
% function

%[rInfoCC, senResampledCC] = compressCoils(rInfo,senResampled,'coilRank',4);

%img = fieldCorrectedRecon(rInfo, senResampled, maskResampled, FMResampled,'Rbeta',20,'dims2penalize',[1,1,1]);



%[rInfoCC, senResampledCC] = compressCoils(rInfo,senResampled,'coilRank',4);
%img = fieldCorrectedRecon(rInfoCC, senResampledCC, maskResampled, FMResampled,'Rbeta',1,'dims2penalize',[1,1,1],'L',0);
