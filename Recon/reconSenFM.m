function [cImages, sen, FM, FMImages] = reconSenFM(varargin)
%reconSenFM - Creates SENSE Maps and Field Maps with default options
%
% Syntax:  [cImages, sen, FM, FMImages] = ReconSenFM()
%
% Inputs:
%    None
%
% Outputs:
%    cImages - Gridded Coil Images of "standard dimensions"
%    sen     - Sense Maps
%    FM      - Field Maps
%    FMImages- SENSE Reconstructed images used to compute the field map
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
% 11-Apr-2017; Last revision: 12-Sep-2017


%% Deal with optional input arguments
% Nothing to add at the current time.
p = inputParser();
TEMask = [];
UseGIRF = false;
PreWhiten = false;
addOptional(p,'UseGIRF',UseGIRF);
addOptional(p,'TEMask',TEMask);
addOptional(p,'PreWhiten',PreWhiten);
parse(p,varargin{:});
UseGIRF = p.Results.UseGIRF;
TEMask = p.Results.TEMask;
PreWhiten = p.Results.PreWhiten;




%% Moving on to actual reconstruction and estimation.
rInfo = recoInfo('','UseGIRF',UseGIRF,'PreWhiten',PreWhiten);

if isempty(TEMask)
    TEMask = logical(ones(1,rInfo.nEchoes));
end

if ~exist('cImages.mat','file')
    cImages = gridCoilImages(rInfo);
    save -v7.3 cImages.mat cImages
else
    load cImages.mat
end
    

%sizeDims = size(cImages);

% if length(sizeDims) == 4
%     cImages = permute(squeeze(cImages),[1 2 5 3 4]);
% end

if ~exist('sen.mat','file')
    if rInfo.nPartitions > 1
        %3D/Multislab style acquisition
        sen = zeros(rInfo.N, rInfo.N, rInfo.nPartitions, rInfo.nSlices,rInfo.nCoils);
        mask = zeros(rInfo.N, rInfo.N, rInfo.nPartitions, rInfo.nSlices,rInfo.nCoils);
        
        for ii = 1:rInfo.nSlices
            senseStack = reshape(cImages(:,:,:,ii,1,1,1,1,1,:),[rInfo.N,rInfo.N,rInfo.nPartitions,rInfo.nCoils]);
            [sen(:,:,:,ii,:), mask(:,:,:,ii,:)] = createSenMap(senseStack, 2);
        end
    else
        %2D Style Acquisition
        senseStack = reshape(cImages(:,:,1,:,1,1,1,1,1,:),[rInfo.N,rInfo.N,rInfo.nSlices,rInfo.nCoils]);
        [sen, mask] = createSenMap(senseStack, 2);
    end
    save -v7.3 sen.mat sen
    save mask.mat mask
else
    load sen.mat
    load mask.mat
end
[x, y] = meshgrid(1:rInfo.N,1:rInfo.N);
maskCirc = sqrt((x-(rInfo.N+1)/2).^2+(y-(rInfo.N+1)/2).^2)<((rInfo.N)/2);
maskCirc = repmat(maskCirc,[1 1 rInfo.nPartitions rInfo.nSlices]);
maskC = squeeze(maskCirc);
save maskC.mat maskC

if ~exist('FM.mat','file')
if rInfo.nPartitions > 1
    % 3D/Multislab style acquisition - To Do: make sure that createFieldMap
    % will process one slab at a time, pass index
    FM = zeros(rInfo.N, rInfo.N, rInfo.nPartitions, rInfo.nSlices);
    FMImages =  zeros(rInfo.N, rInfo.N, rInfo.nPartitions, rInfo.nSlices, rInfo.nEchoes);
    for ii = 1:rInfo.nSlices
        [FM(:,:,:,ii), FMImages(:,:,:,ii,:)] = createFieldMap(rInfo, sen(:,:,:,ii,:), maskCirc(:,:,:,ii,:), 1);
    end
else
    % FM = zeros(rInfo.N, rInfo.N, rInfo.nSlices);
    % FMImages = zeros(rInfo.N, rInfo.N, rInfo.nSlices, rInfo.nEchoes);
    [FM, FMImages,FMmask, FMsmoothed] = createFieldMap(rInfo, sen, squeeze(maskCirc(:,:,1,:)), 2, 'nIterations',1,'TEMask',TEMask);
end 
save FM.mat FM
save FMImages.mat FMImages
save FMsmoothed.mat FMsmoothed
else
    load FM.mat
    load FMImages.mat
    load FMsmoothed.mat
end

end

