function [FM, FMImages, mask,FMsmoothed] = createFieldMap(rInfo, sen, mask, lImagetype, varargin)
%createFieldMap - Creates a Field (Inhomogeneity) Map from multiple echo
%                   time data and a sense map.
%
% Syntax:  [FM] = createFieldMap(rInfo, sen, mask, lImagetype);
%
% Inputs:
%    rInfo - rInfo object that has been initialized with a Siemens
%    VB/VD/VE .dat file that can be read via mapVBVD
%    sen - SENSE Map corresponding to the field map
%    mask - mask corresponding to one from SENSE map generation
%    lImagetype - enumerator value corresponding to type of masking used
%             1 - Default threshold of 0.1 of max with dialate and erodes
%             2 - Default threshold only of 0.1 of max
%             otherwise (-1) - no mask
% Outputs:
%    FMsmoothed - field map of size in rad/s [N,N,NSlices] smoothed outside
%                 of the support mask
%    FMImages - SENSE reconstructed images [N,N,NSlices,NTEs]
%    mask -       The mask used to estimate the field map
%    FM   -       The unsmoothed field maps as returned by
%                 mri_field_map_reg.
%
% Optional Inputs: 
%   nIterations - Control how many times the FMImages are iterated using
%                 field correction. 1 reconstructs the field map images
%                 without using field correction. >= 2 begins to use field
%                 correction to obtain the field map images.
%   TEMask      - Control which echoes are used to reconstruct the field
%                 map. Should be one dimensional array.
%
% Example:
%    rInfo          = rInfo(filename);
%    images         = gridCoilImages(rInfo);
%    [sen, mask]    = createSenMap(images,1);
%    [FM, FMImages] = createFieldMap(rInfo,sen,mask,1);
%
% Other m-files required: rInfo.m, solve_pwls_pcg.m, Robject.m
% Subfunctions: none
% MAT-files required: none
%
% Author:
% University of Illinois at Urbana-Champaign
% email address:
% Website:
% July 2016; Last revision: 1-Jul-2016

%% Deal with input parsing using inputParser Class
p = inputParser;
nIterations = 4;
nIterationsValidationFcn = @(x) isscalar(x) && (x>0);
TEMask = logical(ones(rInfo.nEchoes,1));
p.addOptional('TEMask',TEMask);
p.addOptional('SliceIndex',-1,@isscalar);
p.addOptional('nIterations',nIterations,nIterationsValidationFcn);
p.parse(varargin{:});

inputs = p.Results;

TEMask = p.Results.TEMask;

if isempty(p.Results.SliceIndex)
    SliceIndex = -1;
else
    SliceIndex = inputs.SliceIndex;
end

nIterations = p.Results.nIterations;

%% create field maps
% Calculate a support mask for the field map knowing that the FOV that we
% have in spiral is radial in nature.
[x, y] = meshgrid(1:rInfo.N,1:rInfo.N);
maskCirc = sqrt((x-(rInfo.N+1)/2).^2+(y-(rInfo.N+1)/2).^2)<((rInfo.N)/2);

if rInfo.multibandFactor > 1
    maskCirc = repmat(maskCirc,[1 1 rInfo.multibandFactor rInfo.nSlices]);
else
    maskCirc = repmat(maskCirc,[1 1 rInfo.nPartitions rInfo.nSlices]);
end

FM = zeros(rInfo.N, rInfo.N, rInfo.nPartitions, rInfo.nSlices);

% Coil compress data to speed calculation of field map images
[rInfoCC, senCC] = compressCoils(rInfo,sen,'energyLevel',0.80);

niterImage = 10;
L = ceil(max(rInfo.timingVec(:)/1E-3));
echoes = 1:rInfo.nEchoes;
echoes = echoes(TEMask);
timingVec = repmat(rInfo.timingVec(:,:,:,:,:,:,1,:),[1,1,1,1,1,1,rInfo.nEchoes,1]);
for kk = 1:nIterations % Repetition for field corrected images
    FMImages = fieldCorrectedRecon(rInfoCC,reshape(senCC,rInfoCC.N,rInfoCC.N,1,rInfoCC.nSlices,rInfoCC.nCoils),...
        reshape(maskCirc,rInfoCC.N,rInfoCC.N,1,rInfoCC.nSlices),...
        reshape(FM,rInfoCC.N,rInfoCC.N,1,rInfoCC.nSlices),...
        'Rbeta',1E-3,...
        'L',L,...
        'niter',niterImage,...
        'echoesToRecon',echoes,...
        'timingVec',timingVec);
    FMImages = reshape(FMImages,rInfo.N,rInfo.N,rInfo.nSlices,[]);
    TEs = rInfo.TE;
    
    switch lImagetype
        case 1
            mask = maskCirc;
        case 2
            % Use FSL bet to calculate a nice brain support mask for the
            % field map. Skip the skull completely, only for field map
            % estimation.
            %mask = double(bet(abs(FMImages(:,:,:,1)),'f',.4));
            mask = double(bet(abs(FMImages(:,:,:,end)),'f',.5));
            %mask = (abs(FMImages(:,:,:,1)) > (0.1*max(abs(col(FMImages(:,:,:,1))))));
            %lastTEIndex = find(TEMask,1,'last');
            %maskSmoothing = (col(abs(FMImages(:,:,:,lastTEIndex))) > ...
            %                            (0.15*max(col(abs(FMImages(:,:,:,lastTEIndex))))));
            %maskSmoothing = reshape(maskSmoothing,rInfo.N,rInfo.N,rInfo.nPartitions,rInfo.nSlices);
            %mask = maskSmoothing.*maskCirc;
            
            % Erode the field map mask to eliminate spurious flow voxels
%             nhood = strel('disk',5);
%             for ii = 1:rInfo.nSlices
%                 mask(:,:,ii) = imfill(imerode(mask(:,:,ii),nhood));
%             end
        otherwise
            % Do nothing, use same mask as the one used to reconstruct the
            % field map images.
    end
    
    if rInfo.nSlices < 4
        FM = zeros(rInfo.N, rInfo.N, rInfo.nPartitions, rInfo.nSlices);
        FMConv = zeros(rInfo.N, rInfo.N, rInfo.nPartitions, rInfo.nSlices);
        for sliceIndex = 1:rInfo.nSlices
            %             keyboard
            [FM(:,:,sliceIndex), FMConv(:,:,sliceIndex)] = mri_field_map_reg(squeeze(FMImages(:,:,sliceIndex,TEMask)), TEs(TEMask)*1e-6,'l2b',-3,'mask',(mask(:,:,sliceIndex)>0));
        end
        FM = -1*FM;
        FMConv = -1*FMConv;
    else
        [FM, FMConv] = mri_field_map_reg(squeeze(FMImages(:,:,:,TEMask)), TEs(TEMask)*1e-6,'l2b',-3,'mask',squeeze((mask>0)));
        FM = -1*FM;
        FMConv = -1*FMConv;
    end
end
    %% Smooth field map over entire image, enforcing data fidelity only over the tight mask used in SENSE Map estimation
    beta = 1;
    niter = 1000;
    
    FMsmoothed = zeros(rInfo.N,rInfo.N,rInfo.nSlices);
    
    % Since FM necessarily sucks where there is no long TE signal, let's mask all
    % of that out and smooth over those points.
    lastTEIndex = find(TEMask,1,'last');
    maskSmoothing = (col(abs(FMImages(:,:,:,lastTEIndex))) > (0.1*max(col(abs(FMImages(:,:,:,lastTEIndex))))));
    maskSmoothing = reshape(maskSmoothing,rInfo.N,rInfo.N,rInfo.nSlices);
    data_weight = double(maskSmoothing);
    
    
    %     W = Gdiag(col(sos_img.*mask+~(mask_circ)*mean(sos_img(:))));
    %if rInfo.nSlices < 4
    %%R = Robject(ones(size(FMsmoothed(:,:,1))),'edge_type','tight','order',2,'beta',beta,'type_denom','matlab','potential','quad');
    R = Robj(logical(ones(size(FMsmoothed(:,:,1)))),'beta',beta,'potential','quad');
    A = sensemult(ones(size(FMsmoothed(:,:,1))));
    %init_img = zeros(size(FMsmoothed(:,:,1)));
    %       %init_data = zeros(size(FMsmoothed(:,:,1)));
    %    else
    %       R = Robject(ones(size(FMsmoothed(:,:,:))),'edge_type','tight','order',2,'beta',beta,'type_denom','matlab','potential','quad');
    %       A = sensemult(ones(size(FMsmoothed(:,:,:))));
    %       W = Gdiag(data_weight(:));
    %       %         W = 1;
    %       %     W = Gdiag(col(sos_img));
    %       %init_img = zeros(size(FMsmoothed(:,:,:)));
    %       %init_data = zeros(size(FMsmoothed(:,:,:)));
    %    end
    
    %    if rInfo.nSlices < 4
    for sliceIndex = 1:rInfo.nSlices
        init_img = FM(:,:,sliceIndex).*mask(:,:,sliceIndex);
        init_data = FM(:,:,sliceIndex).*mask(:,:,sliceIndex);
        W = Gdiag(col(data_weight(:,:,sliceIndex)));
        tmp = solve_pwls_pcg(init_img(:), A,W ,init_data(:), R, 'niter',niter,'stepper',{'qs1',1});
        tmp = reshape(tmp,rInfo.N,rInfo.N);
        FMsmoothed(:,:,sliceIndex) = tmp;
    end
    %    else
    %       init_img = FM.*mask;
    %       init_data = FM.*mask;
    %       init_img(isnan(init_img)) = 0;
    %       init_data(isnan(init_data)) = 0;
    %       tmp = solve_pwls_pcg(init_img(:), A,W ,init_data(:), R, 'niter',niter,'stepper',{'qs1',1});
    %       tmp = reshape(tmp,rInfo.N,rInfo.N,rInfo.nSlices);
    %       FMsmoothed(:,:,:) = tmp;
    %    end

end


