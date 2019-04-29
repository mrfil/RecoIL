function [ sen, mask ] = createSenMap(coilImages, lImagetype,bcImage)
%createSenMap - Creates a SENSE (sensitivity) map from coil images and
%                   recoInfo Object
%
% Syntax:  [sen] = createSenMap(recoInfo,coilImages)
%
% Inputs:
%    coilImages - coil Images of dimension [N,N,NSlices,NCoils]
%    lImageType - number representing type of masking to use
% Outputs:
%    sen  - senseMap of size [N,N,NSlices,NCoils]
%    mask - binary object mask of size [N,N,NSlices]
%
% Example:
%    rInfo         = recoInfo(filename);
%    images        = gridCoilImages(rInfo);
%    [sen, mask]   = createSenMap(images,1);
%   
%
% Other m-files required: recoInfo.m, solve_pwls_pcg.m, Robject.m
% Subfunctions: none
% MAT-files required: none
%
% Author:
% University of Illinois at Urbana-Champaign
% email address:
% Website:
% July 2016; Last revision: 1-Jul-2016


%% check input arguments
sizeDims = size(coilImages);

if length(sizeDims) ~= 4
    error('coilImages input must be N x N x NSlices x NCoils in size');
end

if sizeDims(1) ~= sizeDims(2)
    error('Image is not square. createSenMap requires square images inplane.');
end

%% Setup parameters

N = sizeDims(1);
NSlices = sizeDims(3);
NCoils = sizeDims(4);


%% create sum of squares and scaled coil images
if nargin < 3
    sos_img = sqrt(sum(abs(coilImages).^2,4));
else
    sos_img = bcImage;
end
scale_factor = max(abs(sos_img(:)));
sos_img = sos_img/scale_factor;
coil_imgs = coilImages/scale_factor;

[x, y] = meshgrid(1:N,1:N);
mask_circ = sqrt((x-(N+1)/2).^2+(y-(N+1)/2).^2)<((N)/2);
mask_circ = repmat(mask_circ,[1 1 NSlices]);
N_extra = 0;
switch lImagetype
    case 1
        %create a mask for estimate coil sensitivity
        
        %mask objcect based on histogram of values
        [hist_counts, hist_edges]= histcounts(sos_img(:));
        ii = 1;
        while hist_counts(ii)>hist_counts(ii+1)
            ii = ii + 1;
        end
        mask = (sos_img>hist_edges(ii)).*mask_circ;
        
        %make a mask to remove data weighting from voxels with large
        %phase differences
        %TODO: rewrite with linear algebra for simplicity
        [Nx, Ny, NSlices, NCoils] = size(coil_imgs);
        %             Rdiff = Robject(ones(size(mask(:,:,1))),'edge_type','tight','order',1,'beta',1,'type_denom','matlab','potential','quad');
        tmpMask = zeros(size(coil_imgs));
        phaseThreshhold = 0.2;
        for cc = 1:NCoils
            for xx = 1:Nx
                for yy = 1:Ny
                    for zz = 1:NSlices
                        if xx >1
                            tmpMask(xx,yy,zz,cc) = tmpMask(xx,yy,zz,cc) + (abs(angle(exp(1i*(angle(coil_imgs(xx,yy,zz,cc))-angle(coil_imgs(xx-1,yy,zz,cc))))))<phaseThreshhold);
                        end
                        if xx < Nx
                            tmpMask(xx,yy,zz,cc) = tmpMask(xx,yy,zz,cc) + (abs(angle(exp(1i*(angle(coil_imgs(xx,yy,zz,cc))-angle(coil_imgs(xx+1,yy,zz,cc))))))<phaseThreshhold);
                        end
                        if yy >1
                            tmpMask(xx,yy,zz,cc) = tmpMask(xx,yy,zz,cc) + (abs(angle(exp(1i*(angle(coil_imgs(xx,yy,zz,cc))-angle(coil_imgs(xx,yy-1,zz,cc))))))<phaseThreshhold);
                        end
                        if yy < Ny
                            tmpMask(xx,yy,zz,cc) = tmpMask(xx,yy,zz,cc) + (abs(angle(exp(1i*(angle(coil_imgs(xx,yy,zz,cc))-angle(coil_imgs(xx,yy+1,zz,cc))))))<phaseThreshhold);
                        end
                        if zz >1
                            tmpMask(xx,yy,zz,cc) = tmpMask(xx,yy,zz,cc) + (abs(angle(exp(1i*(angle(coil_imgs(xx,yy,zz,cc))-angle(coil_imgs(xx,yy,zz-1,cc))))))<phaseThreshhold);
                        end
                        if zz < NSlices
                            tmpMask(xx,yy,zz,cc) = tmpMask(xx,yy,zz,cc) + (abs(angle(exp(1i*(angle(coil_imgs(xx,yy,zz,cc))-angle(coil_imgs(xx,yy,zz+1,cc))))))<phaseThreshhold);
                        end
                    end
                end
            end
        end
        if NCoils == 4
            if NSlices > 4
                mask = mask.*(mean(tmpMask,4)>3.5);
            else
                mask = mask.*(mean(tmpMask,4)>2);
            end
        elseif NCoils == 32
            if NSlices >4
                mask = mask.*(mean(tmpMask,4)>2.5);
            else
                mask = mask.*(mean(tmpMask,4)>2);
            end
       elseif ((NCoils == 44) || (NCoils == 64))
            if NSlices >4
                mask = mask.*(mean(tmpMask,4)>2.5);
            else
                mask = mask.*(mean(tmpMask,4)>2);
            end
        else
            mask = mask.*(mean(tmpMask,4)>2);
        end
        mask = (mask > 0);
        
        %smoothing terms for SENSE maps
        beta = 1E-6; %AMC
        %beta = 2E6; %
        niter = 100; %AMC
        %niter = 10000; %AMC
        %             init_img = 1e-2*ones(N*N,1);
        %             Cn = C3D_sparse(ones(N,N))*beta;
    case 2
        %create a mask for estimate coil sensitivity
        mask = (abs(sos_img) > (0.1*max(abs(sos_img(:)))));
        
        mask = mask.*mask_circ;
        mask = (mask > 0);
        beta = 1;
        niter = 100;
    case 3
        %create a mask for estimate coil sensitivity
        mask = (abs(sos_img) > (0.1*max(abs(sos_img(:)))));
        
        mask = mask.*mask_circ;
        mask = (mask > 0);
        beta = 1e-3;
        niter = 100;

    otherwise
        %create a mask for estimate coil sensitivity
        mask = true(size(sos_img));
        Cn = 0;
        niter = 100;
        init_img = ones(N*N,1);
        beta = 1;
end

mask = mask_circ; % removing the stupid mask shit for the sense map

%mask coil images with a circular FOV
%for coilIndex= 1:NCoils
%    coil_imgs(:,:,:,coilIndex) = coil_imgs(:,:,:,coilIndex).*mask_circ;
%end

sen = zeros(N,N,NSlices,NCoils);
sen_extra = zeros(N+2*N_extra,N+2*N_extra,NSlices,NCoils);

%weight the edges to force the SENSE map to go to a certain value there
% edge_weight = mean(sos_img(:));
data_weight = zeros(size(sen_extra(:,:,:,1)));
data_weight((1:N)+N_extra,(1:N)+N_extra,:,1) = abs(sos_img).*mask;
% data_weight(1,:,:) = edge_weight;
% data_weight(end,:,:) = edge_weight;
% data_weight(:,1,:) = edge_weight;
% data_weight(:,end,:) = edge_weight;

%     W = Gdiag(col(sos_img.*mask+~(mask_circ)*mean(sos_img(:))));
if NSlices < 4
    %R = Robject(ones(size(sen_extra(:,:,1,1))),'edge_type','tight','order',2,'beta',beta,'type_denom','matlab','potential','quad');
    R = Robj(logical(ones(size(sen_extra(:,:,1,1)))),'beta',beta,'potential','quad');
    A = sensemult(ones(size(sen_extra(:,:,1,1))));
    init_img = zeros(size(sen_extra(:,:,1,1)));
    init_data = zeros(size(sen_extra(:,:,1,1)));
else
    %R = Robject(ones(size(sen_extra(:,:,:,1))),'edge_type','tight','order',2,'beta',beta,'type_denom','matlab','potential','quad');
    R = Robj(logical(ones(size(sen_extra(:,:,:,1)))),'beta',beta,'potential','quad');
    A = sensemult(ones(size(sen_extra(:,:,:,1))));
    W = Gdiag(data_weight(:));
    %         W = 1;
    %     W = Gdiag(col(sos_img));
    init_img = zeros(size(sen_extra(:,:,:,1)));
    init_data = zeros(size(sen_extra(:,:,:,1)));
end
for coilIndex = 1:NCoils
    if NSlices < 4
        for sliceIndex = 1:NSlices
            init_img((1:N)+N_extra,(1:N)+N_extra) = (coil_imgs(:,:,sliceIndex,coilIndex)./sos_img(:,:,sliceIndex).*mask(:,:,sliceIndex));
            init_data((1:N)+N_extra,(1:N)+N_extra) = (coil_imgs(:,:,sliceIndex,coilIndex)./sos_img(:,:,sliceIndex).*mask(:,:,sliceIndex));
            W = Gdiag(col(data_weight(:,:,sliceIndex)));
            tmp = solve_pwls_pcg(init_img(:), A,W ,init_data(:), R, 'niter',niter,'stepper',{'qs1',1});
            tmp = reshape(tmp,N+2*N_extra,N+2*N_extra);
            sen(:,:,sliceIndex,coilIndex) = tmp((1:N)+N_extra,(1:N)+N_extra,:);
        end
    else
        init_img((1:N)+N_extra,(1:N)+N_extra,:,1) = (coil_imgs(:,:,:,coilIndex)./sos_img.*mask);
        init_data((1:N)+N_extra,(1:N)+N_extra,:,1) = (coil_imgs(:,:,:,coilIndex)./sos_img.*mask);
        init_img(isnan(init_img)) = 0;
        init_data(isnan(init_data)) = 0;
        tmp = solve_pwls_pcg(init_img(:), A,W ,init_data(:), R, 'niter',niter,'stepper',{'qs1',1});
        tmp = reshape(tmp,N+2*N_extra,N+2*N_extra,NSlices);
        sen(:,:,:,coilIndex) = tmp((1:N)+N_extra,(1:N)+N_extra,:);
    end
end
%create an object support mask
%     mask = (abs(sos_img) > (0.1*max(abs(sos_img(:)))));
%mask objcect based on histogram of values
% [hist_counts, hist_edges]= histcounts(sos_img(:));
% ii = 1;
% while hist_counts(ii)>hist_counts(ii+1)
%     ii = ii + 1;
% end
%mask = (sos_img>hist_edges(ii)).*mask_circ;
%     mask = mask.*mask_circ;
%     se1 = strel('disk',3);
%     se2 = strel('disk',5);
%     for sliceIndex = 1:NSlices
%         mask(:,:,sliceIndex) = imdilate(mask(:,:,sliceIndex),se1);
%         mask(:,:,sliceIndex) = imerode(mask(:,:,sliceIndex),se2);
%         mask(:,:,sliceIndex) = imdilate(mask(:,:,sliceIndex),se2);
%     end
mask = (mask > 0);

end




