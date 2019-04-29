function [C, wt] = C_sparse(mask,dims2penalize)
%function X = C3D_sparse(mask)
% Creates a roughness penalty matrix for the mask image.
%
% 2015/10/29 JLH: modified for n-dimensional
        
    if ~exist('dims2penalize','var')
        dims2penalize = 1:ndims(mask);
    end

    sz_img = size((mask));
    dims = length(sz_img);
    
    C_ones = ones(sz_img);
    C_onesSp =spdiag(C_ones(:));
    
    C_tmp =[];
    wt = [];
    shift = 1;
    for ii = 1:dims
        C_shiftSp =circshift(C_onesSp,[shift 0]);
        C_tmp = cat(1,C_tmp,C_onesSp-C_shiftSp);
        % remove edge cases and masked voxels
        wt_array_discard = zeros(shift,1);
        wt_array_keep = ones(shift*(sz_img(ii)-1),1);
        wt_array_tmp = cat(1,wt_array_discard,wt_array_keep);
        shift = shift*sz_img(ii);
        wt_array = repmat(wt_array_tmp,numel(mask)/shift,1).*mask(:);
        if (max(ismember(dims2penalize,ii))<1)
            wt = wt*0;
        end
        wt = cat(1,wt,wt_array);
    end
    
    %only calculate for voxels withinth mask
    C = C_tmp(:,logical(mask(:)));
    %set all points that exceed across the max boundary equal to 0
    
    % AMC: This code fragment below is irresponsibly slow...
    %for ii = 1:length(C)
    %    if (double(sum(C(ii,:))~=0))
    %        C(ii,:)=0;
    %    end
    %end
    
    % Hopefully faster version using some logical indexing tricks...
    rowSum = sum(C,2);
    idxToZero = (rowSum ~= 0);
    C(idxToZero,:) = 0;
    
    

end

 function b = spdiag(a)
%function b = spdiag(a, options)
% create a sparse matrix with diagonal given by a
% option:
%	'nowarn'	do not warn about diag_sp
% caution: it may be faster to use my newer diag_sp() object instead.

a = a(:);
a = double(a); % trick: needed because matlab7 doesn't handle single sparse well
n = length(a);
b = sparse(1:n, 1:n, a, n, n, n);
 end
