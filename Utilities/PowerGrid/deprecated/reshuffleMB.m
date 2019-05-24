function [imgND] = reshuffleMB(img)
% 
sizes = size(img);


nX = sizes(1);
nY = sizes(2);
nPartitions = sizes(3);
nSlices = sizes(4);
nAverages = sizes(5);

imgND = zeros(nX,nY,nPartitions*nSlices,nAverages);
%imgFlipped = flipdim(img,3);

for kk = 1:nPartitions
    for ll = 1:nSlices
        for mm = 1:nAverages
            imgND(:,:,kk * nSlices + ll,mm) = img(:,:,kk,ll,mm);
        end
    end
end

%imgND = flipdim(imgND,3);

end
