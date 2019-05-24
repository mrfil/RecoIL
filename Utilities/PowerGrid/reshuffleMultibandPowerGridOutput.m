function [imgND] = reshuffleMultibandPowerGridOutput(img)
% reshuffleMultibandPOwerGridOutput(img)
% img - An N-D array with the raw output from mergePowerGridFileOutput


sizes = size(img);


nX = sizes(1);
nY = sizes(2);
nPartitions = sizes(3);
nSlices = sizes(4);
nRepetitions = sizes(5);
if length(sizes) < 6
    nAverages = 1;
    nEchoes = 1;
    nPhases = 1;
else
    nAverages = sizes(6);
end

if length(sizes) < 7
    nEchoes = 1;
    nPhases = 1;
else
    nEchoes = sizes(7);
end

if length(sizes) < 8
    nPhases = 1;
else 
    nPhases = sizes(8);
end


imgND = zeros(nX,nY,nPartitions*nSlices,nRepetitions, nAverages, nEchoes,nPhases);
%imgFlipped = flipdim(img,3);

for par = 1:nPartitions
    for slc = 1:nSlices
        for rep = 1:nRepetitions
            for avg = 1:nAverages
                for eco = 1:nEchoes
                    for phs = 1:nPhases
                        imgND(:,:,(par-1) * nSlices + slc,rep,avg,eco,phs) = img(:,:,par,slc,rep,avg,eco,phs);
                    end
                end
            end
        end
    end
end



end
