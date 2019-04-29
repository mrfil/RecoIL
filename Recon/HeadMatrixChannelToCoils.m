function [dataCoils,NoiseScan] = HeadMatrixChannelToCoils(ascconv,data,NoiseScan)

%% Figuring out the right order of the channels
reorderArray = zeros(1,12);

for ii = 1:12 %NCoils
    %Cluster Number (1:4 for 12Ch Head Matrix)
    ClusterNumber = str2num(ascconv.asCoilSelectMeas.asList(ii).sCoilElementID.tElement(2));
    switch ascconv.asCoilSelectMeas.asList(ii).sCoilElementID.tElement(3)
        case 'P'
            reorderArray((ClusterNumber-1)*3 + 1) = ii;
        case 'S'
            reorderArray((ClusterNumber-1)*3 + 2) = ii;
        case 'T'
            reorderArray((ClusterNumber-1)*3 + 3) = ii;
        otherwise
            error('CoilElementID in ascconv header does not match expected form. Are you using this function with a coil other than the 12 Ch Head Matrix?');
            keyboard
    end
    
end

data = data(:,:,:,reorderArray);
%% Deal with the noise scan to make sure that it is in coils and not channels after transformation
if (nargin > 2)
    NoiseScan = NoiseScan(reorderArray,:);
    
    for ii = 1:4 %Number of clusters of 3 channels
        %Left Channel
        NoiseScanCoils(:,(ii-1)*3+1) = 0.5.*(NoiseScan(:,(ii-1)*3+1) + sqrt(2).*NoiseScan(:,(ii-1)*3+2) + NoiseScan(:,(ii-1)*3+3));
        %Middle Channel
        NoiseScanCoils(:,(ii-1)*3+2) = 1j.*1./sqrt(2).*(NoiseScan(:,(ii-1)*3+1) - NoiseScan(:,(ii-1)*3+3));
        %Right Channel
        NoiseScanCoils(:,(ii-1)*3+3) = 0.5.*(sqrt(2).*NoiseScan(:,(ii-1)*3+2) - NoiseScan(:,(ii-1)*3+1) - NoiseScan(:,(ii-1)*3+3));
    end
    
end



%% Transforming from the Tim channels to coils
for ii = 1:4 %Number of clusters of 3 channels
    %Left Channel
    dataCoils(:,:,:,(ii-1)*3+1) = 0.5.*(data(:,:,:,(ii-1)*3+1) + sqrt(2).*data(:,:,:,(ii-1)*3+2) + data(:,:,:,(ii-1)*3+3));
    %Middle Channel
    dataCoils(:,:,:,(ii-1)*3+2) = 1j.*1./sqrt(2).*(data(:,:,:,(ii-1)*3+1) - data(:,:,:,(ii-1)*3+3));
    %Right Channel
    dataCoils(:,:,:,(ii-1)*3+3) = 0.5.*(sqrt(2).*data(:,:,:,(ii-1)*3+2) - data(:,:,:,(ii-1)*3+1) - data(:,:,:,(ii-1)*3+3));
end

%dataCoils = permute(dataCoils,[3,2,1,4]);
end