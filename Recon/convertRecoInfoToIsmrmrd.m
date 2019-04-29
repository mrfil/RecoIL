function convertRecoInfoToIsmrmrd(filename,rInfo,sen, FM, PhaseMaps, varargin)
%convertRecoInfoToIsmrmrd Adapter from recoInfo object to ismrmrd file
%
%
%
% for SMS MRE: sense maps need to be [ny nx nz ncoil nslice]

%% Perform input parsing

shotsToRecon = 1:rInfo.nShots;
shotsToReconValidationFcn = @(x) isnumeric(x);

slicesToRecon = 1:rInfo.nSlices;
slicesToReconValidationFcn = @(x) isnumeric(x);

averagesToRecon = 1:rInfo.nAverages;
averagesToReconValidationFcn = @(x) isnumeric(x);

phasesToRecon = 1:rInfo.nPhases;
phasesToReconValidationFcn = @(x) isnumeric(x);

echoesToRecon = 1:rInfo.nEchoes;
echoesToReconValidationFcn = @(x) isnumeric(x);

repetitionsToRecon = 1:rInfo.nRepetitions;
repetitionsToReconValidationFcn = @(x) isnumeric(x);

partitionsToRecon = 1:rInfo.nPartitions;
partitionsToReconValidationFcn = @(x) isnumeric(x);

p = inputParser();
% Create an empty ismrmrd dataset
addOptional(p,'shotsToRecon', shotsToRecon, shotsToReconValidationFcn);
addOptional(p,'slicesToRecon', slicesToRecon, slicesToReconValidationFcn);
addOptional(p,'averagesToRecon', averagesToRecon, averagesToReconValidationFcn);
addOptional(p,'phasesToRecon', phasesToRecon, phasesToReconValidationFcn);
addOptional(p,'echoesToRecon', echoesToRecon, echoesToReconValidationFcn);
addOptional(p,'repetitionsToRecon', repetitionsToRecon, repetitionsToReconValidationFcn);
addOptional(p,'partitionsToRecon', partitionsToRecon, partitionsToReconValidationFcn);

parse(p,varargin{:});

% Override default parameters with their successfully parsed values.
shotsToRecon = p.Results.shotsToRecon;
slicesToRecon = p.Results.slicesToRecon;
averagesToRecon = p.Results.averagesToRecon;
phasesToRecon = p.Results.phasesToRecon;
echoesToRecon = p.Results.echoesToRecon;
repetitionsToRecon = p.Results.repetitionsToRecon;


%% Handle inputs
if exist(filename,'file')p = inputParser();
    error(['File ' filename ' already exists.  Please remove first'])
end
dset = ismrmrd.Dataset(filename);

% load the coil sensitivities

% create NDArray with the data
SENSEMap = ismrmrd.NDArray(col(sen));

% Write to file
dset.appendArray('SENSEMap',SENSEMap);

% load the field map
FM = col(FM);

% create NDArray with the data
FieldMap = ismrmrd.NDArray(FM);

% Write to file
dset.appendArray('FieldMap',FieldMap);

% if exists append PMaps to file
if nargin > 4 && ~isempty(PhaseMaps)
    PMaps = ismrmrd.NDArray(col(PhaseMaps));
    dset.appendArray('PhaseMaps', PMaps);
end

% It is very slow to append one acquisitiliceson at a time, so we're going
% to append a block of acquisitions at a time.
% In this case, we'll do it one repetition at a time to show off this
% feature.  Each block has nY aquisitions
acqblock = ismrmrd.Acquisition(length(shotsToRecon));

% Set the header elemen%delete(filename);ts that don't change
acqblock.head.version(:) = 1;
acqblock.head.number_of_samples(:) = rInfo.shotLength;
acqblock.head.center_sample(:) = 0;
acqblock.head.active_channels(:) = rInfo.nCoils;
acqblock.head.read_dir  = repmat([1 0 0]',[1 length(shotsToRecon)]);
acqblock.head.phase_dir = repmat([0 1 0]',[1 length(shotsToRecon)]);
acqblock.head.slice_dir = repmat([0 0 1]',[1 length(shotsToRecon)]);
acqblock.head.trajectory_dimensions = repmat(4,[1 length(shotsToRecon)]);


% Loop over the acquisitions, set the header, set the data and append
% For things we aren't dealing with for now.
eco = 1;
scanCounter = 0;
for phsIdx = 1:length(phasesToRecon)
    for avgIdx = 1:length(averagesToRecon)
        for repIdx = 1:length(repetitionsToRecon)
            for slcIdx = 1:length(slicesToRecon)
                for parIdx = 1:length(partitionsToRecon)
                    for shotIdx = 1:length(shotsToRecon)
                        
                        phs  = phasesToRecon(phsIdx);
                        avg  = averagesToRecon(avgIdx);
                        rep  = repetitionsToRecon(repIdx);
                        slc  = slicesToRecon(slcIdx);
                        par  = partitionsToRecon(parIdx);
                        shot = shotsToRecon(shotIdx);
                        
                        % Set the header elements that change from acquisition to the next
                        % c-style counting
                        acqblock.head.scan_counter(shotIdx) = scanCounter;
                        scanCounter = scanCounter + 1;
                        acqblock.head.idx.kspace_encode_step_1(shotIdx) = shotIdx-1;
                        acqblock.head.idx.kspace_encode_step_2(shotIdx) = parIdx-1;
                        acqblock.head.idx.repetition(shotIdx) = repIdx - 1;
                        acqblock.head.idx.slice(shotIdx) = slcIdx - 1;
                        acqblock.head.idx.average(shotIdx) = avgIdx - 1;
                        acqblock.head.idx.phase(shotIdx) = phsIdx - 1;
                        
                        % Set the flags
                        acqblock.head.clearAllFlags(shotIdx);
                        
                        if shotIdx == 1
                            acqblock.head.setFlag('ACQ_FIRST_IN_ENCODE_STEP1', shotIdx);
                        end
                        
                        if shotIdx == length(shotsToRecon)
                            acqblock.head.setFlag('ACQ_LAST_IN_ENCODE_STEP1', shotIdx);
                        end
                        
                        if parIdx == 1 && shotIdx == 1
                            acqblock.head.setFlag('ACQ_FIRST_IN_ENCODE_STEP2', shotIdx);
                        end
                        
                        if parIdx == length(partitionsToRecon) && shotIdx == length(shotsToRecon)
                            acqblock.head.setFlag('ACQ_LAST_IN_ENCODE_STEP2', shotIdx);
                        end
                        
                        if slcIdx == 1 && shotIdx == 1 && parIdx == 1
                            acqblock.head.setFlag('ACQ_FIRST_IN_SLICE', shotIdx);
                        end
                        
                        if slcIdx == rInfo.nSlices && parIdx == length(partitionsToRecon) && shotIdx == length(shotsToRecon)
                            acqblock.head.setFlag('ACQ_LAST_IN_SLICE', shotIdx);
                        end
                        
                        if slcIdx == 1 && shotIdx == 1 && parIdx == 1 && repIdx == 1
                            acqblock.head.setFlag('ACQ_FIRST_IN_REPETITION', shotIdx);
                        end
                        
                        
                        if slcIdx == rInfo.nSlices && parIdx == length(partitionsToRecon) && shotIdx == length(shotsToRecon) ...
                                && rep == length(repetitionsToRecon)
                            acqblock.head.setFlag('ACQ_LAST_IN_REPETITION', shotIdx);
                        end
                        
                        if slcIdx == 1 && shotIdx == 1 && parIdx == 1 && repIdx == 1 && avgIdx == 1
                            acqblock.head.setFlag('ACQ_FIRST_IN_AVERAGE', shotIdx);
                        end
                        
                        
                        if slcIdx == rInfo.nSlices && parIdx == length(partitionsToRecon) && shotIdx == length(shotsToRecon) ...
                                && rep == length(repetitionsToRecon) && avg == length(averagesToRecon)
                            acqblock.head.setFlag('ACQ_LAST_IN_AVERAGE', shotIdx);
                        end
                        
                        if slcIdx == 1 && shotIdx == 1 && parIdx == 1 && repIdx == 1 && avgIdx == 1 && phsIdx == 1
                            acqblock.head.setFlag('ACQ_FIRST_IN_PHASE', shotIdx);
                        end
                        
                        
                        if slcIdx == rInfo.nSlices && par == length(partitionsToRecon) && shot == length(shotsToRecon) && ...
                                rep == length(repetitionsToRecon) && avg == length(averagesToRecon) ...
                                && phs == length(phasesToRecon)
                            acqblock.head.setFlag('ACQ_LAST_IN_PHASE', shotIdx);
                        end
                        
                        
                        % fill the data
                        %acqblock.data{shot} = squeeze(data(:,shot,:,rep));
                        acqblock.data{shotIdx} = squeeze(rInfo.dataRead(shot,par,slc,avg,phs,eco,rep));
                        % attach the trajectory
                        %acqblock.traj{shot} = [kx(nR0*(acqno-1)+1:nR0*acqno),ky(nR0*(acqno-1)+1:nR0*acqno),kz(nR0*(acqno-1)+1:nR0*acqno),t(nR0*(acqno-1)+1:nR0*acqno)].';
                        acqblock.traj{shotIdx} = [col(rInfo.kRead(rInfo.dataMask,shot,par,slc,avg,phs,eco,rep)), ...
                            col(rInfo.kPhase(rInfo.dataMask,shot,par,slc,avg,phs,eco,rep)), ...
                            col(rInfo.kSlice(rInfo.dataMask,shot,par,slc,avg,phs,eco,rep)), ...
                            col(rInfo.timingVec(rInfo.dataMask,shot,par,slc,avg,phs,eco,rep))].';
                        
                        
                        
                    end % shot loop
                    % Append the acquisition block
                    dset.appendAcquisition(acqblock);
                end % partition loop
            end % slice loop
        end % rep loop
    end % avg loop
end % phs loop
%% Fill the xml header
% We create a matlab struct and then serialize it to xml.
% Look at the xml schema to see what the field names should be

header = [];

% Experimental Conditions (Required)
header.experimentalConditions.H1resonanceFrequency_Hz = ...
    str2double(rInfo.DataObject.hdr.Dicom.DICOM.lFrequency);

% Acquisition System Information (Optional)
header.acquisitionSystemInformation.systemVendor = rInfo.scannerManufacturer;
header.acquisitionSystemInformation.systemModel = rInfo.scannerModel;
header.acquisitionSystemInformation.receiverChannels = rInfo.nCoils;

% The Encoding (Required)
header.encoding.trajectory = 'spiral'; % Probably a safe bet for our lab
header.encoding.encodedSpace.fieldOfView_mm.x = rInfo.FOV;
header.encoding.encodedSpace.fieldOfView_mm.y = rInfo.FOV;
if rInfo.multibandFactor > 1 % 3D SMS case
    header.encoding.encodedSpace.fieldOfView_mm.z = rInfo.sliceThickness*rInfo.multibandFactor;
elseif rInfo.nPartitions > 1 % 3D case
    header.encoding.encodedSpace.fieldOfView_mm.z = rInfo.sliceThickness*rInfo.nPartitions;
else % 2D case
    header.encoding.encodedSpace.fieldOfView_mm.z = rInfo.sliceThickness;
end

header.encoding.encodedSpace.matrixSize.x = rInfo.N;
header.encoding.encodedSpace.matrixSize.y = rInfo.N;

if rInfo.multibandFactor > 1
    header.encoding.encodedSpace.matrixSize.z = rInfo.multibandFactor;
elseif rInfo.nPartitions > 1
    header.encoding.encodedSpace.matrixSize.z = rInfo.nPartitions;
else
    header.encoding.encodedSpace.matrixSize.z = 1;
end

% Recon Space
% (in this case same as encoding space)
header.encoding.reconSpace = header.encoding.encodedSpace;
% Encoding Limits;

% Step 0 is shots for us
header.encoding.encodingLimits.kspace_encoding_step_1.minimum = 0;
header.encoding.encodingLimits.kspace_encoding_step_1.maximum = length(shotsToRecon)-1;
header.encoding.encodingLimits.kspace_encoding_step_1.center = floor((length(shotsToRecon)-1)/2);
header.encoding.encodingLimits.kspace_encoding_step_2.minimum = 0;

% Step 1 is partitions (kz encodes)
% if rInfo.multibandFactor > 1 % 3D SMS case
%    header.encoding.encodingLimits.kspace_encoding_step_2.maximum = rInfo.multibandFactor-1;
%    header.encoding.encodingLimits.kspace_encoding_step_2.center = floor((rInfo.multibandFactor-1)/2);
% elseif rInfo.nPartitions > 1 % 3D Encoded case
if rInfo.nPartitions > 1 % 3D Encoded case
    header.encoding.encodingLimits.kspace_encoding_step_2.maximum = length(partitionsToRecon)-1;
    header.encoding.encodingLimits.kspace_encoding_step_2.center = floor((length(partitionsToRecon)-1)/2);
else % 2D Case
    header.encoding.encodingLimits.kspace_encoding_step_2.maximum = 0;
    header.encoding.encodingLimits.kspace_encoding_step_2.center = 0;
end

% Deal with Slices
header.encoding.encodingLimits.slice.minimum = 0;
header.encoding.encodingLimits.slice.maximum = length(slicesToRecon)-1;
header.encoding.encodingLimits.slice.center = floor((length(slicesToRecon)-1)/2);

% Deal with Repetitions
header.encoding.encodingLimits.repetition.minimum = 0;
header.encoding.encodingLimits.repetition.maximum = length(repetitionsToRecon)-1;
header.encoding.encodingLimits.repetition.center = floor((length(repetitionsToRecon)-1)/2);

% Deal with Averages
header.encoding.encodingLimits.average.minimum = 0;
header.encoding.encodingLimits.average.maximum = length(averagesToRecon)-1;
header.encoding.encodingLimits.average.center = floor((length(averagesToRecon)-1)/2);

% Deal with Phases
header.encoding.encodingLimits.phase.minimum = 0;
header.encoding.encodingLimits.phase.maximum = length(phasesToRecon)-1;
header.encoding.encodingLimits.phase.center = floor((length(phasesToRecon)-1)/2);

% Add multiband factor for convenience
header.encoding.encodedSpace.multibandFactor = rInfo.multibandFactor;

%% Serialize and write to the data set
xmlstring = ismrmrd.xml.serialize(header);
dset.writexml(xmlstring);

%% Write the dataset
dset.close();

end

