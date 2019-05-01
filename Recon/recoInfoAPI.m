classdef recoInfoAPI
    %recoInfo - Class definining API for reconstruction in RecoIL
    %
    % Syntax:  [rInfo] = recoInfo____
    %          [rInfo] = recoInfoISMRMRD(filename);
    %          [rInfo] = recoInfoTwix(filename);
    %          [rInfo] = recoInfoSim(kspace,data,timingVec));
    %
    % Inputs:
    %    filename - Name (and path) of raw data file
    %
    % Outputs:
    %    rInfo - Initalized recoInfo object created by a concrete base
    %            class such as recoInfoTwix or recoInfoIsmrmrd
    %
    %
    % Example:
    %    rInfo  = recoInfo(filename);
    %    images = gridCoilImages(rInfo);
    %    im(sum(abs(images(:,:,:,:,:,1)).^2,4)); % Sum of Squares image for
    %                                            % first echo time in field map
    %
    %
    % Other m-files required: mapVBVD
    % Subfunctions: none
    % MAT-files required: none
    %
    % Author:
    % University of Illinois at Urbana-Champaign
    % email address:
    % Website:
    % July 2016; Last revision: 27-Apr-2019
    
    properties
        % Siemens raw .dat file related parameters
        filename;
        DataObject;
        

        % Information that helps when using ISMRMRD
        scannerSerial = '';
        scannerManufacturer = '';
        scannerModel = '';
        scannerSoftwareVer = '';
        
        
        % Trajectory design related parameters
        gradAmp;
        gradSlew;
        gradTs;
        adcTs;
        adcTsNav;
        kRead;
        kPhase;
        kSlice;
        kReadNominal;
        kPhaseNominal;
        kSliceNominal;
        kReadNav;
        kPhaseNav;
        kSliceNav;
        kReadNominalNav;
        kPhaseNominalNav;
        kSliceNominalNav;
        nPartitionsNav; %3D partitions for navigator
        ww; % Density compensation vector
        timingVec;
        timingVecNav;
        traj_type;
        nShotsDesigned;
        nShotsUsed; % For sense reconstruction, can use different number of shots for recon
        shotLength; % Length of a single shot
        shotLengthNav;
        dataMask; % A mask that represents which samples to use in the reconstruction
        dataMaskNav;
        
        % Spatial Encoding Volume Parameters
        readShift;
        phaseShift;
        sliceShift;
        rotMatrix = eye(3);
        FOV;
        FOVSlice;
        sliceThickness;
        ExcOrder;
        multibandFactor = 1;
        
        % MDH Related Parameters
        N; % Mapped from lBaseResolution from Header
        NNav;
        nSlices; % Mapped from Slices in MDH
        nCoils; % Mapped from Channels in MDH
        nCoilsActual; % Number of coils in data (independent of coil comp)
        nADCs; % Mapped from Set in MDH
        nPartitions; % Mapped from Partitions to Kz in MDH
        
        nEchoes; % Mapped from Echoes in MDH
        nShots; % Mapped from Lines in MDH
        nRepetitions;
        nAverages; % Mapped from Acquisitions in MDH
        nPhases; % Mapped from Phasese in MDH
        halfVoxelShiftSlice; % HalfVoxelShift in Slice direction
        halfVoxelShiftPhase; % HalfVoxelShift in Phase direction
        halfVoxelShiftRead; % HalfVoxelShift in Read direction
        halfVoxelShiftNavSlice; % HalfVoxelShift for Navigators in Slice direction
        halfVoxelShiftNavPhase; % HalfVoxelShift for Navigators in Phase direction
        halfVoxelShiftNavRead; % HalfVoxelShift for Navigators in Read direction
        % Contrast Related Parameters
        TR;
        TE;
        FAdeg;
        nTimePoints;
        
        % Sampling related parameters
        ptsToDrop;  %Siemens scanners have to discard first points of ADC
        ptsToDropNav;
        nRO;
        nROa;
        lgadc;
        noiseCorr = 1;
        
        % Reconstruction related parameters
        sen = [];
        FM = [];
        mask = []; % Optional object support mask;
        useGIRF;
        nIter = 10;
        L;
        PreWhiten = false;
        noiseDecorr = 1; % No effect
        isCoilCompressed = false;
        coilCompMat  = 1; % Matrix for coil compression such as senseSVD 
        
        GIRF = [1,1,1];
        GIRFts = 2E-6;
    end
    
    methods(Abstract)
        
       dataRead(obj,shotNum,parNum,sliceNum,avgNum,phaseNum,echoNum,repNum,segNum)
       navRead(obj,shotNum,parNum,sliceNum,avgNum,phaseNum,echoNum,repNum,segNum)
       
    end
    
    methods
        
        function rInfo = recoInfo(filename, varargin)
            [filepath,name,ext] = fileparts(filename);
            switch(ext)
                case 'dat'
                    rInfo = recoInfoTwix(filename, varargin);
                case 'h5'
                    rInfo = recoInfoIsmrmrd(filename, varargin);
                otherwise
                    error('Unrecogized filetype. recoInfo currently supports Siemens Twix (.dat) and ISMRMRD HDF5 (.h5) files.')
            end
        end
        
        function voxelSize = voxelSize_mm(obj)
            
            voxelSize = zeros(1,3);
            voxelSize(1) = obj.FOV/(obj.N); % X
            voxelSize(2) = obj.FOV/(obj.N); % Y
            if obj.nPartitions > 1
                voxelSize(3) = obj.FOVSlice/(obj.nPartitions); % Z
            else
                voxelSize(3) = obj.sliceThickness;
            end
        end

    end
    
end


