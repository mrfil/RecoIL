classdef recoInfoTwix < recoInfoAPI
    %recoInfoTwix - Class for reconstruction info object for Siemens VB/VD/VE files
    %
    % Syntax:  [rInfo] = recoInfo(filename)
    %
    % Inputs:
    %    filename - Name (and path) of Siemens VB/VD/VE *.dat file that can be
    %               read with mapVBVD
    %
    % Outputs:
    %    rInfo - Initalized recoInfo object
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
    % April 2019; Last revision: 27-Apr-2019
    
    properties
        % Siemens raw .dat file related parameters
        SequenceString = '';
        
    end
    
    methods
        function obj = recoInfoTwix(filename,varargin)
            
            if nargin < 1 || isempty(filename)
                fileList = dir('*.dat');
                if length(fileList) ~= 1
                    error('Must specify filename if more than one .dat file exists in the current directory.');
                else
                    filename = fileList.name;
                end
            end
            
            UseGIRF = false;
            SequenceString = '';
            PreWhiten = false;
            p = inputParser;
            addOptional(p,'SequenceString',SequenceString);
            addOptional(p,'UseGIRF',UseGIRF);
            addOptional(p,'PreWhiten',PreWhiten);
            
            parse(p,varargin{:});
            
            % Deal with parsed inputs
            obj.useGIRF = p.Results.UseGIRF;
            obj.PreWhiten = p.Results.PreWhiten;
            
            obj.filename = filename;
            DataObjectTemp = mapVBVD(filename); %Open a mapVBVD object matching the file we want to reconstruct
           
            
            % They come as cell matrix most of the time but not always
            if iscell(DataObjectTemp)
                obj.DataObject = DataObjectTemp{1,end};
            else
                obj.DataObject = DataObjectTemp;
            end

            % Let's grab some interesting header information
            obj.scannerManufacturer = obj.DataObject.hdr.Dicom.DICOM.Manufacturer;
            obj.scannerModel = obj.DataObject.hdr.Dicom.DICOM.ManufacturersModelName;
            obj.scannerSoftwareVer = obj.DataObject.hdr.Dicom.DICOM.SoftwareVersions;
            
            % Grab the scanner serial for chosing the right GIRF
            obj.scannerSerial = obj.DataObject.hdr.Dicom.DICOM.DeviceSerialNumber;
            
            % Deal with the Noise Correlation Matrix if it exists
            if isfield(obj.DataObject,'noise')
                noiseData = squeeze(obj.DataObject.noise());
                noiseData = permute(noiseData,[2,1,3]);
                obj.noiseCorr = calcNoiseCorr(noiseData);
            end
            
            obj.gradTs = 10e-06; %Fixed on the siemens platform
            
            sizeDims = size(obj.DataObject.image());
            
            if length(sizeDims) == 3
                obj.nSlices = 1;
                obj.nPartitions = 1;
                obj.nEchoes = 1;
                obj.nAverages = 1;
                obj.nPhases = 1;
            elseif length(sizeDims) == 4
                obj.nSlices = 1;
                obj.nPartitions = sizeDims(4);
                obj.nEchoes = 1;
                obj.nAverages = 1;
                obj.nPhases = 1;
            elseif length(sizeDims) == 5
                obj.nSlices = sizeDims(5);
                obj.nPartitions = sizeDims(4);
                obj.nEchoes = 1;
                obj.nAverages = 1;
                obj.nPhases = 1;
            elseif length(sizeDims) == 6
                obj.nAverages = sizeDims(6);
                obj.nSlices = sizeDims(5);
                obj.nPartitions = sizeDims(4);
                obj.nEchoes = 1;
                obj.nPhases = 1;
            elseif length(sizeDims) == 7
                obj.nAverages = sizeDims(6);
                obj.nSlices = sizeDims(5);
                obj.nPartitions = sizeDims(4);
                obj.nEchoes = 1;
                obj.nPhases = sizeDims(7);
            else
                obj.nPhases = sizeDims(7);
                obj.nAverages = sizeDims(6); % Mapped from Averages (VB) / Acquisitions (VD/VE) in MDH (VB/VD/VE)
                obj.nSlices = sizeDims(5); % Mapped from Channels in MDH
                obj.nPartitions = sizeDims(4); % Mapped from Partitions in MDH
                obj.nEchoes = sizeDims(8); % Mapped from Echoes in MDH
            end
            obj.nCoilsActual = sizeDims(2); % Mapped from Channel in MDH
            obj.nCoils = sizeDims(2);
            obj.nShots = sizeDims(3); % Mapped from Lines in MDH
            if length(sizeDims) < 9
                obj.nRepetitions = 1;
            else
                obj.nRepetitions = sizeDims(9);
            end
            if length(sizeDims) < 10 % If nADCs = 1, there won't be a 10th entry
                obj.nADCs = 1;
            else
                obj.nADCs = sizeDims(10); % Mapped from Set in MDH
            end
            
            % Deal with sequence specific header
            
            
            obj.FOV = obj.DataObject.hdr.Phoenix.sSliceArray.asSlice{1}.dReadoutFOV/10;
            
            if obj.nPartitions > 1
                obj.FOVSlice = obj.DataObject.hdr.Phoenix.sSliceArray.asSlice{1}.dThickness;
                obj.sliceThickness = obj.FOVSlice/obj.nPartitions;
            else
                obj.sliceThickness = obj.DataObject.hdr.Phoenix.sSliceArray.asSlice{1}.dThickness;
            end
            
            obj.sliceThickness = obj.DataObject.hdr.Phoenix.sSliceArray.asSlice{1}.dThickness;
            obj.N = obj.DataObject.hdr.Dicom.DICOM.lBaseResolution;
            if ~strcmp(p.Results.SequenceString,'')
                obj.SequenceString = p.Results.SequenceString;
            else
                obj.SequenceString = obj.DataObject.hdr.Spice.PARC.PIPE.Receiver.source.SequenceString;
                if length(obj.SequenceString) < 7
                    obj.SequenceString = '';
                else
                    obj.SequenceString = obj.SequenceString(1:7);
                end
            end
            
            obj.TE = cell2mat(obj.DataObject.hdr.MeasYaps.alTE);
            obj.TR = cell2mat(obj.DataObject.hdr.MeasYaps.alTR);
            
            [obj.readShift,obj.phaseShift,obj.sliceShift, obj.rotMatrix] = imageShift(obj,1);

            %Deal with Siemens slice order
            if (obj.DataObject.hdr.MeasYaps.sSliceArray.ucMode == 4)%interleaved
                if mod(obj.nSlices,2) %odd number of slices
                    % WIP: AMC
                    %obj.ExcOrder = flip([1:2:obj.nSlices, 2:2:obj.nSlices]);
                    obj.ExcOrder = [1:2:obj.nSlices, 2:2:obj.nSlices];
                else %even number of slices
                    % WIP: AMC
                    %obj.ExcOrder = flip([2:2:obj.nSlices, 1:2:obj.nSlices]);
                    obj.ExcOrder = [2:2:obj.nSlices, 1:2:obj.nSlices];
                end
                
            elseif (obj.DataObject.hdr.MeasYaps.sSliceArray.ucMode == 1)%ascending
                obj.ExcOrder = 1:obj.nSlices;
                %obj.ExcOrder = obj.nSlices:-1:1;
            elseif (obj.DataObject.hdr.MeasYaps.sSliceArray.ucMode == 2)%descending
                obj.ExcOrder = obj.nSlices:-1:1;
                %obj.ExcOrder = 1:obj.nSlices;
            else
                obj.ExcOrder = 1:obj.nSlices;
            end
            
            
            % Calculate the default halfVoxelShift for off center phasing
            
            if obj.nPartitions > 1
                obj.halfVoxelShiftSlice = 1./(2*obj.nPartitions);
            else
                obj.halfVoxelShiftSlice = 0;
            end
            obj.halfVoxelShiftRead = 1./(2*obj.N);
            obj.halfVoxelShiftPhase = 1./(2*obj.N);
            
            
            if (obj.useGIRF == true)
                GIRFFileName = [ obj.scannerModel '_' obj.scannerSerial '.mat'];
                load(['GIRFx_' GIRFFileName], 'GIRFx');
                load(['GIRFy_' GIRFFileName], 'GIRFy');
                load(['GIRFz_' GIRFFileName], 'GIRFz');
                obj.GIRF = horzcat(col(GIRFx),col(GIRFy),col(GIRFz));
            end
            
            % Call switchParser last so parameters can be overridden in the
            % sequence specific parsers
            obj = switchParser(obj,obj.SequenceString);
            % These rely on rInfo.shotLength, so we need o
            obj.L = ceil(obj.shotLength*obj.adcTs/3E-3);
            
            
            obj.dataMask = true(obj.shotLength,1);
            if ~isempty(obj.shotLengthNav) && isempty(obj.dataMaskNav)
                obj.dataMaskNav = true(obj.shotLengthNav,1);
            end
        end
        
        function  obj = switchParser(obj, sequenceString)
            
            % Select correct sequence specific parser
            if (strcmp(sequenceString,'SENFMv1'))
                obj = parseSENFMv1(obj);
            elseif (strcmp(sequenceString,'SENFMv2'))
                obj = parseSENFMv2(obj);
            elseif (strcmp(sequenceString,'SENDIFF'))
                obj = parseSENDIFF(obj);
            elseif (strcmp(sequenceString,'PGSEDTI'))
                obj = parsePGSEDTI(obj);
            elseif (strcmp(sequenceString,'PGSEARB'))
                obj = parsePGSEARB(obj);
            elseif (strcmp(sequenceString,'ExtNav1'))
                obj = parseExtNavTest(obj);
            elseif (strcmp(sequenceString,'ExtNav2'))
                obj = parseExtNav2Test(obj);
            elseif (strcmp(sequenceString,'FlxTPI1'))
                obj = parseFlxTPI1(obj);
            elseif (strcmp(sequenceString,'MRF-Fv1'))
                obj = parseMRFFv1(obj);
            elseif (strcmp(sequenceString,'mre3Dmb'))
                obj = parse_mre3Dmb(obj);
            elseif (strcmp(sequenceString,'iVASOv1'))
                obj = parseiVASOv1(obj);
            elseif (strcmp(sequenceString,'STE_ARB'))
                obj = parseSTE_ARB(obj);
            elseif (strcmp(sequenceString,'SplInFl'))
                obj = parseSplInFl(obj);
            elseif (strcmp(sequenceString,'SpFLASH'))
                obj = parseSpFLASH(obj);
            elseif (strcmp(sequenceString,'MB_STEa'))
                obj = parseMB_STEa(obj);
            elseif (strcmp(sequenceString,'pgseMBa'))
                obj = parse_pgseMBa(obj);
            elseif (strcmp(sequenceString,'MB_STET'))
                obj = parseMB_STET(obj);
            elseif (strcmp(sequenceString,'pgseMBT'))
                obj = parse_pgseMBT(obj);
            elseif (strcmp(sequenceString,'MB_STEb'))
                obj = parseMB_STEb(obj);
            elseif (strcmp(sequenceString,'pgseMBb'))
                obj = parse_pgseMBb(obj);
            elseif (strcmp(sequenceString,'pgseMBc'))
                obj = parse_pgseMBc(obj);
            elseif (strcmp(sequenceString,'pgseMBd'))
                obj = parse_pgseMBd(obj);
            elseif (strcmp(sequenceString,'pgseMBe'))
                obj = parse_pgseMBe(obj);
            elseif (strcmp(sequenceString,'MB_STEe'))
                obj = parseMB_STEe(obj);
            else % fall through to error
                error('Parser not found! Has a parser been specified for the sequence string!?');
            end
        end
        
        function [data] = dataRead(obj,shotNum,parNum,sliceNum,avgNum,phaseNum,echoNum,repNum,segNum)
            
            %% Prepare input arguments
            if nargin < 9
                segNum = 1; % MRFIL doesn't use segment so force segment to be 1 - AMC 2/21/2017
            end
            if isempty(segNum) % segNum = []
                segNum = 1;
            end
            
            if nargin < 8
                repNum = 1:obj.nRepetitions;
            end
            if isempty(repNum)
                repNum = 1:obj.nRepetitions;
            end
            
            if nargin < 7
                echoNum = 1:obj.nEchoes;
            end
            if isempty(echoNum)
                echoNum = 1:obj.nEchoes;
            end
            
            if nargin < 6
                phaseNum = 1:obj.nPhases;
            end
            if isempty(phaseNum)
                phaseNum = 1:obj.nPhases;
            end
            
            if nargin < 5
                avgNum = 1:obj.nAverages;
            end
            if isempty(avgNum)
                avgNum = 1:obj.nAverages;
            end
            
            if nargin < 4
                sliceNum = 1:obj.nSlices;
            end
            if isempty(sliceNum)
                sliceNum = 1:obj.nSlices;
            end
            
            if nargin < 3
                parNum = 1:obj.nPartitions;
            end
            if isempty(parNum)
                parNum = 1:obj.nPartitions;
            end
            
            if nargin < 2
                shotNum = 1:obj.nShots;
            elseif isempty(shotNum)
                shotNum = 1:obj.nShots;
            end
            
            % Reverse look up for slices in case of interleaved or
            % descending slice order
            sliceToMDH(obj.ExcOrder) = 1:obj.nSlices;
            
            %% Read Data via mapVBVD
            data = double(obj.DataObject.image(:,:,shotNum,parNum,sliceToMDH(sliceNum),avgNum,phaseNum,echoNum,repNum,:,segNum,1,1,1,1,1));
            data = permute(data,[1,10,3,4,5,6:8,9,11,2,12,13,14,15,16]); %Assuming that Ida-Ide are not used;
            data = reshape(data,[],length(shotNum),length(parNum),length(sliceNum),length(avgNum),length(phaseNum),length(echoNum),length(repNum),length(segNum),obj.nCoilsActual);
            data = data(1+obj.ptsToDrop:obj.shotLength+obj.ptsToDrop,:,:,:,:,:,:,:,:,:);
            %data = reshape(data,[],obj.nPartitions,obj.nSlices,obj.nEchoes,obj.nRepetitions,obj.nCoils);
            
            if obj.PreWhiten == true
                for ii = 1:length(shotNum)
                    for jj = 1:length(parNum)
                        for kk = 1:length(sliceNum)
                            for ll = 1:length(avgNum)
                                for mm = 1:length(phaseNum)
                                    for nn = 1:length(echoNum)
                                        for oo = 1:length(repNum)
                                            for pp = 1:length(segNum)
                                                data(:,ii,jj,kk,ll,mm,nn,oo,pp,:) = (obj.noiseDecorr*squeeze(data(:,ii,jj,kk,ll,mm,nn,oo,pp,:)).').';
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
            
            % Handle SVD coil compression
            
            if obj.isCoilCompressed == true
                sizeData = size(data);
                dataComp = zeros([sizeData(1:end-1), obj.nCoils]);
                for ii = 1:length(shotNum)
                    for jj = 1:length(parNum)
                        for kk = 1:length(sliceNum)
                            for ll = 1:length(avgNum)
                                for mm = 1:length(phaseNum)
                                    for nn = 1:length(echoNum)
                                        for oo = 1:length(repNum)
                                            for pp = 1:length(segNum)
                                                dataComp(:,ii,jj,kk,ll,mm,nn,oo,pp,:) = squeeze(data(:,ii,jj,kk,ll,mm,nn,oo,pp,:))*obj.coilCompMat;
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
                data = dataComp;
            end
            %             halfVoxelShift = 0;
            %             % Add half voxel shift to the data if 3D encoded
            %             if obj.nPartitions > 1
            %                 halfVoxelShift = 1./(2*obj.nPartitions);
            %             end
            
            % Phase data for off center shifts
            for ii = 1:length(shotNum)
                for jj = 1:length(parNum)
                    for kk = 1:length(sliceNum)
                        for ll = 1:length(avgNum)
                            for mm = 1:length(phaseNum)
                                for nn = 1:length(echoNum)
                                    for oo = 1:length(repNum)
                                        for pp = 1:length(segNum)
                                            for qq = 1:obj.nCoils
                                                curData = squeeze(data(:,ii,jj,kk,ll,mm,nn,oo,pp,qq));
                                                phaseTerm = exp(-1j*( ...
                                                             (obj.phaseShift(obj.ExcOrder(sliceNum(kk))) + obj.halfVoxelShiftPhase) .* obj.kPhase(:,shotNum(ii),parNum(jj),sliceNum(kk),avgNum(ll),phaseNum(mm),echoNum(nn),repNum(oo),segNum(pp)) + ...
                                                             (obj.readShift(obj.ExcOrder(sliceNum(kk))) + obj.halfVoxelShiftRead) .* obj.kRead(:,shotNum(ii),parNum(jj),sliceNum(kk),avgNum(ll),phaseNum(mm),echoNum(nn),repNum(oo),segNum(pp)) + ...
                                                             (obj.sliceShift(obj.ExcOrder(sliceNum(kk))) + obj.halfVoxelShiftSlice) .* obj.kSlice(:,shotNum(ii),parNum(jj),sliceNum(kk),avgNum(ll),phaseNum(mm),echoNum(nn),repNum(oo),segNum(pp))) * 2*pi);
                                                data(:,ii,jj,kk,ll,mm,nn,oo,pp,qq) = curData .* phaseTerm;
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
            
        end
        
        function [data] = navRead(obj,shotNum,parNum,sliceNum,avgNum,phaseNum,echoNum,repNum,segNum)
            
            %% Prepare input arguments
            if nargin < 9
                segNum = 1; % MRFIL doesn't use segment so force segment to be 1 - AMC 2/21/2017
            end
            if isempty(segNum) % segNum = []
                segNum = 1;
            end
            
            if nargin < 8
                repNum = 1:obj.nRepetitions;
            end
            if isempty(repNum)
                repNum = 1:obj.nRepetitions;
            end
            
            if nargin < 7
                echoNum = 1:obj.nEchoes;
            end
            if isempty(echoNum)
                echoNum = 1:obj.nEchoes;
            end
            
            if nargin < 6
                phaseNum = 1:obj.nPhases;
            end
            if isempty(phaseNum)
                phaseNum = 1:obj.nPhases;
            end
            
            if nargin < 5
                avgNum = 1:obj.nAverages;
            end
            if isempty(avgNum)
                avgNum = 1:obj.nAverages;
            end
            
            if nargin < 4
                sliceNum = 1:obj.nSlices;
            end
            if isempty(sliceNum)
                sliceNum = 1:obj.nSlices;
            end
            
            if nargin < 3
                parNum = 1:obj.nPartitions;
            end
            if isempty(parNum)
                parNum = 1:obj.nPartitions;
            end
            
            if nargin < 2
                shotNum = 1:obj.nShots;
            elseif isempty(shotNum)
                shotNum = 1:obj.nShots;
            end
            
            % Reverse look up for slices in case of interleaved or
            % descending slice order
            sliceToMDH(obj.ExcOrder) = 1:obj.nSlices;
            
            %% Read Data via mapVBVD
            data = double(obj.DataObject.RTfeedback(:,:,shotNum,parNum,sliceToMDH(sliceNum),avgNum,phaseNum,echoNum,repNum,:,segNum,1,1,1,1,1));
            data = permute(data,[1,10,3,4,5,6:8,9,11,2,12,13,14,15,16]); %Assuming that Ida-Ide are not used;
            data = reshape(data,[],length(shotNum),length(parNum),length(sliceNum),length(avgNum),length(phaseNum),length(echoNum),length(repNum),length(segNum),obj.nCoilsActual);
            data = data(1+obj.ptsToDropNav:obj.shotLengthNav+obj.ptsToDropNav,:,:,:,:,:,:,:,:,:);
            %data = reshape(data,[],obj.nPartitions,obj.nSlices,obj.nEchoes,obj.nRepetitions,obj.nCoils);
            
            if obj.PreWhiten == true
                for ii = 1:length(shotNum)
                    for jj = 1:length(parNum)
                        for kk = 1:length(sliceNum)
                            for ll = 1:length(avgNum)
                                for mm = 1:length(phaseNum)
                                    for nn = 1:length(echoNum)
                                        for oo = 1:length(repNum)
                                            for pp = 1:length(segNum)
                                                data(:,ii,jj,kk,ll,mm,nn,oo,pp,:) = (obj.noiseDecorr*squeeze(data(:,ii,jj,kk,ll,mm,nn,oo,pp,:)).').';
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
            
            % Handle SVD coil compression
            if obj.isCoilCompressed == true
                dataOut = zeros(obj.shotLengthNav,length(shotNum),length(parNum),length(sliceNum),length(avgNum),length(phaseNum),length(echoNum),length(repNum),length(segNum),obj.nCoils);
                for ii = 1:length(shotNum)
                    for jj = 1:length(parNum)
                        for kk = 1:length(sliceNum)
                            for ll = 1:length(avgNum)
                                for mm = 1:length(phaseNum)
                                    for nn = 1:length(echoNum)
                                        for oo = 1:length(repNum)
                                            for pp = 1:length(segNum)
                                                dataOut(:,ii,jj,kk,ll,mm,nn,oo,pp,:) = (squeeze(data(:,ii,jj,kk,ll,mm,nn,oo,pp,:))*obj.coilCompMat);
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
                data = dataOut;
            end
            
            
            
            % Phase data for off center shifts
            for ii = 1:length(shotNum)
                for jj = 1:length(parNum)
                    for kk = 1:length(sliceNum)
                        for ll = 1:length(avgNum)
                            for mm = 1:length(phaseNum)
                                for nn = 1:length(echoNum)
                                    for oo = 1:length(repNum)
                                        for pp = 1:length(segNum)
                                            for qq = 1:obj.nCoils
                                                curData = squeeze(data(:,ii,jj,kk,ll,mm,nn,oo,pp,qq));
                                                phaseTerm = exp(-1j*( ...
                                                                (obj.phaseShift(obj.ExcOrder(sliceNum(kk))) + obj.halfVoxelShiftNavPhase) .* obj.kPhaseNav(:,shotNum(ii),parNum(jj),sliceNum(kk),avgNum(ll),phaseNum(mm),echoNum(nn),repNum(oo),segNum(pp)) + ...
                                                                (obj.readShift(obj.ExcOrder(sliceNum(kk))) + obj.halfVoxelShiftNavRead) .* obj.kReadNav(:,shotNum(ii),parNum(jj),sliceNum(kk),avgNum(ll),phaseNum(mm),echoNum(nn),repNum(oo),segNum(pp)) + ...
                                                                (obj.sliceShift(obj.ExcOrder(sliceNum(kk))) + obj.halfVoxelShiftNavSlice) .* obj.kSliceNav(:,shotNum(ii),parNum(jj),sliceNum(kk),avgNum(ll),phaseNum(mm),echoNum(nn),repNum(oo),segNum(pp)))*2*pi);
                                                data(:,ii,jj,kk,ll,mm,nn,oo,pp,qq) = curData .* phaseTerm;
                                            end 
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
            
        end
        
    end
    
end


