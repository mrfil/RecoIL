classdef recoInfoIsmrmrd < recoInfoAPI
    %recoInfoIsmrmrd recoInfo Class for ISMRMRD objects
    %   Assumes that your data is prepared and ready for reconstruction
    % Syntax:  [rInfo] = recoInfoIsmrmrd(filename)
    %
    % Inputs:
    %    filename - Name (and path) of Ismrmrd file (*.h5)
    %               
    %
    % Outputs:
    %    rInfo - Initalized recoInfoIsmrmrd object
    %
    %
    % Example:
    %    rInfo  = recoInfoIsmrmrd(filename);
    %    images = gridCoilImages(rInfo);
    %    im(sum(abs(images(:,:,:,:,:,1)).^2,4)); % Sum of Squares image for
    %                                            % first echo time in field map
    %
    %
    % Other m-files required: ISMRMRD library 
    % Subfunctions: none
    % MAT-files required: none
    %
    % Author: Alex Cerjanic
    % University of Illinois at Urbana-Champaign
    % email address:
    % Website:
    % February 2019; Last revision: 2019-02-28

    methods
        function obj = recoInfoIsmrmrd(filename, varargin)
            %recoInfoIsmsmrd Construct an instance of recoInfoIsmrmrd
            %   Detailed explanation goes here
            
            if nargin < 1 || isempty(filename)
                fileList = dir('*.h5');
                if length(fileList) ~= 1
                    error('Must specify filename if more than one .h5 file exists in the current directory.');
                else
                    obj.filename = fileList.name;
                end
            else
                obj.filename = filename;
            end
            
            obj.DataObject = ismrmrd.Dataset(obj.filename);
            hdr = ismrmrd.xml.deserialize(obj.DataObject.readxml);
            
            % Assume square matrix sizes since we deal with spiral mainly.
            obj.N = hdr.encoding.encodedSpace.matrixSize.x;
            
            % Map ISMRMRD extents to the recoInfo format
            if(isfield(hdr.encoding.encodingLimits,'slice'))
            	obj.nSlices = hdr.encoding.encodingLimits.slice.maximum + 1;
            else
                obj.nSlices = 1;
            end
            
            if(isfield(hdr.encoding.encodingLimits,'average'))
                obj.nAverages = hdr.encoding.encodingLimits.average.maximum + 1;
            else
                obj.nAverages = 1;
            end
            
            if(isfield(hdr.encoding.encodingLimits,'repetition'))
                obj.nRepetitions = hdr.encoding.encodingLimits.repetition.maximum + 1;
            else
                obj.nRepetitions = 1;
            end
            
            if(isfield(hdr.encoding.encodingLimits,'contrast'))
                obj.nEchoes = hdr.encoding.encodingLimits.contrast.maximum + 1;
            else
                obj.nEchoes = 1;
            end
            
            if(isfield(hdr.encoding.encodingLimits,'phase'))
                obj.nPhases = hdr.encoding.encodingLimits.phase.maximum + 1;
            else
                obj.nPhases = 1;
            end
            
            if(isfield(hdr.encoding.encodingLimits,'kspace_encoding_step_1'))
                obj.nShots = hdr.encoding.encodingLimits.kspace_encoding_step_1.maximum + 1;
            else
                obj.nShots = 1;
            end
                        
            if(isfield(hdr.encoding.encodingLimits,'kspace_encoding_step_2'))
            	obj.nPartitions = hdr.encoding.encodingLimits.kspace_encoding_step_2.maximum + 1;
            else
                obj.nPartitions = 1;
            end
            
            if(isfield(hdr.encoding.encodingLimits,'slice'))
                obj.nSlices = hdr.encoding.encodingLimits.slice.maximum + 1;
            else
                obj.nSlices = 1;
            end
            
            if(isfield(hdr.encoding.encodedSpace,'multibandFactor'))
                obj.multibandFactor = hdr.encoding.encodedSpace.multibandFactor;
            else
                obj.multibandFactor = 1;
            end


            obj.acqTracking = -1*ones(obj.nShots,obj.nPartitions,obj.nSlices,obj.nRepetitions,obj.nAverages,obj.nEchoes,obj.nPhases);
            
            for ii = 1:obj.DataObject.getNumberOfAcquisitions()
                
                acq = obj.DataObject.readAcquisition(ii);
                
                shot = acq.head.idx.kspace_encode_step_1 + 1;
                par = acq.head.idx.kspace_encode_step_2 + 1;
                slc = acq.head.idx.slice + 1;
                rep = acq.head.idx.repetition + 1;
                avg = acq.head.idx.average + 1;
                eco = acq.head.idx.contrast + 1;
                phs = acq.head.idx.phase + 1;
                
                % Save the index of the acquisition to the acqTracking
                % matrix
                obj.acqTracking(shot,par,slc,rep,avg,eco,phs) = ii;
            end
            
            % Figure out shot length from the first shot
            acq = obj.DataObject.readAcquisition(1);
            obj.shotLength = size(acq.data{1},1);
            obj.dataMask = true(obj.shotLength,1);
            
            obj.nCoils = acq.head.active_channels;
            obj.nCoilsActual = acq.head.active_channels;
            
            % Conslidate the the k space trajectories.
            obj.kRead = zeros(obj.shotLength, obj.nShots, obj.nPartitions, obj.nSlices, obj.nRepetitions, obj.nAverages, obj.nEchoes, obj.nPhases);
            obj.kPhase = zeros(obj.shotLength, obj.nShots, obj.nPartitions, obj.nSlices, obj.nRepetitions, obj.nAverages, obj.nEchoes, obj.nPhases);
            obj.kSlice = zeros(obj.shotLength, obj.nShots, obj.nPartitions, obj.nSlices, obj.nRepetitions, obj.nAverages, obj.nEchoes, obj.nPhases);
            obj.timingVec = zeros(obj.shotLength, obj.nShots, obj.nPartitions, obj.nSlices, obj.nRepetitions, obj.nAverages, obj.nEchoes, obj.nPhases);
            
            for shot = 1:obj.nShots
                for par = 1:obj.nPartitions
                    for slc = 1:obj.nSlices
                        for rep = 1:obj.nRepetitions
                            for avg = 1:obj.nAverages
                                for eco = 1:obj.nEchoes
                                    for phs = 1:obj.nPhases
                                        acq = obj.DataObject.readAcquisition(obj.acqTracking(shot, par, slc, rep, avg, eco, phs));
                                        
                                        trajData = acq.traj{1};
                
                                        obj.kRead(:,shot,par,slc,rep,avg,eco,phs) = trajData(1,:);
                                        obj.kPhase(:,shot,par,slc,rep,avg,eco,phs) = trajData(2,:);
                                        obj.kSlice(:,shot,par,slc,rep,avg,eco,phs) = trajData(3,:);
                                        obj.timingVec(:,shot,par,slc,rep,avg,eco,phs) = trajData(4,:);
                                    end
                                end
                            end
                        end
                    end
                end
            end
            
            
                
        end
        
        function data = dataRead(obj,shotNum,parNum,sliceNum,avgNum,phaseNum,echoNum,repNum,segNum)

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
            
            %% Grab data from ISMRMRD file
            
            data = zeros(obj.shotLength,length(shotNum),length(parNum),length(sliceNum),length(avgNum),length(phaseNum),length(echoNum),length(repNum), obj.nCoilsActual);
            % We need to get a list of the shots.
         
            for shot = 1:length(shotNum)
                for par = 1:length(parNum)
                    for slc = 1:length(sliceNum)
                        for rep = 1:length(repNum)
                            for avg = 1:length(avgNum)
                                for eco = 1:length(echoNum)
                                    for phs = 1:length(phaseNum)
                                        acq = obj.DataObject.readAcquisition(obj.acqTracking(shotNum(shot), parNum(par), sliceNum(slc), repNum(rep), avgNum(avg), echoNum(eco), phaseNum(phs)));
                                        
                                        data(:,shot,par,slc,rep,avg,eco,phs,:) = acq.data{1};
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        
        function data = navRead(obj,shotNum,parNum,sliceNum,avgNum,phaseNum,echoNum,repNum,segNum)
            error 'Not implemented yet for recoInfoISMRMRD'
        end
    end
end

