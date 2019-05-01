classdef recoInfoSim < recoInfoAPI
    %recoInfoSim - Class for reconstruction info object for simulated data
    %
    % Syntax:  [rInfo] = recoInfo(kspace,data, timingVec, N,Nz, varargin)
    %
    % Inputs:
    %    kspace   - k-space trajectory in unitless dimensions [-N/2, N/2]
    %               in an array of dimensions:
    %               [nReadoutPts, nShots, nPartitions, nSlices, nEchoes,
    %               nAverages, nPhases, nRepetitions];
    %   data      - complex data in an array of dimensions matching kspace
    %               with the exception of coils, i.e. : [nReadoutPts, 
    %               nShots, nPartitions, nSlices, nEchoes, nAverages, 
    %               nPhases, nCoils];
    %   timingVec - vector of time indicies in seconds corresponding to
    %               readout. timingVec is relative to a zero at either the
    %               start of a FID (gradient echo) or spin echo. Times
    %               before a spin echo will have a negative index. In
    %               gradient echo, all times will be non-negative.
    %   N         - In-plane Image size
    %   Nz        - Through plane image size (in 3D cases, Nz = 1 for 2D).
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
    % Alex Cerjanic
    % University of Illinois at Urbana-Champaign
    % email address:
    % Website:
    % April 2019; Last revision: 29-Apr-2019
    
    properties
        DataMatrix;
    end
    
    methods
        function obj = recoInfoSim(kRead,kPhase,kSlice,data,timingVec,N,varargin)
            
            if size(kRead) ~= size(kPhase)
                error('Size of kRead must match size of kPhase!');
            end
            
            if size(kRead) ~= size(kSlice)
                error('Size of kRead must match size of kSlice!');
            end
            
            if size(kRead) ~= size(timingVec)
                error('Size of input k-space trajectory and input timing vectors must match!');
            end
            
            UseGIRF = false;
            PreWhiten = false;
            TE = 0;
            p = inputParser;
            addOptional(p,'UseGIRF',UseGIRF);
            addOptional(p,'PreWhiten',PreWhiten);
            addOptional(p,'TE', TE);
            parse(p,varargin{:});
            
            % Deal with parsed inputs
            obj.useGIRF = p.Results.UseGIRF;
            obj.PreWhiten = p.Results.PreWhiten;
            obj.TE = p.Results.TE;
            obj.gradTs = 10e-06; %Fixed on the siemens platform
            
            obj.N = N;
            
            sizeDims = size(data);
            
            obj.shotLength = length(kRead(:,1,1,1,1,1,1,1,1,1));
            
            obj.adcTs = timingVec(2) - timingVec(1);
            
            if length(sizeDims) == 2
                obj.nShots = sizeDims(2);
                obj.nPartitions = 1;
                obj.nSlices = 1;
                obj.nAverages = 1;
                obj.nPhases = 1;
                obj.nEchoes = 1;
                obj.nRepetitions = 1;
                obj.nCoils = 1;
                obj.nCoilsActual = 1;
            elseif length(sizeDims) == 3
                obj.nShots = sizeDims(2);
                obj.nPartitions = sizeDims(3);
                obj.nSlices = 1;
                obj.nAverages = 1;
                obj.nPhases = 1;
                obj.nEchoes = 1;
                obj.nRepetitions = 1;
                obj.nCoils = 1;
                obj.nCoilsActual = 1;
            elseif length(sizeDims) == 4
                obj.nShots = sizeDims(2);
                obj.nPartitions = sizeDims(3);
                obj.nAverages = sizeDims(4);
                obj.nPhases = 1;
                obj.nEchoes = 1;
                obj.nRepetitions = 1;
                obj.nCoils = 1;
                obj.nCoilsActual = 1;
            elseif length(sizeDims) == 5
                obj.nShots = sizeDims(2);
                obj.nPartitions = sizeDims(3);
                obj.nSlices = sizeDims(4);
                obj.nAverages = sizeDims(5);
                obj.nPhases = 1;
                obj.nEchoes = 1;
                obj.nRepetitions = 1;
                obj.nCoils = 1;
                obj.nCoilsActual = 1;
            elseif length(sizeDims) == 6
                obj.nShots = sizeDims(2);
                obj.nPartitions = sizeDims(3);
                obj.nSlices = sizeDims(4);
                obj.nAverages = sizeDims(5);
                obj.nPhases = sizeDims(6);
                obj.nEchoes = 1;
                obj.nRepetitions = 1;
                obj.nCoils = 1;
                obj.nCoilsActual = 1;
            elseif length(sizeDims) == 7
                obj.nShots = sizeDims(2);
                obj.nPartitions = sizeDims(3);
                obj.nSlices = sizeDims(4);
                obj.nAverages = sizeDims(5);
                obj.nPhases = sizeDims(6);
                obj.nEchoes = sizeDims(7);
                obj.nRepetitions = 1;
                obj.nCoils = 1;
                obj.nCoilsActual = 1;
            elseif length(sizeDims) == 8
                obj.nShots = sizeDims(2);
                obj.nSlices = sizeDims(3);
                obj.nPartitions = sizeDims(4);
                obj.nAverages = sizeDims(5);
                obj.nPhases = sizeDims(6);
                obj.nEchoes = sizeDims(7);
                obj.nRepetitions = sizeDims(8);
                obj.nCoils = 1;
                obj.nCoilsActual = 1;
            elseif length(sizeDims) == 9
                obj.nShots = sizeDims(2);
                obj.nSlices = sizeDims(3);
                obj.nPartitions = sizeDims(4);
                obj.nAverages = sizeDims(5);
                obj.nPhases = sizeDims(6);
                obj.nEchoes = sizeDims(7);
                obj.nRepetitions = sizeDims(8);
                obj.nCoils = sizeDims(9);
                obj.nCoilsActual = sizeDims(9);
            else
                error('Too many dimensions to input data matrix!');
            end

            if (obj.nEchoes ~= length(obj.TE))
                error('Number of echoes in data should match number of TEs');
            end
            
         
            % Set readShift and phaseShit, sliceShift, and rotMatrix to
            % nominal values (since this is simulation).
            obj.readShift = zeros(obj.nSlices,1);
            obj.phaseShift = zeros(obj.nSlices,1);
            obj.sliceShift = zeros(obj.nSlices,1);
            obj.rotMatrix = eye(3);
            obj.ExcOrder = 1:obj.nSlices;
            
            
            % Calculate the default halfVoxelShift for off center phasing
           
            obj.halfVoxelShiftSlice = 0;

            obj.halfVoxelShiftRead = 1./(2*obj.N);
            obj.halfVoxelShiftPhase = 1./(2*obj.N);
            
        
            
            % Call switchParser last so parameters can be overridden in the
            % sequence specific parsers

            obj.L = ceil(obj.shotLength*obj.adcTs/2E-3);
            
            
            obj.dataMask = true(obj.shotLength,1);
            
            %Loop through and set the data, timingVec, and kspace arrays.
            obj.DataMatrix = zeros(obj.shotLength, obj.nShots, ...
                                   obj.nPartitions, obj.nSlices, ...
                                   obj.nAverages, obj.nPhases, ...
                                   obj.nEchoes, obj.nRepetitions, ...
                                   obj.nCoils);
            obj.timingVec = zeros(obj.shotLength, obj.nShots, ...
                                   obj.nPartitions, obj.nSlices, ...
                                   obj.nAverages, obj.nPhases, ...
                                   obj.nEchoes, obj.nRepetitions);
            obj.kRead = zeros(obj.shotLength, obj.nShots, ...
                                   obj.nPartitions, obj.nSlices, ...
                                   obj.nAverages, obj.nPhases, ...
                                   obj.nEchoes, obj.nRepetitions);
            obj.kPhase = zeros(obj.shotLength, obj.nShots, ...
                                   obj.nPartitions, obj.nSlices, ...
                                   obj.nAverages, obj.nPhases, ...
                                   obj.nEchoes, obj.nRepetitions);
            obj.kSlice = zeros(obj.shotLength, obj.nShots, ...
                                   obj.nPartitions, obj.nSlices, ...
                                   obj.nAverages, obj.nPhases, ...
                                   obj.nEchoes, obj.nRepetitions);
                               
            for ii = 1:obj.nShots
                for jj = 1:obj.nPartitions
                    for kk = 1:obj.nSlices
                        for ll = 1:obj.nAverages
                            for mm = 1:obj.nPhases
                                for nn = 1:obj.nEchoes
                                    for oo = 1:obj.nRepetitions
                                        obj.DataMatrix(:,ii,jj,kk,ll,mm,nn,oo,:) = data(:,ii,jj,kk,ll,mm,nn,oo,:);
                                        obj.timingVec(:,ii,jj,kk,ll,mm,nn,oo) = timingVec(:,ii,jj,kk,ll,mm,nn,oo);
                                        obj.kRead(:,ii,jj,kk,ll,mm,nn,oo) = kRead(:,ii,jj,kk,ll,mm,nn,oo);
                                        obj.kPhase(:,ii,jj,kk,ll,mm,nn,oo) = kPhase(:,ii,jj,kk,ll,mm,nn,oo);
                                        obj.kSlice(:,ii,jj,kk,ll,mm,nn,oo) = kSlice(:,ii,jj,kk,ll,mm,nn,oo);
                                    end
                                end
                            end
                        end
                    end
                end
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
            
            %% Read Data via mapVBVD
            data = double(obj.DataMatrix(:,shotNum,parNum,sliceNum,avgNum,phaseNum,echoNum,repNum,:));
            data = reshape(data,[],length(shotNum),length(parNum),length(sliceNum),length(avgNum),length(phaseNum),length(echoNum),length(repNum),length(segNum),obj.nCoilsActual);
            
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
            

        end
        
        function [data] = navRead(obj,shotNum,parNum,sliceNum,avgNum,phaseNum,echoNum,repNum,segNum)
            error('Not currently implemented for recoInfoSim!');
        end
    end
    
end


