function obj = parse_pgseMBb(obj)
   %parse_pgseMBb.m - Parsing object for 2D/3D multiband PGSE Diffusion acquisitions.
   %
   % Syntax:  [obj] = parse_pgseMBb(obj)
   %
   % Inputs:
   %    obj - recoInfo object to be initialized.
   %
   % Outputs:
   %    obj - Initalized recoInfo object
   %
   %
   % Example:
   %     % This function should be called by recoInfo(). 
   %     rInfo = recoInfo('filename.dat');
   %
   %
   % Other m-files required: recoInfo.m
   % Subfunctions: none
   % MAT-files required: none
   %
   % Author:
   % Curtis Johnson - University of Delaware
   % Alex Cerjanic - University of Illinois at Urbana-Champaign
   % email address:
   % Website:
   % August 2018; Last revision: 11-Aug-2018
   
   obj.adcTs = obj.DataObject.hdr.Dicom.DICOM.alICEProgramPara(9)*1E-6; %This could and should be added to the alICEProgramPara in the SpiralSBB class
   obj.ptsToDrop = obj.DataObject.hdr.Dicom.DICOM.alICEProgramPara(3);
   ptx = obj.ptsToDrop;
   
   obj.multibandFactor = obj.DataObject.hdr.MeasYaps.sWipMemBlock.alFree{3};
   
   obj.nShotsDesigned = obj.DataObject.hdr.MeasYaps.sWipMemBlock.alFree{1};
   obj.nShotsUsed = obj.DataObject.hdr.MeasYaps.sWipMemBlock.alFree{2};
   
   %% Deal with imaging readout here.
   
   obj = parseSpiralOutReadout(obj);
   % Since Curtis "manually" changed the slice order, we need to override
   % recoInfo's normal (and sane) way of setting up the slice order to fix
   % that.
   
   % if mod(obj.nSlices/obj.multibandFactor,2) %odd number of slices
   %     obj.ExcOrder = [2:2:obj.nSlices, 1:2:obj.nSlices];
   % else %even number of slices
   %     obj.ExcOrder = [1:2:obj.nSlices, 2:2:obj.nSlices];
   % end
   
   %obj.adcTsSliceShiftTmp = obj.sliceShift;
   
   % Ugly hack for isocenter case
   SliceShiftTmp = obj.sliceThickness*(0:(obj.nSlices*obj.multibandFactor-1)) - obj.sliceThickness*((obj.nSlices*obj.multibandFactor-1))/2;
   SliceShiftTmp = SliceShiftTmp + (mean(obj.sliceShift)*obj.sliceThickness); % AMC suggested this, seems right
   for xx = 1:obj.nSlices
      %read_shift(xx) = mean(read_shift_tmp(xx:obj.multibandFactor:end));
      %phase_shift(xx) = mean(phase_shift_tmp(xx:obj.multibandFactor:end));
      %obj.sliceShift(xx) = mean(SliceShiftTmp(xx:obj.nSlices:(obj.nSlices*obj.multibandFactor)))/(obj.nSlices*obj.multibandFactor);
      obj.sliceShift(xx) = mean(SliceShiftTmp(xx:obj.nSlices:(obj.nSlices*obj.multibandFactor)))/(obj.sliceThickness*obj.nSlices*obj.multibandFactor);
   end
   
   % Let's trim the sliceShift vector once we are done to eliminate non-existant
   % voxel positions.
   obj.sliceShift = obj.sliceShift(1:obj.nSlices);
   
   
   obj.halfVoxelShiftRead = 1/(2*obj.N);
   obj.halfVoxelShiftPhase = 1/(2*obj.N);
   if obj.multibandFactor > 1
      obj.halfVoxelShiftSlice = 1/(2*obj.multibandFactor);
   else
      obj.halfVoxelShiftSlice = 0;
   end
   %obj.sliceShift(:) = obj.sliceShift(1);
   
   % obj.adcTs = obj.DataObject.hdr.Dicom.DICOM.alICEProgramPara(9)*1E-6; %This could and should be added to the alICEProgramPara in the SpiralSBB class
   % obj.ptsToDrop = obj.DataObject.hdr.Dicom.DICOM.alICEProgramPara(8); %Empirical factor
   %
   % load kspace_readout.mat
   % obj.ptsToDrop = ptx+((RUP*10e-6)/obj.adcTs);
   % % need to add the RUP points to obj.ptsToDrop
   % kx_vec = kspace(:,:,1);
   % ky_vec = kspace(:,:,2);
   % kz_vec = kspace(:,:,3);
   %
   % obj.nShotsUsed = length(kx_vec(1,:));
   % obj.ShotLength = length(kx_vec(:,1));
   % obj.nShots = obj.nShotsUsed;
   %
   % inplane_shots = obj.nShots/obj.multibandFactor;
   % ww = weight_vor(col(kx_vec(:,1:inplane_shots)),col(ky_vec(:,1:inplane_shots)),inplane_shots,1);
   % ww(ww>2)=2;
   % ww_vec = reshape(repmat(ww,[obj.multibandFactor 1]),[],obj.nShots);
   %
   % % for ii = 1:(obj.nShotsUsed)
   % %     obj.kxNominal(:,ii) = interp1(0:obj.gradTs:obj.gradTs*length(kx_vec(:,ii))-obj.gradTs,kx_vec(:,ii)*obj.FOV,0:obj.adcTs:obj.gradTs*length(kx_vec(:,ii))-obj.adcTs,'linear','extrap')';
   % %     obj.kyNominal(:,ii) = interp1(0:obj.gradTs:obj.gradTs*length(ky_vec(:,ii))-obj.gradTs,ky_vec(:,ii)*obj.FOV,0:obj.adcTs:obj.gradTs*length(ky_vec(:,ii))-obj.adcTs,'linear','extrap')';
   % %     obj.kzNominal(:,ii) = interp1(0:obj.gradTs:obj.gradTs*length(kz_vec(:,ii))-obj.gradTs,kz_vec(:,ii)*obj.FOV/2,0:obj.adcTs:obj.gradTs*length(kz_vec(:,ii))-obj.adcTs,'linear','extrap')';
   % % end
   % obj.kx = repmat(kx_vec,[1,1,obj.nPartitions, obj.nSlices,obj.nAverages, obj.nPhases, obj.nEchoes, obj.nRepetitions,1]);
   % obj.ky = repmat(ky_vec,[1,1,obj.nPartitions, obj.nSlices,obj.nAverages, obj.nPhases, obj.nEchoes, obj.nRepetitions,1]);
   % obj.kz = repmat(kz_vec,[1,1,obj.nPartitions, obj.nSlices,obj.nAverages, obj.nPhases, obj.nEchoes, obj.nRepetitions,1]);
   % obj.ww = repmat(ww_vec,[1,1,obj.nPartitions, obj.nSlices,obj.nAverages, obj.nPhases, obj.nEchoes, obj.nRepetitions,1]);
   %
   % obj.timingVec = 0:obj.adcTs:obj.adcTs*(length(obj.kx(:,1))-1);
   
   %% Deal with navigator trajectory here
   obj.nEchoes = 1; % Ugly hack
   navShotsUsed = 1; % Ugly hack
   FOVz = obj.sliceThickness*obj.nSlices*obj.nSlices*obj.multibandFactor;
   obj.ptsToDropNav = obj.DataObject.hdr.Dicom.DICOM.alICEProgramPara(8);
   %load mre_kspace_nav.mat
   path = fileparts(mfilename('fullpath'));

   [grads_mTm, FOVRef, FOVzRef, RUP] = readExternalGradFile([path '/trajectories/Rxy4Rz1N40Nz4/Rxy4Rz1Nxy40Nz4.xml']);
   RUP
   kspace = calcKSpaceFromGrads(obj, grads_mTm, obj.adcTs, FOVRef, 'FOVz', FOVzRef, 'RUP', RUP, 'UseGIRF', obj.useGIRF );
   obj.ptsToDropNav = obj.DataObject.hdr.Dicom.DICOM.alICEProgramPara(8) + ((RUP*10e-6)/obj.adcTs);
   % need to add the RUP points to obj.ptsToDrop
   kRead_nav = kspace(:,1);
   kPhase_nav = kspace(:,2);
   kSlice_nav = kspace(:,3);
   
   obj.shotLengthNav = length(kRead_nav(:,1)); %Assume kxt and kyt are the same length
   
   obj.kReadNav = repmat(kRead_nav,[1,obj.nShots,obj.nPartitions, obj.nSlices,obj.nAverages, obj.nPhases, obj.nEchoes, obj.nRepetitions,1]);
   obj.kPhaseNav = repmat(kPhase_nav,[1,obj.nShots,obj.nPartitions, obj.nSlices,obj.nAverages, obj.nPhases, obj.nEchoes, obj.nRepetitions,1]);
   obj.kSliceNav = repmat(kSlice_nav,[1,obj.nShots,obj.nPartitions, obj.nSlices,obj.nAverages, obj.nPhases, obj.nEchoes, obj.nRepetitions,1]);
   
   %Hack for now.
   obj.NNav = 40;
   obj.nPartitionsNav = obj.multibandFactor;
   
   obj.timingVecNav = -1*col(flip(0:obj.adcTs:obj.adcTs*(length(obj.kReadNav(:,1))-1)));

   
   %Override default 3D half voxel shift in the z direction
   obj.halfVoxelShiftNavRead = -1/(2*obj.NNav);
   obj.halfVoxelShiftNavPhase = -1/(2*obj.NNav);
   if obj.multibandFactor > 1
      obj.halfVoxelShiftNavSlice = -1/(2*obj.multibandFactor); % CLJ: making this positive, since that is what works
   else
      obj.halfVoxelShiftNavSlice = 0;
   end
   
   % moving these into the parse out of reconMultibandMRE
   obj.dataMask = true(obj.shotLength,1);
   obj.dataMaskNav = true(obj.shotLengthNav,1);
end