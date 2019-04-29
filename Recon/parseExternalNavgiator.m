function obj = parseExternalNavgiator(obj)
%parseExternalNavigator Parses external navigator data into k-space
%trajectory for use with recoInfo.

%obj.nPhases = obj.DataObject.hdr.MeasYaps.sWipMemBlock.alFree{6};
EchoGap = obj.DataObject.hdr.Dicom.DICOM.alICEProgramPara(13); % Echo Gap between navigator and the echo time in usec.

%% Deal with navigator trajectory here
obj.adcTs = obj.DataObject.hdr.Dicom.DICOM.alICEProgramPara(9)*1E-6; %This could and should be added to the alICEProgramPara in the SpiralSBB class
obj.ptsToDrop = obj.DataObject.hdr.Dicom.DICOM.alICEProgramPara(8);
obj.nEchoes = 1; % Ugly hack
navShotsUsed = 1; % Ugly hack
load kspace.mat
%load grads.mat
grads_mTm = readExternalGradFile('CIVIC40x40x20.xml');
% Hard coding for this version of the sequence
obj.NNav = 40;
obj.nPartitionsNav = 20;
%The following code assumes G/cm not mT/m -> convert mT/m kspace.mat file
%to G/cm

grads = grads_mTm/10;

girfTs = 2E-6;
Greadi=interp1(0:obj.gradTs:obj.gradTs*length(grads(:,2))-obj.gradTs,grads(:,2),0:girfTs:obj.gradTs*length(grads(:,2))-girfTs,'previous','extrap')';
Gphasei=interp1(0:obj.gradTs:obj.gradTs*length(grads(:,1))-obj.gradTs,grads(:,1),0:girfTs:obj.gradTs*length(grads(:,1))-girfTs,'previous','extrap')';
Gslicei=interp1(0:obj.gradTs:obj.gradTs*length(grads(:,3))-obj.gradTs,grads(:,3),0:girfTs:obj.gradTs*length(grads(:,3))-girfTs,'previous','extrap')';
Greadi = Greadi(2:end); % Throwaway first point
Gphasei = Gphasei(2:end); % Throwaway first point
Gslicei = Gslicei(2:end); % 2D Sequence
GxN = interp1(0:obj.gradTs:obj.gradTs*length(grads(:,2))-obj.gradTs,grads(:,2),0:obj.adcTs:obj.gradTs*length(grads(:,2))-obj.gradTs)';
GyN = interp1(0:obj.gradTs:obj.gradTs*length(grads(:,1))-obj.gradTs,grads(:,1),0:obj.adcTs:obj.gradTs*length(grads(:,1))-obj.gradTs).';
GzN = interp1(0:obj.gradTs:obj.gradTs*length(grads(:,3))-obj.gradTs,grads(:,3),0:obj.adcTs:obj.gradTs*length(grads(:,3))-obj.gradTs).';
%Gxt=interp1(0:obj.gradTs:obj.gradTs*length(kxi)-obj.gradTs,kxi,0:obj.adcTs:obj.gradTs*length(kxi)-obj.gradTs)';

%sprintf('Performing %d rotations',obj.nShotsDesigned)
%phi = 2*pi/obj.nShotsUsed;

load GIRFx.mat
load GIRFy.mat
load GIRFz.mat

%GIRFx = circshift(GIRFx,1);
%GIRFy = circshift(GIRFy,1);
%GIRFx = 5*decimate(GIRFx,5,'FIR');
%GIRFy = 5*decimate(GIRFy,5,'FIR');
%GIRFx = circshift(GIRFx,-2);
%GIRFy = circshift(GIRFy,-2);



Gphase(:) = Gphasei;
Gread(:) = Greadi;
Gslice(:) = Gslicei;



for jj = 1:length(Greadi)
    temp = obj.rotMatrix*[Gphase(jj);Gread(jj);Gslice(jj)];
    Gx(jj) = temp(1);
    Gy(jj) = temp(2);
    Gz(jj) = temp(3);
end




GreadNominal(:) = GxN;
GphaseNominal(:) = GyN;
GsliceNominal(:) = GzN;


%Convert Gx,Gy to physical gradients
invRotMatrix = inv(obj.rotMatrix);

for jj = 1:length(Greadi)
    temp = invRotMatrix*[Gphase(jj);Gread(jj);Gslice(jj)];
    Gx(jj) = temp(1);
    Gy(jj) = temp(2);
    Gz(jj) = temp(3);
end


%Correct with GIRF

if obj.useGIRF
    Gxc(:) = conv(Gx(:),GIRFx,'same');
    Gyc(:) = conv(Gy(:),GIRFy,'same');
    Gzc(:) = conv(Gz(:),GIRFz,'same');
else
    Gxc(:) = Gx(:);
    Gyc(:) = Gy(:);
    Gzc(:) = Gz(:);
end

%obj.rotMatrix = eye(3);

%Account for apparent negative 1 applied to slice direction
%obj.rotMatrix(3) = -1*obj.rotMatrix(3);

% Return to logical indices
for jj = 1:length(Greadi)
    temp = obj.rotMatrix*[Gxc(jj);Gyc(jj);Gzc(jj)];
    Gphasec(jj) = temp(1);
    Greadc(jj) = temp(2);
    Gslicec(jj) = temp(3);
end

% if obj.nPartitions > 1
%     kzc = -obj.nPartitions/2:(obj.nPartitions/2-1);
% else
%     kzc = 0;
% end
gamma = 2*pi*4.257e3;
gambar = gamma/(2*pi);


kPhaset(:) = cumsum([0; Gphasec(:)])*girfTs*obj.FOV*gambar;
kReadt(:) = cumsum([0; Greadc(:)])*girfTs*obj.FOV*gambar;
kSlicet(:) = cumsum([0; Gslicec(:)])*girfTs*obj.FOV/2*gambar;

obj.kReadNominalNav = interp1(0:obj.gradTs:obj.gradTs*length(kspace(:,1))-obj.gradTs,kspace(:,1)*obj.FOV,0:obj.adcTs:obj.gradTs*length(kspace(:,1))-obj.adcTs,'linear','extrap')';
obj.kPhaseNominalNav = interp1(0:obj.gradTs:obj.gradTs*length(kspace(:,2))-obj.gradTs,kspace(:,2)*obj.FOV,0:obj.adcTs:obj.gradTs*length(kspace(:,2))-obj.adcTs,'linear','extrap')';
obj.kSliceNominalNav = interp1(0:obj.gradTs:obj.gradTs*length(kspace(:,3))-obj.gradTs,kspace(:,3)*obj.FOV/2,0:obj.adcTs:obj.gradTs*length(kspace(:,3))-obj.adcTs,'linear','extrap')';


kReadc = col(interp1(0:girfTs:girfTs*length(kReadt)-girfTs,kPhaset,0:obj.adcTs:girfTs*length(kReadt)-obj.adcTs));
kPhasec = col(interp1(0:girfTs:girfTs*length(kPhaset)-girfTs,kReadt,0:obj.adcTs:girfTs*length(kPhaset)-obj.adcTs));
kSlicec = col(interp1(0:girfTs:girfTs*length(kSlicet)-girfTs,kSlicet,0:obj.adcTs:girfTs*length(kSlicet)-obj.adcTs));


% For navigators, assume single shot design regardless of whether or not,
% the

obj.shotLengthNav = length(kReadc(:)); %Assume kxt and kyt are the same length
obj.dataMaskNav = true(obj.shotLengthNav,1);
obj.dataMaskNav(1:200) = 0; % Masking out the inital center of kspace in favor of the later return to the center of kspace
obj.dataMaskNav = logical(obj.dataMaskNav);

obj.kReadNav = repmat(kReadc,[1,obj.nShots,obj.nPartitions, obj.nSlices,obj.nAverages, obj.nPhases, obj.nEchoes, obj.nRepetitions,1]);
obj.kPhaseNav = repmat(kPhasec,[1,obj.nShots,obj.nPartitions, obj.nSlices,obj.nAverages, obj.nPhases, obj.nEchoes, obj.nRepetitions,1]);
%obj.kz = col(repmat(kzc,length(kxc(:,1))*obj.nShots,1));
obj.kSliceNav = repmat(kSlicec,[1,obj.nShots,obj.nPartitions, obj.nSlices,obj.nAverages, obj.nPhases, obj.nEchoes, obj.nRepetitions,1]);

%Find center of kspace
%[~,loc] = min(sum(sqrt(kspace(10:end,:).^2),2));

timingOffset = EchoGap*1E-6;

obj.timingVecNav = 0:obj.adcTs:obj.adcTs*(length(obj.kReadNav(:,1))-1);

%Make spiral in

obj.timingVecNav = flipdim(obj.timingVecNav,2);

obj.timingVecNav = obj.timingVecNav + timingOffset;

% Deal with half voxel shifts in navigator here

obj.halfVoxelShiftNavRead = 1/(2*obj.NNav);
obj.halfVoxelShiftNavPhase = 1/(2*obj.NNav);
if obj.multibandFactor > 1
    obj.halfVoxelShiftNavSlice = 1/(2*obj.multibandFactor);
else
    obj.halfVoxelShiftNavSlice = 1/(2*obj.nPartitionsNav);
end
end

