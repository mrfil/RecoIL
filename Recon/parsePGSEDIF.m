function obj = parsePGSEDIF(obj)

obj.nShotsDesigned = obj.DataObject.hdr.MeasYaps.sWipMemBlock.alFree{1};
obj.nShotsUsed = obj.DataObject.hdr.MeasYaps.sWipMemBlock.alFree{2};
%% Deal with navigator trajectory here
obj.adcTs = obj.DataObject.hdr.Dicom.DICOM.alICEProgramPara(9)*1E-6; %This could and should be added to the alICEProgramPara in the SpiralSBB class
obj.ptsToDrop = obj.DataObject.hdr.Dicom.DICOM.alICEProgramPara(8); 
obj.nEchoes = 1; % Ugly hack
navShotsUsed = 1; % Ugly hack
load kspace.mat
%load grads.mat
grads_mTm = readExternalGradFile('test.xml');
%The following code assumes G/cm not mT/m -> convert mT/m kspace.mat file
%to G/cm

grads = grads_mTm/10;

girfTs = 2E-6;
Greadi=interp1(0:obj.gradTs:obj.gradTs*length(grads(:,1))-obj.gradTs,grads(:,1),0:girfTs:obj.gradTs*length(grads(:,1))-girfTs,'previous','extrap')';
Gphasei=interp1(0:obj.gradTs:obj.gradTs*length(grads(:,2))-obj.gradTs,grads(:,2),0:girfTs:obj.gradTs*length(grads(:,2))-girfTs,'previous','extrap')';
Gslicei=interp1(0:obj.gradTs:obj.gradTs*length(grads(:,3))-obj.gradTs,grads(:,3),0:girfTs:obj.gradTs*length(grads(:,3))-girfTs,'previous','extrap')';
Greadi = Greadi(2:end); % Throwaway first point
Gphasei = Gphasei(2:end); % Throwaway first point
Gslicei = Gslicei(2:end); % 2D Sequence
GxN = interp1(0:obj.gradTs:obj.gradTs*length(grads(:,1))-obj.gradTs,grads(:,1),0:obj.adcTs:obj.gradTs*length(grads(:,1))-obj.gradTs)';
GyN = interp1(0:obj.gradTs:obj.gradTs*length(grads(:,2))-obj.gradTs,grads(:,2),0:obj.adcTs:obj.gradTs*length(grads(:,2))-obj.gradTs).';
GzN = interp1(0:obj.gradTs:obj.gradTs*length(grads(:,3))-obj.gradTs,grads(:,3),0:obj.adcTs:obj.gradTs*length(grads(:,3))-obj.gradTs).';
%Gxt=interp1(0:obj.gradTs:obj.gradTs*length(kxi)-obj.gradTs,kxi,0:obj.adcTs:obj.gradTs*length(kxi)-obj.gradTs)';

%sprintf('Performing %d rotations',obj.nShotsDesigned)
phi = 2*pi/obj.nShotsUsed;

load GIRFx.mat
load GIRFy.mat
load GIRFz.mat

%GIRFx = circshift(GIRFx,1);
%GIRFy = circshift(GIRFy,1);
%GIRFx = 5*decimate(GIRFx,5,'FIR');
%GIRFy = 5*decimate(GIRFy,5,'FIR');
%GIRFx = circshift(GIRFx,-2);
%GIRFy = circshift(GIRFy,-2);

for ii = 0:(obj.nShotsUsed-1)
    ang_rot = phi*(ii-(navShotsUsed-1)*floor(ii/navShotsUsed));
    Gphase(:,ii+1) = Gphasei*cos(ang_rot) + Greadi*sin(ang_rot);
    Gread(:,ii+1) = Greadi*cos(ang_rot) - Gphasei*sin(ang_rot);
    Gslice(:,ii+1) = Gslicei;
end

for ii = 0:(navShotsUsed-1)
    for jj = 1:length(Greadi)
        temp = obj.rotMatrix*[Gphase(jj,ii+1);Gread(jj,ii+1);Gslice(jj,ii+1)];
        Gx(jj,ii+1) = temp(1);
        Gy(jj,ii+1) = temp(2);
        Gz(jj,ii+1) = temp(3);
    end
end


for ii = 0:(navShotsUsed-1)
    ang_rot = phi*(ii-(obj.nShotsUsed-1)*floor(ii/obj.nShotsUsed));
    GxNominal(:,ii+1) = GxN*cos(ang_rot) + GyN*sin(ang_rot);
    GyNominal(:,ii+1) = GyN*cos(ang_rot) - GxN*sin(ang_rot);
    GzNominal(:,ii+1) = GzN;
end

%Convert Gx,Gy to physical gradients
invRotMatrix = inv(obj.rotMatrix);
for ii = 0:(navShotsUsed-1)
    for jj = 1:length(Greadi)
        temp = invRotMatrix*[Gphase(jj,ii+1);Gread(jj,ii+1);Gslice(jj,ii+1)];
        Gx(jj,ii+1) = temp(1);
        Gy(jj,ii+1) = temp(2);
        Gz(jj,ii+1) = temp(3);
    end
end

%Correct with GIRF
for ii = 0:(navShotsUsed-1)
    if obj.useGIRF
        Gxc(:,ii+1) = conv(Gx(:,ii+1),GIRFx,'same');
        Gyc(:,ii+1) = conv(Gy(:,ii+1),GIRFy,'same');
        Gzc(:,ii+1) = conv(Gz(:,ii+1),GIRFz,'same');
    else
        Gxc(:,ii+1) = Gx(:,ii+1);
        Gyc(:,ii+1) = Gy(:,ii+1);
        Gzc(:,ii+1) = Gz(:,ii+1);
    end
end
obj.rotMatrix = eye(3);

%Account for apparent negative 1 applied to slice direction
%obj.rotMatrix(3) = -1*obj.rotMatrix(3);

% Return to logical indices
for ii = 0:(navShotsUsed-1)
    for jj = 1:length(Greadi)
        
        temp = obj.rotMatrix*[Gxc(jj,ii+1);Gyc(jj,ii+1);Gzc(jj,ii+1)];
        Gphasec(jj,ii+1) = temp(1);
        Greadc(jj,ii+1) = temp(2);
        Gslicec(jj,ii+1) = temp(3);
    end
end

% if obj.nPartitions > 1
%     kzc = -obj.nPartitions/2:(obj.nPartitions/2-1);
% else
%     kzc = 0;
% end
gamma = 2*pi*4.257e3;
gambar = gamma/(2*pi);

for ii = 1:(navShotsUsed)
    kxt(:,ii) = cumsum([0; Greadc(:,ii)])*girfTs*obj.FOV*gambar;
    kyt(:,ii) = cumsum([0; Gphasec(:,ii)])*girfTs*obj.FOV*gambar;
    kzt(:,ii) = -1*cumsum([0; Gslicec(:,ii)])*girfTs*obj.FOV/2*gambar;
end

for ii = 1:(navShotsUsed)
    %obj.kxNominal(:,ii) = interp1(0:girfTs:girfTs*length(kxt(:,ii))-girfTs,kxt(:,ii),0:obj.adcTs:girfTs*length(kxt(:,ii)))';
    %obj.kyNominal(:,ii) = interp1(0:girfTs:girfTs*length(kyt(:,ii))-girfTs,kyt(:,ii),0:obj.adcTs:girfTs*length(kyt(:,ii)))';
    %obj.kzNominal(:,ii) = interp1(0:girfTs:girfTs*length(kzt(:,ii))-girfTs,kzt(:,ii),0:obj.adcTs:girfTs*length(kzt(:,ii)))';
    obj.kxNominalNav(:,ii) = interp1(0:obj.gradTs:obj.gradTs*length(kspace(:,1))-obj.gradTs,kspace(:,1)*obj.FOV,0:obj.adcTs:obj.gradTs*length(kspace(:,1))-obj.adcTs,'linear','extrap')';
    obj.kyNominalNav(:,ii) = interp1(0:obj.gradTs:obj.gradTs*length(kspace(:,2))-obj.gradTs,kspace(:,2)*obj.FOV,0:obj.adcTs:obj.gradTs*length(kspace(:,2))-obj.adcTs,'linear','extrap')';
    obj.kzNominalNav(:,ii) = interp1(0:obj.gradTs:obj.gradTs*length(kspace(:,3))-obj.gradTs,kspace(:,3)*obj.FOV/2,0:obj.adcTs:obj.gradTs*length(kspace(:,3))-obj.adcTs,'linear','extrap')';
end


for ii = 1:(navShotsUsed)
    kxc(:,ii) = interp1(0:girfTs:girfTs*length(kxt(:,ii))-girfTs,kxt(:,ii),0:obj.adcTs:girfTs*length(kxt(:,ii))-obj.adcTs)';
    kyc(:,ii) = interp1(0:girfTs:girfTs*length(kyt(:,ii))-girfTs,kyt(:,ii),0:obj.adcTs:girfTs*length(kyt(:,ii))-obj.adcTs)';
    kzc(:,ii) = interp1(0:girfTs:girfTs*length(kzt(:,ii))-girfTs,kzt(:,ii),0:obj.adcTs:girfTs*length(kzt(:,ii))-obj.adcTs)';
    
end

obj.ShotLengthNav = length(kxc(:,1)); %Assume kxt and kyt are the same length

obj.kxNav = repmat(kxc,[1,1,obj.nPartitions, obj.nSlices,obj.nAverages, obj.nPhases, obj.nEchoes, obj.nRepetitions,1]);
obj.kyNav = repmat(kyc,[1,1,obj.nPartitions, obj.nSlices,obj.nAverages, obj.nPhases, obj.nEchoes, obj.nRepetitions,1]);
%obj.kz = col(repmat(kzc,length(kxc(:,1))*obj.nShots,1));
obj.kzNav = repmat(kzc,[1,1,obj.nPartitions, obj.nSlices,obj.nAverages, obj.nPhases, obj.nEchoes, obj.nRepetitions,1]);

%Find center of kspace
[~,loc] = min(sum(sqrt(kspace(10:end,:).^2),2));

timingOffset = loc*obj.gradTs;

obj.timingVecNav = 0:obj.adcTs:obj.adcTs*(length(obj.kxNominalNav(:,1))-1);

obj.timingVecNav = obj.timingVecNav - timingOffset;

%% Deal with imaging readout here.

obj = parseSpiralOutReadout(obj);

end