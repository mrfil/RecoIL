function obj = parseSpiralInReadoutGC(obj)

girfTs = 2E-6;

obj.ptsToDrop = obj.DataObject.hdr.Dicom.DICOM.alICEProgramPara(3); %Empirical factor

obj.gradAmp = obj.DataObject.hdr.Dicom.DICOM.alICEProgramPara(2)/10; % was always 2.4
obj.gradSlew = obj.DataObject.hdr.Dicom.DICOM.alICEProgramPara(1);  % was 120
obj.adcTs = obj.DataObject.hdr.Dicom.DICOM.alICEProgramPara(4)*1E-6; 
[Greado,Gphaseo,kreadi,kphasei] = genspi(obj.FOV,obj.N,obj.nShotsDesigned,obj.gradAmp,obj.gradSlew,obj.gradTs); %[Gx,Gy,kxi,kyi,sx,sy]
%[Gphaseo,Greado,kphasei,kreadi] = genspi(obj.FOV,obj.N,obj.nShotsDesigned,obj.gradAmp,obj.gradSlew,obj.gradTs); %[Gx,Gy,kxi,kyi,sx,sy]

%% reverse and invert all outputs of genspi
Greado = -1*Greado(end:-1:1);
Gphaseo = -1*Gphaseo(end:-1:1);
kreadi = -1*kreadi(end:-1:1);
kphasei = -1*kphasei(end:-1:1);
    

Greadi  = interp1(0:10:10*length(Greado)-1,Greado,0:2:10*length(kreadi)-1,'next','extrap')';
Gphasei = interp1(0:10:10*length(Gphaseo)-1,Gphaseo,0:2:10*length(kphasei)-1,'next','extrap')';

Greadi  = [zeros(4,1); Greadi ];
Gphasei = [zeros(4,1); Gphasei];

Greadi(isnan(Greadi)) = 0;
Gphasei(isnan(Gphasei)) = 0;
Gslicei = zeros(size(Greadi));

kreadN=interp1(0:obj.gradTs:obj.gradTs*length(kreadi)-obj.gradTs,kreadi,0:obj.adcTs:obj.gradTs*length(kreadi)-obj.adcTs,'linear','extrap')';
kphaseN=interp1(0:obj.gradTs:obj.gradTs*length(kphasei)-obj.gradTs,kphasei,0:obj.adcTs:obj.gradTs*length(kphasei)-obj.adcTs,'linear','extrap')';

%Gxt=interp1(0:obj.gradTs:obj.gradTs*length(kxi)-obj.gradTs,kxi,0:obj.adcTs:obj.gradTs*length(kxi)-obj.adcTs)';
%kxN = kxN(1:end-1);
%kyN = kyN(1:end-1);

%sprintf('Performing %d rotations',obj.nShotsDesigned)
phi = 2*pi/obj.nShotsUsed;

load GIRFx.mat
load GIRFy.mat
load GIRFz.mat

for ii = 0:(obj.nShotsUsed-1)
    ang_rot = phi*(ii-(obj.nShotsUsed-1)*floor(ii/obj.nShotsUsed));
    Gread(:,ii+1) = Greadi*cos(ang_rot) + Gphasei*sin(ang_rot);
    Gphase(:,ii+1) = Gphasei*cos(ang_rot) - Greadi*sin(ang_rot);
    Gslice(:,ii+1) = Gslicei;
end

for ii = 0:(obj.nShotsUsed-1)
    ang_rot = phi*(ii-(obj.nShotsUsed-1)*floor(ii/obj.nShotsUsed));
    obj.ky(:,ii+1) = kreadN*cos(ang_rot) + kphaseN*sin(ang_rot);
    obj.kx(:,ii+1) = kphaseN*cos(ang_rot) - kreadN*sin(ang_rot);
end

%Convert Gx,Gy to physical gradients
invRotMatrix = inv(obj.rotMatrix);
for ii = 0:(obj.nShotsUsed-1)
    for jj = 1:length(Greadi)
        temp = invRotMatrix*[Gphase(jj,ii+1);Gread(jj,ii+1);Gslice(jj,ii+1)];
        Gx(jj,ii+1) = temp(1);
        Gy(jj,ii+1) = temp(2);
        Gz(jj,ii+1) = temp(3);
    end
end

%Correct with GIRF
for ii = 0:(obj.nShotsUsed-1)
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

% Return to logical indices
for ii = 0:(obj.nShotsUsed-1)
    for jj = 1:length(Greadi)
        temp = obj.rotMatrix*[Gxc(jj,ii+1);Gyc(jj,ii+1);Gzc(jj,ii+1)];
        Gphasec(jj,ii+1) = temp(1);
        Greadc(jj,ii+1) = temp(2);
        Gslicec(jj,ii+1) = temp(3);
    end
end

if obj.nPartitions > 1
    kSlicec = -obj.nPartitions/2:(obj.nPartitions/2-1);
else
    kSlicec = 0;
end

% attempting to handle the multibandFactor / kz undersampling here, while
% staying general
% CLJ 9/10/2017
% NOTE: this probably only works for 4/2, 8/4, or 8/2 undersampling cases,
% just like the sequence. Will figure out how to be more flexible at next
% gen
if obj.multibandFactor > 1
    kSlicec = (-obj.multibandFactor/2):(obj.multibandFactor/obj.nPartitions):(obj.multibandFactor/2-1);
end


gamma = 2*pi*4.257e3;
gambar = gamma/(2*pi);

for ii = 1:(obj.nShotsUsed)
    kPhaset(:,ii) = cumsum([Gphase(:,ii)])*girfTs*obj.FOV*gambar;
    kReadt(:,ii) = cumsum([Gread(:,ii)])*girfTs*obj.FOV*gambar;
    kSlicet(:,ii) = cumsum([Gslicec(:,ii)])*girfTs*obj.FOV*gambar;
end
    
for ii = 1:(obj.nShotsUsed)
    kReadc(:,ii) = interp1(0:girfTs:girfTs*length(kReadt(:,ii))-girfTs,kReadt(:,ii),0:obj.adcTs:girfTs*length(kReadt(:,ii))-obj.adcTs)';
    kPhasec(:,ii) = interp1(0:girfTs:girfTs*length(kPhaset(:,ii))-girfTs,kPhaset(:,ii),0:obj.adcTs:girfTs*length(kPhaset(:,ii))-obj.adcTs)';
    %kzc(:,ii) = interp1(0:girfTs:girfTs*length(kzt(:,ii))-girfTs,kzt(:,ii),0:obj.adcTs:girfTs*length(kzt(:,ii))-obj.adcTs)';
end

if obj.useGIRF
    obj.shotLength = length(obj.kx(:,1)); %Assume kxt and kyt are the same length
    %obj.shotLength = 2000; %Assume kxt and kyt are the same length
    obj.kx = repmat(kReadc(1:obj.shotLength,:),[1,1,obj.nPartitions, obj.nSlices,obj.nAverages, obj.nPhases, obj.nEchoes, obj.nRepetitions,1]);
    obj.ky = repmat(kPhasec(1:obj.shotLength,:),[1,1,obj.nPartitions, obj.nSlices,obj.nAverages, obj.nPhases, obj.nEchoes, obj.nRepetitions,1]);
    kSlicec = permute(kSlicec,[3,1,2]);
    obj.kSlice = repmat(kSlicec,[length(kReadc(1:obj.shotLength,1)),obj.nShots,1, obj.nSlices,obj.nAverages, obj.nPhases, obj.nEchoes, obj.nRepetitions,1]);
    %obj.kz = reshape(obj.kz,[],obj.nPartitions);
else % Do not use GIRF, use kxNominal
    obj.shotLength = length(obj.kx(:,1)); %Assume kxt and kyt are the same length
    obj.kx = repmat(obj.kx(1:obj.shotLength,:),[1,1,obj.nPartitions, obj.nSlices,obj.nAverages, obj.nPhases, obj.nEchoes, obj.nRepetitions,1]);
    obj.ky = repmat(obj.ky(1:obj.shotLength,:),[1,1,obj.nPartitions, obj.nSlices,obj.nAverages, obj.nPhases, obj.nEchoes, obj.nRepetitions,1]);
    kSlicec = permute(kSlicec,[3,1,2]);
    obj.kz = repmat(kSlicec,[length(kReadc(1:obj.shotLength,1)),obj.nShots,1, obj.nSlices,obj.nAverages, obj.nPhases, obj.nEchoes, obj.nRepetitions,1]);
    %obj.kx = repmat(col(obj.kxNominal(1:obj.shotLength,:)),1,obj.nPartitions);
    %obj.ky = repmat(col(obj.kyNominal(1:obj.shotLength,:)),1,obj.nPartitions);
end


% Assume that the density compensation function from the first partition
% can be used for all partitions. 
obj.ww = weight_vor(col(obj.kx(:,:,1,1,1,1,1,1,1)),col(obj.ky(:,:,1,1,1,1,1,1,1)),obj.nShotsDesigned,1);
obj.ww(obj.ww > 2) = 2;
obj.ww = reshape(obj.ww,[],obj.nShots);
obj.ww = repmat(obj.ww,[1,1,obj.nPartitions, obj.nSlices,obj.nAverages, obj.nPhases, obj.nEchoes, obj.nRepetitions,1]);
timingOffset = obj.TE(1)*1E-6;

for ii = 1:obj.nEchoes
    obj.timingVec(:,ii) = col(obj.TE(ii)*1E-6:obj.adcTs:(obj.adcTs*obj.shotLength+obj.TE(ii)*1E-6-obj.adcTs));
end

obj.timingVec = obj.timingVec - timingOffset;
obj.timingVec  = reshape(obj.timingVec,[],1,1,1,1,1,obj.nEchoes,1,1);
obj.timingVec = repmat(obj.timingVec,[1,obj.nShots,obj.nPartitions, obj.nSlices,obj.nAverages, obj.nPhases, 1, obj.nRepetitions,1]);


end