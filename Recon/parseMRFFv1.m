function obj = parseMRFFv1(obj)

girfTs = 2E-6;

obj.nShotsDesigned = obj.DataObject.hdr.MeasYaps.sWipMemBlock.alFree{1};
obj.nShotsUsed = 1039; %FISP Parameters

obj.ptsToDrop = obj.DataObject.hdr.Dicom.DICOM.alICEProgramPara(3); %Empirical factor

obj.gradAmp = obj.DataObject.hdr.Dicom.DICOM.alICEProgramPara(2)/10; % was always 2.4
obj.gradSlew = obj.DataObject.hdr.Dicom.DICOM.alICEProgramPara(1);  % was 120
obj.adcTs = obj.DataObject.hdr.Dicom.DICOM.alICEProgramPara(4)*1E-6; 
[Greado,Gphaseo,kxi,kyi] = genspi(obj.FOV,obj.N,obj.nShotsDesigned,obj.gradAmp,obj.gradSlew,obj.gradTs); %[Gx,Gy,kxi,kyi,sx,sy]

Greadi  = interp1(0:10:10*length(Greado)-1,Greado,0:2:10*length(kxi)-1,'next','extrap')';
Gphasei = interp1(0:10:10*length(Gphaseo)-1,Gphaseo,0:2:10*length(kyi)-1,'next','extrap')';

Greadi  = [zeros(4,1); Greadi ];
Gphasei = [zeros(4,1); Gphasei];

Greadi(isnan(Greadi)) = 0;
Gphasei(isnan(Gphasei)) = 0;
Gslicei = zeros(size(Greadi));

kxN=interp1(0:obj.gradTs:obj.gradTs*length(kxi)-obj.gradTs,kxi,0:obj.adcTs:obj.gradTs*length(kxi)-obj.adcTs,'linear','extrap')';
kyN=interp1(0:obj.gradTs:obj.gradTs*length(kyi)-obj.gradTs,kyi,0:obj.adcTs:obj.gradTs*length(kyi)-obj.adcTs,'linear','extrap')';

%Gxt=interp1(0:obj.gradTs:obj.gradTs*length(kxi)-obj.gradTs,kxi,0:obj.adcTs:obj.gradTs*length(kxi)-obj.adcTs)';
%kxN = kxN(1:end-1);
%kyN = kyN(1:end-1);

%sprintf('Performing %d rotations',obj.nShotsDesigned)
phi = 2*pi/obj.nShotsUsed;

load GIRFx.mat
load GIRFy.mat
load GIRFz.mat

for ii = 0:(obj.nShotsUsed-1)
    ang_rot = 7.5*ii/180*pi;
    Gread(:,ii+1) = Greadi*cos(ang_rot) + Gphasei*sin(ang_rot);
    Gphase(:,ii+1) = Gphasei*cos(ang_rot) - Greadi*sin(ang_rot);
    Gslice(:,ii+1) = Gslicei;
end

for ii = 0:(obj.nShotsUsed-1)
    ang_rot = 7.5*ii/180*pi;
    obj.kxNominal(:,ii+1) = kxN*cos(ang_rot) + kyN*sin(ang_rot);
    obj.kyNominal(:,ii+1) = kyN*cos(ang_rot) - kxN*sin(ang_rot);
end

%Convert Gx,Gy to physical gradients
invRotMatrix = inv(obj.rotMatrix);
for ii = 0:(obj.nShotsUsed-1)
    for jj = 1:length(Greadi)
        temp = invRotMatrix*[Gread(jj,ii+1);Gphase(jj,ii+1);Gslice(jj,ii+1)];
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
        Greadc(jj,ii+1) = temp(1);
        Gphasec(jj,ii+1) = temp(2);
        Gslicec(jj,ii+1) = temp(3);
    end
end

if obj.nPartitions > 1
    kzc = -obj.nPartitions/2:(obj.nPartitions/2-1);
else
    kzc = 0;
end

gamma = 2*pi*4.257e3;
gambar = gamma/(2*pi);

for ii = 1:(obj.nShotsUsed)
    kxt(:,ii) = cumsum([Greadc(:,ii)])*girfTs*obj.FOV*gambar;
    kyt(:,ii) = cumsum([Gphasec(:,ii)])*girfTs*obj.FOV*gambar;
    kzt(:,ii) = cumsum([Gslicec(:,ii)])*girfTs*obj.FOV*gambar;
end
    
for ii = 1:(obj.nShotsUsed)
    kxc(:,ii) = interp1(0:girfTs:girfTs*length(kxt(:,ii))-girfTs,kxt(:,ii),0:obj.adcTs:girfTs*length(kxt(:,ii))-obj.adcTs)';
    kyc(:,ii) = interp1(0:girfTs:girfTs*length(kyt(:,ii))-girfTs,kyt(:,ii),0:obj.adcTs:girfTs*length(kyt(:,ii))-obj.adcTs)';
    %kzc(:,ii) = interp1(0:girfTs:girfTs*length(kzt(:,ii))-girfTs,kzt(:,ii),0:obj.adcTs:girfTs*length(kzt(:,ii))-obj.adcTs)';
end

if obj.useGIRF
    obj.ShotLength = length(obj.kxNominal(:,1)); %Assume kxt and kyt are the same length
    %obj.ShotLength = 2000; %Assume kxt and kyt are the same length
    obj.kx = repmat(kxc(1:obj.ShotLength,:),[1,1,obj.nPartitions, obj.nSlices,obj.nAverages, obj.nPhases, obj.nEchoes, obj.nRepetitions,1]);
    obj.ky = repmat(kyc(1:obj.ShotLength,:),[1,1,obj.nPartitions, obj.nSlices,obj.nAverages, obj.nPhases, obj.nEchoes, obj.nRepetitions,1]);
    kzc = permute(kzc,[3,1,2]);
    obj.kz = repmat(kzc,[length(kxc(1:obj.ShotLength,1)),obj.nShots,1, obj.nSlices,obj.nAverages, obj.nPhases, obj.nEchoes, obj.nRepetitions,1]);
    %obj.kz = reshape(obj.kz,[],obj.nPartitions);
else % Do not use GIRF, use kxNominal
    obj.ShotLength = length(obj.kxNominal(:,1)); %Assume kxt and kyt are the same length
    obj.kx = reshape(obj.kxNominal(1:obj.ShotLength,:),[obj.ShotLength,1,obj.nPartitions, obj.nSlices,obj.nAverages, obj.nPhases, obj.nShotsUsed, obj.nRepetitions,1]);
    obj.ky = reshape(obj.kyNominal(1:obj.ShotLength,:),[obj.ShotLength,1,obj.nPartitions, obj.nSlices,obj.nAverages, obj.nPhases, obj.nShotsUsed, obj.nRepetitions,1]);
    kzc = permute(kzc,[3,1,2]);
    obj.kz = repmat(kzc,[length(kxc(1:obj.ShotLength,1)),1,1, obj.nSlices,obj.nAverages, obj.nPhases, obj.nShotsUsed, obj.nRepetitions,1]);
    %obj.kx = repmat(col(obj.kxNominal(1:obj.ShotLength,:)),1,obj.nPartitions);
    %obj.ky = repmat(col(obj.kyNominal(1:obj.ShotLength,:)),1,obj.nPartitions);
end


% Assume that the density compensation function from the first partition
% can be used for all partitions. 
obj.ww = weight_vor(col(obj.kx(:,1,1,1,1,1,1,1,1)),col(obj.ky(:,:,1,1,1,1,1,1,1)),1,1);
obj.ww(obj.ww > 2) = 2;
obj.ww = col(obj.ww);
obj.ww = repmat(obj.ww,[1,1,obj.nPartitions, obj.nSlices,obj.nAverages, obj.nPhases, obj.nShotsUsed, obj.nRepetitions,1]);
timingOffset = obj.TE(1)*1E-6;

%for ii = 1:obj.nEchoes
    obj.timingVec = col(obj.TE(1)*1E-6:obj.adcTs:(obj.adcTs*obj.ShotLength+obj.TE(1)*1E-6-obj.adcTs));
%end

obj.timingVec = obj.timingVec - timingOffset;
%obj.timingVec  = reshape(obj.timingVec,[],1,1,1,1,1,obj.nShotsUsed,1,1);
obj.timingVec = repmat(obj.timingVec,[1,1,obj.nPartitions, obj.nSlices,obj.nAverages, obj.nPhases, obj.nShotsUsed, obj.nRepetitions,1]);


end
