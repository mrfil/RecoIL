function obj = parseFlxTPI1(obj)

% load kspace.mat
% 
% 
% 
% obj.kx = interp1(0:obj.gradTs:obj.gradTs*length(kspace(:,1))-obj.gradTs,kspace(:,1),0:obj.adcTs:obj.gradTs*length(kspace(:,1))-obj.gradTs)';
% obj.ky = interp1(0:obj.gradTs:obj.gradTs*length(kspace(:,2))-obj.gradTs,kspace(:,2),0:obj.adcTs:obj.gradTs*length(kspace(:,2))-obj.gradTs).';
% obj.kz = interp1(0:obj.gradTs:obj.gradTs*length(kspace(:,3))-obj.gradTs,kspace(:,3),0:obj.adcTs:obj.gradTs*length(kspace(:,3))-obj.gradTs).';
% obj.ShotLength = length(obj.kx);
% obj.ww = weight_vor(col(obj.kx),col(obj.ky),obj.nShotsDesigned,1);
% obj.ww(obj.ww > 2) = 2;
% 
% [obj.readShift,obj.phaseShift,obj.sliceShift] = imageShift(obj,1);
% obj.kx = repmat(kxc(:),1,obj.nPartitions);
% obj.ky = repmat(kyc(:),1,obj.nPartitions);
% obj.kz = col(repmat(kzc,length(kxt(:))*obj.nShots,1));
% obj.kz = reshape(obj.kz,[],obj.nPartitions);
% obj.ww = weight_vor(col(obj.kx),col(obj.ky),obj.nShotsDesigned,1);
% obj.ww(obj.ww > 2) = 2;
% obj.ww = reshape(obj.ww,[],obj.nPartitions);
obj.adcTs = obj.DataObject.hdr.Dicom.DICOM.alICEProgramPara(9)*1E-6; %This could and should be added to the alICEProgramPara in the SpiralSBB class
obj.ptsToDrop = obj.DataObject.hdr.Dicom.DICOM.alICEProgramPara(8); %Empirical factor

%load kspace.mat
% %load grads.mat
grads_mTm = readExternalMultishotGradFile('test.xml');

sizeGrads = size(grads_mTm);
obj.nShots = sizeGrads(3); % Shots are along the 3rd direction returned from the read function
obj.FOV = 22; % Ugly Hack, temporary.
for ii=1:obj.nShots
    kspace = calcKSpaceFromGrads(grads_mTm(:,:,ii),obj.adcTs,obj.FOV*10);
    obj.kx(:,ii) = kspace(:,1);
    obj.ky(:,ii) = kspace(:,2);
    obj.kz(:,ii) = kspace(:,3);
end

obj.ShotLength = length(obj.kx(:,1)); %Assume kxt and kyt are the same length
obj.dataMask = ones(obj.ShotLength,1);
% obj.kx = repmat(kxc(:),1,obj.nPartitions);
% obj.ky = repmat(kyc(:),1,obj.nPartitions);
% %obj.kz = col(repmat(kzc,length(kxc(:,1))*obj.nShots,1));
% obj.kz = repmat(kzc,1,obj.nPartitions);
% 
% %Find center of kspace
% [~,loc] = min(sum(sqrt(kspace(10:end,:).^2),2));
% 
% timingOffset = loc*obj.gradTs;

obj.timingVec = 0:obj.adcTs:obj.adcTs*(length(obj.kx)-1);
% 
% obj.timingVec = obj.timingVec - timingOffset;

end