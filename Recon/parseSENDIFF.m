function obj = parseSENFMv1(obj)
            obj.nShotsDesigned = obj.DataObject.hdr.MeasYaps.sWipMemBlock.alFree{1};
            obj.nShotsUsed = obj.DataObject.hdr.MeasYaps.sWipMemBlock.alFree{2};
            obj.ptsToDrop = obj.DataObject.hdr.Dicom.DICOM.alICEProgramPara(3) - 1; %Empirical factor
            
            obj.gradAmp = obj.DataObject.hdr.Dicom.DICOM.alICEProgramPara(2)/10; % was always 2.4
            obj.gradSlew = obj.DataObject.hdr.Dicom.DICOM.alICEProgramPara(1);  % was 120
            obj.adcTs = obj.DataObject.hdr.Dicom.DICOM.alICEProgramPara(4)*1E-6; %This could and should be added to the alICEProgramPara in the SpiralSBB class
            [~,~,kxi,kyi] = genspi(obj.FOV,obj.N,obj.nShotsDesigned,obj.gradAmp,obj.gradSlew,obj.gradTs); %[Gx,Gy,kxi,kyi,sx,sy]
            kxt=interp1(0:obj.gradTs:obj.gradTs*length(kxi)-obj.gradTs,kxi,0:obj.adcTs:obj.gradTs*length(kxi)-obj.gradTs)';
            kyt=interp1(0:obj.gradTs:obj.gradTs*length(kyi)-obj.gradTs,kyi,0:obj.adcTs:obj.gradTs*length(kyi)-obj.gradTs)';
            
            
            obj.ShotLength = length(kxt); %Assume kxt and kyt are the same length
            %sprintf('Performing %d rotations',obj.nShotsDesigned)
            phi = 2*pi/obj.nShotsUsed;
            for ii = 0:(obj.nShotsUsed-1)
                ang_rot = phi*(ii-(obj.nShotsUsed-1)*floor(ii/obj.nShotsUsed));
                kxc(:,ii+1) = kxt*cos(ang_rot) + kyt*sin(ang_rot);
                kyc(:,ii+1) = kyt*cos(ang_rot) - kxt*sin(ang_rot);
            end
            if obj.nPartitions > 1
                kzc = -obj.nPartitions/2:(obj.nPartitions/2-1);
            else
                kzc = 0;
            end
            
            obj.kx = repmat(kxc(:),1,obj.nPartitions);
            obj.ky = repmat(kyc(:),1,obj.nPartitions);
            obj.kz = col(repmat(kzc,length(kxt(:))*obj.nShots,1));
            obj.kz = reshape(obj.kz,[],obj.nPartitions);
            obj.ww = weight_vor(col(obj.kx),col(obj.ky),obj.nShotsDesigned,1);
            obj.ww(obj.ww > 2) = 2;
            obj.ww = reshape(obj.ww,[],obj.nPartitions);

end