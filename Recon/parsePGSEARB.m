function obj = parsePGSEARB(obj)
%% Parse sequence specific data

obj.nShotsDesigned = obj.DataObject.hdr.MeasYaps.sWipMemBlock.alFree{1};
obj.nShotsUsed = obj.DataObject.hdr.MeasYaps.sWipMemBlock.alFree{2};

EchoGap = obj.DataObject.hdr.Dicom.DICOM.alICEProgramPara(13); % Echo Gap between navigator and the echo time in usec.
obj = parseExternalNavgiator(obj);

%% Deal with imaging readout here.

obj = parseSpiralOutReadout(obj);

obj.halfVoxelShiftRead = 1/(2*obj.N);
obj.halfVoxelShiftPhase = 1/(2*obj.N);
if obj.multibandFactor > 1
    obj.halfVoxelShiftSlice = 1/(2*obj.multibandFactor);
else
    obj.halfVoxelShiftSlice = 1/(2*obj.nPartitions);
end

end