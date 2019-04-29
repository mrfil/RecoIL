function obj = parseSpFLASH(obj)

girfTs = 2E-6;

obj.nShotsDesigned = obj.DataObject.hdr.MeasYaps.sWipMemBlock.alFree{1};
obj.nShotsUsed = obj.DataObject.hdr.MeasYaps.sWipMemBlock.alFree{2};

obj = parseSpiralOutReadout(obj);


end
