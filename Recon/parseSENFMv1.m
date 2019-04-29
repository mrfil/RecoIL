function obj = parseSENFMv1(obj)

girfTs = 2E-6;

%Correct the TEs due to an error in the sequence.
for ii = 2:length(obj.TE)
    obj.TE(ii) = obj.TE(1)/2 + 3*obj.TE(ii)/2 - obj.TE(1);
end
obj.dataMask = true(obj.shotLength,1);

obj.nShotsDesigned = obj.DataObject.hdr.MeasYaps.sWipMemBlock.alFree{1};
obj.nShotsUsed = obj.DataObject.hdr.MeasYaps.sWipMemBlock.alFree{2};

obj = parseSpiralOutReadout(obj);


end