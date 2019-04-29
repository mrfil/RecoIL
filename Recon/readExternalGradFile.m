function [grads_mTm, refFOV, refFOVz, RUP, CenterOfKSpace] = readExternalGradFile(filename)
%READEXTERNALGRADFILE Summary of this function goes here
%   Detailed explanation goes here

xmlStruct = xml2struct(filename);

waveformLength = str2double(xmlStruct(1).ExternalGradWaveform.Attributes.WaveformLength);
refGradAmp = str2double(xmlStruct(1).ExternalGradWaveform.Attributes.ReferenceAmplitudemTm);
refFOV = str2double(xmlStruct(1).ExternalGradWaveform.Attributes.ReferenceFieldOfViewmm);
refFOVz = str2double(xmlStruct(1).ExternalGradWaveform.Attributes.ReferenceFieldOfViewmmZ);
grads_mTm = zeros(waveformLength,3);
if(isfield(xmlStruct(1).ExternalGradWaveform.Attributes,'RampUpPts'))
    RUP = str2double(xmlStruct(1).ExternalGradWaveform.Attributes.RampUpPts);
else
    RUP = 0;
end
if(isfield(xmlStruct(1).ExternalGradWaveform.Attributes,'CenterOfKSpaceIndex'))
    CenterOfKSpace = str2double(xmlStruct(1).ExternalGradWaveform.Attributes.CenterOfKSpaceIndex);
else
    CenterOfKSpace = 1;
end
%refGradAmpZ = refGradAmp/FOVz*refFOV;
 
for ii = 1:waveformLength

    grads_mTm(ii,1) = refGradAmp*str2num(xmlStruct(1).ExternalGradWaveform.GradientWaveform{ii}.Attributes.Read);
    grads_mTm(ii,2) = refGradAmp*str2num(xmlStruct(1).ExternalGradWaveform.GradientWaveform{ii}.Attributes.Phase);
    grads_mTm(ii,3) = refGradAmp*str2num(xmlStruct(1).ExternalGradWaveform.GradientWaveform{ii}.Attributes.Slice);
    
end


end

