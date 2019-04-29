function [grads_mTm ] = readExternalMultishotGradFile( filename )
%READEXTERNALGRADFILE Summary of this function goes here
%   Detailed explanation goes here

xmlStruct = xml2struct(filename);

waveformLength = str2num(xmlStruct(1).Attributes(6).Value);
refGradAmp = str2num(xmlStruct(1).Attributes(4).Value);
numberOfShots = str2num(xmlStruct(1).Attributes(2).Value);
numberOfTrajectories = str2num(xmlStruct(1).Attributes(3).Value);
refFOV = str2num(xmlStruct(1).Attributes(5).Value);

% With the multishot setup, we can't just step through waveform length by
% itself any more

raw_grads_mTm = zeros(waveformLength,3,numberOfTrajectories);

% Handle the number of raw trajectories
for jj = 1:numberOfTrajectories
    for ii = 1:waveformLength

        raw_grads_mTm(ii,1,jj) = refGradAmp*str2num(xmlStruct(1).Children(2*jj).Children(2*ii).Attributes(2).Value);
        raw_grads_mTm(ii,2,jj) = refGradAmp*str2num(xmlStruct(1).Children(2*jj).Children(2*ii).Attributes(1).Value);
        raw_grads_mTm(ii,3,jj) = refGradAmp*str2num(xmlStruct(1).Children(2*jj).Children(2*ii).Attributes(3).Value);
    
    end
end

% Now we need to take the raw trajectories and turn them into the actual
% shots that were acquired.
grads_mTm = zeros(waveformLength,3,numberOfShots);

FirstShotIndex = 2*(numberOfTrajectories);

for jj = 1:numberOfShots
    % Get Shot Info 
    
    % Convert from zero index arrays (C/C++) to one index arrays (Matlab)
    trajIndex = str2num(xmlStruct(1).Children(FirstShotIndex+2*jj).Attributes(1).Value) + 1;
    rotAngle = str2num(xmlStruct(1).Children(FirstShotIndex+2*jj).Attributes(2).Value);
    phaseGradScale = str2num(xmlStruct(1).Children(FirstShotIndex+2*jj).Attributes(3).Value);
    readGradScale = str2num(xmlStruct(1).Children(FirstShotIndex+2*jj).Attributes(4).Value);
    sliceGradScale = str2num(xmlStruct(1).Children(FirstShotIndex+2*jj).Attributes(5).Value);
    zDir = str2num(xmlStruct(1).Children(FirstShotIndex+2*jj).Attributes(6).Value);
    % Calculate rotation matrix and scaling vector
    
    rotMatrix = [cos(rotAngle), -sin(rotAngle), 0; sin(rotAngle), cos(rotAngle), 0; 0 0 1];
    scaleVector = [-1*phaseGradScale; -1*readGradScale; sliceGradScale*zDir;];
    
    for ii = 1:waveformLength
        grads_mTm(ii,:,jj) = scaleVector.*(rotMatrix*col(raw_grads_mTm(ii,:,trajIndex)));
    end
end


end

