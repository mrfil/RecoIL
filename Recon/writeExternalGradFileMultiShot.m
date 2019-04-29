function writeExternalGradFileMultiShot(xmlFileName, grads_mTm, shotSchedule, FOV, FOVz, centerOfKSpace)
%writeExternalGradFileMultiShot Write Gradient Waveforms for use with
%Multishot External Trajectory SBB
%
%
%   Inputs:
%   xmlFileName    - Character Array with File Name. (Hint: Append with
%   .xml)
%   grads_mTm      - [WaveformLength x 3 x nTrajectories] Gradient Waveforms 
%   shotSchedule   - [NShots x 2] - Array that corresponds to a shot
%                    schedule, where the first column corresponds to the
%                    Trajectory id (ie the gradient waveform corresponding
%                    to column #3 in grads_mTm), and the second is an
%                    in-plane rotation of the shot in radians.
%   FOV            - Field of View in mm
%   centerOfKSpace - Index ( < WaveformLength) that corresponds to where
%   the echo should form. (I.E. Spiral Out would set this to 0).
%
%
%
%   Outputs:
%   None

%% Write to XML File for moving to scanner
docNode = com.mathworks.xml.XMLUtils.createDocument('ExternalGradWaveform');
docRootNode = docNode.getDocumentElement;
docRootNode.setAttribute('WaveformLength',num2str(length(grads_mTm)));
docRootNode.setAttribute('NumberOfTrajectories',num2str(length(grads_mTm(1,1,:))));
docRootNode.setAttribute('NumberOfShots',num2str(length(shotSchedule(:,1))));
docRootNode.setAttribute('ReferenceAmplitudemTm',num2str(max(col(abs(grads_mTm)))));
docRootNode.setAttribute('ReferenceFieldOfViewmm',num2str(FOV));
docRootNode.setAttribute('ReferenceFieldOfViewmmZ',num2str(FOVz));
docRootNode.setAttribute('CenterOfKSpaceIndex', num2str(centerOfKSpace));

for jj = 1:length(grads_mTm(1,1,:))
    thisElement = docNode.createElement('GradientWaveform'); 
    thisElement.setAttribute('TrajID',num2str(jj-1));
    for ii=1:length(grads_mTm(:,1,1))
        thisPoint = docNode.createElement('GradWaveformPoint'); 
        thisPoint.setAttribute('Read',num2str(grads_mTm(ii,1)/max(col(abs(grads_mTm)))));
        thisPoint.setAttribute('Phase',num2str(grads_mTm(ii,2)/max(col(abs(grads_mTm)))));
        thisPoint.setAttribute('Slice',num2str(grads_mTm(ii,3)/max(col(abs(grads_mTm)))));
        thisElement.appendChild(thisPoint);
    end
    docRootNode.appendChild(thisElement);

end

for ii = 1:length(shotSchedule(:,1))
    thisElement = docNode.createElement('Shot');
    thisElement.setAttribute('ID',num2str(shotSchedule(ii,1)));
    thisElement.setAttribute('RotAngle',num2str(shotSchedule(ii,2)));
    docRootNode.appendChild(thisElement);
end
%xmlFileName = ['Navigator','.xml'];
xmlwrite(xmlFileName,docNode);
%type(xmlFileName);

end