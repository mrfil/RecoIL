function [grads_mTm ] = writeExternalGradFileMultiShot(xmlFileName, grads_mTm,FOV,centerOfKSpace)
%writeExternalGradFile Summary of this function goes here
%   Detailed explanation goes here

%% Write to XML File for moving to scanner
docNode = com.mathworks.xml.XMLUtils.createDocument('ExternalGradWaveform');
docRootNode = docNode.getDocumentElement;
docRootNode.setAttribute('WaveformLength',num2str(length(grads_mTm)));
docRootNode.setAttribute('NumberOfTrajectories',num2str(length(grads_mTm(1,1,:))));

docRootNode.setAttribute('ReferenceAmplitudemTm',num2str(max(col(abs(grads_mTm)))));
docRootNode.setAttribute('ReferenceFieldOfViewmm',num2str(FOV));
docRootNode.setAttribute('CenterOfKSpaceIndex', num2str(centerOfKSpace));
for ii=1:length(grads_mTm)
    thisElement = docNode.createElement('GradientWaveform'); 
    thisElement.setAttribute('Read',num2str(grads_mTm(ii,1)/max(col(abs(grads_mTm)))));
    thisElement.setAttribute('Phase',num2str(grads_mTm(ii,2)/max(col(abs(grads_mTm)))));
    thisElement.setAttribute('Slice',num2str(grads_mTm(ii,3)/max(col(abs(grads_mTm)))));
    docRootNode.appendChild(thisElement);
end
docNode.appendChild(docNode.createComment('this is a comment'));

%xmlFileName = ['Navigator','.xml'];
xmlwrite(xmlFileName,docNode);
%type(xmlFileName);

end