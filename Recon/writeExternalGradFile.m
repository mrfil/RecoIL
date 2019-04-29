function [grads_mTm ] = writeExternalGradFile(xmlFileName, grads_mTm,FOV,FOVz,centerOfKSpace,BaseResolution)
%writeExternalGradFile Summary of this function goes here
%   Detailed explanation goes here

%% Write to XML File for moving to scanner
docNode = com.mathworks.xml.XMLUtils.createDocument('ExternalGradWaveform');
docRootNode = docNode.getDocumentElement;
docRootNode.setAttribute('WaveformLength',num2str(length(grads_mTm)));
docRootNode.setAttribute('ReferenceAmplitudemTm',num2str(max(col(abs(grads_mTm)))));
docRootNode.setAttribute('WaveformLength',num2str(length(grads_mTm)));
docRootNode.setAttribute('ReferenceFieldOfViewmm',num2str(FOV));
docRootNode.setAttribute('ReferenceFieldOfViewmmZ',num2str(FOVz));
docRootNode.setAttribute('CenterOfKSpaceIndex', num2str(centerOfKSpace));
docRootNode.setAttribute('BaseResolution', num2str(BaseResolution));


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