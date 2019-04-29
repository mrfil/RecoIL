function [kspace] = calcKSpaceFromGrads(rInfo, grads_mTm, adcTs, FOV, varargin)
%calcKSpaceFromGrads Calculate k-space coordinates from gradient waveforms
%
%   Calculates a k-space trajectory from a gradient vector. This function
%   handles the interpolation from the gradient raster time based waveform
%   to the sampled trajectory at the ADC rate. This function does not
%   handle application of the Gradient Impulse Response Function (GIRF). If
%   you want your trajectory corrected by the GIRF, apply that to your
%   gradient waveforms before you call this function.
%
%   Input
%   grads_mTm     - [npts x 3] Array gradient waveforms in
%                    units of mT/m.
%   adcTs         - Scalar ADC sampling rate
%
%   Optional Inputs - specify as name - value pairs
%   'gts'         - Gradient Raster time in seconds (s). Default is 10E-6s.
%   'Nuclei'      - Nuclei being imaged as recognized by
%                         getGyromagneticRatio(). Default is '1H'.
%
%
%   Output
%

%% Deal with input arguments

p = inputParser();
RUP = 0;
gts = 10E-6; % Gradient raster time. 10E-6 seconds for Siemens,
% different on other platforms
GIRFts = 2E-6;
UseGIRF = false;

addOptional(p, 'Nuclei', '1H');
addOptional(p, 'gts', gts);
addOptional(p, 'FOVz', FOV);
addOptional(p, 'RUP', RUP);
addOptional(p, 'UseGIRF', UseGIRF);


parse(p,varargin{:});
Nuclei = p.Results.Nuclei;
gts = p.Results.gts;
FOVz = p.Results.FOVz;
RUP = p.Results.RUP;
UseGIRF = p.Results.UseGIRF;

gambar = getGammaBar(Nuclei); %Gamma Bar in units of Hz/T

%% Calculate k-space coordinates

if (UseGIRF == true)
    
    % Interpolate gradients to GIRF time samples (2microseconds)
    gradRead = interp1(0:gts:gts*(length(grads_mTm(:,1))-1),grads_mTm(:,1),0:GIRFts:gts*(length(grads_mTm(:,1))-1),'previous','extrap')';
    gradPhase = interp1(0:gts:gts*(length(grads_mTm(:,2))-1),grads_mTm(:,2),0:GIRFts:gts*(length(grads_mTm(:,2))-1),'previous','extrap')';
    gradSlice = interp1(0:gts:gts*(length(grads_mTm(:,3))-1),grads_mTm(:,3),0:GIRFts:gts*(length(grads_mTm(:,3))-1),'previous','extrap')';
    
    %Convert Gx,Gy to physical gradients
    invRotMatrix = inv(rInfo.rotMatrix);
    
    for jj = 1:length(gradRead)
        temp = invRotMatrix*[gradPhase(jj);gradRead(jj);gradSlice(jj)];
        Gx(jj) = temp(1);
        Gy(jj) = temp(2);
        Gz(jj) = temp(3);
    end
    
    % Debugging Test

    %Correct with GIRF
    Gxc(:) = conv(Gx(:),rInfo.GIRF(:,1),'same');
    Gyc(:) = conv(Gy(:),rInfo.GIRF(:,2),'same');
    Gzc(:) = conv(Gz(:),rInfo.GIRF(:,3),'same');
    
    % Return to logical indices
    
    for jj = 1:length(gradRead)
        temp = rInfo.rotMatrix*[Gxc(jj);Gyc(jj);Gzc(jj)];
        gradPhaseC(jj) = temp(1);
        gradReadC(jj) = temp(2);
        gradSliceC(jj) = temp(3);
    end
    
    
    gradX = interp1(0:GIRFts:GIRFts*(length(gradReadC)-1),gradReadC,0:adcTs:GIRFts*(length(gradReadC)-1),'previous','extrap')';
    gradY = interp1(0:GIRFts:GIRFts*(length(gradPhaseC)-1),gradPhaseC,0:adcTs:GIRFts*(length(gradPhaseC)-1),'previous','extrap')';
    gradZ = interp1(0:GIRFts:GIRFts*(length(gradSliceC)-1),gradSliceC,0:adcTs:GIRFts*(length(gradSliceC)-1),'previous','extrap')';
    
    kspace(:,1) = cumsum([gradX])*adcTs*FOV*gambar/1E6;
    kspace(:,2) = cumsum([gradY])*adcTs*FOV*gambar/1E6;
    kspace(:,3) = cumsum([gradZ])*adcTs*FOVz*gambar/1E6;
    
else
    
    gradX = interp1(0:gts:gts*(length(grads_mTm(:,1))-1),grads_mTm(:,1),0:adcTs:gts*(length(grads_mTm(:,1))-1),'previous','extrap')';
    gradY = interp1(0:gts:gts*(length(grads_mTm(:,2))-1),grads_mTm(:,2),0:adcTs:gts*(length(grads_mTm(:,2))-1),'previous','extrap')';
    gradZ = interp1(0:gts:gts*(length(grads_mTm(:,3))-1),grads_mTm(:,3),0:adcTs:gts*(length(grads_mTm(:,3))-1),'previous','extrap')';
    kspace(:,1) = cumsum([gradX])*adcTs*FOV*gambar/1E6;
    kspace(:,2) = cumsum([gradY])*adcTs*FOV*gambar/1E6;
    kspace(:,3) = cumsum([gradZ])*adcTs*FOVz*gambar/1E6;
    
end
%% Trim by RUP

%RUPPtsToDiscard = ceil((RUP*gts)/adcTs)+2;
RUPPtsToDiscard = ceil((RUP*gts)/adcTs);
if (RUPPtsToDiscard <= 0)
    RUPPtsToDiscard = 1;
end
kspace = kspace(RUPPtsToDiscard:end,:);

end

