function [grads] = calcGradsFromKSpace(kspace, adcTs, FOV, varargin)
%calcGradsFromKSpace Calculate gradient waveform from k-space
%
%   Calculates a gradient waveforms from k-space trajectories. This function does not
%   handle application of the Gradient Impulse Response Function (GIRF). If
%   you want your trajectory corrected by the GIRF, apply that to your
%   gradient waveforms before you call this function.
%
%   Input
%   kspace     - [npts x 3] Array unitless k-space from [-N/2,-N/2]
%   adcTs      - Scalar ADC sampling rate
%   FOV        - Field of View in mm
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
    
    gts = 10E-6; % Gradient raster time. 10E-6 seconds for Siemens, 
    % different on other platforms
    
    addOptional(p, 'Nuclei', '1H');
    addOptional(p, 'gts', gts);
    parse(p,varargin{:});
    Nuclei = p.Results.Nuclei;
    gts = p.Results.gts;
    
    gambar = getGammaBar(Nuclei); %Gamma Bar in units of Hz/T
    
%% Calculate k-space coordinates
    
    gradX = diff([0; col(kspace(:,1))])/(adcTs*FOV*gambar)*1E6;
    gradY = diff([0; col(kspace(:,2))])/(adcTs*FOV*gambar)*1E6;
    gradZ = diff([0; col(kspace(:,3))])/(adcTs*FOV*gambar)*1E6;


    gradX = interp1(0:gts:gts*(length(gradX(:))-1),gradX(:),0:adcTs:gts*(length(gradX(:))-1),'previous','extrap')';
    gradY = interp1(0:gts:gts*(length(gradY(:))-1),gradY(:),0:adcTs:gts*(length(gradY(:))-1),'previous','extrap')';
    gradZ = interp1(0:gts:gts*(length(gradZ(:))-1),gradZ(:),0:adcTs:gts*(length(gradZ(:))-1),'previous','extrap')';
    
    grads = [gradX(:),gradY(:),gradZ(:)];
    
    %kspace(:,1) = cumsum([0; gradX])*adcTs*FOV*gambar/1E6;
    %kspace(:,2) = cumsum([0; gradY])*adcTs*FOV*gambar/1E6;
    %kspace(:,3) = cumsum([0; gradZ])*adcTs*FOV*gambar/1E6;
    
end

