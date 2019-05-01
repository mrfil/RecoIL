% Start Here: an introductory guide to RecoIL!
% Alex Cerjanic
% 2019-04-29

% Tested with MATLAB R2017b

clear
%% Specify some basic factors

% Image Size
N = 64;

% 2D imaging
Nz = 1;

% Single shot imaging
nShots = 1;

% Define number of coils
nCoils = 4;

% Oversampling factor for ground truth data (to avoid the inverse crime)
Nos = 2;

% SNR for simulation
SNR = 50;

% Set sampling/raster time to be the same for simplicity
Ts = 4E-6;

%% Generate a phantom

% Generate a T1 weighted brain slice phantom from data included with matlab
load mri.mat D;
imgGroundTruth = double(D(:,:,15))./double(max(col(D(:,:,15))));


% Generate a sense map
sen_map = generate_birdcage_sensitivities(Nos*N, nCoils, 1.5);

% Reshape sense map into the format we like (single vector image by number
% of coils)
sen_map = reshape(sen_map,[],nCoils);

% Generate a field map with max 100 Hz off-resonance
FM_groundTruth = genSynthFieldMap(Nos*N,100);


%% Let's simulate the calibration scan acquisition (senFM Map in our parlance).

% Generate 18 shot trajectory to limit impacts of field map and distortion

nShotsSenFM = 18;

[kReadSenFM, kPhaseSenFM] = generateSpiral(24, N, 2.2, 120, Ts, nShotsSenFM);

% Generate timing vectors on a per shot basis
timingVecSen_SpinEcho = 0:Ts:(length(kReadSenFM(:,1))-1) * Ts;
timingVecSen_AsymSpinEcho = timingVecSen_SpinEcho + 1E-3;

% Generate first echo (spiral out starting at the spin echo)
G_forwardModel = Gdft(kReadSenFM(:),kPhaseSenFM(:),zeros(size(kReadSenFM)),Nos*N,Nos*N,1*Nz, 2*pi*FM_groundTruth, timingVecSen_SpinEcho);
%G_forwardModel = NUFFT(kReadSenFM(:),kPhaseSenFM(:),zeros(size(kReadSenFM)),Nos*N,Nos*N,1);
S_forwardModel = sense(G_forwardModel,sen_map);

data_firstEcho = S_forwardModel * col(imgGroundTruth);
% Add noise for a good simulation
noisyData_firstEcho = data_firstEcho + 1/SNR.*(randn(size(data_firstEcho)) + 1j*rand(size(data_firstEcho)));

% Generate second echo (spiral out starting at 1ms after the spin echo to
% allow effects of inhomogeneity to develop with little R2* effects)
G_forwardModel = Gdft(kReadSenFM(:),kPhaseSenFM(:),zeros(size(kReadSenFM)),Nos*N,Nos*N,1*Nz, 2*pi*FM_groundTruth, timingVecSen_AsymSpinEcho);
%G_forwardModel = NUFFT(kReadSenFM(:), kPhaseSenFM(:), zeros(size(kReadSenFM)), Nos*N, Nos*N, 1);
S_forwardModel = sense(G_forwardModel,sen_map);

data_secondEcho = S_forwardModel * col(imgGroundTruth);

% Add noise for a good simulation
noisyData_secondEcho = data_secondEcho + 1/SNR.*(randn(size(data_secondEcho)) + 1j*rand(size(data_secondEcho)));

%% Now we plug this data into a recoInfoSim object to use RecoIL recon functions

% Reshape all of the data for recoInfoSim.

% Deal with kspace first
kReadSenFM_rInfo = repmat(kReadSenFM,1,1,1,1,1,1,2,1,1);
kPhaseSenFM_rInfo = repmat(kPhaseSenFM,1,1,1,1,1,1,2,1,1);

% Deal with timing vector
timingVec_Echo1 = repmat(col(timingVecSen_SpinEcho),1,nShotsSenFM,1,1,1,1,1,1,1,1);
timingVec_Echo2 = repmat(col(timingVecSen_AsymSpinEcho),1,nShotsSenFM,1,1,1,1,1,1,1,1);
timingVec_rInfo = cat(7,timingVec_Echo1,timingVec_Echo2);

% Deal with assembling data matrix
noisyData_firstEcho = reshape(noisyData_firstEcho,[],nShotsSenFM,nCoils,1,1,1,1,1,1,1,1);
noisyData_secondEcho = reshape(noisyData_secondEcho,[],nShotsSenFM,nCoils,1,1,1,1,1,1,1,1);

noisyData_firstEcho = permute(noisyData_firstEcho,[1,2,4,5,6,7,8,9,3]);
noisyData_secondEcho = permute(noisyData_secondEcho,[1,2,4,5,6,7,8,9,3]);

noisyData = cat(7,noisyData_firstEcho, noisyData_secondEcho);

% Note that TE is in microseconds (siemens convention)
rInfoSenFM = recoInfoSim(kReadSenFM_rInfo,kPhaseSenFM_rInfo, ...
                         zeros(size(kReadSenFM_rInfo)),noisyData, ...
                         timingVec_rInfo, N, Nz, 'TE', [0, 1000]);

% Calculate density compensation function for multishot spiral as used in
% this example.

rInfoSenFM.ww = calc3DDensityCompensation(col(rInfoSenFM.kRead(:,:,1,1,1,1,1,1,1)),col(rInfoSenFM.kPhase(:,:,1,1,1,1,1,1,1)),col(rInfoSenFM.kSlice(:,:,1,1,1,1,1,1,1)),rInfoSenFM.N,rInfoSenFM.nPartitions);
rInfoSenFM.ww = reshape(rInfoSenFM.ww,size(rInfoSenFM.kRead(:,:,1,1,1,1,1,1,1)));
rInfoSenFM.ww = repmat(rInfoSenFM.ww,[1, 1, rInfoSenFM.nPartitions, rInfoSenFM.nSlices, rInfoSenFM.nAverages, rInfoSenFM.nPhases, rInfoSenFM.nEchoes, rInfoSenFM.nRepetitions,1]);

%% Let's simulate the calibration scan reconstruction

% Start by gridding the images
% Since recoInfo has all of the information in it, as well as a function to
% access our data, gridCoilImages() only needs rInfoSenFM. 
cImages = gridCoilImages(rInfoSenFM);
% We work in a set of standard dimensions for convenience, so coils is
% along the 10th dimension in this data. 

% Form the sum of squares (SOS) image
imSOS = sqrt(sum(abs(cImages).^2,10));

% Create a sense map
% Input has to be N x N x NSlices x NCoils, so we use the first echo
cImages = squeeze(cImages);
[sen, mask] = createSenMap(cImages(:,:,1,:), 2);
sen = reshape(sen,N,N,1,1,nCoils);

% Calculate a field map
[FM, FMImages, mask] = createFieldMap(rInfoSenFM,sen,mask,1,'nIterations',1);

%% Now we can simulate a reconstruction using the results of our calibration scans.


% Let's see what happens in the single shot spiral case. A 64 matrix single
% shot with R=2 is pretty do-able with field correction.

% Generate synthetic data using Gdft (again to avoid the inverse crime).
Rfactor = 2;
nShots = 1;
SNR_gre = 30;

% New trajectory
[kx, ky] = generateSpiral(24, N, 2.2, 120, Ts, Rfactor*nShots);

kRead  = kx(:,1);
kPhase = ky(:,1);

% Simulate a gradient echo acquisition with TE = 15ms.

timingVec_gre = 15E-3:Ts:((length(kRead)-1) * Ts + 15E-3);

% Generate data
G_forwardModel = Gdft(kRead,kPhase,zeros(size(kRead)),Nos*N,Nos*N,1*Nz, 2*pi*FM_groundTruth, col(timingVec_gre));
S_forwardModel = sense(G_forwardModel,sen_map);

data_gre = S_forwardModel * col(imgGroundTruth);
% Add noise for a good simulation
noisyData_gre = data_gre + 1/SNR_gre.*(randn(size(data_gre)) + 1j*rand(size(data_gre)));

noisyData_gre = reshape(noisyData_gre,[],nShots,nCoils,1,1,1,1,1,1,1,1);

noisyData_gre = permute(noisyData_gre,[1,2,4,5,6,7,8,9,3]);

% Generate new recoInfo object with the new single shot trajectory and
% timing vector
rInfo = recoInfoSim(kRead,kPhase, ...
                         zeros(size(kRead)),noisyData_gre, ...
                         col(timingVec_gre), N, Nz, 'TE', [15000]);

% Setup a reconstruction with some parameters
% Beta for the regularization constant
Rbeta = 10;
% Number of CG iterations for reconstruction
Niter = 10;

% Set the number of time segments, one for ever 2ms of readout time
L = ceil(length(kRead)*Ts/2E-3);     

% Setup the NUFFT with the trajectory stored in recoInfo. The transform is
% masked with mask.
G = NUFFT(col(rInfo.kRead(:)), col(rInfo.kPhase(:)), col(rInfo.kSlice(:)), ...
    rInfo.N, rInfo.N, 1, 'mask', logical(mask));

% Setup a time segmentation object
A = TimeSegmentation(G, col(timingVec_gre), col(FM(squeeze(mask))), L);

% Setup a quadratic regularization object with Rbeta regularization
% constant with 2d penalization.
R = Robj(logical(squeeze(mask)), 'edge_type', 'tight', 'order', 2,...
    'beta', Rbeta, 'type_denom', 'matlab', 'potential', 'quad',...
    'dims2penalize', [1,1,0]);

% Reshape the sensemap into an [N*N, nCoils] matrix
sen_tmp = reshape(sen,[],nCoils);

% Create the SENSE object.
S = sense(A,sen_tmp(logical(col(mask)),:));

% Grab the data (only the first data) with all of the shots. We use the []
% empty matrix to signifiy all of the data ala colon operator. Can't pass
% the colon operator to a function.
data = col(rInfo.dataRead([],[],1,1,1,1,1,1));

% Creating an array for the initial value for the image, all zeros.
imginit = zeros(N,N);

% We use a support mask (should be a circle) which matches the support of
% the spiral trajectory. Should be rectangular for a different trajectory.
xinit = col(imginit(squeeze(mask)));

% Run the penalized weighted least squares CG optimizer for 10 iterations.

[ colImage, ~, resid] = solve_pwls_pcg(xinit, S, 1, data, R, 'niter', Niter);

% This function puts the masked data back to the full image.
imgFieldCorrected = embed(colImage, mask);
    
% Now we're going to setup a new recon object for the non-field corrected
% reconstruction.
S = sense(G,sen_tmp(logical(col(mask)),:));

[ colImage, ~, resid] = solve_pwls_pcg(xinit, S, 1, data, R, 'niter', Niter);

imgNoFieldCorrection = embed(colImage, mask);
