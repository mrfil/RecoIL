function reconMultibandDWI(datadir,filename,varargin)
% function reconMultibandDWI(datadir,filename,varargin)
% 
% Inputs (optional):
%
% - datadir:  string pointing to data directory
%             (if not entered, will prep current directory)
%
% - filename: string with filename for recon'ed image file
%             (if not entered, will use "img" as default)
%
%
% Optional Name - Value Pair Inputs:
%    Rbeta   - Add beta for regularization function. Default Beta is zero
%              which is not recommended. (0 is no regularization)
%
%    L       - Overrides the number of time segments used. Otherwise, grabs
%              the value from the recoInfo object rInfo. (0 is no field
%              correction)
%
%   Niter    - Sets the number of CG iterations used to solve for the image
%              by solve_pwls_pcg.m. Default is 10.
%
%   slicesToRecon - Vector of slices to recon. Default is all slices
%
%   repetitionsToRecon - Vector of repetitions to recon. Default is all repetitions  
%
%
% Authors:
% Curtis Johnson - University of Delaware
% Alex Cerjanic - University of Illinois at Urbana-Champaign
% Nov 2017
%
% Note: Built out of the former scratch script "reconMultibandMRE" created
% during MB sequence/recon development. Cleaned and consolidated by CLJ.

curdir = pwd;

if ~isempty(datadir)
    cd(datadir)
end

rInfo = recoInfo();

if ~(exist('sen.mat','file')&&exist('FM.mat','file')&&exist('mask.mat','file')&&exist('PMaps.mat','file'))
    error('Missing prep files; run prepMultibandDWI first')
end
load sen.mat
load FM.mat
load mask.mat
load PMaps.mat


% Deal with Optional Inputs

% Establish Defaults
L = rInfo.L;
LValidationFcn = @(x) isnumeric(x) && isscalar(x) && (x >= 0);

Rbeta = 1; %Use a zero beta, not a good idea in practice.
RBetaValidationFcn = @(x) isnumeric(x) && isscalar(x) && (x >= 0);

Niter = 15;
NIterValidationFcn = @(x) isnumeric(x) && isscalar(x) && (x > 0);

slicesToRecon = 1:rInfo.nSlices;
slicesToReconValidationFcn = @(x) isnumeric(x);

averagesToRecon = 1:rInfo.nAverages;

repetitionsToRecon = 1:rInfo.nRepetitions;
repetitionsToReconValidationFcn = @(x) isnumeric(x);

p = inputParser();
addOptional(p,'Rbeta',Rbeta,RBetaValidationFcn);
addOptional(p,'L',L,LValidationFcn);
addOptional(p,'Niter', Niter, NIterValidationFcn);
addOptional(p,'slicesToRecon', slicesToRecon, slicesToReconValidationFcn);
addOptional(p,'averagesToRecon', averagesToRecon, repetitionsToReconValidationFcn);
addOptional(p,'repetitionsToRecon', repetitionsToRecon, repetitionsToReconValidationFcn);
parse(p,varargin{:});

% Override default parameters with their successfully parsed values.
Rbeta = p.Results.Rbeta;
L = p.Results.L;
Niter = p.Results.Niter;
slicesToRecon = p.Results.slicesToRecon;
repetitionsToRecon = p.Results.repetitionsToRecon;
averagesToRecon = p.Results.averagesToRecon;

[rInfoCC, senMBCC] = compressCoils(rInfo,senMB,'energyLevel',0.9);

% perform recon
img = phaseCorrectedRecon(rInfoCC, senMBCC, maskMB, FMMB, -1*PMaps,'Niter',Niter,'L',L,'Rbeta',Rbeta,'dims2penalize',[1,1,0],'averagesToRecon',averagesToRecon,'repetitionsToRecon',repetitionsToRecon,'slicesToRecon',slicesToRecon);
%img = fieldCorrectedRecon(rInfo, senMB, maskMB, FMMB,'Niter',Niter,'L',L,'Rbeta',Rbeta,'dims2penalize',[1,1,0],'repetitionsToRecon',repetitionsToRecon,'slicesToRecon',slicesToRecon);
% write data out to ISMRMRD format for PowerGrid
if isempty(filename)
    filename = 'img';
end
save(sprintf('%s.mat',filename),'img')

cd(curdir)
  