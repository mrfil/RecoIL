function reconMultibandMRE(datadir,filename,varargin)
% function reconMultibandMRE(datadir,filename,varargin)
% 
% Inputs (optional):
%
% - datadir:  string pointing to data directory
%             (if not entered, will prep current directory)
%
% - filename: string with filename for reconned image file
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

rInfo = recoInfo;

if ~(exist('sen.mat','file')&&exist('FM.mat','file')&&exist('mask.mat','file')&&exist('PMaps.mat','file'))
    error('Missing prep files; run prepMultibandMRE first')
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

repetitionsToRecon = 1:rInfo.nRepetitions;
repetitionsToReconValidationFcn = @(x) isnumeric(x);

p = inputParser();
addOptional(p,'Rbeta',Rbeta,RBetaValidationFcn);
addOptional(p,'L',L,LValidationFcn);
addOptional(p,'Niter', Niter, NIterValidationFcn);
addOptional(p,'slicesToRecon', slicesToRecon, slicesToReconValidationFcn);
addOptional(p,'repetitionsToRecon', repetitionsToRecon, repetitionsToReconValidationFcn);
parse(p,varargin{:});

% Override default parameters with their successfully parsed values.
Rbeta = p.Results.Rbeta;
L = p.Results.L;
Niter = p.Results.Niter;
slicesToRecon = p.Results.slicesToRecon;
repetitionsToRecon = p.Results.repetitionsToRecon;

% perform recon
[rInfoCC,senCC] = compressCoils(rInfo,senMB,'energyLevel',0.90);
% [rInfoCC,senCC] = compressCoils(rInfo,senMB,'coilRank',5);

%if rInfo.nCoils < 32
%img = phaseCorrectedRecon(rInfo, senMB, maskMB, FMMB, PMaps,'Niter',Niter,'L',L,'Rbeta',Rbeta,'dims2penalize',[1,1,0],'repetitionsToRecon',repetitionsToRecon,'slicesToRecon',slicesToRecon);
%else
img = phaseCorrectedRecon(rInfoCC, senCC, maskMB, FMMB, PMaps,'Niter',Niter,'L',L,'Rbeta',Rbeta,'dims2penalize',[1,1,0],'repetitionsToRecon',repetitionsToRecon,'slicesToRecon',slicesToRecon);
%end

% write data out to ISMRMRD format for PowerGrid
if isempty(filename)
    filename = 'img';
end
save(sprintf('%s.mat',filename),'img')

cd(curdir)
  