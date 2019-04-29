function ob = sense_svd(A,sen,varargin)
%function ob = sense_svd(A,sen, varargin)
%	Construct MRI object, which can do Ax and A'y operations using SENSE.
%
%  Allowed Calling syntax
%
%  S_svd = sense_svd(A,sen); Uses default method to determine coil rank
%  S_svd = sense_svd(A,sen,coil_rank); specifies coil_rank as the number of
%                       virtual coils to include
%  S_svd = sense_svd(A,sen,name value pairs); Uses the MATLAB name-value
%                                   pair syntax to take optional arguments
%
%  Calling arguments:
%
%  A         - underlying transform object (such as NUFFT, DFT, or fast_mr)
%  sen       - coil sensitivity matrix (SENSE map) [npts, ncoils]
%  coil_rank - (optional) number of virtual coils to use. (0 < coil_rank <= ncoils)
%
% Optional Name-Value Pairs:
% EnergyLevel - Specifies minimum level of energy on (0,1] measured by absolute
%               value of the singluar values of the coil sensitivity matrix
%               (SENSE map) for the coil rank selection algorithm.

%	default object
ob.A = 0;	% should these be []'s
%ob.sen = 0;
ob.rank_svd = 0;
ob.num_coils = 0;
ob.V = 0;    % This is the transform from original coils to combined
ob.VS = 0;   % This is the combined sensitivity maps
ob.is.empty	= logical(1);
ob.is.transpose = logical(0);
%ob.version = 1.0;

% Deal with input arguments using inputParser
p = inputParser;
p.CaseSensitive = true;
addRequired(p,'A');
senValidationFun = @(x) (isnumeric(x) && ismatrix(x));
addRequired(p,'sen', senValidationFun);
sizeSen = size(sen);
defaultCoilRank = 0;
coilRankValidationFun = @(x) (isnumeric(x) && isscalar(x) && (x > 0) && (x <= sizeSen(2)));
addOptional(p,'CoilRank',defaultCoilRank,coilRankValidationFun);

energyLevelValidationFun = @(x) (isnumeric(x) && isscalar(x) && (x > 0) && (x < 1));
addParameter(p,'EnergyLevel',0.9,energyLevelValidationFun);

parse(p,A,sen,varargin{:});
coil_rank = p.Results.CoilRank;
energyLevel = p.Results.EnergyLevel;

if nargin == 0
    ob = class(ob, 'sense_svd');
    return
end

if isa(A, 'sense_svd')
    ob = A;
    return
end

% if nargin ~= 3  %7
% 	help sense_svd
% 	error nargin
% end

% if size(sen,3) ~= 1
%    sprintf('Sense maps must be num pts x num coils. ')
% end


%	fill object
ob.A = A;
%ob.sen = sen;
[npts num_coils] = size(sen);
ob.num_coils = num_coils;
[~,sig,V] = svd(sen,0);

if (coil_rank == 0)
    sig = diag(sig);
    rank_svd = chooseCoilRank(sig,energyLevel);
    disp(['Sense_svd coil rank selection returned ' num2str(rank_svd) ' virtual coils based on ' num2str(energyLevel) ' energy level.']);
else
    rank_svd = coil_rank;
    disp(['Sense_svd coil rank selection overridden to be ' num2str(rank_svd) ' virtual coils.']);
end


ob.V = V;

ob.rank_svd = rank_svd;
VS = sen*V(:,1:rank_svd);

ob.VS = VS;
%ob.sen = []; % Do not need complete sense map after this point.
ob.is.empty	= logical(0);

%	ob.m = size(we);	% image size
%	ob.n = size(we);

ob = class(ob, 'sense_svd');
end


