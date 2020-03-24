function [ unwrappedPhase ] = prelude( rawPhase, mag, varargin)
%PRELUDE 
%
%   PRELUDE calls FSL Prelude
%
%	unwrappedPhase = PRELUDE( rawPhase, mag )
%	unwrappedPhase = PRELUDE( rawPhase, mag, ... )
%
%    .......................
%
%    rawPhase: 3d array of the wrapped phase values
%    mag:      3d array of the corresponding magnitude data
%
%    .......................
%   
%   The following name-value pairs are supported
%
%        voxelSize
%               default: [1 1 1]
%
%        .mask
%
%
%        isUnwrappingIn2D
%                        default: true
%
%    .......................
% 	topfer@ualberta.ca	2014
%   Modified by Alex Cerjanic, acerja2@illinois.edu University of Illinois, 2019
%		to use name-value pairs for options
%	

	DEFAULT_VOXELSIZE             = [1 1 1] ;
	DEFAULT_ISUNWRAPPINGIN2D      = false ;

	% check inputs

	if nargin < 2 || isempty(rawPhase) || isempty(mag) 
	    error('Function requires at least 2 input arguments.')
	end
	p = inputParser();
	p.addOptional('voxelSize', DEFAULT_VOXELSIZE);
	p.addOptional('isUnwrappingIn2D', false);
	p.addOptional('mask',[])
	p.parse(varargin{:});

	voxelSize = p.Results.voxelSize;
	isUnwrappingIn2D = p.Results.isUnwrappingIn2D;
	mask = p.Results.mask;
	% Create a temp folder for working in.
	tmpFldr = 'tmpPreludeFldr/' ;
	system(['mkdir ' tmpFldr]) ;
	
	save_nii( make_nii( rawPhase, voxelSize ), [tmpFldr 'rawPhase.nii'] ) ;

	save_nii( make_nii( mag, voxelSize ), [tmpFldr 'mag.nii'] ) ;
	
	optionsEtc = '';
	if isUnwrappingIn2D
		optionsEtc = ' -s ' ;
	end

	if ~isempty(mask)
		save_nii(make_nii(mask, voxelSize), [tmpFldr 'mask']) ;

		optionsEtc   = [ optionsEtc ' -m ' tmpFldr 'mask' ] ;
	end

	unwrapCommand = ['$FSLDIR/bin/prelude -p ' tmpFldr 'rawPhase' ...
					 		' -a ' tmpFldr 'mag' ...
					 		' -o ' tmpFldr  'unwrappedPhase.nii.gz' ...
					 		' ' optionsEtc] ; 

	system(unwrapCommand) ;

	unwrappedPhaseNII = load_nii([tmpFldr 'unwrappedPhase.nii.gz']);

	unwrappedPhase = double(unwrappedPhaseNII.img);

	delete( [tmpFldr 'rawPhase.nii'], [tmpFldr 'mag.nii'], ...
		[tmpFldr 'unwrappedPhase.nii'] ) ;
                
    if ~isempty(mask)
        delete( [tmpFldr 'mask.nii'] ) ;
    end
end
