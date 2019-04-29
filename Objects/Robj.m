classdef Robj
% Build roughness penalty regularization "object" based on C = Cdiff() object,
% for regularized solutions to inverse problems.
%
% General form of nonquadratic penalty function:
%	R(x) = \sumk w_k \pot([Cx]_k), where [Cx]_k = \sum_j c_{kj} x_j.
%	
% For quadratic case, \potk(t) = t^2/2 in which case
%	R(x) = x' C' W C x / 2, where W depends on beta and edge_type.
%
% Penalty gradient is C' D C x, where D = diag{\wpot_k([Cx]_k)}.
%
% in
%	kappa	[nx,ny[,nz]]	kappa array, or logical support mask
%
% options
%	'offsets', [?]		offsets to neighboring pixels
%					(see Cdiff for the defaults)
%				use '3d:26' to penalize all 13 pairs of nbrs
%				use '0' for C = I (identity matrix)
%	'beta', ?		global regularization parameter
%					default: 2^0
%	'delta', ?		potential parameter, see potential_func()
%					or {delta, param}.  default: inf
%	'potential', '?'	e.g., 'quad' or 'approxTV'
%				default: 'quad' for quadratic regularization.
%	'user_wt', [?]		User-provided array of penalty weight values
%					size: [nx,ny[,nz],#offsets]
%				of dimension [size(mask) length(offsets)].
%				These are .* the usual wt values for edge_type.
%	'mask'			Override default: mask = (kappa ~= 0)
%
% out
%	R structure has the following inline "methods" 
%	R.penal(R, x)	evaluates R(x)
%	R.cgrad(R, x)	evaluates \cgrad R(x) (column gradient)
%	R.denom(R, x)	evaluates denominator for separable surrogate
%	[pderiv pcurv] = feval(R.dercurv, R, C1*x) derivatives and curvatures
%			for non-separable parabola surrogates
%	R.diag(R)	diagonal of Hessian of R (at x=0), for preconditioners.
%	R.C1		differencing matrix, with entries 1 and -1,
%			almost always should be used in conjunction with R.wt
%
% Typical use:	mask = true([128 128]); % or something more conformal
%		R = Robject(mask, 'beta', 2^7);
%
% Copyright 2004-11-14, Jeff Fessler, The University of Michigan
% Copyright 2019-02-18, Alex Cerjanic, University of Illinois at
% Urbana-Champaign
% Modified to use modern MATLAB class object system for compatibility

    properties
        potential;
        beta;
        delta;
        mask;
        offsets = [];
        dims2penalize = [];
        C1 = [];
        potk;
        wpot;
        dpot;
        wt;
    end
    
    methods
        function obj = Robj(mask, varargin)
            %Robj Construct an instance of this class
            %   Detailed explanation goes here
            
            if ~islogical(mask)
                error 'Mask must be logical'
            end
            
            %% Parse inputs
            p = inputParser();
            p.addOptional('beta', 1, @(x) (x>0));
            p.addOptional('dims2penalize', ones(ndims(mask),1));
            p.addOptional('potential','quad');
            % Deal with inputs here
            p.parse();
            
            obj.beta = p.Results.beta;
            obj.dims2penalize = p.Results.dims2penalize;
            obj.potential = p.Results.potential;
            
            %% Initialize object
            % Set offsets
            if ndims(mask) == 2
                [nx, ~] = size(mask);
                obj.offsets = [1, nx, nx+1, nx-1];
            elseif ndims(mask) == 3
                [nx, ny] = size(mask);
                obj.offsets = [1, nx, nx*ny];
            else
                error 'only 2D and 3D support'
            end
            obj.mask = mask;
            % Not supporting the 26 neigbor 3d regularization due to lack
            % of use from original Robj.m
            
            % Make this play nice with the penaly mex files
            obj.offsets = int32(obj.offsets);
            
            %% Create sparse differencing matrix
            [obj.C1, obj.wt] = C_sparse(obj.mask,obj.dims2penalize);
            obj.wt = obj.beta * obj.wt(:);
            
            
            %% Setup the potential functions
            switch obj.potential
                case 'quad'
                    obj.potk = @(t) (abs(t).^2)./2;
                    obj.wpot = @(t) ones(size(t));
                    obj.dpot = @(t) t;
                case 'approxTV'
                    obj.potk = @(t) (obj.delta.^2).*(sqrt(1.0 + abs(d./obj.delta).^2) - 1.0);
                    obj.wpot = @(t) 1./sqrt(1.0 + abs(d./obj.delta));
                    obj.dpot = @(t) t./sqrt(1.0 + abs(d./obj.delta));
            end
        end
        
        function out = cgrad(obj, x)
            out = obj.C1' * (diag_sp(obj.wt) * obj.dpot(obj.C1 * x));
        end
        
        function out = penal(obj, x)
            out = sum(repmat(obj.wt, ncol(x)) .* obj.potk(obj.C1 * x));
        end
        
        function out = denom(obj, ddir, x)
            Cdir = (obj.C1 * ddir) .* (obj.wt > 0);
            Cx = (obj.C1 * x) .* (obj.wt > 0);
            out = Cdir' * (obj.wpot(Cx) .* Cdir);
        end
    end
end

