classdef TimeSegmentation
    %TimeSegmentation  - Convenience Class that implements field correction
    % for non-Cartesian MRI. Successor object to fast_mr_v2 and variants.
    % Works similarly to the TimeSegmentation class in PowerGrid.
    %
    % Syntax:  A = TimeSegmentation(G, timingVec, fieldMap, nSegs,varargin)
    %
    % Inputs:
    %    A         - An object capable of performing Fourier Transforms
    %                using the interface defined by NUFFT.m
    %    timingVec - A real vector matching the data size (or single shot
    %                length) that details the timing of the relaxation/off resonant signal decay
    %                relative to RF excitation isodelay. It is reccomended
    %                to give the timing vector of a single shot only as
    %                this is faster to compute the temporal interpolator
    %                and replicate than calculate a giant temporal
    %                interpolator.
    %    fieldMap  - A 2D or 3D array matching the image support suppling
    %                the off resonance in radians/s.
    %    nSegs     - Scalar giving number of time segments. For no field
    %                correction, supply a 0 NOT 1.
    %
    % Optional  Name-Value Pair Inputs:
    %    'mask'  - Image domain support mask of logicals to restrict the
    %              transform only to be computed at the true points.
    %
    %    'VoxelBasis' - 'delta'  (Default) Compute the transform assuming
    %                            each voxel is represented in image space
    %                            by a delta function.
    %
    %                   'boxcar' Compute the transform by assuming each
    %                            voxel is represented in image space by a
    %                            boxcar function, requiring a sinc
    %                            weighting to be applied in data/kspace
    %                            domain
    %    'Interpolator' - 'minmax' (Default) Computes the temporal
    %                               interpolator used to correct the off
    %                               resonance using a minmax scheme. The
    %                               computation of this interpolator can be
    %                               signficant.
    %
    %                     'histo'   Computes the minmax temporal interpolator
    %                               using a histogram based approximation.
    %                               Usually is much faster to compute.
    %                    
    %                     'hanning' This option uses a hanning window to
    %                               calculate the temporal interpolator for
    %                               each time segment. While the hanning
    %                               window is quite cheap to compute, the
    %                               accuracy of the hanning temporal
    %                               interpolator is lower than the minmax
    %                               temporal interpolator with the same
    %                               number of time segments. A rule of
    %                               thumb is to double the number of time
    %                               segments used with hanning than for
    %                               minmax .
    %
    % Outputs:
    %    A     - Time segmenation object with overloaded mtimes operator for
    %            computing the field corrected forward and adjoint NUFFT
    %            operations through A*x and A'*x.
    %
    % Example:
    %    % Create a spin-warp kspace
    %    [kx,ky] = ndgrid(-63:64,-63:64);
    %    % Initialize the object for a 2D 128x128 matrix recon.
    %    G =
    %    NUFFT(kx(:),ky(:),zeros(size(kx(:))),128,128,1);
    %
    %    % Compute a data domain estimate
    %    imageModel = phantom(128,128);
    %    dataEstimate = G*imageModel(:);
    %    imageEstimate = G'*dataEstimate;
    %    im(reshape(imageEstimate,128,128));
    %
    %
    % Other m-files required: nufft implementation and private directories
    % Subfunctions: none
    % MAT-files required: none
    %
    % Author: Alex Cerjanic
    % University of Illinois at Urbana-Champaign
    % email address:
    % Website:
    % 4-Sep-2017; Last revision: 4-Sep-2017
    properties
        G = [];  % Transform object, needs to correspond to the interface defined by NUFFT.m object
        fieldMap;     % Field map in radians
        timingVec;     % Timing Vector
        nShots; % Number of shots if a timing vector of only one shot is given, otherwise assume single shot trajectory!
        L; % Number of time segments
        mask;
        
        timeInterp; % Temporal Interpolators array
        
        Interpolator; % Name of the interpolator type
        isMasked = false;
        isEmpty = true;
        isTranspose = false;
        
    end
    
    methods
        function obj = TimeSegmentation(Gobj, tt, we, L, varargin)
            %% Parse inputs including name-value pairs
            obj.Interpolator = 'minmax'; % default
            p = inputParser();
            
            addOptional(p,'Interpolator',obj.Interpolator);
            obj.G = Gobj;
            parse(p,varargin{:});
            obj.timingVec = tt;
            obj.L = L;
            obj.fieldMap = we;
            
            %Use results from inputParser object.
            obj.Interpolator = p.Results.Interpolator;

            if sum(col(we)) == 0
                disp('L overriden to be 0 due to zero field map');
                obj.L = 0;
            end
            %% Prepare temporal interpolator
            % Computing the temporal interpolator can be quite
            % computationally costly. Since most multishot trajectories in
            % MR are made up of identical length shots, we can compute the
            % temporal interpolator for the timing vector of a single shot.
            % This approach is computationally easier, and we can simply
            % repmat back to the full vector.
            obj.nShots = length(obj.G.Kx(:))/length(obj.timingVec(:));
            
            switch(obj.Interpolator)
                case 'minmax'
                    %                 if (sum(abs(we_histo(:,1))) == 0)
                    %                     sprintf('Setting L=0')
                    %                     obj.L=0;
                    %                 end
                    obj.timeInterp = int_tim_seg(obj.timingVec,obj.L,obj.fieldMap(:),1);
                case 'histo'
                    % Calculate histogram of field map to use approximate
                    % minmax temporal interpolator
                    [bin_vals, bin_cens] = histcounts(obj.fieldMap(:),256);
                    % Have to deal with the fact that histcounts does not
                    % work like histc. (histc is being deprecated)
                    bin_cens = 0.5*bin_cens(1:end-1) + 0.5*bin_cens(2:end);
                    we_histo = [col(bin_cens), col(bin_vals)];
                    obj.timeInterp = int_tim_seg(obj.timingVec,obj.L,obj.fieldMap(:),2,we_histo);
                case 'synthetic_histo'
                    % Since actual histogram doesn't matter, just plug in a
                    % synthetic flat histogram with a range set to the
                    % maximum dynamic range of the field map. (Dictated by
                    % echo time difference.)
                    bin_centers = 2*pi*col(linspace(-1/(.5E-3),1/(.5E-3),256));
                    bin_counts = ones(size(bin_centers));
                    we_histo = [bin_centers, bin_counts];
                    % Have to deal with the fact that histcounts does not
                    % work like histc. (histc is being deprecated)
                    %bin_cens = 0.5*bin_cens(1:end-1) + 0.5*bin_cens(2:end);
                    %we_histo = [col(bin_cens), col(bin_vals)];
                    obj.timeInterp = int_tim_seg(obj.timingVec,obj.L,obj.fieldMap(:),2,we_histo);
                case 'hanning'
                    %ob.int = int_tim_seg(tt,L,we(:),int_opt);
                    %obj.timeInterp = int_tim_seghanning(obj.timingVec,obj.L);
                    
                    obj.timeInterp = int_tim_seghanning(obj.timingVec,obj.L);
                otherwise
                        error('TimeSegmentation: Unrecognized temporal interpolator type!');
             end
            
            
            % Mark object as constructed and instantiated
            % Not sure if this is needed with the new class structure.
            % -AMC 2017/08/24
            obj.isEmpty = false;
            
        end
        
        function obj = ctranspose(obj)
            obj.isTranspose = ~obj.isTranspose;
        end
        
        function y = mtimes(obj,x)
            
            if obj.isEmpty
                error('Object is empty and uninitialized!');
            end
            
            % Calculate time span of timing vector
            if (obj.L == 0)
                tau = 0;
            else
                tau = (max(obj.timingVec) - min(obj.timingVec) + eps)/(obj.L);
            end
            
            minTime = min(obj.timingVec);
            sizeG = size(obj.G);
            if ~obj.isTranspose
                % Forward operator
                x = exp(-1j*obj.fieldMap(:)*minTime).*x(:);
                y_temp = zeros(sizeG(1),obj.L+1);
                for ii = 1:(obj.L+1)
                    Wo = exp(-1j*obj.fieldMap(:)*((ii-1)*tau));
                    y_temp(:,ii) = obj.G*(Wo.*x(:));
                end
                aa = repmat(permute(obj.timeInterp, [2,1]), [obj.nShots, 1]);
                y = sum(y_temp.*aa, 2);
                
            else
                % Adjoint operator
                y_temp = zeros(sizeG(2),obj.L+1);
                
                for ii = 1:(obj.L+1)
                    Wo = exp(1j*conj(obj.fieldMap)*((ii-1)*tau));
                    aa = repmat(obj.timeInterp(ii,:)',[obj.nShots 1]);
                    y_temp(:,ii) = Wo(:).*(obj.G'*(aa.*x(:)));
                end
                
                y = sum(y_temp, 2);
                y = exp(1j*conj(obj.fieldMap(:))*minTime).*y;
                
            end
            
        end
        
        function dim = size(obj)
            dim = size(obj.G);
        end
        
        
    end
    
end

