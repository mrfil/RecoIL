classdef Wavelet2
    %WAVELET Simple 2D/3D Convenience class for implementing orthogonal
    %discrete wavelet transform. Inspired by Miki Lustig's @Wavelet class
    %from the sparseMRI package.
    %   Detailed explanation goes here
    
    properties
        isTrans = false;
        Nx;
        Ny;
        mask;
        S;
        Level;
    end
    
    methods
        function obj = Wavelet2(Nx,Ny,mask,Level)
            %WAVELET Construct an instance of this class
            %   Detailed explanation goes here
            
            obj.Nx = Nx;
            obj.Ny = Ny;
            if nargin > 2
                obj.mask = logical(mask);
            else
                obj.mask = true(obj.Nx,obj.Ny);
            end
            if nargin > 3
                obj.Level = Level;
            else
                obj.Level = 4;
            end
            
            % Dummy wavelet transform to get S matrix for studpid matlab;
            [~,obj.S] = wavedec2(zeros(Nx,Ny),obj.Level,'haar');
            
        end
        
        function y = mtimes(obj,x)
            %mtimes Implements mtimes y = A*x operation
            %   Detailed explanation goes here
            if obj.isTrans % Adjoint Transformation
                img = waverec2(real(x),obj.S,'haar') + 1j*waverec2(imag(x),obj.S,'haar');
                y = col(img(obj.mask));
            else % Forward Transformation
                img = embed(x,logical(obj.mask));
                wd =  wavedec2(real(img),obj.Level,'haar') + 1j*wavedec2(imag(img),obj.Level,'haar');
                y = col(wd);
            end
            
        end
        
        function obj = ctranspose(obj)
            obj.isTrans = ~obj.isTrans;
        end
        
        % DO we need a size function? - AMC 11/11/2017
    end
end

