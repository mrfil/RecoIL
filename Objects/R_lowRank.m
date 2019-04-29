classdef R_lowRank
    
    properties
        R; % Robject that deals with a spatial domain image
        N; % Number of temporal images
        Nrank; % Number of ranks
        v; % Temporal basis for low rank
    end
    
    methods 
        function obj = R_lowRank(R, N, Nrank, v)
        
            obj.R = R; 
            obj.N = N;
            obj.v = v;
            obj.Nrank = Nrank;
        end
        
        function y = cgrad(obj, x)
            spatialBasis = reshape(x,[],obj.Nrank);
            imgStack = reshape(spatialBasis*obj.v(:,1:obj.Nrank)',[],obj.N);
            %x = reshape(x,[],obj.Nrank);
            
            for ii = 1:obj.N
               %y(:,ii) = obj.v(:,1:obj.Nrank) * obj.R.cgrad(obj.R,imgStack(:,ii));
               y(:,ii) =  obj.R.cgrad(obj.R,imgStack(:,ii));
            end
            
            y = col(obj.v(:,1:obj.Nrank).' * y.' );
        end
        
        function y = denom(obj, ddir, x)
           spatialBasis = reshape(x,[],obj.Nrank);
           imgStack = reshape(spatialBasis*obj.v(:,1:obj.Nrank)',[],obj.N);
           %x = reshape(x, [], obj.N);
           ddir = reshape(ddir,[],obj.Nrank);
           spatialDdir = reshape(ddir*obj.v(:,1:obj.Nrank)',[],obj.N);
           %ddir = reshape(ddir, [], obj.N);
           temp = zeros(obj.N,1);
           for ii = 1:obj.N
               
%                Cdir = (obj.R.C1 * ddir(:,ii)).*(obj.R.wt > 0);
%                Cx = obj.R.pot.wpot(obj.R.pot,(obj.R.C1 * x(:,ii)).*(obj.R.wt > 0));
%                temp(ii) = Cdir' * (Cx .* Cdir);
                %Rtemp = obj.R;  
                temp(ii) = obj.R.denom(obj.R,spatialDdir(:,ii),imgStack(:,ii));
              
           end
           y = obj.R.beta * sum(temp);
           
        end
        
        function y = penal(obj,x)
           spatialBasis = reshape(x,[],obj.Nrank);
           imgStack = reshape(spatialBasis*obj.v(:,1:obj.Nrank)',[],obj.N);
           for ii = 1:obj.N
               y(ii) = obj.R.penal(obj.R,imgStack(:,ii));
           end
           y = sum(y);
           
        end
        
        
    end
end
