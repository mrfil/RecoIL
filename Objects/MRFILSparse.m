classdef MRFILSparse
   % MRFILSparse Summary of this class goes here
   %   Detailed explanation goes here
   
   properties
      rowIndices = [];
      colIndices = [];
      values = [];
      N;
      M;
      isTranspose = false;
   end
   
   methods
      function obj = MRFILSparse(N,M,rowInd, colInd, values)
         obj.N = uint32(N);
         obj.M = uint32(M);
         obj.rowIndices = uint32(rowInd(:));
         obj.colIndices = uint32(colInd(:));
         obj.values = values(:);
         
      end
      
      
      function obj = ctranspose(obj)
         obj.isTranspose = ~obj.isTranspose;
      end
      
      function sizes = size(obj)
         sizes = [obj.N, obj.M];
      end
      
      function y = mtimes(obj,x)
         
         if ~obj.isTranspose % Non-transpose multiplication
            b = x(obj.colIndices);
            y = accumarray(obj.rowIndices,obj.values.*b);
         else % Transpose multiplication
            b = x(obj.rowIndices);
            y = accumarray(obj.colIndices,conj(obj.values).*b);
         end
         
      end
   end
end


