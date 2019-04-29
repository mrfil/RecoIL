function [ mapOut ] = resampleMapNavUp(mapIn,rInfoOut,rInfoIn)
   %resampleMap Resample maps used in reconstructions using recoInfo.
   %   Detailed explanation goes here
   
   % Determine support of input data in (X,Y,Z) coordinates in mm.
   
   % Determine if data is complex
   isComplex = ~isreal(mapIn);
   
   % In-plane coordinates in and out
   %inputVec  = (((-(rInfoIn.N-1)/2:rInfoIn.N/2)) ) *  rInfoIn.FOV * 10 ./ rInfoIn.N;
   inputVec  = (((0:(rInfoIn.NNav-1)) - rInfoIn.NNav/2) + 0.5) *  rInfoIn.FOV * 10 ./ rInfoIn.NNav;
   outputVec = (((0:(rInfoOut.N-1)) - rInfoOut.N/2) + 0.5) * rInfoOut.FOV * 10 ./ rInfoOut.N;
   
   % Through-plane coordinates in and out
   if rInfoIn.multibandFactor > 1 % 2D SMS Case
      NzIn = rInfoIn.nSlices*rInfoIn.multibandFactor;
      thickness = rInfoIn.sliceThickness;
   elseif rInfoIn.nPartitions > 1% 3D Case
      NzIn = rInfoIn.nSlices*rInfoIn.nPartitionsNav;
      thickness = rInfoIn.sliceThickness/rInfoIn.nPartitionsNav;
   else % 2D Case
      NzIn = rInfoIn.nSlices;
      thickness = rInfoIn.sliceThickness;
   end
   FovZIn = NzIn*thickness;
   
   if rInfoOut.multibandFactor > 1 % 2D SMS Case
      NzOut = rInfoOut.nSlices*rInfoOut.multibandFactor;
      thickness = rInfoOut.sliceThickness;
   elseif rInfoOut.nPartitions > 1 % 3D Case
      NzOut = rInfoOut.nSlices*rInfoOut.nPartitions;
      thickness = rInfoOut.sliceThickness/rInfoOut.nPartitions;
   else % 2D case
      NzOut = rInfoOut.nSlices;
      thickness = rInfoOut.sliceThickness;
   end
   FovZOut = NzOut*thickness;
   
   
   inputVecSlice  = ((( 0:(NzIn-1))  - NzIn /2) + 0.5) * FovZIn   ./NzIn;
   outputVecSlice = ((( 0:(NzOut-1)) - NzOut/2) + 0.5) * FovZOut  ./NzOut;
   
   % Generate the coordinates for the data
   [xInput, yInput, zInput] =  meshgrid(inputVec,inputVec,inputVecSlice);
   [xOutput, yOutput, zOutput] = meshgrid(outputVec,outputVec,outputVecSlice);
   
   if isComplex
      % Perform Interpolation 
      realMapOut = interp3(xInput,yInput,zInput,real(mapIn), xOutput, yOutput, zOutput, 'spline');
      imagMapOut = interp3(xInput,yInput,zInput,imag(mapIn), xOutput, yOutput, zOutput, 'spline');
      mapOut = realMapOut + 1j*imagMapOut;
   else
      mapOut = interp3(xInput, yInput, zInput, mapIn, xOutput, yOutput, zOutput, 'spline');
   end
   
end

