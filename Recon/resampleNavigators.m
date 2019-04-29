function [ mapOut ] = resampleNavigators(mapIn,rInfo)
   %resampleNavigators Resample navigators used in pcSENSE reconstructions 
   % using recoInfo.
   %   Detailed explanation goes here
   
   % Determine support of input data in (X,Y,Z) coordinates in mm.
   if((rInfo.multibandFactor == 1) && (rInfo.nPartitions > 1))
       error('3D navigator case is not supported in this function. Use resampleMapNav instead.');
   end
   
   
   % Determine if data is complex
   isComplex = ~isreal(mapIn);
   
   % In-plane coordinates in and out
   %inputVec  = (((-(rInfoIn.N-1)/2:rInfoIn.N/2)) ) *  rInfoIn.FOV * 10 ./ rInfoIn.N;
   inputVec  = (((0:(rInfo.NNav-1)) - rInfo.NNav/2) + 0.5) *  rInfo.FOV * 10 ./ rInfo.NNav;
   outputVec = (((0:(rInfo.N-1)) - rInfo.N/2) + 0.5) * rInfo.FOV * 10 ./ rInfo.N;
   
   % Through-plane coordinates in and out
   if rInfo.multibandFactor > 1 % 2D SMS Case
      NzIn = rInfo.multibandFactor;
      thickness = rInfo.sliceThickness;
   else % 2D Case
      NzIn = 1;
   end
   FovZIn = NzIn*thickness;
   
   if rInfo.multibandFactor > 1 % 2D SMS Case
      NzOut = rInfo.multibandFactor;
   else % 2D case
      NzOut = 1;
   end
   %FovZOut = NzOut*thickness;
   
   

   
   % Generate the coordinates for the data
   [xInput, yInput] =  meshgrid(inputVec,inputVec);
   [xOutput, yOutput] = meshgrid(outputVec,outputVec);
   mapOut = zeros(rInfo.N,rInfo.N,rInfo.nPartitions);
   for ii = 1:NzOut
       if isComplex
           % Perform Interpolation
           realMapOut = interp2(xInput,yInput,real(mapIn(:,:,ii)), xOutput, yOutput, 'spline');
           imagMapOut = interp2(xInput,yInput,imag(mapIn(:,:,ii)), xOutput, yOutput, 'spline');
           mapOut(:,:,ii) = realMapOut + 1j*imagMapOut;
       else
           mapOut(:,:,ii) = interp2(xInput, yInput, mapIn(:,:,ii), xOutput, yOutput, 'spline');
       end
   end
end

