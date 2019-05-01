function csm = generate_birdcage_sensitivities(matrixSize, nCoils, relativeRadius)
    
    csm = zeros(matrixSize, matrixSize, nCoils);
    
    for ii = 1:nCoils % coil index
        
        coilX = relativeRadius*cos((ii-1)*(2*pi/nCoils));
        coilY = relativeRadius*sin((ii-1)*(2*pi/nCoils));
        
        coilPhase = -(ii-1)*(2*pi/nCoils);
        
        for kk = 1:matrixSize % x
            for jj = 1:matrixSize % y
                yCo = (jj - matrixSize/2)/(matrixSize/2) - coilY;
                xCo = (kk - matrixSize/2)/(matrixSize/2) - coilX;
                rr = sqrt(xCo.^2 + yCo.^2);
                phi = atan2(xCo,-yCo) + coilPhase;
                csm(kk,jj,ii) = 1./rr.*exp(1j*phi);
            end
        end
    end
    csm(isinf(csm)) = 1;
    
end