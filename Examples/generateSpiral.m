function [kRead, kPhase, gRead, gPhase] = generateSpiral(FOV, N, GAmp, GSlew, tSamp, nShots)
    %GENSPIRAL.M Convenience function to use Brian Hargreaves vds function
    %with our particular preferred units and formats.
    
    rmax = 1/(2*FOV/N);
    
    [k,g,s,time,r,theta] = vds(GSlew*100,GAmp/10,tSamp,nShots,FOV,rmax);
    
    % Convert from k-space in units of reciprocal length (1/cm) to unitless
    % k-space for the single shot.
    
    kx = real(k)*FOV;
    ky = imag(k)*FOV;
    
    gx = real(g);
    gy = imag(g);
    
    % Now we need to rotate for all of the shots.
    
    % Handle Rotations to generate complete k-space sampling pattern
    phi = 2*pi/nShots;
    for ii = 0:(nShots-1)
        ang_rot = phi * (ii-(nShots-1) * floor(ii/nShots));
        kRead(:,ii+1) = kx * cos(ang_rot) + ky * sin(ang_rot);
        kPhase(:,ii+1) = ky * cos(ang_rot) - kx * sin(ang_rot);
    end

    phi = 2*pi/nShots;
    for ii = 0:(nShots-1)
        ang_rot = phi * (ii-(nShots-1) * floor(ii/nShots));
        gRead(:,ii+1) = gx * cos(ang_rot) + gy * sin(ang_rot);
        gPhase(:,ii+1) = gy * cos(ang_rot) - gx * sin(ang_rot);
    end
end

