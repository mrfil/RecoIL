function [FM, FMdx, FMdy] = genSynthFieldMap(N,offRes)
   %genSynthFieldMap Generate a synthetic field map for use with simulation and
   %test cases
   %   Inputs:
   %     N  - Matrix Size (Defaults to 256)
   %     offres - Peak Off-resonance in Hz.
   %     
   %  Outputs: 
   %     FM - Field Map with 3 bumps
   
   [im1, im1dx, im1dy] = bump(N,floor(N/20),floor(N/4),floor(N/2),2,5);
   im1 = offRes*im1;
   im1dx = offRes*im1dx;
   im1dy = offRes*im1dy;
   
   [im2, im2dx, im2dy] = bump(N,floor(N/40),floor(5*N/9),floor(N/4),1);
   im2 = offRes/2*im2;
   im2dx = offRes/2*im2dx;
   im2dy = offRes/2*im2dy;
   
   [im3, im3dx, im3dy] = bump(N,floor(N/40),floor(5*N/9),floor(3*N/4),1);
   im3 = offRes/2*im3;
   im3dx = offRes/2*im3dx;
   im3dy = offRes/2*im3dy;
   
   FM = im1+im2+im3;
   FMdx = im1dx+im2dx+im3dx;
   FMdy = im1dy+im2dy+im3dy;
   
end

