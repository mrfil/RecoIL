function mask = masking(image, radius)

% masking .m 
% This m-file will be used to create the masking around the object
%   to set the intensity equal to zero at regions outside the object
% Makes a Circular region

[m n] = size(image);

x = (0:n-1)-n/2+1/2;
y = (0:m-1)-m/2+1/2;

mask = zeros(m,n);

[X Y] = meshgrid(x,y);

mask = (sqrt(X.^2+Y.^2)<radius);

