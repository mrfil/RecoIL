function gammaBar = getGammaBar(nuclei)
%getGammaBar Get gyromanetic ratio for a given MR active nuclei in terms of
%            gamma bar
%   Input
%     nuclei     -    Nuclei name as recognized by getGyromagneticRatio()
%
%   Output
%     gammaBar   -    Gyromagnetic ratio for given nuclei in units of Hz/T. 
%                       If nuclei is an empty string, gamma bar for proton 
%                       will be returned. 

gammaBar = getGyromagneticRatio(nuclei)/(2*pi);

end

