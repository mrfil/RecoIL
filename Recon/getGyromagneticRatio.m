function gamma = getGyromagneticRatio(nuclei)
%getGyromagneticRatio Convenience function to get gyromagnetic ratio for
% X-nuclei applications. For simplicity and sanity, it is reccomended that 
% you use this everywhere you need gamma.
% 
%   Inputs
%   nuclei -    MR Active nuclei in form of:
%                '1H'    : Proton
%                '23Na'  : Sodium-23
%                '31P'   : Phosphorous-31
%
%   Outputs
%   gamma -    Gyromagnetic Ratio in units of rad/s/T.

% Deal with no input case
if nargin < 1
    nuclei = '1H'
end

switch nuclei
    case '1H'
        gamma = 42.57747892*2*pi*1E6;
    case '23Na'
        gamma = 70.61*1E6;
    case '31P'
        gamma = 108.291*1E6;
    otherwise
        error('Unrecognized nuclei provided.');
end

end

