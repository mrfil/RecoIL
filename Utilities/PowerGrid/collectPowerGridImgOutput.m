function [img] = collectPowerGridImgOutput(directory)
    %collectPowerGridImgOutput Convenience function to collect PowerGrid
    %output for use in matlab.
    %   Input:
    %       directory - path to find files. If blank or missing, use
    %                   current directory
    %
    %   Output:
    %       img      - Complex image.
    
    if nargin < 1
        directory = pwd;
    end
    
    [NSlices,NReps,NAvgs,NEchoes,NPhases] = countPowerGridFileOutput(directory);
    
    img = mergePowerGridFileOutput(NSlices,NReps,NAvgs,NEchoes,NPhases,directory);
end

