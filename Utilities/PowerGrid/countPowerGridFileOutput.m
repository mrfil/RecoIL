function [NSlices,NReps,NAvgs,NEchoes,NPhases] = countPowerGridFileOutput(directory)
    %countPowerGridFileOutput Count extents of PowerGrid File output for
    %further processing.

    if nargin < 1
        directory = pwd;
    end
    
    files = dir([directory '/' '*_mag.nii']);
    NFiles = length(files);
    
    filesSorted = natsortfiles({files.name});
    

    extents = sscanf(filesSorted{end}, 'img_Slice%u_Rep%u_Avg%u_Echo%u_Phase%u_mag.nii');
    
    % Convert from zero indexed (Unix/C++ style) to one indexed (MATLAB
    % style)
    NSlices = extents(1) + 1;
    NReps = extents(2) + 1;
    NAvgs = extents(3) + 1;
    NEchoes = extents(4) + 1;
    NPhases = extents(5) + 1;
    
end

