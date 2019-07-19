function rInfo = recoInfo(filename, varargin)
    
    if nargin < 1 || isempty(filename)
        fileList = [dir('*.dat'), dir('*.h5')];
        if length(fileList) ~= 1
            error('Must specify filename if more than one recoInfo compatible file exists in the current directory.');
        else
            filenameSelected = fileList.name;
        end
    else
        filenameSelected = filename;
    end
    
    [filepath,name,ext] = fileparts(filenameSelected);
    switch ext
        case '.dat'
            rInfo = recoInfoTwix(filenameSelected, varargin{:});
        case '.h5'
            rInfo = recoInfoIsmrmrd(filenameSelected, varargin{:});
        otherwise
            error('Unrecogized filetype. recoInfo currently supports Siemens Twix (.dat) and ISMRMRD HDF5 (.h5) files.')
    end
end
