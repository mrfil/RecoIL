% initializeContribPaths.m - Alex Cerjanic 2017/10/15 
% Adds externally developed tools to MATLAB paths. When you add tools to the
% external folder, add the required paths to this file to ensure they are
% available to all users.
   
% Get the root path of the repository
rootPathContrib = fileparts(mfilename('fullpath'));

% Add miscelaneous single files
addpath([ rootPathContrib '/Misc']);

