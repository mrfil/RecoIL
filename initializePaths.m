% Intialize Paths to use Recon Software -- Alex Cerjanic 12-Sep-2017
%
% Instructions: To add necessary paths to use this software to your
% paths. Add the path containing this file to your path in your startup and
% then call this script. As directories are added to the recon repository,
% directories have to be added to this script.

% Get the root path of the repository
rootPath = fileparts(mfilename('fullpath'));

% Replacing the inclusion of the entire IRT with our minimal fork of the
% IRT
addpath([ rootPath '/irt/graph']);
addpath([ rootPath '/irt/mri']);
addpath([ rootPath '/irt/nufft']);
addpath([ rootPath '/irt/nufft/table']);
addpath([ rootPath '/irt/penalty']);
addpath([ rootPath '/irt/systems']);
addpath([ rootPath '/irt/utilities']);
addpath([ rootPath '/irt/mex/v7']);

% Adding lab software.
addpath([ rootPath '']);
addpath([ rootPath '/Recon']);
addpath([ rootPath '/Routines']);
addpath([ rootPath '/Gridding']);
addpath([ rootPath '/Objects']);



% Adding all of the external tools
% Note that we don't use genpath here to control what is added. Exactly what
% paths are added are controlled in initializeContribPaths.m in the /Contrib
% directory.
addpath([ rootPath '/Contrib']);
initializeContribPaths;

