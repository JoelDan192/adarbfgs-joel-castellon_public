%% setup_InvRand Setup file 
%
%    Adds to path the files of the current folder
%    And compiles the mex file libsvmread for use with the ridge regression problems
%    problems

% Get root location of resultsSetup
root = fileparts(which(mfilename)); 

% ----------------------------------------------------------------------
% Add the appropriate subdirs to path.
% ----------------------------------------------------------------------

addpath(genpath(root));
fprintf('%s \n',root);

% ----------------------------------------------------------------------
% Add the appropriate paths for libsvm and UF matrices
% ----------------------------------------------------------------------
LIBSVMpath = './tests/LIBSVM'; %% The local LIBSVM folder,  contains  libsvmread.c, libsvmread.mex and data files (e.g. liver-disorders)
% 
UFgetpath = './tests/UFget';  %% The local UFget folder, that should contain  UFget.m


try
   if exist('libsvmread','file')~=3, mex libsvmread.c; end
   fprintf('Compilation of MEX files for libsvmread was successful!\n');
catch err
    display(err);
    cd(root);
   error('Could not compile MEX files for libsvmread:');
end
fprintf(['\n The above directories and their subdirectories have been'  ,...
         ' temporarily added to your path.\n']);
     
clear;
