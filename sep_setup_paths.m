function sep_setup_paths()
% SEP_SETUP_PATHS Put the project's code on the MATLAB path.
%
% For now this adds only the code root, the folder holding this file. That is
% all today's drivers expect: each of them finds the vendored toolboxes
% (LightSuite, yamlmatlab, matlab_elastix, Bio-Formats) through get_paths and
% adds them itself. The verification tools in tools\ are added by hand when a
% check needs them (see tools\README.md).
%
% Step 4 of docs/REFACTOR_PLAN.md turns this into the explicit list of the
% code folders. It stays a list rather than genpath, so archive\, tests\, the
% Python environments in tools\ and the .git and worktree folders never reach
% the path. No data or atlas folder is added here: those come from get_paths.
%
% A check starts from a fresh session, so nothing an earlier run put on the
% path can shadow the code under test:
%
%   restoredefaultpath; cd('D:\sep_histology\code'); sep_setup_paths

code_dir = fileparts(mfilename('fullpath'));
addpath(code_dir);

end
