function sep_setup_paths()
%SEP_SETUP_PATHS  Put the project's code on the MATLAB path.
%   SEP_SETUP_PATHS() runs once per session, before any driver. A check starts
%   from a fresh session, so nothing an earlier run put on the path can shadow
%   the code under test:
%
%     restoredefaultpath; cd('D:\sep_histology\code'); sep_setup_paths
%
%   It adds an explicit list of our code folders, and the vendored toolboxes
%   in third_party\ (LightSuite, yamlmatlab, matlab_elastix, Bio-Formats) with
%   their subfolders. The order on the path, from the top: our folders, then
%   the toolboxes, then MATLAB's own functions: LightSuite above MATLAB's
%   functions is the precedence every result so far was produced with.
%   tests\sep_test_path checks that no function name is defined twice on this
%   path.
%
%   Never added: the code root with its subfolders (genpath), archive\,
%   tests\, tools\ (a check adds it by hand, see tools\README.md), mapping\
%   (Python), the Python engines' folders, and any data or atlas folder.
%   get_atlas adds the one atlas folder a run needs: LightSuite finds the
%   atlas with which(), so a second atlas folder on the path would be picked
%   up silently.
%
%   Warns about any .m file at the code root other than get_paths.m and this
%   file: an editor tab left open on a file the refactor moved, and saved
%   again, recreates it there, and the old copy would then run in place of the
%   moved one (docs/refactor_name_map.csv says where each file went).

code_dir = fileparts(mfilename('fullpath'));

% the code root holds get_paths.m and this file only; another .m file there is
% most often one the refactor moved, recreated by an editor tab saved after the
% move, and it would run in place of the moved file
listing = dir(fullfile(code_dir, '*.m'));
names = {listing.name};
stray = names(~ismember(lower(names), {'get_paths.m', 'sep_setup_paths.m'}));
for k = 1:numel(stray)
    warning(['sep_setup_paths: %s is at the code root (%s), which holds only ' ...
             'get_paths.m and sep_setup_paths.m. If the refactor moved it, ' ...
             'docs/refactor_name_map.csv says where: carry any edit over to the ' ...
             'moved file, then delete this one.'], stray{k}, code_dir);
end

% the vendored toolboxes first, with their subfolders, so that our folders,
% added after them, end up above them
toolboxes = {'LightSuite', 'yamlmatlab', 'matlab_elastix', 'BioformatsImage'};
for k = 1:numel(toolboxes)
    addpath(genpath(fullfile(code_dir, 'third_party', toolboxes{k})));
end

% our folders, one by one: a new folder is added to this list by hand
folders = {'', 'common', 'atlas', fullfile('atlas', 'qc'), 'preprocessing', ...
    fullfile('preprocessing', 'pipeline'), 'registration', ...
    fullfile('registration', 'pipeline'), fullfile('registration', 'qc'), ...
    fullfile('registration', 'annotation_gui'), ...
    'group_comparison', fullfile('group_comparison', 'pipeline'), 'adult_matlab'};
folders = cellfun(@(f) fullfile(code_dir, f), folders, 'UniformOutput', false);
addpath(folders{:});

end
