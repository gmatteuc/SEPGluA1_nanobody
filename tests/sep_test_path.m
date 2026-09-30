function sep_test_path()
% SEP_TEST_PATH Check that no function name is defined twice on the project path.
%
% Run it in a fresh session, after sep_setup_paths, from the code root:
%
%   restoredefaultpath; cd('D:\sep_histology\code'); sep_setup_paths
%   run(fullfile('tests', 'sep_test_path.m'))
%
% MATLAB takes the first function of a name it finds on the path, silently.
% Two files of the same name, ours or vendored, mean that one of them never
% runs, and which one depends on the order of the path. Both this project and
% the imaging repository (D:\dendrites\code) run in the same MATLAB, so a name
% both define would redirect one project into the other's code.
%
% Lists every function on the project path: the folders under this code root
% that are on the path (third_party\ included; tests\ and tools\, which a check
% may add by hand, left out), with the .m, .p and MEX files and the @class
% folders of each. Fails on a name defined in two folders, compared without
% case (core.ignorecase is true, and Windows does not tell the files apart),
% except the known vendored duplicates of is_known_duplicate below. When
% D:\dendrites\code exists, it also fails on a name that project defines too.
%
% Then runs which -all on each of our names (the folders outside third_party\)
% and fails if one of them shadows a MATLAB or toolbox function, or is shadowed
% by anything. The vendored names that shadow a MATLAB function are listed,
% not failed: the drivers have put LightSuite above MATLAB's functions from
% the start, and sep_setup_paths keeps it there.

code_dir = fileparts(fileparts(mfilename('fullpath')));
dendrites_dir = 'D:\dendrites\code';

% the path must be the one sep_setup_paths of this code folder sets up
if ~strcmpi(fileparts(which('sep_setup_paths')), code_dir)
    error(['sep_test_path: sep_setup_paths resolves to %s, not to this code folder ' ...
           '(%s). Run restoredefaultpath, cd to the code root, then sep_setup_paths.'], ...
           which('sep_setup_paths'), code_dir);
end

% the project's folders on the path
entries = strsplit(path, pathsep);
inside = strcmpi(entries, code_dir) | startsWith(lower(entries), lower([code_dir filesep]));
by_hand = startsWith(lower(entries), lower(fullfile(code_dir, 'tests'))) | ...
    startsWith(lower(entries), lower(fullfile(code_dir, 'tools')));
folders = entries(inside & ~by_hand);
vendored = startsWith(lower(folders), lower([fullfile(code_dir, 'third_party') filesep]));
fprintf('sep_test_path: %d folders of %s on the path (%d vendored)\n', ...
    numel(folders), code_dir, sum(vendored));

% every function name on them, once per folder
names = {};
files = {};
is_ours = false(0, 1);
for k = 1:numel(folders)
    [n, f] = functions_in(folders{k});
    names = [names; n]; %#ok<AGROW>
    files = [files; f]; %#ok<AGROW>
    is_ours = [is_ours; repmat(~vendored(k), numel(n), 1)]; %#ok<AGROW>
end
fprintf('  %d function names (%d ours, %d vendored)\n', numel(names), sum(is_ours), ...
    sum(~is_ours));

% the imaging repository's names, when its checkout is on this computer
other_names = {};
other_files = {};
if isfolder(dendrites_dir)
    other = strsplit(genpath(dendrites_dir), pathsep);
    other = other(~cellfun(@isempty, other));
    for k = 1:numel(other)
        [n, f] = functions_in(other{k});
        other_names = [other_names; n]; %#ok<AGROW>
        other_files = [other_files; f]; %#ok<AGROW>
    end
    fprintf('  %d function names in %s\n', numel(other_names), dendrites_dir);
else
    fprintf('  %s not found: its names are not checked\n', dendrites_dir);
end

n_fail = 0;

% the same name in two folders of the project, without case
[~, ~, group] = unique(lower(names));
counts = accumarray(group, 1);
twice = find(counts > 1);
n_known = 0;
for k = 1:numel(twice)
    clash = files(group == twice(k));
    name = names{find(group == twice(k), 1)};
    if is_known_duplicate(name, clash)
        fprintf('  known %s is defined %d times (see is_known_duplicate):\n', name, ...
            numel(clash));
        n_known = n_known + 1;
    else
        fprintf('  FAIL  %s is defined %d times:\n', name, numel(clash));
        n_fail = n_fail + 1;
    end
    fprintf('          %s\n', clash{:});
end
if numel(twice) == n_known
    fprintf('  PASS  no name is defined twice in the project (%d known duplicates)\n', ...
        n_known);
end

% a name both projects define; duplicates inside the imaging repository are
% its own business (its third-party folders are not all on its path)
n_shared = 0;
for k = 1:numel(names)
    shared = strcmpi(other_names, names{k});
    if any(shared)
        fprintf('  FAIL  %s is also defined in the imaging repository:\n', files{k});
        fprintf('          %s\n', other_files{shared});
        n_fail = n_fail + 1;
        n_shared = n_shared + 1;
    end
end
if isfolder(dendrites_dir) && n_shared == 0
    fprintf('  PASS  no name is shared with %s\n', dendrites_dir);
end

% which -all on each name: ours must neither shadow nor be shadowed; the
% vendored ones that shadow MATLAB's are listed
matlab_dir = lower(matlabroot);
n_shadowing_vendored = 0;
for k = 1:numel(names)
    found = which(names{k}, '-all');
    [own_dir, ~] = fileparts(files{k});
    others = found(~startsWith(lower(found), lower([own_dir filesep])));
    if isempty(others) || is_known_duplicate(names{k}, [files(k); others(:)])
        continue
    end
    from_matlab = startsWith(lower(others), matlab_dir) | startsWith(others, 'built-in');
    if is_ours(k)
        fprintf('  FAIL  %s (%s) shares its name with:\n', names{k}, files{k});
        fprintf('          %s\n', others{:});
        n_fail = n_fail + 1;
    elseif all(from_matlab)
        fprintf('  note  vendored %s shadows %s\n', files{k}, strjoin(others, ', '));
        n_shadowing_vendored = n_shadowing_vendored + 1;
    else
        fprintf('  FAIL  vendored %s shares its name with:\n', files{k});
        fprintf('          %s\n', others{:});
        n_fail = n_fail + 1;
    end
end
if n_shadowing_vendored > 0
    fprintf(['  %d vendored functions shadow a MATLAB function, as they did when the ' ...
             'drivers added them\n'], n_shadowing_vendored);
end

fprintf('sep_test_path: %d failures\n', n_fail);
if n_fail > 0
    error('sep_test_path: %d name clashes on the project path (listed above).', n_fail);
end
end

% ===== Local functions =====

function known = is_known_duplicate(name, clash)
% Duplicates inside a vendored package that nothing calls, left as they came.
% matlab_elastix ships two example scripts named RUN_ALL (one per example
% folder); the drivers have always put both on the path with genpath, and
% neither is called by any code here.

examples = fullfile('third_party', 'matlab_elastix', 'MelastiX_examples');
known = strcmpi(name, 'RUN_ALL') && all(contains(lower(clash), lower(examples)));
end

function [names, files] = functions_in(folder)
% The functions a folder puts on the path: .m, .p and MEX files and @class
% folders, one entry per name (a MEX file beside its .m help file is one
% function).

names = {};
files = {};
listing = dir(folder);
for k = 1:numel(listing)
    item = listing(k);
    [~, base, ext] = fileparts(item.name);
    if item.isdir && startsWith(item.name, '@')
        name = item.name(2:end);
    elseif ~item.isdir && (strcmpi(ext, '.m') || strcmpi(ext, '.p') || ...
            startsWith(lower(ext), '.mex'))
        name = base;
    else
        continue
    end
    if any(strcmpi(names, name))
        continue
    end
    names{end+1, 1} = name; %#ok<AGROW>
    files{end+1, 1} = fullfile(folder, item.name); %#ok<AGROW>
end
end
