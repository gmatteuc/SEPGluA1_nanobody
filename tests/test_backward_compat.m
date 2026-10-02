%% test_backward_compat
% ===== Check that the cohort registry and the paths resolve as before =====
%
% Regression test of the shared definitions: the cohort registry, get_atlas
% and get_paths must resolve exactly the mice, groups and paths of the mouse
% lists and folders the drivers wrote out before them, so the adult results
% stay reproducible while new cohorts are added. Prints PASS or FAIL for each
% check, and stops with an error if any fails.
%
% Run it after any change to get_cohort.m, get_atlas.m, or the cohort
% handling inside the drivers (run_extract_and_center to
% run_group_differences). It checks the code it lives in, on the data root
% get_paths gives (SEP_DATA_ROOT, or the data folder next to the code), and
% reads the lab share, when it is reachable, to find each young mouse's raw
% files. Run sep_setup_paths first, once per MATLAB session.

clear all
close all
clc

%% The code under test

% the code folder this file sits in (tests\ inside it), and the one the path
% resolves: they must be the same, or the test would check other code
code_dir = fileparts(fileparts(mfilename('fullpath')));
paths = get_paths();
if ~strcmpi(paths.code, code_dir)
    error(['test_backward_compat: get_paths resolves to %s, not to this code folder ' ...
           '(%s). Run restoredefaultpath, cd to the code root, then sep_setup_paths.'], ...
           paths.code, code_dir);
end

n_pass = 0;
n_fail = 0;

%% Adult mouse order and groups are frozen

% the mouse list the drivers wrote out before the registry existed
legacy_mice = { ...
    'MG691_Gria1', 'MG692_Gria1', 'MG693_Gria1', 'MG736_Gria1', 'MG737_Gria1', ...
    'CGF027_Gria1', 'CGF028_Gria1', 'CGF033_Gria1', 'CGF034_Gria1', 'CGF035_Gria1', ...
    'MG705_Gria1', 'MG706_Gria1', 'MG709_Gria1', 'MG716_Gria1', 'MG718_Gria1', ...
    'MG725_Gria1', 'MG727_Gria1'};
legacy_types = { ...
    'rws', 'rws', 'rws', 'rws', 'rws', ...
    'naive', 'naive', 'naive', 'naive', 'naive', ...
    'behavior', 'behavior', 'behavior', 'behavior', 'behavior', 'behavior', 'behavior'};

adults = get_cohort('groups', {'rws', 'naive', 'behavior'});
[n_pass, n_fail] = check(isequal({adults.name},  legacy_mice), ...
    'adult mouse list matches the legacy literal, in order', n_pass, n_fail);
[n_pass, n_fail] = check(isequal({adults.group}, legacy_types), ...
    'adult group assignment matches the legacy literal', n_pass, n_fail);

%% Selection by index still picks the same mice

% run_nano_equalisation selects mice_to_process = 1:17; run_residual_correction,
% run_annotate_artifacts and run_register_to_atlas a numeric mouse_idx
all_mice = get_cohort();
[n_pass, n_fail] = check(isequal({all_mice(1:17).name}, legacy_mice), ...
    'indices 1:17 of the full registry are still the 17 adults', n_pass, n_fail);

%% The mouse folders match the legacy paths

% the drivers built [<data root>, '\', mouse_type, '\'] + mouse_name, with the
% data root written out
ok_paths = true;
for i = 1:numel(adults)
    legacy_dir = fullfile(paths.data, legacy_types{i}, legacy_mice{i});
    if ~strcmpi(adults(i).base_dir, legacy_dir)
        fprintf('    MISMATCH %s: %s ~= %s\n', adults(i).name, adults(i).base_dir, ...
            legacy_dir);
        ok_paths = false;
    end
end
[n_pass, n_fail] = check(ok_paths, 'adult base_dir matches legacy path construction', ...
    n_pass, n_fail);

%% The adult folders exist on disk

missing = {};
for i = 1:numel(adults)
    if ~exist(adults(i).base_dir, 'dir')
        missing{end+1} = adults(i).name; %#ok<SAGROW>
    end
end
[n_pass, n_fail] = check(isempty(missing), ...
    sprintf('all adult dirs exist on disk (missing: %s)', strjoin(missing, ', ')), ...
    n_pass, n_fail);

%% The CCF atlas resolves to the folder and files the drivers used

atlas = get_atlas('ccf');
[n_pass, n_fail] = check(strcmpi(atlas.dir, fullfile(paths.data, 'atlas')), ...
    'get_atlas(''ccf'').dir is <data root>\atlas, where the allenDir literal pointed', ...
    n_pass, n_fail);
[n_pass, n_fail] = check(strcmp(atlas.annotation_file, 'annotation_10.nii.gz') && ...
    strcmp(atlas.template_file, 'average_template_10.nii.gz'), ...
    ['get_atlas(''ccf'') filenames match those used by run_residual_correction/' ...
     'run_collect_by_group/alignSliceVolume'], n_pass, n_fail);

%% get_paths gives the folders around this code folder

% the project folders are worked out from where the code sits, which lets the
% tree be copied to another drive; built again here on their own and compared
data_dir = fullfile(fileparts(code_dir), 'data');
if ~isempty(getenv('SEP_DATA_ROOT'))
    data_dir = getenv('SEP_DATA_ROOT');
end

expected = { ...
    'data',       paths.data,       data_dir
    'atlas',      paths.atlas,      fullfile(data_dir, 'atlas')
    'code',       paths.code,       code_dir
    'lightsuite', paths.lightsuite, fullfile(code_dir, 'third_party', 'LightSuite')
    'yaml',       paths.yaml,       fullfile(code_dir, 'third_party', 'yamlmatlab')
    'elastix',    paths.elastix,    fullfile(code_dir, 'third_party', 'matlab_elastix')
    'bioformats', paths.bioformats, fullfile(code_dir, 'third_party', 'BioformatsImage')
    'auto_annotation', paths.auto_annotation, ...
                  fullfile(code_dir, 'registration', 'auto_annotation')};

ok_paths = true;
for i = 1:size(expected, 1)
    if ~strcmpi(expected{i,2}, expected{i,3})
        fprintf('    MISMATCH %s: %s ~= %s\n', expected{i,1}, expected{i,2}, ...
            expected{i,3});
        ok_paths = false;
    end
end
[n_pass, n_fail] = check(ok_paths, ...
    'get_paths() gives the folders laid out around this code folder', n_pass, n_fail);

% a derived path that points nowhere is worse than a literal one
missing_dirs = {};
for i = 1:size(expected, 1)
    if ~exist(expected{i,2}, 'dir')
        missing_dirs{end+1} = expected{i,1}; %#ok<SAGROW>
    end
end
[n_pass, n_fail] = check(isempty(missing_dirs), ...
    sprintf('every derived path exists on disk (missing: %s)', ...
    strjoin(missing_dirs, ', ')), ...
    n_pass, n_fail);

%% The young cohort is consistent

% the size of Sami's list of good-quality brains, as extended on 12 Aug 2026
expected_young = 14;
young = get_cohort('groups', 'young');
[n_pass, n_fail] = check(numel(young) == expected_young, ...
    sprintf('young cohort has %d mice (got %d)', expected_young, numel(young)), ...
    n_pass, n_fail);

% the age in the registry must match the age in the folder name
ok_age = true;
for i = 1:numel(young)
    tok = regexp(young(i).name, '_P(\d+)$', 'tokens', 'once');
    if isempty(tok) || str2double(tok{1}) ~= young(i).age_days
        fprintf('    MISMATCH %s: registry age %g\n', young(i).name, young(i).age_days);
        ok_age = false;
    end
end
[n_pass, n_fail] = check(ok_age, 'young age_days matches the age in each folder name', ...
    n_pass, n_fail);

% share_subdir must point at where the .czi files are on the share (read only)
share_root = 'S:\ElboustaniLab\#SHARE\Data';
if exist(share_root, 'dir')
    ok_share = true;
    for i = 1:numel(young)
        d = fullfile(share_root, young(i).name, young(i).share_subdir);
        if isempty(dir(fullfile(d, '*.czi')))
            fprintf('    MISMATCH %s: no .czi under %s\n', young(i).name, d);
            ok_share = false;
        end
    end
    [n_pass, n_fail] = check(ok_share, ...
        'young share_subdir locates the raw .czi on the share', n_pass, n_fail);
else
    fprintf('  SKIP  share not reachable, cannot check share_subdir\n');
end

%% Summary

fprintf('\n%s\n', repmat('=', [1 60]));
fprintf('backward-compat tests: %d passed, %d failed\n', n_pass, n_fail);
fprintf('%s\n', repmat('=', [1 60]));
if n_fail > 0
    error('Backward-compatibility regression detected.');
end

% ===== Local functions =====

function [n_pass, n_fail] = check(cond, label, n_pass, n_fail)
% Print PASS or FAIL for one check, and add it to the tally.

if cond
    fprintf('  PASS  %s\n', label);
    n_pass = n_pass + 1;
else
    fprintf('  FAIL  %s\n', label);
    n_fail = n_fail + 1;
end

end
