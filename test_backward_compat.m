clear all
close all
clc

% /// Backward-compatibility regression tests ///
% Asserts that the cohort-registry refactor resolves exactly the same mice,
% groups and paths as the pre-refactor hardcoded literals, so adult results
% remain reproducible while new cohorts are added.
%
% Run this after ANY change to get_cohort.m, get_atlas.m, or the cohort
% handling inside P1-P7bis.

%% Add paths

addpath('D:\sep_histology\code')

n_pass = 0;
n_fail = 0;

%% Test 1: adult mouse order and group assignment are frozen

% The exact literal that P1-P7bis carried before the registry existed.
legacy_mice = {'MG691_Gria1', 'MG692_Gria1', 'MG693_Gria1', 'MG736_Gria1', 'MG737_Gria1', 'CGF027_Gria1', 'CGF028_Gria1', 'CGF033_Gria1', 'CGF034_Gria1', 'CGF035_Gria1','MG705_Gria1', 'MG706_Gria1', 'MG709_Gria1', 'MG716_Gria1', 'MG718_Gria1', 'MG725_Gria1', 'MG727_Gria1'};
legacy_types = {'rws','rws','rws','rws','rws','naive','naive','naive','naive','naive','behavior','behavior','behavior','behavior','behavior','behavior','behavior'};

adults = get_cohort('groups', {'rws','naive','behavior'});
[n_pass, n_fail] = check(isequal({adults.name},  legacy_mice), ...
    'adult mouse list matches the legacy literal, in order', n_pass, n_fail);
[n_pass, n_fail] = check(isequal({adults.group}, legacy_types), ...
    'adult group assignment matches the legacy literal', n_pass, n_fail);

%% Test 2: index-based selection still picks the same mice

% P2bis used `mice_to_process = 1:17`; P2/P3/P4 used numeric mouse_idx.
all_mice = get_cohort();
[n_pass, n_fail] = check(isequal({all_mice(1:17).name}, legacy_mice), ...
    'indices 1:17 of the full registry are still the 17 adults', n_pass, n_fail);

%% Test 3: resolved paths match the legacy path construction

% Legacy scripts built: ['D:\sep_histology\data\', mouse_type, '\'] + mouse_name
ok_paths = true;
for i = 1:numel(adults)
    legacy_dir = fullfile('D:\sep_histology\data', legacy_types{i}, legacy_mice{i});
    if ~strcmpi(adults(i).base_dir, legacy_dir)
        fprintf('    MISMATCH %s: %s ~= %s\n', adults(i).name, adults(i).base_dir, legacy_dir);
        ok_paths = false;
    end
end
[n_pass, n_fail] = check(ok_paths, 'adult base_dir matches legacy path construction', n_pass, n_fail);

%% Test 4: adult data directories still exist on disk

missing = {};
for i = 1:numel(adults)
    if ~exist(adults(i).base_dir, 'dir'), missing{end+1} = adults(i).name; end %#ok<SAGROW>
end
[n_pass, n_fail] = check(isempty(missing), ...
    sprintf('all adult dirs exist on disk (missing: %s)', strjoin(missing, ', ')), n_pass, n_fail);

%% Test 5: atlas resolution reproduces the historical hardcoded paths

atlas = get_atlas('ccf');
[n_pass, n_fail] = check(strcmpi(atlas.dir, 'D:\sep_histology\data\atlas'), ...
    'get_atlas(''ccf'').dir matches the hardcoded allenDir literal', n_pass, n_fail);
[n_pass, n_fail] = check(strcmp(atlas.annotation_file, 'annotation_10.nii.gz') && ...
    strcmp(atlas.template_file, 'average_template_10.nii.gz'), ...
    'get_atlas(''ccf'') filenames match those used by P2/P5/alignSliceVolume', n_pass, n_fail);

%% Test 6: young cohort registry is internally consistent

expected_young = 14;   % Sami's good-quality list, extended 2026-08-12
young = get_cohort('groups', 'young');
[n_pass, n_fail] = check(numel(young) == expected_young, ...
    sprintf('young cohort has %d mice (got %d)', expected_young, numel(young)), n_pass, n_fail);

% Age in the registry must match the age embedded in the folder name
ok_age = true;
for i = 1:numel(young)
    tok = regexp(young(i).name, '_P(\d+)$', 'tokens', 'once');
    if isempty(tok) || str2double(tok{1}) ~= young(i).age_days
        fprintf('    MISMATCH %s: registry age %g\n', young(i).name, young(i).age_days);
        ok_age = false;
    end
end
[n_pass, n_fail] = check(ok_age, 'young age_days matches the age in each folder name', n_pass, n_fail);

% share_subdir must point at where the .czi actually live on the share
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
    [n_pass, n_fail] = check(ok_share, 'young share_subdir locates the raw .czi on the share', n_pass, n_fail);
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

%% Local function: report and tally a single assertion

function [n_pass, n_fail] = check(cond, label, n_pass, n_fail)
if cond
    fprintf('  PASS  %s\n', label);
    n_pass = n_pass + 1;
else
    fprintf('  FAIL  %s\n', label);
    n_fail = n_fail + 1;
end
end
