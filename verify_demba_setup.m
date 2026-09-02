clear all
close all
clc

% /// Check that the DeMBA P20 path is actually set up correctly ///
% Registering to the wrong atlas, or to the right atlas in the wrong label
% space, produces a result that looks entirely plausible. Every one of the
% checks below failed silently at some point today, so they are all worth
% asserting before a manual annotation session is spent on top of them.
%
% Run this before annotating a young brain. It touches nothing.

paths = get_paths();
fprintf('=== DeMBA P20 setup check ===\n\n');

n_pass = 0;
n_fail = 0;

%% 1. get_atlas resolves, and only one atlas is on the path

atlas = get_atlas('demba_p20');
resolved = which(atlas.template_file);
[ok, n_pass, n_fail] = report('atlas on path resolves to demba_p20', ...
    strcmpi(fileparts(resolved), atlas.dir), n_pass, n_fail);
fprintf('      %s\n', resolved);

% Whole path ENTRIES, not substrings: 'atlas_demba_p20' contains 'atlas', so a
% substring test reports the adult dir as present whenever the DeMBA one is.
ccf_dir = paths.atlas;
path_entries = strsplit(path, pathsep);
[~, n_pass, n_fail] = report('adult atlas dir is NOT also on the path', ...
    ~any(strcmpi(path_entries, ccf_dir)), n_pass, n_fail);

%% 2. The annotation is in parcellation_index space, not structure IDs

av_y = niftiread(fullfile(atlas.dir, atlas.annotation_file));
av_a = niftiread(fullfile(ccf_dir, 'annotation_10.nii.gz'));

lab_y = unique(av_y(:)); lab_y = lab_y(lab_y ~= 0);
lab_a = unique(av_a(:)); lab_a = lab_a(lab_a ~= 0);
shared = numel(intersect(lab_y, lab_a));

[~, n_pass, n_fail] = report(sprintf(...
    'annotation is in parcellation_index space (%d/%d labels shared with CCF)', ...
    shared, numel(lab_y)), shared / numel(lab_y) > 0.95, n_pass, n_fail);
fprintf('      if this drops to about half, it is the raw BrainGlobe volume\n');
fprintf('      in Allen structure IDs and every region lookup would be wrong\n');

% Caudoputamen is index 662 in this space and structure 672 in the other one
[~, n_pass, n_fail] = report('CP resolves as parcellation_index 662 in both', ...
    any(av_y(:) == 662) && any(av_a(:) == 662), n_pass, n_fail);

%% 3. Resolution and crop

[~, n_pass, n_fail] = report(sprintf('atlas res_um is %g', atlas.res_um), ...
    atlas.res_um == 20, n_pass, n_fail);

crop = atlas.default_aplims;
[~, n_pass, n_fail] = report(sprintf('crop [%d %d] is inside the volume', crop), ...
    crop(1) >= 1 && crop(2) <= size(av_y, 1), n_pass, n_fail);

cropped = av_y(crop(1):crop(2), :, :);
per_plane = squeeze(sum(sum(cropped > 0, 2), 3));
[~, n_pass, n_fail] = report('no empty atlas planes inside the crop', ...
    all(per_plane > 0), n_pass, n_fail);
fprintf('      crop is %d planes = %.2f mm, brain fraction %.4f\n', ...
    size(cropped, 1), size(cropped, 1) * atlas.res_um / 1000, ...
    nnz(cropped > 0) / numel(cropped));

%% 4. Every selected mouse agrees with the atlas

get_cohort('verify');
young = get_cohort('groups', {'young'});

% Only the mice this atlas is FOR. An age-matched atlas is valid for its own
% age and nothing else, so checking a P36 brain against the P20 template would
% be the wrong test -- those brains need their own DeMBA age.
at_age    = [young.age_days] == atlas.age_days;
other_age = ~at_age;

bad_settings = {};
for k = find(at_age)
    f = fullfile(young(k).base_dir, 'local_settings.txt');
    if ~exist(f, 'file'), continue, end
    txt = fileread(f);
    px  = str2double(regexp(txt, 'px_atlas\s*=\s*([\d.]+)', 'tokens', 'once'));
    lim = regexp(txt, 'atlasaplims\s*=\s*\[(\d+)\s+(\d+)\]', 'tokens', 'once');
    if isempty(lim), continue, end
    lim = [str2double(lim{1}) str2double(lim{2})];
    if px ~= atlas.res_um || ~isequal(lim, crop)
        bad_settings{end+1} = sprintf('%s (px_atlas %g, aplims [%d %d])', ...
            young(k).name, px, lim); %#ok<SAGROW>
    end
end
[~, n_pass, n_fail] = report(sprintf(...
    'all %d P%g mice have local_settings matching this atlas', ...
    nnz(at_age), atlas.age_days), isempty(bad_settings), n_pass, n_fail);
for k = 1:numel(bad_settings)
    fprintf('      MISMATCH: %s\n', bad_settings{k});
end

% The rest of the cohort is not a failure, it is unfinished scope
if any(other_age)
    ages = unique([young(other_age).age_days]);
    fprintf('  [note] %d young mice are at other ages (P%s) and each needs its\n', ...
        nnz(other_age), strjoin(arrayfun(@(a) num2str(a), ages, ...
        'UniformOutput', false), ', P'));
    fprintf('         own DeMBA age before it can be registered age-matched.\n');
    fprintf('         BrainGlobe serves every day P4-P56, so it is a download,\n');
    fprintf('         a label remap and a crop measurement each.\n');
end

%% 5. The vendored fixes are in place

src = fileread(fullfile(paths.lightsuite, 'slice_module', 'alignSliceVolume.m'));
[~, n_pass, n_fail] = report('alignSliceVolume takes allenres from px_atlas', ...
    contains(src, 'regopts.allenres     = sliceinfo.px_atlas'), n_pass, n_fail);

src = fileread(fullfile(paths.lightsuite, 'slice_module', 'registerSlicesToAtlas.m'));
[~, n_pass, n_fail] = report('registerSlicesToAtlas control-point placeholder is 0x4', ...
    contains(src, 'zeros(0,4)'), n_pass, n_fail);

src = fileread(fullfile(paths.code, 'P4_register_to_atlas.m'));
[~, n_pass, n_fail] = report('P4 re-parses local_settings instead of MATLAB''s builtin', ...
    contains(src, 'parseSettingsFile(settings_name)'), n_pass, n_fail);
[~, n_pass, n_fail] = report('P4 atlas_key is demba_p20', ...
    contains(src, "atlas_key = 'demba_p20';"), n_pass, n_fail);

%% 6. Which brains are ready to annotate

fprintf('\n--- readiness ---\n');
for k = 1:numel(young)
    d = fullfile(young(k).base_dir, 'lightsuite');
    has_dec = exist(fullfile(d, 'volume_for_ordering_processing_decisions.txt'), 'file') == 2;
    regopts_f = fullfile(d, 'regopts.mat');
    has_reg = exist(regopts_f, 'file') == 2;
    has_ins = exist(fullfile(d, 'volume_for_inspection.tiff'), 'file') == 2;
    has_cps = exist(fullfile(d, 'atlas2histology_tform.mat'), 'file') == 2;

    % An existing regopts.mat is not enough -- it records which atlas the brain
    % was aligned AGAINST, and aligning to one atlas then annotating against
    % another is exactly the silent error this script exists to catch.
    aligned_here = false;
    if has_reg
        R = load(regopts_f, 'allenres', 'atlasaplims');
        aligned_here = isfield(R, 'allenres') && R.allenres == atlas.res_um && ...
                       isfield(R, 'atlasaplims') && isequal(R.atlasaplims(:)', crop);
    end

    if ~ismember(young(k).age_days, atlas.age_days)
        state = sprintf('P%g - needs its own DeMBA age', young(k).age_days);
    elseif ~has_dec
        state = 'needs P1bis ordering';
    elseif ~(has_reg && has_ins)
        state = 'ordered, needs P4 align';
    elseif ~aligned_here
        state = sprintf('STALE - aligned to a different atlas (allenres %g), re-run P4 align', ...
                        R.allenres);
    elseif has_cps
        state = 'HAS control points';
    else
        state = 'ready to annotate';
    end
    fprintf('  %-26s %s\n', young(k).name, state);
end

%% Verdict

fprintf('\n=== %d passed, %d failed ===\n', n_pass, n_fail);
if n_fail > 0
    fprintf('DO NOT annotate until the failures above are fixed.\n');
else
    fprintf('P4 align + annotate is safe to run against demba_p20.\n');
end

fprintf(['\nStill outstanding, and NOT checked by this script: P5 onward\n' ...
         'hardcode the adult atlas and crop. The two cohorts land on\n' ...
         'different grids ([900 800 1140] adult, [994 800 1140] young), so\n' ...
         'those scripts must be made atlas-aware per cohort before any young\n' ...
         'data reaches them.\n']);

%% Local function: one check

function [ok, n_pass, n_fail] = report(name, ok, n_pass, n_fail)
    if ok
        fprintf('  [ ok ] %s\n', name);
        n_pass = n_pass + 1;
    else
        fprintf('  [FAIL] %s\n', name);
        n_fail = n_fail + 1;
    end
end
