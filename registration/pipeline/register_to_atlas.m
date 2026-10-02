function register_to_atlas(run_settings)
%REGISTER_TO_ATLAS  Register the selected mice to the atlas, in one run mode.
%   REGISTER_TO_ATLAS(run_settings) does the work of run_register_to_atlas,
%   which sets the fields of run_settings (paths, groups_to_process,
%   mice_to_process, run_mode, atlas_extent_slices, atlas_key,
%   correction_type, use_equalized_nano, allow_image_only_registration) and
%   says what each one, and each run mode, does.
%
%   The modes 'angle' and 'annotate' open a LightSuite window and return
%   while it is still open: the window keeps what it needs in its own data
%   and saves into the mouse's lightsuite folder, so nothing it does depends
%   on this function's variables.
%
%   The 'align' branch loads the correction outputs with load() and no
%   output, so some of the variables it reads (correctedVol, scaledautoVol,
%   bg_mask_vol, equalized_volume) are assigned nowhere in the code: a
%   function finds them only because no function on the path has their
%   names. The files also hold correction_type, which replaces the setting
%   when they are loaded, as it did in the script.

% settings of run_register_to_atlas, under the names the code below uses
paths = run_settings.paths;
groups_to_process = run_settings.groups_to_process;
mice_to_process = run_settings.mice_to_process;
run_mode = run_settings.run_mode;
atlas_extent_slices = run_settings.atlas_extent_slices;
atlas_key = run_settings.atlas_key;
correction_type = run_settings.correction_type;
use_equalized_nano = run_settings.use_equalized_nano;
allow_image_only_registration = run_settings.allow_image_only_registration;

%% Add paths

% get_atlas puts the chosen atlas folder on the path and takes the others off:
% every atlas folder holds files of the same names, which LightSuite finds with
% which(), so never addpath an atlas folder by hand (the toolboxes are on the
% path from sep_setup_paths)
atlas = get_atlas(atlas_key);

%% Resolve cohort

% check that the registry still lists the adults in their legacy order
get_cohort('verify');

% the mice named, or else every mouse of the groups
if isempty(mice_to_process)
    cohort = get_cohort('groups', groups_to_process);
else
    cohort = get_cohort('names', mice_to_process);
end

% stop on a mistyped mode, which would otherwise fall through to 'align' and redo
% an hour of bridging nobody asked for
if ~ismember(run_mode, {'align', 'angle', 'annotate', 'autoannotate', 'register'})
    error('run_register_to_atlas: unknown run_mode ''%s'' (use ''align'', ''angle'', ''annotate'', ''autoannotate'' or ''register'').', run_mode);
end

fprintf('run_register_to_atlas: %d mouse/mice selected, mode ''%s'', atlas ''%s''.\n', ...
    numel(cohort), run_mode, atlas.key);

%% Check the selected mice

check_atlas_age(cohort, atlas);
check_atlas_resolution(cohort, atlas);

% aligning again rewrites the atlas block that control points, anchors and
% proposals are counted in, so an annotated brain would silently lose its
% annotation (the same rule as 'angle'); checked for every selected mouse before
% any of them is aligned
if strcmp(run_mode, 'align')
    refuse_annotated_mice(cohort);
end

%% Run the selected mode on each mouse

for mouse_idx = 1:numel(cohort)

    % the mouse
    mouse_name = cohort(mouse_idx).name;
    mouse_type = cohort(mouse_idx).group;

    fprintf('\n=== %s (%s) ===\n', mouse_name, mouse_type);

    % its folders
    base_dir = fullfile(paths.data, mouse_type);
    mouse_dir = fullfile(base_dir, mouse_name, '\lightsuite');
    correction_dir = fullfile(base_dir, mouse_name, '\lightsuite', 'correction_output');
    before_correction_dir = fullfile(base_dir, mouse_name, '\lightsuite', ...
        'volume_centered');
    processed_dir = fullfile(base_dir, mouse_name, '\lightsuite', ...
        'volume_centered_processed');
    aligned_dir = fullfile(mouse_dir, 'volume_aligned');
    volorder_dir = fullfile(mouse_dir, 'volume_for_ordering.tiff');

    % the modes after 'align' work off what 'align' wrote, so they skip the
    % expensive bridging below and go straight to their step
    if ismember(run_mode, {'angle', 'annotate', 'autoannotate', 'register'})

        opts = load_regopts(mouse_dir, mouse_name, atlas_extent_slices);

        switch run_mode

            case 'angle'
                set_cutting_angle(cohort, mouse_dir, mouse_name, opts, atlas);

            case 'annotate'
                annotate_control_points(cohort, mouse_dir, opts, atlas);

            case 'autoannotate'
                propose_control_points(mouse_dir, mouse_name);

            case 'register'
                check_control_points(mouse_dir, mouse_name, ...
                    allow_image_only_registration);
                register_slices(opts, mouse_dir, processed_dir, volorder_dir, ...
                    aligned_dir);

        end

        continue

    end

    % (auto) bridge the preprocessing results into LightSuite's own pipeline
    [sliceinfo, correction_type] = bridge_preprocessing(correction_dir, ...
        before_correction_dir, mouse_dir, processed_dir, volorder_dir, aligned_dir, ...
        correction_type, use_equalized_nano);

    % (auto) align the slices and start the registration
    align_slices(sliceinfo, base_dir, mouse_name, mouse_dir, atlas);

end

end

% ===== Local functions: checks before the run =====

function check_atlas_age(cohort, atlas)
% Stop unless every selected mouse is of the age the atlas represents.

% a wrong atlas gives a plausible-looking result, so check the age: an age-matched
% atlas holds for its own age only (a P36 brain belongs on neither P20 nor CCF)
for k = 1:numel(cohort)
    mouse_age = cohort(k).age_days;

    % the adult cohorts are not dated one by one (age_days NaN) and go to 'ccf',
    % whose age_days is 56
    if isnan(mouse_age)
        mouse_age = 56;
    end
    if mouse_age ~= atlas.age_days
        error(['run_register_to_atlas: %s is P%g but atlas ''%s'' represents P%g.\n' ...
               'Register each brain to the atlas for its own age, or add an\n' ...
               'entry for P%g to get_atlas.'], ...
               cohort(k).name, mouse_age, atlas.key, atlas.age_days, mouse_age);
    end
end
fprintf('run_register_to_atlas: all selected mice are P%g, matching atlas ''%s''.\n', ...
    atlas.age_days, atlas.key);
end

function check_atlas_resolution(cohort, atlas)
% Stop unless px_atlas in each mouse's local_settings.txt is the atlas
% resolution (a mouse without the file is skipped).

% px_atlas sets the AP scale of the reconstruction: a young brain left at
% px_atlas = 10 against a 20 um atlas is off by a factor of 2
for k = 1:numel(cohort)
    settings_name = fullfile(cohort(k).base_dir, 'local_settings.txt');
    if ~exist(settings_name, 'file')
        settings_name = fullfile(cohort(k).base_dir, 'lightsuite', 'local_settings.txt');
    end
    if ~exist(settings_name, 'file')
        continue
    end
    txt = fileread(settings_name);
    tok = regexp(txt, 'px_atlas\s*=\s*([\d.]+)', 'tokens', 'once');
    if ~isempty(tok) && str2double(tok{1}) ~= atlas.res_um
        error(['run_register_to_atlas: %s has px_atlas = %s but atlas ''%s'' is %g um.\n' ...
               'Fix px_atlas (and atlasaplims) in\n  %s'], ...
               cohort(k).name, tok{1}, atlas.key, atlas.res_um, settings_name);
    end
end
fprintf('run_register_to_atlas: atlas resolution agrees with local_settings for all selected mice.\n');
end

function refuse_annotated_mice(cohort)
% Stop if any selected mouse already has control points, anchors or a
% proposal, naming the files found.

annotation_files = {'atlas2histology_tform.mat', 'plane_anchors.mat', ...
                    'auto_atlas_planes.mat', 'auto_proposal_controlpoints.mat'};
annotated = {};
for k = 1:numel(cohort)
    found = annotation_files(cellfun(@(f) exist(fullfile(cohort(k).base_dir, ...
        'lightsuite', f), 'file') == 2, annotation_files));
    if ~isempty(found)
        annotated{end+1} = sprintf('  %s: %s', cohort(k).name, strjoin(found, ', ')); %#ok<AGROW>
    end
end
if ~isempty(annotated)
    error(['run_register_to_atlas: these mice already have an annotation, which aligning again would ' ...
           'invalidate:\n%s\nMove those files aside first if you really want to redo ' ...
           'both.'], strjoin(annotated, '\n'));
end
end

% ===== Local functions: the modes that follow align =====

function opts = load_regopts(mouse_dir, mouse_name, atlas_extent_slices)
% The registration options 'align' wrote (regopts.mat), pointed at this
% mouse's folder and given the AP margin, which is written back to the file.

regopts_name = fullfile(mouse_dir, 'regopts.mat');
if ~exist(regopts_name, 'file')
    error(['run_register_to_atlas: no regopts.mat for %s:\n  %s\n' ...
           'Run this script with run_mode = ''align'' for this mouse first.'], ...
           mouse_name, regopts_name);
end
opts = load(regopts_name);

% 'align' writes procpath as an absolute path, drive letter included, which on
% another machine or drive letter does not exist; it is always this folder
opts.procpath = mouse_dir;

% the AP margin (atlas_extent_slices): the GUI reads it from opts and the
% registration from regopts.mat, so it is written back. A mouse with control
% points keeps the margin it was annotated with, since the saved atlas planes are
% counted from the front of that range (the three P20 mice done before this
% setting existed reopen and register again exactly as before).
if opts.extentfactor ~= atlas_extent_slices
    if exist(fullfile(mouse_dir, 'atlas2histology_tform.mat'), 'file')
        fprintf(['  atlas margin around the stack stays at %d slices: this mouse already has ' ...
                 'control points (atlas_extent_slices = %d applies to new mice only)\n'], ...
            opts.extentfactor, atlas_extent_slices);
    else
        fprintf('  atlas margin around the stack: %d -> %d slices (written to regopts.mat)\n', ...
            opts.extentfactor, atlas_extent_slices);
        opts.extentfactor = atlas_extent_slices;
        regopts_disk = load(regopts_name);
        regopts_disk.extentfactor = atlas_extent_slices;
        save(regopts_name, '-struct', 'regopts_disk');
        clear regopts_disk
    end
end
end

function set_cutting_angle(cohort, mouse_dir, mouse_name, opts, atlas)
% Mode 'angle': open LightSuite's cutting-angle GUI on one mouse that has no
% control points yet.

if numel(cohort) > 1
    error('run_register_to_atlas: run_mode ''angle'' opens one GUI at a time; select a single mouse.');
end
if exist(fullfile(mouse_dir, 'atlas2histology_tform.mat'), 'file')
    error(['run_register_to_atlas: %s already has control points. The cutting angle changes the atlas ' ...
           'block those points are counted in, so it has to be set before annotating. ' ...
           'Move atlas2histology_tform.mat aside first if you really want to redo both.'], ...
           mouse_name);
end
angle_name = fullfile(mouse_dir, 'cutting_angle_data.mat');
if exist(angle_name, 'file')
    fprintf('  NOTE: a cutting angle is already saved and will be overwritten on close:\n    %s\n', angle_name);
end
fprintf('  opening the cutting-angle GUI against atlas ''%s''.\n', atlas.key);
fprintf(['  SHIFT+arrows tilt, wheel moves the plane, return saves the plane on the current slice ' ...
         '(3-5 slices, averaged), close the window to write the file.\n']);
determineCuttingAngleGUI(opts);
end

function annotate_control_points(cohort, mouse_dir, opts, atlas)
% Mode 'annotate': open LightSuite's control-point GUI on one mouse, with the
% automatic annotation's keys when it is installed.

% one GUI at a time, or the control points go into the wrong mouse's file
if numel(cohort) > 1
    error('run_register_to_atlas: run_mode ''annotate'' opens one GUI at a time; select a single mouse.');
end
tform_name = fullfile(mouse_dir, 'atlas2histology_tform.mat');
if exist(tform_name, 'file')
    fprintf('  NOTE: control points already exist and will be overwritten on save:\n    %s\n', tform_name);
end
if exist(fullfile(mouse_dir, 'cutting_angle_data.mat'), 'file')
    fprintf('  cutting angle: from cutting_angle_data.mat (set by hand).\n');
else
    fprintf('  cutting angle: the automatic rigid fit''s (no cutting_angle_data.mat; run_mode = ''angle'' to set it by eye).\n');
end

% the automatic annotation's part of the GUI is a plugin, in
% registration\annotation_gui, passed in when auto_annotate finds the engine;
% without it the GUI is LightSuite's with our fixes (PATCHES.md)
engine = auto_annotate('check');
if engine.ok
    opts.plugin = @auto_annotation_plugin;
else
    fprintf(['  the automatic annotation is not installed (%s):\n' ...
             '  the GUI opens without its keys (a, j, k, u).\n'], engine.message);
end
fprintf('  opening the control-point GUI against atlas ''%s''.\n', atlas.key);
fprintf('  place points on every slice, then SAVE and CLOSE, and re-run with run_mode = ''register''.\n');
matchControlPointsInSlices(opts);
end

function propose_control_points(mouse_dir, mouse_name)
% Mode 'autoannotate': propose control points for every slice from the anchor
% planes (auto_proposal_controlpoints.mat).

% a proposal is written next to the annotation, never into it: the GUI loads it
% for review, and only what is accepted there is saved as the annotation
for f = {'plane_anchors.mat', 'auto_atlas_planes.mat'}
    if ~exist(fullfile(mouse_dir, f{1}), 'file')
        error(['run_register_to_atlas: no %s for %s.\n' ...
               'Open run_mode = ''annotate'', set the plane on the suggested anchor\n' ...
               'slices (j jumps between them, a fixes the plane), save with s, then re-run.'], ...
               f{1}, mouse_name);
    end
end
if exist(fullfile(mouse_dir, 'auto_proposal_controlpoints.mat'), 'file')
    fprintf('  NOTE: a previous proposal will be replaced (accepted slices are in the annotation, not there).\n');
end
fprintf('  proposing control points for every slice from the anchor planes...\n');
out = auto_annotate('propose', mouse_dir);
if ~out.ok
    error('run_register_to_atlas: the automatic annotation failed for %s:\n%s', ...
        mouse_name, out.message);
end
fprintf(['  done. Review it: run_mode = ''annotate'' (k accepts a slice, u re-proposes it),\n' ...
         '  then run_mode = ''register''.\n']);
end

function check_control_points(mouse_dir, mouse_name, allow_image_only_registration)
% Mode 'register', before elastix: report the control points and stop on a
% missing file or an unpaired slice (image only when allowed).

% LightSuite does not register well from the images alone, and all 17 adults have
% points on every slice, so a missing annotation is an error, not a fallback
tform_name = fullfile(mouse_dir, 'atlas2histology_tform.mat');
if exist(tform_name, 'file')
    S_cp  = load(tform_name, 'histology_control_points', 'atlas_control_points');
    n_cp  = cellfun(@(c) size(c, 1), S_cp.histology_control_points);
    n_at  = cellfun(@(c) size(c, 1), S_cp.atlas_control_points);
    fprintf('  control points: %d slices, %d-%d points each (median %g)\n', ...
        numel(n_cp), min(n_cp), max(n_cp), median(n_cp));

    % points are paired by row, so a slice whose two lists differ in length has no
    % pairing; elastix would find out only on reaching it, an hour in
    bad = find(n_cp(:) ~= n_at(:))';
    if ~isempty(bad)
        error(['run_register_to_atlas: slice(s) %s have different numbers of histology and atlas points.\n' ...
               'Open run_mode = ''annotate'', go to each, and delete the unpaired point(s)\n' ...
               '(edit mode, d) or press c and re-place them.'], mat2str(bad));
    end

    % registerSlicesToAtlas uses a slice's points only when it has at least five,
    % and silently registers it from the images otherwise
    if any(n_cp < 5)
        fprintf('  NOTE: slice(s) %s have fewer than 5 points and will be registered from image only.\n', ...
            mat2str(find(n_cp(:) < 5)'));
    end

    % a mismatched pair (a point on the wrong structure, a left-right swap) drags
    % the slice's affine, so name the pairs far off a robust fit, by their numbers
    report_suspect_pairs(S_cp.atlas_control_points, S_cp.histology_control_points);
elseif allow_image_only_registration
    fprintf(['  no control points, and allow_image_only_registration is true.\n' ...
             '  Registering from image information alone -- diagnostic only,\n' ...
             '  do not compare the result against the adults.\n']);
else
    error(['run_register_to_atlas: no control points for %s:\n  %s\n' ...
           'LightSuite needs the manual points to register properly, and all 17\n' ...
           'adults have them on every slice. Run this script with\n' ...
           'run_mode = ''annotate'' for this mouse first.\n' ...
           'To register without them anyway (diagnostic only), set\n' ...
           'allow_image_only_registration = true.'], mouse_name, tform_name);
end
end

function report_suspect_pairs(acp, hcp)
% Per slice, print the pairs far off an affine fitted without the outliers, and
% the slices with points on several planes or off as a whole.

annot = annotation_settings();
rule  = annot.outlier_rule;
for k = 1:numel(acp)
    a = acp{k};
    h = hcp{k};
    if size(a,1) < 5 || size(a,1) ~= size(h,1)
        continue
    end

    % points on two atlas planes (a carried-over proposal accepted, a scroll, more
    % points): the registration takes the first point's plane for all of them, so
    % part of the slice is annotated against anatomy 0.2 mm away (seen on MG910)
    planes = unique(a(:,1));
    if numel(planes) > 1
        n_per = arrayfun(@(pl) nnz(a(:,1) == pl), planes);
        fprintf(['  CHECK slice %d: points sit on %d different atlas planes (%s, %s points each). ' ...
                 'Registration uses the first point''s plane for all of them -- clear the slice and ' ...
                 'redo it on one plane.\n'], k, numel(planes), mat2str(planes'), mat2str(n_per'));
    end

    % affine from atlas to histology, refitted up to 8 times without the pairs
    % over twice the median residual, so outliers stop pulling the fit to them
    src = a(:, [3 2]);
    dst = h(:, [3 2]);
    w = ones(size(src,1), 1);
    for it = 1:8
        tf = fitgeotform2d(src(w > 0.999, :), dst(w > 0.999, :), 'affine');
        r  = vecnorm(tf.transformPointsForward(src) - dst, 2, 2);
        s  = median(r) + 1e-6;
        w  = double(r <= 2*s);
        if nnz(w) < 5
            break
        end
    end

    % pairs further than the outlier rule of annotation_settings (3x the median
    % residual and 30 px), the rule the GUI's * mark uses too
    bad = find(r > max(rule.factor*median(r), rule.min_px))';
    if ~isempty(bad)
        fprintf('  CHECK slice %d: pair(s) %s sit %s px off the slice''s own affine (median %.0f). Wrong structure or a left/right swap?\n', ...
            k, mat2str(bad), mat2str(round(r(bad))'), median(r));
    end

    % a slice wrong as a whole, which the relative test above cannot see (MG910's
    % slice 35, median 93 px): a good slice sits at 5-10 px, and past rule.slice_px
    % (20) it is usually the wrong atlas plane, on which elastix tends to fail
    if median(r) > rule.slice_px
        fprintf(['  CHECK slice %d: the WHOLE slice is off -- median residual %.0f px over %d pairs. ' ...
                 'Wrong atlas plane, or points placed on two planes? Clear and redo it.\n'], ...
                 k, median(r), numel(r));
    end
end
end

function register_slices(opts, mouse_dir, processed_dir, volorder_dir, aligned_dir)
% Mode 'register': elastix refinement of every slice, then the registered
% volumes of all five channels.

transformparams = registerSlicesToAtlas(opts); %#ok<NASGU>

transformparams = load(fullfile(mouse_dir, 'transform_params.mat'));
S_slice = load(fullfile(mouse_dir, 'sliceinfo.mat'));
sliceinfo_new = S_slice.sliceinfo;
sliceinfo_new.channames   = {'DAPI','NANO','AUTO','DIFF','MASK'};
sliceinfo_new.slicevol    = processed_dir;
sliceinfo_new.procpath    = mouse_dir;
sliceinfo_new.volorder    = volorder_dir;
sliceinfo_new.slicevolfin = aligned_dir;
generateRegisteredSliceVolume(sliceinfo_new, transformparams);
end

% ===== Local functions: align =====

function [sliceinfo, correction_type] = bridge_preprocessing(correction_dir, ...
    before_correction_dir, mouse_dir, processed_dir, volorder_dir, aligned_dir, ...
    correction_type, use_equalized_nano)
% Mode 'align', first half: the corrected volumes, the nano and the invalid
% mask saved as LightSuite channels; correction_type as the files leave it.

% load the correction outputs and the artifact annotation
matfile0_name = fullfile(correction_dir, sprintf('corrected_volume_%s.mat', ...
    correction_type));
load(matfile0_name);
matfile1_name = fullfile(correction_dir, sprintf('scaled_auto_volume_%s.mat', ...
    correction_type));
load(matfile1_name);
if use_equalized_nano
    matfile0_name = fullfile(correction_dir,'equalized_volume.mat');
    load(matfile0_name);
    nanoVol = equalized_volume;
    clear equalized_volume
else
end
matfile2_name = fullfile(correction_dir, sprintf('artifact_mask_volume_%s.mat', ...
    correction_type));
if not(exist(matfile2_name))
    artifact_mask_vol = false(size(bg_mask_vol),'like',bg_mask_vol);
else
    load(matfile2_name);
end

% combine the background and artifact masks
[H, W, Z] = size(nanoVol);
invalid_mask = single(or(artifact_mask_vol,bg_mask_vol));

% load the original DAPI channel, the one registered
matfile3_name = fullfile(before_correction_dir, sprintf('chan01_DAPI.tiff'));
dapiVol = single(loadVolume({matfile3_name}, 1));

% load sliceinfo
sliceinfo_name = fullfile(mouse_dir, sprintf('sliceinfo.mat'));
load(sliceinfo_name);

% the channels and the metadata to save
slicevol_new = uint16(permute(cat(4, dapiVol, nanoVol, scaledautoVol, correctedVol, ...
    invalid_mask), [1 2 4 3]));
sliceinfo_new = sliceinfo;
sliceinfo_new.channames = {'DAPI','NANO','AUTO','DIFF','MASK'};
sliceinfo_new.slicevol = processed_dir;
sliceinfo_new.procpath = mouse_dir;
sliceinfo_new.volorder = volorder_dir;
sliceinfo_new.slicevolfin = aligned_dir;
sliceinfo_new.backvalues = recompute_backvalues(slicevol_new);

% save the processed data as LightSuite channels
saveLargeSliceVolume(slicevol_new, sliceinfo_new.channames, sliceinfo_new.slicevol);

% the sliceinfo with those channels
sliceinfo = sliceinfo_new;
end

function align_slices(sliceinfo, base_dir, mouse_name, mouse_dir, atlas)
% Mode 'align', second half: re-read the mouse's local settings, then align
% the slices and fit the atlas rigidly (writes regopts.mat).

% re-read local_settings.txt (the mouse folder first), so an edit to it takes
% effect: sliceinfo.mat holds the values run_extract_and_center wrote
settings_name = fullfile(base_dir, mouse_name, 'local_settings.txt');
if ~exist(settings_name, 'file')
    settings_name = fullfile(mouse_dir, 'local_settings.txt');
end
mouse_settings = parseSettingsFile(settings_name);
sliceinfo = copyStructBtoA(sliceinfo, mouse_settings);

% the atlas in force wins over the file, so the two cannot disagree about the
% resolution or the crop whatever a stale settings file says
sliceinfo.px_atlas    = atlas.res_um;
sliceinfo.atlasaplims = atlas.default_aplims;

fprintf('  settings: px_atlas %g um, atlasaplims %s, slicethickness %g\n', ...
    sliceinfo.px_atlas, mat2str(sliceinfo.atlasaplims), sliceinfo.slicethickness);

alignedvol = alignSliceVolume(sliceinfo.slicevol, sliceinfo);

% next, the cutting angle by eye (run_mode 'angle', optional but every adult had
% it), then 'annotate'
fprintf(['  aligned. regopts.mat and volume_for_inspection.tiff are written, so this\n' ...
         '  mouse is ready for run_mode = ''annotate''.\n']);
end
