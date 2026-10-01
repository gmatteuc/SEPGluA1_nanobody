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

% The settings of run_register_to_atlas, under the names the code below uses
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

% get_atlas puts the chosen atlas dir on the path and takes the other one off,
% which matters because every atlas dir holds files with identical names and
% LightSuite finds them with which(). Do not addpath an atlas dir by hand.
% The toolboxes are on the path from sep_setup_paths.
atlas = get_atlas(atlas_key);

%% Resolve cohort

get_cohort('verify');
if isempty(mice_to_process)
    cohort = get_cohort('groups', groups_to_process);
else
    cohort = get_cohort('names', mice_to_process);
end
% Catch a mistyped mode here rather than letting it fall through to 'align' and
% quietly redo an hour of bridging nobody asked for.
if ~ismember(run_mode, {'align', 'angle', 'annotate', 'autoannotate', 'register'})
    error('run_register_to_atlas: unknown run_mode ''%s'' (use ''align'', ''angle'', ''annotate'', ''autoannotate'' or ''register'').', run_mode);
end

fprintf('run_register_to_atlas: %d mouse/mice selected, mode ''%s'', atlas ''%s''.\n', ...
    numel(cohort), run_mode, atlas.key);

% Registering against the wrong atlas produces a perfectly plausible-looking
% result, so check the age rather than trusting the operator. An age-matched
% atlas is valid for its own age and nothing else: a P36 brain does not belong
% on the P20 template any more than it belongs on the adult one. Adults carry
% age_days = NaN in the registry and go to 'ccf' (age_days 56).
for k = 1:numel(cohort)
    mouse_age = cohort(k).age_days;
    if isnan(mouse_age)
        mouse_age = 56;     % the adult cohorts are not dated individually
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

% The atlas resolution has to agree with what each mouse's local_settings.txt
% says, because px_atlas is what sets the AP scale of the reconstruction. A
% young brain left at px_atlas = 10 against a 20 um atlas is off by a factor 2.
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

% Aligning again rewrites the atlas block that control points, anchors and
% proposals are counted in, so an annotated brain would silently lose its
% annotation (same rule as 'angle' below). Checked for every selected mouse
% before any of them is aligned.
if strcmp(run_mode, 'align')
    annotation_files = {'atlas2histology_tform.mat', 'plane_anchors.mat', ...
                        'auto_atlas_planes.mat', 'auto_proposal_controlpoints.mat'};
    annotated = {};
    for k = 1:numel(cohort)
        found = annotation_files(cellfun(@(f) exist(fullfile(cohort(k).base_dir, ...
            'lightsuite', f), 'file') == 2, annotation_files));
        if ~isempty(found)
            annotated{end+1} = sprintf('  %s: %s', cohort(k).name, strjoin(found, ', ')); %#ok<SAGROW>
        end
    end
    if ~isempty(annotated)
        error(['run_register_to_atlas: these mice already have an annotation, which aligning again would ' ...
               'invalidate:\n%s\nMove those files aside first if you really want to redo ' ...
               'both.'], strjoin(annotated, '\n'));
    end
end

%% Loop over mice

for mouse_idx = 1:numel(cohort)

    %% Fetch data and bridge preprocessing results to the original Lightsuite registartion pipeline

    % Get current mouse name and type
    mouse_name = cohort(mouse_idx).name;
    mouse_type = cohort(mouse_idx).group;

    fprintf('\n=== %s (%s) ===\n', mouse_name, mouse_type);

    % Get dirs
    base_dir = fullfile(paths.data, mouse_type);
    mouse_dir = fullfile(base_dir, mouse_name, '\lightsuite');
    correction_dir = fullfile(base_dir, mouse_name, '\lightsuite', 'correction_output');
    before_correction_dir = fullfile(base_dir, mouse_name, '\lightsuite', 'volume_centered');
    processed_dir = fullfile(base_dir, mouse_name, '\lightsuite', 'volume_centered_processed');
    aligned_dir = fullfile(mouse_dir,  'volume_aligned');
    volorder_dir = fullfile(mouse_dir, 'volume_for_ordering.tiff');

    % 'annotate' and 'register' both work off what 'align' already wrote, so
    % they skip the expensive bridging below and go straight to their step.
    if ismember(run_mode, {'angle', 'annotate', 'autoannotate', 'register'})

        regopts_name = fullfile(mouse_dir, 'regopts.mat');
        if ~exist(regopts_name, 'file')
            error(['run_register_to_atlas: no regopts.mat for %s:\n  %s\n' ...
                   'Run this script with run_mode = ''align'' for this mouse first.'], ...
                   mouse_name, regopts_name);
        end
        opts = load(regopts_name);
        % regopts.procpath is written by 'align' as an absolute path, drive
        % letter included, so a mouse aligned on another machine (or the same
        % disk under a different letter) points the GUI and the registration
        % at a folder that does not exist here. The folder is always this
        % mouse's lightsuite directory, so say so from where the code sits.
        opts.procpath = mouse_dir;

        % Apply the AP margin (see atlas_extent_slices above). The GUI reads it
        % from opts and the registration re-reads regopts.mat from disk, so
        % the two must agree: write it back. A mouse that already has control
        % points keeps whatever margin it was annotated with, because the
        % saved atlas planes are counted from the front of that range --
        % the three P20 mice done before this parameter existed reopen and
        % re-register exactly as before, with no setting to remember.
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

        switch run_mode

            case 'angle'
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

            case 'annotate'
                % One GUI at a time, or the control points get placed in the
                % wrong mouse's file.
                if numel(cohort) > 1
                    error('run_register_to_atlas: run_mode ''annotate'' opens one GUI at a time; select a single mouse.');
                end
                tform_name = fullfile(mouse_dir, 'atlas2histology_tform.mat');
                if exist(tform_name, 'file')
                    fprintf('  NOTE: control points already exist and will be overwritten on save:\n    %s\n', tform_name);
                end
                % The proposal keys (r, t) call the image matcher through a
                % persistent Python worker; start it now so the first press is
                % fast. Optional: if Python is missing the GUI still opens and
                % r reports why. See landmark_refine/README.md.
                landmark_refine_worker('start');
                if exist(fullfile(mouse_dir, 'cutting_angle_data.mat'), 'file')
                    fprintf('  cutting angle: from cutting_angle_data.mat (set by hand).\n');
                else
                    fprintf('  cutting angle: the automatic rigid fit''s (no cutting_angle_data.mat; run_mode = ''angle'' to set it by eye).\n');
                end
                fprintf('  opening the control-point GUI against atlas ''%s''.\n', atlas.key);
                fprintf('  place points on every slice, then SAVE and CLOSE, and re-run with run_mode = ''register''.\n');
                matchControlPointsInSlices(opts);

            case 'autoannotate'
                % The automatic control points, proposed from the anchor planes
                % set in the GUI. Writes a PROPOSAL next to the annotation; the
                % GUI loads it for review, and only what is accepted there is
                % ever saved as atlas2histology_tform.mat.
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
                    error('run_register_to_atlas: the automatic annotation failed for %s:\n%s', mouse_name, out.message);
                end
                fprintf(['  done. Review it: run_mode = ''annotate'' (k accepts a slice, u re-proposes it),\n' ...
                         '  then run_mode = ''register''.\n']);

            case 'register'
                % LightSuite does not register well from image information
                % alone -- the manual control points are what make it work, and
                % all 17 adults have them on every slice. So a missing
                % atlas2histology_tform.mat is an error, not a fallback.
                % Registering without them would quietly produce a volume that
                % looks fine in the folder and is not comparable to the adults.
                tform_name = fullfile(mouse_dir, 'atlas2histology_tform.mat');
                if exist(tform_name, 'file')
                    S_cp  = load(tform_name, 'histology_control_points', 'atlas_control_points');
                    n_cp  = cellfun(@(c) size(c, 1), S_cp.histology_control_points);
                    n_at  = cellfun(@(c) size(c, 1), S_cp.atlas_control_points);
                    fprintf('  control points: %d slices, %d-%d points each (median %g)\n', ...
                        numel(n_cp), min(n_cp), max(n_cp), median(n_cp));
                    % Points are paired by row, so a slice whose two lists differ
                    % in length has no valid pairing at all. elastix would only
                    % find out when it reached that slice, an hour in. Say so now.
                    bad = find(n_cp(:) ~= n_at(:))';
                    if ~isempty(bad)
                        error(['run_register_to_atlas: slice(s) %s have different numbers of histology and atlas points.\n' ...
                               'Open run_mode = ''annotate'', go to each, and delete the unpaired point(s)\n' ...
                               '(edit mode, d) or press c and re-place them.'], mat2str(bad));
                    end
                    % registerSlicesToAtlas uses the points only when a slice has
                    % at least five; with fewer it silently registers that slice
                    % from image information alone.
                    if any(n_cp < 5)
                        fprintf('  NOTE: slice(s) %s have fewer than 5 points and will be registered from image only.\n', ...
                            mat2str(find(n_cp(:) < 5)'));
                    end
                    % A mismatched pair -- a point on the wrong structure, or a
                    % left/right swap -- drags the whole slice's affine and only
                    % shows up as a bad overlay after the run. Fit each slice
                    % robustly and name the pairs that stand far off it, by the
                    % numbers the GUI draws next to them.
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

        continue

    end

    % Load the saved artifact annotation and correction outputs
    matfile0_name = fullfile(correction_dir, sprintf('corrected_volume_%s.mat', correction_type)); %#ok<UNRCH>
    load(matfile0_name);
    matfile1_name = fullfile(correction_dir, sprintf('scaled_auto_volume_%s.mat', correction_type));
    load(matfile1_name);
    if use_equalized_nano
        matfile0_name = fullfile(correction_dir,'equalized_volume.mat');
        load(matfile0_name);
        nanoVol = equalized_volume;
        clear equalized_volume
    else
    end
    matfile2_name = fullfile(correction_dir, sprintf('artifact_mask_volume_%s.mat', correction_type));
    if not(exist(matfile2_name))
        artifact_mask_vol = false(size(bg_mask_vol),'like',bg_mask_vol);
    else
        load(matfile2_name);
    end

    % Combine background and artifact mask
    [H, W, Z] = size(nanoVol);
    invalid_mask = single(or(artifact_mask_vol,bg_mask_vol));

    % Load original dapi channel for registration
    matfile3_name = fullfile(before_correction_dir, sprintf('chan01_DAPI.tiff'));
    dapiVol = single(loadVolume({matfile3_name}, 1));

    % Load sliceinfo
    sliceinfo_name = fullfile(mouse_dir, sprintf('sliceinfo.mat'));
    load(sliceinfo_name);

    % Prepare data and metadata for re-saving
    slicevol_new = uint16(permute(cat(4, dapiVol, nanoVol, scaledautoVol, correctedVol, invalid_mask), [1 2 4 3]));
    sliceinfo_new = sliceinfo;
    sliceinfo_new.channames = {'DAPI','NANO','AUTO','DIFF','MASK'};
    sliceinfo_new.slicevol = processed_dir;
    sliceinfo_new.procpath = mouse_dir;
    sliceinfo_new.volorder = volorder_dir;
    sliceinfo_new.slicevolfin = aligned_dir;
    sliceinfo_new.backvalues = recompute_backvalues(slicevol_new);

    % Re-save processed data as new channels
    saveLargeSliceVolume(slicevol_new, sliceinfo_new.channames, sliceinfo_new.slicevol);

    %% (auto) Align slices and initialize registration

    % and slicevol channelsnames
    sliceinfo = sliceinfo_new;

    % This used to read copyStructBtoA(sliceinfo, settings), but `settings` is
    % never defined in this script -- it resolved to MATLAB's own builtin, so
    % the line quietly copied the matlab/database/parallel setting groups into
    % sliceinfo and refreshed nothing. The values actually used for alignment
    % stayed frozen at whatever run_extract_and_center baked into sliceinfo.mat,
    % which means editing local_settings.txt had no effect at all. Re-read the
    % file properly.
    settings_name = fullfile(base_dir, mouse_name, 'local_settings.txt');
    if ~exist(settings_name, 'file')
        settings_name = fullfile(mouse_dir, 'local_settings.txt');
    end
    mouse_settings = parseSettingsFile(settings_name);
    sliceinfo = copyStructBtoA(sliceinfo, mouse_settings);

    % The atlas in force wins over the file, so the two can never disagree
    % about resolution or crop no matter what a stale settings file says.
    sliceinfo.px_atlas    = atlas.res_um;
    sliceinfo.atlasaplims = atlas.default_aplims;

    fprintf('  settings: px_atlas %g um, atlasaplims %s, slicethickness %g\n', ...
        sliceinfo.px_atlas, mat2str(sliceinfo.atlasaplims), sliceinfo.slicethickness);

    alignedvol = alignSliceVolume(sliceinfo.slicevol, sliceinfo);

    % The cutting angle can be set by eye next: run_mode = 'angle' (optional,
    % but every adult had it), then 'annotate'.

    fprintf(['  aligned. regopts.mat and volume_for_inspection.tiff are written, so this\n' ...
             '  mouse is ready for run_mode = ''annotate''.\n']);

end

end


function report_suspect_pairs(acp, hcp)
% Per slice: affine from atlas to histology, reweighted a few times so that
% outliers stop pulling the fit toward themselves, then list any pair further
% than 3x the slice's median residual (and at least 30 px) from it.
%
% Two failures that test cannot see, both found on MG910:
%
%   Points spread over two atlas planes. The GUI pins the plane to the median
%   of the points already placed, but a slice can still end up with points from
%   two planes (accept a carried-over proposal, scroll, place more). The
%   registration then takes the plane of the FIRST point and silently treats
%   the rest as if they were on it, so part of the slice is annotated against
%   anatomy 0.2 mm away. Nothing downstream complains.
%
%   A slice that is uniformly wrong. The outlier test is relative to the
%   slice's own median, so when most pairs are bad none of them stands out --
%   MG910 slice 35 had a median residual of 93 px and was never flagged.
for k = 1:numel(acp)
    a = acp{k}; h = hcp{k};
    if size(a,1) < 5 || size(a,1) ~= size(h,1), continue, end
    planes = unique(a(:,1));
    if numel(planes) > 1
        n_per = arrayfun(@(pl) nnz(a(:,1) == pl), planes);
        fprintf(['  CHECK slice %d: points sit on %d different atlas planes (%s, %s points each). ' ...
                 'Registration uses the first point''s plane for all of them -- clear the slice and ' ...
                 'redo it on one plane.\n'], k, numel(planes), mat2str(planes'), mat2str(n_per'));
    end
    src = a(:, [3 2]); dst = h(:, [3 2]);
    w = ones(size(src,1), 1);
    for it = 1:8
        tf = fitgeotform2d(src(w > 0.999, :), dst(w > 0.999, :), 'affine');
        r  = vecnorm(tf.transformPointsForward(src) - dst, 2, 2);
        s  = median(r) + 1e-6;
        w  = double(r <= 2*s);                      % drop, refit, repeat
        if nnz(w) < 5, break, end
    end
    bad = find(r > max(3*median(r), 30))';
    if ~isempty(bad)
        fprintf('  CHECK slice %d: pair(s) %s sit %s px off the slice''s own affine (median %.0f). Wrong structure or a left/right swap?\n', ...
            k, mat2str(bad), mat2str(round(r(bad))'), median(r));
    end
    % A well annotated slice sits at 5-10 px. Past 20 the slice as a whole does
    % not fit an affine, which usually means the wrong atlas plane rather than
    % a few wrong points -- and elastix tends to fail outright on those.
    if median(r) > 20
        fprintf(['  CHECK slice %d: the WHOLE slice is off -- median residual %.0f px over %d pairs. ' ...
                 'Wrong atlas plane, or points placed on two planes? Clear and redo it.\n'], ...
                 k, median(r), numel(r));
    end
end
end
