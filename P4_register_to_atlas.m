clear all
close all
clc

% /// Pipeline script #4: feed residual correction anlaysis outputs in orgiginal atlas registration pipeline ///
%
% This runs in three modes, because the manual control-point step sits in the
% middle of it and the two automatic halves have to be run either side:
%
%   run_mode = 'align'     (auto)   bridge the corrected volumes into LightSuite,
%                                   align the slices and fit the atlas rigidly.
%                                   Writes regopts.mat and volume_for_inspection.tiff.
%   run_mode = 'annotate'  (MANUAL) open the control-point GUI on one mouse.
%                                   Writes atlas2histology_tform.mat.
%   run_mode = 'register'  (auto)   elastix refinement and the registered volumes.
%                                   Picks up the control points if they exist.
%
% The adults were done this way, one mouse at a time, with the GUI lines
% uncommented by hand. Every one of them has control points on every slice, so
% a young brain registered without them is not being treated the same way --
% see the note on 'annotate' below.

%% User-defined parameters

% Where the project lives. Derived from the location of the code rather than
% written out, so the tree can be moved or copied to another drive as is.
paths = get_paths();

% Cohort selection (mice come from the shared registry get_cohort.m).
% Set mice_to_process to {} to process every mouse in groups_to_process.
groups_to_process = {'young'};                  % 'rws' | 'naive' | 'behavior' | 'young'
mice_to_process   = {'MG903_SepGluA_P20', 'MG897_SepGluA_P20'};   % {} = all mice in groups_to_process
                                                % P20 first: curated and the age Sami wants prioritised

% Which half of the script to run. 'annotate' takes one mouse at a time.
run_mode = 'align';                             % 'align' | 'annotate' | 'register'

% Reference atlas.
%
% DECIDED 2026-09-02: BOTH cohorts register to the adult Allen CCF. The point
% is that young and adult get identical treatment -- registering the two groups
% to two differently built templates would put a methodological difference
% exactly along the axis the whole comparison rests on.
%
% The age-matched alternative was built and measured before this was settled,
% and it did not earn its costs: region centroids sit 0.22 mm apart between the
% two atlases (2.4% of brain length, inside DeMBA's own 0.10-0.19 mm build
% accuracy), cortical thickness matches at ratio 1.03, and MG903 registered to
% the adult CCF better than a typical adult does. Meanwhile the adult template
% carries 1.4x the global CV and 1.6x the local smoothed gradient, which helps
% both the image-driven B-spline and the placing of control points against it.
% See registration_qc\atlas_region_comparison.png.
%
% 'demba_p20' remains built, label-remapped and crop-verified, so this is a
% one-line switch if it is ever wanted. Do not mix atlases within a run.
atlas_key = 'ccf';                              % 'ccf' | 'demba_p20'

% Choose correction type
correction_type = 'slicewise';

% Set if to use equalized nano volumes
use_equalized_nano = 1;

% Register without manual control points. Off, and it should stay off for
% anything that ends up in a figure: LightSuite relies on the control points to
% register well, so an image-only run is a diagnostic, not a result.
allow_image_only_registration = false;

%% Add paths

% get_atlas puts the chosen atlas dir on the path and takes the other one off,
% which matters because every atlas dir holds files with identical names and
% LightSuite finds them with which(). Do not addpath an atlas dir by hand.
atlas = get_atlas(atlas_key);
lightsuiteDir = paths.lightsuite;
yamlDir = paths.yaml;
elastixDir = paths.elastix;
addpath(genpath(lightsuiteDir))
addpath(genpath(yamlDir))
addpath(genpath(elastixDir))

%% Resolve cohort

get_cohort('verify');
if isempty(mice_to_process)
    cohort = get_cohort('groups', groups_to_process);
else
    cohort = get_cohort('names', mice_to_process);
end
% Catch a mistyped mode here rather than letting it fall through to 'align' and
% quietly redo an hour of bridging nobody asked for.
if ~ismember(run_mode, {'align', 'annotate', 'register'})
    error('P4: unknown run_mode ''%s'' (use ''align'', ''annotate'' or ''register'').', run_mode);
end

fprintf('P4: %d mouse/mice selected, mode ''%s'', atlas ''%s''.\n', ...
    numel(cohort), run_mode, atlas.key);

% Registering against the wrong atlas produces a perfectly plausible-looking
% result, so say out loud which one is in force whenever it is not the one the
% project settled on. Every cohort goes to 'ccf'; anything else is a deliberate
% experiment and should be treated as one.
if ~strcmp(atlas.key, 'ccf')
    warning(['P4: registering to ''%s'', not the adult CCF that both cohorts ' ...
             'use. Outputs from this run are NOT comparable with the rest of ' ...
             'the dataset -- keep them out of P5.'], atlas.key);
end

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
        error(['P4: %s has px_atlas = %s but atlas ''%s'' is %g um.\n' ...
               'Fix px_atlas (and atlasaplims) in\n  %s'], ...
               cohort(k).name, tok{1}, atlas.key, atlas.res_um, settings_name);
    end
end
fprintf('P4: atlas resolution agrees with local_settings for all selected mice.\n');

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
    if ismember(run_mode, {'annotate', 'register'})

        regopts_name = fullfile(mouse_dir, 'regopts.mat');
        if ~exist(regopts_name, 'file')
            error(['P4: no regopts.mat for %s:\n  %s\n' ...
                   'Run this script with run_mode = ''align'' for this mouse first.'], ...
                   mouse_name, regopts_name);
        end
        opts = load(regopts_name);

        switch run_mode

            case 'annotate'
                % One GUI at a time, or the control points get placed in the
                % wrong mouse's file.
                if numel(cohort) > 1
                    error('P4: run_mode ''annotate'' opens one GUI at a time; select a single mouse.');
                end
                tform_name = fullfile(mouse_dir, 'atlas2histology_tform.mat');
                if exist(tform_name, 'file')
                    fprintf('  NOTE: control points already exist and will be overwritten on save:\n    %s\n', tform_name);
                end
                fprintf('  opening the control-point GUI against atlas ''%s''.\n', atlas.key);
                fprintf('  place points on every slice, then SAVE and CLOSE, and re-run with run_mode = ''register''.\n');
                matchControlPointsInSlices(opts);

            case 'register'
                % LightSuite does not register well from image information
                % alone -- the manual control points are what make it work, and
                % all 17 adults have them on every slice. So a missing
                % atlas2histology_tform.mat is an error, not a fallback.
                % Registering without them would quietly produce a volume that
                % looks fine in the folder and is not comparable to the adults.
                tform_name = fullfile(mouse_dir, 'atlas2histology_tform.mat');
                if exist(tform_name, 'file')
                    S_cp  = load(tform_name, 'histology_control_points');
                    n_cp  = cellfun(@(c) size(c, 1), S_cp.histology_control_points);
                    fprintf('  control points: %d slices, %d-%d points each (median %g)\n', ...
                        numel(n_cp), min(n_cp), max(n_cp), median(n_cp));
                    if any(n_cp == 0)
                        fprintf('  NOTE: %d slice(s) carry no points; those fall back to image only.\n', ...
                            nnz(n_cp == 0));
                    end
                elseif allow_image_only_registration
                    fprintf(['  no control points, and allow_image_only_registration is true.\n' ...
                             '  Registering from image information alone -- diagnostic only,\n' ...
                             '  do not compare the result against the adults.\n']);
                else
                    error(['P4: no control points for %s:\n  %s\n' ...
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
    % stayed frozen at whatever P1 baked into sliceinfo.mat, which means editing
    % local_settings.txt had no effect at all. Re-read the file properly.
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

    % %% (manual) Determine cutting angle gui if you are not happy with the original estimation
    % opts = load(fullfile(sliceinfo.procpath, "regopts.mat"));
    % determineCuttingAngleGUI(opts)

    fprintf(['  aligned. regopts.mat and volume_for_inspection.tiff are written, so this\n' ...
             '  mouse is ready for run_mode = ''annotate''.\n']);

end