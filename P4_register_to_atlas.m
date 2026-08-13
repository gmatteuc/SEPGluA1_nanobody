clear all
close all
clc

% /// Pipeline script #4: feed residual correction anlaysis outputs in orgiginal atlas registration pipeline ///

%% User-defined parameters

% Where the project lives. Derived from the location of the code rather than
% written out, so the tree can be moved or copied to another drive as is.
paths = get_paths();

% Cohort selection (mice come from the shared registry get_cohort.m).
% Set mice_to_process to {} to process every mouse in groups_to_process.
groups_to_process = {'young'};                  % 'rws' | 'naive' | 'behavior' | 'young'
mice_to_process   = {'MG903_SepGluA_P20'};      % {} = all mice in groups_to_process

% Reference atlas
atlas_key = 'ccf';

% Choose correction type
correction_type = 'slicewise';

% Set if to use equalized nano volumes
use_equalized_nano = 1;

%% Add paths

atlas = get_atlas(atlas_key);
allenDir = atlas.dir;
lightsuiteDir = paths.lightsuite;
yamlDir = paths.yaml;
elastixDir = paths.elastix;
addpath(allenDir)
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
fprintf('P4: %d mouse/mice selected.\n', numel(cohort));

%% Loop over mice

for mouse_idx = 1:numel(cohort)

    %% Fetch data and bridge preprocessing results to the original Lightsuite registartion pipeline

    % Get current mouse name and type
    mouse_name = cohort(mouse_idx).name;
    mouse_type = cohort(mouse_idx).group;

    % Get dirs
    base_dir = fullfile(paths.data, mouse_type);
    mouse_dir = fullfile(base_dir, mouse_name, '\lightsuite');
    correction_dir = fullfile(base_dir, mouse_name, '\lightsuite', 'correction_output');
    before_correction_dir = fullfile(base_dir, mouse_name, '\lightsuite', 'volume_centered');
    processed_dir = fullfile(base_dir, mouse_name, '\lightsuite', 'volume_centered_processed');
    aligned_dir = fullfile(mouse_dir,  'volume_aligned');
    volorder_dir = fullfile(mouse_dir, 'volume_for_ordering.tiff');

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
    sliceinfo          = sliceinfo_new;
    sliceinfo          = copyStructBtoA(sliceinfo, settings);
    alignedvol         = alignSliceVolume(sliceinfo.slicevol, sliceinfo);

    % %% (manual) Determine cutting angle gui if you are not happy with the original estimation
    % opts = load(fullfile(sliceinfo.procpath, "regopts.mat"));
    % determineCuttingAngleGUI(opts)

    %% (manual) Match control points to determine cutting angle and gaps
    % opts = load(fullfile(sliceinfo_new.procpath, "regopts.mat"));
    % matchControlPointsInSlices(opts)

    %% (auto) Refine registation with control points and elastix
    opts            = load(fullfile(sliceinfo_new.procpath, "regopts.mat"));
    transformparams = registerSlicesToAtlas(opts); %#ok<NASGU>

    %% (auto) Apply registration to all color channels to generate registered volumes
    transformparams = load(fullfile(sliceinfo_new.procpath, "transform_params.mat"));
    generateRegisteredSliceVolume(sliceinfo_new, transformparams);

end