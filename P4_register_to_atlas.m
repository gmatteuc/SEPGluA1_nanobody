clear all
close all
clc

% /// Pipeline script #4: feed residual correction anlaysis outputs in orgiginal atlas registration pipeline ///

%% User-defined parameters

mice = {'MG691_Gria1', 'MG692_Gria1', 'MG693_Gria1', 'MG736_Gria1', 'MG737_Gria1', 'CGF027_Gria1', 'CGF028_Gria1', 'CGF033_Gria1', 'CGF034_Gria1', 'CGF035_Gria1','MG705_Gria1', 'MG706_Gria1', 'MG709_Gria1', 'MG716_Gria1', 'MG718_Gria1', 'MG725_Gria1', 'MG727_Gria1'};
mousetypes = {'rws','rws','rws','rws','rws','naive','naive','naive','naive','naive','behavior','behavior','behavior','behavior','behavior','behavior','behavior'};

% Choose correction type
correction_type = 'slicewise';

% Set if to use equalized nano volumes
use_equalized_nano = 1;

%% Add paths

allenDir = 'D:\sep_histology\data\atlas';
lightsuiteDir = 'D:\sep_histology\code\LightSuite-main';
yamlDir = 'D:\sep_histology\code\yamlmatlab';
elastixDir = 'D:\sep_histology\code\matlab_elastix-master';
addpath(allenDir)
addpath(genpath(lightsuiteDir))
addpath(genpath(yamlDir))
addpath(genpath(elastixDir))

%% Loop over mice

for mouse_idx = 6:numel(mice) %[11,12,13,15,17]% [6,9,10,14,16]%1:numel(mice)

    %% Fetch data and bridge preprocessing results to the original Lightsuite registartion pipeline

    % Get current mouse name and type
    mouse_name = mice{mouse_idx};
    mouse_type = mousetypes{mouse_idx};

    % Get dirs
    base_dir = ['D:\sep_histology\data\', mouse_type, '\'];
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