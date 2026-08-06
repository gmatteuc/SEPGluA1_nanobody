clear all
close all
clc

% /// Pipeline script #3: display outputs of residual correction anlaysis in a GUI where user can annotate artifact for future removal /// 

%% User-defined parameters

mice = {'MG691_Gria1', 'MG692_Gria1', 'MG693_Gria1', 'MG736_Gria1', 'MG737_Gria1', 'CGF027_Gria1', 'CGF028_Gria1', 'CGF033_Gria1', 'CGF034_Gria1', 'CGF035_Gria1','MG705_Gria1', 'MG706_Gria1', 'MG709_Gria1', 'MG716_Gria1', 'MG718_Gria1', 'MG725_Gria1', 'MG727_Gria1'};
mousetypes = {'rws','rws','rws','rws','rws','naive','naive','naive','naive','naive','behavior','behavior','behavior','behavior','behavior','behavior','behavior'};

% Choose correction type
correction_type = 'slicewise';

%% Loop over mice

for mouse_idx = 3%:numel(mice)

    % Get current mouse name and type
    mouse_name = mice{mouse_idx};
    mouse_type = mousetypes{mouse_idx};

    % Get dirs
    base_dir = ['D:\sep_histology\data\', mouse_type, '\'];
    output_dir = fullfile(base_dir, mouse_name, 'correction_output');
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    % Load the saved matfile
    matfile_name = fullfile(output_dir, sprintf('scaled_auto_volume_%s.mat', correction_type));
    if ~exist(matfile_name, 'file')
        fprintf('File not found for %s: %s\n', mouse_name, matfile_name);
        continue;
    end
    load(matfile_name);  % Loads scaledautoVol, nanoVol, bg_mask_vol, slice_data, average_slope, average_intercept, correction_type
    [H, W, Z] = size(nanoVol);

    % Invoke the artifact annotation GUI
    ArtifactAnnotator(scaledautoVol, nanoVol, bg_mask_vol, slice_data, mouse_name, output_dir, correction_type);

    % Clear per-mouse variables to save memory
    clear scaledautoVol nanoVol bg_mask_vol slice_data average_slope average_intercept

end