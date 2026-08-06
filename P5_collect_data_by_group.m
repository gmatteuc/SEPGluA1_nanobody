clear all
close all
clc

% /// Pipeline script #5: load output of registartion for all mice and assemble 4D volumes stording data from all mice of each group /// 

%% User-defined parameters

% Set mice and mousetypes
mice = {'MG691_Gria1', 'MG692_Gria1', 'MG693_Gria1', 'MG736_Gria1', 'MG737_Gria1', 'CGF027_Gria1', 'CGF028_Gria1', 'CGF033_Gria1', 'CGF034_Gria1', 'CGF035_Gria1','MG705_Gria1', 'MG706_Gria1', 'MG709_Gria1', 'MG716_Gria1', 'MG718_Gria1', 'MG725_Gria1', 'MG727_Gria1'};
mousetypes = {'rws','rws','rws','rws','rws','naive','naive','naive','naive','naive','behavior','behavior','behavior','behavior','behavior','behavior','behavior'};
mousetypes_list = {'rws','naive','behavior'};
% Choose correction type
correction_type = 'slicewise';

%% Add paths

% Define paths
allenDir = 'D:\sep_histology\data\atlas';
lightsuiteDir = 'D:\sep_histology\code\LightSuite-main';
yamlDir = 'D:\sep_histology\code\yamlmatlab';
elastixDir = 'D:\sep_histology\code\matlab_elastix-master';
% Add defined paths
addpath(allenDir)
addpath(genpath(lightsuiteDir))
addpath(genpath(yamlDir))
addpath(genpath(elastixDir))

%% Get atlas

% Load Allen atlas annotations
AllenFile = [allenDir,'\annotation_10.nii.gz'];
AllenVol = niftiread(AllenFile);
limits = [180 1079];
AllenCrop = AllenVol(limits(1):limits(2),:,:);
brainMask = AllenCrop > 0;

%% Process per mouse type

for mousetype_idx = 2:3%1:numel(mousetypes_list)

    % Get current mouse type and select mice
    current_mouse_type = mousetypes_list{mousetype_idx};
    current_mice_idx = strcmp(mousetypes, current_mouse_type);
    current_mice = mice(current_mice_idx);
    num_current = numel(current_mice);

    % Initialize cells for volumes of this type
    nanoVols_type = cell(1, num_current);
    autoVols_type = cell(1, num_current);
    maskVols_type = cell(1, num_current);

    for i = 1:num_current

        % Get current mouse
        mouse_name = current_mice{i};

        % Get dirs (base_dir is per type)
        base_dir = ['D:\sep_histology\data\', current_mouse_type, ''];
        mouse_dir = fullfile(base_dir, mouse_name, 'lightsuite');
        correction_dir = fullfile(base_dir, mouse_name, 'lightsuite', 'correction_output');
        before_correction_dir = fullfile(base_dir, mouse_name, 'lightsuite', 'volume_centered');
        processed_dir = fullfile(base_dir, mouse_name, 'lightsuite', 'volume_centered_processed');
        aligned_dir = fullfile(mouse_dir, 'volume_aligned');
        registered_dir = fullfile(mouse_dir, 'volume_registered');

        % Load needed data for current mouse
        file1_name = fullfile(registered_dir, sprintf('chan02_NANO.tiff'));
        file2_name = fullfile(registered_dir, sprintf('chan03_AUTO.tiff'));
        file4_name = fullfile(registered_dir, sprintf('chan05_MASK.tiff'));
        nanoVol = single(loadVolume({file1_name}, 1));
        autoVol = single(loadVolume({file2_name}, 1));
        maskVol = single(loadVolume({file4_name}, 1));

        % Store in cells
        nanoVols_type{i} = nanoVol;
        autoVols_type{i} = autoVol;
        maskVols_type{i} = maskVol;

    end

    % Concatenate volumes along the 4th dimension
    nano_4d = cat(4, nanoVols_type{:});
    auto_4d = cat(4, autoVols_type{:});
    mask_4d = cat(4, maskVols_type{:});

    % Clear temporary cells to free memory
    clear nanoVols_type autoVols_type maskVols_type

    % Compute averages along the 4th dimension
    avg_nano = nanmean(nano_4d, 4); %#ok<NANMEAN>
    avg_auto = nanmean(auto_4d, 4); %#ok<NANMEAN>
    sum_mask_4d = nansum(mask_4d, 4); %#ok<NANSUM>

    % Form derived mask for display
    min_num_contrib=3;
    avg_mask=sum_mask_4d<=min_num_contrib;

    % Compute derived diff volumes
    diff_4d_new = (nano_4d-auto_4d)./auto_4d;
    avg_diff_new = nanmean(diff_4d_new, 4); %#ok<NANMEAN>

    % Save 4D volumes and averages in base_dir (use -v7.3 for large matrices)
    save(fullfile(base_dir, 'nano_4d.mat'), 'nano_4d', '-v7.3');
    % save(fullfile(base_dir, 'auto_4d.mat'), 'auto_4d', '-v7.3');
    % save(fullfile(base_dir, 'mask_4d.mat'), 'mask_4d', '-v7.3');
    % save(fullfile(base_dir, 'diff_4d_new.mat'), 'diff_4d_new', '-v7.3');
    % save(fullfile(base_dir, 'avg_nano.mat'), 'avg_nano', '-v7.3');
    % save(fullfile(base_dir, 'avg_auto.mat'), 'avg_auto', '-v7.3');
    % save(fullfile(base_dir, 'avg_mask.mat'), 'avg_mask', '-v7.3');
    % save(fullfile(base_dir, 'sum_mask_4d.mat'), 'sum_mask_4d', '-v7.3');
    % save(fullfile(base_dir, 'avg_diff_new.mat'), 'avg_diff_new', '-v7.3');

    % Clear variables to free memory before next type
    clear nano_4d auto_4d mask_4d diff_4d_new avg_nano avg_auto avg_mask sum_mask_4d avg_diff_new

end