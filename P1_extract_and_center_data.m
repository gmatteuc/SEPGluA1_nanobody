close all
clear all
clc

% /// Pipeline script #1: extracts data from raw .czi files and save centered volumes for further processing  /// 

% Set mice and mousetypes
mice = {'MG691_Gria1', 'MG692_Gria1', 'MG693_Gria1', 'MG736_Gria1', 'MG737_Gria1', 'CGF027_Gria1', 'CGF028_Gria1', 'CGF033_Gria1', 'CGF034_Gria1', 'CGF035_Gria1','MG705_Gria1', 'MG706_Gria1', 'MG709_Gria1', 'MG716_Gria1', 'MG718_Gria1', 'MG725_Gria1', 'MG727_Gria1'};
mousetypes = {'rws','rws','rws','rws','rws','naive','naive','naive','naive','naive','behavior','behavior','behavior','behavior','behavior','behavior','behavior'};
base_root = 'D:\sep_histology';

for mouse_idx=1:numel(mice)

    %% Initialize extraction via Lightsuite

    % Get current mouse metadata
    mousename      = mice{mouse_idx};
    dp = fullfile(base_root, mousetypes{mouse_idx}, sprintf('*%s*', mousename));
    dp = dir(dp);
    dp = fullfile(dp.folder, dp.name);
    sliceinfo = parseSettingsFile(fullfile(dp, 'local_settings.txt'));
    sliceinfo.mousename   = mousename;
    filelistcheck         = dir(fullfile(dp, '*.czi'));
    filepaths             = fullfile({filelistcheck(:).folder}', {filelistcheck(:).name}');
    sliceinfo.filepaths   = filepaths;
    sliceinfo             = getSliceInfo(sliceinfo);

    %% Generate the slice volume (auto)

    % Generate centered volume from raw data
    slicevol = generateSliceVolume(sliceinfo, sliceinfo.regchan);

    %% Reorder, flip and discard slices if needed (manual)

    % Edit slice ordering annotation if needed
    generateReordedVolume(sliceinfo);
    % SliceOrderEditor(sliceinfo.volorder)

end