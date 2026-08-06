clear all
close all
clc

% /// Pipeline script #6: perform final normalization across mice in each group and plot averages /// 

%%  Set user-defined parameters

% Set mice and mousetypes
mice = {'MG691_Gria1', 'MG692_Gria1', 'MG693_Gria1', 'MG736_Gria1', 'MG737_Gria1', 'CGF027_Gria1', 'CGF028_Gria1', 'CGF033_Gria1', 'CGF034_Gria1', 'CGF035_Gria1','MG705_Gria1', 'MG706_Gria1', 'MG709_Gria1', 'MG716_Gria1', 'MG718_Gria1', 'MG725_Gria1', 'MG727_Gria1'};
mousetypes = {'rws','rws','rws','rws','rws','naive','naive','naive','naive','naive','behavior','behavior','behavior','behavior','behavior','behavior','behavior'};
mousetypes_list = {'rws','naive','behavior'};

% Set list of selected mice
selected_mice_idx_list{1} = 1:5; % rws
selected_mice_idx_list{2} = 1:5; % naive
selected_mice_idx_list{3} = [1,3,4,5]; % behavior (excluding 3 bad mice)

% Choose correction type
correction_type = 'slicewise';

% Set plotting and saving
plot_sigle_mouse_videos = false;
plot_diagnostic_plots_normalized = true;
plot_diagnostic_plots_raw = true;
save_scaled_diff_volume = true;
plot_average_video_normalized = true;
plot_average_video_raw = true;

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

%% Allen atlas setup

AllenFile = [allenDir,'\annotation_10.nii.gz'];
AllenVol = niftiread(AllenFile);
limits = [180 1079];
AllenCrop = AllenVol(limits(1):limits(2),:,:);
brainMask = AllenCrop > 0;

%% Plot videos per mouse type

for mousetype_idx = 1% 1:numel(mousetypes_list)
    tic

    % --- Fetch single mouse data for current type ---

    % Get current mouse type and select mice
    current_mouse_type = mousetypes_list{mousetype_idx};
    current_mice_idx = strcmp(mousetypes, current_mouse_type);
    current_mice = mice(current_mice_idx);
    num_current = numel(current_mice);
    % Get base_dirs
    base_dir = ['D:\sep_histology\data\', current_mouse_type, '\'];
    global_diagnostics_dir = ['D:\sep_histology\data\', current_mouse_type, '\global_diagnostics'];
    if not(exist(global_diagnostics_dir)) %#ok<EXIST>
        mkdir(global_diagnostics_dir)
    end
    % Load needed 4D data for individual mouse videos
    load(fullfile(base_dir, 'diff_4d_new.mat'));
    load(fullfile(base_dir, 'mask_4d.mat'));
    if plot_diagnostic_plots_raw
        load(fullfile(base_dir, 'nano_4d.mat')); %#ok<UNRCH>
        load(fullfile(base_dir, 'auto_4d.mat'));
    end
    % Set common parameters for plotting
    clim_to_use = 1;
    z_label_to_use = 'Relative difference (Nano - Auto)';

    % --- Plot video for each individual mouse from the 4D matrix (diff)
    if plot_sigle_mouse_videos

        for i = 1:num_current
            tic
            mouse_name = current_mice{i};
            diff_ind = diff_4d_new(:,:,:,i);
            mask_ind = not(mask_4d(:,:,:,i));
            write_diff_video(diff_ind, mask_ind, AllenCrop, brainMask, global_diagnostics_dir, ['diff_video_' mouse_name '.mp4'], clim_to_use, z_label_to_use);
            clear diff_ind mask_ind
            toc
        end

    end


    if plot_diagnostic_plots_raw

        % --- Plot example frame from single mouse (diff)
        for i = 1:num_current
            mouse_name = current_mice{i};
            diff_ind = diff_4d_new(:,:,:,i);
            mask_ind = not(mask_4d(:,:,:,i));
            atlas_vol = AllenCrop;
            brain_mask = brainMask;
            data_vol = diff_ind;
            mask_vol = mask_ind;
            clim_value = clim_to_use;
            cb_label = z_label_to_use;
            slice_num = 500;
            titlestring = ['Example - Diff - ',strrep(mouse_name,'_',' ')];
            fh = plot_diff_slice(data_vol, mask_vol, atlas_vol, brain_mask, slice_num, clim_value, cb_label, titlestring);
            saveas(fh, fullfile(global_diagnostics_dir, [mouse_name '_slice_' num2str(slice_num) '_diff.fig']));
            exportgraphics(fh, fullfile(global_diagnostics_dir, [mouse_name '_slice_' num2str(slice_num) '_diff.png']), 'Resolution', 300);
            close(fh);
        end

        % --- Plot example frame from single mouse (nano)
        for i = 1:num_current
            mouse_name = current_mice{i};
            nano_ind = nano_4d(:,:,:,i);
            mask_ind = not(mask_4d(:,:,:,i));
            atlas_vol = AllenCrop;
            brain_mask = brainMask;
            data_vol = nano_ind;
            mask_vol = mask_ind;
            cb_label = 'Absolute Intensity (Nano)';
            slice_num = 500;
            titlestring = ['Example - Nano - ',strrep(mouse_name,'_',' ')];
            min_clim = quantile(data_vol(:), 0.01);
            max_clim = quantile(data_vol(:), 0.99);
            clim_values = [min_clim,max_clim]';
            fh = plot_abs_slice(data_vol, mask_vol, atlas_vol, brain_mask, slice_num, clim_values, cb_label, titlestring);
            saveas(fh, fullfile(global_diagnostics_dir, [mouse_name '_slice_' num2str(slice_num) '_nano.fig']));
            exportgraphics(fh, fullfile(global_diagnostics_dir, [mouse_name '_slice_' num2str(slice_num) '_nano.png']), 'Resolution', 300);
            close(fh);
        end

        % --- Plot example frame from single mouse (auto)
        for i = 1:num_current
            mouse_name = current_mice{i};
            auto_ind = auto_4d(:,:,:,i);
            nano_ind = nano_4d(:,:,:,i);
            mask_ind = not(mask_4d(:,:,:,i));
            atlas_vol = AllenCrop;
            brain_mask = brainMask;
            data_vol = auto_ind;
            mask_vol = mask_ind;
            cb_label = 'Absolute Intensity (Auto)';
            slice_num = 500;
            titlestring = ['Example - Auto - ',strrep(mouse_name,'_',' ')];
            min_clim = quantile(nano_ind(:), 0.01);
            max_clim = quantile(nano_ind(:), 0.99);
            clim_values = [min_clim,max_clim]';
            fh = plot_abs_slice(data_vol, mask_vol, atlas_vol, brain_mask, slice_num, clim_values, cb_label, titlestring);
            saveas(fh, fullfile(global_diagnostics_dir, [mouse_name '_slice_' num2str(slice_num) '_auto_reg.fig']));
            exportgraphics(fh, fullfile(global_diagnostics_dir, [mouse_name '_slice_' num2str(slice_num) '_auto_reg.png']), 'Resolution', 300);
            close(fh);
        end

    end

    % --- Normalize diff volumes

    normalization_method = 'atlas_region';  % 'global_quantile' or 'atlas_region'
    norm_factors = NaN(1, num_current);
    switch normalization_method
        case 'global_quantile'

            % --- Global Quantile ---
            fprintf('Computing normalization factors based on Global Quantile (0.99999)...\n');
            target_quantile = 0.99999;
            for i = 1:num_current
                diff_ind = diff_4d_new(:,:,:,i);
                mask_ind = not(mask_4d(:,:,:,i));
                valid_mask = mask_ind & brainMask;
                values = diff_ind(valid_mask);
                values(isinf(values)) = [];
                norm_factors(i) = quantile(values(:), target_quantile);
            end

        case 'atlas_region'

            % --- Atlas based ---
            fprintf('Computing normalization factors based on Atlas Regions ...\n');
            target_regions = {'Striatum', 'Hippocampal formation'};
            norm_mask_3d = get_allen_region_mask(allenDir, AllenCrop, target_regions, brainMask);
            for i = 1:num_current
                diff_ind = diff_4d_new(:,:,:,i);
                mask_ind = not(mask_4d(:,:,:,i));
                valid_roi_mask = norm_mask_3d & mask_ind;
                roi_values = diff_ind(valid_roi_mask);
                roi_values(isinf(roi_values) | isnan(roi_values)) = [];
                norm_factors(i) = quantile(roi_values,0.99); 
            end

    end

    % --- Apply Normalization ---

    scaled_diff_4d_new = diff_4d_new;
    for i = 1:num_current
        scale_factor = norm_factors(i);
        if scale_factor == 0 || isnan(scale_factor)
            scale_factor = 1;
        end
        scaled_diff_4d_new(:,:,:,i) = diff_4d_new(:,:,:,i) ./ scale_factor;
        fprintf('Mouse: %s | Method: %s | Factor: %.4f\n', ...
            current_mice{i}, normalization_method, scale_factor);
    end

    % Save scaled volumes
    if save_scaled_diff_volume
        current_selected_mice_idx = selected_mice_idx_list{mousetype_idx};
        current_mouse_type = mousetypes_list{mousetype_idx};
        current_mice = mice(current_mice_idx);
        save(fullfile(base_dir, 'scaled_diff_4d_new.mat'),...
            'scaled_diff_4d_new',...
            'current_selected_mice_idx',...
            'current_mouse_type',...
            'current_mice',...
            '-v7.3');
    end

    % --- Plot example frame from single mouse (normalized)
    if plot_diagnostic_plots_normalized
        for i = 1:num_current
            mouse_name = current_mice{i};
            diff_ind = scaled_diff_4d_new(:,:,:,i);
            mask_ind = not(mask_4d(:,:,:,i));
            atlas_vol = AllenCrop;
            brain_mask = brainMask;
            data_vol = diff_ind;
            mask_vol = mask_ind;
            clim_value = 1;  % Now symmetric [-1, 1]
            cb_label = 'Normalized Relative difference (Nano - Auto)';
            slice_num = 500;
            titlestring = ['Normalized Example - Diff - ', strrep(mouse_name, '_', ' ')];
            fh = plot_diff_slice(data_vol, mask_vol, atlas_vol, brain_mask, slice_num, clim_value, cb_label, titlestring);
            saveas(fh, fullfile(global_diagnostics_dir, [mouse_name '_normalized_slice_' num2str(slice_num) '_diff.fig']));
            exportgraphics(fh, fullfile(global_diagnostics_dir, [mouse_name '_normalized_slice_' num2str(slice_num) '_diff.png']), 'Resolution', 300);
            close(fh);
        end
    end

    % --- Plot diff histograms
    if plot_diagnostic_plots_normalized
        edges = linspace(0, 5, 1001); 
        num = num_current;
        intensity = linspace(1, 0.2, num)';
        colors = [intensity, zeros(num,1), zeros(num,1)];
        h_handles = gobjects(num, 1);
        mouse_labels = cell(num, 1);
        fh = figure('visible', 'on', 'units', 'normalized', 'outerposition', [0 0 1 1]);
        hold on;
        for i = 1:num_current
            mouse_name = current_mice{i};
            diff_ind = diff_4d_new(:,:,:,i);
            mask_ind = not(mask_4d(:,:,:,i));
            data_vol = diff_ind;
            mask_vol = mask_ind;
            brain_mask = brainMask;
            valid_mask = mask_vol & brain_mask;
            values = data_vol(valid_mask);
            subsample_size = min(1e5, length(values));
            sub_idx = randsample(length(values), subsample_size);
            sub_values = values(sub_idx);
            h_handles(i) = histogram(sub_values, 'BinEdges', edges, 'Normalization', 'probability', ...
                'DisplayStyle', 'stairs', 'EdgeColor', colors(i,:), 'LineWidth', 2);
            mouse_labels{i} = strrep(mouse_name, '_', ' ');
        end
        hold off;
        xlabel('Value');
        ylabel('Probability');
        title(['Overlaid histograms - ' current_mouse_type]);
        legend(h_handles, mouse_labels, 'Location', 'best');
        grid on;
        xlim([0 5]);
        saveas(fh, fullfile(global_diagnostics_dir, ['overlaid_histograms_' current_mouse_type '.fig']));
        exportgraphics(fh, fullfile(global_diagnostics_dir, ['overlaid_histograms_' current_mouse_type '.png']), 'Resolution', 300);
        close(fh);
    end

    % Plot video for the average (normalized)
    if plot_average_video_normalized
        load(fullfile(base_dir, 'avg_mask.mat'));
        avg_diff_new_scaled = nanmean(scaled_diff_4d_new(:,:,:,selected_mice_idx_list{mousetype_idx}),4); %#ok<NANMEAN>
        write_diff_video(avg_diff_new_scaled, avg_mask, AllenCrop, brainMask, base_dir, ['avg_diff_scaled_video_' current_mouse_type '.mp4'], 0.033*clim_to_use, z_label_to_use);
    end

    % Plot video for the average (raw)
    if plot_average_video_raw
        load(fullfile(base_dir, 'avg_mask.mat'));
        avg_diff_new = nanmean(diff_4d_new(:,:,:,selected_mice_idx_list{mousetype_idx}),4); %#ok<NANMEAN>
        write_diff_video(avg_diff_new, avg_mask, AllenCrop, brainMask, base_dir, ['avg_diff_video_' current_mouse_type '.mp4'], clim_to_use, z_label_to_use);
    end

    % Clear 4D volumes to free memory before loading averages
    clear diff_4d_new mask_4d

    toc

end
