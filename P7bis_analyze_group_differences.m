clear all
close all
clc

% /// Pipeline script #7: analyze hemispheric differences and compare groups using normalized signal channels (nano or auto) ///

%%  Set user-defined parameters

% Set mice and mousetypes
mice = {'MG691_Gria1', 'MG692_Gria1', 'MG693_Gria1', 'MG736_Gria1', 'MG737_Gria1',...
    'CGF027_Gria1', 'CGF028_Gria1', 'CGF033_Gria1', 'CGF034_Gria1', 'CGF035_Gria1',...
    'MG705_Gria1', 'MG706_Gria1', 'MG709_Gria1', 'MG716_Gria1', 'MG718_Gria1', 'MG725_Gria1', 'MG727_Gria1'};
mousetypes = {'rws','rws','rws','rws','rws', ...
    'naive','naive','naive','naive','naive', ...
    'behavior','behavior','behavior','behavior','behavior','behavior','behavior'};
mousetypes_list = {'rws','naive','behavior'};

% Define your two groups
ctrl_type = 'naive';
exp_type  = 'behavior';

% Selection Filters
selected_mice_idx_list{1} = 1:5;        % rws
selected_mice_idx_list{2} = 1:5;        % naive
selected_mice_idx_list{3} = [1,3,4,5];    % behavior

% Get mousenames
ctrl_group_idx = find(strcmp(mousetypes_list, ctrl_type));
all_ctrl_indices = find(strcmp(mousetypes, ctrl_type));
final_ctrl_indices = all_ctrl_indices(selected_mice_idx_list{ctrl_group_idx});
ctrl_mousenames = mice(final_ctrl_indices);
exp_group_idx = find(strcmp(mousetypes_list, exp_type));
all_exp_indices = find(strcmp(mousetypes, exp_type));
final_exp_indices = all_exp_indices(selected_mice_idx_list{exp_group_idx});
exp_mousenames = mice(final_exp_indices);

% Set whether to generate difference videos
generate_diff_videos = true;
generate_individual_diff_videos = true;
generate_t_scored_videos = true;
generate_surprise_videos = true;
generate_rolling_videos = true;
generate_signed_diff_videos = true;

perform_area_based_analysis_fine= false;
perform_area_based_analysis_coarse = false;

% Smoothing Settings
apply_smoothing = true;
smooth_sigma = 5.0;

% Channel to analyze: 'nano' (default, surface GluA1) or 'auto' (autofluorescence control).
% Loads <channel>_4d_normalized.mat from each cohort folder. The channel is rolled
% into comp_tag so nano and auto outputs go to separate comparisons/ subfolders.
channel = 'nano';

% Comparison tag for filenames (includes channel so nano/auto runs don't collide)
comp_tag = [ctrl_type '_vs_' exp_type '_' channel];

% Base directory (common part)
base_root = 'D:\sep_histology\data\';

% Construct full paths
ctrl_dir = fullfile(base_root, ctrl_type);
exp_dir  = fullfile(base_root, exp_type);

% Define a specific directory for comparison results to avoid clutter
comp_out_dir = fullfile(base_root, 'comparisons', comp_tag);
if ~exist(comp_out_dir, 'dir'), mkdir(comp_out_dir); end

%% Allen atlas setup

% Prepare atlas
allenDir = 'D:\sep_histology\data\atlas';
addpath(allenDir);
AllenFile = fullfile(allenDir, 'annotation_10.nii.gz');
AllenVol = niftiread(AllenFile);
limits = [180 1079];
AllenCrop = AllenVol(limits(1):limits(2),:,:);
brainMask = AllenCrop > 0;
half_atlas = AllenCrop(:,:,1:end);
clear bg_L_c bg_R_c bg_L_e bg_R_e AllenVol

%% Load data for both groups

% Build channel-aware filenames + dynamic field name once
norm_var_name      = [channel '_4d_normalized'];                                  % field inside the .mat
norm_filename      = [channel '_4d_normalized.mat'];                              % e.g. 'nano_4d_normalized.mat'
bkgmask_filename   = [channel '_4d_normalized_bkgmask.mat'];

% Load control group data
fprintf('Loading Control group data (%s, channel=%s)...\n', ctrl_type, channel);
S_ctrl_vol  = load(fullfile(ctrl_dir, norm_filename),    norm_var_name);
S_ctrl_mask = load(fullfile(ctrl_dir, bkgmask_filename), 'recomputed_bkg_mask_4d');
data_4d_new_ctrl         = S_ctrl_vol.(norm_var_name); %(:,:,:,selected_mice_idx_list{ctrl_group_idx});
data_4d_new_ctrl_bkgmask = S_ctrl_mask.recomputed_bkg_mask_4d;
clear S_ctrl_vol S_ctrl_mask
fprintf('  Kept %d mice based on selection.\n', size(data_4d_new_ctrl, 4));

% Load experimental group data
fprintf('Loading Experimental group data (%s, channel=%s)...\n', exp_type, channel);
S_exp_vol  = load(fullfile(exp_dir, norm_filename),    norm_var_name);
S_exp_mask = load(fullfile(exp_dir, bkgmask_filename), 'recomputed_bkg_mask_4d');
if strcmp(exp_type,'rws')
    data_4d_new_exp         = S_exp_vol.(norm_var_name); %(:,:,:,selected_mice_idx_list{exp_group_idx});
    data_4d_new_exp_bkgmask = S_exp_mask.recomputed_bkg_mask_4d;
elseif strcmp(exp_type,'behavior')
    data_4d_new_exp         = S_exp_vol.(norm_var_name)(:,:,:,[1,2,3,4]); % NB: for behavior subselect 3 of the originally saved and selected 4 ('MG705_Gria1' 'MG716_Gria1' 'MG718_Gria1')
    data_4d_new_exp_bkgmask = S_exp_mask.recomputed_bkg_mask_4d(:,:,:,[1,2,3,4]);
end
clear S_exp_vol S_exp_mask
fprintf('  Kept %d mice based on selection.\n', size(data_4d_new_exp, 4));

%% Recompute 4D masks for entire registered volume

% Mask computation for control
fprintf('Recomputing background masks controls (%s)...\n', ctrl_type);
med_data_4d_ctrl = nan(size(data_4d_new_ctrl,1), size(data_4d_new_ctrl,4));
total_slices = size(data_4d_new_ctrl, 1);
for iii = 1:size(data_4d_new_ctrl,4)
    fprintf('  Processing Mouse %d ...\n', iii);
    for slice_idx_loop = 1:total_slices
        img_data = squeeze(data_4d_new_ctrl(slice_idx_loop,:,:,iii));
        bg_mask = squeeze(data_4d_new_ctrl_bkgmask(slice_idx_loop,:,:,iii));
        med_data_4d_ctrl(slice_idx_loop,iii) = nanmean(img_data(~bg_mask));
    end
end
% Just get the mask as recomputed one
recomputed_bkg_mask_4d_ctrl = data_4d_new_ctrl_bkgmask;
clear data_4d_new_ctrl_bkgmask

% Mask computation for experimental
fprintf('Recomputing background masks experimentals (%s)...\n', exp_type);
med_data_4d_exp = nan(size(data_4d_new_exp,1), size(data_4d_new_exp,4));
total_slices = size(data_4d_new_exp, 1);
for iii = 1:size(data_4d_new_exp,4)
    fprintf('  Processing Mouse %d ...\n', iii);
    for slice_idx_loop = 1:total_slices
        img_data = squeeze(data_4d_new_exp(slice_idx_loop,:,:,iii));
        bg_mask = squeeze(data_4d_new_exp_bkgmask(slice_idx_loop,:,:,iii));
        med_data_4d_exp(slice_idx_loop,iii) = nanmean(img_data(~bg_mask));
    end
end
% Just get the mask as recomputed one
recomputed_bkg_mask_4d_exp = data_4d_new_exp_bkgmask;
clear data_4d_new_exp_bkgmask

%% Apply NaN-tolerant Spatial Smoothing

if apply_smoothing
    fprintf('Applying NaN-Robust 3D Gaussian Smoothing (Sigma = %.1f)...\n', smooth_sigma);

    % Define helper function for NaN tolerant smoothing with imgaussfilt3
    nan_smooth_3d = @(v, sig) imgaussfilt3(fillmissing(v, 'constant', 0), sig) ./ ...
        imgaussfilt3(double(~isnan(v)), sig);

    % --- Control Group ---
    fprintf('  Processing Control group...\n');
    for i = 1:size(data_4d_new_ctrl, 4)
        tic
        vol = data_4d_new_ctrl(:,:,:,i);
        bg_mask = recomputed_bkg_mask_4d_ctrl(:,:,:,i);
        is_valid_tissue = brainMask & ~bg_mask;
        vol(~is_valid_tissue) = NaN;
        vol_smoothed = nan_smooth_3d(vol, smooth_sigma);
        vol_smoothed(~is_valid_tissue) = 0;
        vol_smoothed(isnan(vol_smoothed)) = 0;
        data_4d_new_ctrl(:,:,:,i) = vol_smoothed;
        toc
    end

    % --- Experimental Group ---
    fprintf('  Processing Experimental group...\n');
    for i = 1:size(data_4d_new_exp, 4)
        tic
        vol = data_4d_new_exp(:,:,:,i);
        bg_mask = recomputed_bkg_mask_4d_exp(:,:,:,i);
        is_valid_tissue = brainMask & ~bg_mask;
        vol(~is_valid_tissue) = NaN;
        vol_smoothed = nan_smooth_3d(vol, smooth_sigma);
        vol_smoothed(~is_valid_tissue) = 0;
        vol_smoothed(isnan(vol_smoothed)) = 0;
        data_4d_new_exp(:,:,:,i) = vol_smoothed;
        toc
    end

    fprintf('  Smoothing complete.\n');
end

%% Compute final normalization

% Set region of interest for calculating the fit
interest_region = 200:700;

% Calculate Group Mean Profiles
mean_profile_ctrl = nanmean(med_data_4d_ctrl, 2);
mean_profile_exp  = nanmean(med_data_4d_exp, 2);
% Extract Data in Interest Region for Fitting
y_target = mean_profile_ctrl(interest_region);
x_source = mean_profile_exp(interest_region);
% Remove NaNs if any exist in the mean profile
valid_idx = ~isnan(x_source) & ~isnan(y_target);
x_source = x_source(valid_idx);
y_target = y_target(valid_idx);
% We want to transform Exp to match Ctrl
p = polyfit(x_source, y_target, 1);
slope = p(1);
intercept = p(2);

fprintf('Alignment Parameters (Exp -> Ctrl): Slope = %.4f, Intercept = %.4f\n', slope, intercept);

% Control stays as is (Reference)
norm_ctrl = med_data_4d_ctrl;
% Experimental is transformed: New = Slope * Old + Intercept
norm_exp = (med_data_4d_exp .* slope) + intercept;

% Get numeric scaling factors
interest_region_bis = 300:500;
norm_ctrl_med_fact = nanmean(nanmean(med_data_4d_ctrl(interest_region_bis,:),1));
norm_exp_med_fact = nanmean(nanmean(norm_exp(interest_region_bis,:),1));
avg_med_fact = (norm_ctrl_med_fact + norm_exp_med_fact)./2;
norm_ctrl_med_fact = avg_med_fact;
norm_exp_med_fact = avg_med_fact;

% --- Plotting ---
fig_norm = figure('Visible', 'off', 'Name', ['Normalization_Profile_LR_' comp_tag], 'Position', [100, 100, 1200, 500], 'Color', 'w');
cmap_ctrl = [linspace(0.6, 0, 5)', linspace(0.7, 0.2, 5)', linspace(1, 0.4, 5)'];
cmap_exp  = [linspace(1, 0.5, 5)', linspace(0.6, 0.1, 5)', linspace(0.6, 0.1, 5)'];
slices = 1:size(med_data_4d_ctrl, 1);

% Subplot 1: Raw Data
subplot(1,2,1);
hold on; box on; grid on;
% Plot Individual Mice
for i = 1:size(med_data_4d_ctrl, 2), plot(slices, med_data_4d_ctrl(:,i), 'Color', cmap_ctrl(i,:), 'LineWidth', 1.2); end
for i = 1:size(med_data_4d_exp, 2),  plot(slices, med_data_4d_exp(:,i), 'Color', cmap_exp(i,:), 'LineWidth',1.2); end
% Plot Means & SEM
mean_c = nanmean(med_data_4d_ctrl, 2); sem_c = nanstd(med_data_4d_ctrl, [], 2) / sqrt(size(med_data_4d_ctrl,2));
mean_e = nanmean(med_data_4d_exp, 2);  sem_e = nanstd(med_data_4d_exp, [], 2) / sqrt(size(med_data_4d_exp,2));
fill([slices fliplr(slices)], [mean_c-sem_c; flipud(mean_c+sem_c)], [0 0.45 0.74], 'FaceAlpha',0.3, 'EdgeColor','none');
plot(slices, mean_c, 'Color', [0 0.2 0.5], 'LineWidth', 3.5);
fill([slices fliplr(slices)], [mean_e-sem_e; flipud(mean_e+sem_e)], [0.85 0.33 0.10], 'FaceAlpha',0.3, 'EdgeColor','none');
plot(slices, mean_e, 'Color', [0.64 0.08 0.18], 'LineWidth', 3.5);
% Formatting
title('Raw Nanobody intensity', 'FontSize',12);
xlabel('Coronal index', 'FontSize',12); ylabel('Intensity', 'FontSize',12);
xline(interest_region(1), 'k--'); xline(interest_region(end), 'k--'); % Show region
set(gca,'FontSize',12);

% Subplot 2: Normalized Data
subplot(1,2,2);
hold on; box on; grid on;
% Plot Individual Mice
for i = 1:size(norm_ctrl, 2), plot(slices, norm_ctrl(:,i), 'Color', cmap_ctrl(i,:), 'LineWidth', 1.2); end
for i = 1:size(norm_exp, 2),  plot(slices, norm_exp(:,i), 'Color', cmap_exp(i,:), 'LineWidth',1.2); end
% Plot Means & SEM
mean_nc = nanmean(norm_ctrl, 2); sem_nc = nanstd(norm_ctrl, [], 2) / sqrt(size(norm_ctrl,2));
mean_ne = nanmean(norm_exp, 2);  sem_ne = nanstd(norm_exp, [], 2) / sqrt(size(norm_exp,2));
fill([slices fliplr(slices)], [mean_nc-sem_nc; flipud(mean_nc+sem_nc)], [0 0.45 0.74], 'FaceAlpha',0.3, 'EdgeColor','none');
plot(slices, mean_nc, 'Color', [0 0.2 0.5], 'LineWidth', 3.5);
fill([slices fliplr(slices)], [mean_ne-sem_ne; flipud(mean_ne+sem_ne)], [0.85 0.33 0.10], 'FaceAlpha',0.3, 'EdgeColor','none');
plot(slices, mean_ne, 'Color', [0.64 0.08 0.18], 'LineWidth', 3.5);
% Formatting
title(sprintf('Linearly Aligned (Slope=%.2f, Int=%.0f)', slope, intercept), 'FontSize',12);
xlabel('Coronal index', 'FontSize',12); ylabel('Aligned Intensity', 'FontSize',12);
xline(interest_region(1), 'k--'); xline(interest_region(end), 'k--');
set(gca,'FontSize',12);
sgtitle(['Profile Alignment: ' ctrl_type ' (Ref) vs ' exp_type ' (Aligned)'], ...
    'FontSize',14, 'Color','k', 'Interpreter', 'none');

% Save plot
set(fig_norm, 'InvertHardcopy', 'off');
saveas(fig_norm, fullfile(comp_out_dir, ['Normalization_Profiles_LR_' comp_tag '.fig']));
exportgraphics(fig_norm, fullfile(comp_out_dir, ['Normalization_Profiles_LR_' comp_tag '.png']), 'Resolution', 300);
fprintf('Saved LR Normalization Profile plot to: %s\n', comp_out_dir);

%% Compute LR differences and group differences on normalized volumes

% Compute LR difference (Using Linear Regression params: Exp -> Ctrl)
[lr_diff_ctrl, lr_sum_ctrl] = compute_lr_stats(data_4d_new_ctrl./norm_ctrl_med_fact);
[lr_diff_exp,  lr_sum_exp]  = compute_lr_stats(((data_4d_new_exp .* slope) + intercept)./norm_exp_med_fact);
avg_lr_diff_ctrl = nanmean(abs(lr_diff_ctrl), 4); %#ok<*NANMEAN>
avg_lr_sum_ctrl  = nanmean(abs(lr_sum_ctrl),  4);
avg_lr_diff_exp  = nanmean(abs(lr_diff_exp),  4);
avg_lr_sum_exp   = nanmean(abs(lr_sum_exp),   4);
avg_lr_diff_groupdiff = avg_lr_diff_exp - avg_lr_diff_ctrl;
avg_lr_sum_groupdiff  = avg_lr_sum_exp  - avg_lr_sum_ctrl;

% Crop mask and atlas to actual data size
brainMask_cropped = brainMask(:, :, 1:size(avg_lr_diff_ctrl, 3));

% Setup Dimensions & Masks
n_half = size(lr_diff_ctrl, 3);
n_full = size(recomputed_bkg_mask_4d_ctrl, 3);
bg_L_c = recomputed_bkg_mask_4d_ctrl(:, :, 1:n_half, :);
bg_R_c = recomputed_bkg_mask_4d_ctrl(:, :, (n_full - n_half + 1):end, :);
mask_bg_ctrl = logical(bg_L_c | flip(bg_R_c, 3));                          % NB: hemisphere mask for individual mice (ctrl)
bg_L_e = recomputed_bkg_mask_4d_exp(:, :, 1:n_half, :);
bg_R_e = recomputed_bkg_mask_4d_exp(:, :, (n_full - n_half + 1):end, :);
mask_bg_exp  = logical(bg_L_e | flip(bg_R_e, 3));                          % NB: hemisphere mask for individual mice (exp)

% Find pixels that are tissue in at least one mouse
tissue_3d_ctrl = any(~mask_bg_ctrl, 4);
tissue_3d_exp  = any(~mask_bg_exp, 4);
brainMask_cropped_no_bkg_ctrl = brainMask_cropped & tissue_3d_ctrl;        % NB: average hemisphere mask (ctrl)
brainMask_cropped_no_bkg_exp  = brainMask_cropped & tissue_3d_exp;         % NB: average hemisphere mask (exp)
brainMask_group_diff = brainMask_cropped_no_bkg_ctrl & brainMask_cropped_no_bkg_exp; % NB: average hemisphere mask (diff)

if generate_diff_videos
    % Generate average ctrl group video
    write_lr_video(avg_lr_diff_ctrl, avg_lr_sum_ctrl, half_atlas, brainMask_cropped_no_bkg_ctrl, ...
        comp_out_dir, [['lr_diff_sum_' channel '_'] ctrl_type '.mp4'], [0, 2.5], ctrl_type, ...
        ['LR abs difference (' channel ') average - ' ctrl_type], ['LR abs sum (' channel ') average - ' ctrl_type]);
    % Generate average exp group video
    write_lr_video(avg_lr_diff_exp, avg_lr_sum_exp, half_atlas, brainMask_cropped_no_bkg_exp, ...
        comp_out_dir, [['lr_diff_sum_' channel '_'] exp_type '.mp4'], [0, 2.5], exp_type, ...
        ['LR abs difference (' channel ') average - ' exp_type], ['LR abs sum (' channel ') average - ' exp_type]);
    % Generate average group difference video
    write_lr_video(avg_lr_diff_groupdiff, avg_lr_sum_groupdiff, half_atlas, brainMask_group_diff, ...
        comp_out_dir, ['lr_diff_sum_' channel '_groupdiff_' comp_tag '.mp4'], [-2.5, 2.5], ...
        [exp_type ' - ' ctrl_type], ...
        ['LR abs diff groupdiff (' comp_tag ')'], ['LR abs sum groupdiff (' comp_tag ')']);
end
clear data_4d_new_ctrl data_4d_new_exp recomputed_bkg_mask_4d_ctrl recomputed_bkg_mask_4d_exp

%% Produce individual LR-asymmetry videos

if generate_individual_diff_videos

    fprintf('Generating Individual LR-Asymmetry Videos...\n');

    % Visualization Params
    clim_diff = [0 1.5];
    clim_sum  = [0 10];

    % Generate Control Group Video
    write_lr_indiv_video( ...
        abs(lr_diff_ctrl), ...                                  % Diff Data (4D) - ensuring absolute value
        lr_sum_ctrl, ...                                        % Sum Data (4D)
        mask_bg_ctrl, ...                                       % Background Masks (4D)
        AllenCrop, ...                                          % Atlas (Full 3D, function will slice it)
        comp_out_dir, ...                                       % Save Directory
        ['Individual_LR_Diff_Sum_' ctrl_type '.mp4'], ...       % Filename
        clim_diff, clim_sum, ...                                % Limits
        ctrl_type, ...                                          % Group Name
        ctrl_mousenames ...                                     % Mouse Names
        );
    % Generate Experimental Group Video
    write_lr_indiv_video( ...
        abs(lr_diff_exp), ...
        lr_sum_exp, ...
        mask_bg_exp, ...
        AllenCrop, ...
        comp_out_dir, ...
        ['Individual_LR_Diff_Sum_' exp_type '.mp4'], ...
        clim_diff, clim_sum, ...
        exp_type, ...
        exp_mousenames ...
        );

    fprintf('Individual videos generation complete.\n');

end

%% Produce individual LR-asymmetry videos (signed)

if generate_signed_diff_videos

    fprintf('Generating Signed (Directional) LR-Asymmetry Videos...\n');

    % Visualization Params
    clim_diff_signed = [-0.75 0.75];
    clim_sum  = [0 10];

    % Generate Control Group Video
    write_lr_indiv_video_signed( ...
        lr_diff_ctrl, ...
        lr_sum_ctrl, ...
        mask_bg_ctrl, ...
        AllenCrop, ...
        comp_out_dir, ...
        ['Individual_LR_Directional_' ctrl_type '.mp4'], ...
        clim_diff_signed, clim_sum, ...
        ctrl_type, ...
        ctrl_mousenames ...
        );

    % Generate Experimental Group Video
    write_lr_indiv_video_signed( ...
        lr_diff_exp, ...
        lr_sum_exp, ...
        mask_bg_exp, ...
        AllenCrop, ...
        comp_out_dir, ...
        ['Individual_LR_Directional_' exp_type '.mp4'], ...
        clim_diff_signed, clim_sum, ...
        exp_type, ...
        exp_mousenames ...
        );

    fprintf('Individual directional videos generation complete.\n');
end

%% Plot t-scored videos

sem_lr_diff_ctrl = nanstd(abs(lr_diff_ctrl), [], 4) ./ sqrt(sum(~isnan(lr_diff_ctrl),4)); %#ok<*NANSTD>
sem_lr_sum_ctrl  = nanstd(abs(lr_sum_ctrl),  [], 4) ./ sqrt(sum(~isnan(lr_sum_ctrl), 4));
sem_lr_diff_exp  = nanstd(abs(lr_diff_exp),  [], 4) ./ sqrt(sum(~isnan(lr_diff_exp),  4));
sem_lr_sum_exp   = nanstd(abs(lr_sum_exp),   [], 4) ./ sqrt(sum(~isnan(lr_sum_exp),   4));

% Avoid division by zero
sem_lr_diff_ctrl(sem_lr_diff_ctrl==0) = NaN;
sem_lr_sum_ctrl(sem_lr_sum_ctrl==0)   = NaN;
sem_lr_diff_exp(sem_lr_diff_exp==0)   = NaN;
sem_lr_sum_exp(sem_lr_sum_exp==0)     = NaN;

% Within-group t-scores
t_lr_diff_ctrl = avg_lr_diff_ctrl ./ sem_lr_diff_ctrl;
t_lr_sum_ctrl  = avg_lr_sum_ctrl  ./ sem_lr_sum_ctrl;
t_lr_diff_exp  = avg_lr_diff_exp  ./ sem_lr_diff_exp;
t_lr_sum_exp   = avg_lr_sum_exp   ./ sem_lr_sum_exp;

% Pooled SEM for between-group difference
sem_diff_lr_diff = sqrt(sem_lr_diff_ctrl.^2 + sem_lr_diff_exp.^2);
sem_diff_lr_sum  = sqrt(sem_lr_sum_ctrl.^2  + sem_lr_sum_exp.^2);
sem_diff_lr_diff(isnan(sem_diff_lr_diff) | sem_diff_lr_diff==0) = NaN;
sem_diff_lr_sum(isnan(sem_diff_lr_sum)  | sem_diff_lr_sum==0)  = NaN;

% Between-group t-scores (Welch)
t_lr_diff_groupdiff = avg_lr_diff_groupdiff ./ sem_diff_lr_diff;
t_lr_sum_groupdiff  = avg_lr_sum_groupdiff  ./ sem_diff_lr_sum;

if generate_t_scored_videos
    % Generate t-score groupdiff video
    t_lim = [-6 6];
    write_lr_video(t_lr_diff_groupdiff, t_lr_sum_groupdiff, half_atlas, brainMask_group_diff, ... % brainMask_cropped
        comp_out_dir, [['t_lr_diff_sum_' channel '_groupdiff_'] comp_tag '.mp4'], t_lim, ...
        [exp_type ' - ' ctrl_type ' (t-score)'], ...
        ['LR abs diff t-score (' comp_tag ')'], ['LR abs sum t-score (' comp_tag ')']);
end

%% Plot surprise and surprise-masked t-scored videos

% Get number of mice per group
n_ctrl = max(max(max(sum(~isnan(lr_diff_ctrl), 4))));
n_exp  = max(max(max(sum(~isnan(lr_diff_exp),  4))));

% Compute Welch–Satterthwaite degrees of freedom
var1_diff = sem_lr_diff_ctrl.^2;   var2_diff = sem_lr_diff_exp.^2;
var1_sum  = sem_lr_sum_ctrl.^2;    var2_sum  = sem_lr_sum_exp.^2;
df_diff = (var1_diff + var2_diff).^2 ./ (var1_diff.^2./(n_ctrl-1) + var2_diff.^2./(n_exp-1));
df_sum  = (var1_sum  + var2_sum ).^2 ./ (var1_sum.^2 ./ (n_ctrl-1) + var2_sum.^2 ./ (n_exp-1));
df_diff(n_ctrl < 2 | n_exp < 2) = NaN;
df_sum(n_ctrl  < 2 | n_exp < 2) = NaN;

% Compute t-zest p-values and surprise
p_diff = 2 * tcdf(-abs(t_lr_diff_groupdiff), df_diff);
p_sum  = 2 * tcdf(-abs(t_lr_sum_groupdiff),  df_sum);
surp_diff = -log10(p_diff);
surp_sum  = -log10(p_sum);

if generate_surprise_videos
    % Generate surprise video
    write_lr_video(surp_diff, surp_sum, half_atlas, brainMask_group_diff, ... % brainMask_cropped
        comp_out_dir, ['surp_lr_diff_sum_' channel '_groupdiff_' comp_tag '.mp4'], [0 8], ...
        ['-log_{10}(p) | ' comp_tag], ...
        ['LR abs diff surprise (' comp_tag ')'], ['LR abs sum surprise (' comp_tag ')']);
    % Generate surprise-masked t-score video
    surp_thresh = -log10(0.05);
    write_lr_video_surpmask( ...
        t_lr_diff_groupdiff, t_lr_sum_groupdiff, ...
        half_atlas, brainMask_group_diff, ... % brainMask_cropped
        comp_out_dir, [['t_lr_diff_sum_' channel '_groupdiff_'] comp_tag '_surpmask.mp4'], ...
        t_lim, [exp_type ' - ' ctrl_type ' (t-score, p<0.05)'], ...
        ['LR abs diff t-score masked (' comp_tag ')'], ['LR abs sum t-score masked (' comp_tag ')'], ...
        surp_diff, surp_thresh);
end

% Display output
disp('Generalized LR analysis completed successfully!');
fprintf('All comparison results saved to: %s\n', comp_out_dir);

%% Whole-Brain Ontology Analysis (montage)

if perform_area_based_analysis_fine

    fprintf('Starting Whole-Brain Multi-Metric Analysis (Leaf Nodes)...\n');

    % 1. Setup Atlas & ID Mapping
    [~, ~, n_width] = size(AllenCrop);
    half_width = floor(n_width / 2);
    atlas_left = AllenCrop(:, :, 1:half_width);

    % Find valid pixels (Atlas > 0)
    valid_mask = atlas_left > 0;
    all_ids_vec = atlas_left(valid_mask);

    % Map large Allen IDs to sequential indices 1..N
    [unique_ids, ~, id_indices] = unique(all_ids_vec);
    n_unique_regions = length(unique_ids);

    fprintf('  Found %d unique leaf regions.\n', n_unique_regions);

    % Define metrics
    metric_names = {'Mean', 'Median', 'Q1', 'Q3', 'P99'};
    n_metrics = length(metric_names);

    % Helper to compute quantiles ignoring NaNs
    nan_quant = @(x, p) quantile(x(~isnan(x)), p);
    funcs = {
        @(x) median(x, 'omitnan'), ...
        @(x) nan_quant(x, 0.25), ...
        @(x) nan_quant(x, 0.75), ...
        @(x) nan_quant(x, 0.99) ...
        };

    % Pre-slice masks (Common for both Diff and Sum)
    bg_mask_ctrl_left = mask_bg_ctrl;
    bg_mask_exp_left  = mask_bg_exp;

    % Analysis Loop: Diff and Sum
    analysis_types = {'Diff', 'Sum'};
    for a_idx = 1:length(analysis_types)

        res_type = analysis_types{a_idx}; % 'Diff' or 'Sum'
        fprintf('--- Whole-Brain Analysis: Processing %s Data ---\n', res_type);

        % Select Data Input
        if strcmp(res_type, 'Diff')
            input_4d_ctrl = lr_diff_ctrl;
            input_4d_exp  = lr_diff_exp;
        else
            input_4d_ctrl = lr_sum_ctrl;
            input_4d_exp  = lr_sum_exp;
        end

        % Preallocate Storage
        leaf_stats_ctrl = nan(n_unique_regions, n_ctrl, n_metrics);
        leaf_stats_exp  = nan(n_unique_regions, n_exp, n_metrics);

        % --- Extract Data & Compute Stats (Control) ---
        fprintf('  Computing stats for control group (%s)...\n', res_type);
        for m = 1:n_ctrl
            % Calculate abs() for one mouse at a time
            vol_data = abs(input_4d_ctrl(:,:,:,m));
            vol_bg   = bg_mask_ctrl_left(:,:,:,m);

            % Flatten to 1D vectors
            data_vec = vol_data(valid_mask);
            bg_vec   = vol_bg(valid_mask);

            % Set Background to NaN
            data_vec(bg_vec) = NaN;

            % Metric 1: Mean
            sums   = accumarray(id_indices, data_vec, [n_unique_regions 1], @nansum); %#ok<*NANSUM>
            counts = accumarray(id_indices, ~isnan(data_vec), [n_unique_regions 1], @sum);
            leaf_stats_ctrl(:, m, 1) = sums ./ counts;

            % Metrics 2-5: Robust Quantiles
            for f = 1:length(funcs)
                leaf_stats_ctrl(:, m, f+1) = accumarray(id_indices, data_vec, [n_unique_regions 1], funcs{f});
            end
        end

        % --- Extract Data & Compute Stats (Experimental) ---
        fprintf('  Computing stats for experimental group (%s)...\n', res_type);
        for m = 1:n_exp
            vol_data = abs(input_4d_exp(:,:,:,m));
            vol_bg   = bg_mask_exp_left(:,:,:,m);

            data_vec = vol_data(valid_mask);
            bg_vec   = vol_bg(valid_mask);
            data_vec(bg_vec) = NaN;

            % Metric 1: Mean
            sums   = accumarray(id_indices, data_vec, [n_unique_regions 1], @nansum);
            counts = accumarray(id_indices, ~isnan(data_vec), [n_unique_regions 1], @sum);
            leaf_stats_exp(:, m, 1) = sums ./ counts;

            % Metrics 2-5: Robust Quantiles
            for f = 1:length(funcs)
                leaf_stats_exp(:, m, f+1) = accumarray(id_indices, data_vec, [n_unique_regions 1], funcs{f});
            end
        end

        % Visualization Loop

        % Setup Visualization Parameters
        slices_to_show = 150:50:size(atlas_left,1)-150;
        n_cols = 5;
        n_rows = ceil(length(slices_to_show)/n_cols);

        % Custom Red-White-Blue Colormap
        n_steps = 256;
        c_blue = [0 0 1]; c_white = [1 1 1]; c_red = [1 0 0];
        cmap_neg = [linspace(c_blue(1), c_white(1), n_steps/2)', linspace(c_blue(2), c_white(2), n_steps/2)', linspace(c_blue(3), c_white(3), n_steps/2)'];
        cmap_pos = [linspace(c_white(1), c_red(1), n_steps/2)', linspace(c_white(2), c_red(2), n_steps/2)', linspace(c_white(3), c_red(3), n_steps/2)'];
        custom_cmap = [cmap_neg; cmap_pos];

        fprintf('  Generating T-Maps for %s metrics...\n', res_type);

        for i_met = 1:n_metrics
            metric_name = metric_names{i_met};

            % Get Data for this metric
            data_c = squeeze(leaf_stats_ctrl(:, :, i_met));
            data_e = squeeze(leaf_stats_exp(:, :, i_met));

            % Compute T-Scores ... NB: T = (Mean_Exp - Mean_Ctrl) / Pooled_SEM
            mu_c = nanmean(data_c, 2);
            mu_e = nanmean(data_e, 2);
            diff_mu = mu_e - mu_c;

            sem_c = nanstd(data_c, [], 2) ./ sqrt(n_ctrl);
            sem_e = nanstd(data_e, [], 2) ./ sqrt(n_exp);
            pooled_sem = sqrt(sem_c.^2 + sem_e.^2);

            t_scores_vec = diff_mu ./ pooled_sem;

            % Filter invalid T-scores
            t_scores_vec(isnan(t_scores_vec) | isinf(t_scores_vec)) = 0;

            % Determine range
            t_vals = t_scores_vec(t_scores_vec ~= 0);
            if isempty(t_vals)
                max_t = 1;
                fprintf('    [Warning] Metric %s yielded all zero T-scores.\n', metric_name);
            else
                max_t = quantile(abs(t_vals), 0.95); % 95th percentile of absolute T
                if max_t < 0.1, max_t = 0.1; end % Avoid tiny ranges
            end
            t_lims = [-max_t, max_t];

            % Reconstruct 3D volume
            t_score_vol = zeros(size(atlas_left), 'single');
            % Map indices back to volume
            t_score_valid_pixels = t_scores_vec(id_indices);
            t_score_vol(valid_mask) = t_score_valid_pixels;

            % Plot Montage
            fig_h = figure('Visible', 'off', 'Name', ['WholeBrain_TMap_' res_type '_' metric_name '_' comp_tag], 'Color', 'k', 'Position', [50 50 1200 900]);

            for k = 1:length(slices_to_show)
                s_idx = slices_to_show(k);
                subplot(n_rows, n_cols, k);
                im_slice = squeeze(t_score_vol(s_idx, :, :));
                mask_slice = squeeze(atlas_left(s_idx, :, :));
                alpha_data = double(mask_slice > 0);
                imagesc(im_slice);
                axis image; axis off;
                set(gca, 'Color', 'k');
                set(findobj(gca, 'Type', 'image'), 'AlphaData', alpha_data);
                colormap(gca, custom_cmap);
                clim(t_lims);
                title(['Slice ' num2str(s_idx)], 'Color', 'w', 'FontSize', 8);
            end
            c = colorbar;
            c.Position = [0.92 0.1 0.02 0.8]; c.Color = 'w';
            c.Label.String = sprintf('T-score (%s %s)', res_type, metric_name);
            colormap(c, custom_cmap);
            clim(t_lims);

            sgtitle(['Whole-Brain T-Scores (' res_type '): ' metric_name], 'Color', 'w', 'FontSize', 14);

            % Save
            set(fig_h, 'InvertHardcopy', 'off');
            save_base = fullfile(comp_out_dir, ['WholeBrain_TMap_' res_type '_' metric_name '_' comp_tag]);
            saveas(fig_h, [save_base '.fig']);
            exportgraphics(fig_h, [save_base '.png'], 'Resolution', 300, 'BackgroundColor', 'current');

            % Save volume data
            save([save_base '.mat'], 't_score_vol', 't_scores_vec', 'unique_ids', 'metric_name', 'res_type');
        end
    end

    % Clean up large masks
    clear bg_mask_ctrl_left bg_mask_exp_left t_score_vol t_score_valid_pixels input_4d_ctrl input_4d_exp

    fprintf('Whole-brain multi-metric analysis (Diff & Sum) saved to: %s\n', comp_out_dir);

end

%% Area-based Analysis: Vectorized Target Region Quantification (Diff & Sum)

if perform_area_based_analysis_coarse

    fprintf('Starting Area-based quantification (Vectorized)...\n');

    % 0. Define metric to plot
    choosen_metric = 'P99';

    % Define Areas for Detailed Inspection
    rois_to_inspect = {};

    % 1. Define Target Regions
    roi_list = { ...
        'Primary somatosensory area, barrel field', ...
        'Primary somatosensory area, trunk', ...
        'Primary somatosensory area, upper limb', ...
        'Primary somatosensory area, lower limb', ...
        'Supplemental somatosensory area', ...
        'Primary motor area', ...
        'Secondary motor area', ...
        'Primary visual area', ...
        'Lateral visual area', ...
        'Anterolateral visual area', ...
        'Anteromedial visual area', ...
        'Primary auditory area', ...
        'Dorsal auditory area', ...
        'Ventral auditory area', ...
        'Anterior cingulate area', ...
        'Olfactory tubercle', ...
        'Prelimbic area', ...
        'Infralimbic area', ...
        'Visceral area',...
        'Gustatory areas',...
        'Piriform area',...
        'Subiculum',...
        'Orbital area', ...
        'Claustrum', ...
        'Agranular insular area', ...
        'Anterior area', ...
        'Rostrolateral visual area', ...
        'Temporal association areas', ...
        'Perirhinal area', ...
        'Ectorhinal area', ...
        'Retrosplenial area', ...
        'Nucleus accumbens', ...
        'Caudoputamen', ...
        'Globus pallidus, external segment', ...
        'Subthalamic nucleus', ...
        'Ventral posteromedial nucleus of the thalamus', ...
        'Ventral posterolateral nucleus of the thalamus', ...
        'Ventral medial nucleus of the thalamus',...
        'Zona incerta',...
        'Posterior complex of the thalamus', ...
        'Lateral posterior nucleus of the thalamus', ...
        'Lateral dorsal nucleus of thalamus', ...
        'Ventral anterior-lateral complex of the thalamus', ...
        'Mediodorsal nucleus of the thalamus', ...
        'Parafascicular nucleus', ...
        'Nucleus of reuniens', ...
        'Central lateral nucleus of the thalamus', ...
        'Reticular nucleus of the thalamus', ...
        'Geniculate group, dorsal thalamus', ...
        'Midbrain, motor related', ...
        'Superior colliculus, motor related', ...
        'Superior colliculus, sensory related', ...
        'Hippocampal formation', ...
        'Basolateral amygdalar nucleus', ...
        'Parafascicular nucleus', ...
        'Hypothalamus' ...
        };
    n_rois = length(roi_list);

    % 2. Setup Atlas Vectors
    fprintf('Mapping Atlas Volume...\n');
    atlas_left = AllenCrop(:, :, 1:half_width);
    valid_pixels = atlas_left > 0;
    pixel_ids = atlas_left(valid_pixels);

    % 3. Pre-calculate Masks
    fprintf('Building masks for %d regions (Handling overlaps)...\n', n_rois);

    roi_masks = false(length(pixel_ids), n_rois);
    roi_pixel_counts = zeros(n_rois, 1);

    for r = 1:n_rois
        region_name = roi_list{r};
        try
            % Get hierarchical IDs for this region name
            mask_temp = get_allen_region_mask(allenDir, AllenCrop, {region_name}, brainMask, '');
            if size(mask_temp, 3) >= half_width
                mask_temp = mask_temp(:, :, 1:half_width);
            end
            % Get unique IDs belonging to this region
            ids_in_region = unique(atlas_left(mask_temp));
            ids_in_region(ids_in_region == 0) = [];
            % Create boolean mask for the flattened pixel vector
            roi_masks(:, r) = ismember(pixel_ids, ids_in_region);
            roi_pixel_counts(r) = sum(roi_masks(:, r));
        catch
            warning('Could not map region: %s', region_name);
        end
    end

    % Warning for empty regions
    empty_rois = find(roi_pixel_counts == 0);
    if ~isempty(empty_rois)
        fprintf('Warning: The following regions have 0 pixels:\n');
        disp(roi_list(empty_rois)');
    end

    % 4. Initialize Stats Storage
    metric_names = {'Mean', 'Median', 'Q1', 'Q3', 'P99'};
    n_metrics = length(metric_names);
    roi_stats_ctrl = nan(n_rois, n_ctrl, n_metrics);
    roi_stats_exp  = nan(n_rois, n_exp, n_metrics);

    % Pre-slice background masks
    bg_mask_ctrl = mask_bg_ctrl;
    bg_mask_exp  = mask_bg_exp;

    % Inspection Setup
    inspect_indices = find(ismember(roi_list, rois_to_inspect));
    dist_out_dir = fullfile(comp_out_dir, 'distribution_checks');
    if ~exist(dist_out_dir, 'dir') && ~isempty(inspect_indices), mkdir(dist_out_dir); end

    % 5. Analysis Loop
    analysis_types = {'Diff', 'Sum'};

    for a_idx = 1:length(analysis_types)
        res_type = analysis_types{a_idx};
        fprintf('--- Processing %s Data for ROIs ---\n', res_type);

        % Define Limits based on Type
        if strcmp(res_type, 'Diff')
            input_vol_ctrl = abs(lr_diff_ctrl);
            input_vol_exp  = abs(lr_diff_exp);
            hist_xlim = [0 2.5];
        else
            input_vol_ctrl = abs(lr_sum_ctrl);
            input_vol_exp  = abs(lr_sum_exp);
            hist_xlim = [0 5];
        end

        % Initialize Inspection Figures for this Analysis Type
        inspect_figs = gobjects(length(inspect_indices), 1);
        for k = 1:length(inspect_indices)
            r_idx = inspect_indices(k);
            inspect_figs(k) = figure('Visible', 'off', 'Name', ['Dist_' roi_list{r_idx} '_' res_type], ...
                'Color', 'w', 'Visible', 'off', 'Position', [100 100 300*max(n_ctrl,n_exp) 600]);
        end

        % --- Control Group ---
        for m = 1:n_ctrl
            tic
            % Flatten Data for this mouse
            vol_data = input_vol_ctrl(:,:,:,m);
            vol_bg   = bg_mask_ctrl(:,:,:,m);

            data_vec = vol_data(valid_pixels);
            bg_vec   = vol_bg(valid_pixels);
            data_vec(bg_vec) = NaN;

            % Loop over ROIs
            for r = 1:n_rois
                if roi_pixel_counts(r) == 0, continue; end
                % Extract pixels for this ROI
                roi_vals = data_vec(roi_masks(:, r));

                % Compute Stats
                roi_stats_ctrl(r, m, 1) = mean(roi_vals, 'omitnan'); % Mean

                % Compute other metrics
                roi_vals_clean = roi_vals(~isnan(roi_vals));
                if ~isempty(roi_vals_clean)
                    roi_stats_ctrl(r, m, 2) = median(roi_vals_clean);
                    roi_stats_ctrl(r, m, 3) = quantile(roi_vals_clean, 0.25);
                    roi_stats_ctrl(r, m, 4) = quantile(roi_vals_clean, 0.75);
                    roi_stats_ctrl(r, m, 5) = quantile(roi_vals_clean, 0.99);
                    % --- Inspection Plot (Control) ---
                    if ismember(r, inspect_indices)
                        k_fig = find(inspect_indices == r);
                        set(0, 'CurrentFigure', inspect_figs(k_fig));
                        subplot(2, max(n_ctrl, n_exp), m);
                        histogram(roi_vals_clean, 100, 'EdgeColor', 'none', 'FaceColor', [0 0.45 0.74]); hold on;
                        xline(roi_stats_ctrl(r, m, 1), 'k-', 'LineWidth', 1.5);
                        xline(roi_stats_ctrl(r, m, 2), 'k--', 'LineWidth', 1.5);
                        xline(roi_stats_ctrl(r, m, 5), '-', 'LineWidth', 1.5, 'Color',[1,0,1]);
                        xlim(hist_xlim);
                        txt_str = [sprintf('Mean: %.2f\nP99: %.2f', roi_stats_ctrl(r, m, 1), roi_stats_ctrl(r, m, 5)),' (N=',num2str(numel(roi_vals_clean)),')'];
                        text(0.95, 0.9, txt_str, 'Units', 'normalized', 'HorizontalAlignment', 'right', ...
                            'FontSize', 8, 'BackgroundColor', 'w', 'EdgeColor', 'k');
                        if exist('ctrl_mousenames','var'), t_str=ctrl_mousenames{m}; else, t_str=sprintf('M%d',m); end
                        title(t_str, 'Interpreter', 'none', 'FontSize', 8);
                        if m==1, ylabel(['Control (' res_type ')'], 'FontWeight', 'bold'); end
                        grid on;
                    end
                end
            end
            toc
        end

        % --- Experimental Group ---
        for m = 1:n_exp
            tic
            vol_data = input_vol_exp(:,:,:,m);
            vol_bg   = bg_mask_exp(:,:,:,m);
            data_vec = vol_data(valid_pixels);
            bg_vec   = vol_bg(valid_pixels);
            data_vec(bg_vec) = NaN;
            for r = 1:n_rois
                if roi_pixel_counts(r) == 0, continue; end
                roi_vals = data_vec(roi_masks(:, r));
                roi_stats_exp(r, m, 1) = mean(roi_vals, 'omitnan');

                roi_vals_clean = roi_vals(~isnan(roi_vals));
                if ~isempty(roi_vals_clean)
                    roi_stats_exp(r, m, 2) = median(roi_vals_clean);
                    roi_stats_exp(r, m, 3) = quantile(roi_vals_clean, 0.25);
                    roi_stats_exp(r, m, 4) = quantile(roi_vals_clean, 0.75);
                    roi_stats_exp(r, m, 5) = quantile(roi_vals_clean, 0.99);
                    % --- Inspection Plot (Exp) ---
                    if ismember(r, inspect_indices)
                        k_fig = find(inspect_indices == r);
                        set(0, 'CurrentFigure', inspect_figs(k_fig));
                        subplot(2, max(n_ctrl, n_exp), m + max(n_ctrl, n_exp));
                        histogram(roi_vals_clean, 100, 'EdgeColor', 'none', 'FaceColor', [0.85 0.33 0.10]); hold on;
                        xline(roi_stats_exp(r, m, 1), 'k-', 'LineWidth', 1.5);
                        xline(roi_stats_exp(r, m, 2), 'k--', 'LineWidth', 1.5);
                        xline(roi_stats_exp(r, m, 5), '-', 'LineWidth', 1.5, 'Color',[1,0,1]);
                        xlim(hist_xlim);
                        txt_str = [sprintf('Mean: %.2f\nP99: %.2f', roi_stats_exp(r, m, 1), roi_stats_exp(r, m, 5)),' (N=',num2str(numel(roi_vals_clean)),')'];
                        text(0.95, 0.9, txt_str, 'Units', 'normalized', 'HorizontalAlignment', 'right', ...
                            'FontSize', 8, 'BackgroundColor', 'w', 'EdgeColor', 'k');
                        if exist('exp_mousenames','var'), t_str=exp_mousenames{m}; else, t_str=sprintf('M%d',m); end
                        title(t_str, 'Interpreter', 'none', 'FontSize', 8);
                        if m==1, ylabel(['Exp (' res_type ')'], 'FontWeight', 'bold'); end
                        grid on;
                    end
                end
            end
            toc
        end

        % --- Save Inspection Figures ---
        for k = 1:length(inspect_figs)
            if isvalid(inspect_figs(k))
                r_idx = inspect_indices(k);
                r_name = roi_list{r_idx};
                clean_name = regexprep(r_name, '[^a-zA-Z0-9]', '_');
                sgtitle(inspect_figs(k), ['Distribution: ' r_name ' (' res_type ') | Solid=Mean, Dash=Med, Purple=P99'], 'Interpreter', 'none');

                saveas(inspect_figs(k), fullfile(dist_out_dir, ['Dist_' res_type '_' clean_name '.fig']));
                exportgraphics(inspect_figs(k), fullfile(dist_out_dir, ['Dist_' res_type '_' clean_name '.png']));
                close(inspect_figs(k));
            end
        end

        % --- Plotting: Bar Chart ---
        target_metric_idx = find(strcmp(metric_names, choosen_metric));
        if isempty(target_metric_idx), target_metric_idx = 1; end

        roi_means_ctrl = roi_stats_ctrl(:, :, target_metric_idx);
        roi_means_exp  = roi_stats_exp(:, :, target_metric_idx);

        mean_ctrl_roi = nanmean(roi_means_ctrl, 2);
        mean_exp_roi  = nanmean(roi_means_exp, 2);
        diff_means    = mean_exp_roi - mean_ctrl_roi;

        sem_c = nanstd(roi_means_ctrl,[],2) ./ sqrt(n_ctrl);
        sem_e = nanstd(roi_means_exp, [],2) ./ sqrt(n_exp);
        pooled_sem = sqrt(sem_c.^2 + sem_e.^2);
        t_score_roi = diff_means ./ pooled_sem;

        valid_rows = ~isnan(t_score_roi);
        t_score_roi = t_score_roi(valid_rows);
        current_rois = roi_list(valid_rows);

        [sorted_t, sort_idx] = sort(t_score_roi, 'ascend');
        sorted_rois = current_rois(sort_idx);

        fig_bars = figure('Visible', 'off', 'Name', ['Region_Analysis_BarChart_' res_type '_' comp_tag], ...
            'Color', 'w', 'Units', 'Normalized', 'Position', [0 0 0.9 0.9]);

        b = barh(sorted_t); b.FaceColor = 'flat';
        for k = 1:length(sorted_t)
            if sorted_t(k) > 0, b.CData(k,:) = [0.85 0.33 0.10]; else, b.CData(k,:) = [0 0.45 0.74]; end
        end

        yticks(1:length(sorted_rois)); yticklabels(sorted_rois);
        xlabel(['t-score (' exp_type ' - ' ctrl_type ')']);
        title(['Regional LR - ' res_type ' differences - ', strrep(comp_tag,'_',' '), ' - ', choosen_metric]);
        grid on; set(gca, 'FontSize', 10); ylim([0 length(sorted_rois)+1]);

        % Highlighting
        switch exp_type
            case 'rws', highlighted_areas = {'Primary somatosensory area, barrel field', 'Ventral posteromedial nucleus of the thalamus', 'Posterior complex of the thalamus', 'Supplemental somatosensory area', 'Zona incerta'};
            case 'behavior', highlighted_areas = {'Primary somatosensory area, barrel field', 'Ventral posteromedial nucleus of the thalamus', 'Posterior complex of the thalamus', 'Supplemental somatosensory area', 'Zona incerta'};
            otherwise, highlighted_areas = {};
        end

        ax = gca;
        ytl = ax.YTickLabel; colored_labels = repmat({''}, size(ytl));
        for i = 1:length(ytl)
            if ismember(ytl{i}, highlighted_areas), colored_labels{i} = ['\color{magenta} \bf ' strrep(ytl{i}, '_', ' ') ]; else, colored_labels{i} = ['\color{black} ' strrep(ytl{i}, '_', ' ') ]; end
        end
        ax = gca; ax.YTickLabel = colored_labels;

        % Significance lines
        n_regions_tested = length(sorted_t);
        df = n_ctrl + n_exp - 2;
        t_crit_bonf = tinv(1 - 0.05/(2*n_regions_tested), df);
        hold on; xline(t_crit_bonf, 'k:', 'LineWidth', 2); xline(-t_crit_bonf, 'k:', 'LineWidth', 2); hold off;

        set(fig_bars, 'InvertHardcopy', 'off');
        saveas(fig_bars, fullfile(comp_out_dir, ['Region_Stats_Bar_' res_type '_' comp_tag '.fig']));
        exportgraphics(fig_bars, fullfile(comp_out_dir, ['Region_Stats_Bar_' res_type '_' comp_tag '.png']), 'Resolution', 300);
    end

    clear data_vec roi_masks atlas_left input_vol_ctrl input_vol_exp vol vol_smoothed bg_vec roi_vals roi_vals_clean
    fprintf(' Area-based analysis saved to: %s\n', comp_out_dir);

end

%% Whole-Brain Ontology Analysis (annotated video)

if perform_area_based_analysis_fine

    fprintf('Generating Annotated T-Score Video (Direct Mapping)...\n');

    % Select metric to plot
    metric_to_plot = 'P99';

    % 1. Load Mapping Table (Pixel Index -> Acronym)
    map_file = fullfile(allenDir, 'parcellation_to_parcellation_term_membership.csv');
    id2acronym = containers.Map('KeyType', 'double', 'ValueType', 'char');

    if exist(map_file, 'file')
        % Read table safely
        opts = detectImportOptions(map_file);
        opts.VariableTypes = repmat({'string'}, 1, length(opts.VariableTypes));
        T_map = readtable(map_file, opts);
        % Check for required columns
        if ismember('parcellation_index', T_map.Properties.VariableNames) && ...
                ismember('parcellation_term_acronym', T_map.Properties.VariableNames)
            % Extract Columns
            raw_pixels = str2double(T_map.parcellation_index);
            raw_acros  = T_map.parcellation_term_acronym;
            % Hierarchy filtering
            if ismember('parcellation_term_set_name', T_map.Properties.VariableNames)
                term_sets = lower(T_map.parcellation_term_set_name);
                is_leaf = contains(term_sets, 'structure');
                if sum(is_leaf) == 0
                    warning('No "structure" level found. Using all terms (may get parent categories).');
                    is_leaf = true(size(raw_pixels));
                end
            else
                is_leaf = true(size(raw_pixels));
            end
            % Apply Filter
            leaf_pixels = raw_pixels(is_leaf);
            leaf_acros  = raw_acros(is_leaf);
            % Remove NaNs
            valid_idx = ~isnan(leaf_pixels);
            leaf_pixels = leaf_pixels(valid_idx);
            leaf_acros = cellstr(leaf_acros(valid_idx));
            % Handle Duplicates
            [unique_pixels, idx] = unique(leaf_pixels, 'last');
            unique_acros = leaf_acros(idx);
            % Build Map
            if ~isempty(unique_pixels)
                id2acronym = containers.Map(unique_pixels, unique_acros);
                fprintf('  Mapping created successfully: %d leaf regions mapped.\n', length(unique_pixels));
            else
                warning('Mapping failed. No valid pixel-acronym pairs found.');
            end
        else
            warning('Required columns (parcellation_index, parcellation_term_acronym) missing.');
        end
    else
        warning('Mapping CSV not found. Labels will be skipped.');
    end

    % 2. Load T-Score Data & Atlas
    file_diff = fullfile(comp_out_dir, ['WholeBrain_TMap_Diff_' metric_to_plot '_' comp_tag '.mat']);
    file_sum  = fullfile(comp_out_dir, ['WholeBrain_TMap_Sum_'  metric_to_plot '_' comp_tag '.mat']);

    if exist(file_diff, 'file') && exist(file_sum, 'file')
        data_diff = load(file_diff, 't_score_vol');
        vol_diff = data_diff.t_score_vol;

        data_sum = load(file_sum, 't_score_vol');
        vol_sum = data_sum.t_score_vol;

        % Atlas preparation
        atlas_hi_res = double(AllenCrop);
        atlas_left_hires = atlas_hi_res(:, :, 1:size(vol_diff,3));

        % 3. Setup Video Writer
        video_filename = ['Annotated_TMap_' metric_to_plot '_' comp_tag '.mp4'];
        full_video_path = fullfile(comp_out_dir, video_filename);

        vidObj = VideoWriter(full_video_path, 'MPEG-4');
        vidObj.FrameRate = 10;
        vidObj.Quality = 95;
        open(vidObj);
        n_slices = size(vol_diff, 1);

        % Colormap
        n_steps = 256;
        c_blue = [0 0 1]; c_white = [1 1 1]; c_red = [1 0 0];
        cmap_neg = [linspace(c_blue(1), c_white(1), n_steps/2)', linspace(c_blue(2), c_white(2), n_steps/2)', linspace(c_blue(3), c_white(3), n_steps/2)'];
        cmap_pos = [linspace(c_white(1), c_red(1), n_steps/2)', linspace(c_white(2), c_red(2), n_steps/2)', linspace(c_white(3), c_red(3), n_steps/2)'];
        custom_cmap = [cmap_neg; cmap_pos];

        fh = figure('visible', 'off', 'units', 'normalized', 'outerposition', [0 0 1 1], 'Color', 'k');
        set(fh, 'InvertHardcopy', 'off');

        fprintf('Writing annotated video: %s\n', video_filename);

        for j = 1:n_slices

            mask_slice_hires = squeeze(atlas_left_hires(j, :, :));
            if sum(mask_slice_hires(:) > 0) == 0, continue; end %#ok<LOGSUM>

            alpha_data = double(mask_slice_hires > 0);

            % --- Pre-calculate Labels ---
            regions_in_slice = unique(mask_slice_hires(mask_slice_hires > 0));

            lbl_data = struct('x', {}, 'y', {}, 'str', {});
            idx_lbl = 1;

            if ~isempty(id2acronym)
                for k = 1:length(regions_in_slice)
                    r_id = regions_in_slice(k);

                    bin_mask = (mask_slice_hires == r_id);
                    if sum(bin_mask(:)) < 80, continue; end

                    if isKey(id2acronym, r_id)
                        [py, px] = find(bin_mask);
                        lbl_data(idx_lbl).x = mean(px);
                        lbl_data(idx_lbl).y = mean(py);
                        lbl_data(idx_lbl).str = id2acronym(r_id);
                        idx_lbl = idx_lbl + 1;
                    end
                end
            end

            % --- Plot Diff ---
            subplot(1, 2, 1);
            imagesc(squeeze(vol_diff(j, :, :)));
            clim([-4 4]); colormap(gca, custom_cmap);
            c = colorbar;
            c.Color = 'w';
            c.Label.String = sprintf('T-score (%s %s)', 'Diff', metric_to_plot);
            set(findobj(gca, 'Type', 'image'), 'AlphaData', alpha_data);
            axis image; axis off; set(gca, 'Color', 'k');
            title(['Diff T-Score (' metric_to_plot ')'], 'Color', 'w', 'FontSize', 12);

            if ~isempty(lbl_data)
                text([lbl_data.x], [lbl_data.y], {lbl_data.str}, ...
                    'Color', 'k', 'FontSize', 6, 'FontWeight', 'bold', ...
                    'HorizontalAlignment', 'center', 'Interpreter', 'none', ...
                    'BackgroundColor', 'w', 'Margin', 0.5, 'EdgeColor', 'none');
            end

            % --- Plot Sum ---
            subplot(1, 2, 2);
            imagesc(squeeze(vol_sum(j, :, :)));
            clim([-4 4]); colormap(gca, custom_cmap);
            c = colorbar;
            c.Color = 'w';
            c.Label.String = sprintf('T-score (%s %s)', 'Sum', metric_to_plot);
            set(findobj(gca, 'Type', 'image'), 'AlphaData', alpha_data);
            axis image; axis off; set(gca, 'Color', 'k');
            title(['Sum T-Score (' metric_to_plot ')'], 'Color', 'w', 'FontSize', 12);

            if ~isempty(lbl_data)
                text([lbl_data.x], [lbl_data.y], {lbl_data.str}, ...
                    'Color', 'k', 'FontSize', 6, 'FontWeight', 'bold', ...
                    'HorizontalAlignment', 'center', 'Interpreter', 'none', ...
                    'BackgroundColor', 'w', 'Margin', 0.5, 'EdgeColor', 'none');
            end

            sgtitle(['Slice # ' num2str(j)], 'Color', 'w', 'FontSize', 14);

            frame = getframe(fh);
            writeVideo(vidObj, frame);

            if mod(j, 50) == 0, fprintf('  Frame %d written...\n', j); end
            clf(fh);
        end

        close(vidObj);
        close(fh);
        fprintf('Annotated video saved: %s\n', full_video_path);

    else
        warning('T-Map .mat files not found. Run the Whole-Brain Analysis section first.');
    end
    clear atlas_hi_res atlas_left_hires

end

%% Plot Surprise-Masked Slab (Detailed inspection of rolling median around target slice)

% Set parameters
target_slice = 565;         % 565;
slab_range   = 10;          % +/- slices
save_fig     = true;

% Define slab indices
z_start = max(1, target_slice - slab_range);
z_end   = min(size(t_lr_diff_groupdiff, 1), target_slice + slab_range);
z_indices = z_start:z_end;

fprintf('Averaging signal across slices %d to %d (Target: %d)...\n', z_start, z_end, target_slice);

% Extract masks
mask_slab_3d = brainMask_group_diff(z_indices, :, :);

% Diff data
raw_diff = t_lr_diff_groupdiff(z_indices, :, :);
raw_diff(~mask_slab_3d) = NaN;
slab_diff = squeeze(nanmedian(raw_diff, 1));

% Sum data
raw_sum = t_lr_sum_groupdiff(z_indices, :, :);
raw_sum(~mask_slab_3d) = NaN;
slab_sum = squeeze(nanmedian(raw_sum, 1));

% Surprise data (diff)
raw_surp_d = surp_diff(z_indices, :, :);
raw_surp_d(~mask_slab_3d) = NaN;
slab_surp_diff = squeeze(nanmedian(raw_surp_d, 1));

% Surprise data (sum)
raw_surp_s = surp_sum(z_indices, :, :);
raw_surp_s(~mask_slab_3d) = NaN;
slab_surp_sum = squeeze(nanmedian(raw_surp_s, 1));

% Prepare alpha masks
p_thresh = 0.01;
surp_thresh = -log10(p_thresh);

% Helper to calculate alpha
calc_alpha = @(vol) min(1, max(0, vol ./ surp_thresh));
% Alpha for diff
alpha_diff = calc_alpha(slab_surp_diff);
alpha_diff(isnan(alpha_diff)) = 0;
% Alpha for sum
alpha_sum  = calc_alpha(slab_surp_sum);
alpha_sum(isnan(alpha_sum)) = 0;

% Combine with structural mask
slab_mask_2d = squeeze(max(mask_slab_3d, [], 1));
alpha_mask_diff = alpha_diff .* double(slab_mask_2d);
alpha_mask_sum  = alpha_sum  .* double(slab_mask_2d);

% Prepare atlas boundaries (from central slice)
atlas_slice = squeeze(half_atlas(target_slice, :, 1:size(slab_diff,2)));
[gy, gx] = gradient(single(atlas_slice));
boundaries = (abs(gx) + abs(gy)) > 0 & (atlas_slice > 0);
[b_row, b_col] = find(boundaries);

% Plotting
fig_slab = figure('Visible', 'off', 'Name', sprintf('Slab_Avg_%d', target_slice), 'Color', 'k', 'Position', [100 100 1200 600]);
set(fig_slab, 'InvertHardcopy', 'off');

% Left: diff T-Score
subplot(1, 2, 1);
h1 = imagesc(slab_diff);
set(h1, 'AlphaData', alpha_mask_diff);
axis image; axis off; set(gca, 'Color', 'k');
hold on;
plot(b_col, b_row, '.', 'Color', [0.7 0.7 0.7], 'MarkerSize', 0.25);
clim([-6 6]); % T-score limits
colormap(gca, get_color2color_colormap([0 0 1], [1 0 0]));
cb1 = colorbar; cb1.Color = 'w'; cb1.Label.String = 'T-Score (Diff)'; cb1.Label.Color = 'w';
title(['LR Diff (T-Score) - Slab Avg ' num2str(target_slice) '\pm' num2str(slab_range)], 'Color', 'w');

% Right: sum T-Score
subplot(1, 2, 2);
h2 = imagesc(slab_sum);
set(h2, 'AlphaData', alpha_mask_sum);
axis image; axis off; set(gca, 'Color', 'k');
hold on;
plot(b_col, b_row, '.', 'Color', [0.7 0.7 0.7], 'MarkerSize', 0.25);
clim([-6 6]);
colormap(gca, get_color2color_colormap([0 0 1], [1 0 0]));
cb2 = colorbar; cb2.Color = 'w'; cb2.Label.String = 'T-Score (Sum)'; cb2.Label.Color = 'w';
title(['LR Sum (T-Score) - Slab Avg ' num2str(target_slice) '\pm' num2str(slab_range)], 'Color', 'w');

sgtitle(['Slab average - ',exp_type ' vs ' ctrl_type ' - surprise masked (p<',num2str(p_thresh),')'], 'Color', 'w', 'FontSize', 14);

% Save
if save_fig
    if ~exist(comp_out_dir, 'dir'), mkdir(comp_out_dir); end
    filename = sprintf('Slab_Avg_%d_%s_surpmask', target_slice, comp_tag);
    saveas(fig_slab, fullfile(comp_out_dir, [filename '.fig']));
    exportgraphics(fig_slab, fullfile(comp_out_dir, [filename '.png']), 'Resolution', 300, 'BackgroundColor', 'current');
    fprintf('Saved slab average figure to: %s\n', fullfile(comp_out_dir, filename));
end

%% Plot Individual Mice Slab Average (Detailed inspection of rolling median around target slice)

fprintf('Generating Individual Mice Slab Average Plots...\n');

% Set parameters
target_slice = 565;                 % 565;
slab_range   = 10;                  % +/- slices
clim_diff_indiv = [0 1.5];
clim_sum_indiv  = [0 10];
bool_overlay_atlas = false;

% Define slab indices
z_start = max(1, target_slice - slab_range);
z_end   = min(size(lr_diff_ctrl, 1), target_slice + slab_range);
z_indices = z_start:z_end;

% Prepare atlas boundaries (from central slice)
atlas_slice = squeeze(half_atlas(target_slice, :, 1:size(lr_diff_ctrl,3)));
[gy, gx] = gradient(single(atlas_slice));
boundaries = (abs(gx) + abs(gy)) > 0 & (atlas_slice > 0);
[b_row, b_col] = find(boundaries);

% Loop over groups
group_names     = {ctrl_type, exp_type};
group_data_diff = {abs(lr_diff_ctrl), abs(lr_diff_exp)};
group_data_sum  = {lr_sum_ctrl, lr_sum_exp};
group_masks     = {mask_bg_ctrl, mask_bg_exp};
group_mice_list = {ctrl_mousenames, exp_mousenames};

for g_idx = 1:2
    curr_grp  = group_names{g_idx};
    curr_diff = group_data_diff{g_idx};
    curr_sum  = group_data_sum{g_idx};
    curr_mask = group_masks{g_idx};
    curr_mice = group_mice_list{g_idx};

    n_mice = size(curr_diff, 4);

    % Initialize Figure
    fig_indiv = figure('Visible', 'off', 'Name', ['Individual_Slab_Avg_' curr_grp], ...
        'Color', 'k', 'Units', 'normalized', 'Position', [-0.05 -0.05 0.9 0.9]);
    set(fig_indiv, 'InvertHardcopy', 'off');

    for k = 1:n_mice
        mouse_name = strrep(curr_mice{k}, '_', ' ');

        % Extract slab data
        raw_slab_d = curr_diff(z_indices, :, :, k);
        raw_slab_s = curr_sum(z_indices, :, :, k);

        % Extract mask slab
        mask_slab = curr_mask(z_indices, :, :, k);
        mask_slab_logical = logical(mask_slab);

        % Apply mask
        raw_slab_d(mask_slab_logical) = NaN;
        raw_slab_s(mask_slab_logical) = NaN;

        % Compute rolling median
        slab_diff_m = squeeze(nanmedian(raw_slab_d, 1)); %#ok<*NANMEDIAN>
        slab_sum_m  = squeeze(nanmedian(raw_slab_s, 1));

        % Prepare alpha mask for plotting
        slab_bg_m = squeeze(min(mask_slab, [], 1));

        % Combine with atlas mask
        valid_pixels = (atlas_slice > 0) & (~slab_bg_m);
        alpha_data = double(valid_pixels);

        % Row 1: diff plot
        subplot(2, n_mice, k);
        imagesc(slab_diff_m);
        set(findobj(gca, 'Type', 'image'), 'AlphaData', alpha_data);
        clim(clim_diff_indiv);
        colormap(gca, hot);
        axis image; axis off; set(gca, 'Color', 'k'); hold on;
        % Overlay atlas
        if bool_overlay_atlas
            plot(b_col, b_row, '.', 'Color', [0.5 0.5 0.5], 'MarkerSize', 0.1);
        end

        title(mouse_name, 'Color', 'w', 'FontSize', 10, 'Interpreter', 'none');
        % Colorbar for every plot
        cb = colorbar;
        cb.Label.String = '|L - R|';
        cb.Color = 'w'; cb.Label.Color = 'w';

        % Row 2: sum plot
        subplot(2, n_mice, k + n_mice);
        imagesc(slab_sum_m);
        set(findobj(gca, 'Type', 'image'), 'AlphaData', alpha_data);
        clim(clim_sum_indiv);
        colormap(gca, hot);
        axis image; axis off; set(gca, 'Color', 'k'); hold on;
        % Overlay atlas
        if bool_overlay_atlas
            plot(b_col, b_row, '.', 'Color', [0.5 0.5 0.5], 'MarkerSize', 0.1);
        end
        % Colorbar for every plot
        cb = colorbar;
        cb.Label.String = 'L + R';
        cb.Color = 'w'; cb.Label.Color = 'w';

    end
    sgtitle(['Individual Slab Median (' curr_grp ') - Slice ' num2str(target_slice) '\pm' num2str(slab_range)], ...
        'Color', 'w', 'FontSize', 16);

    % Save
    if exist('comp_out_dir', 'var')
        filename = sprintf('Indiv_Slab_Avg_%d_%s', target_slice, curr_grp);
        saveas(fig_indiv, fullfile(comp_out_dir, [filename '.fig']));
        exportgraphics(fig_indiv, fullfile(comp_out_dir, [filename '.png']), 'Resolution', 300, 'BackgroundColor', 'current');
        fprintf('Saved individual slab figure for %s.\n', curr_grp);
    end

end

clear group_names group_data_diff group_data_sum group_masks group_mice_list curr_grp curr_diff curr_sum curr_mask curr_mice

%% Generate rolling surprise-masked video (spatially smoothed by rolling median)

if generate_rolling_videos

    % Parameters
    slab_range = 10;
    surp_thresh = -log10(0.01);
    t_lim = [-6 6];

    % Generate filename
    vid_name = [['t_lr_diff_sum_' channel '_groupdiff_'] comp_tag '_surpmask_rolling_slab.mp4'];

    fprintf('Generating Rolling Slab Surpmask Video (Window: +/- %d)...\n', slab_range);

    write_lr_video_surpmask_rolling( ...
        t_lr_diff_groupdiff, ...                    % Diff Data
        t_lr_sum_groupdiff, ...                     % Sum Data
        half_atlas, ...                             % Atlas
        brainMask_group_diff, ...                   % Background Masks (1=Bkg)
        comp_out_dir, ...                           % Save Directory
        vid_name, ...                               % Filename
        t_lim, ...                                  % Color limits
        [exp_type ' - ' ctrl_type ' (p<0.01)'], ... % Title
        'T-Score Diff (Slab)', ...                  % Label Left
        'T-Score Sum (Slab)', ...                   % Label Right
        surp_diff, ...                              % Surprise Diff (for Left Alpha)
        surp_sum, ...                               % Surprise Sum (for Right Alpha)
        surp_thresh, ...                            % Threshold
        slab_range ...                              % Window
        );

end

%% Generate individual mice LR rolling videos (spatially smoothed by rolling median)

if generate_rolling_videos

    % Parameters
    slab_range = 10; % +/- 10 slices
    clim_diff = [0 2];
    clim_sum  = [0 10];

    % Use left-hemi atlas for boundaries
    if ~exist('atlas_left', 'var')
        atlas_left = double(AllenCrop(:, :, 1:size(lr_diff_ctrl,3)));
    end

    fprintf('Generating Individual LR Rolling Slab Videos...\n');

    % Control group
    write_lr_indiv_rolling_video( ...
        abs(lr_diff_ctrl), ...                              % Diff Data
        lr_sum_ctrl, ...                                    % Sum Data
        mask_bg_ctrl, ...                                   % Background Masks (1=Bkg)
        atlas_left, ...                                     % Atlas Vol
        comp_out_dir, ...                                   % Save Dir
        ['Individual_Rolling_Slab_' ctrl_type '.mp4'], ...  % Filename
        ctrl_type, ...                                      % Group Name
        ctrl_mousenames, ...                                % Mouse Names
        slab_range, ...                                     % Slab Radius
        clim_diff, clim_sum ...                             % Color Limits
        );

    % Experimental group
    write_lr_indiv_rolling_video( ...
        abs(lr_diff_exp), ...
        lr_sum_exp, ...
        mask_bg_exp, ...
        atlas_left, ...
        comp_out_dir, ...
        ['Individual_Rolling_Slab_' exp_type '.mp4'], ...
        exp_type, ...
        exp_mousenames, ...
        slab_range, ...
        clim_diff, clim_sum ...
        );

    fprintf('Individual rolling videos generation complete.\n');

end

%% Area-based Analysis: Regional Surprise Quantification (Rolling Median - Diff & Sum)

fprintf('Starting Regional Surprise Analysis (Rolling Median - Diff & Sum)...\n');

% Set parameters
slab_range = 10;                            % +/- slices for rolling median
p_thresh_agg = 0.01;                        % Threshold for aggregation count
surp_thresh_val = -log10(p_thresh_agg);

% Define target regions
roi_list_surp = { ...
    'Primary somatosensory area, barrel field', ...
    'Primary somatosensory area, trunk', ...
    'Primary somatosensory area, upper limb', ...
    'Primary somatosensory area, lower limb', ...
    'Supplemental somatosensory area', ...
    'Primary motor area', ...
    'Secondary motor area', ...
    'Primary visual area', ...
    'Lateral visual area', ...
    'Anterolateral visual area', ...
    'Anteromedial visual area', ...
    'Primary auditory area', ...
    'Dorsal auditory area', ...
    'Ventral auditory area', ...
    'Anterior cingulate area', ...
    'Olfactory tubercle', ...
    'Prelimbic area', ...
    'Infralimbic area', ...
    'Visceral area',...
    'Gustatory areas',...
    'Piriform area',...
    'Subiculum',...
    'Orbital area', ...
    'Claustrum', ...
    'Agranular insular area', ...
    'Anterior area', ...
    'Rostrolateral visual area', ...
    'Temporal association areas', ...
    'Perirhinal area', ...
    'Ectorhinal area', ...
    'Retrosplenial area', ...
    'Nucleus accumbens', ...
    'Caudoputamen', ...
    'Globus pallidus, external segment', ...
    'Subthalamic nucleus', ...
    'Ventral posteromedial nucleus of the thalamus', ...
    'Ventral posterolateral nucleus of the thalamus', ...
    'Ventral medial nucleus of the thalamus',...
    'Zona incerta',...
    'Posterior complex of the thalamus', ...
    'Lateral posterior nucleus of the thalamus', ...
    'Lateral dorsal nucleus of thalamus', ...
    'Ventral anterior-lateral complex of the thalamus', ...
    'Mediodorsal nucleus of the thalamus', ...
    'Parafascicular nucleus', ...
    'Nucleus of reuniens', ...
    'Central lateral nucleus of the thalamus', ...
    'Reticular nucleus of the thalamus', ...
    'Geniculate group, dorsal thalamus', ...
    'Midbrain, motor related', ...
    'Superior colliculus, motor related', ...
    'Superior colliculus, sensory related', ...
    'Hippocampal formation', ...
    'Basolateral amygdalar nucleus', ...
    'Parafascicular nucleus', ...
    'Hypothalamus' ...
    };
n_rois_surp = length(roi_list_surp);

% Setup atlas vectors
fprintf('  Mapping Atlas Volume...\n');
[~, ~, n_width] = size(AllenCrop);
half_width = floor(n_width / 2);
atlas_left = AllenCrop(:, :, 1:half_width);
valid_pixels = atlas_left > 0;
pixel_ids = atlas_left(valid_pixels);

% Pre-calculate masks
fprintf('  Building masks for %d regions...\n', n_rois_surp);
roi_masks_surp = false(length(pixel_ids), n_rois_surp);
roi_pixel_counts_surp = zeros(n_rois_surp, 1);

for r = 1:n_rois_surp
    region_name = roi_list_surp{r};
    try
        mask_temp = get_allen_region_mask(allenDir, AllenCrop, {region_name}, brainMask, '');
        if size(mask_temp, 3) >= half_width, mask_temp = mask_temp(:, :, 1:half_width); end
        ids_in_region = unique(atlas_left(mask_temp));
        ids_in_region(ids_in_region == 0) = [];
        roi_masks_surp(:, r) = ismember(pixel_ids, ids_in_region);
        roi_pixel_counts_surp(r) = sum(roi_masks_surp(:, r));
    catch
        warning('Could not map region: %s', region_name);
    end
end

% Initialize figure
fig_surp = figure('Visible', 'off', 'Name', ['Region_Surprise_BarChart_DiffSum_' comp_tag], ...
    'Color', 'w', 'Units', 'Normalized', 'Position', [0 0 0.9 0.9]);

% Loop for diff and sum
modes = {'Diff', 'Sum'};

for m_idx = 1:2
    mode_name = modes{m_idx};

    % Select input data
    if strcmp(mode_name, 'Diff')
        if ~exist('surp_diff', 'var'), error('surp_diff not found'); end
        raw_surp_vol = surp_diff;
    else
        if ~exist('surp_sum', 'var'), error('surp_sum not found'); end
        raw_surp_vol = surp_sum;
    end

    % Apply rolling median with masking
    fprintf('  [%s] Calculating rolling median (slab +/- %d)...\n', mode_name, slab_range);
    vol_surp = zeros(size(raw_surp_vol), 'single');
    n_slices = size(raw_surp_vol, 1);
    for z = 1:n_slices
        z_start = max(1, z - slab_range);
        z_end   = min(n_slices, z + slab_range);
        slab_data = raw_surp_vol(z_start:z_end, :, :);
        % Apply background mask
        if exist('brainMask_group_diff', 'var')
            slab_mask = brainMask_group_diff(z_start:z_end, :, :);
            slab_data(~slab_mask) = NaN;
        end
        vol_surp(z, :, :) = nanmedian(slab_data, 1);
    end

    % Flatten & aggregate
    surp_vec = vol_surp(valid_pixels);
    surp_vec(isnan(surp_vec)) = 0;
    roi_surp_agg = nan(n_rois_surp, 1);
    for r = 1:n_rois_surp
        if roi_pixel_counts_surp(r) == 0, continue; end
        vals = surp_vec(roi_masks_surp(:, r));
        if ~isempty(vals)
            roi_surp_agg(r) = sum(vals(vals > surp_thresh_val));
        end
    end

    % Sorting & plotting prep
    sort_metric = roi_surp_agg;
    [sorted_surp, sort_idx] = sort(sort_metric, 'ascend');
    sorted_rois_surp = roi_list_surp(sort_idx);

    % Filter valid
    valid_k = sorted_surp > 0 & ~isnan(sorted_surp);
    sorted_surp = sorted_surp(valid_k);
    sorted_rois_surp = sorted_rois_surp(valid_k);

    % Plot subplot
    subplot(1, 2, m_idx);

    b = barh(sorted_surp);
    b.FaceColor = 'flat';

    % Color logic
    c_map_surp = gray(256);
    c_map_surp = flipud(c_map_surp(1:200, :));
    if ~isempty(sorted_surp)
        c_vals = round( (sorted_surp / max(sorted_surp)) * size(c_map_surp,1) );
        c_vals(c_vals < 1) = 1; c_vals(isnan(c_vals)) = 1;
        for k = 1:length(sorted_surp), b.CData(k,:) = c_map_surp(c_vals(k), :); end
    end

    % Axis labels
    yticks(1:length(sorted_rois_surp));
    yticklabels(sorted_rois_surp);
    xlabel(['Aggregated surprise ( > ' num2str(surp_thresh_val, '%.1f') ')']);
    title(['Regional ' lower(mode_name) ' significance']);
    grid on; set(gca, 'FontSize', 10);
    ylim([0 length(sorted_rois_surp)+1]);

    % Highlighting
    switch exp_type
        case 'rws', highlighted_areas = {'Primary somatosensory area, barrel field', 'Ventral posteromedial nucleus of the thalamus', 'Posterior complex of the thalamus', 'Supplemental somatosensory area', 'Zona incerta', 'Rostrolateral visual area'};
        case 'behavior', highlighted_areas = {'Primary somatosensory area, barrel field', 'Ventral posteromedial nucleus of the thalamus', 'Posterior complex of the thalamus', 'Supplemental somatosensory area', 'Zona incerta', 'Rostrolateral visual area'};
        otherwise, highlighted_areas = {};
    end

    ax = gca;
    ytl = ax.YTickLabel;
    colored_labels = repmat({''}, size(ytl));
    for i = 1:length(ytl)
        if ismember(ytl{i}, highlighted_areas)
            colored_labels{i} = ['\color{magenta} \bf ' strrep(ytl{i}, '_', ' ') ];
        else
            colored_labels{i} = ['\color{black} ' strrep(ytl{i}, '_', ' ') ];
        end
    end
    ax.YTickLabel = colored_labels;
end

sgtitle(['Regional integrated significance (rolling median) - ' strrep(comp_tag,'_',' ')], 'FontSize', 14, 'FontWeight', 'bold');

% Save
set(fig_surp, 'InvertHardcopy', 'off');
saveas(fig_surp, fullfile(comp_out_dir, ['Region_Surprise_Bar_DiffSum_' comp_tag '.fig']));
exportgraphics(fig_surp, fullfile(comp_out_dir, ['Region_Surprise_Bar_DiffSum_' comp_tag '.png']), 'Resolution', 300);

fprintf('Regional surprise analysis (Diff & Sum) saved to: %s\n', comp_out_dir);
% clear roi_masks_surp surp_vec valid_pixels pixel_ids vol_surp