function group_differences(run_settings)
%GROUP_DIFFERENCES  Left-right differences of two groups, and their comparison.
%   GROUP_DIFFERENCES(run_settings) does the work of run_group_differences,
%   which sets the fields of run_settings and says what each one does.
%
%   Reads <channel>_4d_normalized.mat and <channel>_4d_normalized_bkgmask.mat
%   (run_normalise_groups) from the two groups' folders, and writes into
%   comp_out_dir, each as .fig and .png: Normalization_Profiles_LR_<comp_tag>
%   (the profile alignment), Slab_Avg_565_<comp_tag>_surpmask and
%   Indiv_Slab_Avg_565_<group> (the slab figures), and
%   Region_Surprise_Bar_DiffSum_<comp_tag> (the regional surprise bars); then
%   the videos switched on, and the region analyses if switched on.
%
%   The volumes are on the adults' CCF crop (planes 180 to 1079 of the 10 um
%   annotation). Each mouse's values outside its tissue (the atlas brain,
%   outside its background mask) are NaN, and its tissue is smoothed, when
%   apply_smoothing is set, by a normalised convolution: the Gaussian-smoothed
%   values over the Gaussian-smoothed tissue mask. The experimental group is put
%   on the control group's scale by the line between the two groups' mean plane
%   profiles, and both are divided by one common factor. Each mouse is folded
%   onto the left hemisphere: L - R and L + R of every voxel and its mirror
%   image, NaN where either side is outside the mouse's tissue. The group maps
%   are the means of the absolute values over the mice with a value, the t maps
%   Welch t of their difference,
%       t = (mean_exp - mean_ctrl) / sqrt(sem_ctrl^2 + sem_exp^2)
%   and the surprise -log10 of its two-sided p, with Welch-Satterthwaite
%   degrees of freedom. Every mean, SEM and degree of freedom counts, voxel by
%   voxel, only the mice with a value there.

% settings of run_group_differences, under the names the code below uses
paths = run_settings.paths;
ctrl_type = run_settings.ctrl_type;
exp_type = run_settings.exp_type;
ctrl_mousenames = run_settings.ctrl_mousenames;
exp_mousenames = run_settings.exp_mousenames;
behavior_mice = run_settings.behavior_mice;
generate_diff_videos = run_settings.generate_diff_videos;
generate_individual_diff_videos = run_settings.generate_individual_diff_videos;
generate_t_scored_videos = run_settings.generate_t_scored_videos;
generate_surprise_videos = run_settings.generate_surprise_videos;
generate_rolling_videos = run_settings.generate_rolling_videos;
generate_signed_diff_videos = run_settings.generate_signed_diff_videos;
perform_area_based_analysis_fine = run_settings.perform_area_based_analysis_fine;
perform_area_based_analysis_coarse = run_settings.perform_area_based_analysis_coarse;
apply_smoothing = run_settings.apply_smoothing;
smooth_sigma = run_settings.smooth_sigma;
channel = run_settings.channel;
comp_tag = run_settings.comp_tag;
ctrl_dir = run_settings.ctrl_dir;
exp_dir = run_settings.exp_dir;
comp_out_dir = run_settings.comp_out_dir;

%% Atlas

[allenDir, AllenCrop, brainMask, half_atlas] = load_allen_atlas(paths);

%% Load both groups

[data_4d_new_ctrl, data_4d_new_exp, med_data_4d_ctrl, recomputed_bkg_mask_4d_ctrl, ...
    med_data_4d_exp, recomputed_bkg_mask_4d_exp, exp_mousenames] = load_groups(channel, ...
    ctrl_type, exp_type, ctrl_dir, exp_dir, behavior_mice, exp_mousenames);

%% Each mouse's tissue, smoothed if asked

% NaN outside each mouse's tissue, so that no voxel outside it enters a mean, an
% SEM or a count of mice below
[data_4d_new_ctrl, data_4d_new_exp] = tissue_groups(data_4d_new_ctrl, ...
    data_4d_new_exp, recomputed_bkg_mask_4d_ctrl, recomputed_bkg_mask_4d_exp, ...
    brainMask, apply_smoothing, smooth_sigma);

%% Align the experimental group onto the control group

[interest_region, norm_ctrl, norm_exp, slope, intercept, norm_ctrl_med_fact, ...
    norm_exp_med_fact] = align_exp_to_ctrl(med_data_4d_ctrl, med_data_4d_exp);
plot_alignment_profiles(med_data_4d_ctrl, med_data_4d_exp, norm_ctrl, norm_exp, ...
    interest_region, slope, intercept, ctrl_type, exp_type, comp_tag, comp_out_dir);

%% Left-right maps of each mouse and group

[lr_diff_ctrl, lr_sum_ctrl, lr_diff_exp, lr_sum_exp, avg_lr_diff_ctrl, ...
    avg_lr_sum_ctrl, avg_lr_diff_exp, avg_lr_sum_exp, avg_lr_diff_groupdiff, ...
    avg_lr_sum_groupdiff] = compute_group_lr(data_4d_new_ctrl, data_4d_new_exp, ...
    norm_ctrl_med_fact, norm_exp_med_fact, slope, intercept);
[mask_bg_ctrl, mask_bg_exp, brainMask_cropped_no_bkg_ctrl, ...
    brainMask_cropped_no_bkg_exp, brainMask_group_diff] = hemisphere_masks(brainMask, ...
    lr_diff_ctrl, lr_diff_exp);

% videos of the group means and of their difference
if generate_diff_videos
    write_group_videos(avg_lr_diff_ctrl, avg_lr_sum_ctrl, avg_lr_diff_exp, ...
        avg_lr_sum_exp, avg_lr_diff_groupdiff, avg_lr_sum_groupdiff, half_atlas, ...
        brainMask_cropped_no_bkg_ctrl, brainMask_cropped_no_bkg_exp, ...
        brainMask_group_diff, comp_out_dir, channel, ctrl_type, exp_type, comp_tag);
end
clear data_4d_new_ctrl data_4d_new_exp recomputed_bkg_mask_4d_ctrl ...
    recomputed_bkg_mask_4d_exp

%% Videos of every mouse

if generate_individual_diff_videos
    write_individual_videos(lr_diff_ctrl, lr_sum_ctrl, mask_bg_ctrl, lr_diff_exp, ...
        lr_sum_exp, mask_bg_exp, AllenCrop, comp_out_dir, ctrl_type, exp_type, ...
        ctrl_mousenames, exp_mousenames);
end

%% Videos of every mouse, signed

if generate_signed_diff_videos
    write_signed_videos(lr_diff_ctrl, lr_sum_ctrl, mask_bg_ctrl, lr_diff_exp, ...
        lr_sum_exp, mask_bg_exp, AllenCrop, comp_out_dir, ctrl_type, exp_type, ...
        ctrl_mousenames, exp_mousenames);
end

%% Group t and surprise maps

[t_lr_diff_groupdiff, t_lr_sum_groupdiff, n_ctrl, n_exp, surp_diff, surp_sum] = ...
    group_t_and_surprise(lr_diff_ctrl, lr_sum_ctrl, lr_diff_exp, lr_sum_exp, ...
    avg_lr_diff_ctrl, avg_lr_sum_ctrl, avg_lr_diff_exp, avg_lr_sum_exp, ...
    avg_lr_diff_groupdiff, avg_lr_sum_groupdiff, half_atlas, brainMask_group_diff, ...
    comp_out_dir, channel, comp_tag, ctrl_type, exp_type, generate_t_scored_videos, ...
    generate_surprise_videos);

disp('Generalized LR analysis completed successfully!');
fprintf('All comparison results saved to: %s\n', comp_out_dir);

%% Region t-scores, every leaf region (off in production)

if perform_area_based_analysis_fine
    [atlas_left, half_width] = wholebrain_tmaps(AllenCrop, mask_bg_ctrl, mask_bg_exp, ...
        lr_diff_ctrl, lr_diff_exp, lr_sum_ctrl, lr_sum_exp, n_ctrl, n_exp, comp_tag, ...
        comp_out_dir);
end

%% Region t-scores, a fixed list of regions (off in production)

if perform_area_based_analysis_coarse

    % the width of the left hemisphere along ML, which the fine analysis returns;
    % without it, worked out here as wholebrain_tmaps does
    if ~perform_area_based_analysis_fine
        [~, ~, n_width] = size(AllenCrop);
        half_width = floor(n_width / 2);
    end

    coarse_region_tstats(AllenCrop, allenDir, brainMask, half_width, lr_diff_ctrl, ...
        lr_diff_exp, lr_sum_ctrl, lr_sum_exp, mask_bg_ctrl, mask_bg_exp, n_ctrl, ...
        n_exp, ctrl_mousenames, exp_mousenames, ctrl_type, exp_type, comp_tag, ...
        comp_out_dir);

    % cleared in this workspace for atlas_left, which the montage above returns
    % and the rolling videos below rebuild when it is missing
    clear data_vec roi_masks atlas_left input_vol_ctrl input_vol_exp vol vol_smoothed ...
        bg_vec roi_vals roi_vals_clean
    fprintf(' Area-based analysis saved to: %s\n', comp_out_dir);

end

%% Annotated video of the region t-scores (off in production)

if perform_area_based_analysis_fine

    write_annotated_tmap_video(AllenCrop, allenDir, comp_out_dir, comp_tag);

end

%% Slab figure of the group t-maps

% median over the planes around plane 565, shown where the surprise is high
plot_group_slab(t_lr_diff_groupdiff, t_lr_sum_groupdiff, surp_diff, surp_sum, ...
    brainMask_group_diff, half_atlas, exp_type, ctrl_type, comp_out_dir, comp_tag);

%% Slab figures of every mouse

plot_individual_slabs(lr_diff_ctrl, lr_diff_exp, lr_sum_ctrl, lr_sum_exp, ...
    mask_bg_ctrl, mask_bg_exp, ctrl_mousenames, exp_mousenames, ctrl_type, exp_type, ...
    half_atlas, comp_out_dir);

%% Rolling slab video of the group t-maps

if generate_rolling_videos

    write_rolling_tscore_video(t_lr_diff_groupdiff, t_lr_sum_groupdiff, half_atlas, ...
        brainMask_group_diff, surp_diff, surp_sum, comp_out_dir, channel, comp_tag, ...
        exp_type, ctrl_type);

end

%% Rolling slab videos of every mouse

if generate_rolling_videos

    % median over +/- 10 planes; colour limits of the difference and the sum
    slab_range = 10;
    clim_diff = [0 2];
    clim_sum = [0 10];

    % the left half of the atlas, for the boundaries, unless the leaf region
    % analysis left it
    if ~exist('atlas_left', 'var')
        atlas_left = double(AllenCrop(:, :, 1:size(lr_diff_ctrl, 3)));
    end

    write_individual_rolling_videos(lr_diff_ctrl, lr_sum_ctrl, mask_bg_ctrl, ...
        lr_diff_exp, lr_sum_exp, mask_bg_exp, atlas_left, comp_out_dir, ctrl_type, ...
        exp_type, ctrl_mousenames, exp_mousenames, slab_range, clim_diff, clim_sum);

end

%% Regional surprise bars

regional_surprise_bars(surp_diff, surp_sum, brainMask_group_diff, AllenCrop, allenDir, ...
    brainMask, exp_type, comp_tag, comp_out_dir);

end

% ===== Local functions: loading and alignment =====

function [allenDir, AllenCrop, brainMask, half_atlas] = load_allen_atlas(paths)
% The Allen annotation cropped to the AP range of the registered volumes, its
% brain mask, and the atlas the videos draw on.

% the 10 um annotation, cropped to the adults' AP planes 180 to 1079
allenDir = paths.atlas;
addpath(allenDir);
AllenFile = fullfile(allenDir, 'annotation_10.nii.gz');
AllenVol = niftiread(AllenFile);
limits = [180 1079];
AllenCrop = AllenVol(limits(1):limits(2), :, :);
brainMask = AllenCrop > 0;

% the atlas the videos draw their boundaries from: the whole width, despite the
% name (the videos show the left half)
half_atlas = AllenCrop(:, :, 1:end);
clear bg_L_c bg_R_c bg_L_e bg_R_e AllenVol
end

function [data_4d_new_ctrl, data_4d_new_exp, med_data_4d_ctrl, ...
    recomputed_bkg_mask_4d_ctrl, med_data_4d_exp, recomputed_bkg_mask_4d_exp, ...
    exp_mousenames] = load_groups(channel, ctrl_type, exp_type, ctrl_dir, exp_dir, ...
    behavior_mice, exp_mousenames)
% Both groups' normalised volumes and background masks, and each mouse's mean
% tissue intensity per plane; for behavior, the mice named in behavior_mice, which
% then name the experimental mice.

% the files of the channel, and the name of the volume inside them
norm_var_name = [channel '_4d_normalized'];
norm_filename = [channel '_4d_normalized.mat'];
bkgmask_filename = [channel '_4d_normalized_bkgmask.mat'];

% the control group: every mouse saved, since run_normalise_groups kept the
% selected ones (the selection here is off:
% (:,:,:,selected_mice_idx_list{ctrl_group_idx}))
fprintf('Loading Control group data (%s, channel=%s)...\n', ctrl_type, channel);
S_ctrl_vol = load(fullfile(ctrl_dir, norm_filename), norm_var_name);
S_ctrl_mask = load(fullfile(ctrl_dir, bkgmask_filename), 'recomputed_bkg_mask_4d');
data_4d_new_ctrl = S_ctrl_vol.(norm_var_name);
data_4d_new_ctrl_bkgmask = S_ctrl_mask.recomputed_bkg_mask_4d;
clear S_ctrl_vol S_ctrl_mask
fprintf('  Kept %d mice based on selection.\n', size(data_4d_new_ctrl, 4));

% the experimental group: every mouse saved for rws (the selection is off, as
% above), the mice behavior_mice names for behavior
fprintf('Loading Experimental group data (%s, channel=%s)...\n', exp_type, channel);
S_exp_vol = load(fullfile(exp_dir, norm_filename), norm_var_name);
S_exp_mask = load(fullfile(exp_dir, bkgmask_filename), 'recomputed_bkg_mask_4d');
if strcmp(exp_type, 'rws')
    data_4d_new_exp = S_exp_vol.(norm_var_name);
    data_4d_new_exp_bkgmask = S_exp_mask.recomputed_bkg_mask_4d;
elseif strcmp(exp_type, 'behavior')
    behavior_idx = saved_positions(behavior_mice, fullfile(exp_dir, norm_filename));
    data_4d_new_exp = S_exp_vol.(norm_var_name)(:, :, :, behavior_idx);
    data_4d_new_exp_bkgmask = S_exp_mask.recomputed_bkg_mask_4d(:, :, :, behavior_idx);
    exp_mousenames = behavior_mice;
else
    error('run_group_differences: unknown exp_type ''%s'' (use ''rws'' or ''behavior'').', ...
        exp_type);
end
clear S_exp_vol S_exp_mask
fprintf('  Kept %d mice based on selection.\n', size(data_4d_new_exp, 4));

% each control mouse's mean tissue intensity per plane, on the saved masks
fprintf(['Plane means of the tissue, on the saved background masks, controls ' ...
         '(%s)...\n'], ctrl_type);
med_data_4d_ctrl = plane_tissue_means(data_4d_new_ctrl, data_4d_new_ctrl_bkgmask);

% the saved masks are the ones used from here on
recomputed_bkg_mask_4d_ctrl = data_4d_new_ctrl_bkgmask;
clear data_4d_new_ctrl_bkgmask

% the same for each experimental mouse
fprintf(['Plane means of the tissue, on the saved background masks, experimentals ' ...
         '(%s)...\n'], exp_type);
med_data_4d_exp = plane_tissue_means(data_4d_new_exp, data_4d_new_exp_bkgmask);
recomputed_bkg_mask_4d_exp = data_4d_new_exp_bkgmask;
clear data_4d_new_exp_bkgmask
end

function mouse_idx = saved_positions(mouse_names, norm_file)
% The positions of the named mice along the fourth dimension of a normalised
% volume, from the mouse names run_normalise_groups saved with it.

S_saved = load(norm_file, 'current_mice');
[is_saved, mouse_idx] = ismember(mouse_names, S_saved.current_mice);
if ~all(is_saved)
    error(['run_group_differences: %s not among the mice run_normalise_groups ' ...
           'saved in %s (%s). Name some of those in behavior_mice, or select the ' ...
           'mice in run_normalise_groups and rerun it.'], ...
          strjoin(mouse_names(~is_saved), ', '), norm_file, ...
          strjoin(S_saved.current_mice, ', '));
end
end

function med_data_4d = plane_tissue_means(data_4d, bkgmask_4d)
% Each mouse's mean tissue intensity per plane (planes x mice), outside its
% background mask.

med_data_4d = nan(size(data_4d, 1), size(data_4d, 4));
total_slices = size(data_4d, 1);
for iii = 1:size(data_4d, 4)
    fprintf('  Processing Mouse %d ...\n', iii);
    for slice_idx_loop = 1:total_slices
        img_data = squeeze(data_4d(slice_idx_loop, :, :, iii));
        bg_mask = squeeze(bkgmask_4d(slice_idx_loop, :, :, iii));
        med_data_4d(slice_idx_loop, iii) = nanmean(img_data(~bg_mask));
    end
end
end

function [data_4d_new_ctrl, data_4d_new_exp] = tissue_groups(data_4d_new_ctrl, ...
    data_4d_new_exp, recomputed_bkg_mask_4d_ctrl, recomputed_bkg_mask_4d_exp, ...
    brainMask, apply_smoothing, smooth_sigma)
% Each mouse's volume NaN outside its tissue, and its tissue smoothed in 3D when
% apply_smoothing is set. The volumes come back under the names they came in
% with, so MATLAB updates them in place.

fprintf('Setting each mouse''s voxels outside its tissue to NaN...\n');
if apply_smoothing
    fprintf('Applying NaN-Robust 3D Gaussian Smoothing (Sigma = %.1f)...\n', smooth_sigma);
end

% the control group
fprintf('  Processing Control group...\n');
for i = 1:size(data_4d_new_ctrl, 4)
    tic
    data_4d_new_ctrl(:, :, :, i) = tissue_only(data_4d_new_ctrl(:, :, :, i), ...
        recomputed_bkg_mask_4d_ctrl(:, :, :, i), brainMask, apply_smoothing, ...
        smooth_sigma);
    toc
end

% the same for the experimental group
fprintf('  Processing Experimental group...\n');
for i = 1:size(data_4d_new_exp, 4)
    tic
    data_4d_new_exp(:, :, :, i) = tissue_only(data_4d_new_exp(:, :, :, i), ...
        recomputed_bkg_mask_4d_exp(:, :, :, i), brainMask, apply_smoothing, ...
        smooth_sigma);
    toc
end

fprintf('  Done.\n');
end

function vol = tissue_only(vol, bg_mask, brainMask, apply_smoothing, smooth_sigma)
% One mouse's volume, NaN outside its tissue (in the atlas brain, outside its
% background, with a value), the tissue smoothed when apply_smoothing is set.

% the tissue: in the atlas brain, outside the background, reached by a section
% (run_normalise_groups left NaN where none was)
is_tissue = brainMask & ~bg_mask & ~isnan(vol);
vol(~is_tissue) = NaN;

% normalised convolution, both with the same Gaussian: the smoothed values (NaN
% as 0) over the smoothed tissue mask, so the voxels outside the tissue neither
% count as zero nor spread into it; outside the tissue it stays NaN, since a 0
% there would enter the group means and the counts of mice as a measured value
if apply_smoothing
    smoothed_values = imgaussfilt3(fillmissing(vol, 'constant', 0), smooth_sigma);
    smoothed_mask = imgaussfilt3(double(is_tissue), smooth_sigma);
    vol = smoothed_values ./ smoothed_mask;
    vol(~is_tissue) = NaN;
end
end

function [interest_region, norm_ctrl, norm_exp, slope, intercept, norm_ctrl_med_fact, ...
    norm_exp_med_fact] = align_exp_to_ctrl(med_data_4d_ctrl, med_data_4d_exp)
% The line that maps the experimental group's mean plane profile onto the
% control group's, and the common scale of both groups.

% planes the line is fitted on
interest_region = 200:700;

% the groups' mean plane profiles
mean_profile_ctrl = nanmean(med_data_4d_ctrl, 2);
mean_profile_exp = nanmean(med_data_4d_exp, 2);

% the line from the experimental profile to the control one, on the fit planes
% where both have a value
y_target = mean_profile_ctrl(interest_region);
x_source = mean_profile_exp(interest_region);
valid_idx = ~isnan(x_source) & ~isnan(y_target);
x_source = x_source(valid_idx);
y_target = y_target(valid_idx);
p = polyfit(x_source, y_target, 1);
slope = p(1);
intercept = p(2);

fprintf('Alignment Parameters (Exp -> Ctrl): Slope = %.4f, Intercept = %.4f\n', slope, ...
    intercept);

% the control profiles stay as they are, the experimental ones go through the line
norm_ctrl = med_data_4d_ctrl;
norm_exp = (med_data_4d_exp .* slope) + intercept;

% one common scale for both groups: the average of the two groups' mean
% intensity over planes 300-500, after the alignment
interest_region_bis = 300:500;
norm_ctrl_med_fact = nanmean(nanmean(med_data_4d_ctrl(interest_region_bis, :), 1));
norm_exp_med_fact = nanmean(nanmean(norm_exp(interest_region_bis, :), 1));
avg_med_fact = (norm_ctrl_med_fact + norm_exp_med_fact)./2;
norm_ctrl_med_fact = avg_med_fact;
norm_exp_med_fact = avg_med_fact;
end

function plot_alignment_profiles(med_data_4d_ctrl, med_data_4d_exp, norm_ctrl, ...
    norm_exp, interest_region, slope, intercept, ctrl_type, exp_type, comp_tag, ...
    comp_out_dir)
% The figure of both groups' plane profiles, before and after the alignment.

fig_norm = figure('Visible', 'off', 'Name', ['Normalization_Profile_LR_' comp_tag], ...
    'Position', [100, 100, 1200, 500], 'Color', 'w');

% five shades of blue for the control mice, five of red for the experimental ones
cmap_ctrl = [linspace(0.6, 0, 5)', linspace(0.7, 0.2, 5)', linspace(1, 0.4, 5)'];
cmap_exp = [linspace(1, 0.5, 5)', linspace(0.6, 0.1, 5)', linspace(0.6, 0.1, 5)'];

% a group of more than five mice takes the shades again from the first
n_rows = max([5, size(med_data_4d_ctrl, 2), size(med_data_4d_exp, 2)]);
cmap_ctrl = cmap_ctrl(mod(0:n_rows - 1, 5) + 1, :);
cmap_exp = cmap_exp(mod(0:n_rows - 1, 5) + 1, :);
slices = 1:size(med_data_4d_ctrl, 1);

% left: the raw profiles
draw_raw_profiles(med_data_4d_ctrl, slices, cmap_ctrl, med_data_4d_exp, cmap_exp);

% the fit planes
xline(interest_region(1), 'k--');
xline(interest_region(end), 'k--');
set(gca, 'FontSize', 12);

% right: the aligned profiles
draw_aligned_profiles(norm_ctrl, slices, cmap_ctrl, norm_exp, cmap_exp, slope, ...
    intercept);
xline(interest_region(1), 'k--');
xline(interest_region(end), 'k--');
set(gca, 'FontSize', 12);
sgtitle(['Profile Alignment: ' ctrl_type ' (Ref) vs ' exp_type ' (Aligned)'], ...
    'FontSize', 14, 'Color', 'k', 'Interpreter', 'none');

% save it
set(fig_norm, 'InvertHardcopy', 'off');
saveas(fig_norm, ...
    fullfile(comp_out_dir, ['Normalization_Profiles_LR_' comp_tag '.fig']));
exportgraphics(fig_norm, ...
    fullfile(comp_out_dir, ['Normalization_Profiles_LR_' comp_tag '.png']), ...
    'Resolution', 300);
fprintf('Saved LR Normalization Profile plot to: %s\n', comp_out_dir);
end

function draw_raw_profiles(med_data_4d_ctrl, slices, cmap_ctrl, med_data_4d_exp, ...
    cmap_exp)
% Every mouse's raw plane profile, and each group's mean and SEM.

% left: the raw profiles of every mouse, and each group's mean and SEM
subplot(1, 2, 1);
hold on;
box on;
grid on;
for i = 1:size(med_data_4d_ctrl, 2)
    plot(slices, med_data_4d_ctrl(:, i), 'Color', cmap_ctrl(i, :), 'LineWidth', 1.2);
end
for i = 1:size(med_data_4d_exp, 2)
    plot(slices, med_data_4d_exp(:, i), 'Color', cmap_exp(i, :), 'LineWidth', 1.2);
end

% each group's mean and SEM, over the mice with tissue on the plane
mean_c = nanmean(med_data_4d_ctrl, 2);
sem_c = nanstd(med_data_4d_ctrl, [], 2) ./ sqrt(sum(~isnan(med_data_4d_ctrl), 2));
mean_e = nanmean(med_data_4d_exp, 2);
sem_e = nanstd(med_data_4d_exp, [], 2) ./ sqrt(sum(~isnan(med_data_4d_exp), 2));
fill([slices fliplr(slices)], [mean_c-sem_c; flipud(mean_c+sem_c)], ...
    sep_palette('control'), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(slices, mean_c, 'Color', sep_palette('control_mean'), 'LineWidth', 3.5);
fill([slices fliplr(slices)], [mean_e-sem_e; flipud(mean_e+sem_e)], ...
    sep_palette('experimental'), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(slices, mean_e, 'Color', sep_palette('experimental_mean'), 'LineWidth', 3.5);
title('Raw Nanobody intensity', 'FontSize', 12);
xlabel('Coronal index', 'FontSize', 12);
ylabel('Intensity', 'FontSize', 12);
end

function draw_aligned_profiles(norm_ctrl, slices, cmap_ctrl, norm_exp, cmap_exp, ...
    slope, intercept)
% The same after the alignment, titled with the line.

% right: the same after the alignment
subplot(1, 2, 2);
hold on;
box on;
grid on;
for i = 1:size(norm_ctrl, 2)
    plot(slices, norm_ctrl(:, i), 'Color', cmap_ctrl(i, :), 'LineWidth', 1.2);
end
for i = 1:size(norm_exp, 2)
    plot(slices, norm_exp(:, i), 'Color', cmap_exp(i, :), 'LineWidth', 1.2);
end

% each group's mean and SEM, over the mice with tissue on the plane
mean_nc = nanmean(norm_ctrl, 2);
sem_nc = nanstd(norm_ctrl, [], 2) ./ sqrt(sum(~isnan(norm_ctrl), 2));
mean_ne = nanmean(norm_exp, 2);
sem_ne = nanstd(norm_exp, [], 2) ./ sqrt(sum(~isnan(norm_exp), 2));
fill([slices fliplr(slices)], [mean_nc-sem_nc; flipud(mean_nc+sem_nc)], ...
    sep_palette('control'), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(slices, mean_nc, 'Color', sep_palette('control_mean'), 'LineWidth', 3.5);
fill([slices fliplr(slices)], [mean_ne-sem_ne; flipud(mean_ne+sem_ne)], ...
    sep_palette('experimental'), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(slices, mean_ne, 'Color', sep_palette('experimental_mean'), 'LineWidth', 3.5);
title(sprintf('Linearly Aligned (Slope=%.2f, Int=%.0f)', slope, intercept), ...
    'FontSize', 12);
xlabel('Coronal index', 'FontSize', 12);
ylabel('Aligned Intensity', 'FontSize', 12);
end

% ===== Local functions: left-right maps and videos =====

function [lr_diff_ctrl, lr_sum_ctrl, lr_diff_exp, lr_sum_exp, avg_lr_diff_ctrl, ...
    avg_lr_sum_ctrl, avg_lr_diff_exp, avg_lr_sum_exp, avg_lr_diff_groupdiff, ...
    avg_lr_sum_groupdiff] = compute_group_lr(data_4d_new_ctrl, data_4d_new_exp, ...
    norm_ctrl_med_fact, norm_exp_med_fact, slope, intercept)
% Each mouse's left-right difference and sum on the aligned scale, the group
% means of their absolute values, and the difference of those means.

% each mouse's L - R and L + R on the common scale, the experimental group
% through the alignment line first
[lr_diff_ctrl, lr_sum_ctrl] = compute_lr_stats(data_4d_new_ctrl./norm_ctrl_med_fact);
[lr_diff_exp, lr_sum_exp] = compute_lr_stats( ...
    ((data_4d_new_exp .* slope) + intercept)./norm_exp_med_fact);

% the group means of the absolute values, and experimental minus control
avg_lr_diff_ctrl = nanmean(abs(lr_diff_ctrl), 4); %#ok<*NANMEAN>
avg_lr_sum_ctrl = nanmean(abs(lr_sum_ctrl), 4);
avg_lr_diff_exp = nanmean(abs(lr_diff_exp), 4);
avg_lr_sum_exp = nanmean(abs(lr_sum_exp), 4);
avg_lr_diff_groupdiff = avg_lr_diff_exp - avg_lr_diff_ctrl;
avg_lr_sum_groupdiff = avg_lr_sum_exp - avg_lr_sum_ctrl;
end

function [mask_bg_ctrl, mask_bg_exp, brainMask_cropped_no_bkg_ctrl, ...
    brainMask_cropped_no_bkg_exp, brainMask_group_diff] = hemisphere_masks(brainMask, ...
    lr_diff_ctrl, lr_diff_exp)
% The folded voxels without a left-right value: per mouse, and the voxels shown
% for each group and for the group difference.

% the brain mask over the folded width
brainMask_cropped = brainMask(:, :, 1:size(lr_diff_ctrl, 3));

% each mouse's background, folded: the voxels without a value, outside its
% tissue on either side (so a figure never draws a NaN as a colour)
mask_bg_ctrl = isnan(lr_diff_ctrl);
mask_bg_exp = isnan(lr_diff_exp);

% each group's mask: the brain where at least one mouse has tissue; the group
% difference's: where both groups do
tissue_3d_ctrl = any(~mask_bg_ctrl, 4);
tissue_3d_exp = any(~mask_bg_exp, 4);
brainMask_cropped_no_bkg_ctrl = brainMask_cropped & tissue_3d_ctrl;
brainMask_cropped_no_bkg_exp = brainMask_cropped & tissue_3d_exp;
brainMask_group_diff = brainMask_cropped_no_bkg_ctrl & brainMask_cropped_no_bkg_exp;
end

function write_group_videos(avg_lr_diff_ctrl, avg_lr_sum_ctrl, avg_lr_diff_exp, ...
    avg_lr_sum_exp, avg_lr_diff_groupdiff, avg_lr_sum_groupdiff, half_atlas, ...
    brainMask_cropped_no_bkg_ctrl, brainMask_cropped_no_bkg_exp, brainMask_group_diff, ...
    comp_out_dir, channel, ctrl_type, exp_type, comp_tag)
% Videos of each group's mean left-right difference and sum, and of the group
% difference.

% the control group's means
write_lr_video(avg_lr_diff_ctrl, avg_lr_sum_ctrl, half_atlas, ...
    brainMask_cropped_no_bkg_ctrl, comp_out_dir, ...
    [['lr_diff_sum_' channel '_'] ctrl_type '.mp4'], [0, 2.5], ctrl_type, ...
    ['LR abs difference (' channel ') average - ' ctrl_type], ...
    ['LR abs sum (' channel ') average - ' ctrl_type]);

% the experimental group's means
write_lr_video(avg_lr_diff_exp, avg_lr_sum_exp, half_atlas, ...
    brainMask_cropped_no_bkg_exp, comp_out_dir, ...
    [['lr_diff_sum_' channel '_'] exp_type '.mp4'], [0, 2.5], exp_type, ...
    ['LR abs difference (' channel ') average - ' exp_type], ...
    ['LR abs sum (' channel ') average - ' exp_type]);

% their difference
write_lr_video(avg_lr_diff_groupdiff, avg_lr_sum_groupdiff, half_atlas, ...
    brainMask_group_diff, comp_out_dir, ...
    ['lr_diff_sum_' channel '_groupdiff_' comp_tag '.mp4'], [-2.5, 2.5], ...
    [exp_type ' - ' ctrl_type], ...
    ['LR abs diff groupdiff (' comp_tag ')'], ['LR abs sum groupdiff (' comp_tag ')']);
end

function write_individual_videos(lr_diff_ctrl, lr_sum_ctrl, mask_bg_ctrl, lr_diff_exp, ...
    lr_sum_exp, mask_bg_exp, AllenCrop, comp_out_dir, ctrl_type, exp_type, ...
    ctrl_mousenames, exp_mousenames)
% Videos of every mouse's absolute left-right difference and sum, one per group.

fprintf('Generating Individual LR-Asymmetry Videos...\n');

% colour limits of the difference and the sum
clim_diff = [0 1.5];
clim_sum = [0 10];

% the control group: the absolute difference, the sum, the background masks and
% the whole atlas (the video crops it to the hemisphere)
write_lr_indiv_video( ...
    abs(lr_diff_ctrl), ...
    lr_sum_ctrl, ...
    mask_bg_ctrl, ...
    AllenCrop, ...
    comp_out_dir, ...
    ['Individual_LR_Diff_Sum_' ctrl_type '.mp4'], ...
    clim_diff, clim_sum, ...
    ctrl_type, ...
    ctrl_mousenames ...
    );

% the experimental group
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

function write_signed_videos(lr_diff_ctrl, lr_sum_ctrl, mask_bg_ctrl, lr_diff_exp, ...
    lr_sum_exp, mask_bg_exp, AllenCrop, comp_out_dir, ctrl_type, exp_type, ...
    ctrl_mousenames, exp_mousenames)
% Videos of every mouse's signed left-right difference and sum, one per group.

fprintf('Generating Signed (Directional) LR-Asymmetry Videos...\n');

% colour limits of the signed difference and the sum
clim_diff_signed = [-0.75 0.75];
clim_sum = [0 10];

% the control group
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

% the experimental group
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

function [t_lr_diff_groupdiff, t_lr_sum_groupdiff, n_ctrl, n_exp, surp_diff, ...
    surp_sum] = group_t_and_surprise(lr_diff_ctrl, lr_sum_ctrl, lr_diff_exp, ...
    lr_sum_exp, avg_lr_diff_ctrl, avg_lr_sum_ctrl, avg_lr_diff_exp, avg_lr_sum_exp, ...
    avg_lr_diff_groupdiff, avg_lr_sum_groupdiff, half_atlas, brainMask_group_diff, ...
    comp_out_dir, channel, comp_tag, ctrl_type, exp_type, generate_t_scored_videos, ...
    generate_surprise_videos)
% Welch t and surprise (-log10 p) maps of the group difference, and their videos.
% One function, since the surprise video also draws the t maps, within t_lim.

% the number of mice of each group, for the region analyses
n_ctrl = size(lr_diff_ctrl, 4);
n_exp = size(lr_diff_exp, 4);

% each voxel's number of mice with a value in each group (in single, to spare
% memory); the sum has the same missing voxels as the difference
n_vox_ctrl = single(sum(~isnan(lr_diff_ctrl), 4));
n_vox_exp = single(sum(~isnan(lr_diff_exp), 4));

% the voxels with a t: shown for both groups, with at least two mice with a value
% in each (elsewhere the t and its surprise are NaN, which a video would draw in
% its colormap's first colour)
has_t = brainMask_group_diff & n_vox_ctrl >= 2 & n_vox_exp >= 2;

% the SEMs and the Welch t of the group difference
[sem_lr_diff_ctrl, sem_lr_sum_ctrl, sem_lr_diff_exp, sem_lr_sum_exp, ...
    t_lr_diff_groupdiff, t_lr_sum_groupdiff] = group_welch_t(lr_diff_ctrl, ...
    lr_sum_ctrl, lr_diff_exp, lr_sum_exp, n_vox_ctrl, n_vox_exp, avg_lr_diff_ctrl, ...
    avg_lr_sum_ctrl, avg_lr_diff_exp, avg_lr_sum_exp, avg_lr_diff_groupdiff, ...
    avg_lr_sum_groupdiff);

% the colour limits of the t maps, in both videos below
t_lim = [-6 6];

% its video, shown on the voxels with a t (not on the whole brain mask,
% brainMask_cropped)
if generate_t_scored_videos
    write_lr_video(t_lr_diff_groupdiff, t_lr_sum_groupdiff, half_atlas, ...
        has_t, comp_out_dir, ...
        [['t_lr_diff_sum_' channel '_groupdiff_'] comp_tag '.mp4'], t_lim, ...
        [exp_type ' - ' ctrl_type ' (t-score)'], ...
        ['LR abs diff t-score (' comp_tag ')'], ['LR abs sum t-score (' comp_tag ')']);
end

% the surprise of the t maps
[surp_diff, surp_sum] = welch_surprise(n_vox_ctrl, n_vox_exp, sem_lr_diff_ctrl, ...
    sem_lr_diff_exp, sem_lr_sum_ctrl, sem_lr_sum_exp, t_lr_diff_groupdiff, ...
    t_lr_sum_groupdiff);

% the surprise video, and the t video shown where p < 0.05 (as above, on the
% voxels with a t)
if generate_surprise_videos
    write_lr_video(surp_diff, surp_sum, half_atlas, has_t, ...
        comp_out_dir, ['surp_lr_diff_sum_' channel '_groupdiff_' comp_tag '.mp4'], ...
        [0 8], ['-log_{10}(p) | ' comp_tag], ...
        ['LR abs diff surprise (' comp_tag ')'], ['LR abs sum surprise (' comp_tag ')']);
    surp_thresh = -log10(0.05);
    write_lr_video_surpmask( ...
        t_lr_diff_groupdiff, t_lr_sum_groupdiff, ...
        half_atlas, has_t, ...
        comp_out_dir, ...
        [['t_lr_diff_sum_' channel '_groupdiff_'] comp_tag '_surpmask.mp4'], ...
        t_lim, [exp_type ' - ' ctrl_type ' (t-score, p<0.05)'], ...
        ['LR abs diff t-score masked (' comp_tag ')'], ...
        ['LR abs sum t-score masked (' comp_tag ')'], ...
        surp_diff, surp_thresh);
end
end

function [sem_lr_diff_ctrl, sem_lr_sum_ctrl, sem_lr_diff_exp, sem_lr_sum_exp, ...
    t_lr_diff_groupdiff, t_lr_sum_groupdiff] = group_welch_t(lr_diff_ctrl, ...
    lr_sum_ctrl, lr_diff_exp, lr_sum_exp, n_vox_ctrl, n_vox_exp, avg_lr_diff_ctrl, ...
    avg_lr_sum_ctrl, avg_lr_diff_exp, avg_lr_sum_exp, avg_lr_diff_groupdiff, ...
    avg_lr_sum_groupdiff)
% Each group's SEM of the absolute values over each voxel's mice with a value
% (zero to NaN), and the Welch t of the group difference.

% each group's SEM of the absolute values, over the mice with a value (the
% counts in double, as before, so that a voxel every mouse has keeps its old
% value to the last bit)
sem_lr_diff_ctrl = nanstd(abs(lr_diff_ctrl), [], 4) ./ ...
    sqrt(double(n_vox_ctrl)); %#ok<*NANSTD>
sem_lr_sum_ctrl = nanstd(abs(lr_sum_ctrl), [], 4) ./ sqrt(double(n_vox_ctrl));
sem_lr_diff_exp = nanstd(abs(lr_diff_exp), [], 4) ./ sqrt(double(n_vox_exp));
sem_lr_sum_exp = nanstd(abs(lr_sum_exp), [], 4) ./ sqrt(double(n_vox_exp));

% a zero SEM to NaN, so the t is NaN rather than infinite (one mouse gives a zero
% SEM, so a voxel needs two in each group)
sem_lr_diff_ctrl(sem_lr_diff_ctrl==0) = NaN;
sem_lr_sum_ctrl(sem_lr_sum_ctrl==0) = NaN;
sem_lr_diff_exp(sem_lr_diff_exp==0) = NaN;
sem_lr_sum_exp(sem_lr_sum_exp==0) = NaN;

% the t of each group's mean against zero
t_lr_diff_ctrl = avg_lr_diff_ctrl ./ sem_lr_diff_ctrl;
t_lr_sum_ctrl = avg_lr_sum_ctrl ./ sem_lr_sum_ctrl;
t_lr_diff_exp = avg_lr_diff_exp ./ sem_lr_diff_exp;
t_lr_sum_exp = avg_lr_sum_exp ./ sem_lr_sum_exp;

% the SEM of the group difference
sem_diff_lr_diff = sqrt(sem_lr_diff_ctrl.^2 + sem_lr_diff_exp.^2);
sem_diff_lr_sum = sqrt(sem_lr_sum_ctrl.^2 + sem_lr_sum_exp.^2);
sem_diff_lr_diff(isnan(sem_diff_lr_diff) | sem_diff_lr_diff==0) = NaN;
sem_diff_lr_sum(isnan(sem_diff_lr_sum) | sem_diff_lr_sum==0) = NaN;

% the Welch t of the group difference
t_lr_diff_groupdiff = avg_lr_diff_groupdiff ./ sem_diff_lr_diff;
t_lr_sum_groupdiff = avg_lr_sum_groupdiff ./ sem_diff_lr_sum;
end

function [surp_diff, surp_sum] = welch_surprise(n_vox_ctrl, n_vox_exp, ...
    sem_lr_diff_ctrl, sem_lr_diff_exp, sem_lr_sum_ctrl, sem_lr_sum_exp, ...
    t_lr_diff_groupdiff, t_lr_sum_groupdiff)
% The Welch-Satterthwaite degrees of freedom from each voxel's numbers of mice
% with a value, and the surprise -log10 p of the t maps.

% Welch-Satterthwaite degrees of freedom, voxel by voxel (the counts in double,
% as the SEMs')
var1_diff = sem_lr_diff_ctrl.^2;
var2_diff = sem_lr_diff_exp.^2;
var1_sum = sem_lr_sum_ctrl.^2;
var2_sum = sem_lr_sum_exp.^2;
df_diff = (var1_diff + var2_diff).^2 ./ ...
    (var1_diff.^2 ./ (double(n_vox_ctrl) - 1) + var2_diff.^2 ./ (double(n_vox_exp) - 1));
df_sum = (var1_sum + var2_sum).^2 ./ ...
    (var1_sum.^2 ./ (double(n_vox_ctrl) - 1) + var2_sum.^2 ./ (double(n_vox_exp) - 1));

% NaN where a group has fewer than two mice with a value
too_few_mice = n_vox_ctrl < 2 | n_vox_exp < 2;
df_diff(too_few_mice) = NaN;
df_sum(too_few_mice) = NaN;

% the two-sided p of the Welch t, and the surprise -log10 p
p_diff = 2 * tcdf(-abs(t_lr_diff_groupdiff), df_diff);
p_sum = 2 * tcdf(-abs(t_lr_sum_groupdiff), df_sum);
surp_diff = -log10(p_diff);
surp_sum = -log10(p_sum);
end

% ===== Local functions: region analyses (off in production) =====

function [atlas_left, half_width] = wholebrain_tmaps(AllenCrop, mask_bg_ctrl, ...
    mask_bg_exp, lr_diff_ctrl, lr_diff_exp, lr_sum_ctrl, lr_sum_exp, n_ctrl, n_exp, ...
    comp_tag, comp_out_dir)
% Group t-scores per leaf region (five statistics of the left-right difference
% and sum), as montages and volumes. Returns the half atlas and its width.

fprintf('Starting Whole-Brain Multi-Metric Analysis (Leaf Nodes)...\n');

% the left half of the atlas
[~, ~, n_width] = size(AllenCrop);
half_width = floor(n_width / 2);
atlas_left = AllenCrop(:, :, 1:half_width);

% its voxels in a region, and their region numbered 1 to N in place of the Allen
% ids
valid_mask = atlas_left > 0;
all_ids_vec = atlas_left(valid_mask);
[unique_ids, ~, id_indices] = unique(all_ids_vec);
n_unique_regions = length(unique_ids);

fprintf('  Found %d unique leaf regions.\n', n_unique_regions);

% the statistics: the mean, then the median and three quantiles, NaN left out
metric_names = {'Mean', 'Median', 'Q1', 'Q3', 'P99'};
n_metrics = length(metric_names);
nan_quant = @(x, p) quantile(x(~isnan(x)), p);
funcs = {
    @(x) median(x, 'omitnan'), ...
    @(x) nan_quant(x, 0.25), ...
    @(x) nan_quant(x, 0.75), ...
    @(x) nan_quant(x, 0.99) ...
    };

% the background masks, the same for the difference and the sum
bg_mask_ctrl_left = mask_bg_ctrl;
bg_mask_exp_left = mask_bg_exp;

% the difference, then the sum
analysis_types = {'Diff', 'Sum'};
for a_idx = 1:length(analysis_types)

    res_type = analysis_types{a_idx};
    fprintf('--- Whole-Brain Analysis: Processing %s Data ---\n', res_type);

    if strcmp(res_type, 'Diff')
        input_4d_ctrl = lr_diff_ctrl;
        input_4d_exp = lr_diff_exp;
    else
        input_4d_ctrl = lr_sum_ctrl;
        input_4d_exp = lr_sum_exp;
    end

    [leaf_stats_ctrl, leaf_stats_exp] = leaf_region_stats(input_4d_ctrl, input_4d_exp, ...
        bg_mask_ctrl_left, bg_mask_exp_left, valid_mask, id_indices, n_unique_regions, ...
        n_ctrl, n_exp, n_metrics, funcs, res_type);

    plot_tmap_montages(leaf_stats_ctrl, leaf_stats_exp, atlas_left, valid_mask, ...
        id_indices, unique_ids, metric_names, n_metrics, res_type, comp_tag, ...
        comp_out_dir);
end

% free the large arrays
clear bg_mask_ctrl_left bg_mask_exp_left t_score_vol t_score_valid_pixels ...
    input_4d_ctrl input_4d_exp

fprintf('Whole-brain multi-metric analysis (Diff & Sum) saved to: %s\n', comp_out_dir);
end

function [leaf_stats_ctrl, leaf_stats_exp] = leaf_region_stats(input_4d_ctrl, ...
    input_4d_exp, bg_mask_ctrl_left, bg_mask_exp_left, valid_mask, id_indices, ...
    n_unique_regions, n_ctrl, n_exp, n_metrics, funcs, res_type)
% Per mouse and leaf region, the mean and four quantiles of the absolute
% left-right values in the tissue.

% the control group, one mouse at a time
fprintf('  Computing stats for control group (%s)...\n', res_type);
leaf_stats_ctrl = group_leaf_stats(input_4d_ctrl, bg_mask_ctrl_left, valid_mask, ...
    id_indices, n_unique_regions, n_ctrl, n_metrics, funcs);

% the same for the experimental group
fprintf('  Computing stats for experimental group (%s)...\n', res_type);
leaf_stats_exp = group_leaf_stats(input_4d_exp, bg_mask_exp_left, valid_mask, ...
    id_indices, n_unique_regions, n_exp, n_metrics, funcs);
end

function leaf_stats = group_leaf_stats(input_4d, bg_mask_left, valid_mask, id_indices, ...
    n_unique_regions, n_mice, n_metrics, funcs)
% One group's statistics per leaf region and mouse (region x mouse x statistic).

leaf_stats = nan(n_unique_regions, n_mice, n_metrics);
for m = 1:n_mice

    % the mouse's absolute values in the atlas, NaN on its background
    vol_data = abs(input_4d(:, :, :, m));
    vol_bg = bg_mask_left(:, :, :, m);
    data_vec = vol_data(valid_mask);
    bg_vec = vol_bg(valid_mask);
    data_vec(bg_vec) = NaN;

    % the mean of each region
    sums = accumarray(id_indices, data_vec, [n_unique_regions 1], @nansum); %#ok<*NANSUM>
    counts = accumarray(id_indices, ~isnan(data_vec), [n_unique_regions 1], @sum);
    leaf_stats(:, m, 1) = sums ./ counts;

    % the median and quantiles of each region
    for f = 1:length(funcs)
        leaf_stats(:, m, f+1) = accumarray(id_indices, data_vec, ...
            [n_unique_regions 1], funcs{f});
    end
end
end

function plot_tmap_montages(leaf_stats_ctrl, leaf_stats_exp, atlas_left, valid_mask, ...
    id_indices, unique_ids, metric_names, n_metrics, res_type, comp_tag, comp_out_dir)
% One montage of group t-scores per statistic, each saved with its volume.

% every 50th plane from 150 to 150 before the end, five per row
slices_to_show = 150:50:size(atlas_left, 1)-150;
n_cols = 5;
n_rows = ceil(length(slices_to_show)/n_cols);

% blue-white-red colormap
custom_cmap = sep_palette('difference');

fprintf('  Generating T-Maps for %s metrics...\n', res_type);

for i_met = 1:n_metrics
    metric_name = metric_names{i_met};

    % the t-score of each region, and the colour limits
    [t_scores_vec, t_lims] = region_t_scores(leaf_stats_ctrl, i_met, leaf_stats_exp, ...
        metric_name);

    % the region t-scores back into a volume
    t_score_vol = zeros(size(atlas_left), 'single');
    t_score_valid_pixels = t_scores_vec(id_indices);
    t_score_vol(valid_mask) = t_score_valid_pixels;

    % the montage
    fig_h = draw_tmap_montage(res_type, metric_name, comp_tag, slices_to_show, n_rows, ...
        n_cols, t_score_vol, atlas_left, custom_cmap, t_lims);

    % save the figure, and the volume
    set(fig_h, 'InvertHardcopy', 'off');
    save_base = fullfile(comp_out_dir, ...
        ['WholeBrain_TMap_' res_type '_' metric_name '_' comp_tag]);
    saveas(fig_h, [save_base '.fig']);
    exportgraphics(fig_h, [save_base '.png'], 'Resolution', 300, 'BackgroundColor', ...
        'current');
    save([save_base '.mat'], 't_score_vol', 't_scores_vec', 'unique_ids', ...
        'metric_name', 'res_type');
end
end

function [t_scores_vec, t_lims] = region_t_scores(leaf_stats_ctrl, i_met, ...
    leaf_stats_exp, metric_name)
% The group t-score of each region for one statistic (NaN or infinite to 0), and
% symmetric colour limits at the 95th percentile of |t|, at least 0.1.

% the statistic of every mouse
data_c = squeeze(leaf_stats_ctrl(:, :, i_met));
data_e = squeeze(leaf_stats_exp(:, :, i_met));

% t = (mean_exp - mean_ctrl) / pooled SEM, per region
mu_c = nanmean(data_c, 2);
mu_e = nanmean(data_e, 2);
diff_mu = mu_e - mu_c;

sem_c = sem_over_mice(data_c);
sem_e = sem_over_mice(data_e);
pooled_sem = sqrt(sem_c.^2 + sem_e.^2);

t_scores_vec = diff_mu ./ pooled_sem;

% a NaN or infinite t to 0
t_scores_vec(isnan(t_scores_vec) | isinf(t_scores_vec)) = 0;

% colour limits: the 95th percentile of |t|, at least 0.1
t_vals = t_scores_vec(t_scores_vec ~= 0);
if isempty(t_vals)
    max_t = 1;
    fprintf('    [Warning] Metric %s yielded all zero T-scores.\n', metric_name);
else
    max_t = quantile(abs(t_vals), 0.95);
    if max_t < 0.1
        max_t = 0.1;
    end
end
t_lims = [-max_t, max_t];
end

function sem = sem_over_mice(data)
% The SEM of each row over the mice (columns) with a value; NaN with fewer than
% two, where a spread cannot be measured.

n_mice = sum(~isnan(data), 2);
sem = nanstd(data, [], 2) ./ sqrt(n_mice);
sem(n_mice < 2) = NaN;
end

function fig_h = draw_tmap_montage(res_type, metric_name, comp_tag, slices_to_show, ...
    n_rows, n_cols, t_score_vol, atlas_left, custom_cmap, t_lims)
% The montage of one statistic's t-score volume, every 50th plane.

% the montage
fig_h = figure('Visible', 'off', 'Name', ...
    ['WholeBrain_TMap_' res_type '_' metric_name '_' comp_tag], 'Color', 'k', ...
    'Position', [50 50 1200 900]);

for k = 1:length(slices_to_show)
    s_idx = slices_to_show(k);
    subplot(n_rows, n_cols, k);
    im_slice = squeeze(t_score_vol(s_idx, :, :));
    mask_slice = squeeze(atlas_left(s_idx, :, :));
    alpha_data = double(mask_slice > 0);
    imagesc(im_slice);
    axis image;
    axis off;
    set(gca, 'Color', 'k');
    set(findobj(gca, 'Type', 'image'), 'AlphaData', alpha_data);
    colormap(gca, custom_cmap);
    clim(t_lims);
    title(['Slice ' num2str(s_idx)], 'Color', 'w', 'FontSize', 8);
end
c = colorbar;
c.Position = [0.92 0.1 0.02 0.8];
c.Color = 'w';
c.Label.String = sprintf('T-score (%s %s)', res_type, metric_name);
colormap(c, custom_cmap);
clim(t_lims);

sgtitle(['Whole-Brain T-Scores (' res_type '): ' metric_name], 'Color', 'w', ...
    'FontSize', 14);
end

function coarse_region_tstats(AllenCrop, allenDir, brainMask, half_width, ...
    lr_diff_ctrl, lr_diff_exp, lr_sum_ctrl, lr_sum_exp, mask_bg_ctrl, mask_bg_exp, ...
    n_ctrl, n_exp, ctrl_mousenames, exp_mousenames, ctrl_type, exp_type, comp_tag, ...
    comp_out_dir)
% Group t-scores of the left-right difference and sum in a fixed list of
% regions, as bar charts, with optional per-region distributions.

fprintf('Starting Area-based quantification (Vectorized)...\n');

% the statistic of the bar chart
choosen_metric = 'P99';

% regions whose distributions are drawn (none)
rois_to_inspect = {};

% the regions, and their voxels in the left hemisphere
roi_list = coarse_roi_list();
n_rois = length(roi_list);

[valid_pixels, roi_masks, roi_pixel_counts] = coarse_roi_masks(AllenCrop, half_width, ...
    allenDir, brainMask, roi_list, n_rois);

% region x mouse x statistic
metric_names = {'Mean', 'Median', 'Q1', 'Q3', 'P99'};
n_metrics = length(metric_names);
roi_stats_ctrl = nan(n_rois, n_ctrl, n_metrics);
roi_stats_exp = nan(n_rois, n_exp, n_metrics);

% the background masks
bg_mask_ctrl = mask_bg_ctrl;
bg_mask_exp = mask_bg_exp;

% the regions to inspect, and their folder
inspect_indices = find(ismember(roi_list, rois_to_inspect));
dist_out_dir = fullfile(comp_out_dir, 'distribution_checks');
if ~exist(dist_out_dir, 'dir') && ~isempty(inspect_indices)
    mkdir(dist_out_dir);
end

% the difference, then the sum
analysis_types = {'Diff', 'Sum'};

for a_idx = 1:length(analysis_types)
    res_type = analysis_types{a_idx};
    fprintf('--- Processing %s Data for ROIs ---\n', res_type);

    % the absolute values, and the histograms' range
    if strcmp(res_type, 'Diff')
        input_vol_ctrl = abs(lr_diff_ctrl);
        input_vol_exp = abs(lr_diff_exp);
        hist_xlim = [0 2.5];
    else
        input_vol_ctrl = abs(lr_sum_ctrl);
        input_vol_exp = abs(lr_sum_exp);
        hist_xlim = [0 5];
    end

    % one inspection figure per region to inspect
    inspect_figs = inspection_figures(inspect_indices, roi_list, res_type, n_ctrl, ...
        n_exp);

    % the control group, in the top row of the inspection figures
    roi_stats_ctrl = coarse_group_stats(roi_stats_ctrl, input_vol_ctrl, bg_mask_ctrl, ...
        valid_pixels, roi_masks, roi_pixel_counts, n_rois, n_ctrl, n_ctrl, n_exp, ...
        inspect_indices, inspect_figs, hist_xlim, ctrl_mousenames, res_type, 0, ...
        sep_palette('control'), 'Control');

    % the experimental group, in the bottom row
    roi_stats_exp = coarse_group_stats(roi_stats_exp, input_vol_exp, bg_mask_exp, ...
        valid_pixels, roi_masks, roi_pixel_counts, n_rois, n_exp, n_ctrl, n_exp, ...
        inspect_indices, inspect_figs, hist_xlim, exp_mousenames, res_type, ...
        max(n_ctrl, n_exp), sep_palette('experimental'), 'Exp');

    % save the inspection figures
    save_inspection_figures(inspect_figs, inspect_indices, roi_list, res_type, ...
        dist_out_dir);

    % the bar chart
    plot_coarse_bars(roi_stats_ctrl, roi_stats_exp, metric_names, choosen_metric, ...
        roi_list, n_ctrl, n_exp, res_type, exp_type, ctrl_type, comp_tag, comp_out_dir);
end
end

function inspect_figs = inspection_figures(inspect_indices, roi_list, res_type, ...
    n_ctrl, n_exp)
% One hidden figure per region to inspect.

% one inspection figure per region to inspect
inspect_figs = gobjects(length(inspect_indices), 1);
for k = 1:length(inspect_indices)
    r_idx = inspect_indices(k);
    inspect_figs(k) = figure('Visible', 'off', 'Name', ...
        ['Dist_' roi_list{r_idx} '_' res_type], 'Color', 'w', 'Visible', 'off', ...
        'Position', [100 100 300*max(n_ctrl, n_exp) 600]);
end
end

function save_inspection_figures(inspect_figs, inspect_indices, roi_list, res_type, ...
    dist_out_dir)
% Title, save and close the inspection figures, in dist_out_dir.

% save the inspection figures
for k = 1:length(inspect_figs)
    if isvalid(inspect_figs(k))
        r_idx = inspect_indices(k);
        r_name = roi_list{r_idx};
        clean_name = regexprep(r_name, '[^a-zA-Z0-9]', '_');
        sgtitle(inspect_figs(k), ['Distribution: ' r_name ' (' res_type ...
            ') | Solid=Mean, Dash=Med, Purple=P99'], 'Interpreter', 'none');

        saveas(inspect_figs(k), ...
            fullfile(dist_out_dir, ['Dist_' res_type '_' clean_name '.fig']));
        exportgraphics(inspect_figs(k), ...
            fullfile(dist_out_dir, ['Dist_' res_type '_' clean_name '.png']));
        close(inspect_figs(k));
    end
end
end

function roi_list = coarse_roi_list()
% The regions of the coarse analysis.

% each name as the atlas spells it: get_allen_region_mask resolves a name the atlas
% lacks to the first name containing it, which can be another region
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
    'Visceral area', ...
    'Gustatory areas', ...
    'Piriform area', ...
    'Subiculum', ...
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
    'Ventral medial nucleus of the thalamus', ...
    'Zona incerta', ...
    'Posterior complex of the thalamus', ...
    'Lateral posterior nucleus of the thalamus', ...
    'Lateral dorsal nucleus of thalamus', ...
    'Ventral anterior-lateral complex of the thalamus', ...
    'Mediodorsal nucleus of thalamus', ...
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
    'Hypothalamus' ...
    };
end

function [valid_pixels, roi_masks, roi_pixel_counts] = coarse_roi_masks(AllenCrop, ...
    half_width, allenDir, brainMask, roi_list, n_rois)
% The atlas voxels of the left hemisphere, and which of them fall in each region.

% the left hemisphere's atlas voxels, as a vector
fprintf('Mapping Atlas Volume...\n');
atlas_left = AllenCrop(:, :, 1:half_width);
valid_pixels = atlas_left > 0;
pixel_ids = atlas_left(valid_pixels);

% for each region, which of those voxels carry one of its ids (regions may overlap)
fprintf('Building masks for %d regions (Handling overlaps)...\n', n_rois);

roi_masks = false(length(pixel_ids), n_rois);
roi_pixel_counts = zeros(n_rois, 1);

for r = 1:n_rois
    region_name = roi_list{r};
    try
        % the region and its descendants, in the left hemisphere
        mask_temp = get_allen_region_mask(allenDir, AllenCrop, {region_name}, ...
            brainMask, '');
        if size(mask_temp, 3) >= half_width
            mask_temp = mask_temp(:, :, 1:half_width);
        end

        % the atlas ids inside it, and the voxels that carry them
        ids_in_region = unique(atlas_left(mask_temp));
        ids_in_region(ids_in_region == 0) = [];
        roi_masks(:, r) = ismember(pixel_ids, ids_in_region);
        roi_pixel_counts(r) = sum(roi_masks(:, r));
    catch
        warning('Could not map region: %s', region_name);
    end
end

% the regions with no voxel
empty_rois = find(roi_pixel_counts == 0);
if ~isempty(empty_rois)
    fprintf('Warning: The following regions have 0 pixels:\n');
    disp(roi_list(empty_rois)');
end
end

function roi_stats = coarse_group_stats(roi_stats, input_vol, bg_mask, valid_pixels, ...
    roi_masks, roi_pixel_counts, n_rois, n_mice, n_ctrl, n_exp, inspect_indices, ...
    inspect_figs, hist_xlim, mousenames, res_type, row_offset, face_color, group_label)
% Per mouse of one group and region, the mean and four quantiles of the absolute
% values, with the distributions drawn into the inspection figures: the control
% group in the top row (row_offset 0), the experimental one below.

for m = 1:n_mice
    tic

    % the mouse's values in the atlas, NaN on its background
    vol_data = input_vol(:, :, :, m);
    vol_bg = bg_mask(:, :, :, m);

    data_vec = vol_data(valid_pixels);
    bg_vec = vol_bg(valid_pixels);
    data_vec(bg_vec) = NaN;

    for r = 1:n_rois
        if roi_pixel_counts(r) == 0
            continue;
        end

        % the region's values, their mean
        roi_vals = data_vec(roi_masks(:, r));
        roi_stats(r, m, 1) = mean(roi_vals, 'omitnan');

        % the median and quantiles, without the NaN
        roi_vals_clean = roi_vals(~isnan(roi_vals));
        if ~isempty(roi_vals_clean)
            roi_stats(r, m, 2) = median(roi_vals_clean);
            roi_stats(r, m, 3) = quantile(roi_vals_clean, 0.25);
            roi_stats(r, m, 4) = quantile(roi_vals_clean, 0.75);
            roi_stats(r, m, 5) = quantile(roi_vals_clean, 0.99);

            % the histogram, in the group's row of the region's inspection figure
            if ismember(r, inspect_indices)
                k_fig = find(inspect_indices == r);
                set(0, 'CurrentFigure', inspect_figs(k_fig));
                subplot(2, max(n_ctrl, n_exp), m + row_offset);
                histogram(roi_vals_clean, 100, 'EdgeColor', 'none', 'FaceColor', ...
                    face_color);
                hold on;
                xline(roi_stats(r, m, 1), 'k-', 'LineWidth', 1.5);
                xline(roi_stats(r, m, 2), 'k--', 'LineWidth', 1.5);
                xline(roi_stats(r, m, 5), '-', 'LineWidth', 1.5, 'Color', [1, 0, 1]);
                xlim(hist_xlim);
                txt_str = [sprintf('Mean: %.2f\nP99: %.2f', roi_stats(r, m, 1), ...
                    roi_stats(r, m, 5)), ' (N=', num2str(numel(roi_vals_clean)), ')'];
                text(0.95, 0.9, txt_str, 'Units', 'normalized', 'HorizontalAlignment', ...
                    'right', 'FontSize', 8, 'BackgroundColor', 'w', 'EdgeColor', 'k');
                if exist('mousenames', 'var')
                    t_str = mousenames{m};
                else
                    t_str = sprintf('M%d', m);
                end
                title(t_str, 'Interpreter', 'none', 'FontSize', 8);
                if m==1
                    ylabel([group_label ' (' res_type ')'], 'FontWeight', 'bold');
                end
                grid on;
            end
        end
    end
    toc
end
end

function plot_coarse_bars(roi_stats_ctrl, roi_stats_exp, metric_names, choosen_metric, ...
    roi_list, n_ctrl, n_exp, res_type, exp_type, ctrl_type, comp_tag, comp_out_dir)
% The bar chart of region t-scores for the chosen statistic, with Bonferroni
% thresholds.

% the t-score of each region, sorted
[sorted_t, sorted_rois] = coarse_t_scores(metric_names, choosen_metric, ...
    roi_stats_ctrl, roi_stats_exp, roi_list);

% the bars: red above zero, blue below
fig_bars = figure('Visible', 'off', 'Name', ['Region_Analysis_BarChart_' res_type '_' ...
    comp_tag], 'Color', 'w', 'Units', 'Normalized', 'Position', [0 0 0.9 0.9]);

b = barh(sorted_t);
b.FaceColor = 'flat';
for k = 1:length(sorted_t)
    if sorted_t(k) > 0
        b.CData(k, :) = sep_palette('experimental');
    else
        b.CData(k, :) = sep_palette('control');
    end
end

yticks(1:length(sorted_rois));
yticklabels(sorted_rois);
xlabel(['t-score (' exp_type ' - ' ctrl_type ')']);
title(['Regional LR - ' res_type ' differences - ', strrep(comp_tag, '_', ' '), ' - ', ...
    choosen_metric]);
grid on;
set(gca, 'FontSize', 10);
ylim([0 length(sorted_rois)+1]);

% the regions expected to change, labelled in bold magenta
highlight_coarse_regions(exp_type);

% the Bonferroni thresholds over the regions tested, p < 0.05 two-sided
n_regions_tested = length(sorted_t);
df = n_ctrl + n_exp - 2;
t_crit_bonf = tinv(1 - 0.05/(2*n_regions_tested), df);
hold on;
xline(t_crit_bonf, 'k:', 'LineWidth', 2);
xline(-t_crit_bonf, 'k:', 'LineWidth', 2);
hold off;

% save it
set(fig_bars, 'InvertHardcopy', 'off');
saveas(fig_bars, ...
    fullfile(comp_out_dir, ['Region_Stats_Bar_' res_type '_' comp_tag '.fig']));
exportgraphics(fig_bars, ...
    fullfile(comp_out_dir, ['Region_Stats_Bar_' res_type '_' comp_tag '.png']), ...
    'Resolution', 300);
end

function [sorted_t, sorted_rois] = coarse_t_scores(metric_names, choosen_metric, ...
    roi_stats_ctrl, roi_stats_exp, roi_list)
% The t-score of each region for the chosen statistic, the regions with one,
% sorted.

% the chosen statistic (the mean if it is not in the list)
target_metric_idx = find(strcmp(metric_names, choosen_metric));
if isempty(target_metric_idx)
    target_metric_idx = 1;
end

roi_means_ctrl = roi_stats_ctrl(:, :, target_metric_idx);
roi_means_exp = roi_stats_exp(:, :, target_metric_idx);

% t = (mean_exp - mean_ctrl) / pooled SEM, per region
mean_ctrl_roi = nanmean(roi_means_ctrl, 2);
mean_exp_roi = nanmean(roi_means_exp, 2);
diff_means = mean_exp_roi - mean_ctrl_roi;

sem_c = sem_over_mice(roi_means_ctrl);
sem_e = sem_over_mice(roi_means_exp);
pooled_sem = sqrt(sem_c.^2 + sem_e.^2);
t_score_roi = diff_means ./ pooled_sem;

% the regions with a t, sorted by it
valid_rows = ~isnan(t_score_roi);
t_score_roi = t_score_roi(valid_rows);
current_rois = roi_list(valid_rows);

[sorted_t, sort_idx] = sort(t_score_roi, 'ascend');
sorted_rois = current_rois(sort_idx);
end

function highlight_coarse_regions(exp_type)
% The tick labels of the regions expected to change, in bold magenta.

% the regions expected to change, labelled in bold magenta
switch exp_type
    case {'rws', 'behavior'}
        highlighted_areas = {'Primary somatosensory area, barrel field', ...
            'Ventral posteromedial nucleus of the thalamus', ...
            'Posterior complex of the thalamus', 'Supplemental somatosensory area', ...
            'Zona incerta'};
    otherwise
        highlighted_areas = {};
end

ax = gca;
ytl = ax.YTickLabel;
colored_labels = repmat({''}, size(ytl));
for i = 1:length(ytl)
    if ismember(ytl{i}, highlighted_areas)
        colored_labels{i} = ['\color{magenta} \bf ' strrep(ytl{i}, '_', ' ')];
    else
        colored_labels{i} = ['\color{black} ' strrep(ytl{i}, '_', ' ')];
    end
end
ax = gca;
ax.YTickLabel = colored_labels;
end

function write_annotated_tmap_video(AllenCrop, allenDir, comp_out_dir, comp_tag)
% Video of the P99 t-score volumes the montage step saved, every plane labelled
% with its region acronyms.

fprintf('Generating Annotated T-Score Video (Direct Mapping)...\n');

% the statistic shown
metric_to_plot = 'P99';

% Allen id to acronym
id2acronym = load_acronym_map(allenDir);

% the t-score volumes of the difference and the sum, saved by the montage step
file_diff = fullfile(comp_out_dir, ...
    ['WholeBrain_TMap_Diff_' metric_to_plot '_' comp_tag '.mat']);
file_sum = fullfile(comp_out_dir, ...
    ['WholeBrain_TMap_Sum_' metric_to_plot '_' comp_tag '.mat']);

if exist(file_diff, 'file') && exist(file_sum, 'file')
    data_diff = load(file_diff, 't_score_vol');
    vol_diff = data_diff.t_score_vol;

    data_sum = load(file_sum, 't_score_vol');
    vol_sum = data_sum.t_score_vol;

    % the atlas over the volumes' width
    atlas_hi_res = double(AllenCrop);
    atlas_left_hires = atlas_hi_res(:, :, 1:size(vol_diff, 3));

    % 10 frames per second, quality 95
    video_filename = ['Annotated_TMap_' metric_to_plot '_' comp_tag '.mp4'];
    full_video_path = fullfile(comp_out_dir, video_filename);

    vidObj = VideoWriter(full_video_path, 'MPEG-4');
    vidObj.FrameRate = 10;
    vidObj.Quality = 95;
    open(vidObj);
    n_slices = size(vol_diff, 1);

    % blue-white-red colormap
    custom_cmap = sep_palette('difference');

    % one figure, cleared after each frame
    fh = figure('visible', 'off', 'units', 'normalized', 'outerposition', [0 0 1 1], ...
        'Color', 'k');
    set(fh, 'InvertHardcopy', 'off');

    fprintf('Writing annotated video: %s\n', video_filename);

    % every plane with atlas voxels, labelled
    write_annotated_frames(n_slices, atlas_left_hires, id2acronym, vol_diff, vol_sum, ...
        custom_cmap, metric_to_plot, fh, vidObj);

    close(vidObj);
    close(fh);
    fprintf('Annotated video saved: %s\n', full_video_path);

else
    warning('T-Map .mat files not found. Run the Whole-Brain Analysis section first.');
end
clear atlas_hi_res atlas_left_hires
end

function write_annotated_frames(n_slices, atlas_left_hires, id2acronym, vol_diff, ...
    vol_sum, custom_cmap, metric_to_plot, fh, vidObj)
% One labelled frame per plane with atlas voxels, written to vidObj.

for j = 1:n_slices

    % skip the planes with no atlas voxel
    mask_slice_hires = squeeze(atlas_left_hires(j, :, :));
    if sum(mask_slice_hires(:) > 0) == 0 %#ok<LOGSUM>
        continue;
    end

    draw_annotated_frame(mask_slice_hires, id2acronym, vol_diff, vol_sum, j, ...
        custom_cmap, metric_to_plot);

    frame = getframe(fh);
    writeVideo(vidObj, frame);

    if mod(j, 50) == 0
        fprintf('  Frame %d written...\n', j);
    end
    clf(fh);
end
end

function id2acronym = load_acronym_map(allenDir)
% Allen region id to acronym, from the parcellation table when it is there
% (an empty map otherwise).

map_file = fullfile(allenDir, 'parcellation_to_parcellation_term_membership.csv');
id2acronym = containers.Map('KeyType', 'double', 'ValueType', 'char');

if exist(map_file, 'file')

    % every column read as text
    opts = detectImportOptions(map_file);
    opts.VariableTypes = repmat({'string'}, 1, length(opts.VariableTypes));
    T_map = readtable(map_file, opts);

    % the id and acronym columns must be there
    if ismember('parcellation_index', T_map.Properties.VariableNames) && ...
            ismember('parcellation_term_acronym', T_map.Properties.VariableNames)
        raw_pixels = str2double(T_map.parcellation_index);
        raw_acros = T_map.parcellation_term_acronym;

        % keep the structure level of the hierarchy when the table has one
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
        leaf_pixels = raw_pixels(is_leaf);
        leaf_acros = raw_acros(is_leaf);

        % without the ids that are not numbers
        valid_idx = ~isnan(leaf_pixels);
        leaf_pixels = leaf_pixels(valid_idx);
        leaf_acros = cellstr(leaf_acros(valid_idx));

        % one acronym per id, the last one listed
        [unique_pixels, idx] = unique(leaf_pixels, 'last');
        unique_acros = leaf_acros(idx);

        if ~isempty(unique_pixels)
            id2acronym = containers.Map(unique_pixels, unique_acros);
            fprintf('  Mapping created successfully: %d leaf regions mapped.\n', ...
                length(unique_pixels));
        else
            warning('Mapping failed. No valid pixel-acronym pairs found.');
        end
    else
        warning('Required columns (parcellation_index, parcellation_term_acronym) missing.');
    end
else
    warning('Mapping CSV not found. Labels will be skipped.');
end
end

function draw_annotated_frame(mask_slice_hires, id2acronym, vol_diff, vol_sum, j, ...
    custom_cmap, metric_to_plot)
% One frame of the annotated video, drawn into the current figure.

alpha_data = double(mask_slice_hires > 0);

% a label at the centre of each region of at least 80 voxels
lbl_data = region_labels(mask_slice_hires, id2acronym);

% left: the difference
subplot(1, 2, 1);
imagesc(squeeze(vol_diff(j, :, :)));
clim([-4 4]);
colormap(gca, custom_cmap);
c = colorbar;
c.Color = 'w';
c.Label.String = sprintf('T-score (%s %s)', 'Diff', metric_to_plot);
set(findobj(gca, 'Type', 'image'), 'AlphaData', alpha_data);
axis image;
axis off;
set(gca, 'Color', 'k');
title(['Diff T-Score (' metric_to_plot ')'], 'Color', 'w', 'FontSize', 12);

if ~isempty(lbl_data)
    text([lbl_data.x], [lbl_data.y], {lbl_data.str}, ...
        'Color', 'k', 'FontSize', 6, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'Interpreter', 'none', ...
        'BackgroundColor', 'w', 'Margin', 0.5, 'EdgeColor', 'none');
end

% right: the sum
subplot(1, 2, 2);
imagesc(squeeze(vol_sum(j, :, :)));
clim([-4 4]);
colormap(gca, custom_cmap);
c = colorbar;
c.Color = 'w';
c.Label.String = sprintf('T-score (%s %s)', 'Sum', metric_to_plot);
set(findobj(gca, 'Type', 'image'), 'AlphaData', alpha_data);
axis image;
axis off;
set(gca, 'Color', 'k');
title(['Sum T-Score (' metric_to_plot ')'], 'Color', 'w', 'FontSize', 12);

if ~isempty(lbl_data)
    text([lbl_data.x], [lbl_data.y], {lbl_data.str}, ...
        'Color', 'k', 'FontSize', 6, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'Interpreter', 'none', ...
        'BackgroundColor', 'w', 'Margin', 0.5, 'EdgeColor', 'none');
end

sgtitle(['Slice # ' num2str(j)], 'Color', 'w', 'FontSize', 14);
end

function lbl_data = region_labels(mask_slice_hires, id2acronym)
% The acronym and centre of each region of at least 80 voxels on the plane.

% a label at the centre of each region of at least 80 voxels with a known acronym
regions_in_slice = unique(mask_slice_hires(mask_slice_hires > 0));

lbl_data = struct('x', {}, 'y', {}, 'str', {});
idx_lbl = 1;

if ~isempty(id2acronym)
    for k = 1:length(regions_in_slice)
        r_id = regions_in_slice(k);

        bin_mask = (mask_slice_hires == r_id);
        if sum(bin_mask(:)) < 80
            continue;
        end

        if isKey(id2acronym, r_id)
            [py, px] = find(bin_mask);
            lbl_data(idx_lbl).x = mean(px);
            lbl_data(idx_lbl).y = mean(py);
            lbl_data(idx_lbl).str = id2acronym(r_id);
            idx_lbl = idx_lbl + 1;
        end
    end
end
end

% ===== Local functions: slab figures, rolling videos, regional surprise =====

function plot_group_slab(t_lr_diff_groupdiff, t_lr_sum_groupdiff, surp_diff, surp_sum, ...
    brainMask_group_diff, half_atlas, exp_type, ctrl_type, comp_out_dir, comp_tag)
% The group t-maps around plane 565, median over the slab, shown where the
% surprise reaches p < 0.01.

% the slab: plane 565 and 10 planes on each side
target_slice = 565;
slab_range = 10;
save_fig = true;

z_start = max(1, target_slice - slab_range);
z_end = min(size(t_lr_diff_groupdiff, 1), target_slice + slab_range);
z_indices = z_start:z_end;

fprintf('Averaging signal across slices %d to %d (Target: %d)...\n', z_start, z_end, ...
    target_slice);

% the medians over the slab, and the opacities
[slab_diff, slab_sum, p_thresh, alpha_mask_diff, alpha_mask_sum] = group_slab_medians( ...
    brainMask_group_diff, z_indices, t_lr_diff_groupdiff, t_lr_sum_groupdiff, ...
    surp_diff, surp_sum);

% the atlas boundaries of the central plane
atlas_slice = squeeze(half_atlas(target_slice, :, 1:size(slab_diff, 2)));
[gy, gx] = gradient(single(atlas_slice));
boundaries = (abs(gx) + abs(gy)) > 0 & (atlas_slice > 0);
[b_row, b_col] = find(boundaries);

% the figure
fig_slab = draw_group_slab(target_slice, slab_diff, alpha_mask_diff, b_col, b_row, ...
    slab_range, slab_sum, alpha_mask_sum, exp_type, ctrl_type, p_thresh);

% save it
if save_fig
    if ~exist(comp_out_dir, 'dir')
        mkdir(comp_out_dir);
    end
    filename = sprintf('Slab_Avg_%d_%s_surpmask', target_slice, comp_tag);
    saveas(fig_slab, fullfile(comp_out_dir, [filename '.fig']));
    exportgraphics(fig_slab, fullfile(comp_out_dir, [filename '.png']), ...
        'Resolution', 300, 'BackgroundColor', 'current');
    fprintf('Saved slab average figure to: %s\n', fullfile(comp_out_dir, filename));
end
end

function [slab_diff, slab_sum, p_thresh, alpha_mask_diff, alpha_mask_sum] = ...
    group_slab_medians(brainMask_group_diff, z_indices, t_lr_diff_groupdiff, ...
    t_lr_sum_groupdiff, surp_diff, surp_sum)
% The medians over the slab of the t maps and their surprise, and each panel's
% opacity: the surprise over -log10(0.01), clipped, where the slab has voxels.

% the medians over the slab, inside the voxels both groups have: the t of the
% difference and of the sum, and their surprise
mask_slab_3d = brainMask_group_diff(z_indices, :, :);

raw_diff = t_lr_diff_groupdiff(z_indices, :, :);
raw_diff(~mask_slab_3d) = NaN;
slab_diff = squeeze(nanmedian(raw_diff, 1));

raw_sum = t_lr_sum_groupdiff(z_indices, :, :);
raw_sum(~mask_slab_3d) = NaN;
slab_sum = squeeze(nanmedian(raw_sum, 1));

raw_surp_d = surp_diff(z_indices, :, :);
raw_surp_d(~mask_slab_3d) = NaN;
slab_surp_diff = squeeze(nanmedian(raw_surp_d, 1));

raw_surp_s = surp_sum(z_indices, :, :);
raw_surp_s(~mask_slab_3d) = NaN;
slab_surp_sum = squeeze(nanmedian(raw_surp_s, 1));

% opacity: the surprise over -log10(0.01), clipped to 0-1
p_thresh = 0.01;
surp_thresh = -log10(p_thresh);
calc_alpha = @(vol) min(1, max(0, vol ./ surp_thresh));
alpha_diff = calc_alpha(slab_surp_diff);
alpha_diff(isnan(alpha_diff)) = 0;
alpha_sum = calc_alpha(slab_surp_sum);
alpha_sum(isnan(alpha_sum)) = 0;

% and only where some plane of the slab has the voxel
slab_mask_2d = squeeze(max(mask_slab_3d, [], 1));
alpha_mask_diff = alpha_diff .* double(slab_mask_2d);
alpha_mask_sum = alpha_sum .* double(slab_mask_2d);
end

function fig_slab = draw_group_slab(target_slice, slab_diff, alpha_mask_diff, b_col, ...
    b_row, slab_range, slab_sum, alpha_mask_sum, exp_type, ctrl_type, p_thresh)
% The slab figure: the t of the difference and of the sum, at their opacities.

fig_slab = figure('Visible', 'off', 'Name', sprintf('Slab_Avg_%d', target_slice), ...
    'Color', 'k', 'Position', [100 100 1200 600]);
set(fig_slab, 'InvertHardcopy', 'off');

% left: the t of the difference, on t limits of +/- 6
subplot(1, 2, 1);
h1 = imagesc(slab_diff);
set(h1, 'AlphaData', alpha_mask_diff);
axis image;
axis off;
set(gca, 'Color', 'k');
hold on;
plot(b_col, b_row, '.', 'Color', [0.7 0.7 0.7], 'MarkerSize', 0.25);
clim([-6 6]);
colormap(gca, sep_palette('difference'));
cb1 = colorbar;
cb1.Color = 'w';
cb1.Label.String = 'T-Score (Diff)';
cb1.Label.Color = 'w';
title(['LR Diff (T-Score) - Slab Avg ' num2str(target_slice) '\pm' ...
    num2str(slab_range)], 'Color', 'w');

% right: the t of the sum
subplot(1, 2, 2);
h2 = imagesc(slab_sum);
set(h2, 'AlphaData', alpha_mask_sum);
axis image;
axis off;
set(gca, 'Color', 'k');
hold on;
plot(b_col, b_row, '.', 'Color', [0.7 0.7 0.7], 'MarkerSize', 0.25);
clim([-6 6]);
colormap(gca, sep_palette('difference'));
cb2 = colorbar;
cb2.Color = 'w';
cb2.Label.String = 'T-Score (Sum)';
cb2.Label.Color = 'w';
title(['LR Sum (T-Score) - Slab Avg ' num2str(target_slice) '\pm' ...
    num2str(slab_range)], 'Color', 'w');

sgtitle(['Slab average - ', exp_type ' vs ' ctrl_type ' - surprise masked (p<', ...
    num2str(p_thresh), ')'], 'Color', 'w', 'FontSize', 14);
end

function plot_individual_slabs(lr_diff_ctrl, lr_diff_exp, lr_sum_ctrl, lr_sum_exp, ...
    mask_bg_ctrl, mask_bg_exp, ctrl_mousenames, exp_mousenames, ctrl_type, exp_type, ...
    half_atlas, comp_out_dir)
% Every mouse's left-right difference and sum around plane 565, median over the
% slab, one figure per group.

fprintf('Generating Individual Mice Slab Average Plots...\n');

% the slab: plane 565 and 10 planes on each side; colour limits of the
% difference and the sum; the atlas boundaries off
target_slice = 565;
slab_range = 10;
clim_diff_indiv = [0 1.5];
clim_sum_indiv = [0 10];
bool_overlay_atlas = false;

z_start = max(1, target_slice - slab_range);
z_end = min(size(lr_diff_ctrl, 1), target_slice + slab_range);
z_indices = z_start:z_end;

% the atlas boundaries of the central plane
atlas_slice = squeeze(half_atlas(target_slice, :, 1:size(lr_diff_ctrl, 3)));
[gy, gx] = gradient(single(atlas_slice));
boundaries = (abs(gx) + abs(gy)) > 0 & (atlas_slice > 0);
[b_row, b_col] = find(boundaries);

% the control group, then the experimental one
group_names = {ctrl_type, exp_type};
group_data_diff = {abs(lr_diff_ctrl), abs(lr_diff_exp)};
group_data_sum = {lr_sum_ctrl, lr_sum_exp};
group_masks = {mask_bg_ctrl, mask_bg_exp};
group_mice_list = {ctrl_mousenames, exp_mousenames};

for g_idx = 1:2
    curr_grp = group_names{g_idx};
    curr_diff = group_data_diff{g_idx};
    curr_sum = group_data_sum{g_idx};
    curr_mask = group_masks{g_idx};
    curr_mice = group_mice_list{g_idx};

    n_mice = size(curr_diff, 4);

    fig_indiv = figure('Visible', 'off', 'Name', ['Individual_Slab_Avg_' curr_grp], ...
        'Color', 'k', 'Units', 'normalized', 'Position', [-0.05 -0.05 0.9 0.9]);
    set(fig_indiv, 'InvertHardcopy', 'off');

    % every mouse: the difference above, the sum below
    draw_mouse_slab_medians(n_mice, curr_mice, curr_diff, z_indices, curr_sum, ...
        curr_mask, atlas_slice, clim_diff_indiv, bool_overlay_atlas, b_col, b_row, ...
        clim_sum_indiv);
    sgtitle(['Individual Slab Median (' curr_grp ') - Slice ' num2str(target_slice) ...
        '\pm' num2str(slab_range)], 'Color', 'w', 'FontSize', 16);

    % save it
    if exist('comp_out_dir', 'var')
        filename = sprintf('Indiv_Slab_Avg_%d_%s', target_slice, curr_grp);
        saveas(fig_indiv, fullfile(comp_out_dir, [filename '.fig']));
        exportgraphics(fig_indiv, fullfile(comp_out_dir, [filename '.png']), ...
            'Resolution', 300, 'BackgroundColor', 'current');
        fprintf('Saved individual slab figure for %s.\n', curr_grp);
    end

end

clear group_names group_data_diff group_data_sum group_masks group_mice_list curr_grp ...
    curr_diff curr_sum curr_mask curr_mice
end

function draw_mouse_slab_medians(n_mice, curr_mice, curr_diff, z_indices, curr_sum, ...
    curr_mask, atlas_slice, clim_diff_indiv, bool_overlay_atlas, b_col, b_row, ...
    clim_sum_indiv)
% Every mouse's slab median: the difference in the top row, the sum below.

for k = 1:n_mice
    mouse_name = strrep(curr_mice{k}, '_', ' ');

    % the mouse's slab, NaN on its background
    raw_slab_d = curr_diff(z_indices, :, :, k);
    raw_slab_s = curr_sum(z_indices, :, :, k);
    mask_slab = curr_mask(z_indices, :, :, k);
    mask_slab_logical = logical(mask_slab);
    raw_slab_d(mask_slab_logical) = NaN;
    raw_slab_s(mask_slab_logical) = NaN;

    % the median over the slab
    slab_diff_m = squeeze(nanmedian(raw_slab_d, 1)); %#ok<*NANMEDIAN>
    slab_sum_m = squeeze(nanmedian(raw_slab_s, 1));

    % shown inside the atlas where any plane of the slab is tissue
    slab_bg_m = squeeze(min(mask_slab, [], 1));
    valid_pixels = (atlas_slice > 0) & (~slab_bg_m);
    alpha_data = double(valid_pixels);

    % top row: the difference
    subplot(2, n_mice, k);
    imagesc(slab_diff_m);
    set(findobj(gca, 'Type', 'image'), 'AlphaData', alpha_data);
    clim(clim_diff_indiv);
    colormap(gca, sep_palette('intensity'));
    axis image;
    axis off;
    set(gca, 'Color', 'k');
    hold on;
    if bool_overlay_atlas
        plot(b_col, b_row, '.', 'Color', [0.5 0.5 0.5], 'MarkerSize', 0.1);
    end

    title(mouse_name, 'Color', 'w', 'FontSize', 10, 'Interpreter', 'none');
    cb = colorbar;
    cb.Label.String = '|L - R|';
    cb.Color = 'w';
    cb.Label.Color = 'w';

    % bottom row: the sum
    subplot(2, n_mice, k + n_mice);
    imagesc(slab_sum_m);
    set(findobj(gca, 'Type', 'image'), 'AlphaData', alpha_data);
    clim(clim_sum_indiv);
    colormap(gca, sep_palette('intensity'));
    axis image;
    axis off;
    set(gca, 'Color', 'k');
    hold on;
    if bool_overlay_atlas
        plot(b_col, b_row, '.', 'Color', [0.5 0.5 0.5], 'MarkerSize', 0.1);
    end
    cb = colorbar;
    cb.Label.String = 'L + R';
    cb.Color = 'w';
    cb.Label.Color = 'w';

end
end

function write_rolling_tscore_video(t_lr_diff_groupdiff, t_lr_sum_groupdiff, ...
    half_atlas, brainMask_group_diff, surp_diff, surp_sum, comp_out_dir, channel, ...
    comp_tag, exp_type, ctrl_type)
% Video of the group t-maps, median over a rolling slab, masked by surprise.

% median over +/- 10 planes, opaque from p < 0.01, t limits of +/- 6
slab_range = 10;
surp_thresh = -log10(0.01);
t_lim = [-6 6];

vid_name = [['t_lr_diff_sum_' channel '_groupdiff_'] comp_tag ...
    '_surpmask_rolling_slab.mp4'];

fprintf('Generating Rolling Slab Surpmask Video (Window: +/- %d)...\n', slab_range);

% the t of the difference and the sum, each panel's opacity from its own surprise
write_lr_video_surpmask_rolling( ...
    t_lr_diff_groupdiff, ...
    t_lr_sum_groupdiff, ...
    half_atlas, ...
    brainMask_group_diff, ...
    comp_out_dir, ...
    vid_name, ...
    t_lim, ...
    [exp_type ' - ' ctrl_type ' (p<0.01)'], ...
    'T-Score Diff (Slab)', ...
    'T-Score Sum (Slab)', ...
    surp_diff, ...
    surp_sum, ...
    surp_thresh, ...
    slab_range ...
    );
end

function write_individual_rolling_videos(lr_diff_ctrl, lr_sum_ctrl, mask_bg_ctrl, ...
    lr_diff_exp, lr_sum_exp, mask_bg_exp, atlas_left, comp_out_dir, ctrl_type, ...
    exp_type, ctrl_mousenames, exp_mousenames, slab_range, clim_diff, clim_sum)
% Videos of every mouse's left-right difference and sum, median over a rolling
% slab, one per group.

fprintf('Generating Individual LR Rolling Slab Videos...\n');

% the control group: the absolute difference, the sum, the background masks and
% the left half of the atlas
write_lr_indiv_rolling_video( ...
    abs(lr_diff_ctrl), ...
    lr_sum_ctrl, ...
    mask_bg_ctrl, ...
    atlas_left, ...
    comp_out_dir, ...
    ['Individual_Rolling_Slab_' ctrl_type '.mp4'], ...
    ctrl_type, ...
    ctrl_mousenames, ...
    slab_range, ...
    clim_diff, clim_sum ...
    );

% the experimental group
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

function regional_surprise_bars(surp_diff, surp_sum, brainMask_group_diff, AllenCrop, ...
    allenDir, brainMask, exp_type, comp_tag, comp_out_dir)
% Bar charts of the surprise summed over each region of a fixed list, after a
% rolling median over planes, for the difference and the sum.

fprintf('Starting Regional Surprise Analysis (Rolling Median - Diff & Sum)...\n');

% median over +/- 10 planes; a voxel counts in its region's sum from p < 0.01
slab_range = 10;
p_thresh_agg = 0.01;
surp_thresh_val = -log10(p_thresh_agg);

[roi_list_surp, n_rois_surp, valid_pixels, roi_masks_surp, roi_pixel_counts_surp] = ...
    surprise_roi_masks(AllenCrop, allenDir, brainMask);

plot_regional_surprise(surp_diff, surp_sum, brainMask_group_diff, roi_list_surp, ...
    n_rois_surp, valid_pixels, roi_masks_surp, roi_pixel_counts_surp, slab_range, ...
    surp_thresh_val, exp_type, comp_tag, comp_out_dir);
end

function [roi_list_surp, n_rois_surp, valid_pixels, roi_masks_surp, ...
    roi_pixel_counts_surp] = surprise_roi_masks(AllenCrop, allenDir, brainMask)
% The regions of the surprise bars, the atlas voxels of the left hemisphere,
% and which of them fall in each region.

% the regions, as in the coarse analysis
roi_list_surp = coarse_roi_list();
n_rois_surp = length(roi_list_surp);

% the left hemisphere's atlas voxels, as a vector
fprintf('  Mapping Atlas Volume...\n');
[~, ~, n_width] = size(AllenCrop);
half_width = floor(n_width / 2);
atlas_left = AllenCrop(:, :, 1:half_width);
valid_pixels = atlas_left > 0;
pixel_ids = atlas_left(valid_pixels);

% for each region, which of those voxels carry one of its ids
fprintf('  Building masks for %d regions...\n', n_rois_surp);
roi_masks_surp = false(length(pixel_ids), n_rois_surp);
roi_pixel_counts_surp = zeros(n_rois_surp, 1);

for r = 1:n_rois_surp
    region_name = roi_list_surp{r};
    try
        mask_temp = get_allen_region_mask(allenDir, AllenCrop, {region_name}, ...
            brainMask, '');
        if size(mask_temp, 3) >= half_width
            mask_temp = mask_temp(:, :, 1:half_width);
        end
        ids_in_region = unique(atlas_left(mask_temp));
        ids_in_region(ids_in_region == 0) = [];
        roi_masks_surp(:, r) = ismember(pixel_ids, ids_in_region);
        roi_pixel_counts_surp(r) = sum(roi_masks_surp(:, r));
    catch
        warning('Could not map region: %s', region_name);
    end
end
end

function plot_regional_surprise(surp_diff, surp_sum, brainMask_group_diff, ...
    roi_list_surp, n_rois_surp, valid_pixels, roi_masks_surp, roi_pixel_counts_surp, ...
    slab_range, surp_thresh_val, exp_type, comp_tag, comp_out_dir)
% The surprise bar charts, difference and sum side by side.

fig_surp = figure('Visible', 'off', 'Name', ...
    ['Region_Surprise_BarChart_DiffSum_' comp_tag], 'Color', 'w', ...
    'Units', 'Normalized', 'Position', [0 0 0.9 0.9]);

% the difference, then the sum
modes = {'Diff', 'Sum'};

for m_idx = 1:2
    mode_name = modes{m_idx};

    if strcmp(mode_name, 'Diff')
        if ~exist('surp_diff', 'var')
            error('surp_diff not found');
        end
        raw_surp_vol = surp_diff;
    else
        if ~exist('surp_sum', 'var')
            error('surp_sum not found');
        end
        raw_surp_vol = surp_sum;
    end

    vol_surp = rolling_surprise_median(raw_surp_vol, slab_range, brainMask_group_diff, ...
        mode_name);

    % each region's sum of the surprise above the threshold
    roi_surp_agg = summed_surprise(vol_surp, valid_pixels, n_rois_surp, ...
        roi_pixel_counts_surp, roi_masks_surp, surp_thresh_val);

    % sorted, and only the regions with a sum above zero
    sort_metric = roi_surp_agg;
    [sorted_surp, sort_idx] = sort(sort_metric, 'ascend');
    sorted_rois_surp = roi_list_surp(sort_idx);

    valid_k = sorted_surp > 0 & ~isnan(sorted_surp);
    sorted_surp = sorted_surp(valid_k);
    sorted_rois_surp = sorted_rois_surp(valid_k);

    subplot(1, 2, m_idx);

    % the bars
    draw_surprise_bars(sorted_surp, sorted_rois_surp, surp_thresh_val, mode_name);

    % the regions expected to change, labelled in bold magenta
    highlight_surprise_regions(exp_type);
end

sgtitle(['Regional integrated significance (rolling median) - ' ...
    strrep(comp_tag, '_', ' ')], 'FontSize', 14, 'FontWeight', 'bold');

% save it
set(fig_surp, 'InvertHardcopy', 'off');
saveas(fig_surp, ...
    fullfile(comp_out_dir, ['Region_Surprise_Bar_DiffSum_' comp_tag '.fig']));
exportgraphics(fig_surp, ...
    fullfile(comp_out_dir, ['Region_Surprise_Bar_DiffSum_' comp_tag '.png']), ...
    'Resolution', 300);

fprintf('Regional surprise analysis (Diff & Sum) saved to: %s\n', comp_out_dir);

% clearing these is not needed, since the workspace goes when the function returns
% clear roi_masks_surp surp_vec valid_pixels pixel_ids vol_surp
end

function roi_surp_agg = summed_surprise(vol_surp, valid_pixels, n_rois_surp, ...
    roi_pixel_counts_surp, roi_masks_surp, surp_thresh_val)
% Each region's sum of the surprise above surp_thresh_val (NaN counts as 0); NaN
% for a region with no voxel.

% each region's sum of the surprise above the threshold (NaN counts as 0)
surp_vec = vol_surp(valid_pixels);
surp_vec(isnan(surp_vec)) = 0;
roi_surp_agg = nan(n_rois_surp, 1);
for r = 1:n_rois_surp
    if roi_pixel_counts_surp(r) == 0
        continue;
    end
    vals = surp_vec(roi_masks_surp(:, r));
    if ~isempty(vals)
        roi_surp_agg(r) = sum(vals(vals > surp_thresh_val));
    end
end
end

function draw_surprise_bars(sorted_surp, sorted_rois_surp, surp_thresh_val, mode_name)
% The bars of one panel, grey by their sum relative to the largest.

b = barh(sorted_surp);
b.FaceColor = 'flat';

% bar colour by its sum relative to the largest: light grey (0.78) to black
c_map_surp = sep_palette('bars');
if ~isempty(sorted_surp)
    c_vals = round((sorted_surp / max(sorted_surp)) * size(c_map_surp, 1));
    c_vals(c_vals < 1) = 1;
    c_vals(isnan(c_vals)) = 1;
    for k = 1:length(sorted_surp)
        b.CData(k, :) = c_map_surp(c_vals(k), :);
    end
end

yticks(1:length(sorted_rois_surp));
yticklabels(sorted_rois_surp);
xlabel(['Aggregated surprise ( > ' num2str(surp_thresh_val, '%.1f') ')']);
title(['Regional ' lower(mode_name) ' significance']);
grid on;
set(gca, 'FontSize', 10);
ylim([0 length(sorted_rois_surp)+1]);
end

function highlight_surprise_regions(exp_type)
% The tick labels of the regions expected to change, in bold magenta.

% the regions expected to change, labelled in bold magenta
switch exp_type
    case {'rws', 'behavior'}
        highlighted_areas = {'Primary somatosensory area, barrel field', ...
            'Ventral posteromedial nucleus of the thalamus', ...
            'Posterior complex of the thalamus', ...
            'Supplemental somatosensory area', 'Zona incerta', ...
            'Rostrolateral visual area'};
    otherwise
        highlighted_areas = {};
end

ax = gca;
ytl = ax.YTickLabel;
colored_labels = repmat({''}, size(ytl));
for i = 1:length(ytl)
    if ismember(ytl{i}, highlighted_areas)
        colored_labels{i} = ['\color{magenta} \bf ' strrep(ytl{i}, '_', ' ')];
    else
        colored_labels{i} = ['\color{black} ' strrep(ytl{i}, '_', ' ')];
    end
end
ax.YTickLabel = colored_labels;
end

function vol_surp = rolling_surprise_median(raw_surp_vol, slab_range, ...
    brainMask_group_diff, mode_name)
% The surprise, median over a rolling slab of planes, inside the tissue of both
% groups.

fprintf('  [%s] Calculating rolling median (slab +/- %d)...\n', mode_name, slab_range);
vol_surp = zeros(size(raw_surp_vol), 'single');
n_slices = size(raw_surp_vol, 1);
for z = 1:n_slices
    z_start = max(1, z - slab_range);
    z_end = min(n_slices, z + slab_range);
    slab_data = raw_surp_vol(z_start:z_end, :, :);

    % NaN outside the voxels both groups have
    if exist('brainMask_group_diff', 'var')
        slab_mask = brainMask_group_diff(z_start:z_end, :, :);
        slab_data(~slab_mask) = NaN;
    end
    vol_surp(z, :, :) = nanmedian(slab_data, 1);
end
end
