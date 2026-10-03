%% run_group_differences
% ===== Left-right differences, control against experimental group =====
%
% Plasticity comparison, step 3 of 3:
%   1. run_collect_by_group    stack each group's registered volumes
%   2. run_normalise_groups    bring the mice of a group onto one intensity scale
%   3. run_group_differences   left-right differences, control against experimental group  <- this script
%
% Compares the hemispheric asymmetry of the signal between a control group
% (naive) and an experimental one (rws or behavior), voxel by voxel, from the
% normalised stacks of run_normalise_groups. Each mouse's mean tissue intensity
% per plane is taken, and its tissue smoothed in 3D (Gaussian, NaN-tolerant) if
% apply_smoothing is set. The experimental group is aligned onto the control one
% by a line fitted between the two groups' mean profiles over planes 200-700,
% then both are put on one common scale from planes 300-500. The hemispheres are
% folded: for every mouse |L - R| and L + R, then the group means, their
% difference, and Welch t and surprise (-log10 p) maps of it. Saves, in
% data\comparisons\<ctrl>_vs_<exp>_<channel>\, the profile alignment figure, the
% slab figures around plane 565 (group t-maps masked by surprise, individual
% mice), the regional surprise bars and the videos switched on below.
%
% Setup: naive against behavior (the four behavior mice), nano channel, smoothed
% with sigma 5, every video on, the region analyses off. Run sep_setup_paths
% first, once per MATLAB session; the code is in pipeline\group_differences.m.

clear; clc; close all;

%% Settings

% project folders, worked out from where the code sits, so the tree can be moved
% or copied to another drive as is (SEP_DATA_ROOT points the data elsewhere)
paths = get_paths();

% the adult mice, the group of each, and the groups in the order of
% selected_mice_idx_list
mice = {'MG691_Gria1', 'MG692_Gria1', 'MG693_Gria1', 'MG736_Gria1', 'MG737_Gria1', ...
    'CGF027_Gria1', 'CGF028_Gria1', 'CGF033_Gria1', 'CGF034_Gria1', 'CGF035_Gria1', ...
    'MG705_Gria1', 'MG706_Gria1', 'MG709_Gria1', 'MG716_Gria1', 'MG718_Gria1', ...
    'MG725_Gria1', 'MG727_Gria1'};
mousetypes = {'rws', 'rws', 'rws', 'rws', 'rws', ...
    'naive', 'naive', 'naive', 'naive', 'naive', ...
    'behavior', 'behavior', 'behavior', 'behavior', 'behavior', 'behavior', 'behavior'};
mousetypes_list = {'rws', 'naive', 'behavior'};

% the two groups compared: the control group, and the experimental one ('rws' or
% 'behavior')
ctrl_type = 'naive';
exp_type = 'behavior';

% the mice run_normalise_groups kept in each group, by position within the group,
% to name the mice in the figures and videos (must match its selection)
selected_mice_idx_list{1} = 1:5;
selected_mice_idx_list{2} = 1:5;
selected_mice_idx_list{3} = [1,3,4,5];

% mice of the behavior volume to analyse, by position along its fourth
% dimension: run_normalise_groups saves the four behavior mice selected above
% (MG705, MG709, MG716, MG718), and all four are used. Applied only when
% exp_type is 'behavior'; the other groups use every mouse saved.
behavior_subset = [1, 2, 3, 4];

% the names of the mice of the two groups
ctrl_group_idx = find(strcmp(mousetypes_list, ctrl_type));
all_ctrl_indices = find(strcmp(mousetypes, ctrl_type));
final_ctrl_indices = all_ctrl_indices(selected_mice_idx_list{ctrl_group_idx});
ctrl_mousenames = mice(final_ctrl_indices);
exp_group_idx = find(strcmp(mousetypes_list, exp_type));
all_exp_indices = find(strcmp(mousetypes, exp_type));
final_exp_indices = all_exp_indices(selected_mice_idx_list{exp_group_idx});
exp_mousenames = mice(final_exp_indices);

% videos to write: the group means, every mouse, the group t-scores, the surprise,
% the rolling slab, the signed differences
generate_diff_videos = true;
generate_individual_diff_videos = true;
generate_t_scored_videos = true;
generate_surprise_videos = true;
generate_rolling_videos = true;
generate_signed_diff_videos = true;

% region analyses, off in production: fine (every leaf region, montages and an
% annotated video) and coarse (a fixed list of regions, bar charts); either runs
% without the other
perform_area_based_analysis_fine = false;
perform_area_based_analysis_coarse = false;

% smooth each mouse's tissue in 3D before the comparison, with a NaN-tolerant
% Gaussian of sigma smooth_sigma, in 10 um voxels
apply_smoothing = true;
smooth_sigma = 5.0;

% channel to compare ('nano', surface GluA1, or 'auto', the autofluorescence
% control): loads <channel>_4d_normalized.mat from each cohort folder
channel = 'nano';

% comparison tag for the file names; it holds the channel, so the nano and auto
% outputs go to separate folders under comparisons\
comp_tag = [ctrl_type '_vs_' exp_type '_' channel];

% the two groups' folders, and the comparison's own folder for the results
base_root = paths.data;
ctrl_dir = fullfile(base_root, ctrl_type);
exp_dir = fullfile(base_root, exp_type);
comp_out_dir = fullfile(base_root, 'comparisons', comp_tag);
if ~exist(comp_out_dir, 'dir')
    mkdir(comp_out_dir);
end

%% Run

% pass the settings to the code, under the same names
run_settings = struct();
run_settings.paths = paths;
run_settings.ctrl_type = ctrl_type;
run_settings.exp_type = exp_type;
run_settings.ctrl_mousenames = ctrl_mousenames;
run_settings.exp_mousenames = exp_mousenames;
run_settings.behavior_subset = behavior_subset;
run_settings.generate_diff_videos = generate_diff_videos;
run_settings.generate_individual_diff_videos = generate_individual_diff_videos;
run_settings.generate_t_scored_videos = generate_t_scored_videos;
run_settings.generate_surprise_videos = generate_surprise_videos;
run_settings.generate_rolling_videos = generate_rolling_videos;
run_settings.generate_signed_diff_videos = generate_signed_diff_videos;
run_settings.perform_area_based_analysis_fine = perform_area_based_analysis_fine;
run_settings.perform_area_based_analysis_coarse = perform_area_based_analysis_coarse;
run_settings.apply_smoothing = apply_smoothing;
run_settings.smooth_sigma = smooth_sigma;
run_settings.channel = channel;
run_settings.comp_tag = comp_tag;
run_settings.ctrl_dir = ctrl_dir;
run_settings.exp_dir = exp_dir;
run_settings.comp_out_dir = comp_out_dir;
group_differences(run_settings);
