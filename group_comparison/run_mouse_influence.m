%% run_mouse_influence
% ===== Each mouse's part in the barrel-field cluster =====
%
% Plasticity comparison, after its step 3:
%   1. run_collect_by_group    stack each group's registered volumes
%   2. run_normalise_groups    bring the mice of a group onto one intensity scale
%   3. run_group_differences   left-right differences, control against experimental group
%   4. run_per_mouse_values    each mouse's values in regions named in advance
%   5. run_mouse_influence     each mouse's part in the barrel-field cluster  <- this script
%
% Which mice the barrel-field cluster of step 3's test rests on, and why a few
% mice are enough for its test while their means over the cluster do not
% separate the groups. Five parts, from the caches of run_per_mouse_values:
%   1. each mouse left out in turn, the comparison redone without it as step 3
%      does it (the alignment refitted, the t with min_mice_per_group mice per
%      group, the surprise, its rolling median, p < cluster_p, connected within
%      the region): the heaviest cluster where |L - R| is higher in the
%      experimental group, its voxels, mass and overlap with the cluster of all
%      the mice, and its p over every split of the remaining mice (126 for five
%      against four), for that sign and for either sign as step 3 takes it;
%      then the same with the alignment held at that of all the mice, which
%      tells a mouse's own values from its part in the alignment's slope;
%   2. the comparison of all the mice redone at other slopes of the alignment
%      (sweep_slopes): the slope multiplies every experimental mouse's |L - R|
%      on the test's maps, and the splits keep it fixed;
%   3. each mouse's asymmetry, |L - R| / (L + R) on its raw stack less its
%      off-tissue level (the reading to look at), on its autofluorescence and,
%      as |L - R|, on the test's maps, in a window around the cluster: a coronal
%      view over the cluster's planes and a view from above over the cluster's
%      depth below the pia, the cluster outlined; and along AP through the
%      cluster, with the edges between the mouse's sections;
%   4. at every voxel of the region, how many mice of the experimental group are
%      more asymmetric than every mouse of the control group (0 to 5): how many
%      voxels reach each count, and the largest connected patch at each count,
%      against the same under every split of the mice; on the raw stack, the
%      test's maps, the test's maps at slope 1 and the autofluorescence;
%   5. the cluster's size in micrometres against one barrel column, its depth
%      below the pia, and where each mouse's own peak asymmetry sits.
% Saves, in data\comparisons\<ctrl>_vs_<exp>_<channel>\, Influence_Folds_<tag>
% (.csv and the figure of the folds and the slopes), Influence_Slopes_<tag>.csv,
% Influence_Maps_<tag>_<reading> (one figure per reading),
% Influence_Profiles_<tag> and Influence_Consistency_<tag> (.csv and the
% figure), Influence_Mice_<tag>.csv, Influence_Cluster_<tag>.csv and the cache
% of the splits, Influence_<tag>.mat, <tag> being the comparison and the
% smoothing.
%
% Setup: naive against rws, nano channel, smoothed with sigma 5, the barrel
% field, the cluster settings of step 3; run_per_mouse_values must have run with
% the same settings (its maps and autofluorescence caches are read). Run
% sep_setup_paths first, once per MATLAB session; the code is in
% pipeline\mouse_influence.m (for ten mice and no pool, about 25 minutes for the
% folds under every split, as much again with the alignment held, 20 for the
% slopes and 4 for the counts; a few minutes from the cache).

clear; clc; close all;

%% Settings

% project folders, worked out from where the code sits, so the tree can be moved
% or copied to another drive as is (SEP_DATA_ROOT points the data elsewhere)
paths = get_paths();

% the two groups compared: the control group, and the experimental one ('rws' or
% 'behavior')
ctrl_type = 'naive';
exp_type = 'rws';

% behavior mice to analyse, by name, as in run_per_mouse_values; applied only
% when exp_type is 'behavior'
behavior_mice = {'MG705_Gria1', 'MG709_Gria1', 'MG716_Gria1', 'MG718_Gria1'};

% the regions and the clusters' region of run_per_mouse_values, whose caches are
% read: they must be the same
regions = {'SSp-bfd', 'VISp'};
cluster_region = 'SSp-bfd';

% smoothing of each mouse's tissue, as in run_per_mouse_values
apply_smoothing = true;
smooth_sigma = 5.0;

% the clusters' test, as run_group_differences and run_per_mouse_values do it: a
% t with at least min_mice_per_group mice per group, the surprise's median over
% +/- slab_range planes, voxels at p < cluster_p, touching by a face or an edge;
% the patches of the counts are connected the same way
min_mice_per_group = 3;
slab_range = 10;
cluster_p = 0.01;
cluster_connectivity = 18;

% the slopes of the alignment the comparison of all the mice is redone at,
% besides the fitted one: from 1, each group on its own scale of step 2, to 1.7,
% past the slopes refitted without one mouse (1.11 to 1.67 for naive against
% rws); the alignment's intercept and common factor change no t
sweep_slopes = 1.0:0.1:1.7;

% the window of the maps: the cluster's bounding box and this much on each side,
% in micrometres (600 um: about two barrel columns of 300 um, Lefort et al.
% 2009, beyond the cluster's edge on each side)
window_um = 600;

% workers of the region test's thread pool under every split (0: no pool)
n_workers = 0;

% redo the splits rather than read them from the cache (the cache is read only
% when made from these maps and settings)
force_recompute = false;

% channel compared, as in run_per_mouse_values; the autofluorescence is read
% beside it
channel = 'nano';

% comparison tag and folders, as in run_per_mouse_values: the results go into
% the comparison's folder, beside its caches
comp_tag = [ctrl_type '_vs_' exp_type '_' channel];
base_root = paths.data;
ctrl_dir = fullfile(base_root, ctrl_type);
exp_dir = fullfile(base_root, exp_type);
comp_out_dir = fullfile(base_root, 'comparisons', comp_tag);

%% Run

% pass the settings to the code, under the same names
run_settings = struct();
run_settings.paths = paths;
run_settings.ctrl_type = ctrl_type;
run_settings.exp_type = exp_type;
run_settings.behavior_mice = behavior_mice;
run_settings.regions = regions;
run_settings.cluster_region = cluster_region;
run_settings.apply_smoothing = apply_smoothing;
run_settings.smooth_sigma = smooth_sigma;
run_settings.min_mice_per_group = min_mice_per_group;
run_settings.slab_range = slab_range;
run_settings.cluster_p = cluster_p;
run_settings.cluster_connectivity = cluster_connectivity;
run_settings.sweep_slopes = sweep_slopes;
run_settings.window_um = window_um;
run_settings.n_workers = n_workers;
run_settings.force_recompute = force_recompute;
run_settings.channel = channel;
run_settings.comp_tag = comp_tag;
run_settings.ctrl_dir = ctrl_dir;
run_settings.exp_dir = exp_dir;
run_settings.comp_out_dir = comp_out_dir;
mouse_influence(run_settings);
