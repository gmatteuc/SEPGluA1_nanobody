%% run_per_mouse_values
% ===== Each mouse's values in the regions named in advance =====
%
% Plasticity comparison, after its step 3:
%   1. run_collect_by_group    stack each group's registered volumes
%   2. run_normalise_groups    bring the mice of a group onto one intensity scale
%   3. run_group_differences   left-right differences, control against experimental group
%   4. run_per_mouse_values    each mouse's values in regions named in advance  <- this script
%
% The voxel maps of step 3 give a group result; this gives one value per mouse
% in regions named before looking (the barrel field, and the primary visual
% area as a control), to show which mice carry an effect and whether every
% mouse of a group does. Each mouse is taken as step 3 takes it (its tissue
% smoothed, folded onto the left hemisphere, the experimental group aligned
% onto the control group), and in each region gets its asymmetry index, mean
% |L - R| over mean (L + R), on those maps and on its collected stack less its
% off-tissue level (no normalisation, no alignment); the signed index, for the
% side of the asymmetry; and its L + R relative to its own isocortex. Then each
% mouse is left out in turn, the comparison is done again without it as step
% 3 does it, and the mouse's index is read in the heaviest positive cluster of
% the barrel field's L - R map that the others give, so no mouse is read in a
% cluster its own data helped define. Every value is compared between the
% groups by an exact permutation of the mice (252 splits for 5 against 5), with
% the Welch t and Hedges' g beside it. Saves, in
% data\comparisons\<ctrl>_vs_<exp>_<channel>\, Per_Mouse_Values_<tag> (.csv
% and the figure, .fig and .png), Per_Mouse_Stats_<tag>.csv,
% Per_Mouse_LOO_<tag>.csv and the cache of the mice's maps,
% Per_Mouse_Maps_<tag>.mat, <tag> being the comparison and the smoothing.
%
% Setup: naive against rws, nano channel, smoothed with sigma 5 as step 3,
% the barrel field and the primary visual area, the cluster settings of step
% 3. For naive against behavior set exp_type = 'behavior'. Run sep_setup_paths
% first, once per MATLAB session; the code is in
% pipeline\per_mouse_region_values.m (about 45 minutes for ten mice, a minute
% from the cache).

clear; clc; close all;

%% Settings

% project folders, worked out from where the code sits, so the tree can be moved
% or copied to another drive as is (SEP_DATA_ROOT points the data elsewhere)
paths = get_paths();

% the two groups compared: the control group, and the experimental one ('rws' or
% 'behavior')
ctrl_type = 'naive';
exp_type = 'rws';

% behavior mice to analyse, by name, among the four run_normalise_groups saved,
% as in run_group_differences; applied only when exp_type is 'behavior'
behavior_mice = {'MG705_Gria1', 'MG709_Gria1', 'MG716_Gria1', 'MG718_Gria1'};

% regions named in advance, by atlas acronym (regions of the bars of step 3): the
% barrel field, since RWS stimulates the whiskers and the behavior task uses one
% whisker, and the primary visual area as a cortical control
regions = {'SSp-bfd', 'VISp'};

% the region of the leave-one-out clusters, one of regions
loo_region = 'SSp-bfd';

% smoothing of each mouse's tissue, as in run_group_differences (the test's maps)
apply_smoothing = true;
smooth_sigma = 5.0;

% the leave-one-out's test, as run_group_differences and its bars do it: a t
% with at least min_mice_per_group mice per group, the surprise's median over
% +/- slab_range planes, voxels at p < cluster_p, touching by a face or an edge
min_mice_per_group = 3;
slab_range = 10;
cluster_p = 0.01;
cluster_connectivity = 18;

% redo every mouse's maps rather than read them from the cache, after a change
% to the stacks or to how a mouse is processed
force_recompute_mice = false;

% channel to compare ('nano' or 'auto'), as in run_group_differences
channel = 'nano';

% comparison tag and folders, as in run_group_differences: the results go into
% the comparison's folder, beside those of step 3
comp_tag = [ctrl_type '_vs_' exp_type '_' channel];
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
run_settings.behavior_mice = behavior_mice;
run_settings.regions = regions;
run_settings.loo_region = loo_region;
run_settings.apply_smoothing = apply_smoothing;
run_settings.smooth_sigma = smooth_sigma;
run_settings.min_mice_per_group = min_mice_per_group;
run_settings.slab_range = slab_range;
run_settings.cluster_p = cluster_p;
run_settings.cluster_connectivity = cluster_connectivity;
run_settings.force_recompute_mice = force_recompute_mice;
run_settings.channel = channel;
run_settings.comp_tag = comp_tag;
run_settings.ctrl_dir = ctrl_dir;
run_settings.exp_dir = exp_dir;
run_settings.comp_out_dir = comp_out_dir;
per_mouse_region_values(run_settings);
