%% run_per_mouse_values
% ===== Each mouse's values in the regions named in advance =====
%
% Plasticity comparison, after its step 3:
%   1. run_collect_by_group    stack each group's registered volumes
%   2. run_normalise_groups    bring the mice of a group onto one intensity scale
%   3. run_group_differences   left-right differences, control against experimental group
%   4. run_per_mouse_values    each mouse's values in regions named in advance  <- this script
%
% The voxel maps of step 3 give a group result; this gives one value per mouse,
% to show which mice carry an effect and whether every mouse of a group does.
% Each mouse is taken as step 3 takes it (its tissue smoothed, folded onto the
% left hemisphere, the experimental group aligned onto the control group). Its
% asymmetry index, mean |L - R| over mean (L + R), is read first in the
% barrel-field cluster of step 3's test, the heaviest where |L - R| is higher
% in the experimental group (a positive t), found with all the mice: on its
% collected stack less its off-tissue level (no normalisation, no alignment),
% on step 3's maps and on its autofluorescence. That is step 3's test seen
% mouse by mouse, and its p redoes the cluster search under every split of the
% mice (the selection-matched test). Then over whole regions named before
% looking (the barrel field, and the primary visual area as a control): the
% index and the L + R relative to the mouse's own isocortex, their p the
% exact permutation of the values. No signed value: left and right are not
% certain for every brain, so the stimulated side is not known mouse by mouse.
% Last, stricter, each mouse is left out in turn, the comparison is done again
% without it, and its index is read in the cluster the others give, so no
% mouse is read in a cluster its own data helped define; its p redoes the
% leave-one-out under every split. Every p is one-sided for the experimental
% group higher and two-sided. Saves, in data\comparisons\<ctrl>_vs_<exp>_<channel>\,
% Per_Mouse_Values_<tag> (.csv and the figure, .fig and .png),
% Per_Mouse_Stats_<tag>.csv, Per_Mouse_LOO_<tag> (the leave-one-out's figure,
% .csv, and .mat with the folds' clusters) and the caches of the mice's maps,
% Per_Mouse_Maps_<tag>.mat, of their autofluorescence, Per_Mouse_Auto_<tag>.mat,
% of the selection-matched test, Per_Mouse_Selection_<tag>.mat, and of the
% redone leave-one-out, Per_Mouse_LOO_Relabelled_<tag>.mat, <tag> being the
% comparison and the smoothing.
%
% Setup: naive against rws, nano channel, smoothed with sigma 5 as step 3,
% the barrel field and the primary visual area, the cluster settings of step
% 3, the autofluorescence stacked by run_collect_by_group with channels =
% {'auto'}. For naive against behavior set exp_type = 'behavior'. Run
% sep_setup_paths first, once per MATLAB session; the code is in
% pipeline\per_mouse_region_values.m (for ten mice, about 25 minutes for the
% maps, 25 for the autofluorescence, 6 for the selection-matched test and an
% hour for the leave-one-out under every split, 25 minutes for 5 against 4; a
% few minutes from the caches).

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

% the region of the clusters, of the selection-matched test and of the
% leave-one-out, one of regions
cluster_region = 'SSp-bfd';

% smoothing of each mouse's tissue, as in run_group_differences (the test's maps)
apply_smoothing = true;
smooth_sigma = 5.0;

% the clusters' test, as run_group_differences and its bars do it: a t
% with at least min_mice_per_group mice per group, the surprise's median over
% +/- slab_range planes, voxels at p < cluster_p, touching by a face or an edge
min_mice_per_group = 3;
slab_range = 10;
cluster_p = 0.01;
cluster_connectivity = 18;

% whether the one-sided direction, the experimental group higher, was named
% before any number: for RWS it was (Giulio, 9 October 2026: RWS potentiates
% the stimulated barrels' synapses and brings AMPA receptors to their surface,
% Gambino et al. 2014), and its one-sided p comes first; for behavior it is
% carried over from RWS, not named on its own, and the two-sided p comes first
direction_named = strcmp(exp_type, 'rws');

% the selection-matched test, the main per-mouse test (fixed on 9 October 2026
% before any of its numbers; the README says why it is valid): step 3's test of
% the barrel field seen mouse by mouse, as permissive as it. With the true
% groups, the heaviest cluster of the region where |L - R| is higher in the
% experimental group (a positive t), found with all the mice as
% region_permutation_test finds it; every mouse's AI read in it, on the raw
% stack (the value to read), on the test's maps and on the autofluorescence;
% the statistic, the experimental group's mean minus the control group's. Its
% null: under every split of the mice (252 for 5 and 5, 126 for 5 and 4) the
% cluster found again with the split's groups, the split's experimental group
% in the experimental group's place, every mouse read in that split's cluster,
% the same statistic. One-sided p: the share of the splits reaching the
% observed statistic, the observed split included. Two-sided: the search in
% either direction, as the region test's score takes the larger of its two
% signs (per split, the larger of the statistic in its positive cluster and of
% control minus experimental in its negative one). A split without a cluster,
% or with a group without a mouse with a value, gives 0; the p with those
% splits left out is given beside it. Its first split must give step 3's
% cluster, and the cluster mass over its splits step 3's p
selection_matched = true;

% the autofluorescence read in the same clusters and tested the same way,
% against a misregistration of the surface: each group's auto_4d.mat, stacked
% from the registered volumes by run_collect_by_group with channels = {'auto'}
% (one without the files it was read from, as production's older auto_4d.mat,
% is refused), checked against the nano stack voxel for voxel (the voxels a
% section reached), less its off-tissue level and smoothed in the nano
% channel's tissue, as the raw stack
auto_control = true;

% redo the leave-one-out under every split of the mice (each fold's cluster
% found again with the split's groups), for the full p of its two values; read
% from its cache after the first run. Stricter than the selection-matched test:
% each mouse is read in a cluster found without it, four mice against five
loo_relabel = true;

% redo every mouse's maps rather than read them from the cache, after a change
% to how a mouse is processed (a change to the normalised stacks stops the run)
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
run_settings.cluster_region = cluster_region;
run_settings.apply_smoothing = apply_smoothing;
run_settings.smooth_sigma = smooth_sigma;
run_settings.min_mice_per_group = min_mice_per_group;
run_settings.slab_range = slab_range;
run_settings.cluster_p = cluster_p;
run_settings.cluster_connectivity = cluster_connectivity;
run_settings.direction_named = direction_named;
run_settings.selection_matched = selection_matched;
run_settings.auto_control = auto_control;
run_settings.loo_relabel = loo_relabel;
run_settings.force_recompute_mice = force_recompute_mice;
run_settings.channel = channel;
run_settings.comp_tag = comp_tag;
run_settings.ctrl_dir = ctrl_dir;
run_settings.exp_dir = exp_dir;
run_settings.comp_out_dir = comp_out_dir;
per_mouse_region_values(run_settings);
