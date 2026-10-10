%% run_scale_free_test
% ===== Step 3's region test on the raw asymmetry index, which no scale changes =====
%
% Plasticity comparison, after its step 3:
%   1. run_collect_by_group    stack each group's registered volumes
%   2. run_normalise_groups    bring the mice of a group onto one intensity scale
%   3. run_group_differences   left-right differences, control against experimental group
%   4. run_per_mouse_values    each mouse's values in regions named in advance
%   5. run_mouse_influence     each mouse's part in the barrel-field cluster
%   6. run_scale_free_test     step 3's region test on a reading no scale changes  <- this script
%
% Step 3 puts the experimental group on the control group's scale by a line,
% whose slope multiplies every experimental mouse's |L - R| on its maps, and
% its splits keep that slope. run_mouse_influence found that the barrel-field
% result after RWS depends on it: at the fitted slope (1.49, set mostly by one
% bright naive mouse, CGF033; jackknife standard error 0.45) the RWS-higher
% cluster is the heavier sign (p 0.040), at slope 1 the naive-higher one is
% (its p 0.004), and the signs cross between 1.3 and 1.4. This runs step 3's
% region test, unchanged, on a reading no scale can change: per mouse and
% voxel, the asymmetry index |L - R| / (L + R) on the collected stack less the
% mouse's off-tissue level, smoothed as step 3 smooths before the ratio. No
% normalisation of step 2 and no alignment of the groups: a factor on a mouse
% or on a group divides out of the ratio.
%
% The design was fixed on 10 October 2026 before any of its numbers, and is
% reported whatever comes out; no other variant is tried (README, the
% scale-free test). The test named in advance: the barrel field's heaviest
% cluster where the index is higher after RWS, its mass against that of every
% one of the 252 splits (one-sided, Gambino et al. 2014), the p of either sign
% beside it; and every region's score, p and corrected p as step 3 gives them.
% After behavior the same over its 126 splits, the p of either sign first.
% Saves, in data\comparisons\<ctrl>_vs_<exp>_<channel>\,
% Scale_Free_Regions_<tag>.csv, Scale_Free_Clusters_<tag>.csv,
% Scale_Free_Bars_<tag> and Scale_Free_TMap_<tag> (.fig and .png), and the
% cache of the test, Scale_Free_<tag>.mat, <tag> being the comparison and the
% smoothing.
%
% Setup: naive against rws, nano channel, smoothed with sigma 5 as step 3, the
% statistics of step 3's region test, the barrel field named in advance; for
% naive against behavior set exp_type = 'behavior'. run_per_mouse_values must
% have run with the same comparison and smoothing: its leave-one-out file holds
% step 3's cluster voxel for voxel, and its maps check the reading. Run
% sep_setup_paths first, once per MATLAB session; the code is in
% pipeline\scale_free_test.m.

clear; clc; close all;

%% Settings

% project folders, worked out from where the code sits, so the tree can be moved
% or copied to another drive as is (SEP_DATA_ROOT points the data elsewhere)
paths = get_paths();

% the two groups compared: the control group, and the experimental one ('rws' or
% 'behavior')
ctrl_type = 'naive';
exp_type = 'rws';

% behavior mice to analyse, by name, as in run_group_differences; applied only
% when exp_type is 'behavior'
behavior_mice = {'MG705_Gria1', 'MG709_Gria1', 'MG716_Gria1', 'MG718_Gria1'};

% the reading, fixed on 10 October 2026 before any number: per mouse and voxel
% the asymmetry index |L - R| / (L + R) of the collected stack, the mouse's
% off-tissue level subtracted (as run_per_mouse_values' raw stack), its tissue
% smoothed in 3D before the folding exactly as step 3 smooths it (NaN-aware,
% sigma 5 in 10 um voxels), missing where L + R is not above 0 or outside the
% mouse's tissue; the folded grid, the regions and the masks are step 3's
apply_smoothing = true;
smooth_sigma = 5.0;

% the statistics, step 3's region test as it is (fixed with the reading): the
% Welch t where each group has at least min_mice_per_group mice with a value,
% its surprise -log10 p, the median over +/- slab_range planes, voxels at
% p < cluster_p touching by a face or an edge, clusters within each region's
% own voxels; each region scored by its heaviest cluster's mass with its sign,
% the heavier of the two signs, its p over every split of the mice
% (n_permutations = 'all': 252 for 5 against 5, 126 for 5 against 4), its
% corrected p the largest score over the regions
min_mice_per_group = 3;
slab_range = 10;
cluster_p = 0.01;
cluster_connectivity = 18;
n_permutations = 'all';

% the region named in advance: the barrel field. For RWS its direction too,
% the index higher after RWS (RWS potentiates the stimulated barrels' synapses
% and brings AMPA receptors to their surface, Gambino et al. 2014): the
% one-sided p, the share of the splits whose cluster of that sign is as heavy,
% is the test, the p of either sign beside it. For behavior the direction is
% carried over from RWS, not named on its own, and the p of either sign comes
% first
a_priori_region = 'SSp-bfd';
direction_named = strcmp(exp_type, 'rws');

% the colour limits of the t map, those of step 3's t videos
t_limit = 6;

% thread workers of the permutation test, as in run_group_differences
permutation_workers = 16;

% redo the maps and the test rather than read them from the cache (the cache is
% read only when made from these stacks and settings)
force_recompute = false;

% channel compared, as in run_group_differences
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
run_settings.apply_smoothing = apply_smoothing;
run_settings.smooth_sigma = smooth_sigma;
run_settings.min_mice_per_group = min_mice_per_group;
run_settings.slab_range = slab_range;
run_settings.cluster_p = cluster_p;
run_settings.cluster_connectivity = cluster_connectivity;
run_settings.n_permutations = n_permutations;
run_settings.a_priori_region = a_priori_region;
run_settings.direction_named = direction_named;
run_settings.t_limit = t_limit;
run_settings.permutation_workers = permutation_workers;
run_settings.force_recompute = force_recompute;
run_settings.channel = channel;
run_settings.comp_tag = comp_tag;
run_settings.ctrl_dir = ctrl_dir;
run_settings.exp_dir = exp_dir;
run_settings.comp_out_dir = comp_out_dir;
scale_free_test(run_settings);
