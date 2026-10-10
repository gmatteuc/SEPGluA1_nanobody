%% run_rescaling_figure
% ===== The barrel-field cluster test at each rescaling of the RWS mice =====
%
% Plasticity comparison, after its step 3:
%   1. run_collect_by_group    stack each group's registered volumes
%   2. run_normalise_groups    bring the mice of a group onto one intensity scale
%   3. run_group_differences   left-right differences, control against experimental group
%   4. run_per_mouse_values    each mouse's values in regions named in advance
%   5. run_mouse_influence     each mouse's part in the barrel-field cluster
%   6. run_scale_free_test     step 3's region test on a reading no scale changes
%   7. run_rescaling_figure    a figure of the cluster test at each rescaling  <- this script
%
% Step 3 puts the experimental mice on the control mice's scale with a line
% whose slope multiplies every experimental mouse's |L - R|, and its splits
% keep that slope. One figure of what steps 5 and 6 found about it, read from
% their tables, nothing computed again: against the slope, the one-sided p of
% the region's heaviest cluster of each sign (the experimental-higher one the
% test the scale-free test named in advance) and each cluster's size, with the
% fitted slope, its jackknife standard error and the slope refitted without
% the mouse that moves it most; under the title, the experimental-higher
% cluster's p at the fitted slope, one-sided and of either sign (step 3's p,
% the one quoted), and on the asymmetry index |L - R| / (L + R), where no
% scale enters. Saves, in data\comparisons\<ctrl>_vs_<exp>_<channel>\,
% Rescaling_Test_<tag> (.fig and .png), <tag> being the comparison and the
% smoothing.
%
% Setup: naive against rws, nano channel, smoothed with sigma 5, the barrel
% field; run_mouse_influence and run_scale_free_test must have run with the
% same comparison and smoothing (Influence_Slopes_<tag>.csv,
% Influence_Folds_<tag>.csv, Scale_Free_Clusters_<tag>.csv are read). Run
% sep_setup_paths first, once per MATLAB session; the code is in
% pipeline\plot_rescaling_test.m (a few seconds).

clear; clc; close all;

%% Settings

% project folders, worked out from where the code sits, so the tree can be moved
% or copied to another drive as is (SEP_DATA_ROOT points the data elsewhere)
paths = get_paths();

% the two groups compared: the control group, and the experimental one; the
% influence of each mouse was run for 'rws' only
ctrl_type = 'naive';
exp_type = 'rws';

% the names the figure's text gives the two groups, and the region of the
% clusters (run_mouse_influence's cluster_region and run_scale_free_test's
% a_priori_region, SSp-bfd)
ctrl_name = 'naive';
exp_name = 'RWS';
region_name = 'barrel field';

% smoothing of the tables read, as in run_mouse_influence and run_scale_free_test
apply_smoothing = true;
smooth_sigma = 5.0;

% channel compared, as in run_mouse_influence
channel = 'nano';

% comparison tag and folder, as in run_mouse_influence: the figure goes into the
% comparison's folder, beside the tables it is drawn from
comp_tag = [ctrl_type '_vs_' exp_type '_' channel];
comp_out_dir = fullfile(paths.data, 'comparisons', comp_tag);

%% Run

% pass the settings to the code, under the same names
run_settings = struct();
run_settings.exp_type = exp_type;
run_settings.ctrl_name = ctrl_name;
run_settings.exp_name = exp_name;
run_settings.region_name = region_name;
run_settings.apply_smoothing = apply_smoothing;
run_settings.smooth_sigma = smooth_sigma;
run_settings.comp_tag = comp_tag;
run_settings.comp_out_dir = comp_out_dir;
plot_rescaling_test(run_settings);
