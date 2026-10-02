%% run_normalise_groups
% ===== Bring the mice of each group onto one intensity scale =====
%
% Plasticity comparison, step 2 of 3:
%   1. run_collect_by_group    stack each group's registered volumes
%   2. run_normalise_groups    bring the mice of a group onto one intensity scale  <- this script
%   3. run_group_differences   left-right differences, control against experimental group
%
% Loads a group's stack from run_collect_by_group (<channel>_4d<tag>.mat), keeps
% the selected mice, and recomputes a background mask for every plane of every
% mouse, with percentile bounds that change along the AP axis. The cortical
% tissue of every fifth plane (isocortex layers 1 to 5, without the
% retrosplenial, anterior cingulate and prelimbic areas) is pooled, its median
% across mice taken as the reference, and each mouse fitted onto it with a
% robust line (slope and intercept). (raw - intercept) / slope is then applied
% to the whole volume, the voxels no section reached are set to NaN, and the
% result is saved for run_group_differences as <channel>_4d_normalized<tag>.mat,
% with the background masks in <channel>_4d_normalized_bkgmask<tag>.mat.
% Diagnostic figures go to the group's folder (the background trace) and to
% global_diagnostics\normalization_checks_<channel>\, with a video if asked for.
%
% Setup: the plasticity comparison, the rws and naive groups, nano channel
% (behavior and the young cohorts through SEP_COHORT_SPECS). Run sep_setup_paths
% first, once per MATLAB session; the code is in pipeline\normalise_groups.m.

clear; clc; close all;

%% Settings

% project folders, worked out from where the code sits, so the tree can be moved
% or copied to another drive as is (SEP_DATA_ROOT points the data elsewhere)
paths = get_paths();

% the adult mice, in the order of their groups' stacks, the group of each, and
% the groups in the order of selected_mice_idx_list
mice = {'MG691_Gria1', 'MG692_Gria1', 'MG693_Gria1', 'MG736_Gria1', 'MG737_Gria1', ...
    'CGF027_Gria1', 'CGF028_Gria1', 'CGF033_Gria1', 'CGF034_Gria1', 'CGF035_Gria1', ...
    'MG705_Gria1', 'MG706_Gria1', 'MG709_Gria1', 'MG716_Gria1', 'MG718_Gria1', ...
    'MG725_Gria1', 'MG727_Gria1'};
mousetypes = {'rws', 'rws', 'rws', 'rws', 'rws', ...
    'naive', 'naive', 'naive', 'naive', 'naive', ...
    'behavior', 'behavior', 'behavior', 'behavior', 'behavior', 'behavior', 'behavior'};
mousetypes_list = {'rws', 'naive', 'behavior'};

% mice kept in each adult group, by position within the group: rws and naive all
% five; behavior four of seven (MG705, MG709, MG716, MG718), three excluded as bad
selected_mice_idx_list{1} = 1:5;
selected_mice_idx_list{2} = 1:5;
selected_mice_idx_list{3} = [1,3,4,5];

% draw a video of every plane, normalised against raw (slow)
plot_verification_video = false;

% channel to normalise ('nano', surface GluA1, or 'auto', the autofluorescence
% control): loads <channel>_4d.mat and saves <channel>_4d_normalized.mat in each
% cohort folder, for run_group_differences, P8 and P9
channel = 'nano';

% cohorts to normalise, as specs: an adult group by name ('rws', 'naive',
% 'behavior'), with the mouse lists and selections above, or the young group with
% an age ('young_P20'), with all the mice run_collect_by_group stacked
% (collected_mice<tag>.mat) and its own atlas; see get_cohort_spec
cohort_specs = {'rws', 'naive'};

% batch runs pick the cohorts without editing this file:
% set SEP_COHORT_SPECS=young_P20 (comma-separated for several)
if ~isempty(getenv('SEP_COHORT_SPECS'))
    cohort_specs = strtrim(strsplit(getenv('SEP_COHORT_SPECS'), ','));
end

%% Run

% pass the settings to the code, under the same names
run_settings = struct();
run_settings.paths = paths;
run_settings.mice = mice;
run_settings.mousetypes = mousetypes;
run_settings.mousetypes_list = mousetypes_list;
run_settings.selected_mice_idx_list = selected_mice_idx_list;
run_settings.plot_verification_video = plot_verification_video;
run_settings.channel = channel;
run_settings.cohort_specs = cohort_specs;
normalise_groups(run_settings);
