clear all
close all
clc

% /// Plasticity comparison, step 2 of 3: normalise the mice of each group ///
% Brings the mice of a group onto a common intensity scale, so that they can
% be averaged and compared, and saves the normalised stack for
% run_group_differences:
%   (1) loads the group's stack from run_collect_by_group
%       (<channel>_4d<tag>.mat) and keeps the selected mice
%   (2) recomputes a background mask for every plane of every mouse, with
%       percentile bounds that change along the AP axis
%   (3) pools the cortical tissue pixels (isocortex layers 1 to 5, without
%       the retrosplenial, anterior cingulate and prelimbic areas) of every
%       fifth plane, takes their median across mice as the reference, and
%       fits each mouse onto it with a robust line (slope and intercept)
%   (4) draws diagnostic figures (the background trace in the group's
%       folder, the rest in global_diagnostics\normalization_checks_<channel>\),
%       and optionally a video
%   (5) applies (raw - intercept) / slope to the whole volume, sets the
%       voxels no section reached to NaN, and saves
%       <channel>_4d_normalized<tag>.mat and the background masks,
%       <channel>_4d_normalized_bkgmask<tag>.mat
% Run sep_setup_paths first, once per MATLAB session. The settings are below,
% the code is in pipeline\normalise_groups.m.

% Where the project lives. Derived from the location of the code rather than
% written out, so the tree can be moved or copied to another drive as is.
paths = get_paths();

%%  Set user-defined parameters

% Set mice and mousetypes
mice = {'MG691_Gria1', 'MG692_Gria1', 'MG693_Gria1', 'MG736_Gria1', 'MG737_Gria1', 'CGF027_Gria1', 'CGF028_Gria1', 'CGF033_Gria1', 'CGF034_Gria1', 'CGF035_Gria1','MG705_Gria1', 'MG706_Gria1', 'MG709_Gria1', 'MG716_Gria1', 'MG718_Gria1', 'MG725_Gria1', 'MG727_Gria1'};
mousetypes = {'rws','rws','rws','rws','rws','naive','naive','naive','naive','naive','behavior','behavior','behavior','behavior','behavior','behavior','behavior'};
mousetypes_list = {'rws','naive','behavior'};

% Set list of selected mice
selected_mice_idx_list{1} = 1:5; % rws
selected_mice_idx_list{2} = 1:5; % naive
selected_mice_idx_list{3} = [1,3,4,5]; % behavior (excluding 3 bad mice)

% Set plotting and saving
plot_verification_video = false;

% Channel to normalize: 'nano' (default, surface GluA1) or 'auto' (autofluorescence control).
% Loads <channel>_4d.mat and saves <channel>_4d_normalized.mat in each cohort folder.
% Pipeline scripts downstream (run_group_differences, P8, P9) read whichever
% channel they're configured for.
channel = 'nano'; % 'auto'

% Cohorts to normalise, as specs: an adult group by name ('rws', 'naive',
% 'behavior'), or the young group with an age ('young_P20'). An adult group
% keeps the mouse lists and selections above, exactly as before. A young cohort
% takes the mice run_collect_by_group stacked (collected_mice<tag>.mat), all of
% them, and its own atlas. See get_cohort_spec.
cohort_specs = {'rws', 'naive'};
% Batch runs can pick the cohorts without editing this file:
%   set SEP_COHORT_SPECS=young_P20   (comma-separated for several)
if ~isempty(getenv('SEP_COHORT_SPECS'))
    cohort_specs = strtrim(strsplit(getenv('SEP_COHORT_SPECS'), ','));
end

%% Run

% The settings above go to the code under the same names
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
