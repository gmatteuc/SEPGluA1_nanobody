%% run_collect_by_group
% ===== Collect the registered volumes of each group =====
%
% Plasticity comparison, step 1 of 3:
%   1. run_collect_by_group    stack each group's registered volumes  <- this script
%   2. run_normalise_groups    bring the mice of a group onto one intensity scale
%   3. run_group_differences   left-right differences, control against experimental group
%
% Stacks the registered nano volumes of a group's mice into one 4D array
% (AP x DV x ML x mouse), saved in the group's folder for run_normalise_groups.
% The mice come from the cohort registry (get_cohort), in registry order,
% narrowed to some ages if age_filter is set; mice whose registration has not
% finished are skipped with a note, so a stack can be rebuilt as more brains
% finish. Saves <group>\nano_4d<tag>.mat and collected_mice<tag>.mat, the
% order of the mice along the fourth dimension (<tag> is empty for a whole
% group, _P20 for age_filter = [20]). Only the nano stack is saved; the help
% of collect_by_group says why the other files on disk are not to be read.
%
% Setup: the young cohort, P20 brains only, for the comparison with the
% adults. For the plasticity comparison, the groups are rws, naive and
% behavior, with no age filter. Run sep_setup_paths first, once per MATLAB
% session; the code is in pipeline\collect_by_group.m.

clear all
close all
clc

%% Settings

% project folders, worked out from where the code sits, so the tree can be moved
% or copied to another drive as is (SEP_DATA_ROOT points the data elsewhere)
paths = get_paths();

% groups to collect ('rws', 'naive', 'behavior', 'young'); the mice of each come
% from the cohort registry, get_cohort
mousetypes_list = {'young'};

% ages to keep in each group ([] = whole group, or e.g. [20], [16 20 22]); P20 alone,
% since the P32 and P36 brains would dilute the youngest against the adults
age_filter = [20];

% skip, with a note, the mice whose registration has not finished, so the stack
% holds whatever is ready
skip_missing = true;

% preprocessing correction ('slicewise'); not passed to the code below, so it has
% no effect on this step
correction_type = 'slicewise';

%% Run

% pass the settings to the code, under the same names
run_settings = struct();
run_settings.paths = paths;
run_settings.mousetypes_list = mousetypes_list;
run_settings.age_filter = age_filter;
run_settings.skip_missing = skip_missing;
collect_by_group(run_settings);
