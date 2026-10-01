clear all
close all
clc

% /// Plasticity comparison, step 1 of 3: collect the registered volumes by group ///
% For each group, stacks the registered nano volumes of its mice into one 4D
% array (AP x DV x ML x mouse) and saves it in the group's folder, where
% run_normalise_groups reads it:
%   (1) takes the group's mice from the cohort registry (get_cohort), in
%       registry order, narrowed to some ages if age_filter is set
%   (2) skips, with a note, the mice whose registration has not finished
%   (3) loads each mouse's registered nano, auto and mask volumes
%       (lightsuite\volume_registered\) and stacks them
%   (4) saves <group>\nano_4d<tag>.mat, and collected_mice<tag>.mat with the
%       order of the mice along the fourth dimension (<tag> is empty for a
%       whole group, _P20 for age_filter = [20])
% Only the nano stack is saved; the other saves are commented out (see the
% note above them, in pipeline\collect_by_group.m).
% Run sep_setup_paths first, once per MATLAB session. The settings are below,
% the code is in pipeline\collect_by_group.m.

%% User-defined parameters

% Where the project lives. Derived from the location of the code rather than
% written out, so the tree can be moved or copied to another drive as is.
paths = get_paths();

% Cohort selection (mice come from the shared registry get_cohort.m).
mousetypes_list = {'young'};                    % 'rws' | 'naive' | 'behavior' | 'young'

% Optional age filter, applied within each group above. Leave empty to take the
% whole group. The young cohort spans P16-P36, but the comparison Sami wants
% first is the youngest ages against the adults, so a P20-only aggregate is
% assembled separately rather than diluting it with the P32/P36 brains.
age_filter = [20];                              % [] = whole group, e.g. [20] or [16 20 22]

% Only mice that actually reached the end of run_register_to_atlas can be
% collected here; the rest are skipped with a warning rather than killing the
% run, so the aggregate can be rebuilt as more brains finish registering.
skip_missing = true;

% Choose correction type
correction_type = 'slicewise';

%% Run

% The settings above go to the code under the same names
run_settings = struct();
run_settings.paths = paths;
run_settings.mousetypes_list = mousetypes_list;
run_settings.age_filter = age_filter;
run_settings.skip_missing = skip_missing;
collect_by_group(run_settings);
