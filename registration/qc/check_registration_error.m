%% check_registration_error
% ===== Atlas-alignment error of each aligned mouse, read off disk =====
%
% Every mouse through the 'align' mode of run_register_to_atlas has a
% regopts.mat holding errall, the 3D fit error of the sample against the atlas
% at each optimisation step. Nothing is recomputed: this prints the last value
% for every aligned mouse beside its section count, and how far it sits from a
% fit of the error against the section count.
%
% Two limits to the numbers:
%
%   errall does not compare across atlases. The adults go to the adult CCF and the
%   young brains to DeMBA P20, and the two templates disagree on the length of the
%   brain in AP by about 11% (registration_qc\ATLAS_PARAMETERS.md), so young
%   against adult measures the difference between the atlases, not the quality of
%   the registration.
%
%   Within one atlas the section count dominates it: across the 17 adults, errall
%   against the section count gives r = -0.97, so a short brain scores badly for
%   reasons that have nothing to do with the registration. The fair reading is
%   against brains on the same atlas with a similar count, which the fit gives.
%
% Run sep_setup_paths first, once per MATLAB session.

clear all
close all
clc

%% Settings

% project folders, worked out from where the code sits, so the tree can be moved
% or copied to another drive as is
paths = get_paths();

% groups to report ('rws', 'naive', 'behavior', 'young'); keep to groups on one
% atlas for the fit at the bottom to mean anything ('young' lists them too)
groups_to_report = {'rws', 'naive', 'behavior'};

% groups the fit is computed over, all on the same atlas
reference_groups = {'rws', 'naive', 'behavior'};

%% Collect the stored alignment error

get_cohort('verify');
cohort = get_cohort();

mouse_names = {};
mouse_group = {};
err_end = [];
n_slices = [];

for k = 1:numel(cohort)

    if ~ismember(cohort(k).group, groups_to_report)
        continue
    end

    % skip a mouse not aligned yet
    regopts_name = fullfile(cohort(k).base_dir, 'lightsuite', 'regopts.mat');
    if ~exist(regopts_name, 'file')
        continue
    end

    S = load(regopts_name, 'errall');
    if ~isfield(S, 'errall') || isempty(S.errall)
        continue
    end

    % the sections that went in, after the removals in run_order_slices, to tell a
    % bad fit from a short brain
    decisions_name = fullfile(cohort(k).base_dir, 'lightsuite', ...
        'volume_for_ordering_processing_decisions.txt');
    if exist(decisions_name, 'file')
        T = readtable(decisions_name);
        kept = sum(T.FlipState ~= -1);
    else
        kept = NaN;
    end

    mouse_names{end+1} = cohort(k).name; %#ok<SAGROW>
    mouse_group{end+1} = cohort(k).group; %#ok<SAGROW>
    err_end(end+1) = S.errall(end); %#ok<SAGROW>
    n_slices(end+1) = kept; %#ok<SAGROW>

end

if isempty(err_end)
    error(['No aligned mice found in %s. regopts.mat is written by the alignment ' ...
           'stage of run_register_to_atlas, so run it with run_mode = ''align'' first.'], ...
           strjoin(groups_to_report, ', '));
end

%% Print

fprintf('\n%-26s %-9s %8s %8s %10s\n', 'mouse', 'group', 'err', 'slices', 'vs fit');
fprintf('%s\n', repmat('-', 1, 66));

is_ref = ismember(mouse_group, reference_groups) & ~isnan(n_slices);

% the error a brain with this many sections usually gets, within one atlas: a
% straight line through the reference groups, and the SD of their residuals
if nnz(is_ref) >= 4
    coef = polyfit(n_slices(is_ref), err_end(is_ref), 1);
    resid_ref = err_end(is_ref) - polyval(coef, n_slices(is_ref));
    resid_sd = std(resid_ref);
else
    coef = [];
end

% each mouse's error, and its distance from the fit in residual SDs
for k = 1:numel(mouse_names)
    if ~isempty(coef) && ~isnan(n_slices(k)) && ismember(mouse_group{k}, reference_groups)
        z = (err_end(k) - polyval(coef, n_slices(k))) / resid_sd;
        z_str = sprintf('%+.2f sd', z);
    else
        z_str = '-';
    end
    fprintf('%-26s %-9s %8.2f %8g %10s\n', ...
        mouse_names{k}, mouse_group{k}, err_end(k), n_slices(k), z_str);
end

if ~isempty(coef)
    r = corr(n_slices(is_ref)', err_end(is_ref)');
    fprintf(['\nwithin %s (n = %d): err = %.3f * sections + %.2f, ' ...
             'r = %.3f, residual sd %.2f\n'], ...
             strjoin(reference_groups, '/'), nnz(is_ref), coef(1), coef(2), r, resid_sd);
    fprintf(['Section count explains most of the spread, so judge a brain by the ' ...
             '"vs fit" column\nrather than by the raw error -- and only against ' ...
             'brains on the same atlas.\n']);
end
