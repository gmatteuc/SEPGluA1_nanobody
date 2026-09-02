clear all
close all
clc

% /// QC helper: atlas-alignment error, read straight off disk ///
% Every mouse that has been through the alignment stage of P4 has a
% regopts.mat holding `errall`, the 3D fit error of the sample against the
% atlas at each optimization step. Nothing has to be recomputed to read it, so
% this prints it for every aligned mouse alongside the section count.
%
% READ THIS BEFORE USING THE NUMBERS:
%
%   errall is NOT comparable between mice registered to different atlases.
%   The adults go to the adult CCF and the young brains to DeMBA P20, and the
%   two templates disagree about how long the brain is in AP by about 11%
%   (see registration_qc\ATLAS_PARAMETERS.md). A young-versus-adult comparison
%   of these numbers therefore measures the difference between the atlases,
%   not the quality of the registration. An earlier version of this script drew
%   exactly that plot; it was wrong and has been removed.
%
%   Within one atlas the number is dominated by how many sections a brain has.
%   Across the 17 adults, errall against section count gives r = -0.97: a short
%   brain scores badly for reasons that have nothing to do with registration.
%   So the only fair reading is against other brains on the same atlas WITH a
%   similar section count, which is what the fit below provides.

%% User-defined parameters

% Where the project lives. Derived from the location of the code rather than
% written out, so the tree can be moved or copied to another drive as is.
paths = get_paths();

% Which groups to report. Keep to groups sharing one atlas if you want the fit
% at the bottom to mean anything.
groups_to_report = {'rws', 'naive', 'behavior'};   % add 'young' to list them too

% Groups the fit is computed over (must all be on the same atlas)
reference_groups = {'rws', 'naive', 'behavior'};

%% Collect the stored alignment error

get_cohort('verify');
cohort = get_cohort();

mouse_names = {};
mouse_group = {};
err_end     = [];
n_slices    = [];

for k = 1:numel(cohort)

    if ~ismember(cohort(k).group, groups_to_report)
        continue
    end

    regopts_name = fullfile(cohort(k).base_dir, 'lightsuite', 'regopts.mat');
    if ~exist(regopts_name, 'file')
        continue    % not aligned yet
    end

    S = load(regopts_name, 'errall');
    if ~isfield(S, 'errall') || isempty(S.errall)
        continue
    end

    % How many sections actually went in, after the manual removals in P1bis.
    % Needed to tell a genuinely bad fit from a merely short brain.
    decisions_name = fullfile(cohort(k).base_dir, 'lightsuite', ...
                     'volume_for_ordering_processing_decisions.txt');
    if exist(decisions_name, 'file')
        T = readtable(decisions_name);
        kept = sum(T.FlipState ~= -1);
    else
        kept = NaN;
    end

    mouse_names{end+1} = cohort(k).name;      %#ok<SAGROW>
    mouse_group{end+1} = cohort(k).group;     %#ok<SAGROW>
    err_end(end+1)     = S.errall(end);       %#ok<SAGROW>
    n_slices(end+1)    = kept;                %#ok<SAGROW>

end

if isempty(err_end)
    error(['No aligned mice found in %s. regopts.mat is written by the alignment ' ...
           'stage of P4, so run P4 with run_mode = ''align'' first.'], ...
           strjoin(groups_to_report, ', '));
end

%% Print

fprintf('\n%-26s %-9s %8s %8s %10s\n', 'mouse', 'group', 'err', 'slices', 'vs fit');
fprintf('%s\n', repmat('-', 1, 66));

is_ref = ismember(mouse_group, reference_groups) & ~isnan(n_slices);

% The within-atlas expectation: what error does a brain with this many
% sections usually get? Anything else is not a fair comparison.
if nnz(is_ref) >= 4
    coef = polyfit(n_slices(is_ref), err_end(is_ref), 1);
    resid_ref = err_end(is_ref) - polyval(coef, n_slices(is_ref));
    resid_sd = std(resid_ref);
else
    coef = [];
end

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
