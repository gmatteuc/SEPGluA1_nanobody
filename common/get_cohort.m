function cohort = get_cohort(varargin)
%GET_COHORT  The mouse cohort registry, the one list of the project's mice.
%   cohort = GET_COHORT() returns every registered mouse.
%   cohort = GET_COHORT('groups', {'naive', 'rws'}) returns only those groups.
%   cohort = GET_COHORT('names', {'CGF027_Gria1'}) returns only those mice.
%   GET_COHORT('verify') runs the legacy-order self-check.
%   Options combine; the registry's order is kept whatever the filter.
%
%   cohort is a struct array with the fields:
%     name          mouse folder name under <data>\<group>\
%     group         'rws', 'naive', 'behavior' or 'young'
%     age_days      postnatal age in days (NaN where not recorded)
%     share_subdir  subfolder of the mouse folder on the lab share that holds
%                   the raw .czi ('' when they sit at the mouse folder's root).
%                   Only the copy step uses it: the local copies always put
%                   the .czi at the mouse folder's root, so every later script
%                   sees one layout.
%     base_dir      full local path, <data>\<group>\<name>
%
%   The adult entries (rws, naive, behavior) are listed in the exact order of
%   the mouse lists the drivers from run_extract_and_center to
%   run_group_differences wrote out before this registry. Several scripts
%   select mice by index (run_nano_equalisation, mice_to_process = 1:17), so
%   reordering these entries would silently change which mouse is processed.
%   GET_COHORT('verify') checks that the order still matches that list.

%% Parse inputs

group_filter = {};
name_filter  = {};
do_verify    = false;

k = 1;
while k <= numel(varargin)
    arg = varargin{k};
    if strcmpi(arg, 'verify')
        do_verify = true;
        k = k + 1;
    elseif strcmpi(arg, 'groups')
        group_filter = cellstr(varargin{k+1});
        k = k + 2;
    elseif strcmpi(arg, 'names')
        name_filter = cellstr(varargin{k+1});
        k = k + 2;
    else
        error('get_cohort: unknown option "%s" (use ''groups'', ''names'' or ''verify'').', ...
            string(arg));
    end
end

%% Registry

paths = get_paths();
base_root = paths.data;

% name, group, age_days, share_subdir; the adults in their legacy order, never
% to be reordered (see the help)
registry_rows = { ...
    'MG691_Gria1',        'rws',       NaN, ''
    'MG692_Gria1',        'rws',       NaN, ''
    'MG693_Gria1',        'rws',       NaN, ''
    'MG736_Gria1',        'rws',       NaN, ''
    'MG737_Gria1',        'rws',       NaN, ''
    'CGF027_Gria1',       'naive',     NaN, ''
    'CGF028_Gria1',       'naive',     NaN, ''
    'CGF033_Gria1',       'naive',     NaN, ''
    'CGF034_Gria1',       'naive',     NaN, ''
    'CGF035_Gria1',       'naive',     NaN, ''
    'MG705_Gria1',        'behavior',  NaN, ''
    'MG706_Gria1',        'behavior',  NaN, ''
    'MG709_Gria1',        'behavior',  NaN, ''
    'MG716_Gria1',        'behavior',  NaN, ''
    'MG718_Gria1',        'behavior',  NaN, ''
    'MG725_Gria1',        'behavior',  NaN, ''
    'MG727_Gria1',        'behavior',  NaN, ''
    % the young cohort: Sami's list of good-quality brains, 3 Aug 2026
    'MG895_SepGluA_P36',  'young',      36, fullfile('Anatomy','Axioscan')
    'MG896_SepGluA_P28',  'young',      28, fullfile('Anatomy','Axioscan')
    'MG897_SepGluA_P20',  'young',      20, fullfile('Anatomy','Axioscan')
    'MG903_SepGluA_P20',  'young',      20, ''
    'MG904_SepGluA_P22',  'young',      22, ''
    'MG906_SepGluA_P32',  'young',      32, ''
    'MG907_SepGluA_P36',  'young',      36, ''
    'MG908_SepGluA_P32',  'young',      32, ''
    % added to that list on 12 Aug 2026; MG911, at P16, is the youngest brain
    'MG909_SepGluA_P20',  'young',      20, ''
    'MG910_SepGluA_P20',  'young',      20, ''
    'MG911_SepGluA_P16',  'young',      16, ''
    'MG912_SepGluA_P20',  'young',      20, ''
    'MG913_SepGluA_P20',  'young',      20, ''
    'MG914_SepGluA_P28',  'young',      28, ''
    };

%% Build the struct array

cohort = struct('name', {}, 'group', {}, 'age_days', {}, 'share_subdir', {}, ...
    'base_dir', {});
for i = 1:size(registry_rows, 1)
    cohort(i).name         = registry_rows{i, 1}; %#ok<AGROW>
    cohort(i).group        = registry_rows{i, 2}; %#ok<AGROW>
    cohort(i).age_days     = registry_rows{i, 3}; %#ok<AGROW>
    cohort(i).share_subdir = registry_rows{i, 4}; %#ok<AGROW>
    cohort(i).base_dir     = fullfile(base_root, registry_rows{i, 2}, ...
        registry_rows{i, 1}); %#ok<AGROW>
end

%% Check the legacy order

if do_verify
    verify_legacy_order(cohort);
end

%% Apply the filters

% the registry's order is kept
if ~isempty(group_filter)
    keep   = ismember({cohort.group}, group_filter);
    cohort = cohort(keep);
end

if ~isempty(name_filter)
    keep   = ismember({cohort.name}, name_filter);
    missing = setdiff(name_filter, {cohort.name});
    if ~isempty(missing)
        error('get_cohort: name(s) not in registry: %s', strjoin(missing, ', '));
    end
    cohort = cohort(keep);
end

end

% ===== Local functions =====

function verify_legacy_order(cohort)
% Check the adults of the registry against the mouse list the drivers wrote out
% before the registry existed, which selection by index depends on.

legacy_mice = { ...
    'MG691_Gria1', 'MG692_Gria1', 'MG693_Gria1', 'MG736_Gria1', 'MG737_Gria1', ...
    'CGF027_Gria1', 'CGF028_Gria1', 'CGF033_Gria1', 'CGF034_Gria1', 'CGF035_Gria1', ...
    'MG705_Gria1', 'MG706_Gria1', 'MG709_Gria1', 'MG716_Gria1', 'MG718_Gria1', ...
    'MG725_Gria1', 'MG727_Gria1'};
legacy_types = {'rws', 'rws', 'rws', 'rws', 'rws', ...
    'naive', 'naive', 'naive', 'naive', 'naive', ...
    'behavior', 'behavior', 'behavior', 'behavior', 'behavior', 'behavior', 'behavior'};

n = numel(legacy_mice);
assert(numel(cohort) >= n, 'get_cohort: registry has fewer than the %d legacy mice.', n);

got_mice  = {cohort(1:n).name};
got_types = {cohort(1:n).group};

assert(isequal(got_mice, legacy_mice), ...
    'get_cohort: adult mouse order changed. Index-based selection (e.g. run_nano_equalisation mice_to_process = 1:17) would break.');
assert(isequal(got_types, legacy_types), ...
    'get_cohort: adult group assignment changed relative to the legacy literal.');

fprintf('get_cohort: legacy order OK (%d adults frozen, %d mice total).\n', n, ...
    numel(cohort));

end
