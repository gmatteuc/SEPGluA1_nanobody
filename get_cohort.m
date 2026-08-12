function cohort = get_cohort(varargin)
% GET_COHORT Single source of truth for the mouse cohort registry.
%
%   cohort = get_cohort()                          returns every registered mouse
%   cohort = get_cohort('groups', {'naive','rws'}) returns only those groups
%   cohort = get_cohort('names',  {'CGF027_Gria1'}) returns only those mice
%   get_cohort('verify')                           runs the legacy-order self-check
%
%   Returned struct array fields:
%     name          mouse folder name under <base_root>\<group>\
%     group         'rws' | 'naive' | 'behavior' | 'young'
%     age_days      postnatal age in days (NaN where not recorded)
%     share_subdir  subfolder under the mouse dir ON THE LAB SHARE holding the
%                   raw .czi ('' when they sit at the mouse root). Only used by
%                   the ingest/copy step; local copies always put .czi at the
%                   mouse root so every downstream script sees one layout.
%     base_dir      full local path: <base_root>\<group>\<name>
%
%   IMPORTANT - the adult entries (rws, naive, behavior) are listed in the exact
%   legacy order used by the hardcoded `mice`/`mousetypes` literals that P1-P7bis
%   carried before this registry existed. Several scripts select mice by numeric
%   index (e.g. P2bis `mice_to_process = 1:17`), so reordering these entries would
%   silently change which mouse is processed. get_cohort('verify') asserts the
%   order still matches the frozen literal.

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
        error('get_cohort: unknown option "%s" (use ''groups'', ''names'' or ''verify'').', string(arg));
    end
end

%% Registry

base_root = 'D:\sep_histology\data';

% name, group, age_days, share_subdir
% --- adult cohorts: LEGACY ORDER, DO NOT REORDER (see header) ---
reg = { ...
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
    % --- developmental cohort (Sami's good-quality list, 2026-08-03) ---
    'MG895_SepGluA_P36',  'young',      36, fullfile('Anatomy','Axioscan')
    'MG896_SepGluA_P28',  'young',      28, fullfile('Anatomy','Axioscan')
    'MG897_SepGluA_P20',  'young',      20, fullfile('Anatomy','Axioscan')
    'MG903_SepGluA_P20',  'young',      20, ''
    'MG904_SepGluA_P22',  'young',      22, ''
    'MG906_SepGluA_P32',  'young',      32, ''
    'MG907_SepGluA_P36',  'young',      36, ''
    'MG908_SepGluA_P32',  'young',      32, ''
    % --- added to the good-quality list 2026-08-12. MG911 at P16 extends the
    %     range below the P20 floor the cohort had until now.
    'MG909_SepGluA_P20',  'young',      20, ''
    'MG910_SepGluA_P20',  'young',      20, ''
    'MG911_SepGluA_P16',  'young',      16, ''
    'MG912_SepGluA_P20',  'young',      20, ''
    'MG913_SepGluA_P20',  'young',      20, ''
    'MG914_SepGluA_P28',  'young',      28, ''
    };

%% Build struct array

cohort = struct('name', {}, 'group', {}, 'age_days', {}, 'share_subdir', {}, 'base_dir', {});
for i = 1:size(reg, 1)
    cohort(i).name         = reg{i, 1}; %#ok<AGROW>
    cohort(i).group        = reg{i, 2}; %#ok<AGROW>
    cohort(i).age_days     = reg{i, 3}; %#ok<AGROW>
    cohort(i).share_subdir = reg{i, 4}; %#ok<AGROW>
    cohort(i).base_dir     = fullfile(base_root, reg{i, 2}, reg{i, 1}); %#ok<AGROW>
end

%% Legacy-order self-check

if do_verify
    verify_legacy_order(cohort);
end

%% Apply filters (order within the registry is always preserved)

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

%% Local function: assert the adult entries still match the frozen literal

function verify_legacy_order(cohort)
% Compares the adult portion of the registry against the exact literal that
% P1-P7bis hardcoded before the registry existed. Index-based mouse selection
% depends on this order.

legacy_mice = {'MG691_Gria1', 'MG692_Gria1', 'MG693_Gria1', 'MG736_Gria1', 'MG737_Gria1', ...
    'CGF027_Gria1', 'CGF028_Gria1', 'CGF033_Gria1', 'CGF034_Gria1', 'CGF035_Gria1', ...
    'MG705_Gria1', 'MG706_Gria1', 'MG709_Gria1', 'MG716_Gria1', 'MG718_Gria1', ...
    'MG725_Gria1', 'MG727_Gria1'};
legacy_types = {'rws','rws','rws','rws','rws', ...
    'naive','naive','naive','naive','naive', ...
    'behavior','behavior','behavior','behavior','behavior','behavior','behavior'};

n = numel(legacy_mice);
assert(numel(cohort) >= n, 'get_cohort: registry has fewer than the %d legacy mice.', n);

got_mice  = {cohort(1:n).name};
got_types = {cohort(1:n).group};

assert(isequal(got_mice, legacy_mice), ...
    'get_cohort: adult mouse order changed. Index-based selection (e.g. P2bis mice_to_process = 1:17) would break.');
assert(isequal(got_types, legacy_types), ...
    'get_cohort: adult group assignment changed relative to the legacy literal.');

fprintf('get_cohort: legacy order OK (%d adults frozen, %d mice total).\n', n, numel(cohort));
end
