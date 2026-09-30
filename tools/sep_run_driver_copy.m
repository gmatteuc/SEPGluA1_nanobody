function sep_run_driver_copy(driver, data_root, assign, rewrites)
%SEP_RUN_DRIVER_COPY  Run a driver script on a copy of the data tree.
%   SEP_RUN_DRIVER_COPY(driver, data_root, assign, rewrites) makes a copy of
%   the driver script (by name, e.g. 'P5_collect_data_by_group') with the
%   settings changed as in sep_make_script_copy, points the whole data tree,
%   inputs and outputs, to data_root through the environment variable
%   SEP_DATA_ROOT, which get_paths honours, and runs it. The variable's
%   previous value is restored when the run ends, also after an error.
%
%   Refused: the production data (D:\sep_histology\data) or any folder
%   inside it, and any folder inside the snapshot on G:. The copy runs only
%   if get_paths then resolves data_root, so code that ignores the variable
%   cannot write into the production tree.
%
%   Run it in a fresh MATLAB session, through tools\run_matlab_detached.ps1 or
%   from the command line (one line):
%     matlab -batch "restoredefaultpath; cd('G:\sep_refactor\check\code');
%       sep_setup_paths; addpath('tools');
%       sep_run_driver_copy('P5_collect_data_by_group', 'G:\sep_refactor\check\data')"
%   A session that already ran other code can give different figures (see
%   tools\README.md).
%
%   The copy is written next to the data tree, as
%   <data_root>\..\driver_copies\<driver>\<driver>_copy.m, and kept as a
%   record of what ran. Its folder goes on the path for the run only, so the
%   current folder stays the one the driver would run from.
%
%   See also SEP_MAKE_SCRIPT_COPY, SEP_COMPARE_OUTPUTS.

if nargin < 3 || isempty(assign)
    assign = struct();
end
if nargin < 4
    rewrites = cell(0, 2);
end

% a check never writes into the production data or into the snapshot; both
% are refused from the path alone, before anything there is touched
data_root = canonical(absolute_path(data_root));
production_data = 'D:\sep_histology\data';
if strcmpi(data_root, production_data) || ...
        startsWith(lower(data_root), lower([production_data '\']))
    error(['sep_run_driver_copy: %s is the production data. A check runs on a copy ' ...
        '(docs/REFACTOR_PLAN.md, Verification design).'], data_root);
end
if startsWith(lower(data_root), lower('G:\sep_histology_snapshot'))
    error('sep_run_driver_copy: %s is inside the snapshot on G:, which is a backup.', ...
        data_root);
end
if ~isfolder(data_root)
    error('sep_run_driver_copy: data root not found: %s', data_root);
end

src = which([driver '.m']);
if isempty(src)
    error(['sep_run_driver_copy: driver %s not found on the path. Run ' ...
        'sep_setup_paths first.'], driver);
end

% copy of the driver, next to the data tree
copy_folder = fullfile(fileparts(data_root), 'driver_copies', driver);
copy_name = [driver '_copy'];
copy_file = fullfile(copy_folder, [copy_name '.m']);
sep_make_script_copy(src, copy_file, assign, rewrites);

% redirect the data tree, and put the variable back however the run ends
previous_root = getenv('SEP_DATA_ROOT');
setenv('SEP_DATA_ROOT', data_root);
restore_root = onCleanup(@() setenv('SEP_DATA_ROOT', previous_root));

% the code must honour the redirect before anything runs
p = get_paths();
if ~strcmpi(canonical(p.data), data_root)
    error(['sep_run_driver_copy: get_paths (%s) resolves the data root to %s, ' ...
        'not %s. Nothing was run.'], which('get_paths'), p.data, data_root);
end
fprintf(['running %s\n  driver     %s\n  copy       %s\n  data root  %s\n' ...
    '  get_paths  %s\n'], driver, src, copy_file, p.data, which('get_paths'));

% the copy's folder goes first on the path, for this run only
addpath(copy_folder);
restore_path = onCleanup(@() rmpath(copy_folder));
run_script(copy_name);
end

% ===== Local functions =====

function run_script(name)
% Run a script in this function's own workspace. Drivers start with clear
% all, which clears the workspace they run in: here that is this one, not
% the caller's, where the cleanup that restores SEP_DATA_ROOT lives.

eval([name ';']);
end

function p = absolute_path(p)
% Absolute path, worked out from the text alone: a relative path is taken
% from the current folder, and . and .. are resolved.

p = char(p);
if isempty(regexp(p, '^([A-Za-z]:[\\/]|[\\/][\\/])', 'once'))
    p = fullfile(pwd, p);
end
p = char(java.io.File(p).toPath().normalize().toString());
end

function s = canonical(s)
% Backslashes and no trailing separator, the spelling get_paths compares.

s = strrep(char(s), '/', '\');
while numel(s) > 3 && s(end) == '\'
    s(end) = [];
end
end
