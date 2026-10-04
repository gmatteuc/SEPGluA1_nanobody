function [T, ok] = sep_check_code_identity(new_dir, ref_dir, name_map)
%SEP_CHECK_CODE_IDENTITY  Check that .m files differ only in comments and layout.
%   T = SEP_CHECK_CODE_IDENTITY(new_dir, ref_dir) compares every .m file under
%   new_dir with the file at the same relative path under ref_dir, after
%   removing comments and blank lines (the code is parsed with mtree and
%   printed back with tree2str). Use it after an edit that should change
%   comments, headers or spacing only.
%
%   T = SEP_CHECK_CODE_IDENTITY(new_dir, ref_dir, name_map) first looks each
%   file up in an old-to-new name map, so a file that was moved or renamed is
%   compared with its old version. name_map is a two-column cell array
%   {old_path, new_path; ...} or a CSV file whose first line is
%   old_path,new_path, with paths relative to ref_dir and new_dir. An empty
%   old_path lists a file added on purpose, an empty new_path one removed on
%   purpose. Rows about other files than .m are ignored, so one table can
%   serve the Python check too. A map entry naming a .m file that does not
%   exist is an error.
%
%   Every file must be accounted for: a new file with no reference file, and
%   a reference file that no new file was compared with, are failures, so a
%   moved file cannot pass uncompared. Two folders with no .m file at all are
%   an error, since an empty comparison would pass without checking anything.
%
%   Skipped on both sides: hidden folders (git's folder, worktrees, the
%   environments), __pycache__, and venv*.
%
%   T is a table with one row per file: file, ref_file, status. Statuses:
%     'same code'          comments and layout only
%     'CODE CHANGED'       failure
%     'NO COUNTERPART'     new file with no reference file; failure
%     'ONLY IN REF'        reference file no new file was compared with; failure
%     'SYNTAX ERROR'       either file does not parse; failure
%     'added (listed)'     no reference file, added on purpose in the map
%     'removed (listed)'   no new file, removed on purpose in the map
%   ok is true when there is no failure. A summary is printed, then every
%   file that is not 'same code'.

% no name map: every file is compared at its own path
if nargin < 3
    name_map = cell(0, 2);
end
name_map = read_name_map(name_map);

%% List the files

% relative paths with forward slashes
new_files = list_m_files(new_dir);
ref_files = list_m_files(ref_dir);

% an empty comparison would pass without checking anything
if isempty(new_files) && isempty(ref_files)
    error('sep_check_code_identity: no .m files to compare in\n  %s\n  %s', ...
        new_dir, ref_dir);
end

%% Check the name map

% every map entry must name a file that exists: a typo would otherwise leave
% the file it meant uncompared
for i = 1:size(name_map, 1)
    if ~isempty(name_map{i, 1}) && ~ismember(name_map{i, 1}, ref_files)
        error('sep_check_code_identity: the name map lists %s, which is not in %s', ...
            name_map{i, 1}, ref_dir);
    end
    if ~isempty(name_map{i, 2}) && ~ismember(name_map{i, 2}, new_files)
        error('sep_check_code_identity: the name map lists %s, which is not in %s', ...
            name_map{i, 2}, new_dir);
    end
end

% each new file can come from one old file only
listed_new = name_map(~cellfun(@isempty, name_map(:, 2)), 2);
if numel(unique(listed_new)) < numel(listed_new)
    error('sep_check_code_identity: a new path appears twice in the name map');
end

%% Compare each file

% one row per new file
n = numel(new_files);
file = new_files(:);
ref_file = cell(n, 1);
status = cell(n, 1);

for k = 1:n

    % the reference file: from the map, or at the same path
    row = find(strcmp(name_map(:, 2), file{k}), 1);
    if ~isempty(row)
        ref_file{k} = name_map{row, 1};
    elseif ismember(file{k}, ref_files)
        ref_file{k} = file{k};
    else
        ref_file{k} = '';
    end

    % no reference file: added on purpose, or a failure
    if isempty(ref_file{k}) && ~isempty(row)
        status{k} = 'added (listed)';
        continue
    end
    if isempty(ref_file{k})
        status{k} = 'NO COUNTERPART';
        continue
    end

    % code without comments and whitespace
    [new_code, new_ok] = code_of(fullfile(new_dir, file{k}));
    [ref_code, ref_ok] = code_of(fullfile(ref_dir, ref_file{k}));
    if ~new_ok || ~ref_ok
        status{k} = 'SYNTAX ERROR';
        continue
    end

    % compare
    if strcmp(new_code, ref_code)
        status{k} = 'same code';
    else
        status{k} = 'CODE CHANGED';
    end
end

%% Reference files left over

% reference files no new file was compared with: removed on purpose, or a failure
listed_removed = name_map(cellfun(@isempty, name_map(:, 2)), 1);
unused = setdiff(ref_files, ref_file);
for i = 1:numel(unused)
    file{end+1, 1} = ''; %#ok<AGROW>
    ref_file{end+1, 1} = unused{i}; %#ok<AGROW>
    if ismember(unused{i}, listed_removed)
        status{end+1, 1} = 'removed (listed)'; %#ok<AGROW>
    else
        status{end+1, 1} = 'ONLY IN REF'; %#ok<AGROW>
    end
end

%% Summary

% the table, and ok when no file failed
T = table(file, ref_file, status);
failures = {'CODE CHANGED', 'NO COUNTERPART', 'ONLY IN REF', 'SYNTAX ERROR'};
ok = ~any(ismember(status, failures));

% the count of each status
fprintf( ...
    ['%d new and %d reference files: %d same code, %d changed, %d with no ' ...
     'counterpart, %d only in ref, %d with syntax errors, %d added and %d removed ' ...
     'as listed\n'], numel(new_files), numel(ref_files), ...
    sum(strcmp(status, 'same code')), sum(strcmp(status, 'CODE CHANGED')), ...
    sum(strcmp(status, 'NO COUNTERPART')), sum(strcmp(status, 'ONLY IN REF')), ...
    sum(strcmp(status, 'SYNTAX ERROR')), sum(strcmp(status, 'added (listed)')), ...
    sum(strcmp(status, 'removed (listed)')));

% list the files that are not the same code
flagged = ~strcmp(status, 'same code');
if any(flagged)
    disp(T(flagged, :));
end
end

% ===== Local functions =====

function map = read_name_map(map)
% The old-to-new name map as a two-column cell array of relative paths with
% forward slashes, from a cell array or from a CSV file.

if ischar(map) || isstring(map)
    csv_file = char(map);
    lines = splitlines(fileread(csv_file));

    % a file saved by a spreadsheet can start with a byte order mark
    lines{1} = strrep(lines{1}, char(65279), '');
    if ~strcmp(strtrim(lines{1}), 'old_path,new_path')
        error(['sep_check_code_identity: the first line of %s must be ' ...
               'old_path,new_path'], csv_file);
    end

    % one old,new pair per line; blank lines are skipped
    map = cell(0, 2);
    for i = 2:numel(lines)
        line = strtrim(lines{i});
        if isempty(line)
            continue
        end
        parts = strsplit(line, ',', 'CollapseDelimiters', false);
        if numel(parts) ~= 2
            error(['sep_check_code_identity: line %d of %s is not ' ...
                   'old_path,new_path: %s'], i, csv_file, line);
        end
        map(end+1, :) = strtrim(parts); %#ok<AGROW>
    end
end

% an empty map of any shape becomes 0 x 2; anything else must have two columns
if isempty(map)
    map = cell(0, 2);
end
if ~iscell(map) || size(map, 2) ~= 2
    error('sep_check_code_identity: the name map must be {old_path, new_path; ...}');
end

% forward slashes, as in the file lists
map = strrep(map, '\', '/');

% one table can serve every language: keep the rows about .m files
is_m = endsWith(map(:, 1), '.m') | endsWith(map(:, 2), '.m');
map = map(is_m, :);
end

function [code, ok] = code_of(path)
% Code of a file without comments and layout; ok is false if it does not
% parse (mtree then returns a single error node).

% parse, then print the tree back as code, without comments
tree = mtree(fileread(path));
ok = ~(tree.count == 1 && strcmp(tree.root.kind, 'ERR'));
code = '';
if ok
    code = tree2str(tree);
end
end

function files = list_m_files(folder)
% The .m files under folder, as sorted relative paths with forward slashes,
% without descending into the skipped folders.

if ~isfolder(folder)
    error('sep_check_code_identity: folder not found: %s', folder);
end

% walk the tree, keeping a list of the folders still to read
files = {};
pending = {''};
while ~isempty(pending)
    rel = pending{end};
    pending(end) = [];
    entries = dir(fullfile(folder, rel));
    for k = 1:numel(entries)
        name = entries(k).name;
        if entries(k).isdir && ~is_skipped(name)
            pending{end+1} = join_relative(rel, name); %#ok<AGROW>
        elseif ~entries(k).isdir && endsWith(name, '.m')
            files{end+1} = join_relative(rel, name); %#ok<AGROW>
        end
    end
end
files = sort(files)';
end

function skip = is_skipped(name)
% Folders that hold no code of ours: hidden ones (. and .., git's folder, worktrees,
% environments), Python caches and environments (venv*).

skip = startsWith(name, '.') || strcmp(name, '__pycache__') || startsWith(name, 'venv');
end

function p = join_relative(folder, name)
% A relative path with forward slashes.

if isempty(folder)
    p = name;
else
    p = [folder '/' name];
end
end
