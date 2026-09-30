function T = sep_compare_outputs(ref_dir, new_dir, opts)
%SEP_COMPARE_OUTPUTS  Compare two output folders file by file.
%   T = SEP_COMPARE_OUTPUTS(ref_dir, new_dir) compares every file under
%   ref_dir (the reference run) with the file at the same relative path under
%   new_dir (the run after a change), and lists files present in only one of
%   them.
%
%   How each type is compared:
%     .mat        contents, via sep_summarise_value (large arrays by hash,
%                 objects such as transforms by their values)
%     .fig        the plotted data: data, limits and text of every graphics
%                 object (the file itself holds its creation date)
%     .csv .tsv   table contents, column by column
%     .png .jpg   bytes, then pixels (with the colour map of an indexed
%                 image), then sep_compare_figure_images: 'same render' if
%                 every changed pixel, channel by channel, is a mix of the
%                 original's colours within 1 px (anti-aliasing), else
%                 'DIFFERENT', as is an image of another size. The detail
%                 gives the number of changed pixels
%     .eps        text, without the creation date and title lines
%     .tif .tiff  pixels, every page of a stack
%     .avi .mp4   per-frame hashes
%     .m .txt .md .json
%                 text
%     other       bytes
%   The Python route's arrays (.npy, .npz) are compared by bytes here, which
%   says whether they differ but not where: tools/compare_outputs.py compares
%   them array by array.
%
%   T is a table with one row per file: file, result, detail. Results are
%   'same', 'same render', 'DIFFERENT', 'only in ref', 'only in new',
%   'NOT REWRITTEN', 'compare failed'. A summary is printed, then every file
%   that is not 'same', 'same render' included.
%
%   Options (opts):
%     ignore_files   regexps of relative paths to skip
%                    (default {'\.log$', 'log\.txt$'})
%     ignore_fields  field names and table columns left out of .mat and
%                    .csv comparisons (run times)               (default {})
%     replace_text   {old, new; ...}: text replaced in the reference's values
%                    and text files before comparing, for the paths a run
%                    saves, which name its own tree, e.g.
%                    {'G:\sep_refactor\ref\data', 'G:\sep_refactor\check\data'}
%                                                               (default {})
%     newer_than     datetime: a file of new_dir last written before it is
%                    'NOT REWRITTEN', so a run that wrote nothing cannot
%                    pass; compare output folders only, since inputs are
%                    never rewritten                            (default [])
%     max_detail     differences listed per file                (default 5)

if nargin < 3
    opts = struct();
end
opts = sep_struct_defaults(opts, struct('ignore_files', {{'\.log$', 'log\.txt$'}}, ...
    'ignore_fields', {{}}, 'replace_text', {cell(0, 2)}, 'newer_than', [], ...
    'max_detail', 5));

% absolute paths, so relative inputs work too
ref_dir = absolute_path(ref_dir);
new_dir = absolute_path(new_dir);

% files of both folders, as relative paths, with the new files' dates
ref_files = relative_files(ref_dir);
[new_files, new_dates] = relative_files(new_dir);

% skip ignored files
is_ignored = @(list) cellfun( ...
    @(f) any(~cellfun(@isempty, regexp(f, opts.ignore_files, 'once'))), list);
ref_files = ref_files(~is_ignored(ref_files));
kept = ~is_ignored(new_files);
new_files = new_files(kept);
new_dates = new_dates(kept);

% an empty comparison would pass without checking anything
if isempty(ref_files) && isempty(new_files)
    error('sep_compare_outputs: no files to compare in\n  %s\n  %s', ref_dir, new_dir);
end

% one row per file of either folder
all_files = union(ref_files, new_files);
all_files = all_files(:);
n = numel(all_files);
result = cell(n, 1);
detail = cell(n, 1);

for k = 1:n

    % files in only one of the folders
    f = all_files{k};
    detail{k} = '';
    if ~ismember(f, new_files)
        result{k} = 'only in ref';
        continue
    end
    if ~ismember(f, ref_files)
        result{k} = 'only in new';
        continue
    end

    % compare by type; a failure is reported, not raised
    try
        [result{k}, detail{k}] = compare_file(fullfile(ref_dir, f), ...
            fullfile(new_dir, f), opts);
    catch err
        result{k} = 'compare failed';
        detail{k} = sprintf('%s (%s, line %d)', err.message, err.stack(1).name, ...
            err.stack(1).line);
    end

    % a file the run did not write can only repeat what was there before
    if ~isempty(opts.newer_than)
        written = datetime(new_dates(strcmp(new_files, f)), 'ConvertFrom', 'datenum');
        if written < opts.newer_than
            detail{k} = sprintf('written %s, compared: %s', ...
                char(written, 'yyyy-MM-dd HH:mm:ss'), result{k});
            result{k} = 'NOT REWRITTEN';
        end
    end
end

T = table(all_files, result, detail, 'VariableNames', {'file', 'result', 'detail'});

% summary
kinds = {'same', 'same render', 'DIFFERENT', 'only in ref', 'only in new', ...
    'NOT REWRITTEN', 'compare failed'};
fprintf('%d files:', n);
for i = 1:numel(kinds)
    fprintf('  %s %d', kinds{i}, sum(strcmp(result, kinds{i})));
end
fprintf('\n');

% list the files that are not the same: 'same render' too, since its count of
% changed pixels is worth a look
flagged = ~strcmp(result, 'same');
if any(flagged)
    disp(T(flagged, :));
end
end

% ===== Local functions =====

function [result, detail] = compare_file(ref, new, opts)
% Compare one pair of files according to their type. Replacements apply to
% the reference side only.

detail = '';
no_replacements = cell(0, 2);
[~, ~, ext] = fileparts(ref);
switch lower(ext)
    case '.mat'
        % one file in memory at a time: some cohort files are several GB
        a = sep_summarise_value(clean_value(load(ref), opts.ignore_fields, ...
            opts.replace_text));
        b = sep_summarise_value(clean_value(load(new), opts.ignore_fields, ...
            no_replacements));
        d = walk_differences(a, b, '', {});
    case '.fig'
        d = value_differences( ...
            clean_value(figure_contents(ref), {}, opts.replace_text), ...
            figure_contents(new), 'figure', {});
    case {'.csv', '.tsv'}
        d = table_differences( ...
            clean_table(read_table(ref), opts.ignore_fields, opts.replace_text), ...
            clean_table(read_table(new), opts.ignore_fields, no_replacements));
    case {'.png', '.jpg', '.jpeg'}
        [result, detail] = compare_images(ref, new);
        return
    case '.eps'
        d = text_differences( ...
            clean_value(strip_eps_dates(fileread(ref)), {}, opts.replace_text), ...
            strip_eps_dates(fileread(new)));
    case {'.tif', '.tiff'}
        d = value_differences(tiff_hashes(ref), tiff_hashes(new), 'pages', {});
    case {'.avi', '.mp4'}
        d = value_differences(video_hashes(ref), video_hashes(new), 'frames', {});
    case {'.m', '.txt', '.md', '.json'}
        d = text_differences(clean_value(fileread(ref), {}, opts.replace_text), ...
            fileread(new));
    otherwise
        d = value_differences(sep_hash_array(file_bytes(ref)), ...
            sep_hash_array(file_bytes(new)), 'bytes', {});
end

% list the first differences
if isempty(d)
    result = 'same';
else
    result = 'DIFFERENT';
    detail = strjoin(d(1:min(end, opts.max_detail)), ' | ');
end
end

function [result, detail] = compare_images(ref, new)
% Same bytes, same pixels, or same figure up to rendering.

detail = '';
if isequal(file_bytes(ref), file_bytes(new))
    result = 'same';
    return
end

% same pixels, of the same class, and for an indexed image the same colour
% map (the pixels are then only indices into it)
[image_ref, map_ref] = imread(ref);
[image_new, map_new] = imread(new);
if isequal(image_ref, image_new) && strcmp(class(image_ref), class(image_new)) ...
        && isequal(map_ref, map_new)
    result = 'same';
    return
end

% otherwise, compare as figures
R = sep_compare_figure_images(ref, new);
if R.ok
    result = 'same render';
else
    result = 'DIFFERENT';
end
if ~strcmp(R.size_a, R.size_b)
    detail = sprintf('size %s vs %s', R.size_a, R.size_b);
else
    detail = sprintf(['%d pixels changed, %d of them not a mix of the original''s ' ...
        'colours within 1 px (largest distance from a mix %.0f of 255)'], ...
        R.n_changed, R.n_real, R.max_dist);
end
end

function d = value_differences(a, b, path, d)
% Paths at which two values differ (values summarised first).

a = sep_summarise_value(a);
b = sep_summarise_value(b);
d = walk_differences(a, b, path, d);
end

function d = walk_differences(a, b, path, d)
% Recursive part of value_differences, on summarised values.

% equal values of different classes (single against double) still differ
if isequaln(a, b) && strcmp(class(a), class(b))
    return
end

% summarised struct: field lists are reported item by item below, so only
% compare size and items here
if isstruct(a) && isstruct(b) && isscalar(a) && isscalar(b) && ...
        all(isfield(a, {'fields', 'size', 'items'})) && ...
        all(isfield(b, {'fields', 'size', 'items'}))
    d = walk_differences(a.size, b.size, [path '.size'], d);
    if isequal(a.size, [1 1]) && isequal(b.size, [1 1])
        d = walk_differences(a.items{1}, b.items{1}, path, d);
    else
        d = walk_differences(a.items, b.items, path, d);
    end
    return
end

% structs: compare field by field
if isstruct(a) && isstruct(b) && isequal(size(a), size(b)) && isscalar(a)
    fa = fieldnames(a);
    fb = fieldnames(b);

    % loop over indices: a for loop over an empty cell can still run once
    only_ref = setdiff(fa, fb);
    only_new = setdiff(fb, fa);
    common = intersect(fa, fb);
    for i = 1:numel(only_ref)
        d{end+1} = [path '.' only_ref{i} ' only in ref']; %#ok<AGROW>
    end
    for i = 1:numel(only_new)
        d{end+1} = [path '.' only_new{i} ' only in new']; %#ok<AGROW>
    end
    for i = 1:numel(common)
        d = walk_differences(a.(common{i}), b.(common{i}), ...
            [path '.' common{i}], d);
    end
    return
end

% cells: compare element by element
if iscell(a) && iscell(b) && isequal(size(a), size(b))
    for i = 1:numel(a)
        d = walk_differences(a{i}, b{i}, sprintf('%s{%d}', path, i), d);
    end
    return
end

% small arrays of the same class and size: the first element that differs
same_shape = strcmp(class(a), class(b)) && isequal(size(a), size(b));
if same_shape && (isnumeric(a) || islogical(a))
    i = find(~(a == b | (isnan(double(a)) & isnan(double(b)))), 1);
    d{end+1} = sprintf('%s(%d): %s vs %s', path, i, short_text(a(i)), short_text(b(i)));
    return
end

% anything else: report the two values, with their classes when they differ
if strcmp(class(a), class(b))
    d{end+1} = sprintf('%s: %s vs %s', path, short_text(a), short_text(b));
else
    d{end+1} = sprintf('%s: %s %s vs %s %s', path, class(a), short_text(a), ...
        class(b), short_text(b));
end
end

function s = clean_value(s, names, replacements)
% Remove the named fields at any depth of a struct, and apply the text
% replacements ({old, new; ...}) to every text value.

if isstruct(s)
    s = rmfield(s, intersect(fieldnames(s), names));
    f = fieldnames(s);
    for i = 1:numel(s)
        for k = 1:numel(f)
            s(i).(f{k}) = clean_value(s(i).(f{k}), names, replacements);
        end
    end
elseif iscell(s)
    for i = 1:numel(s)
        s{i} = clean_value(s{i}, names, replacements);
    end
elseif ischar(s) || isstring(s)
    for i = 1:size(replacements, 1)
        s = strrep(s, replacements{i, 1}, replacements{i, 2});
    end
end
end

function t = read_table(f)
% A .csv or .tsv file as a table. Column names that are not valid MATLAB
% names are changed the same way in both files, so the warning is not shown.

warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');
restore = onCleanup(@() warning('on', 'MATLAB:table:ModifiedAndSavedVarnames'));
[~, ~, ext] = fileparts(f);
if strcmpi(ext, '.tsv')
    t = readtable(f, 'FileType', 'text', 'Delimiter', '\t');
else
    t = readtable(f);
end
end

function t = clean_table(t, names, replacements)
% A table without the ignored columns, with the text replacements applied to
% its text columns.

t(:, intersect(t.Properties.VariableNames, names)) = [];
columns = t.Properties.VariableNames;
for k = 1:numel(columns)
    column = t.(columns{k});
    if iscell(column) || ischar(column) || isstring(column)
        t.(columns{k}) = clean_value(column, {}, replacements);
    end
end
end

function d = table_differences(a, b)
% Differing columns of two tables, each with its first differing row.

d = {};
names_a = a.Properties.VariableNames;
names_b = b.Properties.VariableNames;
if ~isequal(names_a, names_b)
    d{end+1} = sprintf('columns: %s vs %s', strjoin(names_a, ' '), strjoin(names_b, ' '));
    return
end
if height(a) ~= height(b)
    d{end+1} = sprintf('%d rows vs %d', height(a), height(b));
    return
end

for k = 1:numel(names_a)
    x = a.(names_a{k});
    y = b.(names_a{k});
    if isequaln(x, y) && strcmp(class(x), class(y))
        continue
    end

    % the first row that differs, with both values
    for r = 1:height(a)
        if ~isequaln(x(r, :), y(r, :)) || ~strcmp(class(x), class(y))
            d{end+1} = sprintf('%s, row %d: %s vs %s', names_a{k}, r, ...
                short_text(x(r, :)), short_text(y(r, :))); %#ok<AGROW>
            break
        end
    end
end
end

function s = figure_contents(f)
% What a saved figure shows: the data, limits and text of every graphics
% object, in the order the figure holds them. Window positions are left out,
% since they depend on the screen.

fig = openfig(f, 'invisible');
closer = onCleanup(@() delete(fig));
objects = findall(fig);
names = {'Type', 'XData', 'YData', 'ZData', 'CData', 'String', 'XLim', 'YLim', ...
    'ZLim', 'CLim', 'Colormap', 'Visible'};
s = cell(numel(objects), 1);
for k = 1:numel(objects)
    item = struct();
    for i = 1:numel(names)
        if isprop(objects(k), names{i})
            item.(names{i}) = get(objects(k), names{i});
        end
    end
    s{k} = item;
end
end

function d = text_differences(a, b)
% First differing lines of two texts.

d = {};
if strcmp(a, b)
    return
end

% compare line by line, ignoring leading and trailing blanks; list the first 3
la = splitlines(a);
lb = splitlines(b);
n_diff = 0;
for i = 1:max(numel(la), numel(lb))
    x = '';
    y = '';
    if i <= numel(la)
        x = strtrim(la{i});
    end
    if i <= numel(lb)
        y = strtrim(lb{i});
    end
    if ~strcmp(x, y)
        n_diff = n_diff + 1;
        if n_diff <= 3
            d{end+1} = sprintf('line %d: "%s" vs "%s"', i, x, y); %#ok<AGROW>
        end
    end
end

% 0 differing lines: the texts differ only in line endings or blanks at the
% ends of lines, still reported
d{end+1} = sprintf('%d differing lines', n_diff);
end

function t = strip_eps_dates(t)
% EPS text without the lines that change at every export.

t = regexprep(t, '%%(CreationDate|Title):[^\n]*', '');
end

function h = tiff_hashes(f)
% One hash per page of a TIFF file, so a stack is compared whole and not
% only its first page.

n_pages = numel(imfinfo(f));
h = cell(1, n_pages);
for k = 1:n_pages
    h{k} = sep_hash_array(imread(f, k));
end
end

function h = video_hashes(f)
% One hash per video frame.

v = VideoReader(f);
h = {};
while hasFrame(v)
    h{end+1} = sep_hash_array(readFrame(v)); %#ok<AGROW>
end
end

function b = file_bytes(f)
% Raw bytes of a file.

fid = fopen(f, 'r');
closer = onCleanup(@() fclose(fid));
b = fread(fid, Inf, '*uint8');
end

function [files, dates] = relative_files(folder)
% All files under folder, as relative paths with forward slashes, and the
% date each was last written (datenum).

d = dir(fullfile(folder, '**', '*'));
d = d(~[d.isdir]);
files = cell(numel(d), 1);
dates = zeros(numel(d), 1);
for k = 1:numel(d)
    full_path = fullfile(d(k).folder, d(k).name);
    files{k} = strrep(full_path(numel(folder) + 2:end), '\', '/');
    dates(k) = d(k).datenum;
end
end

function p = absolute_path(p)
% Absolute path of an existing folder.

if ~isfolder(p)
    error('sep_compare_outputs: folder not found: %s', p);
end
info = dir(p);
p = info(1).folder;
end

function s = short_text(v)
% Short printable version of a value.

if iscell(v) && isscalar(v)
    s = short_text(v{1});
elseif iscell(v)
    s = sprintf('cell %s', mat2str(size(v)));
elseif ischar(v)
    s = v;
elseif isnumeric(v) || islogical(v)
    s = mat2str(v(1:min(end, 6)));
else
    s = class(v);
end
if numel(s) > 80
    s = [s(1:80) '...'];
end
end
