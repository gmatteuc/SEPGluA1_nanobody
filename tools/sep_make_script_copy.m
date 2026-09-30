function out = sep_make_script_copy(src, out, assign, rewrites)
%SEP_MAKE_SCRIPT_COPY  Copy of a driver script with some settings changed.
%   SEP_MAKE_SCRIPT_COPY(src, out, assign, rewrites) writes a copy of the
%   script src to out, with:
%     assign    struct: for each field, the assignment to that variable is
%               replaced by field = value (the value as MATLAB text), e.g.
%               struct('mice_to_process', '{''MG914''}')
%     rewrites  {pattern, replacement; ...}: regular expressions for any other
%               line, e.g.
%               {'produce_videos\s*=\s*true;', 'produce_videos = false;'}
%   Every assignment and pattern must match, or it errors: a silent miss would
%   run the driver with its real settings.
%
%   Used to run a driver on a copy of the data tree without touching the
%   driver itself. Field names cannot contain a dot, so settings like
%   opts.smooth go through rewrites. The copy is read and written as UTF-8,
%   like the drivers, so a micro sign in a figure label survives the copy.
%
%   See also SEP_RUN_DRIVER_COPY.

if nargin < 3 || isempty(assign)
    assign = struct();
end
if nargin < 4
    rewrites = cell(0, 2);
end
txt = fileread(src, 'Encoding', 'UTF-8');

names = fieldnames(assign);
for k = 1:numel(names)
    name = names{k};

    % a multi-line cell assignment, or a one-line assignment
    pattern_cell = ['^[ \t]*' name '\s*=\s*\{.*?\};'];
    pattern_line = ['^[ \t]*' name '\s*=\s*[^;\n]*;'];
    if ~isempty(regexp(txt, pattern_cell, 'once', 'dotall', 'lineanchors'))
        pattern = pattern_cell;
        options = {'once', 'dotall', 'lineanchors'};
    else
        pattern = pattern_line;
        options = {'once', 'lineanchors'};
    end

    % replace it, or stop if there is none
    assert(~isempty(regexp(txt, pattern, options{:})), ...
        'sep_make_script_copy: no assignment to %s in %s', name, src);
    replacement = sprintf('%s = %s;', name, assign.(name));
    txt = regexprep(txt, pattern, regexptranslate('escape', replacement), options{:});
end

% apply the other rewrites, each of which must match
for k = 1:size(rewrites, 1)
    assert(~isempty(regexp(txt, rewrites{k,1}, 'once')), ...
        'sep_make_script_copy: rewrite "%s" not found in %s', rewrites{k,1}, src);
    txt = regexprep(txt, rewrites{k,1}, regexptranslate('escape', rewrites{k,2}));
end

% write the copy, in the source's encoding
folder = fileparts(out);
if ~exist(folder, 'dir')
    mkdir(folder);
end
fid = fopen(out, 'w', 'n', 'UTF-8');
fprintf(fid, '%s', txt);
fclose(fid);
end
