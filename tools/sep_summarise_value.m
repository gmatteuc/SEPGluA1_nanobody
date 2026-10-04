function s = sep_summarise_value(v)
%SEP_SUMMARISE_VALUE  Comparable summary of any value.
%   s = SEP_SUMMARISE_VALUE(v) walks structs, cells, tables and objects and
%   replaces numeric arrays of more than 64 elements by their hash
%   (SEP_HASH_ARRAY), so two large results can be compared with isequaln and
%   the differences located. Struct fields are sorted, so field order does
%   not count as a difference.
%
%   Objects are summarised by their values, since their class name alone
%   would make any two of them look equal: datetimes and durations as
%   numbers, categoricals as labels, containers.Map by keys and values, and
%   other value objects (the registration transforms) by their public
%   properties. Other handle objects are summarised by their class only.

if isstruct(v)
    % field by field, the fields sorted, for each element
    fields = sort(fieldnames(v))';
    s = struct('fields', {fields}, 'size', size(v), 'items', {{}});
    for i = 1:numel(v)
        item = struct();
        for k = 1:numel(fields)
            item.(fields{k}) = sep_summarise_value(v(i).(fields{k}));
        end
        s.items{end+1} = item;
    end

elseif iscell(v)
    s = cellfun(@sep_summarise_value, v, 'UniformOutput', false);

elseif isstring(v)
    s = cellstr(v);

elseif ischar(v)
    s = v;

elseif isnumeric(v) || islogical(v)
    % small arrays are kept as they are, so a difference is readable
    if numel(v) <= 64
        s = v;
    else
        s = sep_hash_array(v);
    end

elseif isa(v, 'function_handle')
    s = func2str(v);

elseif istable(v)
    s = sep_summarise_value(table2struct(v));

elseif isdatetime(v) || isduration(v)
    % exact values (datetimes as seconds since 1970) and the display format
    if isdatetime(v)
        numbers = posixtime(v);
    else
        numbers = seconds(v);
    end
    s = struct('class', class(v), 'format', v.Format, ...
        'values', sep_summarise_value(numbers));

elseif isa(v, 'categorical')
    % categorical arrays have no public properties: labels and categories
    s = struct('class', class(v), 'categories', {categories(v)}, ...
        'values', {sep_summarise_value(cellstr(v))});

elseif isa(v, 'containers.Map')
    s = struct('class', class(v), 'keys', {keys(v)}, ...
        'values', {sep_summarise_value(values(v))});

elseif isobject(v) && ~isa(v, 'handle')
    s = summarise_object(v);

else
    s = class(v);
end
end

% ===== Local functions =====

function s = summarise_object(v)
% Class and public properties of each element of a value object array; a class with
% none is read through struct(), which shows the private ones too, so its values count.

names = properties(v);
items = cell(1, numel(v));
for i = 1:numel(v)
    if isempty(names)
        % struct() on an object warns that it shows private properties: on purpose
        warning('off', 'MATLAB:structOnObject');
        items{i} = sep_summarise_value(struct(v(i)));
        warning('on', 'MATLAB:structOnObject');
        continue
    end
    item = struct();
    for k = 1:numel(names)
        item.(names{k}) = sep_summarise_value(v(i).(names{k}));
    end
    items{i} = item;
end
s = struct('class', class(v), 'size', size(v), 'items', {items});
end
