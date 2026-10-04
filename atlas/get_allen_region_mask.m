function mask_3d = get_allen_region_mask(csv_dir, atlas_vol, target_roots, brain_mask, ...
    name_filter)
%GET_ALLEN_REGION_MASK  Mask of the voxels of named regions and all their descendants.
%   mask_3d = GET_ALLEN_REGION_MASK(csv_dir, atlas_vol, target_roots)
%   returns a logical mask, the size of atlas_vol, of the voxels that belong to
%   the regions named in target_roots or to any region below them in the
%   ontology.
%
%   mask_3d = GET_ALLEN_REGION_MASK(..., brain_mask) also requires brain_mask.
%   mask_3d = GET_ALLEN_REGION_MASK(..., brain_mask, name_filter) keeps only
%   the regions whose name contains one of the terms of name_filter.
%
%   Inputs:
%     csv_dir      folder of the ontology CSVs (parcellation_term.csv and
%                   parcellation_to_parcellation_term_membership.csv)
%     atlas_vol     annotation volume in Allen parcellation_index values
%     target_roots  cell array of region names; each is matched exactly, else
%                   as the first name containing it, without case
%     brain_mask    voxels allowed in the mask (default: all)
%     name_filter   a term or a cell array of terms, matched without case
%                   anywhere in the name (default: no filter)
%
%   The regions are resolved by name through the ontology and mapped to the
%   atlas values through the membership table. A root that is not found is a
%   warning; none found is an error. An empty mask is a warning.

if nargin < 5
    name_filter = {};
end

% the name filter as a cell array
if ischar(name_filter) || (isstring(name_filter) && numel(name_filter) == 1)
    if strlength(string(name_filter)) == 0
        name_filter = {};
    else
        name_filter = {name_filter};
    end
end

% no brain mask: every voxel is allowed
if nargin < 4 || isempty(brain_mask)
    brain_mask = true(size(atlas_vol));
end

%% Load the ontology

term_file = fullfile(csv_dir, 'parcellation_term.csv');
map_file  = fullfile(csv_dir, 'parcellation_to_parcellation_term_membership.csv');

if ~exist(term_file, 'file') || ~exist(map_file, 'file')
    error('Allen Atlas CSV files not found in: %s', csv_dir);
end

terms = readtable(term_file);
mapping = readtable(map_file);

%% Find the root regions

root_ids = find_root_ids(target_roots, terms);

%% Add all their descendants

candidate_names = region_and_descendants(root_ids, terms);

%% Apply the name filter

target_names = apply_name_filter(name_filter, candidate_names);

%% Map the names to atlas values

% the rows of the membership table with these names
if ismember('parcellation_term_name', mapping.Properties.VariableNames)
    valid_map_rows = ismember(mapping.parcellation_term_name, target_names);
else
    error('Mapping file does not contain "parcellation_term_name" column.');
end

% their atlas values (parcellation_index)
target_pixel_vals = mapping.parcellation_index(valid_map_rows);

%% Make the mask

mask_3d = ismember(atlas_vol, target_pixel_vals) & brain_mask;

if sum(mask_3d(:)) == 0
    warning('Mask Gen: Resulting mask is empty. Check region names or atlas alignment.');
end

end

% ===== Local functions =====

function root_ids = find_root_ids(target_roots, terms)
% The ontology identifiers of the named roots: an exact name first, else the
% first name containing it; stops if none is found.

root_ids = {};
for i = 1:numel(target_roots)

    % an exact match first, else the first name that contains it
    idx = find(strcmpi(terms.name, target_roots{i}), 1);
    if isempty(idx)
        idx = find(contains(lower(terms.name), lower(target_roots{i})), 1);
    end

    if ~isempty(idx)
        root_ids = [root_ids; terms.identifier(idx)]; %#ok<AGROW>
        fprintf('Mask Gen: Found root region "%s"\n', terms.name{idx});
    else
        warning('Mask Gen: Root region "%s" not found in ontology.', target_roots{i});
    end
end

if isempty(root_ids)
    error('No valid root regions found. Mask cannot be generated.');
end

end

function candidate_names = region_and_descendants(root_ids, terms)
% The names of the roots and of every region below them.

% walk down the hierarchy one level at a time
final_term_ids = root_ids;
current_parents = root_ids;
while ~isempty(current_parents)
    is_child = ismember(terms.parent_identifier, current_parents);
    new_children = terms.identifier(is_child);

    % only children not seen yet, so a cycle cannot loop forever
    new_children = setdiff(new_children, final_term_ids);

    if isempty(new_children)
        break;
    end

    final_term_ids = [final_term_ids; new_children]; %#ok<AGROW>
    current_parents = new_children;
end

% the names of every region found, roots, inner nodes and leaves
found_rows_logical = ismember(terms.identifier, final_term_ids);
candidate_names = terms.name(found_rows_logical);

end

function target_names = apply_name_filter(name_filter, candidate_names)
% The names that contain any of the filter's terms, or all of them without one.

if ~isempty(name_filter)
    fprintf('Mask Gen: Filtering %d regions with %d filter term(s)...\n', ...
        numel(candidate_names), numel(name_filter));

    % keep a name that contains any of the terms
    keep_idx = false(size(candidate_names));
    for k = 1:numel(name_filter)
        current_term = name_filter{k};
        matches = contains(lower(candidate_names), lower(current_term));
        keep_idx = keep_idx | matches;
    end

    target_names = candidate_names(keep_idx);

    if isempty(target_names)
        warning('Mask Gen: Filters resulted in 0 regions (original selection had %d).', ...
            numel(candidate_names));
    else
        fprintf('Mask Gen: Filter retained %d regions.\n', numel(target_names));
    end
else

    % no filter: keep every region found
    target_names = candidate_names;
    fprintf('Mask Gen: Total ontology terms found (Root + Descendants): %d\n', ...
        numel(target_names));
end

end
