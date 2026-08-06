function mask_3d = get_allen_region_mask(allenDir, atlas_vol, target_roots, brain_mask, name_filter)

if nargin < 5
    name_filter = {};
end

% Ensure name_filter is a cell array for consistent processing
if ischar(name_filter) || (isstring(name_filter) && numel(name_filter) == 1)
    if strlength(string(name_filter)) == 0
        name_filter = {};
    else
        name_filter = {name_filter};
    end
end

if nargin < 4 || isempty(brain_mask)
    brain_mask = true(size(atlas_vol)); % If no brain mask provided, assume all valid
end

% Load the CSVs
termFile = fullfile(allenDir, 'parcellation_term.csv');
mapFile  = fullfile(allenDir, 'parcellation_to_parcellation_term_membership.csv');

if ~exist(termFile, 'file') || ~exist(mapFile, 'file')
    error('Allen Atlas CSV files not found in: %s', allenDir);
end

terms = readtable(termFile);
mapping = readtable(mapFile);

% --- Step 1: Find Root IDs ---
root_ids = {};
for i = 1:numel(target_roots)
    % Find exact match first
    idx = find(strcmpi(terms.name, target_roots{i}), 1);

    % If not found, try partial match
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

% --- Step 2: Iteratively find ALL descendants (Hierarchy Traversal) ---
final_term_ids = root_ids;
current_parents = root_ids;

while ~isempty(current_parents)
    is_child = ismember(terms.parent_identifier, current_parents);
    new_children = terms.identifier(is_child);

    % Avoid infinite loops or duplicates
    new_children = setdiff(new_children, final_term_ids);

    if isempty(new_children)
        break;
    end

    final_term_ids = [final_term_ids; new_children]; %#ok<AGROW>
    current_parents = new_children;
end

% Get the list of names corresponding to ALL found IDs (Tree + Leaves)
found_rows_logical = ismember(terms.identifier, final_term_ids);
candidate_names = terms.name(found_rows_logical);

% --- Step 3: Apply Name Filter (List of terms) ---
if ~isempty(name_filter)
    fprintf('Mask Gen: Filtering %d regions with %d filter term(s)...\n', numel(candidate_names), numel(name_filter));

    % Initialize a logical array for indices to keep
    keep_idx = false(size(candidate_names));

    % Loop through each filter term and accumulate matches (Logical OR)
    for k = 1:numel(name_filter)
        current_term = name_filter{k};
        matches = contains(lower(candidate_names), lower(current_term));
        keep_idx = keep_idx | matches;
    end

    target_names = candidate_names(keep_idx);

    if isempty(target_names)
        warning('Mask Gen: Filters resulted in 0 regions (original selection had %d).', numel(candidate_names));
    else
        fprintf('Mask Gen: Filter retained %d regions.\n', numel(target_names));
    end
else
    % No filter, keep everything
    target_names = candidate_names;
    fprintf('Mask Gen: Total ontology terms found (Root + Descendants): %d\n', numel(target_names));
end

% --- Step 4: Map Names to Atlas Pixel Values ---

% Find these names in the mapping file
if ismember('parcellation_term_name', mapping.Properties.VariableNames)
    valid_map_rows = ismember(mapping.parcellation_term_name, target_names);
else
    error('Mapping file does not contain "parcellation_term_name" column.');
end

% Extract the pixel values (parcellation_index)
target_pixel_vals = mapping.parcellation_index(valid_map_rows);

% --- Step 5: Create the 3D Mask ---
mask_3d = ismember(atlas_vol, target_pixel_vals) & brain_mask;

if sum(mask_3d(:)) == 0
    warning('Mask Gen: Resulting mask is empty. Check region names or atlas alignment.');
end

end