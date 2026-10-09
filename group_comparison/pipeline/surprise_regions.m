function [T_regions, valid_pixels, region_of_voxel] = surprise_regions(AllenCrop, ...
    allenDir)
%SURPRISE_REGIONS  The regions of the surprise bars, and their voxels.
%   [T_regions, valid_pixels, region_of_voxel] = SURPRISE_REGIONS(AllenCrop,
%   allenDir) returns the table of the regions (surprise_region_acronyms), one
%   row each (acronym, name, group, label, voxels in the left hemisphere, the
%   regions of the list taken out of it), the atlas voxels of the left
%   hemisphere of AllenCrop (the width of the folded maps), and the region
%   each of them falls in (its row, 0 for none). A region of the list that
%   holds another gives up that one's voxels, so no voxel is in two regions.
%   allenDir is the folder of the ontology tables. Used by group_differences
%   (the region bars and their permutation test) and per_mouse_region_values.

% the ontology: every term with its parent, and each atlas value's term at the
% atlas's five levels (organ, category, division, structure, substructure)
term_file = fullfile(allenDir, 'parcellation_term.csv');
T_terms = readtable(term_file);
T_members = readtable(fullfile(allenDir, ...
    'parcellation_to_parcellation_term_membership.csv'));

% the isocortical areas, taken from the atlas, then the subcortical regions,
% declared
acronyms = surprise_region_acronyms(T_members);

% each region's name and atlas values, its own and those of every term below it
[names, region_values] = region_atlas_values(acronyms, T_terms, T_members, term_file);

% a region of the list that holds another gives up that one's values, so no voxel
% counts in two bars
[region_values, removed_idx] = remove_nested_regions(acronyms, region_values);

% each region's group: its division in the atlas
groups = region_divisions(region_values, T_members);

% the left hemisphere's atlas voxels, and the region of each
[valid_pixels, region_of_voxel] = label_left_hemisphere(AllenCrop, region_values);

% the table of the regions, printed
T_regions = region_table(acronyms, names, groups, removed_idx, region_of_voxel);
print_region_table(T_regions);
end

% ===== Local functions =====

function [names, region_values] = region_atlas_values(acronyms, T_terms, T_members, ...
    term_file)
% Each region's name, and its atlas values (parcellation_index): those of its term
% and of every term below it in the ontology (T_terms, read from term_file).

n_regions = numel(acronyms);
names = cell(n_regions, 1);
region_values = cell(n_regions, 1);
for r = 1:n_regions

    % the term, by its exact acronym: an acronym the atlas lacks stops the run,
    % rather than resolving to another region
    term_row = find(strcmp(T_terms.acronym, acronyms{r}), 1);
    if isempty(term_row)
        error(['surprise_regions: the surprise-bar region %s is not an ' ...
               'acronym of %s. Spell it as the atlas does.'], acronyms{r}, term_file);
    end
    names{r} = T_terms.name{term_row};

    % the term and every term below it, one level at a time
    term_ids = T_terms.identifier(term_row);
    parent_ids = term_ids;
    while ~isempty(parent_ids)
        is_child = ismember(T_terms.parent_identifier, parent_ids);
        child_ids = setdiff(unique(T_terms.identifier(is_child)), term_ids);
        term_ids = [term_ids; child_ids(:)]; %#ok<AGROW>
        parent_ids = child_ids;
    end

    % the atlas values whose term, at any level, is one of these
    term_acronyms = unique(T_terms.acronym(ismember(T_terms.identifier, term_ids)));
    is_member = ismember(T_members.parcellation_term_acronym, term_acronyms);
    region_values{r} = unique(T_members.parcellation_index(is_member));

    % a region with no atlas value would have no voxel
    if isempty(region_values{r})
        error(['surprise_regions: the surprise-bar region %s has no value in ' ...
               'the atlas annotation. Take it out of the list.'], acronyms{r});
    end
end
end

function [own_values, removed_idx] = remove_nested_regions(acronyms, region_values)
% Each region's atlas values without those of the regions of the list inside it
% (HPF without SUB), and the positions of those regions in the list. Stops if a
% region has no value left, or if two regions still share one.

n_regions = numel(acronyms);
own_values = region_values;
removed_idx = cell(n_regions, 1);
for r = 1:n_regions
    for k = 1:n_regions

        % region k inside region r: every value of k is one of r's, and r has more
        is_inside = k ~= r && all(ismember(region_values{k}, region_values{r})) && ...
            numel(region_values{k}) < numel(region_values{r});
        if is_inside
            own_values{r} = setdiff(own_values{r}, region_values{k});
            removed_idx{r} = [removed_idx{r} k];
        end
    end

    % a region the regions inside it cover entirely would have no voxel of its own
    if isempty(own_values{r})
        error(['surprise_regions: the surprise-bar region %s is covered by the ' ...
               'regions of the list inside it (%s). Take it out of the list.'], ...
              acronyms{r}, strjoin(acronyms(removed_idx{r}), ', '));
    end
end

% no value in two regions: two that overlap without one holding the other, or
% that have the same values, would count voxels in two bars
for r = 1:n_regions
    for k = r + 1:n_regions
        shared_values = intersect(own_values{r}, own_values{k});
        if ~isempty(shared_values)
            error(['surprise_regions: the surprise-bar regions %s and %s share ' ...
                   '%d atlas values, and neither holds the other. Keep one of them.'], ...
                  acronyms{r}, acronyms{k}, numel(shared_values));
        end
    end
end
end

function groups = region_divisions(region_values, T_members)
% Each region's division in the atlas (Isocortex, OLF, HPF, CTXsp, STR, PAL, TH,
% HY, MB, ...), from its atlas values; the divisions joined if it spans several.

% each atlas value's division
is_division = strcmp(T_members.parcellation_term_set_name, 'division');
division_values = T_members.parcellation_index(is_division);
division_acronyms = T_members.parcellation_term_acronym(is_division);

% the divisions of each region's values
groups = cell(numel(region_values), 1);
for r = 1:numel(region_values)
    in_region = ismember(division_values, region_values{r});
    groups{r} = strjoin(unique(division_acronyms(in_region)), ', ');
end
end

function [valid_pixels, region_of_voxel] = label_left_hemisphere(AllenCrop, ...
    region_values)
% The atlas voxels of the left hemisphere (the width of the folded maps), and the
% region each falls in (uint16, 0 for none), from regions that do not overlap.

% the left hemisphere's atlas voxels, as a vector
[~, ~, n_width] = size(AllenCrop);
half_width = floor(n_width / 2);
atlas_left = AllenCrop(:, :, 1:half_width);
valid_pixels = atlas_left > 0;
pixel_ids = atlas_left(valid_pixels);

% a table from atlas value to region, value v at row v + 1, then looked up for
% every voxel (with the annotation's integer values as indices, which spares a
% copy of the voxels in double)
all_values = vertcat(region_values{:});
n_table = max(double(max(pixel_ids)), max(all_values)) + 1;
value_to_region = zeros(n_table, 1, 'uint16');
for r = 1:numel(region_values)
    value_to_region(region_values{r} + 1) = r;
end
region_of_voxel = value_to_region(pixel_ids + 1);
end

function T_regions = region_table(acronyms, names, groups, removed_idx, ...
    region_of_voxel)
% One row per region: acronym, name, group (its atlas division), the label of its
% bar (the name, and the regions taken out of it), its voxels in the left
% hemisphere, and the regions taken out of it with their voxels.

n_regions = numel(acronyms);

% each region's voxels
n_voxels = zeros(n_regions, 1);
for r = 1:n_regions
    n_voxels(r) = nnz(region_of_voxel == r);
end

% the regions taken out of each, and their voxels: those regions' own voxels,
% which hold any region inside them in turn
removed = repmat({''}, n_regions, 1);
n_voxels_removed = zeros(n_regions, 1);
labels = names;
for r = 1:n_regions
    if ~isempty(removed_idx{r})
        removed{r} = strjoin(acronyms(removed_idx{r}), ', ');
        n_voxels_removed(r) = sum(n_voxels(removed_idx{r}));
        labels{r} = [names{r} ' (without ' removed{r} ')'];
    end
end

T_regions = table(acronyms, names, groups, labels, n_voxels, removed, ...
    n_voxels_removed, 'VariableNames', {'acronym', 'name', 'group', 'label', ...
    'n_voxels', 'removed', 'n_voxels_removed'});
end

function print_region_table(T_regions)
% The regions, one line each: acronym, group, voxels, name, and the regions taken
% out of it; a warning for the regions with no voxel.

fprintf('  %d regions, voxels in the left hemisphere:\n', height(T_regions));
for r = 1:height(T_regions)
    fprintf('    %-8s %-9s %10d  %s', T_regions.acronym{r}, T_regions.group{r}, ...
        T_regions.n_voxels(r), T_regions.name{r});
    if ~isempty(T_regions.removed{r})
        fprintf(' (without %s: %d voxels)', T_regions.removed{r}, ...
            T_regions.n_voxels_removed(r));
    end
    fprintf('\n');
end

% a region outside the volumes' AP range has no voxel, and no bar
has_no_voxel = T_regions.n_voxels == 0;
if any(has_no_voxel)
    warning('surprise_regions: no voxel in the volumes for %s; no bar.', ...
        strjoin(T_regions.acronym(has_no_voxel), ', '));
end
end
