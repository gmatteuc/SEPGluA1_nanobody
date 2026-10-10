function raw_idx = collected_index(names, group, M_raw, raw_var_name, raw_file)
%COLLECTED_INDEX  The mice's places in a group's collected stack.
%   raw_idx = COLLECTED_INDEX(names, group, M_raw, raw_var_name, raw_file)
%   returns the place of each mouse of names in the collected stack
%   raw_var_name of raw_file (open as the matfile M_raw), which
%   run_collect_by_group fills in the cohort table's order; stops if the stack
%   does not hold the group's mice of the table. Used by
%   per_mouse_region_values and mouse_influence.

cohort = get_cohort('groups', {group});
raw_names = {cohort.name};
raw_size = size(M_raw, raw_var_name);
if numel(raw_size) < 4
    raw_size(4) = 1;
end
if raw_size(4) ~= numel(raw_names)
    error(['collected_index: %s holds %d mice, the cohort table %d for %s. ' ...
           'Rerun run_collect_by_group.'], raw_file, raw_size(4), numel(raw_names), ...
           group);
end
[~, raw_idx] = ismember(names, raw_names);
end
