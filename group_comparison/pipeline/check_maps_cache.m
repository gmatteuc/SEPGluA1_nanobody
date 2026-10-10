function check_maps_cache(G, group_dir, channel, cache_file)
%CHECK_MAPS_CACHE  Stop unless cached maps were made from the group's normalised stack.
%   CHECK_MAPS_CACHE(G, group_dir, channel, cache_file) stops unless the
%   normalised stack of the group in group_dir still holds the mice and the
%   lines of run_normalise_groups that the cached maps G (a group of
%   Per_Mouse_Maps_<tag>.mat, in cache_file) were made from: after step 2 is
%   run again, the cache would otherwise give the old maps without a word.
%   Used by per_mouse_region_values and mouse_influence.

norm_file = fullfile(group_dir, [channel '_4d_normalized.mat']);
S_norm = load(norm_file, 'current_mice', 'norm_params');
if ~isfield(G, 'stack_norm_params') || ~isequal(G.stack_mice, S_norm.current_mice) ...
        || ~isequal(G.stack_norm_params, S_norm.norm_params)
    error(['check_maps_cache: the maps of %s in %s were not made from the ' ...
           'normalised stack now in %s (its mice or their lines differ, or the ' ...
           'cache predates the check). Rerun run_per_mouse_values with ' ...
           'force_recompute_mice = true.'], G.group, cache_file, norm_file);
end
end
