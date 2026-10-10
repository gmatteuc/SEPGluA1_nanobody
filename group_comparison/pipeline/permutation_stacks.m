function [stacks, geom] = permutation_stacks(maps, valid_pixels, region_of_voxel, ...
    n_regions, min_mice_per_group)
%PERMUTATION_STACKS  Every mouse's maps on the voxels a split can give a t, for the region test.
%   [stacks, geom] = PERMUTATION_STACKS(maps, valid_pixels, region_of_voxel,
%   n_regions, min_mice_per_group) takes the folded maps the region test
%   compares, maps, a cell with one {control, experimental} pair per map, each
%   AP x DV x ML x mice on the folded grid (NaN where a mouse has no value),
%   and returns for region_permutation_test each map's absolute values (the
%   values the group means and SEMs take) on the candidate voxels, n_cand x
%   mice in single, control mice first, and the candidates' geometry. The
%   candidates are the atlas voxels of the left hemisphere (valid_pixels)
%   where at least twice min_mice_per_group mice have a value in the first
%   map, the fewest with which some split gives a t; region_of_voxel gives
%   the region of each of valid_pixels' voxels, 0 for none, of n_regions.
%   Used by group_differences, with L - R and L + R (the sum has the same
%   missing voxels as the difference), and by scale_free_test, with the
%   asymmetry index.

% the candidates: atlas voxels of the left hemisphere where at least twice
% min_mice_per_group mice have a value, the fewest with which some split gives a
% t
n_with_value = zeros(size(valid_pixels), 'uint8');
for g = 1:2
    for i = 1:size(maps{1}{g}, 4)
        n_with_value = n_with_value + uint8(~isnan(maps{1}{g}(:, :, :, i)));
    end
end
candidates = valid_pixels & n_with_value >= 2 * min_mice_per_group;
clear n_with_value

% their region, through the left hemisphere's voxels
region_vol = zeros(size(valid_pixels), 'like', region_of_voxel);
region_vol(valid_pixels) = region_of_voxel;

% the candidates, with the grid of the folded maps: the 10 um voxels of the
% registered volumes (get_atlas_crop)
geom = struct();
geom.grid_size = size(valid_pixels);
geom.cand_lin = uint32(find(candidates));
geom.cand_region = region_vol(candidates);
geom.n_regions = n_regions;
geom.voxel_mm = 0.01;
clear region_vol candidates

% each map's values, mouse by mouse, by linear index into the 4D maps
n_cand = numel(geom.cand_lin);
n_voxels = prod(geom.grid_size);
stacks = cell(1, numel(maps));
for m = 1:numel(maps)
    n_mice = size(maps{m}{1}, 4) + size(maps{m}{2}, 4);
    stacks{m} = zeros(n_cand, n_mice, 'single');
    column = 0;
    for g = 1:2
        for i = 1:size(maps{m}{g}, 4)
            column = column + 1;
            values = maps{m}{g}(double(geom.cand_lin) + (i - 1) * n_voxels);
            stacks{m}(:, column) = single(abs(values));
        end
    end
end
fprintf('  %d candidate voxels; stacks of %d mice, %d maps, %.1f GB.\n', n_cand, ...
    size(stacks{1}, 2), numel(maps), numel(maps) * n_cand * size(stacks{1}, 2) * 4 / 1e9);
end
