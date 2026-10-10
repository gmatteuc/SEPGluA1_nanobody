function [geom, is_cand] = band_geometry(stack, masks, perm_settings)
%BAND_GEOMETRY  The candidate voxels of the band, as region_permutation_test takes them.
%   [geom, is_cand] = BAND_GEOMETRY(stack, masks, perm_settings) returns the
%   geometry region_permutation_test needs for a stack of band_stack, and
%   which of the band's voxels are candidates: those where at least twice
%   min_mice_per_group mice have a value, the fewest with which a t is
%   possible. The region's voxels are labelled 1, the band around it 0, which
%   the rolling median reads but no cluster takes. Which mice have a value does
%   not depend on their groups, so one set of candidates serves every split of
%   the mice. Used by per_mouse_region_values and mouse_influence.

is_cand = sum(~isnan(stack), 2) >= 2 * perm_settings.min_mice_per_group;
geom = struct();
geom.grid_size = masks.box_size;
geom.cand_lin = masks.band_lin(is_cand);
geom.cand_region = uint16(masks.band_in_region(is_cand));
geom.n_regions = 1;
geom.voxel_mm = 0.01;
end
