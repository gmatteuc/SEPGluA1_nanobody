function masks = region_box_masks(T_regions, valid_pixels, region_of_voxel, regions, ...
    cluster_region, slab_range)
%REGION_BOX_MASKS  The voxels of the regions named, and the box the clusters are read in.
%   masks = REGION_BOX_MASKS(T_regions, valid_pixels, region_of_voxel, regions,
%   cluster_region, slab_range) takes the regions of the bars (surprise_regions)
%   and returns, as linear indices into the folded grid (the left hemisphere),
%   the voxels of each region named in regions (region_lin, one cell each) and
%   of the isocortex, every isocortical area of the bars (iso_lin); and the box
%   around cluster_region that the clusters are read in: its planes and the
%   slab_range planes on each side, which the rolling median of a voxel in it
%   reads (the band).
%
%   The box: its planes (box_ap), rows (box_dv), columns of the left hemisphere
%   (box_ml_left) and, after them, their mirror images in the right hemisphere
%   (box_ml), so compute_lr_stats pairs each column of a mouse's box with its
%   mirror image as it pairs them across the whole width; box_size, the size of
%   the folded box; band_lin, the band's voxels in it; band_in_region, which of
%   them are in cluster_region. Used by per_mouse_region_values and
%   mouse_influence.

% the region of every voxel of the folded grid, 0 for none
region_vol = zeros(size(valid_pixels), 'uint16');
region_vol(valid_pixels) = region_of_voxel;

% each region named in advance, by its row of the bars' table
masks = struct();
masks.folded_size = size(valid_pixels);
masks.region_lin = cell(numel(regions), 1);
for r = 1:numel(regions)
    row = find(strcmp(T_regions.acronym, regions{r}));
    if isempty(row)
        error(['region_box_masks: the region %s is not one of the bars'' regions. ' ...
               'Use their atlas acronyms: %s.'], regions{r}, ...
               strjoin(T_regions.acronym, ', '));
    end
    masks.region_lin{r} = uint32(find(region_vol == row));
    fprintf('  %s: %d voxels in the left hemisphere.\n', regions{r}, ...
        numel(masks.region_lin{r}));
end

% the isocortex: every isocortical area of the bars
iso_rows = find(strcmp(T_regions.group, 'Isocortex'));
masks.iso_lin = uint32(find(ismember(region_vol, iso_rows)));
fprintf('  isocortex: %d areas, %d voxels in the left hemisphere.\n', ...
    numel(iso_rows), numel(masks.iso_lin));

% the band of the clusters: the region and slab_range planes on each side along
% AP, within the atlas
in_region = region_vol == find(strcmp(T_regions.acronym, cluster_region));
band = imdilate(in_region, true(2 * slab_range + 1, 1)) & valid_pixels;
clear region_vol

% the band's bounding box: its planes and rows, its columns of the left
% hemisphere and, after them, their mirror images in the right hemisphere
[ap, dv, ml] = ind2sub(size(band), find(band));
masks.box_ap = min(ap):max(ap);
masks.box_dv = min(dv):max(dv);
masks.box_ml_left = min(ml):max(ml);
masks.box_ml = [masks.box_ml_left, flip(2 * size(band, 3) + 1 - masks.box_ml_left)];
clear ap dv ml

% the band and the region in the box's folded grid
band_box = band(masks.box_ap, masks.box_dv, masks.box_ml_left);
region_box = in_region(masks.box_ap, masks.box_dv, masks.box_ml_left);
masks.box_size = size(band_box);
masks.band_lin = uint32(find(band_box));
masks.band_in_region = region_box(masks.band_lin);
fprintf(['  box of the clusters around %s: %d x %d x %d voxels, %d in the band, %d in ' ...
         'the region.\n'], cluster_region, masks.box_size, numel(masks.band_lin), ...
        nnz(masks.band_in_region));
end
