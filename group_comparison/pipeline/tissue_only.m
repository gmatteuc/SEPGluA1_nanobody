function vol = tissue_only(vol, bg_mask, brainMask, apply_smoothing, smooth_sigma)
%TISSUE_ONLY  One mouse's volume, NaN outside its tissue, smoothed if asked.
%   vol = TISSUE_ONLY(vol, bg_mask, brainMask, apply_smoothing, smooth_sigma)
%   sets to NaN every voxel of vol (AP x DV x ML) outside the mouse's tissue:
%   outside the atlas brain (brainMask), in its background mask (bg_mask,
%   true for background), or without a value. With apply_smoothing, the
%   tissue is then smoothed in 3D by a normalised convolution with a
%   Gaussian of sigma smooth_sigma, in 10 um voxels. Used by
%   group_differences and per_mouse_region_values, so both smooth alike.

% the tissue: in the atlas brain, outside the background, reached by a section
% (run_normalise_groups left NaN where none was)
is_tissue = brainMask & ~bg_mask & ~isnan(vol);
vol(~is_tissue) = NaN;

% normalised convolution, both with the same Gaussian: the smoothed values (NaN
% as 0) over the smoothed tissue mask, so the voxels outside the tissue neither
% count as zero nor spread into it; outside the tissue it stays NaN, since a 0
% there would enter the group means and the counts of mice as a measured value
if apply_smoothing
    smoothed_values = imgaussfilt3(fillmissing(vol, 'constant', 0), smooth_sigma);
    smoothed_mask = imgaussfilt3(double(is_tissue), smooth_sigma);
    vol = smoothed_values ./ smoothed_mask;
    vol(~is_tissue) = NaN;
end
end
