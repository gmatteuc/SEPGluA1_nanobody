function [level, plane_levels] = off_tissue_level(raw, bg_mask, brainMask)
%OFF_TISSUE_LEVEL  A mouse's off-tissue level in its collected stack.
%   [level, plane_levels] = OFF_TISSUE_LEVEL(raw, bg_mask, brainMask) returns
%   the median of the collected volume raw (AP x DV x ML) over the mouse's
%   background voxels (bg_mask, run_normalise_groups' mask) outside the atlas
%   brain (brainMask) that a section reached: the slide around the section.
%   Over the whole stack (level), and plane by plane (plane_levels, NaN for a
%   plane with fewer than 1000 of them). Used by per_mouse_region_values and
%   scale_free_test.

% run_normalise_groups' mask marks, plane by plane, the voxels below the knee
% between background and tissue; outside the atlas brain they are the slide
% around the section, and a raw 0 is a voxel no section reached
is_off = bg_mask & ~brainMask & raw > 0;
level = double(median(raw(is_off)));

% plane by plane, for the spread along AP; 1000 voxels give a stable median
n_planes = size(raw, 1);
plane_levels = nan(n_planes, 1);
for z = 1:n_planes
    plane_raw = raw(z, :, :);
    plane_off = is_off(z, :, :);
    if nnz(plane_off) >= 1000
        plane_levels(z) = double(median(plane_raw(plane_off)));
    end
end
end
