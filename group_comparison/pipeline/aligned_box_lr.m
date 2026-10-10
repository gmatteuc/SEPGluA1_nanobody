function [lr_diff, lr_sum] = aligned_box_lr(box, is_exp, slope, intercept, common_factor)
%ALIGNED_BOX_LR  A mouse's box on the common scale of the test, folded.
%   [lr_diff, lr_sum] = ALIGNED_BOX_LR(box, is_exp, slope, intercept,
%   common_factor) puts a mouse's box (its normalised volume in the clusters'
%   box, region_box_masks) on the common scale as run_group_differences puts a
%   whole volume on it, a mouse of the experimental group through the line of
%   align_exp_to_ctrl first, then folds it (compute_lr_stats). Used by
%   per_mouse_region_values, band_stack and mouse_influence.

if is_exp
    vol = ((box .* slope) + intercept) ./ common_factor;
else
    vol = box ./ common_factor;
end
[lr_diff, lr_sum] = compute_lr_stats(vol);
end
