function stack = band_stack(boxes, is_exp, slope, intercept, common_factor, masks)
%BAND_STACK  The mice's |L - R| on the band, as the stacks of the region test.
%   stack = BAND_STACK(boxes, is_exp, slope, intercept, common_factor, masks)
%   returns one column per mouse, in the order of boxes (each mouse's box,
%   region_box_masks), of its |L - R| at the band's voxels (masks.band_lin), on
%   the common scale (aligned_box_lr: the mice with is_exp through the line):
%   n_band x n_mice single, NaN where the mouse has no value, the stack
%   region_permutation_test takes. Used by per_mouse_region_values and
%   mouse_influence.

stack = zeros(numel(masks.band_lin), numel(boxes), 'single');
for k = 1:numel(boxes)
    lr_diff = aligned_box_lr(boxes{k}, is_exp(k), slope, intercept, common_factor);
    stack(:, k) = abs(lr_diff(masks.band_lin));
end
end
