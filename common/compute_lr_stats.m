function [lr_diff, lr_sum] = compute_lr_stats(volume)
%COMPUTE_LR_STATS  Left minus right, and left plus right, of a volume.
%   [lr_diff, lr_sum] = COMPUTE_LR_STATS(volume) mirrors the right half of
%   volume (AP x DV x ML, or AP x DV x ML x mouse) onto the left along ML and
%   returns their difference and sum, each AP x DV x floor(ML/2) (x mouse).
%   For an odd ML width the last plane is left out, so plane i is paired with
%   plane ML - i rather than its mirror (the registered grids are 1140 wide).
%   A 3D input gives 3D outputs.

[~, ~, n_width, ~] = size(volume);
half_width = floor(n_width / 2);

% left and right halves; for an odd width the right half stops before the last plane
left = volume(:, :, 1:half_width, :);
right_start = half_width + 1;
right_end = min(half_width * 2, n_width);
right = volume(:, :, right_start:right_end, :);

% mirror the right half along ML (dim 3)
flipped_right = flip(right, 3);

% crop the mirrored half to half_width (right_end above already keeps it there)
if size(flipped_right, 3) > half_width
    flipped_right = flipped_right(:, :, 1:half_width, :);
end

lr_diff = left - flipped_right;
lr_sum = left + flipped_right;

% a 3D input has no mouse dimension to keep
if ndims(volume) == 3
    lr_diff = squeeze(lr_diff);
    lr_sum = squeeze(lr_sum);
end

end
