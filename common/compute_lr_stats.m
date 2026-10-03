function [lr_diff, lr_sum] = compute_lr_stats(volume)
%COMPUTE_LR_STATS  Left minus right, and left plus right, of a volume.
%   [lr_diff, lr_sum] = COMPUTE_LR_STATS(volume) mirrors the right half of
%   volume (AP x DV x ML, or AP x DV x ML x mouse) onto the left along ML and
%   returns their difference and sum, each AP x DV x floor(ML/2) (x mouse).
%   Plane i is paired with its mirror, plane ML + 1 - i; for an odd ML width
%   the middle plane, the midline, is left out (the registered grids are 1140
%   wide, so even). A 3D input gives 3D outputs.

[~, ~, n_width, ~] = size(volume);
half_width = floor(n_width / 2);

% left and right halves, half_width planes each; for an odd width the middle
% plane belongs to neither
left = volume(:, :, 1:half_width, :);
right_start = n_width - half_width + 1;
right = volume(:, :, right_start:n_width, :);

% mirror the right half along ML (dim 3)
flipped_right = flip(right, 3);

lr_diff = left - flipped_right;
lr_sum = left + flipped_right;

% a 3D input has no mouse dimension to keep
if ndims(volume) == 3
    lr_diff = squeeze(lr_diff);
    lr_sum = squeeze(lr_sum);
end

end
