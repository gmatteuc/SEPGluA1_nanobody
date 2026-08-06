function [lr_diff, lr_sum] = compute_lr_stats(volume)

    [~, ~, n_width, ~] = size(volume);
    half_width = floor(n_width / 2);
    
    % Extract left and right halves
    left = volume(:, :, 1:half_width, :);
    right_start = half_width + 1;
    right_end = min(half_width * 2, n_width);  % Handle odd width
    right = volume(:, :, right_start:right_end, :);
    
    % Flip right along width (dim 3)
    flipped_right = flip(right, 3);
    
    % Crop flipped_right to half_width if necessary (for odd width)
    if size(flipped_right, 3) > half_width
        flipped_right = flipped_right(:, :, 1:half_width, :);
    end
    
    lr_diff = left - flipped_right;
    lr_sum = left + flipped_right;
    
    if ndims(volume) == 3  % Original was 3D, squeeze out singleton mouse dim
        lr_diff = squeeze(lr_diff);
        lr_sum = squeeze(lr_sum);
    end
    
end