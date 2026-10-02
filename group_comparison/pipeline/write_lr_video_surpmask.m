function write_lr_video_surpmask(lr_diff_vol, lr_sum_vol, atlas_vol, brain_mask, ...
    save_dir, video_filename, clim_values, ...
    group_name, label_string_diff, label_string_sum, ...
    surp_vol, surp_thresh)
%WRITE_LR_VIDEO_SURPMASK  Video of left-right maps, shown where the surprise is high.
%   WRITE_LR_VIDEO_SURPMASK(lr_diff_vol, lr_sum_vol, atlas_vol, brain_mask,
%   save_dir, video_filename, clim_values, group_name, label_string_diff,
%   label_string_sum, surp_vol, surp_thresh) writes save_dir\video_filename
%   (MPEG-4, 15 frames per second) as WRITE_LR_VIDEO does, with both panels on
%   the limits clim_values and each voxel's opacity surp_vol / surp_thresh,
%   clipped to 0-1: a voxel at the threshold surprise (-log10 p) or above is
%   opaque. The one surprise volume sets the opacity of both panels.
%
%   Run by group_differences.
%
%   See also WRITE_LR_VIDEO, WRITE_LR_VIDEO_SURPMASK_ROLLING.

% create the folder if needed
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

% open the video
full_video_path = fullfile(save_dir, video_filename);
vidObj = VideoWriter(full_video_path, 'MPEG-4');
vidObj.FrameRate = 15;
vidObj.Quality = 95;
open(vidObj);
n_slices = size(lr_diff_vol, 1);

fprintf('Writing video: %s\n', video_filename);

for j = 1:n_slices

    % skip the planes with no brain
    if sum(sum(brain_mask(j, :, :))) == 0
        continue;
    end

    % atlas boundaries: where the annotation changes along ML
    atlasim = squeeze(atlas_vol(j, :, :));
    atlasim = single(atlasim);
    av_warp_boundaries = gradient(atlasim) ~= 0 & (atlasim > 1);
    [row, col] = ind2sub(size(atlasim), find(av_warp_boundaries));

    % opacity from the surprise, clipped to 0-1, and transparent where it is NaN
    surp_slice = squeeze(surp_vol(j, :, :));
    alpha_surp = surp_slice ./ surp_thresh;
    alpha_surp(alpha_surp < 0) = 0;
    alpha_surp(alpha_surp > 1) = 1;
    alpha_surp(isnan(alpha_surp)) = 0;

    brain_slice = squeeze(brain_mask(j, :, 1:size(lr_diff_vol, 3)));
    brain_slice = double(brain_slice);

    % combine it with the brain mask
    alpha_mask = alpha_surp .* brain_slice;

    fh = figure('visible', 'off', 'units', 'normalized', 'outerposition', [0 0 1 1], ...
        'Color', 'k');
    set(fh, 'InvertHardcopy', 'off');

    % left: the difference
    subplot(1, 2, 1);
    h1 = imagesc(squeeze(lr_diff_vol(j, :, :)));
    clim(clim_values);

    set(h1, 'AlphaData', alpha_mask);

    % blue-red for symmetric limits (jet if the colormap function is missing)
    if abs(clim_values(1)) == abs(clim_values(2))
        try
            colormap(gca, sep_palette('difference'));
        catch
            colormap(gca, jet);
        end
    else
        colormap(gca, sep_palette('intensity'));
    end

    ax1 = gca;
    ax1.Color = 'k';
    axis equal;
    axis off;
    hold on;

    line(col, row, 'Marker', '.', 'LineStyle', 'none', ...
        'Color', [0.66 0.66 0.66], 'MarkerSize', 0.5);

    xlim([0, size(lr_diff_vol, 3)]);

    title([group_name ' - ', label_string_diff], 'Color', 'w', 'FontSize', 12);
    cb1 = colorbar;
    cb1.Label.String = label_string_diff;
    cb1.Label.FontSize = 10;
    cb1.Color = 'w';
    cb1.Label.Color = 'w';

    % right: the sum
    subplot(1, 2, 2);
    h2 = imagesc(squeeze(lr_sum_vol(j, :, :)));
    clim(clim_values);

    set(h2, 'AlphaData', alpha_mask);

    if abs(clim_values(1)) == abs(clim_values(2))
        try
            colormap(gca, sep_palette('difference'));
        catch
            colormap(gca, jet);
        end
    else
        colormap(gca, sep_palette('intensity'));
    end

    ax2 = gca;
    ax2.Color = 'k';
    axis equal;
    axis off;
    hold on;

    line(col, row, 'Marker', '.', 'LineStyle', 'none', ...
        'Color', [0.66 0.66 0.66], 'MarkerSize', 0.5);

    xlim([0, size(lr_diff_vol, 3)]);

    title([group_name ' - ', label_string_sum], 'Color', 'w', 'FontSize', 12);
    cb2 = colorbar;
    cb2.Label.String = label_string_sum;
    cb2.Label.FontSize = 10;
    cb2.Color = 'w';
    cb2.Label.Color = 'w';

    sgtitle(['Slice # ' num2str(j)], 'Color', 'w', 'FontSize', 14);

    % write the frame
    frame = getframe(fh);
    writeVideo(vidObj, frame);
    close(fh);

    if mod(j, 50) == 0
        fprintf('  Frame %d written...\n', j);
    end
end
close(vidObj);
fprintf('Video saved: %s\n', full_video_path);
end
