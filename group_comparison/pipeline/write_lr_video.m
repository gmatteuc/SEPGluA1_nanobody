function write_lr_video(lr_diff_vol, lr_sum_vol, atlas_vol, brain_mask, save_dir, ...
    video_filename, clim_values, group_name, label_string_diff, label_string_sum)
%WRITE_LR_VIDEO  Video of a left-right difference map and its sum map, plane by plane.
%   WRITE_LR_VIDEO(lr_diff_vol, lr_sum_vol, atlas_vol, brain_mask, save_dir,
%   video_filename, clim_values, group_name, label_string_diff,
%   label_string_sum) writes save_dir\video_filename (MPEG-4, 15 frames per
%   second), one frame per AP plane that has a voxel in brain_mask: the
%   difference on the left, with colour limits clim_values, the sum on the
%   right, with limits 2 * clim_values, both over the atlas boundaries.
%   Symmetric limits give the blue-red colormap, others hot.
%
%   lr_diff_vol, lr_sum_vol   AP x DV x ML, one hemisphere
%   atlas_vol                 the annotation the boundaries are drawn from
%   brain_mask                AP x DV x ML, the voxels shown (1) or hidden (0)
%   group_name                the start of each panel's title
%   label_string_diff, label_string_sum
%                             the end of each panel's title, and its colorbar label
%
%   Run by group_differences, and by P8.

% open the video, in a folder made if needed
[vidObj, full_video_path] = open_lr_video(save_dir, video_filename);
n_slices = size(lr_diff_vol, 1);

fprintf('Writing video: %s\n', video_filename);

for j = 1:n_slices

    % skip the planes with no brain
    if sum(sum(brain_mask(j, :, :))) == 0
        continue;
    end

    % atlas boundaries: where the annotation changes along ML
    [row, col] = lr_atlas_boundaries(atlas_vol, j);

    fh = figure('visible', 'off', 'units', 'normalized', 'outerposition', [0 0 1 1], ...
        'Color', 'k');

    set(fh, 'InvertHardcopy', 'off');

    % left: the difference
    subplot(1, 2, 1);
    h1 = imagesc(squeeze(lr_diff_vol(j, :, :)));
    clim(clim_values);

    set(h1, 'AlphaData', squeeze(brain_mask(j, :, 1:size(lr_diff_vol, 3))));

    % blue-red for symmetric limits (jet if the colormap function is missing)
    set_lr_colormap(clim_values);

    ax1 = gca;
    ax1.Color = 'k';
    axis equal;
    axis off;
    hold on;

    line(col, row, 'Marker', '.', 'LineStyle', 'none', 'Color', [0.66 0.66 0.66], ...
        'MarkerSize', 0.5);

    xlim([0, size(lr_diff_vol, 3)]);

    title([group_name ' - ' label_string_diff], 'Color', 'w', 'FontSize', 12);

    cb1 = colorbar;
    cb1.Label.String = label_string_diff;
    cb1.Label.FontSize = 10;
    cb1.Color = 'w';
    cb1.Label.Color = 'w';

    % right: the sum, on twice the limits
    subplot(1, 2, 2);
    h2 = imagesc(squeeze(lr_sum_vol(j, :, :)));
    clim(2*clim_values);

    set(h2, 'AlphaData', squeeze(brain_mask(j, :, 1:size(lr_sum_vol, 3))));

    set_lr_colormap(clim_values);

    ax2 = gca;
    ax2.Color = 'k';
    axis equal;
    axis off;
    hold on;

    line(col, row, 'Marker', '.', 'LineStyle', 'none', 'Color', [0.66 0.66 0.66], ...
        'MarkerSize', 0.5);

    xlim([0, size(lr_diff_vol, 3)]);
    title([group_name ' - ' label_string_sum], 'Color', 'w', 'FontSize', 12);

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
