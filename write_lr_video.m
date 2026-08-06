function write_lr_video(lr_diff_vol, lr_sum_vol, atlas_vol, brain_mask, save_dir, video_filename, clim_values, group_name, label_string_diff, label_string_sum)

% Ensure directory exists
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

full_video_path = fullfile(save_dir, video_filename);
vidObj = VideoWriter(full_video_path, 'MPEG-4');
vidObj.FrameRate = 15;
vidObj.Quality = 95;
open(vidObj);
n_slices = size(lr_diff_vol, 1);

fprintf('Writing video: %s\n', video_filename);

for j = 1:n_slices

    if sum(sum(brain_mask(j,:,:))) == 0
        continue;
    end

    atlasim = squeeze(atlas_vol(j, :, :));
    atlasim = single(atlasim);
    av_warp_boundaries = gradient(atlasim) ~= 0 & (atlasim > 1);
    [row, col] = ind2sub(size(atlasim), find(av_warp_boundaries));

    fh = figure('visible', 'off', 'units', 'normalized', 'outerposition', [0 0 1 1], 'Color', 'k');

    set(fh, 'InvertHardcopy', 'off');

    subplot(1, 2, 1);
    h1 = imagesc(squeeze(lr_diff_vol(j, :, :)));
    clim(clim_values);

    set(h1, 'AlphaData', squeeze(brain_mask(j, :, 1:size(lr_diff_vol, 3))));

    if abs(clim_values(1)) == abs(clim_values(2))
        try
            colormap(gca, get_color2color_colormap([0, 0, 1], [1, 0, 0]));
        catch
            colormap(gca, jet); 
        end
    else
        colormap(gca, hot);
    end

    ax1 = gca;
    ax1.Color = 'k';
    axis equal; axis off;
    hold on;

    line(col, row, 'Marker', '.', 'LineStyle', 'none', 'Color', [0.66 0.66 0.66], 'MarkerSize', 0.5);

    xlim([0, size(lr_diff_vol,3)]);

    title([group_name ' - ' label_string_diff], 'Color', 'w', 'FontSize', 12);

    cb1 = colorbar;
    cb1.Label.String = label_string_diff;
    cb1.Label.FontSize = 10;
    cb1.Color = 'w';
    cb1.Label.Color = 'w';

    subplot(1, 2, 2);
    h2 = imagesc(squeeze(lr_sum_vol(j, :, :)));
    clim(2*clim_values);

    set(h2, 'AlphaData', squeeze(brain_mask(j, :, 1:size(lr_sum_vol, 3))));

    if abs(clim_values(1)) == abs(clim_values(2))
        try
            colormap(gca, get_color2color_colormap([0, 0, 1], [1, 0, 0]));
        catch
            colormap(gca, jet);
        end
    else
        colormap(gca, hot);
    end

    ax2 = gca;
    ax2.Color = 'k';
    axis equal; axis off;
    hold on;

    line(col, row, 'Marker', '.', 'LineStyle', 'none', 'Color', [0.66 0.66 0.66], 'MarkerSize', 0.5);

    xlim([0, size(lr_diff_vol,3)]);
    title([group_name ' - ' label_string_sum], 'Color', 'w', 'FontSize', 12);

    cb2 = colorbar;
    cb2.Label.String = label_string_sum;
    cb2.Label.FontSize = 10;
    cb2.Color = 'w';
    cb2.Label.Color = 'w';

    sgtitle(['Slice # ' num2str(j)], 'Color', 'w', 'FontSize', 14);

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