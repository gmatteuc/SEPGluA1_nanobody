function write_diff_video(data_vol, mask_vol, atlas_vol, brain_mask, save_dir, video_filename, clim_value, cb_label)

full_video_path = fullfile(save_dir, video_filename);
vidObj = VideoWriter(full_video_path, 'MPEG-4');
vidObj.FrameRate = 1;
vidObj.Quality = 100;
open(vidObj);

for j = 1:size(data_vol, 1)
    atlasim = squeeze(atlas_vol(j, :, :));
    atlasim = single(atlasim);
    av_warp_boundaries = gradient(atlasim) ~= 0 & (atlasim > 1);
    [row, col] = ind2sub(size(atlasim), find(av_warp_boundaries));
    mask_borders = squeeze(brain_mask(j, :, :));
    figure('visible', 'off', 'units', 'normalized', 'outerposition', [0 0 1 1]);
    h = imagesc(squeeze(data_vol(j, :, :)));
    set(h, 'AlphaData', squeeze(mask_borders .* squeeze(mask_vol(j, :, :))));
    clim([-clim_value, clim_value]);
    colormap(get_color2color_colormap([0, 0, 1], [1, 0, 0]));
    cb = colorbar;
    cb.Label.String = cb_label;
    cb.Label.FontSize = 12;
    ax = gca;
    ax.Color = 'k';
    axis equal;
    grid on;
    line(col, row, 'Marker', '.', 'LineStyle', 'none', 'Color', [0.66 0.66 0.66], 'MarkerSize', 0.5)
    ylim([0, size(data_vol, 2)]);
    xlim([0, size(data_vol, 3)]);
    title(sprintf('Average Diff - Slice # %d', j));
    frame = getframe(gcf);
    writeVideo(vidObj, frame);
    close(gcf);
    if mod(j,10)==0
        disp(['Frame ',num2str(j), ' written', ' - ',video_filename])
    end
end
close(vidObj);

end