function fh=plot_abs_slice(data_vol, mask_vol, atlas_vol, brain_mask, slice_num, clim_values, cb_label, titlestring)

    if nargin < 6 || isempty(clim_values)
        min_clim = quantile(data_vol(:), 0.01);
        max_clim = quantile(data_vol(:), 0.99);
        clim_values = [min_clim,max_clim];
    end
    if nargin < 7 || isempty(cb_label)
        cb_label = 'Absolute intensity';
    end
    
    atlasim = squeeze(atlas_vol(slice_num, :, :));
    atlasim = single(atlasim);
    av_warp_boundaries = gradient(atlasim) ~= 0 & (atlasim > 1);
    [row, col] = ind2sub(size(atlasim), find(av_warp_boundaries));
    mask_borders = squeeze(brain_mask(slice_num, :, :));
    
    fh=figure('visible', 'on', 'units', 'normalized', 'outerposition', [0 0 1 1]);
    h = imagesc(squeeze(data_vol(slice_num, :, :)));
    set(h, 'AlphaData', squeeze(mask_borders .* squeeze(mask_vol(slice_num, :, :))));
    clim(clim_values);
    colormap('gray');
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
    title([titlestring,sprintf(' - Slice # %d', slice_num)]);
    
end