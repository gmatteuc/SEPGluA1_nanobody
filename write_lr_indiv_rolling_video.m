function write_lr_indiv_rolling_video(diff_4d, sum_4d, bg_mask_4d, atlas_vol, ...
    save_dir, video_filename, group_name, mouse_names, slab_range, clim_diff, clim_sum)

% Ensure directory exists
if ~exist(save_dir, 'dir'), mkdir(save_dir); end

full_video_path = fullfile(save_dir, video_filename);
vidObj = VideoWriter(full_video_path, 'MPEG-4');
vidObj.FrameRate = 15;
vidObj.Quality = 95;
open(vidObj);

[n_slices, ~, ~, n_mice] = size(diff_4d);

fprintf('Writing Rolling Video: %s (Slab +/- %d)\n', video_filename, slab_range);

% Initialize figure
fh = figure('visible', 'off', 'units', 'normalized', 'outerposition', [0 0 1 1], 'Color', 'k');
set(fh, 'InvertHardcopy', 'off');

for j = 1:n_slices

    % Define Slab Indices
    z_start = max(1, j - slab_range);
    z_end   = min(n_slices, j + slab_range);
    z_indices = z_start:z_end;

    % Skip empty atlas slices (optimization)
    atlas_slice = squeeze(atlas_vol(j, :, :));
    if sum(atlas_slice(:) > 0) == 0
        if mod(j, 100) == 0, fprintf('  Skipping slice %d\n', j); end
        continue;
    end

    % Prepare atlas overlay
    [gy, gx] = gradient(single(atlas_slice));
    boundaries = (abs(gx) + abs(gy)) > 0 & (atlas_slice > 0);
    [b_row, b_col] = find(boundaries);

    for k = 1:n_mice
        mouse_name = strrep(mouse_names{k}, '_', ' ');

        % Extract slabs
        raw_slab_d = diff_4d(z_indices, :, :, k);
        raw_slab_s = sum_4d(z_indices, :, :, k);
        mask_slab  = bg_mask_4d(z_indices, :, :, k);

        % Apply NaN to background
        raw_slab_d(logical(mask_slab)) = NaN;
        raw_slab_s(logical(mask_slab)) = NaN;

        % Compute rolling median
        slab_diff = squeeze(nanmedian(raw_slab_d, 1)); %#ok<NANMEDIAN>
        slab_sum  = squeeze(nanmedian(raw_slab_s, 1));

        % Prepare alpha mask
        slab_bg = squeeze(min(mask_slab, [], 1));
        valid_pixels = (atlas_slice > 0) & (~slab_bg);
        alpha_data = double(valid_pixels);

        % --- Plot diff (row 1) ---
        subplot(2, n_mice, k);
        imagesc(slab_diff);
        set(findobj(gca, 'Type', 'image'), 'AlphaData', alpha_data);
        clim(clim_diff); colormap(gca, hot);
        axis image; axis off; set(gca, 'Color', 'k'); hold on;
        plot(b_col, b_row, '.', 'Color', [0.5 0.5 0.5], 'MarkerSize', 0.1);
        title(mouse_name, 'Color', 'w', 'FontSize', 10, 'Interpreter', 'none');
        ylabel('LR Diff (Slab)', 'Color','w','FontSize',12,'FontWeight','bold');
        cb = colorbar('Location', 'westoutside');
        cb.Label.String = '|L - R|'; cb.Color = 'w';

        % --- Plot sum (row 2) ---
        subplot(2, n_mice, k + n_mice);
        imagesc(slab_sum);
        set(findobj(gca, 'Type', 'image'), 'AlphaData', alpha_data);
        clim(clim_sum); colormap(gca, hot);
        axis image; axis off; set(gca, 'Color', 'k'); hold on;
        plot(b_col, b_row, '.', 'Color', [0.5 0.5 0.5], 'MarkerSize', 0.1);
        ylabel('LR Sum (Slab)', 'Color','w','FontSize',12,'FontWeight','bold');
        cb = colorbar('Location', 'westoutside');
        cb.Label.String = 'L + R'; cb.Color = 'w';
    end

    sgtitle(['Slice ' num2str(j) ' (Slab \pm' num2str(slab_range) ') - ' group_name], 'Color', 'w', 'FontSize', 14);

    frame = getframe(fh);
    writeVideo(vidObj, frame);
    clf(fh);

    if mod(j, 50) == 0
        fprintf('  Frame %d written...\n', j);
    end
end

close(vidObj);
close(fh);
fprintf('Video saved: %s\n', full_video_path);
end