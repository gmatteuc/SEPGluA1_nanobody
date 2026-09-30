function write_lr_indiv_video_signed(lr_diff_4d, lr_sum_4d, mask_bg_4d, atlas_vol, ...
    save_dir, video_filename, clim_diff, clim_sum, ...
    group_name, mouse_names)

% Ensure directory exists
if ~exist(save_dir, 'dir'), mkdir(save_dir); end

full_video_path = fullfile(save_dir, video_filename);
vidObj = VideoWriter(full_video_path, 'MPEG-4');
vidObj.FrameRate = 15;
vidObj.Quality = 95;
open(vidObj)

[n_slices, ~, ~, n_mice] = size(lr_diff_4d);

% Handle atlas slicing
n_depth_data = size(lr_diff_4d, 3);
atlas_subset = double(atlas_vol(:,:,1:n_depth_data));

% Get colormap
crwb = get_color2color_colormap([1 0 0],[0 0 1]);

fprintf('Writing signed individual video: %s\n', video_filename);

% Initialize figure
fh = figure('visible', 'off', 'units', 'normalized', 'outerposition', [0 0 1 1], 'Color', 'k');
set(fh, 'InvertHardcopy', 'off');

for j = 1:n_slices

    % Skip empty atlas slices
    atlas_slice = squeeze(atlas_subset(j, :, :));
    if sum(atlas_slice(:) > 0) == 0 %#ok<LOGSUM>
        if mod(j, 50) == 0, fprintf('  Skipping slice %d\n', j); end
        continue;
    end

    % Compute atlas boundaries (optional overlay)
    [gy, gx] = gradient(atlas_slice);
    boundaries = (abs(gx) + abs(gy)) > 0 & (atlas_slice > 0);
    [b_row, b_col] = find(boundaries); %#ok<ASGLU>

    for k = 1:n_mice
        if k > length(mouse_names)
            m_name = sprintf('Mouse %d', k);
        else
            m_name = strrep(mouse_names{k}, '_', ' ');
        end

        % Get data for this mouse
        data_diff = squeeze(lr_diff_4d(j, :, :, k));
        data_sum  = squeeze(lr_sum_4d(j, :, :, k));
        mask_bg   = squeeze(mask_bg_4d(j, :, :, k));

        % Valid pixels
        valid_pixels = (atlas_slice > 0) & (~mask_bg);
        alpha_data = double(valid_pixels);

        % --- Row 1: LR signed difference ---
        subplot(2, n_mice, k);
        imagesc(data_diff);
        clim(clim_diff);
        colormap(gca, crwb);

        set(findobj(gca, 'Type', 'image'), 'AlphaData', alpha_data);
        axis image; axis off; set(gca, 'Color', 'k');
        hold on;

        % Overlay atlas (Optional)
        % plot(b_col, b_row, '.', 'Color', [0.5 0.5 0.5], 'MarkerSize', 0.1);

        title(m_name, 'Color', 'w', 'FontSize', 10, 'Interpreter', 'none');

        ylabel('L - R (Signed)', 'Color','w','FontSize',12,'FontWeight','bold');
        cb = colorbar('Location', 'eastoutside');
        cb.Label.String = 'R > L  |  L > R';
        cb.Color = 'w'; cb.Label.Color = 'w';

        % --- Row 2: LR sum ---
        subplot(2, n_mice, k + n_mice);
        imagesc(data_sum);
        clim(clim_sum);
        colormap(gca, hot); % Keep 'hot' for intensity

        set(findobj(gca, 'Type', 'image'), 'AlphaData', alpha_data);
        axis image; axis off; set(gca, 'Color', 'k');
        hold on;

        % Overlay atlas (Optional)
        % plot(b_col, b_row, '.', 'Color', [0.5 0.5 0.5], 'MarkerSize', 0.1);

        ylabel('L + R (Sum)', 'Color','w','FontSize',12,'FontWeight','bold');
        cb = colorbar('Location', 'eastoutside');
        cb.Label.String = 'Intensity';
        cb.Color = 'w'; cb.Label.Color = 'w';
    end

    sgtitle(['Slice # ' num2str(j) ' - Directional Asymmetry (' group_name,')'], 'Color', 'w', 'FontSize', 14);

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