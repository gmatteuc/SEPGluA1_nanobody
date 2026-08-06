function write_lr_indiv_video(lr_diff_4d, lr_sum_4d, mask_bg_4d, atlas_vol, ...
    save_dir, video_filename, clim_diff, clim_sum, ...
    group_name, mouse_names)

% Ensure directory exists
if ~exist(save_dir, 'dir'), mkdir(save_dir); end

full_video_path = fullfile(save_dir, video_filename);
vidObj = VideoWriter(full_video_path, 'MPEG-4');
vidObj.FrameRate = 15;
vidObj.Quality = 95;
open(vidObj);

[n_slices, ~, ~, n_mice] = size(lr_diff_4d);

% Handle atlas slicing
n_depth_data = size(lr_diff_4d, 3);
atlas_subset = double(atlas_vol(:,:,1:n_depth_data));

fprintf('Writing individual video: %s\n', video_filename);

% Initialize figure
fh = figure('visible', 'off', 'units', 'normalized', 'outerposition', [0 0 1 1], 'Color', 'k');
set(fh, 'InvertHardcopy', 'off');

for j = 1:n_slices

    % Skip empty atlas slices
    atlas_slice = squeeze(atlas_subset(j, :, :));
    if sum(atlas_slice(:) > 0) == 0
        if mod(j, 50) == 0, fprintf('  Skipping slice %d\n', j); end
        continue;
    end

    % Compute atlas boundaries
    [gy, gx] = gradient(atlas_slice);
    boundaries = (abs(gx) + abs(gy)) > 0 & (atlas_slice > 0);
    [b_row, b_col] = find(boundaries); %#ok<ASGLU>

    for k = 1:n_mice
        mouse_name = strrep(mouse_names{k}, '_', ' ');

        % Get data for this mouse
        data_diff = squeeze(lr_diff_4d(j, :, :, k));
        data_sum  = squeeze(lr_sum_4d(j, :, :, k));
        mask_bg   = squeeze(mask_bg_4d(j, :, :, k));

        % Valid pixels
        valid_pixels = (atlas_slice > 0) & (~mask_bg);
        alpha_data = double(valid_pixels);

        % Row 1: LR difference
        subplot(2, n_mice, k);
        imagesc(data_diff);
        clim(clim_diff);
        colormap(gca, hot);
        set(findobj(gca, 'Type', 'image'), 'AlphaData', alpha_data);
        axis image; axis off; set(gca, 'Color', 'k');
        hold on;

        % Overlay atlas
        % plot(b_col, b_row, '.', 'Color', [0.5 0.5 0.5], 'MarkerSize', 0.1);

        title(mouse_name, 'Color', 'w', 'FontSize', 10, 'Interpreter', 'none');

        % Labels
        ylabel('LR Diff', 'Color','w','FontSize',12,'FontWeight','bold');
        cb = colorbar('Location', 'westoutside');
        cb.Label.String = '|L - R|'; cb.Color = 'w'; cb.Label.Color = 'w';

        % row 2: LR sum
        subplot(2, n_mice, k + n_mice);
        imagesc(data_sum);
        clim(clim_sum);
        colormap(gca, hot);
        set(findobj(gca, 'Type', 'image'), 'AlphaData', alpha_data);
        axis image; axis off; set(gca, 'Color', 'k');
        hold on;

        % Overlay atlas
        % plot(b_col, b_row, '.', 'Color', [0.5 0.5 0.5], 'MarkerSize', 0.1);

        ylabel('LR Sum', 'Color','w','FontSize',12,'FontWeight','bold');
        cb = colorbar('Location', 'westoutside');
        cb.Label.String = 'L + R'; cb.Color = 'w'; cb.Label.Color = 'w';
    end

    sgtitle(['Slice # ' num2str(j) ' - Left Hemishpere asymmetry (' group_name,')'], 'Color', 'w', 'FontSize', 14);

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