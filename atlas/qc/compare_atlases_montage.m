clear all
close all
clc

% /// QC helper: adult CCF against DeMBA P20, coronal, at matched levels ///
% Puts the two reference atlases side by side so the actual anatomical
% differences between a P56 and a P20 brain can be judged by eye:
%   (1) Crops each atlas to its own atlasaplims, which are the spans that
%       correspond to the same anatomy in the two volumes
%   (2) Brings both to 20 um in plane, so a pixel means the same thing in
%       either image and the two are literally the same scale
%   (3) Draws a montage, one column per matched level, adult on top and P20
%       below, with the two brain outlines superimposed underneath
%   (4) Prints the measured width and height of each brain at each level
%
% The two atlases are NOT the same modality of intensity (both are serial
% two-photon, but scaled differently), so each is contrast-normalised on its
% own percentiles. Compare shape and size here, not brightness.
%
% Matched levels come from the fractional position within each crop rather
% than from a plane index, because the two volumes have different AP extents
% and origins. That is only as good as the crops themselves -- see the note on
% default_aplims in get_atlas.m.

%% User-defined parameters

% Where the project lives. Derived from the location of the code rather than
% written out, so the tree can be moved or copied to another drive as is.
paths = get_paths();

% The two atlases to compare, as get_atlas keys
atlas_key_a = 'ccf';                % top row
atlas_key_b = 'demba_p20';          % bottom row

% Where along each crop to cut, as a fraction from the anterior end
level_fracs = [0.05 0.20 0.35 0.50 0.65 0.80 0.95];

% Common display resolution. Both atlases are brought to this in plane, so the
% montage shows the two brains at the same physical scale.
display_res_um = 20;

% Contrast normalisation percentiles for the template images
clim_pcts = [1 99.5];

% Draw the third row with the two outlines superimposed
show_outline_overlay = true;

% Draw a fourth row of sagittal outlines. The coronal rows are matched level by
% level, which by construction hides the difference in AP length; the sagittal
% view is where that difference is visible, so it is worth its own row.
show_sagittal_row = true;

% Where to cut the sagittal planes, as a fraction of the ML half-width from the
% midline outward. Both atlases have the same ML extent, so the same fraction
% is the same place in both.
sagittal_ml_fracs = [0.00 0.08 0.16 0.24 0.32 0.40 0.48];

% Length of the scale bar, in mm
scalebar_mm = 1;

% Output
out_dir = fullfile(paths.data, 'young', 'registration_qc');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
save_figure = true;

% Color palette
outline_a_color = [0.35 0.35 0.35];   % gray, adult (atlas gray, as elsewhere)
outline_b_color = [0.95 0.55 0.10];   % strong orange, P20
squeezed_color  = [0.20 0.55 0.85];   % blue, P20 squeezed to adult length
outline_width   = 1.3;                % outlines are drawn as lines, not painted pixels
label_color     = [0.15 0.15 0.15];

%% Load both atlases

% Read by explicit full path rather than through which(), because this is the
% one script that legitimately wants both atlases in memory at once and
% get_atlas deliberately keeps only one of them on the MATLAB path.
atlas_a = get_atlas(atlas_key_a);
atlas_b = get_atlas(atlas_key_b);

fprintf('A: %s\n   %s\n', atlas_a.description, atlas_a.dir);
fprintf('B: %s\n   %s\n', atlas_b.description, atlas_b.dir);

tv_a = niftiread(fullfile(atlas_a.dir, atlas_a.template_file));
av_a = niftiread(fullfile(atlas_a.dir, atlas_a.annotation_file));
tv_b = niftiread(fullfile(atlas_b.dir, atlas_b.template_file));
av_b = niftiread(fullfile(atlas_b.dir, atlas_b.annotation_file));

% Crop each to its own AP span
lim_a = atlas_a.default_aplims;
lim_b = atlas_b.default_aplims;
tv_a = tv_a(lim_a(1):lim_a(2), :, :);
av_a = av_a(lim_a(1):lim_a(2), :, :);
tv_b = tv_b(lim_b(1):lim_b(2), :, :);
av_b = av_b(lim_b(1):lim_b(2), :, :);

n_ap_a = size(tv_a, 1);
n_ap_b = size(tv_b, 1);

fprintf('\n%s: %d AP planes at %g um = %.2f mm\n', ...
    atlas_a.key, n_ap_a, atlas_a.res_um, n_ap_a * atlas_a.res_um / 1000);
fprintf('%s: %d AP planes at %g um = %.2f mm\n', ...
    atlas_b.key, n_ap_b, atlas_b.res_um, n_ap_b * atlas_b.res_um / 1000);

% Scale factors that bring each atlas to the common display resolution
scale_a = atlas_a.res_um / display_res_um;
scale_b = atlas_b.res_um / display_res_um;

%% Build the matched planes

n_levels = numel(level_fracs);
plane_a  = cell(1, n_levels);
plane_b  = cell(1, n_levels);
mask_a   = cell(1, n_levels);
mask_b   = cell(1, n_levels);

fprintf('\n%-8s %-14s %-14s %-18s %-18s\n', ...
    'level', 'adult plane', 'P20 plane', 'adult w x h (mm)', 'P20 w x h (mm)');
fprintf('%s\n', repmat('-', 1, 76));

for k = 1:n_levels

    f = level_fracs(k);

    % Fractional position, not plane index: the two crops hold different
    % numbers of planes even though they span the same anatomy.
    i_a = max(1, min(n_ap_a, round(1 + f * (n_ap_a - 1))));
    i_b = max(1, min(n_ap_b, round(1 + f * (n_ap_b - 1))));

    % Take the coronal plane first, then rescale in 2D. Resizing the whole
    % volume in 3D would interpolate along AP as well, which would blur the
    % very levels we are trying to compare.
    img_a = single(squeeze(tv_a(i_a, :, :)));
    img_b = single(squeeze(tv_b(i_b, :, :)));
    lbl_a = squeeze(av_a(i_a, :, :));
    lbl_b = squeeze(av_b(i_b, :, :));

    if scale_a ~= 1
        img_a = imresize(img_a, scale_a, 'bilinear');
        lbl_a = imresize(lbl_a, scale_a, 'nearest');
    end
    if scale_b ~= 1
        img_b = imresize(img_b, scale_b, 'bilinear');
        lbl_b = imresize(lbl_b, scale_b, 'nearest');
    end

    plane_a{k} = img_a;
    plane_b{k} = img_b;
    mask_a{k}  = lbl_a > 0;
    mask_b{k}  = lbl_b > 0;

    % Measured extent of the brain at this level, so the size difference is a
    % number and not just an impression
    [wa, ha] = mask_extent_mm(mask_a{k}, display_res_um);
    [wb, hb] = mask_extent_mm(mask_b{k}, display_res_um);

    fprintf('%-8.2f %-14d %-14d %-18s %-18s\n', f, i_a, i_b, ...
        sprintf('%.2f x %.2f', wa, ha), sprintf('%.2f x %.2f', wb, hb));

end

%% Put every panel on one common canvas

% So that a brain that is genuinely smaller looks smaller, rather than each
% panel being cropped to its own content.
all_h = cellfun(@(x) size(x, 1), [plane_a plane_b]);
all_w = cellfun(@(x) size(x, 2), [plane_a plane_b]);
canvas_h = max(all_h);
canvas_w = max(all_w);

for k = 1:n_levels
    plane_a{k} = pad_to_canvas(plane_a{k}, canvas_h, canvas_w, 0);
    plane_b{k} = pad_to_canvas(plane_b{k}, canvas_h, canvas_w, 0);
    mask_a{k}  = pad_to_canvas(mask_a{k},  canvas_h, canvas_w, false);
    mask_b{k}  = pad_to_canvas(mask_b{k},  canvas_h, canvas_w, false);
end

% Trim the empty border that both atlases carry around the tissue. One common
% box for every panel, so the panels stay comparable -- it only removes
% background, and stops the montage being mostly black.
any_content = false(canvas_h, canvas_w);
for k = 1:n_levels
    any_content = any_content | mask_a{k} | mask_b{k};
end
rows = find(any(any_content, 2));
cols = find(any(any_content, 1));
margin_px = round(0.04 * canvas_w);
r0 = max(1, rows(1) - margin_px);   r1 = min(canvas_h, rows(end) + margin_px);
c0 = max(1, cols(1) - margin_px);   c1 = min(canvas_w, cols(end) + margin_px);

for k = 1:n_levels
    plane_a{k} = plane_a{k}(r0:r1, c0:c1);
    plane_b{k} = plane_b{k}(r0:r1, c0:c1);
    mask_a{k}  = mask_a{k}(r0:r1, c0:c1);
    mask_b{k}  = mask_b{k}(r0:r1, c0:c1);
end
canvas_h = r1 - r0 + 1;
canvas_w = c1 - c0 + 1;

%% Sagittal planes, for the length comparison

% Built on their own common canvas, because a sagittal plane is AP x DV and so
% has nothing to do with the coronal canvas above. Both are brought to the same
% micrometres per pixel, and neither is stretched to fit, so the difference in
% AP length is shown at true scale.
if show_sagittal_row

    n_sag = numel(sagittal_ml_fracs);
    sag_a = cell(1, n_sag);
    sag_b = cell(1, n_sag);
    sag_c = cell(1, n_sag);

    ml_a = size(av_a, 3);
    ml_b = size(av_b, 3);

    fprintf('\nSagittal outlines, AP length of the labelled brain at each ML level:\n');
    fprintf('%-10s %-18s %-18s %-8s\n', 'ML frac', 'adult AP (mm)', 'P20 AP (mm)', 'ratio');
    fprintf('%s\n', repmat('-', 1, 58));

    for k = 1:n_sag

        f = sagittal_ml_fracs(k);

        % Measured outward from the ML midline, same physical place in both
        i_a = round(ml_a / 2 + f * ml_a / 2);
        i_b = round(ml_b / 2 + f * ml_b / 2);
        i_a = min(max(i_a, 1), ml_a);
        i_b = min(max(i_b, 1), ml_b);

        m_a = squeeze(av_a(:, :, i_a)) > 0;
        m_b = squeeze(av_b(:, :, i_b)) > 0;

        if scale_a ~= 1
            m_a = imresize(m_a, scale_a, 'nearest');
        end
        if scale_b ~= 1
            m_b = imresize(m_b, scale_b, 'nearest');
        end

        % The volumes are AP x DV, which would draw a sagittal view on its side.
        % Transposing puts AP along x and DV down y; the atlases are stored
        % anterior-first and superior-first, so this comes out anterior-left and
        % dorsal-up, the way a sagittal section is normally read.
        ap_a = ap_extent_mm(m_a, display_res_um);
        ap_b = ap_extent_mm(m_b, display_res_um);

        sag_a{k} = m_a';
        sag_b{k} = m_b';

        % A third outline: the young one squeezed uniformly along AP until it is
        % the adult's length. If that lands on top of the adult outline, the
        % difference between the two atlases is a plain proportional stretch
        % rather than a change of shape -- which is not something anchoring the
        % anterior ends together can tell you on its own.
        sq = imresize(m_b', [size(m_b, 2), round(size(m_b, 1) * ap_a / ap_b)], 'nearest');
        sag_c{k} = sq; %#ok<SAGROW>
        fprintf('%-10.2f %-18.2f %-18.2f %-8.3f\n', f, ap_a, ap_b, ap_b / ap_a);

    end

    sag_h = max(cellfun(@(x) size(x, 1), [sag_a sag_b sag_c]));
    sag_w = max(cellfun(@(x) size(x, 2), [sag_a sag_b sag_c]));

    % Anchored at the anterior end rather than centred, so the extra length
    % accumulates at one end and can actually be read off the picture.
    for k = 1:n_sag
        sag_a{k} = pad_anterior(sag_a{k}, sag_h, sag_w);
        sag_b{k} = pad_anterior(sag_b{k}, sag_h, sag_w);
        sag_c{k} = pad_anterior(sag_c{k}, sag_h, sag_w);
    end

end

%% Figure

n_rows_grid = 2 + double(show_outline_overlay) + double(show_sagittal_row);

% Size the figure to the grid rather than to the screen. With 'axis image' on
% wide, short coronal planes, a full-screen figure leaves most of each tile
% empty and the montage comes out mostly white space.
tile_w_px = 300;
tile_h_px = tile_w_px * canvas_h / canvas_w;
fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'pixels', ...
             'Position', [50 50 ...
                          round(n_levels * tile_w_px + 120), ...
                          round(n_rows_grid * tile_h_px + 110)]);
tl = tiledlayout(fig, n_rows_grid, n_levels, ...
                 'TileSpacing', 'compact', 'Padding', 'compact');

% Row 1 - adult
for k = 1:n_levels
    nexttile(tl, k);
    show_plane(plane_a{k}, clim_pcts);
    title(sprintf('%.0f%%', 100 * level_fracs(k)), ...
          'FontWeight', 'normal', 'Color', label_color);
    if k == 1
        ylabel(sprintf('%s  (%g um)', atlas_a.key, atlas_a.res_um), ...
               'Interpreter', 'none', 'Visible', 'on');
    end
end

% Row 2 - young
for k = 1:n_levels
    nexttile(tl, n_levels + k);
    show_plane(plane_b{k}, clim_pcts);
    if k == 1
        ylabel(sprintf('%s  (%g um)', atlas_b.key, atlas_b.res_um), ...
               'Interpreter', 'none', 'Visible', 'on');
    end
    if k == n_levels
        draw_scalebar(gca, scalebar_mm, display_res_um, canvas_h, canvas_w);
    end
end

% Row 3 - the two outlines on top of each other, which is where a difference
% in size or shape actually becomes visible
if show_outline_overlay
    for k = 1:n_levels
        nexttile(tl, 2 * n_levels + k);
        draw_outline(mask_b{k}, outline_b_color, outline_width);
        draw_outline(mask_a{k}, outline_a_color, outline_width);
        xlim([1 canvas_w]); ylim([1 canvas_h]);
        set(gca, 'YDir', 'reverse'); axis image off
        if k == 1
            ylabel('outlines', 'Visible', 'on');
        end
    end
end

% Row 4 - sagittal outlines. The coronal rows are matched level by level, which
% deliberately cancels the AP length difference; here it is left in.
if show_sagittal_row
    row_offset = (2 + double(show_outline_overlay)) * n_levels;
    for k = 1:min(n_sag, n_levels)
        nexttile(tl, row_offset + k);
        draw_outline(sag_b{k}, outline_b_color, outline_width);
        draw_outline(sag_c{k}, squeezed_color,  outline_width);
        draw_outline(sag_a{k}, outline_a_color, outline_width);
        xlim([1 sag_w]); ylim([1 sag_h]);
        set(gca, 'YDir', 'reverse'); axis image off
        title(sprintf('ML %.0f%%', 100 * sagittal_ml_fracs(k)), ...
              'FontWeight', 'normal', 'Color', label_color);
        if k == 1
            ylabel('sagittal', 'Visible', 'on');
            draw_scalebar(gca, scalebar_mm, display_res_um, sag_h, sag_w);
        end
    end
end

title(tl, sprintf(['Adult CCF vs DeMBA P20, both at %g um, same scale throughout' ...
                   '   |   gray = %s, orange = %s\n' ...
                   'sagittal: anterior ends anchored together, so the AP difference ' ...
                   'accumulates rightwards.\nBlue = P20 squeezed uniformly to the ' ...
                   'adult length -- where it lands on gray, the difference is pure scale'], ...
                   display_res_um, atlas_a.key, atlas_b.key), ...
      'Interpreter', 'none');

if save_figure
    png_name = fullfile(out_dir, 'atlas_comparison_adult_vs_p20.png');
    fig_name = fullfile(out_dir, 'atlas_comparison_adult_vs_p20.fig');
    exportgraphics(fig, png_name, 'Resolution', 160);
    savefig(fig, fig_name);
    fprintf('\nsaved:\n  %s\n  %s\n', png_name, fig_name);
end

%% Local function: physical width and height of a brain mask

function [w_mm, h_mm] = mask_extent_mm(mask, res_um)
    rows = find(any(mask, 2));
    cols = find(any(mask, 1));
    if isempty(rows)
        w_mm = 0; h_mm = 0;
        return
    end
    h_mm = (rows(end) - rows(1) + 1) * res_um / 1000;
    w_mm = (cols(end) - cols(1) + 1) * res_um / 1000;
end

%% Local function: AP extent of a sagittal mask, in mm

function ap_mm = ap_extent_mm(mask, res_um)
    rows = find(any(mask, 2));
    if isempty(rows)
        ap_mm = 0;
        return
    end
    ap_mm = (rows(end) - rows(1) + 1) * res_um / 1000;
end

%% Local function: put a sagittal mask on a canvas, anchored at the anterior end

function out = pad_anterior(img, canvas_h, canvas_w)
    out = false(canvas_h, canvas_w);
    cols = find(any(img, 1));
    if isempty(cols)
        return
    end
    % Trim to the brain, then place its anterior tip at a fixed column, so the
    % two atlases start together at the left and any extra length runs off the
    % posterior end where it can be read off directly.
    img = img(:, cols(1):cols(end));
    [h, w] = size(img);
    w = min(w, canvas_w);
    r0 = floor((canvas_h - h) / 2) + 1;
    r0 = max(r0, 1);
    h = min(h, canvas_h - r0 + 1);
    out(r0:r0 + h - 1, 1:w) = img(1:h, 1:w);
end

%% Local function: centre an image on a common canvas

function out = pad_to_canvas(img, canvas_h, canvas_w, fill_value)
    out = repmat(cast(fill_value, 'like', img), canvas_h, canvas_w);
    [h, w] = size(img);
    r0 = floor((canvas_h - h) / 2) + 1;
    c0 = floor((canvas_w - w) / 2) + 1;
    out(r0:r0 + h - 1, c0:c0 + w - 1) = img;
end

%% Local function: draw one template plane in gray

function show_plane(img, clim_pcts)
    inside = img(img > 0);
    if isempty(inside)
        lo = 0; hi = 1;
    else
        lo = prctile(inside, clim_pcts(1));
        hi = prctile(inside, clim_pcts(2));
    end
    if hi <= lo, hi = lo + 1; end
    imagesc(img, [lo hi]);
    colormap(gca, gray);
    axis image off
end

%% Local function: draw a mask's boundary as continuous lines

function draw_outline(mask, color, lw)
    % Painting bwperim into an RGB image gave a one-pixel outline that broke up
    % into dashes once the panel was scaled down. Tracing the boundary and
    % drawing it as a line stays continuous at any size.
    hold on
    B = bwboundaries(mask, 'noholes');
    for k = 1:numel(B)
        if size(B{k}, 1) < 20
            continue    % specks, not anatomy
        end
        plot(B{k}(:, 2), B{k}(:, 1), '-', 'Color', color, 'LineWidth', lw);
    end
end

%% Local function: scale bar, so the montage carries its own ruler

function draw_scalebar(ax, bar_mm, res_um, canvas_h, canvas_w)
    bar_px = bar_mm * 1000 / res_um;
    x0 = canvas_w - bar_px - round(0.05 * canvas_w);
    y0 = canvas_h - round(0.08 * canvas_h);
    hold(ax, 'on')
    plot(ax, [x0 x0 + bar_px], [y0 y0], 'w-', 'LineWidth', 3);
    text(ax, x0 + bar_px / 2, y0 - round(0.04 * canvas_h), ...
         sprintf('%g mm', bar_mm), 'Color', 'w', ...
         'HorizontalAlignment', 'center', 'FontSize', 9);
    hold(ax, 'off')
end
