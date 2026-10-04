%% compare_atlases_montage
% ===== The adult CCF and DeMBA P20 side by side, coronal, at matched levels =====
%
% Atlas QC, run by hand, so the anatomical differences between a P56 and a
% P20 brain can be judged by eye:
%   1. crops each atlas to its own atlasaplims, the spans that hold the same
%      anatomy in the two volumes
%   2. brings both to 20 um in plane, so a pixel means the same in either
%      image and the two are at the same scale
%   3. draws a montage, one column per matched level, adult on top and P20
%      below, then the two brain outlines on top of each other, then sagittal
%      outlines
%   4. prints the measured width and height of each brain at each level
% Saves atlas_comparison_adult_vs_p20.png and .fig in young\registration_qc
% under the data root.
%
% The two templates are both serial two-photon images but scaled differently,
% so each is contrast-normalised on its own percentiles: compare shape and
% size here, not brightness. Matched levels are fractional positions within
% each crop rather than plane indices, because the two volumes differ in AP
% extent and origin; that is only as good as the crops themselves (see
% default_aplims in get_atlas).
%
% Setup: the adult CCF against DeMBA P20. Run sep_setup_paths first, once
% per MATLAB session.

clear; clc; close all;

%% Settings

% project folders, worked out from where the code sits, so the tree can be moved
% or copied to another drive as is (SEP_DATA_ROOT points the data elsewhere)
paths = get_paths();

% the two atlases to compare, as get_atlas keys: the top row, then the bottom
atlas_key_a = 'ccf';
atlas_key_b = 'demba_p20';

% where along each crop to cut, as a fraction from the anterior end
level_fracs = [0.05 0.20 0.35 0.50 0.65 0.80 0.95];

% common display resolution in plane, in um, so the montage shows the two
% brains at the same physical scale
display_res_um = 20;

% contrast normalisation percentiles of the template images
clim_pcts = [1 99.5];

% draw a third row with the two outlines on top of each other
show_outline_overlay = true;

% draw a fourth row of sagittal outlines: the coronal rows are matched level by
% level, which hides the difference in AP length; the sagittal view shows it
show_sagittal_row = true;

% where to cut the sagittal planes, as a fraction of the ML half-width from the
% midline outward; both atlases have the same ML extent, so the same fraction
% is the same place in both
sagittal_ml_fracs = [0.00 0.08 0.16 0.24 0.32 0.40 0.48];

% length of the scale bar, in mm
scalebar_mm = 1;

% output folder, and whether to save the figure
out_dir = fullfile(paths.data, 'young', 'registration_qc');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end
save_figure = true;

% colours: the adult outline grey (atlas grey, as elsewhere), P20 orange, P20
% squeezed to the adult length blue, labels dark grey; outlines are drawn as
% lines of this width, not painted pixels
outline_a_color = [0.35 0.35 0.35];
outline_b_color = [0.95 0.55 0.10];
squeezed_color  = [0.20 0.55 0.85];
outline_width   = 1.3;
label_color     = [0.15 0.15 0.15];

%% Load both atlases

% read by full path rather than through which(): this script wants both
% atlases in memory at once, and get_atlas keeps only one of them on the path
atlas_a = get_atlas(atlas_key_a);
atlas_b = get_atlas(atlas_key_b);

fprintf('A: %s\n   %s\n', atlas_a.description, atlas_a.dir);
fprintf('B: %s\n   %s\n', atlas_b.description, atlas_b.dir);

tv_a = niftiread(fullfile(atlas_a.dir, atlas_a.template_file));
av_a = niftiread(fullfile(atlas_a.dir, atlas_a.annotation_file));
tv_b = niftiread(fullfile(atlas_b.dir, atlas_b.template_file));
av_b = niftiread(fullfile(atlas_b.dir, atlas_b.annotation_file));

% crop each to its own AP span
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

% scale factors that bring each atlas to the common display resolution
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

    % fractional position, not plane index: the two crops hold different
    % numbers of planes although they span the same anatomy
    i_a = max(1, min(n_ap_a, round(1 + f * (n_ap_a - 1))));
    i_b = max(1, min(n_ap_b, round(1 + f * (n_ap_b - 1))));

    % the coronal plane first, then rescaled in 2D: resizing the whole volume
    % would interpolate along AP too, blurring the very levels compared
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

    % the measured extent of the brain at this level, so the size difference
    % is a number and not an impression
    [wa, ha] = mask_extent_mm(mask_a{k}, display_res_um);
    [wb, hb] = mask_extent_mm(mask_b{k}, display_res_um);

    fprintf('%-8.2f %-14d %-14d %-18s %-18s\n', f, i_a, i_b, ...
        sprintf('%.2f x %.2f', wa, ha), sprintf('%.2f x %.2f', wb, hb));

end

%% Put every panel on one common canvas

% so a brain that is smaller looks smaller, rather than each panel being
% cropped to its own content
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

% trim the empty border both atlases carry around the tissue, with one box for
% every panel, so the panels stay comparable and the montage is not mostly black
any_content = false(canvas_h, canvas_w);
for k = 1:n_levels
    any_content = any_content | mask_a{k} | mask_b{k};
end
rows = find(any(any_content, 2));
cols = find(any(any_content, 1));
margin_px = round(0.04 * canvas_w);
r0 = max(1, rows(1) - margin_px);
r1 = min(canvas_h, rows(end) + margin_px);
c0 = max(1, cols(1) - margin_px);
c1 = min(canvas_w, cols(end) + margin_px);

for k = 1:n_levels
    plane_a{k} = plane_a{k}(r0:r1, c0:c1);
    plane_b{k} = plane_b{k}(r0:r1, c0:c1);
    mask_a{k}  = mask_a{k}(r0:r1, c0:c1);
    mask_b{k}  = mask_b{k}(r0:r1, c0:c1);
end
canvas_h = r1 - r0 + 1;
canvas_w = c1 - c0 + 1;

%% Sagittal planes, for the length comparison

% on a canvas of their own (a sagittal plane is AP x DV), both at the same um
% per pixel and neither stretched to fit, so the AP difference is at true scale
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

        % measured outward from the ML midline, the same place in both
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

        % transposed, AP along x and DV down y: the atlases are stored anterior
        % first and superior first, so anterior-left and dorsal-up, as usual
        ap_a = ap_extent_mm(m_a, display_res_um);
        ap_b = ap_extent_mm(m_b, display_res_um);

        sag_a{k} = m_a';
        sag_b{k} = m_b';

        % a third outline, the young one squeezed along AP to the adult's
        % length: on top of the adult outline, the difference is a plain stretch
        sq = imresize(m_b', [size(m_b, 2), round(size(m_b, 1) * ap_a / ap_b)], 'nearest');
        sag_c{k} = sq; %#ok<SAGROW>
        fprintf('%-10.2f %-18.2f %-18.2f %-8.3f\n', f, ap_a, ap_b, ap_b / ap_a);

    end

    sag_h = max(cellfun(@(x) size(x, 1), [sag_a sag_b sag_c]));
    sag_w = max(cellfun(@(x) size(x, 2), [sag_a sag_b sag_c]));

    % anchored at the anterior end rather than centred, so the extra length
    % gathers at one end, where it can be read off
    for k = 1:n_sag
        sag_a{k} = pad_anterior(sag_a{k}, sag_h, sag_w);
        sag_b{k} = pad_anterior(sag_b{k}, sag_h, sag_w);
        sag_c{k} = pad_anterior(sag_c{k}, sag_h, sag_w);
    end

end

%% Figure

n_rows_grid = 2 + double(show_outline_overlay) + double(show_sagittal_row);

% the figure sized to the grid rather than the screen: with 'axis image' on
% wide, short coronal planes, a full-screen figure is mostly white space
tile_w_px = 300;
tile_h_px = tile_w_px * canvas_h / canvas_w;
fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'pixels', ...
             'Position', [50 50 ...
                          round(n_levels * tile_w_px + 120), ...
                          round(n_rows_grid * tile_h_px + 110)]);
tl = tiledlayout(fig, n_rows_grid, n_levels, ...
                 'TileSpacing', 'compact', 'Padding', 'compact');

% row 1: adult
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

% row 2: young
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

% row 3: the two outlines on top of each other, where a difference in size or
% shape shows
if show_outline_overlay
    for k = 1:n_levels
        nexttile(tl, 2 * n_levels + k);
        draw_outline(mask_b{k}, outline_b_color, outline_width);
        draw_outline(mask_a{k}, outline_a_color, outline_width);
        xlim([1 canvas_w]);
        ylim([1 canvas_h]);
        set(gca, 'YDir', 'reverse');
        axis image off
        if k == 1
            ylabel('outlines', 'Visible', 'on');
        end
    end
end

% row 4: sagittal outlines; the coronal rows cancel the AP length difference
% level by level, here it is left in
if show_sagittal_row
    row_offset = (2 + double(show_outline_overlay)) * n_levels;
    for k = 1:min(n_sag, n_levels)
        nexttile(tl, row_offset + k);
        draw_outline(sag_b{k}, outline_b_color, outline_width);
        draw_outline(sag_c{k}, squeezed_color,  outline_width);
        draw_outline(sag_a{k}, outline_a_color, outline_width);
        xlim([1 sag_w]);
        ylim([1 sag_h]);
        set(gca, 'YDir', 'reverse');
        axis image off
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

% ===== Local functions =====

function [w_mm, h_mm] = mask_extent_mm(mask, res_um)
% Width and height of a brain mask in mm (0 for an empty mask).

rows = find(any(mask, 2));
cols = find(any(mask, 1));
if isempty(rows)
    w_mm = 0;
    h_mm = 0;
    return
end
h_mm = (rows(end) - rows(1) + 1) * res_um / 1000;
w_mm = (cols(end) - cols(1) + 1) * res_um / 1000;

end

function ap_mm = ap_extent_mm(mask, res_um)
% AP extent of a sagittal mask in mm (0 for an empty mask).

rows = find(any(mask, 2));
if isempty(rows)
    ap_mm = 0;
    return
end
ap_mm = (rows(end) - rows(1) + 1) * res_um / 1000;

end

function out = pad_anterior(img, canvas_h, canvas_w)
% A sagittal mask on a canvas, its anterior tip at the left edge, centred in DV.

out = false(canvas_h, canvas_w);
cols = find(any(img, 1));
if isempty(cols)
    return
end

% trim to the brain, then place its anterior tip at a fixed column, so the two
% atlases start together and any extra length runs off the posterior end
img = img(:, cols(1):cols(end));
[h, w] = size(img);
w = min(w, canvas_w);
r0 = floor((canvas_h - h) / 2) + 1;
r0 = max(r0, 1);
h = min(h, canvas_h - r0 + 1);
out(r0:r0 + h - 1, 1:w) = img(1:h, 1:w);

end

function show_plane(img, clim_pcts)
% Draw one template plane in gray, scaled to the percentiles clim_pcts of its
% tissue.

inside = img(img > 0);
if isempty(inside)
    lo = 0;
    hi = 1;
else
    lo = prctile(inside, clim_pcts(1));
    hi = prctile(inside, clim_pcts(2));
end
if hi <= lo
    hi = lo + 1;
end
imagesc(img, [lo hi]);
colormap(gca, sep_palette('anatomy'));
axis image off

end

function draw_outline(mask, color, lw)
% Draw a mask's boundary as continuous lines, specks under 20 points left out.

% a boundary traced and drawn as a line stays continuous at any size, where a
% painted one-pixel bwperim breaks into dashes once the panel is scaled down
hold on
B = bwboundaries(mask, 'noholes');
for k = 1:numel(B)

    % specks, not anatomy
    if size(B{k}, 1) < 20
        continue
    end
    plot(B{k}(:, 2), B{k}(:, 1), '-', 'Color', color, 'LineWidth', lw);
end

end

function draw_scalebar(ax, bar_mm, res_um, canvas_h, canvas_w)
% A white scale bar of bar_mm in the bottom right corner, so the montage
% carries its own ruler.

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
