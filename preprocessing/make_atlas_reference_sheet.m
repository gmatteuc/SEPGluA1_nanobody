%% make_atlas_reference_sheet
% ===== Coronal plates of the atlas, spaced as the sections =====
%
% A tool, run by hand. Draws one sheet of coronal plates of the atlas template,
% slice_spacing_um apart, so plate N and plate N+1 are one section apart in our
% data: to keep open next to SliceOrderEditor when deciding the order of the
% slices. Each plate is labelled with its number and its distance from the
% anterior end of the cropped range, in mm.
%
% Setup: the CCF template, a plate every 150 um, written to
% young\atlas_coronal_reference.png. Run sep_setup_paths first, once per MATLAB
% session.

close all
clear all
clc

%% Settings

% project folders, worked out from where the code sits, so the tree can be moved
% or copied to another drive as is (SEP_DATA_ROOT points the data elsewhere)
paths = get_paths();

% atlas ('ccf')
atlas_key = 'ccf';

% spacing between plates, in um; equal to slicethickness, so one plate is one
% section
slice_spacing_um = 150;

% the sheet
out_file = fullfile(paths.data, 'young', 'atlas_coronal_reference.png');

% plates per row
n_cols = 8;

%% Load the atlas template

atlas = get_atlas(atlas_key);
template_path = fullfile(atlas.dir, atlas.template_file);

fprintf('Loading %s ...\n', template_path);
vol = niftiread(template_path);

% the volume is AP x DV x ML; crop AP to the limits the pipeline uses, so the
% plates cover the region analysed
ap_limits = atlas.default_aplims;
vol = vol(ap_limits(1):ap_limits(2), :, :);

fprintf('Cropped volume is %d x %d x %d (AP x DV x ML)\n', size(vol, 1), size(vol, 2), ...
    size(vol, 3));

%% Pick the AP positions to show

step = round(slice_spacing_um / atlas.res_um);
ap_positions = 1:step:size(vol, 1);
n_plates = numel(ap_positions);

fprintf('Sampling every %d voxels (%d um): %d plates\n', step, slice_spacing_um, ...
    n_plates);

%% Build the montage

% plates at 0.35 of their size, so the sheet stays manageable
scale = 0.35;

plates = cell(n_plates, 1);
for k = 1:n_plates
    plane = squeeze(vol(ap_positions(k), :, :));
    plates{k} = imresize(plane, scale);
end

plate_h = size(plates{1}, 1);
plate_w = size(plates{1}, 2);

% the plates in rows of n_cols on one canvas
n_rows = ceil(n_plates / n_cols);
canvas = zeros(n_rows * plate_h, n_cols * plate_w, 'like', plates{1});

for k = 1:n_plates
    row = floor((k - 1) / n_cols);
    col = mod(k - 1, n_cols);
    canvas(row * plate_h + (1:plate_h), col * plate_w + (1:plate_w)) = plates{k};
end

%% Draw it and label each plate

fig = figure('Visible', 'off', 'Color', 'k', 'Units', 'pixels', ...
    'Position', [100 100 min(2200, n_cols * plate_w) min(2200, n_rows * plate_h)]);

ax = axes('Parent', fig, 'Position', [0 0 1 1]);
imshow(canvas, [], 'Parent', ax);
colormap(ax, gray);

for k = 1:n_plates
    row = floor((k - 1) / n_cols);
    col = mod(k - 1, n_cols);

    % distance from the anterior end of the cropped range, enough to order the
    % sections against each other
    depth_mm = (ap_positions(k) - 1) * atlas.res_um / 1000;

    label = sprintf('%d   %.2f mm', k, depth_mm);
    text(ax, col * plate_w + 6, row * plate_h + 16, label, ...
        'Color', 'y', 'FontSize', 10, 'FontWeight', 'bold', 'Interpreter', 'none');
end

exportgraphics(fig, out_file, 'Resolution', 130);
close(fig);

fprintf('Wrote %s\n', out_file);
fprintf('Plate 1 is the anterior end of the cropped range (atlas AP index %d).\n', ...
    ap_limits(1));
