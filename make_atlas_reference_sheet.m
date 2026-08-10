close all
clear all
clc

% Makes a coronal reference sheet from the Allen template, sampled at the same
% AP spacing as our sections. Handy to keep open next to SliceOrderEditor when
% deciding what order the slices should go in.
%
% The plates are spaced slice_spacing_um apart, so plate N and plate N+1 are one
% section apart in our data. That makes it easy to walk through your sections and
% the reference side by side.

%% User-defined parameters

% Reference atlas
atlas_key = 'ccf';

% Spacing between plates, in micrometres. Keep this equal to slicethickness
% so one plate corresponds to one section.
slice_spacing_um = 150;

% Where to write the sheet
out_file = 'D:\sep_histology\data\young\_atlas_coronal_reference.png';

% How many plates per row
n_cols = 8;

%% Load the atlas template

atlas = get_atlas(atlas_key);
template_path = fullfile(atlas.dir, atlas.template_file);

fprintf('Loading %s ...\n', template_path);
vol = niftiread(template_path);

% The volume is AP x DV x ML. We crop the AP range to the same limits the
% pipeline uses, so the plates cover exactly the region we analyse.
ap_limits = atlas.default_aplims;
vol = vol(ap_limits(1):ap_limits(2), :, :);

fprintf('Cropped volume is %d x %d x %d (AP x DV x ML)\n', size(vol,1), size(vol,2), size(vol,3));

%% Pick the AP positions to show

step = round(slice_spacing_um / atlas.res_um);
ap_positions = 1:step:size(vol,1);
n_plates = numel(ap_positions);

fprintf('Sampling every %d voxels (%d um): %d plates\n', step, slice_spacing_um, n_plates);

%% Build the montage

% Scale the plates down a bit, otherwise the sheet gets unwieldy
scale = 0.35;

plates = cell(n_plates, 1);
for k = 1:n_plates
    plane = squeeze(vol(ap_positions(k), :, :));
    plates{k} = imresize(plane, scale);
end

plate_h = size(plates{1}, 1);
plate_w = size(plates{1}, 2);

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

    % Distance from the front of the cropped range, which is all we need to
    % order sections relative to each other
    depth_mm = (ap_positions(k) - 1) * atlas.res_um / 1000;

    label = sprintf('%d   %.2f mm', k, depth_mm);
    text(ax, col * plate_w + 6, row * plate_h + 16, label, ...
        'Color', 'y', 'FontSize', 10, 'FontWeight', 'bold', 'Interpreter', 'none');
end

exportgraphics(fig, out_file, 'Resolution', 130);
close(fig);

fprintf('Wrote %s\n', out_file);
fprintf('Plate 1 is the anterior end of the cropped range (atlas AP index %d).\n', ap_limits(1));
