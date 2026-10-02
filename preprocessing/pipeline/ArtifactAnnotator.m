function ArtifactAnnotator(scaledautoVol, nanoVol, bg_mask_vol, slice_data, ...
    mouse_name, output_dir, correction_type)
%ARTIFACTANNOTATOR  Outline artifacts by hand, slice by slice, as polygons.
%   ARTIFACTANNOTATOR(scaledautoVol, nanoVol, bg_mask_vol, slice_data,
%   mouse_name, output_dir, correction_type) shows each slice with the
%   scaled autofluorescence in red and the nano in green, both scaled to the
%   slice's used_clim in slice_data, and the background in grey. The left and
%   right arrows move between slices, d draws a new polygon, delete removes
%   the slice's last one, s saves and escape closes; a right click on a
%   polygon deletes it, and its vertices can be dragged. The window holds
%   MATLAB (uiwait) until it is closed; closing asks whether to save.
%
%   Saves output_dir\artifact_mask_volume_<correction_type>.mat with
%   artifact_mask_vol (H x W x Z logical) and poly_vertices (the polygons of
%   each slice, so the work can go on), and reloads that file when it exists.

[H, W, Z] = size(scaledautoVol);
mask_file = fullfile(output_dir, sprintf('artifact_mask_volume_%s.mat', correction_type));

% the window's data: the volumes, and no polygon on any slice yet
gui_data = struct();
gui_data.scaledautoVol = scaledautoVol;
gui_data.nanoVol = nanoVol;
gui_data.bg_mask_vol = bg_mask_vol;
gui_data.slice_data = slice_data;
gui_data.artifact_mask = false(H, W, Z);
gui_data.poly_vertices = cell(Z, 1);
gui_data.rois = cell(Z, 1);
for zz = 1:Z
    gui_data.poly_vertices{zz} = {};
    gui_data.rois{zz} = gobjects(0);
end

% carry on from the saved annotations, if any
if exist(mask_file, 'file')
    loaded_data = load(mask_file);
    if isfield(loaded_data, 'artifact_mask_vol')
        gui_data.artifact_mask = loaded_data.artifact_mask_vol;
    end
    if isfield(loaded_data, 'poly_vertices')
        gui_data.poly_vertices = loaded_data.poly_vertices;
    else
        % without poly_vertices no polygon can be edited; the mask is kept
    end
    disp(sprintf('Loaded existing annotations for %s from %s', mouse_name, mask_file)); %#ok<DSPSP>
end
gui_data.current_z = 1;
gui_data.Z = Z;
gui_data.mouse_name = mouse_name;
gui_data.output_dir = output_dir;
gui_data.correction_type = correction_type;
gui_data.H = H;
gui_data.W = W;

% the window, with its key, click and close handlers
gui_fig = figure('Name', sprintf('Artifact Annotator - %s', mouse_name), ...
    'NumberTitle', 'off', ...
    'Toolbar', 'none', 'Menubar', 'none', 'Color', 'w', ...
    'WindowState', 'maximized', ...
    'CloseRequestFcn', @(src, evt) callback_close_gui_request(src, evt), ...
    'KeyPressFcn', @callback_keypress, ...
    'Resize', 'on', ...
    'WindowButtonDownFcn', @on_mouse_click);
gui_data.imageAxes = axes('Parent', gui_fig, 'Units', 'normalized', ...
    'Position', [0.05 0.05 0.9 0.9]);

% an RGB image, three channels from the start
gui_data.imageHandle = imshow(zeros(H, W, 3), 'Parent', gui_data.imageAxes);
gui_data.titleHandle = title(gui_data.imageAxes, '', 'FontSize', 12);
guidata(gui_fig, gui_data);

% show the first slice and wait until the window is closed
display_current_slice(gui_fig);
if Z > 0
    uiwait(gui_fig);
else
    if ishandle(gui_fig)
        delete(gui_fig);
    end
    disp('No slices to annotate. GUI closed.');
    return;
end
end

% ===== Local functions =====

function display_current_slice(fig)
% Draw the current slice as an RGB overlay, with its polygons, editable.

gui_data = guidata(fig);
z = gui_data.current_z;
if z < 1 || z > gui_data.Z
    return;
end

% the slice's display limits, 0 to 1 when it has none
used_clim = gui_data.slice_data(z).used_clim;
if any(isnan(used_clim))
    used_clim = [0, 1];
end
clim_J = double(used_clim);
auto_slice = double(gui_data.scaledautoVol(:, :, z));
nano_slice = double(gui_data.nanoVol(:, :, z));
bg_slice = logical(gui_data.bg_mask_vol(:, :, z));

% autofluorescence in red, nano in green, both with the same limits
red_ch = mat2gray(auto_slice, clim_J);
green_ch = mat2gray(nano_slice, clim_J);
blue_ch = ones(size(auto_slice));

% the background in grey (0.5 in every channel), the tissue as it is
mask_fg = ~bg_slice;
red_ch = red_ch .* double(mask_fg) + 0.5*(blue_ch .* double(bg_slice));
green_ch = green_ch .* double(mask_fg) + 0.5*(blue_ch .* double(bg_slice));
blue_ch = 0.5*(blue_ch .* double(bg_slice));
rgb = cat(3, red_ch, green_ch, blue_ch);

% show it
set(gui_data.imageHandle, 'CData', rgb);
axis(gui_data.imageAxes, 'image', 'off');

% remove the previous slice's polygons, and draw this slice's from their vertices
delete(findobj(gui_data.imageAxes, 'Type', 'images.roi.Polygon'));
gui_data.rois{z} = gobjects(length(gui_data.poly_vertices{z}), 1);
hold(gui_data.imageAxes, 'on');
for i = 1:length(gui_data.poly_vertices{z})
    verts = gui_data.poly_vertices{z}{i};
    if size(verts, 1) >= 3
        h = drawpolygon(gui_data.imageAxes, 'Position', verts, ...
            'Color', 'yellow', 'LineWidth', 2, 'FaceAlpha', 0.3);
        addlistener(h, 'ROIMoved', @(src, evt) on_roi_changed(fig, src, z));
        gui_data.rois{z}(i) = h;
    end
end
hold(gui_data.imageAxes, 'off');

% the title, with the keys
title_str = {sprintf('Slice %d/%d - Red: Scaled Autofluo, Green: Nano | Keys: <-/-> nav, d=draw new, del=delete last, s=save', ...
    z, gui_data.Z), ...
    'Edit polygons by dragging vertices. Right-click on polygon to delete. Do not press ESC during drawing.'};
set(gui_data.titleHandle, 'String', title_str);
guidata(fig, gui_data);
drawnow;
end

function callback_keypress(~, eventdata)
% Keys: the arrows move between slices, d draws a polygon, delete removes the
% last one, s saves, escape closes.

fig = gcbf;
gui_data = guidata(fig);
z = gui_data.current_z;
key = lower(eventdata.Key);
switch key
    case 'leftarrow'
        gui_data.current_z = max(1, z - 1);
    case 'rightarrow'
        gui_data.current_z = min(gui_data.Z, z + 1);
    case 'd'

        % draw a new polygon, with the title saying how
        normal_title = {sprintf('Slice %d/%d - Red: Scaled Autofluo, Green: Nano | Keys: <-/-> nav, d=draw new, del=delete last, s=save', ...
            z, gui_data.Z), ...
            'Edit polygons by dragging vertices. Right-click on polygon to delete. Do not press ESC during drawing.'};
        set(gui_data.titleHandle, 'String', ...
            'Drawing new polygon: Click points, double-click or right-click to close. Press ESC to cancel drawing only.');
        drawnow;

        % while drawing, a right click must not delete and escape must not close
        temp_window_button_down = get(fig, 'WindowButtonDownFcn');
        set(fig, 'WindowButtonDownFcn', '');
        temp_key_press = get(fig, 'KeyPressFcn');
        set(fig, 'KeyPressFcn', '');
        h = drawpolygon(gui_data.imageAxes, 'Color', 'yellow', 'LineWidth', 2, ...
            'FaceAlpha', 0.3);
        set(fig, 'KeyPressFcn', temp_key_press);
        set(fig, 'WindowButtonDownFcn', temp_window_button_down);

        % keep the polygon unless the drawing was cancelled
        if ~isempty(h) && isvalid(h)
            addlistener(h, 'ROIMoved', @(src, evt) on_roi_changed(fig, src, z));
            idx = length(gui_data.poly_vertices{z}) + 1;
            gui_data.poly_vertices{z}{idx} = h.Position;
            gui_data.rois{z}(end+1) = h;
            gui_data.artifact_mask(:, :, z) = compute_mask_from_rois(gui_data, z);
            guidata(fig, gui_data);
        end
        set(gui_data.titleHandle, 'String', normal_title);
        drawnow;
        return;
    case 'delete'

        % delete the slice's last polygon
        if ~isempty(gui_data.rois{z}) && isvalid(gui_data.rois{z}(end))
            delete(gui_data.rois{z}(end));
            gui_data.rois{z}(end) = [];
            if ~isempty(gui_data.poly_vertices{z})
                gui_data.poly_vertices{z}(end) = [];
            end
            gui_data.artifact_mask(:, :, z) = compute_mask_from_rois(gui_data, z);
            guidata(fig, gui_data);
            display_current_slice(fig);
        end
        return;
    case 's'
        save_artifact_masks(gui_data);
        return;
    case 'escape'
        callback_close_gui_request(fig, []);
        return;
end
guidata(fig, gui_data);
display_current_slice(fig);
end

function on_mouse_click(~, ~)
% A right click on a polygon deletes it (the last drawn, where polygons overlap).

sel_type = get(gcbf, 'SelectionType');

% 'alt' is a right click
if strcmp(sel_type, 'alt')
    fig = gcbf;
    gui_data = guidata(fig);
    z = gui_data.current_z;
    cp = get(gui_data.imageAxes, 'CurrentPoint');
    if ~isempty(cp)
        cp = cp(1, 1:2);
        for i = length(gui_data.rois{z}):-1:1
            h = gui_data.rois{z}(i);
            if isvalid(h) && inpolygon(cp(1), cp(2), h.Position(:, 1), h.Position(:, 2))
                delete(h);
                gui_data.rois{z}(i) = [];
                gui_data.poly_vertices{z}(i) = [];
                gui_data.artifact_mask(:, :, z) = compute_mask_from_rois(gui_data, z);
                guidata(fig, gui_data);
                display_current_slice(fig);
                break;
            end
        end
    end
end
end

function on_roi_changed(fig, src, z)
% Store a dragged polygon's new vertices and update the slice's mask.

gui_data = guidata(fig);
idx = find(gui_data.rois{z} == src, 1);
if ~isempty(idx)
    gui_data.poly_vertices{z}{idx} = src.Position;
    gui_data.artifact_mask(:, :, z) = compute_mask_from_rois(gui_data, z);
    guidata(fig, gui_data);
end
end

function mask = compute_mask_from_rois(gui_data, z)
% The mask of slice z: the union of its polygons of three vertices or more.

mask = false(gui_data.H, gui_data.W);
for i = 1:length(gui_data.poly_vertices{z})
    verts = gui_data.poly_vertices{z}{i};
    if size(verts, 1) >= 3
        pm = poly2mask(verts(:, 1), verts(:, 2), gui_data.H, gui_data.W);
        mask = mask | pm;
    end
end
end

function save_artifact_masks(gui_data)
% Save every slice's mask, recomputed from its polygons, and the polygons.

for z = 1:gui_data.Z
    gui_data.artifact_mask(:, :, z) = compute_mask_from_rois(gui_data, z);
end
mask_file = fullfile(gui_data.output_dir, sprintf('artifact_mask_volume_%s.mat', ...
    gui_data.correction_type));
artifact_mask_vol = gui_data.artifact_mask;
poly_vertices = gui_data.poly_vertices;
save(mask_file, 'artifact_mask_vol', 'poly_vertices', '-v7.3');
disp(sprintf('Saved artifact masks for %s in %s', gui_data.mouse_name, mask_file)); %#ok<DSPSP>
end

function callback_close_gui_request(src, ~)
% Ask whether to save, then close the window and let MATLAB go on, unless
% cancelled.

choice = questdlg('Save annotations before closing?', 'Confirm Close', 'Save & Close', ...
    'Discard & Close', 'Cancel', 'Save & Close');
if strcmp(choice, 'Save & Close')
    gui_data = guidata(src);
    save_artifact_masks(gui_data);
end
if strcmp(choice, 'Save & Close') || strcmp(choice, 'Discard & Close')
    uiresume(src);
    delete(src);
end
end
