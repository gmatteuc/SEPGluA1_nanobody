function ArtifactAnnotator(scaledautoVol, nanoVol, bg_mask_vol, slice_data, mouse_name, output_dir, correction_type)
% Interactively annotate artifacts by drawing polygons on RGB overlays (red: scaled autofluorescence, green: nanobody).
% Navigation: arrow keys; Draw: 'd' (add new polygon); Save: 's'; Close: ESC.
% Polygons are interactive and editable. Edits update mask automatically.
% Right-click on a polygon to delete it. 'delete' key deletes last polygon.
% Saves artifact_mask_vol as logical 3D array and poly_vertices for resuming.
[H, W, Z] = size(scaledautoVol);
mask_file = fullfile(output_dir, sprintf('artifact_mask_volume_%s.mat', correction_type));
% Initialize GUI data
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
% Load existing annotations if file exists
if exist(mask_file, 'file')
    loaded_data = load(mask_file);
    if isfield(loaded_data, 'artifact_mask_vol')
        gui_data.artifact_mask = loaded_data.artifact_mask_vol;
    end
    if isfield(loaded_data, 'poly_vertices')
        gui_data.poly_vertices = loaded_data.poly_vertices;
    else
        % If no poly_vertices, initialize empty (mask is non-editable but visible via recreation if needed)
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
% GUI setup
gui_fig = figure('Name', sprintf('Artifact Annotator - %s', mouse_name), 'NumberTitle', 'off', ...
    'Toolbar', 'none', 'Menubar', 'none', 'Color', 'w', ...
    'WindowState', 'maximized', ...
    'CloseRequestFcn', @(src,evt) callback_close_gui_request(src, evt), ...
    'KeyPressFcn', @callback_keypress, ...
    'Resize', 'on', ...
    'WindowButtonDownFcn', @on_mouse_click);
gui_data.imageAxes = axes('Parent', gui_fig, 'Units', 'normalized', 'Position', [0.05 0.05 0.9 0.9]);
gui_data.imageHandle = imshow(zeros(H,W,3), 'Parent', gui_data.imageAxes); % Ensure 3-channel
gui_data.titleHandle = title(gui_data.imageAxes, '', 'FontSize', 12);
guidata(gui_fig, gui_data);
% Display initial slice
display_current_slice(gui_fig);
if Z > 0
    uiwait(gui_fig);
else
    if ishandle(gui_fig), delete(gui_fig); end
    disp('No slices to annotate. GUI closed.');
    return;
end
end
function display_current_slice(fig)
gui_data = guidata(fig);
z = gui_data.current_z;
if z < 1 || z > gui_data.Z, return; end
% Get slice data
used_clim = gui_data.slice_data(z).used_clim;
if any(isnan(used_clim))
    used_clim = [0, 1]; % Fallback
end
clim_J = double(used_clim);
auto_slice = double(gui_data.scaledautoVol(:,:,z));
nano_slice = double(gui_data.nanoVol(:,:,z));
bg_slice = logical(gui_data.bg_mask_vol(:,:,z)); % Ensure logical for masking
% Normalize and create RGB (same clim for both)
red_ch = mat2gray(auto_slice, clim_J);
green_ch = mat2gray(nano_slice, clim_J);
blue_ch = ones(size(auto_slice));
% Apply mask to make background visible (black out background areas)
mask_fg = ~bg_slice; % Foreground mask (inverse of background)
red_ch = red_ch .* double(mask_fg) + 0.5*(blue_ch .* double(bg_slice));
green_ch = green_ch .* double(mask_fg) + 0.5*(blue_ch .* double(bg_slice));
blue_ch = 0.5*(blue_ch .* double(bg_slice));
rgb = cat(3, red_ch, green_ch, blue_ch);
% Display image
set(gui_data.imageHandle, 'CData', rgb);
axis(gui_data.imageAxes, 'image', 'off');
% Clear existing ROIs on the axes (previous slice's)
delete(findobj(gui_data.imageAxes, 'Type', 'images.roi.Polygon'));
% Re-create ROIs for current slice from stored vertices
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
% Update title
title_str = {sprintf('Slice %d/%d - Red: Scaled Autofluo, Green: Nano | Keys: <-/-> nav, d=draw new, del=delete last, s=save', z, gui_data.Z), ...
    'Edit polygons by dragging vertices. Right-click on polygon to delete. Do not press ESC during drawing.'};
set(gui_data.titleHandle, 'String', title_str);
guidata(fig, gui_data);
drawnow;
end
function callback_keypress(~, eventdata)
fig = gcbf; % Get current figure
gui_data = guidata(fig);
z = gui_data.current_z;
key = lower(eventdata.Key);
switch key
    case 'leftarrow'
        gui_data.current_z = max(1, z - 1);
    case 'rightarrow'
        gui_data.current_z = min(gui_data.Z, z + 1);
    case 'd'
        % Draw new polygon
        normal_title = {sprintf('Slice %d/%d - Red: Scaled Autofluo, Green: Nano | Keys: <-/-> nav, d=draw new, del=delete last, s=save', z, gui_data.Z), ...
            'Edit polygons by dragging vertices. Right-click on polygon to delete. Do not press ESC during drawing.'};
        set(gui_data.titleHandle, 'String', 'Drawing new polygon: Click points, double-click or right-click to close. Press ESC to cancel drawing only.');
        drawnow;
        temp_window_button_down = get(fig, 'WindowButtonDownFcn'); % Save original
        set(fig, 'WindowButtonDownFcn', ''); % Temporarily disable right-click delete during drawing
        temp_key_press = get(fig, 'KeyPressFcn');
        set(fig, 'KeyPressFcn', ''); % Temporarily disable keypress during drawing to prevent ESC close
        h = drawpolygon(gui_data.imageAxes, 'Color', 'yellow', 'LineWidth', 2, 'FaceAlpha', 0.3);
        set(fig, 'KeyPressFcn', temp_key_press); % Restore keypress
        set(fig, 'WindowButtonDownFcn', temp_window_button_down); % Restore
        if ~isempty(h) && isvalid(h)
            addlistener(h, 'ROIMoved', @(src, evt) on_roi_changed(fig, src, z));
            idx = length(gui_data.poly_vertices{z}) + 1;
            gui_data.poly_vertices{z}{idx} = h.Position;
            gui_data.rois{z}(end+1) = h;
            gui_data.artifact_mask(:,:,z) = compute_mask_from_rois(gui_data, z);
            guidata(fig, gui_data);
        end
        set(gui_data.titleHandle, 'String', normal_title);
        drawnow;
        return;
    case 'delete'
        % Delete last ROI on current slice
        if ~isempty(gui_data.rois{z}) && isvalid(gui_data.rois{z}(end))
            delete(gui_data.rois{z}(end));
            gui_data.rois{z}(end) = [];
            if ~isempty(gui_data.poly_vertices{z})
                gui_data.poly_vertices{z}(end) = [];
            end
            gui_data.artifact_mask(:,:,z) = compute_mask_from_rois(gui_data, z);
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
% For right-click to delete ROI
sel_type = get(gcbf, 'SelectionType');
if strcmp(sel_type, 'alt') % Right-click
    fig = gcbf;
    gui_data = guidata(fig);
    z = gui_data.current_z;
    cp = get(gui_data.imageAxes, 'CurrentPoint');
    if ~isempty(cp)
        cp = cp(1,1:2);
        for i = length(gui_data.rois{z}):-1:1
            h = gui_data.rois{z}(i);
            if isvalid(h) && inpolygon(cp(1), cp(2), h.Position(:,1), h.Position(:,2))
                delete(h);
                gui_data.rois{z}(i) = [];
                gui_data.poly_vertices{z}(i) = [];
                gui_data.artifact_mask(:,:,z) = compute_mask_from_rois(gui_data, z);
                guidata(fig, gui_data);
                display_current_slice(fig);
                break;
            end
        end
    end
end
end
function on_roi_changed(fig, src, z)
gui_data = guidata(fig);
idx = find(gui_data.rois{z} == src, 1);
if ~isempty(idx)
    gui_data.poly_vertices{z}{idx} = src.Position;
    gui_data.artifact_mask(:,:,z) = compute_mask_from_rois(gui_data, z);
    guidata(fig, gui_data);
end
end
function mask = compute_mask_from_rois(gui_data, z)
mask = false(gui_data.H, gui_data.W);
for i = 1:length(gui_data.poly_vertices{z})
    verts = gui_data.poly_vertices{z}{i};
    if size(verts, 1) >= 3
        pm = poly2mask(verts(:,1), verts(:,2), gui_data.H, gui_data.W);
        mask = mask | pm;
    end
end
end
function save_artifact_masks(gui_data)
for z = 1:gui_data.Z
    gui_data.artifact_mask(:,:,z) = compute_mask_from_rois(gui_data, z);
end
mask_file = fullfile(gui_data.output_dir, sprintf('artifact_mask_volume_%s.mat', gui_data.correction_type));
artifact_mask_vol = gui_data.artifact_mask;
poly_vertices = gui_data.poly_vertices;
save(mask_file, 'artifact_mask_vol', 'poly_vertices', '-v7.3');
disp(sprintf('Saved artifact masks for %s in %s', gui_data.mouse_name, mask_file)); %#ok<DSPSP>
end
function callback_close_gui_request(src, ~)
choice = questdlg('Save annotations before closing?', 'Confirm Close', 'Save & Close', 'Discard & Close', 'Cancel', 'Save & Close');
if strcmp(choice, 'Save & Close')
    gui_data = guidata(src);
    save_artifact_masks(gui_data);
end
if strcmp(choice, 'Save & Close') || strcmp(choice, 'Discard & Close')
    uiresume(src);
    delete(src);
end
end