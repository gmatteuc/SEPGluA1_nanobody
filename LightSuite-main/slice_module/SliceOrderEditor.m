function SliceOrderEditor(optionalVolumePath)
% Interactively reorders slices using a direct "Move to Position" workflow.
% This is the final, redesigned version for clarity and ease of use.

    % --- Initial Setup & GUI Creation (No changes here) ---
    inputFileFullPath = '';
    if nargin > 0 && ~isempty(optionalVolumePath) && exist(optionalVolumePath, 'file')
        [~, ~, ext] = fileparts(optionalVolumePath);
        if strcmpi(ext, '.tif') || strcmpi(ext, '.tiff')
            inputFileFullPath = optionalVolumePath;
        else
            warning('Provided path is not a TIFF file. Please select a file manually.');
        end
    end
    if isempty(inputFileFullPath)
        [fileName, pathName] = uigetfile({'*.tif;*.tiff', 'TIFF Files (*.tif, *.tiff)'}, ...
                                         'Select Multi-Page TIFF File');
        if isequal(fileName, 0) || isequal(pathName, 0)
            disp('No file selected. Exiting.');
            return;
        end
        inputFileFullPath = fullfile(pathName, fileName);
    end
    try
        tiffInfo = imfinfo(inputFileFullPath);
        numSlices = numel(tiffInfo);
        if numSlices == 0, errordlg('The selected TIFF file contains no images.', 'File Error'); return; end
    catch ME
        errordlg(['Error reading TIFF file info: ' ME.message], 'File Error');
        return;
    end
    gui_data.originalSliceImages = cell(numSlices, 1);
    hWaitBar = waitbar(0, 'Loading slices...');
    try
        for i = 1:numSlices
            waitbar(i/numSlices, hWaitBar);
            gui_data.originalSliceImages{i} = imread(inputFileFullPath, i, 'Info', tiffInfo);
        end
    catch ME
        if ishandle(hWaitBar); close(hWaitBar); end
        errordlg(['Error loading slices from TIFF: ' ME.message], 'Image Loading Error');
        return;
    end
    if ishandle(hWaitBar); close(hWaitBar); end
    [filePathStr, baseName, ~] = fileparts(inputFileFullPath);
    gui_data.processingDecisionsFilename = fullfile(filePathStr, strcat(baseName, '_processing_decisions.txt'));
    screenSize = get(0, 'ScreenSize');
    gui_aspect_ratio = 1.6; gui_width_fraction = 0.5;
    gui_width_px = screenSize(3) * gui_width_fraction;
    gui_position = [(screenSize(3)-gui_width_px)/2, (screenSize(4)-gui_width_px/gui_aspect_ratio)/2, gui_width_px, gui_width_px/gui_aspect_ratio];
    gui_fig = figure('Name', 'Slice Order Editor - Close-up', 'NumberTitle', 'off', 'Toolbar', 'none', 'Menubar', 'none', 'Color', 'w', 'Units', 'pixels', 'Position', gui_position, 'CloseRequestFcn', @(src,evt) callback_close_gui_request(src), 'KeyPressFcn', @callback_keypress);
    gui_data.numSlices = numSlices;
    gui_data.currentDisplayPosition = 1;
    gui_data.flipState = zeros(numSlices, 1);
    gui_data.displaySequenceOriginalIndices = (1:numSlices)';
    % Montage companion window. Everything it needs lives under these fields, and
    % all of its code sits in the block at the bottom of this file, so it can be
    % removed again by deleting that block plus the few calls marked "montage".
    gui_data.showMontage      = true;   % set false for the original single-window editor
    gui_data.montageColumns   = 8;
    gui_data.montageThumbWidth = 150;   % pixels; drives how big the montage window is
    gui_data.montageFig       = [];
    gui_data.montageAxes      = [];
    gui_data.montageImage     = [];
    gui_data.montageHighlight = [];
    gui_data.montageTarget    = [];     % dashed frame shown while dragging
    gui_data.montageLabels    = [];
    gui_data.montageMarks     = [];     % red crosses on excluded slices
    gui_data.montageThumbs    = {};
    gui_data.dragFromPosition = [];
    load_processing_decisions();
    gui_data.imageAxes = axes('Parent', gui_fig, 'Units', 'normalized', 'Position', [0.05 0.05 0.9 0.9]);
    gui_data.imageAxes.Colormap = colormap('gray'); axis(gui_data.imageAxes, 'image', 'off'); 
    gui_data.imageHandle = image(gui_data.imageAxes, []); 
    gui_data.titleHandle = title(gui_data.imageAxes, '', 'FontSize', 10);
    gui_data.excludeMarkerHandle = [];
    gui_data.orderTextHandle = text(gui_data.imageAxes, 0, 0, '', 'FontSize', 24, 'Color', 'yellow', 'FontWeight', 'bold', 'BackgroundColor', [0 0 0 0.5], 'VerticalAlignment', 'top');
    guidata(gui_fig, gui_data);
    display_current_slice(gui_fig);
    build_montage_window(gui_fig);   % montage
    if gui_data.numSlices > 0, uiwait(gui_fig);
    elseif ishandle(gui_fig), delete(gui_fig); disp('No slices to display. GUI closed.'); end

    function load_processing_decisions()
        if exist(gui_data.processingDecisionsFilename, 'file') && gui_data.numSlices > 0
            try
                loadedTable = readtable(gui_data.processingDecisionsFilename, 'Delimiter', '\t', 'ReadVariableNames', true);
                if height(loadedTable) ~= gui_data.numSlices || ~all(ismember({'OriginalIndex', 'FlipState', 'NewOrderOriginalIndex'}, loadedTable.Properties.VariableNames))
                    warning('Decision file has mismatching dimensions or columns. Using defaults.'); return;
                end
                flipData = loadedTable.FlipState;
                if isnumeric(flipData)
                    gui_data.flipState = flipData;
                else 
                    tempFlipState = zeros(gui_data.numSlices,1);
                    for k_f = 1:gui_data.numSlices
                        val_str = strtrim(string(flipData{k_f}));
                        if strcmpi(val_str, "1"), tempFlipState(k_f) = 1;
                        elseif strcmpi(val_str, "-1"), tempFlipState(k_f) = -1;
                        end
                    end
                    gui_data.flipState = tempFlipState;
                end
                loadedOrder = loadedTable.NewOrderOriginalIndex;
                isValidPermutation = ~any(isnan(loadedOrder)) && numel(unique(loadedOrder)) == gui_data.numSlices;
                if isValidPermutation
                    gui_data.displaySequenceOriginalIndices = loadedOrder;
                    disp('Loaded saved order from decisions file.');
                else
                    disp('Saved order is incomplete or invalid. Using default 1-N order.');
                end
            catch ME_load
                warning('Error loading/parsing decisions file: %s. Using defaults.', ME_load.message);
            end
        end
    end
end

function display_current_slice(fig)
    gui_data = guidata(fig);
    if gui_data.numSlices == 0; return; end
    if ishandle(gui_data.excludeMarkerHandle); delete(gui_data.excludeMarkerHandle); gui_data.excludeMarkerHandle = []; end
    originalSliceIdx = gui_data.displaySequenceOriginalIndices(gui_data.currentDisplayPosition);
    imgData = gui_data.originalSliceImages{originalSliceIdx};
    currentFlipState = gui_data.flipState(originalSliceIdx);
    status_str = 'Normal';
    if currentFlipState == 1, imgData = fliplr(imgData); status_str = 'FLIPPED';
    elseif currentFlipState == -1, status_str = 'EXCLUDED'; end
    set(gui_data.imageHandle, 'CData', imgData); axis(gui_data.imageAxes, 'image', 'off');   
    if currentFlipState == -1
        [h, w, ~] = size(imgData);
        gui_data.excludeMarkerHandle = text(gui_data.imageAxes, w/2, h/2, 'X', 'FontSize', h/4, 'FontWeight', 'bold', 'Color', 'r', 'HorizontalAlignment', 'center');
    end
    orderNumStr = sprintf('%d', gui_data.currentDisplayPosition);
    set(gui_data.orderTextHandle, 'Position', [0.02*size(imgData,2), 0.02*size(imgData,1)], 'String', orderNumStr);
    title_str = {sprintf('Slice at Order Position: %d/%d (Original Index: %d) - %s', ...
                         gui_data.currentDisplayPosition, gui_data.numSlices, originalSliceIdx, status_str), ...
                 'Keys: Nav = Left/Right | Flip = f | Exclude = o | Reorder = Enter | Save = s | Quit = Esc'};
    set(gui_data.titleHandle, 'String', title_str);
    guidata(fig, gui_data);
    update_montage_highlight(fig);   % montage
end

function callback_keypress(fig, eventdata)
    gui_data = guidata(fig);
    if gui_data.numSlices == 0; return; end
    originalSliceIdx = gui_data.displaySequenceOriginalIndices(gui_data.currentDisplayPosition);
    refresh_tiles = false;   % montage: flipping or excluding changes how a tile looks
    switch eventdata.Key
        case 'leftarrow'
            gui_data.currentDisplayPosition = max(1, gui_data.currentDisplayPosition - 1);
        case 'rightarrow'
            gui_data.currentDisplayPosition = min(gui_data.numSlices, gui_data.currentDisplayPosition + 1);
        case 'f'
            if gui_data.flipState(originalSliceIdx) ~= -1, gui_data.flipState(originalSliceIdx) = ~gui_data.flipState(originalSliceIdx); end
            refresh_tiles = true;
        case 'o'
            if gui_data.flipState(originalSliceIdx) == -1, gui_data.flipState(originalSliceIdx) = 0;
            else, gui_data.flipState(originalSliceIdx) = -1; end
            refresh_tiles = true;
        case {'return', 'enter'}
            reorder_slice_callback(fig); return;
        case 's', guidata(fig, gui_data); save_processing_decisions(gui_data);return;
        case 'escape', callback_close_gui_request(fig); return;
    end
    guidata(fig, gui_data);
    display_current_slice(fig);
    if refresh_tiles, refresh_montage(fig); end   % montage
end

function reorder_slice_callback(fig)
    gui_data = guidata(fig);
    current_pos = gui_data.currentDisplayPosition;
    slice_to_move_idx = gui_data.displaySequenceOriginalIndices(current_pos);
    
    % --- REDESIGNED WORKFLOW ---
    prompt = {sprintf('Move slice (Orig. Idx: %d) to new order position (1-%d):', slice_to_move_idx, gui_data.numSlices)};
    dlg_title = 'Move to Position';
    answer = inputdlg(prompt, dlg_title, [1 60], {num2str(current_pos)});
    
    if isempty(answer), return; end 
    
    new_pos = round(str2double(answer{1}));
    
    if isnan(new_pos) || new_pos < 1 || new_pos > gui_data.numSlices
        warndlg('Invalid input. Please enter a valid position number.', 'Input Error');
        return;
    end
    
    % --- Call the new, simpler, direct logic ---
    new_order_vec = perform_move_to_position(gui_data.displaySequenceOriginalIndices, slice_to_move_idx, new_pos);
    
    gui_data.displaySequenceOriginalIndices = new_order_vec;
    gui_data.currentDisplayPosition = new_pos; % The new position is exactly what the user entered

    guidata(fig, gui_data);
    display_current_slice(fig);
    refresh_montage(fig);   % montage
end

function new_order_vec = perform_move_to_position(order_vec, slice_to_move_idx, new_pos)
    % This is the new, simpler logic: remove a slice and insert it at a specific index.
    
    % Remove the slice from the list
    temp_order_vec = order_vec(order_vec ~= slice_to_move_idx);
    
    if new_pos == 1
        % Insert at the beginning
        new_order_vec = [slice_to_move_idx; temp_order_vec];
    elseif new_pos == numel(order_vec)
        % Insert at the end
        new_order_vec = [temp_order_vec; slice_to_move_idx];
    else
        % Insert in the middle
        new_order_vec = [temp_order_vec(1:new_pos-1); ...
                         slice_to_move_idx; ...
                         temp_order_vec(new_pos:end)];
    end
    disp(sprintf('Moved slice %d to position %d.', slice_to_move_idx, new_pos));
end

function callback_close_gui_request(fig)
    choice = questdlg('Save changes to order and states before closing?', 'Confirm Close', 'Save & Close', 'Discard & Close', 'Cancel', 'Save & Close');
    if strcmp(choice, 'Save & Close'), save_and_close_confirmed(fig);
    elseif strcmp(choice, 'Discard & Close')
        close_montage_window(fig);   % montage
        if ishandle(fig); delete(fig); end
    end
end

function save_and_close_confirmed(fig)
    gui_data = guidata(fig);
    try
        save_processing_decisions(gui_data);
        disp('Order and state decisions saved.');
    catch ME
        errordlg(['Error saving decisions: ', ME.message], 'Save Error');
        return;
    end
    close_montage_window(fig);   % montage
    if ishandle(fig); delete(fig); end
end

function save_processing_decisions(gui_data)
    decisionsMatrix = zeros(gui_data.numSlices, 3);
    decisionsMatrix(:,1) = (1:gui_data.numSlices)'; 
    decisionsMatrix(:,2) = gui_data.flipState;
    decisionsMatrix(:,3) = gui_data.displaySequenceOriginalIndices;
    try
        T = array2table(decisionsMatrix, 'VariableNames', {'OriginalIndex', 'FlipState', 'NewOrderOriginalIndex'});
        writetable(T, gui_data.processingDecisionsFilename, 'WriteVariableNames', true, 'Delimiter', '\t');
        disp(['Processing decisions saved to: ', gui_data.processingDecisionsFilename]);
    catch ME_write
        rethrow(ME_write);
    end
end

% =========================================================================
%  Montage window
% -------------------------------------------------------------------------
%  A companion window showing every slice at once, in the current order. You
%  can click a tile to jump the close-up window to it, or drag a tile onto
%  another position to reorder. Arrow keys work from either window and the
%  selection stays in step.
%
%  Everything below is self contained: it only reads and writes the montage*
%  fields of gui_data, plus displaySequenceOriginalIndices and flipState which
%  the rest of the editor already owns. Delete this block and the handful of
%  calls marked "% montage" to get the original single-window editor back.
% =========================================================================

function build_montage_window(fig)
gui_data = guidata(fig);
if ~gui_data.showMontage || gui_data.numSlices == 0, return; end

% Thumbnails are made once. The slices themselves never change, only their
% order and their flip state, so there is no reason to rescale them again on
% every redraw.
gui_data.montageThumbs = make_thumbnails(gui_data.originalSliceImages, gui_data.montageThumbWidth);

[th, tw, ~] = size(gui_data.montageThumbs{1});
n_cols = gui_data.montageColumns;
n_rows = ceil(gui_data.numSlices / n_cols);

% Size the window to the montage, but keep it on screen
screen = get(0, 'ScreenSize');
win_w = min(n_cols * tw + 40, screen(3) * 0.9);
win_h = min(n_rows * th + 60, screen(4) * 0.85);

gui_data.montageFig = figure( ...
    'Name', 'Slice Order Editor - Montage', 'NumberTitle', 'off', ...
    'Toolbar', 'none', 'Menubar', 'none', 'Color', 'k', 'Units', 'pixels', ...
    'Position', [40, max(60, screen(4) - win_h - 100), win_w, win_h], ...
    'CloseRequestFcn', @(src, evt) montage_close_request(src), ...
    'KeyPressFcn', @(src, evt) forward_key_to_main(fig, evt), ...
    'WindowButtonDownFcn', @(src, evt) montage_button_down(fig));

gui_data.montageAxes = axes('Parent', gui_data.montageFig, ...
    'Units', 'normalized', 'Position', [0 0 1 1], 'Color', 'k');
gui_data.montageImage = image(gui_data.montageAxes, zeros(n_rows*th, n_cols*tw, 3, 'uint8'));
axis(gui_data.montageAxes, 'image', 'off');
hold(gui_data.montageAxes, 'on');

% Frame around the slice currently open in the close-up window
gui_data.montageHighlight = rectangle('Parent', gui_data.montageAxes, ...
    'Position', [0.5 0.5 tw th], 'EdgeColor', [1 1 0], 'LineWidth', 3);

% Dashed frame showing where a dragged slice would land, hidden until needed
gui_data.montageTarget = rectangle('Parent', gui_data.montageAxes, ...
    'Position', [0.5 0.5 tw th], 'EdgeColor', [0 1 1], 'LineWidth', 2, ...
    'LineStyle', '--', 'Visible', 'off');

% One label per tile, reused rather than recreated on every redraw
gui_data.montageLabels = gobjects(gui_data.numSlices, 1);
for k = 1:gui_data.numSlices
    [x0, y0] = tile_corner(k, n_cols, th, tw);
    gui_data.montageLabels(k) = text(gui_data.montageAxes, x0 + 4, y0 + 10, '', ...
        'Color', [1 1 0], 'FontSize', 9, 'FontWeight', 'bold', ...
        'VerticalAlignment', 'top', 'Interpreter', 'none', 'HitTest', 'off');
end

guidata(fig, gui_data);
refresh_montage(fig);

% Keep the keyboard on the close-up window, which is where most work happens
if ishandle(fig), figure(fig); end
end


function refresh_montage(fig)
% Rebuild the composite from the current order and flip states.
gui_data = guidata(fig);
if ~gui_data.showMontage || isempty(gui_data.montageFig) || ~ishandle(gui_data.montageFig)
    return
end

[th, tw, ~] = size(gui_data.montageThumbs{1});
n_cols = gui_data.montageColumns;
n_rows = ceil(gui_data.numSlices / n_cols);

canvas = zeros(n_rows * th, n_cols * tw, 3, 'uint8');

% Red crosses are drawn per excluded slice, so clear the previous ones
delete(gui_data.montageMarks(ishandle(gui_data.montageMarks)));
gui_data.montageMarks = gobjects(0);

for pos = 1:gui_data.numSlices
    orig_idx = gui_data.displaySequenceOriginalIndices(pos);
    thumb = gui_data.montageThumbs{orig_idx};
    state = gui_data.flipState(orig_idx);

    if state == 1
        thumb = fliplr(thumb);
    elseif state == -1
        thumb = uint8(single(thumb) * 0.35);   % dim what is marked for removal
    end

    [x0, y0] = tile_corner(pos, n_cols, th, tw);
    canvas(y0:(y0+th-1), x0:(x0+tw-1), :) = thumb;

    if state == -1
        gui_data.montageMarks(end+1) = text(gui_data.montageAxes, ...
            x0 + tw/2, y0 + th/2, 'X', 'Color', [1 0.2 0.2], 'FontSize', 22, ...
            'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', 'HitTest', 'off'); %#ok<AGROW>
    end

    label = sprintf('%d [%d]', pos, orig_idx);
    if state == 1, label = [label ' F']; end
    set(gui_data.montageLabels(pos), 'String', label);
end

set(gui_data.montageImage, 'CData', canvas);
guidata(fig, gui_data);
update_montage_highlight(fig);
end


function update_montage_highlight(fig)
% Move the yellow frame to whichever slice the close-up window is showing.
gui_data = guidata(fig);
if ~gui_data.showMontage || isempty(gui_data.montageHighlight) || ~ishandle(gui_data.montageHighlight)
    return
end

[th, tw, ~] = size(gui_data.montageThumbs{1});
[x0, y0] = tile_corner(gui_data.currentDisplayPosition, gui_data.montageColumns, th, tw);
set(gui_data.montageHighlight, 'Position', [x0 - 0.5, y0 - 0.5, tw, th]);
end


function montage_button_down(fig)
% A press starts a possible drag. Whether it turns out to be a click or a drag
% is decided on release, by comparing where the mouse went up with where it
% went down.
gui_data = guidata(fig);
if gui_data.numSlices == 0, return; end

pos = position_under_cursor(gui_data);
if isempty(pos), return; end

gui_data.dragFromPosition = pos;
guidata(fig, gui_data);

set(gui_data.montageFig, ...
    'WindowButtonMotionFcn', @(src, evt) montage_drag(fig), ...
    'WindowButtonUpFcn',     @(src, evt) montage_release(fig));
end


function montage_drag(fig)
gui_data = guidata(fig);
if isempty(gui_data.dragFromPosition), return; end

pos = position_under_cursor(gui_data);
if isempty(pos) || pos == gui_data.dragFromPosition
    set(gui_data.montageTarget, 'Visible', 'off');
    return
end

[th, tw, ~] = size(gui_data.montageThumbs{1});
[x0, y0] = tile_corner(pos, gui_data.montageColumns, th, tw);
set(gui_data.montageTarget, 'Position', [x0 - 0.5, y0 - 0.5, tw, th], 'Visible', 'on');
end


function montage_release(fig)
gui_data = guidata(fig);
from_pos = gui_data.dragFromPosition;
if isempty(from_pos), return; end

set(gui_data.montageFig, 'WindowButtonMotionFcn', '', 'WindowButtonUpFcn', '');
set(gui_data.montageTarget, 'Visible', 'off');
gui_data.dragFromPosition = [];
guidata(fig, gui_data);

to_pos = position_under_cursor(gui_data);

if isempty(to_pos) || to_pos == from_pos
    % Released where it started: treat as a click and just jump there
    gui_data.currentDisplayPosition = from_pos;
    guidata(fig, gui_data);
    display_current_slice(fig);
    return
end

% Dragged onto another tile: move the slice to that position
slice_to_move = gui_data.displaySequenceOriginalIndices(from_pos);
gui_data.displaySequenceOriginalIndices = ...
    perform_move_to_position(gui_data.displaySequenceOriginalIndices, slice_to_move, to_pos);
gui_data.currentDisplayPosition = to_pos;
guidata(fig, gui_data);

display_current_slice(fig);
refresh_montage(fig);
end


function pos = position_under_cursor(gui_data)
% Which order position the mouse is over, or empty if it is off the grid.
pos = [];
if isempty(gui_data.montageAxes) || ~ishandle(gui_data.montageAxes), return; end

cp = get(gui_data.montageAxes, 'CurrentPoint');
x = cp(1, 1);
y = cp(1, 2);

[th, tw, ~] = size(gui_data.montageThumbs{1});
col = floor((x - 0.5) / tw);
row = floor((y - 0.5) / th);
if col < 0 || col >= gui_data.montageColumns || row < 0, return; end

candidate = row * gui_data.montageColumns + col + 1;
if candidate >= 1 && candidate <= gui_data.numSlices
    pos = candidate;
end
end


function [x0, y0] = tile_corner(pos, n_cols, th, tw)
% Top-left pixel of the tile at a given order position, 1-based.
row = floor((pos - 1) / n_cols);
col = mod(pos - 1, n_cols);
x0 = col * tw + 1;
y0 = row * th + 1;
end


function thumbs = make_thumbnails(images, target_width)
% Scale every slice down to a common tile size. Slices can differ slightly in
% size, so each is resized to the same width and then padded to the tallest
% height, which keeps the montage on a clean grid.
n = numel(images);
scaled = cell(n, 1);

for k = 1:n
    img = images{k};
    if size(img, 3) == 1
        img = repmat(img, 1, 1, 3);   % the editor also accepts grayscale stacks
    end
    scale = target_width / size(img, 2);
    scaled{k} = imresize(img, scale);
end

max_h = max(cellfun(@(x) size(x, 1), scaled));
thumbs = cell(n, 1);
for k = 1:n
    pad = max_h - size(scaled{k}, 1);
    thumbs{k} = padarray(scaled{k}, [pad 0 0], 0, 'post');
end
end


function forward_key_to_main(fig, evt)
% Arrow keys and shortcuts should work whichever window has focus.
if ishandle(fig)
    callback_keypress(fig, evt);
end
end


function montage_close_request(montage_fig)
% Closing the montage on its own just hides that view; the close-up window
% stays in charge of saving and quitting.
if ishandle(montage_fig), delete(montage_fig); end
end


function close_montage_window(fig)
if ~ishandle(fig), return; end
gui_data = guidata(fig);
if isfield(gui_data, 'montageFig') && ~isempty(gui_data.montageFig) && ishandle(gui_data.montageFig)
    delete(gui_data.montageFig);
end
end