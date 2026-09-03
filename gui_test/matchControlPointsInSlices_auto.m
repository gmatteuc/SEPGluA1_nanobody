function matchControlPointsInSlices_auto(opts)
% Manually align histology slices and matched CCF slices

% Initialize guidata
gui_data = struct;
opts.downfac_reg = opts.allenres/opts.registres;

% Load atlas
allen_atlas_path = fileparts(which('average_template_10.nii.gz'));
if isempty(allen_atlas_path)
    error('No CCF atlas found (add CCF atlas to path)')
end
disp('Loading Allen CCF atlas...')
gui_data.tv      = niftiread(fullfile(allen_atlas_path,'average_template_10.nii.gz'));
factv            = 255/single(max(gui_data.tv,[],"all"));
gui_data.tv      = uint8(single(gui_data.tv(opts.atlasaplims(1):opts.atlasaplims(2), :, :))*factv);
gui_data.tv      = imresize3(gui_data.tv,opts.downfac_reg);
gui_data.av      = niftiread(fullfile(allen_atlas_path,'annotation_10.nii.gz'));
gui_data.av      = imresize3(gui_data.av(opts.atlasaplims(1):opts.atlasaplims(2), :, :),...
    opts.downfac_reg, "Method","nearest");
disp('Done.')

gui_data.save_path = opts.procpath;

% --- Load sample volume ---
volume_dir       = dir(fullfile(opts.procpath,'*inspection.tif*'));

volpath          = fullfile(volume_dir.folder, volume_dir.name);
volload          = readDownStack(volpath);
volload          = permute(volload, [1 2 4 3]);
volload          = permute(volload, [opts.howtoperm 4]);
volload          = single(volload);
minvals          = single(quantile(volload, 0.01, [2 3]));
maxvals          = single(quantile(volload, 0.999, [2 3]));
volload          = 255 * (volload - minvals)./(maxvals - minvals);
gui_data.volume  = uint8(volload);
gui_data.Nslices = size(gui_data.volume, 1);
gui_data.colsuse    = 1:size(volload, 4);

nfac    = ceil(opts.extentfactor * 7.5/opts.pxsizes(1));

% --- Setup spatial referencing ---
Ratlas  = imref3d(size(gui_data.tv));
Rvolume = imref3d(size(gui_data.volume, 1:3), 1, opts.pxsizes(1), 1);
yworld  = [Rvolume.YWorldLimits(1)-opts.pxsizes(1)*nfac, Rvolume.YWorldLimits(2)+nfac*opts.pxsizes(1)];
ypix    = ceil(range(yworld));
Rout    = imref3d([ypix, size(gui_data.tv, [2 3])], Rvolume.XWorldLimits, yworld,Rvolume.ZWorldLimits);
tformuse =  opts.tformrigid_allen_to_samp_20um;
%-------------------------------------------------------------------------
% Define path to the saved cutting angle data
angle_data_path = fullfile(opts.procpath, 'cutting_angle_data.mat');

% Check if the cutting angle data file exists
if exist(angle_data_path, 'file')
    fprintf('Found saved cutting angle data. Overwriting transformation rotation...\n');
    
    % --- 1. Load and process angle data ---
    data = load(angle_data_path);
    cutting_angle_data = data.cutting_angle_data;
    if ~isfield(cutting_angle_data, 'saved_slice_vectors')
        error('The GUI must be updated to save "saved_slice_vectors" for this to work.');
    end
    data = load(angle_data_path);
    cutting_angle_data = data.cutting_angle_data;
    tformuse = applyAngleToTransform(tformuse, Ratlas, cutting_angle_data);    
end
gui_data.tv = imwarp(gui_data.tv, Ratlas, tformuse, 'linear',  'OutputView', Rout);
gui_data.av = imwarp(gui_data.av, Ratlas, tformuse, 'nearest', 'OutputView', Rout);
%-------------------------------------------------------------------------
% This section remains the same
yatlasvals     = linspace(yworld(1), yworld(2), ypix + 1);
yatlasvals     = yatlasvals(1:end-1) + median(diff(yatlasvals))/2;
ysamplevals    = linspace(Rvolume.YWorldLimits(1), Rvolume.YWorldLimits(2), gui_data.Nslices+1);
ysamplevals    = ysamplevals(1:end-1) + median(diff(ysamplevals))/2;
[~, atlasinds] = min(pdist2(ysamplevals',yatlasvals'), [],2);
gui_data.atlasinds   = atlasinds;
gui_data.atlasindsuse  = atlasinds;
gui_data.slicewidth = median(diff(yatlasvals))* opts.pxsizes(1);
gui_data.atlasvals  = yatlasvals;
gui_data.Rmoving    = imref2d(size(gui_data.av, [2 3]));
gui_data.Rfixed     = imref2d(size(volload, [2 3]));

%--------------------------------------------------------------------------
[rows, columns] = size(gui_data.tv, [2 3]);
linespacing     = round(min([rows, columns]/12));

xspaces = 1:linespacing:columns;
xlinesx = [repmat(xspaces, [2 1]); nan(1, numel(xspaces))];
xlinesy = [repmat([1; rows], [1 numel(xspaces)]) ; nan(1, numel(xspaces))];

yspaces    = 1:linespacing:rows;
ylinesy   = [repmat(yspaces, [2 1]); nan(1, numel(yspaces))];
ylinesx   = [repmat([1; columns], [1 numel(yspaces)]) ; nan(1, numel(yspaces))];
xlinesall = [xlinesx(:); ylinesx(:)];
ylinesall = [xlinesy(:); ylinesy(:)];
%--------------------------------------------------------------------------
% Load automated alignment
auto_ccf_alignment_fn = fullfile(gui_data.save_path,'atlas2histology_tform.mat');
if exist(auto_ccf_alignment_fn,'file')
    oldtform = load(auto_ccf_alignment_fn);
    % gui_data.histology_ccf_auto_alignment = oldtform.atlas2histology_tform;
    gui_data.histology_control_points     = oldtform.histology_control_points;
    gui_data.atlas_control_points         = oldtform.atlas_control_points;
else
    % Initialize alignment control points and tform matricies
    gui_data.histology_control_points = repmat({zeros(0,3)}, gui_data.Nslices, 1);
    gui_data.atlas_control_points     = repmat({zeros(0,3)}, gui_data.Nslices,1);
end

% Editing state. Points are paired by ROW INDEX between the two lists -- row i
% of the atlas points goes with row i of the histology points in fitgeotform2d
% -- so anything that removes a point has to remove the same row from both, or
% every later pair silently shifts onto the wrong partner.
gui_data.edit_mode     = false;    % 'e' toggles; click selects instead of adding
gui_data.carry_forward = false;    % 'p' toggles; new slice inherits the last one's points
gui_data.sel_side      = '';       % 'histology' | 'atlas' | '' when nothing selected
gui_data.sel_idx       = 0;
gui_data.select_radius_frac = 0.04;   % of image width, for click-to-select

% Points carried onto a slice are PROVISIONAL until touched: drawn hollow in a
% different colour, and they do not pin the atlas plane, so the wheel still
% scrolls freely and they follow the view. The moment one is added, dragged or
% deleted the slice commits, the normal colours come back, and the plane locks
% the way it always did -- that lock is a safety feature, not an obstacle, and
% it should still protect anything actually placed by hand.
gui_data.provisional = false(gui_data.Nslices, 1);
gui_data.order_problem = '';   % set by update_atlas_slice, painted purple below

% Landmark proposals (landmark_refine.m). Nothing runs on its own: you find
% the best-matching atlas plane with the wheel, then press r, and the
% neighbouring slice's points are refined through the image matcher against
% THAT plane and placed here as provisional points. A plain copy was measured
% at 15 px from the final answer -- no better than starting empty -- the
% refined one at 10. Matching uses the autofluorescence channel, the closest
% modality to the template.
gui_data.refine_chan  = min(3, size(gui_data.volume, 4));
gui_data.hist_labels  = gobjects(0);   % the little index numbers next to each point
gui_data.atlas_labels = gobjects(0);

% Create figure, set button functions
screen_size_px = get(0,'screensize');
gui_aspect_ratio = 1.7; % width/length
gui_width_fraction = 0.6; % fraction of screen width to occupy
gui_width_px = screen_size_px(3).*gui_width_fraction;
gui_position = [...
    (screen_size_px(3)-gui_width_px)/2, ... % left x
    (screen_size_px(4)-gui_width_px/gui_aspect_ratio)/2, ... % bottom y
    gui_width_px,gui_width_px/gui_aspect_ratio]; % width, height

% The mouse this session belongs to, for the window title. procpath ends in
% \lightsuite, so the folder above it is the mouse.
[mouse_path, ~] = fileparts(opts.procpath);
[~, gui_data.mouse_name] = fileparts(mouse_path);

gui_fig = figure('KeyPressFcn',@keypress, ...
    'WindowScrollWheelFcn',@scroll_atlas_slice,...
    'Toolbar','none','Menubar','none','color','w', ...
    'Units','pixels','Position',gui_position, ...
    'Name', gui_data.mouse_name, 'NumberTitle', 'off', ...
    'CloseRequestFcn',@close_gui);


% gui_data.curr_slice = randperm(numel(chooselist), 1);
% curr_image = volumeIdtoImage(gui_data.volume, gui_data.volindids(1, :));
gui_data.curr_slice = 1;
curr_image = squeeze(gui_data.volume(gui_data.curr_slice, :, :, :));

% curr_image = adapthisteq(curr_image);

% Set up axis for histology image
gui_data.histology_ax = subplot(1,2,1,'YDir','reverse'); 
set(gui_data.histology_ax,'Position',[0,0,0.5,0.9]);
hold on; colormap(gray); axis image off;
gui_data.histology_im_h = image(curr_image,...
    'Parent',gui_data.histology_ax,'ButtonDownFcn',@mouseclick_histology);
gui_data.histology_grid = line(xlinesall(:), ylinesall(:), 'Color', 'w', 'LineWidth', 0.5, 'PickableParts', 'none');
set(gui_data.histology_grid,'Visible','off')

% Set up histology-aligned atlas overlay
% (and make it invisible to mouse clicks)
% histology_aligned_atlas_boundaries_init = ...
%     zeros(size(curr_image));
% gui_data.histology_aligned_atlas_boundaries = ...
%     imagesc(histology_aligned_atlas_boundaries_init,'Parent',gui_data.histology_ax, ...
%     'AlphaData',histology_aligned_atlas_boundaries_init,'PickableParts','none');

histology_aligned_atlas_boundaries_init = ...
    zeros(size(curr_image));
gui_data.histology_aligned_atlas_boundaries = ...
    plot(histology_aligned_atlas_boundaries_init(:,1), histology_aligned_atlas_boundaries_init(:,2),...
    'w.','MarkerSize',2, 'Parent',gui_data.histology_ax, 'PickableParts','none');

gui_data.atlas_slice = gui_data.atlasinds(1);
curr_atlas = squeeze(gui_data.tv(gui_data.atlas_slice, :, :));

% Set up axis for atlas slice
gui_data.atlas_ax = subplot(1,2,2,'YDir','reverse'); 
set(gui_data.atlas_ax,'Position',[0.5,0,0.5,0.9]);
hold on; axis image off; colormap(gray); clim([0,255]);
gui_data.atlas_im_h = imagesc(curr_atlas, ...
    'Parent',gui_data.atlas_ax,'ButtonDownFcn',@mouseclick_atlas);
gui_data.atlas_grid = line(xlinesall(:), ylinesall(:), 'Color', 'w', 'LineWidth', 0.5, 'PickableParts', 'none');
set(gui_data.atlas_grid,'Visible','off')

title(gui_data.atlas_ax, sprintf('Atlas slice = %2.2f h-slice widths', ...
        gui_data.atlas_slice/gui_data.slicewidth));

% PickableParts none, or the markers swallow the click. Both axes carry their
% ButtonDownFcn on the IMAGE, and a marker drawn on top of it intercepts the
% press without doing anything -- so clicking straight at a point, which is
% exactly what selecting one means, did nothing at all.
gui_data.histology_control_points_plot = plot(gui_data.histology_ax,nan,nan,...
    'o','MarkerSize', 7, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'g', ...
    'PickableParts', 'none');
gui_data.atlas_control_points_plot = plot(gui_data.atlas_ax,nan,nan,...
    'o','MarkerSize', 7, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'r', ...
    'PickableParts', 'none');

% Highlight for the point currently picked up in edit mode
gui_data.histology_sel_plot = plot(gui_data.histology_ax,nan,nan,'o', ...
    'MarkerSize', 13, 'LineWidth', 2, 'MarkerEdgeColor', [1 0.85 0.1], ...
    'MarkerFaceColor', 'none', 'PickableParts', 'none');
gui_data.atlas_sel_plot = plot(gui_data.atlas_ax,nan,nan,'o', ...
    'MarkerSize', 13, 'LineWidth', 2, 'MarkerEdgeColor', [1 0.85 0.1], ...
    'MarkerFaceColor', 'none', 'PickableParts', 'none');

% If there was previously auto-alignment, intitialize with that
if isfield(gui_data,'histology_ccf_auto_alignment')
    gui_data.histology_ccf_manual_alignment = gui_data.histology_ccf_auto_alignment;
end

% Upload gui data
guidata(gui_fig,gui_data);

% Initialize alignment - FIX!!!
align_ccf_to_histology(gui_fig);

% Print controls
show_controls(gui_fig);

end


function show_controls(gui_fig)
% The controls list, in its own non-modal window.
%
% This lives in a function rather than inline at startup so that 'h' can bring
% it back. It is an ordinary window, easy to close by accident and easy to lose
% behind the main one, and once it is gone nothing on screen says what the keys
% do.

gui_data = guidata(gui_fig);

% Already open: raise it instead of stacking up copies
if isfield(gui_data, 'controls_fig') && ~isempty(gui_data.controls_fig) && ...
        isvalid(gui_data.controls_fig)
    figure(gui_data.controls_fig);
    return
end

CreateStruct.Interpreter = 'tex';
CreateStruct.WindowStyle = 'non-modal';

% Keys in one column and meanings in another, with every line short enough to
% fit: msgbox sizes itself to the longest line, so one long sentence stretches
% the whole window and then wraps anyway.
gui_data.controls_fig = msgbox( ...
    {'\fontsize{11}' ...
    '\bf NAVIGATE\rm', ...
    '   left / right      switch slice', ...
    '   scroll wheel      change atlas plane', ...
    '   1 2 3 / 0         channels  (0 = all)', ...
    '   space             overlay on / off', ...
    '   g                 grid on / off', ...
    '   h                 show this window', ...
    ' ', ...
    '\bf PLACE\rm', ...
    '   click             add a point', ...
    '                     (3 per slice minimum)', ...
    '   ctrl + z          undo the last one', ...
    ' ', ...
    '\bf EDIT   [ e ]\rm', ...
    '   click + drag      grab nearest point, move it', ...
    '   d                 delete it, and its pair', ...
    '   esc               let go', ...
    ' ', ...
    '   a held point also rings its partner', ...
    '   on the other panel. A ring on one', ...
    '   side only means it has no partner.', ...
    ' ', ...
    '\bf CARRY FORWARD   [ p ]\rm', ...
    '   an empty slice inherits a copy of', ...
    '   the one before it, as hollow orange', ...
    '   circles. The wheel still scrolls them', ...
    '   while orange. Touching one commits', ...
    '   the slice and locks the plane.', ...
    '   Slices left orange are NOT saved.', ...
    ' ', ...
    '\bf PROPOSE   [ r ]\rm', ...
    '   scroll to the best-matching atlas', ...
    '   plane, then r: the neighbouring', ...
    '   slice''s points are refined against', ...
    '   this plane by the image matcher and', ...
    '   placed here, orange. A few seconds.', ...
    '   Never overwrites hand-placed points.', ...
    '   t                 same, but a plain copy:', ...
    '                     no matcher, instant', ...
    ' ', ...
    '\bf FINISH\rm', ...
    '   c                 clear every point here', ...
    '   s                 save'}, ...
    'Controls',CreateStruct);

guidata(gui_fig, gui_data);

end


function keypress(gui_fig,eventdata)

% Get guidata
gui_data = guidata(gui_fig);

switch eventdata.Key
    
    % left/right arrows: move slice
    case 'leftarrow'
        gui_data = step_slice(gui_data, -1);
        guidata(gui_fig,gui_data);
        update_slice(gui_fig);

    case 'rightarrow'
        gui_data = step_slice(gui_data, +1);
        guidata(gui_fig,gui_data);
        update_slice(gui_fig);

    % e: edit mode. Clicks grab an existing point instead of adding a new one.
    case 'e'
        gui_data.edit_mode = ~gui_data.edit_mode;
        gui_data.sel_side  = '';
        gui_data.sel_idx   = 0;
        if gui_data.edit_mode
            disp('EDIT mode ON: click near a point to grab it, drag to move, d deletes.');
        else
            disp('EDIT mode off: clicks add points again.');
        end
        guidata(gui_fig,gui_data);
        update_window_title(gui_fig);
        update_slice(gui_fig);

    % p: carry the current slice's points onto the next one
    case 'p'
        gui_data.carry_forward = ~gui_data.carry_forward;
        if gui_data.carry_forward
            disp('CARRY FORWARD on: an empty slice inherits a copy of the one you came from.');
        else
            disp('CARRY FORWARD off.');
        end
        guidata(gui_fig,gui_data);
        update_window_title(gui_fig);

    % d / delete: remove the selected point AND its pair
    case {'d', 'delete', 'backspace'}
        if isempty(gui_data.sel_side)
            disp('Nothing selected. Turn on edit mode with e, then click a point.');
        else
            idx = gui_data.sel_idx;
            h = gui_data.histology_control_points{gui_data.curr_slice};
            a = gui_data.atlas_control_points{gui_data.curr_slice};
            % Both sides, same row: the two lists are matched by index, so
            % dropping one alone would re-pair every later point.
            if size(h,1) >= idx, h(idx,:) = []; end
            if size(a,1) >= idx, a(idx,:) = []; end
            gui_data.histology_control_points{gui_data.curr_slice} = h;
            gui_data.atlas_control_points{gui_data.curr_slice}     = a;
            gui_data.sel_side = '';
            gui_data.sel_idx  = 0;
            gui_data.provisional(gui_data.curr_slice) = false;
            fprintf('Deleted pair %d (%d histology / %d atlas points left).\n', ...
                idx, size(h,1), size(a,1));
            guidata(gui_fig,gui_data);
            update_slice(gui_fig);
        end

    case 'escape'
        gui_data.sel_side = '';
        gui_data.sel_idx  = 0;
        guidata(gui_fig,gui_data);
        update_slice(gui_fig);
        
    % space: toggle overlay visibility
    case 'space'
        curr_visibility = ...
            get(gui_data.histology_aligned_atlas_boundaries,'Visible');
        set(gui_data.histology_aligned_atlas_boundaries,'Visible', ...
            cell2mat(setdiff({'on','off'},curr_visibility)))
         % space: toggle overlay visibility
    case 'g'
        curr_visibility = ...
            get(gui_data.histology_grid,'Visible');
        set(gui_data.histology_grid,'Visible', ...
            cell2mat(setdiff({'on','off'},curr_visibility)))
        set(gui_data.atlas_grid,'Visible', ...
            cell2mat(setdiff({'on','off'},curr_visibility)))
    % bring the controls window back after it has been closed
    case 'h'
        show_controls(gui_fig);

    case {'1', '2', '3'}
        gui_data.colsuse = str2double(eventdata.Key);
        guidata(gui_fig,gui_data);
        update_slice(gui_fig, true);
    case '0'
        gui_data.colsuse = [1 2 3];
        guidata(gui_fig,gui_data);
        update_slice(gui_fig, true);
        
    % r: propose points for THIS slice at the atlas plane on screen, by
    % refining the neighbouring slice's points through the image matcher.
    % Deliberately not automatic -- find the best-matching plane with the
    % wheel first, then ask. Hand-placed points are never overwritten.
    case 'r'
        sl = gui_data.curr_slice;
        if ~isempty(gui_data.histology_control_points{sl}) && ~gui_data.provisional(sl)
            disp('This slice has hand-placed points. Press c to clear them first if you want a proposal.');
        else
            src = source_slice(gui_data, sl);
            if src == 0
                disp('Nothing to propose from: neither neighbouring slice has points.');
            else
                plane = round(gui_data.atlas_slice);
                h = gui_data.histology_control_points{src};
                a = gui_data.atlas_control_points{src};
                fprintf('Proposing %d point(s) for slice %d from slice %d at atlas plane %d...\n', ...
                    size(a,1), sl, src, plane);
                % Whatever goes wrong on the Python side -- no interpreter, a
                % broken install, a corrupt reply -- ends here as a message.
                % A key press must never take the GUI down.
                try
                    [h2, a2, ok, msg] = auto_refine_points(gui_data, sl, plane, h, a);
                catch err
                    ok = false; msg = err.message;
                end
                if ok
                    now_stamp = convertTo(datetime('now'), 'datenum');
                    h2(:, 1) = sl;      h2(:, 4) = now_stamp;
                    a2(:, 1) = plane;   a2(:, 4) = now_stamp;
                    gui_data.histology_control_points{sl} = h2;
                    gui_data.atlas_control_points{sl}     = a2;
                    gui_data.provisional(sl) = true;
                    gui_data.sel_side = '';
                    gui_data.sel_idx  = 0;
                    fprintf('Proposed %d point(s). They are provisional: touch one to keep them, c to discard.\n', ...
                        size(h2,1));
                    guidata(gui_fig, gui_data);
                    update_slice(gui_fig);
                else
                    fprintf('Proposal failed: %s\n', msg);
                end
            end
        end

    % t: take the neighbouring slice's points exactly as they are -- carry-
    % forward on demand, at the atlas plane on screen, no matcher involved.
    % The cheap start when the section barely changed, and the fallback when
    % Python is not there. Same rules as r: provisional, never overwrites.
    case 't'
        sl = gui_data.curr_slice;
        if ~isempty(gui_data.histology_control_points{sl}) && ~gui_data.provisional(sl)
            disp('This slice has hand-placed points. Press c to clear them first.');
        else
            src = source_slice(gui_data, sl);
            if src == 0
                disp('Nothing to take: neither neighbouring slice has points.');
            else
                plane = round(gui_data.atlas_slice);
                h = gui_data.histology_control_points{src};
                a = gui_data.atlas_control_points{src};
                now_stamp = convertTo(datetime('now'), 'datenum');
                h(:, 1) = sl;       h(:, 4) = now_stamp;
                a(:, 1) = plane;    a(:, 4) = now_stamp;
                gui_data.histology_control_points{sl} = h;
                gui_data.atlas_control_points{sl}     = a;
                gui_data.provisional(sl) = true;
                gui_data.sel_side = '';
                gui_data.sel_idx  = 0;
                fprintf('Took %d point(s) from slice %d as they are, at atlas plane %d. Provisional: touch one to keep, c to discard.\n', ...
                    size(h,1), src, plane);
                guidata(gui_fig, gui_data);
                update_slice(gui_fig);
            end
        end

    % ctrl+z: take back the last point placed on this slice
    case 'z'
        if any(strcmp(eventdata.Modifier, 'control'))
            gui_data = undo_last_point(gui_data);
            guidata(gui_fig,gui_data);
            update_slice(gui_fig);
        end

    % c: clear current points
    case 'c'
        gui_data.histology_control_points{gui_data.curr_slice} = zeros(0,3);
        gui_data.atlas_control_points{gui_data.curr_slice} = zeros(0,3);

        guidata(gui_fig,gui_data);
        update_slice(gui_fig);
        
    % s: save
    case 's'
        [histology_control_points, atlas_control_points] = points_for_saving(gui_data);
        save_fn = fullfile(gui_data.save_path,'atlas2histology_tform.mat');
        save(save_fn,'atlas_control_points', 'histology_control_points');
        fprintf('Saved %d annotated slice(s) to %s\n', ...
            nnz(~cellfun(@isempty, histology_control_points)), save_fn);
        
end

end


function mouseclick_histology(gui_fig,eventdata)
% Draw new point for alignment, or in edit mode grab an existing one

% Get guidata
gui_data = guidata(gui_fig);
toplot   = [2 3];

if gui_data.edit_mode
    grab_point(gui_fig, 'histology', flip(eventdata.IntersectionPoint(1:2)));
    return
end

cpt(1)   = gui_data.curr_slice;
cpt(2:3) = flip(eventdata.IntersectionPoint(1:2));

% Placing a point commits the slice: it stops being a carried guess
gui_data.provisional(gui_data.curr_slice) = false;

% Add clicked location to control points
gui_data.histology_control_points{gui_data.curr_slice} = ...
    vertcat(gui_data.histology_control_points{gui_data.curr_slice}, ...
     [cpt convertTo(datetime('now'), 'datenum')]);

set(gui_data.histology_control_points_plot, ...
    'XData',gui_data.histology_control_points{gui_data.curr_slice}(:,toplot(2)), ...
    'YData',gui_data.histology_control_points{gui_data.curr_slice}(:,toplot(1)));

% Upload gui data
guidata(gui_fig, gui_data);

% If equal number of histology/atlas control points > 3, draw boundaries
if size(gui_data.histology_control_points{gui_data.curr_slice},1) == ...
        size(gui_data.atlas_control_points{gui_data.curr_slice},1) || ...
        (size(gui_data.histology_control_points{gui_data.curr_slice},1) > 3 && ...
        size(gui_data.atlas_control_points{gui_data.curr_slice},1) > 3)
    align_ccf_to_histology(gui_fig)
end

end


function mouseclick_atlas(gui_fig,eventdata)
% Draw new point for alignment, or in edit mode grab an existing one

% Get guidata
gui_data = guidata(gui_fig);

if gui_data.edit_mode
    grab_point(gui_fig, 'atlas', flip(eventdata.IntersectionPoint(1:2)));
    return
end

cpt      = zeros(1,3);
cpt(1)   = gui_data.atlas_slice;
cpt(2:3) =  flip(eventdata.IntersectionPoint(1:2));
toplot   = [2 3];

% Placing a point commits the slice: it stops being a carried guess
gui_data.provisional(gui_data.curr_slice) = false;

% Add clicked location to control points
gui_data.atlas_control_points{gui_data.curr_slice} = ...
    vertcat(gui_data.atlas_control_points{gui_data.curr_slice}, ...
    [cpt convertTo(datetime('now'), 'datenum')]);

set(gui_data.atlas_control_points_plot, ...
    'XData',gui_data.atlas_control_points{gui_data.curr_slice}(:,toplot(2)), ...
    'YData',gui_data.atlas_control_points{gui_data.curr_slice}(:,toplot(1)));

% Upload gui data
guidata(gui_fig, gui_data);

% If equal number of histology/atlas control points > 3, draw boundaries
if size(gui_data.histology_control_points{gui_data.curr_slice},1) == ...
        size(gui_data.atlas_control_points{gui_data.curr_slice},1) || ...
        (size(gui_data.histology_control_points{gui_data.curr_slice},1) > 3 && ...
        size(gui_data.atlas_control_points{gui_data.curr_slice},1) > 3)
    align_ccf_to_histology(gui_fig)
end

end


function align_ccf_to_histology(gui_fig)

% Get guidata 
gui_data = guidata(gui_fig);


Nmin = 3;
cptsatlas     = gui_data.atlas_control_points{gui_data.curr_slice};
idatlas       = gui_data.atlas_slice;
cptsatlas     = cptsatlas(:, [3 2]); 
cptshistology = gui_data.histology_control_points{gui_data.curr_slice};
cptshistology = cptshistology(:, [3 2]);
currim        = squeeze(gui_data.av(idatlas, :, :));

tstrcurr      = sprintf('Slice %d/%d', gui_data.curr_slice, gui_data.Nslices);

% Second title line: which atlas plane this slice sits on, and whether that is
% its own or borrowed. A slice carrying points is ANCHORED there and stays put;
% one without is only where the interpolation between the surrounding anchors
% currently puts it, and will move as neighbours get annotated. Same units as
% the atlas panel, so the two read against each other directly.
anchorpts = gui_data.atlas_control_points{gui_data.curr_slice};
if ~isempty(anchorpts)
    anchorstr = sprintf('anchored at atlas %2.2f h-slice widths', ...
        median(anchorpts(:,1))/gui_data.slicewidth);
elseif isfield(gui_data, 'atlasindsuse')
    anchorstr = sprintf('not anchored, predicted atlas %2.2f h-slice widths', ...
        gui_data.atlasindsuse(gui_data.curr_slice)/gui_data.slicewidth);
else
    anchorstr = 'not anchored';
end

if size(cptshistology,1) == size(cptsatlas,1) && ...
        (size(cptshistology,1) >= Nmin && size(cptsatlas,1) >= Nmin)
    
    tform = fitgeotform2d(cptsatlas, cptshistology, "affine");
    mse = mean(sqrt(sum((cptshistology-tform.transformPointsForward(cptsatlas)).^2, 2)));

    % tform = estgeotform3d(cptsatlas, cptshistology, 'similarity');
    currim = imwarp(currim, gui_data.Rmoving, tform, 'OutputView',gui_data.Rfixed);

    tstrcurr = sprintf('%s, Npoints = %d, mse = %2.2f ', tstrcurr, size(cptshistology,1), mse);

elseif size(gui_data.histology_control_points{gui_data.curr_slice},1) >= 1 ||  ...
        size(gui_data.atlas_control_points{gui_data.curr_slice},1) >= 1
    % If less than 3 or nonmatching points, use auto but don't draw
    tstrcurr = sprintf('%s, New alignment ', tstrcurr);

    % This branch returns before the title call below, so set it here too --
    % otherwise a slice part-way through a pair keeps the previous slice's
    % title, anchor line included.
    title(gui_data.histology_ax, {tstrcurr, anchorstr})

    % Upload gui data
    guidata(gui_fig, gui_data);
    return

elseif isfield(gui_data,'histology_ccf_auto_alignment')
    % If no points, use automated outline if available
    tform = gui_data.histology_ccf_auto_alignment;
    tform = affinetform2d(tform);
    tstrcurr = sprintf('%s, Previous alignment ', tstrcurr);
    currim = imwarp(currim, gui_data.Rmoving, tform, 'OutputView',gui_data.Rfixed);
else
    tform = affinetform2d;
end

title(gui_data.histology_ax, {tstrcurr, anchorstr})

av_warp_boundaries = round(conv2(currim,ones(3)./9,'same')) ~= currim;

[row,col] = ind2sub(size(currim), find(av_warp_boundaries));


 % set(gui_data.histology_aligned_atlas_boundaries, ...
 %    'CData',av_warp_boundaries, ...
 %    'AlphaData',av_warp_boundaries*0.6);
set(gui_data.histology_aligned_atlas_boundaries, ...
'XData', col, 'YData', row);

% Update transform matrix
gui_data.histology_ccf_manual_alignment = tform.A;

% Upload gui data
guidata(gui_fig, gui_data);

end


function update_slice(gui_fig, varargin)
% Draw histology and CCF slice

if nargin < 2 
    sliceonly = false;
else
    sliceonly = varargin{1};
end

% Get guidata
gui_data = guidata(gui_fig);

atlas_cpoints    = gui_data.atlas_control_points;
hascp            = ~cellfun(@isempty,   atlas_cpoints);
useratlasinds    = cellfun(@(x) x(1,1), atlas_cpoints(hascp));
replaceinds      = gui_data.atlasinds;

% How many atlas planes one sample slice is worth. This is fixed by the
% section spacing, so it is known before any point is placed and does not have
% to be inferred from the annotations.
known_step = (gui_data.atlasinds(end) - gui_data.atlasinds(1)) / ...
             max(1, gui_data.Nslices - 1);

if nnz(hascp) > 0
    meanrem     = mean(gui_data.atlasinds(hascp));
    newinds     = gui_data.atlasinds - meanrem + mean(useratlasinds);
    replaceinds = round(newinds);
end
if nnz(hascp) > 1
    % refine remaining
    pfit = polyfit(gui_data.atlasinds(hascp), useratlasinds, 1);
    % A fit that runs backwards means later slices would sit on earlier atlas
    % planes, which cannot be right. Keep the offset, drop the bad slope.
    if pfit(1) > 0
        replaceinds = round(polyval(pfit, gui_data.atlasinds));
    else
        replaceinds = round(gui_data.atlasinds - mean(gui_data.atlasinds(hascp)) ...
                            + mean(useratlasinds));
    end
end
if nnz(hascp) > 3
    % Interpolate BETWEEN the annotated slices, but step out beyond them at the
    % known rate instead of continuing the slope of the last pair.
    %
    % The original used interp1(..., 'extrap'), which takes its slope from the
    % final two annotated slices alone. If those are adjacent and the second
    % sits even one plane behind the first, every slice past them marches
    % backwards, runs off the front of the atlas and gets clamped to plane 1 --
    % which is the atlas appearing to travel the wrong way as you advance.
    anchors = find(hascp);
    allsl   = (1:gui_data.Nslices)';
    out     = nan(gui_data.Nslices, 1);

    inside        = allsl >= anchors(1) & allsl <= anchors(end);
    out(inside)   = interp1(anchors, useratlasinds, allsl(inside), 'linear');

    before        = allsl < anchors(1);
    out(before)   = useratlasinds(1)   + (allsl(before) - anchors(1))   * known_step;

    after         = allsl > anchors(end);
    out(after)    = useratlasinds(end) + (allsl(after)  - anchors(end)) * known_step;

    replaceinds = round(out);
end

newinds = replaceinds;

% Say so when the mapping has gone somewhere impossible, rather than silently
% pinning slices to the first or last plane and leaving it looking like the
% atlas is running backwards. Record what is wrong rather than printing it once
% into a console nobody is watching. update_slice paints the title purple off this,
% so the problem is visible in the window for as long as it lasts.
gui_data.order_problem = '';

if nnz(hascp) > 1
    bad = find(diff(useratlasinds) < 0);
    if ~isempty(bad)
        anchors_all = find(hascp);
        gui_data.order_problem = sprintf(...
            'ORDER CONFLICT: slice %d (plane %d) is behind slice %d (plane %d)', ...
            anchors_all(bad(1)+1), useratlasinds(bad(1)+1), ...
            anchors_all(bad(1)),   useratlasinds(bad(1)));
    end
end

Natlas = size(gui_data.tv, 1);
n_clamped = nnz(newinds < 1 | newinds > Natlas);
if n_clamped > 0 && isempty(gui_data.order_problem)
    gui_data.order_problem = sprintf('%d slice(s) map off the end of the atlas', n_clamped);
end

newinds(newinds<1)      = 1;
newinds(newinds>Natlas) = Natlas;

gui_data.atlasindsuse = newinds;


% Unchanged from the original: the displayed plane always comes from the
% points when a slice has any, otherwise from the prediction. Arriving at a
% carried slice therefore shows the plane step_slice stamped into the copies,
% not whatever the previous slice happened to be on. Whether the wheel can
% MOVE that plane is decided in scroll_atlas_slice, not here.
cpointsatlas = gui_data.atlas_control_points{gui_data.curr_slice};
if ~isempty(cpointsatlas)
    gui_data.atlas_slice = median(cpointsatlas(:,1));
else
    gui_data.atlas_slice = gui_data.atlasindsuse(gui_data.curr_slice);
end


toplot     = [2 3];

% Set next histology slice

curr_image = squeeze(gui_data.volume(gui_data.curr_slice, :, :, gui_data.colsuse));
set(gui_data.histology_im_h,'CData', curr_image)

% Plot control points for slice
set(gui_data.histology_control_points_plot, ...
    'XData',gui_data.histology_control_points{gui_data.curr_slice}(:,toplot(2)), ...
    'YData',gui_data.histology_control_points{gui_data.curr_slice}(:,toplot(1)));

% Carried-but-untouched points are drawn as hollow orange circles, so it is
% obvious at a glance that they are a guess inherited from the previous slice
% and that the atlas plane is still free to scroll. They become solid green and
% red the moment the slice is committed.
if gui_data.provisional(gui_data.curr_slice)
    set(gui_data.histology_control_points_plot, ...
        'MarkerFaceColor', 'none', 'MarkerEdgeColor', [0.95 0.55 0.10], 'LineWidth', 1.5);
    set(gui_data.atlas_control_points_plot, ...
        'MarkerFaceColor', 'none', 'MarkerEdgeColor', [0.95 0.55 0.10], 'LineWidth', 1.5);
else
    set(gui_data.histology_control_points_plot, ...
        'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'g', 'LineWidth', 0.5);
    set(gui_data.atlas_control_points_plot, ...
        'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'r', 'LineWidth', 0.5);
end

histology_aligned_atlas_boundaries_init = nan(1,2);
set(gui_data.histology_aligned_atlas_boundaries, ...
    'XData',histology_aligned_atlas_boundaries_init(:,1), 'YData',histology_aligned_atlas_boundaries_init(:,2));

% Upload gui data
guidata(gui_fig, gui_data);

% Update atlas boundaries
align_ccf_to_histology(gui_fig)

% Keep the edit-mode highlight in step with whatever is on screen
draw_selection(gui_fig)
draw_point_numbers(gui_fig)

if ~sliceonly

    % update atlas slice
    update_atlas_slice(gui_fig)
    % update title
    % set_histology_title(gui_fig)

    % % clear points that are not used
    % cptsatlas     = gui_data.atlas_control_points;
    % cptshistology = gui_data.histology_control_points;
    % Nptsatlas     = cellfun(@(x) size(x, 1), cptsatlas);
    % Nptshisto     = cellfun(@(x) size(x, 1), cptshistology);
    % badpts        = find(~(Nptsatlas == Nptshisto));
    % for ii = 1:numel(badpts)
    %     % curratlas = cptsatlas{badpts(ii)};
    %     % currhisto = cptshistology{badpts(ii)};
    %     % temp fix, clearing bad stuff, use timestamps later
    %     gui_data.atlas_control_points{badpts(ii)} = {zeros(0,4)};
    %     gui_data.histology_control_points{badpts(ii)} = {zeros(0,4)};
    % end
end

end


function scroll_atlas_slice(gui_fig,eventdata)
% Move point to draw atlas slice perpendicular to the camera

% Get guidata
gui_data = guidata(gui_fig);

% Move slice point
gui_data.atlas_slice = ...
    gui_data.atlas_slice + eventdata.VerticalScrollCount;

gui_data.atlas_slice = max(gui_data.atlas_slice, 1);
gui_data.atlas_slice = min(gui_data.atlas_slice, size(gui_data.tv, 1));

% Carried-but-untouched points follow the wheel; committed ones do not.
%
% update_atlas_slice pins the displayed plane to median(points(:,1)) once a
% slice has atlas points, which is deliberate -- it stops an accidental scroll
% silently relocating work already placed by hand. That protection stays. Only
% provisional points, copied from the previous slice and not yet touched, are
% re-stamped so the wheel can still be used to find the right plane.
if gui_data.provisional(gui_data.curr_slice)
    cpts = gui_data.atlas_control_points{gui_data.curr_slice};
    if ~isempty(cpts)
        gui_data.atlas_control_points{gui_data.curr_slice}(:,1) = gui_data.atlas_slice;
    end
end

% Upload gui data
guidata(gui_fig, gui_data);

% Update slice
update_atlas_slice(gui_fig)

end

function close_gui(gui_fig,~)

% Get guidata
gui_data = guidata(gui_fig);

opts.Default = 'Yes';
opts.Interpreter = 'tex';
user_confirm = questdlg('\fontsize{14} Save?','Confirm exit',opts);
switch user_confirm
    case 'Yes'
        % Save and close
        atlas2histology_tform = ...
            gui_data.histology_ccf_manual_alignment;
        [histology_control_points, atlas_control_points] = points_for_saving(gui_data);
        save_fn = fullfile(gui_data.save_path,'atlas2histology_tform.mat');

        save(save_fn,'atlas2histology_tform', 'atlas_control_points', 'histology_control_points');
        fprintf('Saved %d annotated slice(s) to %s\n', ...
            nnz(~cellfun(@isempty, histology_control_points)), save_fn);
        delete(gui_fig);

    case 'No'
        % Close without saving
        delete(gui_fig);

    case 'Cancel'
        % Do nothing

end   

end


function update_atlas_slice(gui_fig)
% Draw atlas slice through plane perpendicular to camera through set point

% Get guidata
gui_data = guidata(gui_fig);

toplot     = [2 3];
sluse      = gui_data.atlas_slice;
curr_atlas = squeeze(gui_data.tv(sluse, :, :));
curr_atlas = adapthisteq(curr_atlas);

set(gui_data.atlas_im_h,'CData', curr_atlas);
set(gui_data.atlas_control_points_plot, ...
    'XData',gui_data.atlas_control_points{gui_data.curr_slice}(:,toplot(2)), ...
    'YData',gui_data.atlas_control_points{gui_data.curr_slice}(:,toplot(1)));

% Reset histology-aligned atlas boundaries if not
% histology_aligned_atlas_boundaries_init = zeros(size(curr_image));
% set(gui_data.histology_aligned_atlas_boundaries, ...
%     'CData',histology_aligned_atlas_boundaries_init, ...
%     'AlphaData',histology_aligned_atlas_boundaries_init);


% Second line: where the neighbouring slices are anchored, in the same units,
% so scrolling for the best plane has a floor and a ceiling on screen instead
% of in memory -- going further back than the previous slice is a mistake
% that is easy to make and hard to notice.
sl = gui_data.curr_slice;
nb = {};
if sl > 1 && ~isempty(gui_data.atlas_control_points{sl-1})
    nb{end+1} = sprintf('slice %d is at %2.2f', sl-1, ...
        median(gui_data.atlas_control_points{sl-1}(:,1))/gui_data.slicewidth);
end
if sl < gui_data.Nslices && ~isempty(gui_data.atlas_control_points{sl+1})
    nb{end+1} = sprintf('slice %d is at %2.2f', sl+1, ...
        median(gui_data.atlas_control_points{sl+1}(:,1))/gui_data.slicewidth);
end
if isempty(nb)
    nbline = 'no annotated neighbour';
else
    nbline = strjoin(nb, '    |    ');
end
lines = { sprintf('Atlas slice = %2.2f h-slice widths', sluse/gui_data.slicewidth), nbline };

if isempty(gui_data.order_problem)
    title(gui_data.atlas_ax, lines, 'Color', 'k', 'FontWeight', 'normal');
else
    % Purple, in the window, naming the slices. The registration code is left
    % exactly as it was, so a contradictory order is not caught downstream --
    % it has to be caught here, while it can still be corrected.
    title(gui_data.atlas_ax, [lines, {gui_data.order_problem}], ...
        'Color', [0.55 0.10 0.75], 'FontWeight', 'bold');
end

% Upload gui_data
guidata(gui_fig, gui_data);

end


function gui_data = step_slice(gui_data, dir)
% Move one slice, optionally carrying the current points along.
%
% Copying keeps the in-plane positions but re-stamps the slice index: the
% histology point gets the new sample slice, and the atlas point gets the atlas
% plane predicted for it, so the copy lands on the right section rather than
% dragging the previous slice's atlas plane along with it.

from_slice = gui_data.curr_slice;
to_slice   = min(max(from_slice + dir, 1), gui_data.Nslices);

if gui_data.carry_forward && to_slice ~= from_slice
    src_h = gui_data.histology_control_points{from_slice};
    src_a = gui_data.atlas_control_points{from_slice};
    dst_h = gui_data.histology_control_points{to_slice};
    dst_a = gui_data.atlas_control_points{to_slice};

    % Only seed an empty slice. Never overwrite work already done there.
    if isempty(dst_h) && isempty(dst_a) && ~isempty(src_h)
        now_stamp   = convertTo(datetime('now'), 'datenum');
        target_plane = gui_data.atlasindsuse(to_slice);

        new_h = src_h;
        new_h(:, 1) = to_slice;
        new_h(:, 4) = now_stamp;

        new_a = src_a;
        if ~isempty(new_a)
            new_a(:, 1) = target_plane;
            new_a(:, 4) = now_stamp;
        end

        gui_data.histology_control_points{to_slice} = new_h;
        gui_data.atlas_control_points{to_slice}     = new_a;
        gui_data.provisional(to_slice)              = true;
        fprintf('Carried %d point(s) from slice %d to %d (atlas plane %d). Press r for a refined proposal.\n', ...
            size(new_h,1), from_slice, to_slice, target_plane);
    end
end

gui_data.curr_slice = to_slice;
gui_data.sel_side   = '';
gui_data.sel_idx    = 0;
end


function grab_point(gui_fig, side, yx)
% Pick the nearest existing point on this panel and start dragging it.

gui_data = guidata(gui_fig);

if strcmp(side, 'histology')
    pts = gui_data.histology_control_points{gui_data.curr_slice};
    width = size(gui_data.volume, 3);
else
    pts = gui_data.atlas_control_points{gui_data.curr_slice};
    width = size(gui_data.av, 3);
end

if isempty(pts)
    disp('No points on this slice to grab.');
    return
end

% Columns 2 and 3 hold y and x
d = hypot(pts(:,2) - yx(1), pts(:,3) - yx(2));
[dmin, idx] = min(d);

if dmin > gui_data.select_radius_frac * width
    fprintf('No point within reach (nearest is %.0f px away).\n', dmin);
    return
end

gui_data.sel_side = side;
gui_data.sel_idx  = idx;
gui_data.provisional(gui_data.curr_slice) = false;   % touching it commits the slice
guidata(gui_fig, gui_data);
draw_selection(gui_fig);

% Follow the mouse until the button comes back up.
%
% ButtonDownFcn hands the callback the object that was CLICKED, so the handle
% arriving here is the image, not the figure. guidata does not care -- it walks
% up to the figure either way, which is why selecting and deleting worked -- but
% WindowButtonMotionFcn only exists on a figure, so setting it on the image
% quietly failed and no drag was ever installed.
fig = ancestor(gui_fig, 'figure');
set(fig, 'WindowButtonMotionFcn', @drag_motion, ...
         'WindowButtonUpFcn',     @drag_stop);
end


function drag_motion(gui_fig, ~)
% Move the held point to wherever the pointer is now

gui_fig  = ancestor(gui_fig, 'figure');
gui_data = guidata(gui_fig);
if isempty(gui_data.sel_side)
    return
end

if strcmp(gui_data.sel_side, 'histology')
    ax = gui_data.histology_ax;
else
    ax = gui_data.atlas_ax;
end
cp = get(ax, 'CurrentPoint');
x  = cp(1,1);
y  = cp(1,2);

% Ignore excursions outside the panel rather than flinging the point away
xl = xlim(ax); yl = ylim(ax);
if x < xl(1) || x > xl(2) || y < yl(1) || y > yl(2)
    return
end

if strcmp(gui_data.sel_side, 'histology')
    gui_data.histology_control_points{gui_data.curr_slice}(gui_data.sel_idx, 2:3) = [y x];
else
    gui_data.atlas_control_points{gui_data.curr_slice}(gui_data.sel_idx, 2:3) = [y x];
end

guidata(gui_fig, gui_data);
redraw_points(gui_fig);
draw_selection(gui_fig);
end


function drag_stop(gui_fig, ~)
% Let go, and refresh the alignment with the point in its new place

fig = ancestor(gui_fig, 'figure');
set(fig, 'WindowButtonMotionFcn', '', 'WindowButtonUpFcn', '');
update_slice(fig);
end


function redraw_points(gui_fig)
% Refresh both point plots from the stored coordinates

gui_data = guidata(gui_fig);
h = gui_data.histology_control_points{gui_data.curr_slice};
a = gui_data.atlas_control_points{gui_data.curr_slice};
set(gui_data.histology_control_points_plot, ...
    'XData', h(:,3), 'YData', h(:,2));
set(gui_data.atlas_control_points_plot, ...
    'XData', a(:,3), 'YData', a(:,2));
draw_point_numbers(gui_fig)     % the numbers follow a dragged point
end


function draw_selection(gui_fig)
% Ring the point being held AND its partner on the other panel.
%
% The two lists are paired by row -- row i of the histology points belongs
% with row i of the atlas points, which is exactly the pairing the affine
% fit later relies on -- so the partner is just the same index on the other
% side. Same colour on both, so they read as one pair; the held one gets the
% thicker ring so it stays obvious which of the two the mouse has hold of.
%
% A ring that appears on one side only means that index has no partner, which
% makes an unpaired point visible instead of something you find out about at
% registration time.

gui_data = guidata(gui_fig);
set(gui_data.histology_sel_plot, 'XData', nan, 'YData', nan);
set(gui_data.atlas_sel_plot,     'XData', nan, 'YData', nan);

if isempty(gui_data.sel_side)
    return
end

idx      = gui_data.sel_idx;
hpts     = gui_data.histology_control_points{gui_data.curr_slice};
apts     = gui_data.atlas_control_points{gui_data.curr_slice};
heldhist = strcmp(gui_data.sel_side, 'histology');

set_sel_ring(gui_data.histology_sel_plot, hpts, idx,  heldhist);
set_sel_ring(gui_data.atlas_sel_plot,     apts, idx, ~heldhist);
end


function [hpts, apts] = points_for_saving(gui_data)
% The points as they should reach disk: everything placed by hand, and nothing
% that is still a carried guess.
%
% A provisional slice holds a copy of its neighbour that nobody has looked at.
% Those are dropped rather than written, because the provisional flag itself is
% NOT saved -- once on disk a carried copy is indistinguishable from a real
% placement, so writing one is a door that does not open again, and the
% registration ends up anchored on guesses that look like decisions. An empty
% slice is the honest record, and registerSlicesToAtlas already interpolates a
% plane for slices that have no points of their own.
%
% Touching a slice at all -- placing, grabbing or deleting a point -- clears its
% provisional flag. So accepting a carried set exactly as it stands is a matter
% of grabbing any one of its points and letting go.

hpts = gui_data.histology_control_points;
apts = gui_data.atlas_control_points;

prov = find(gui_data.provisional(:)');
for k = prov
    hpts{k} = zeros(0,4);
    apts{k} = zeros(0,4);
end

if ~isempty(prov)
    fprintf(['NOTE: %d slice(s) hold carried, untouched points and were NOT saved:\n' ...
             '      %s\n' ...
             '      Grab any point on one to accept it as it stands.\n'], ...
        numel(prov), mat2str(prov));
end

end


function gui_data = undo_last_point(gui_data)
% Take back the last control point placed on this slice.
%
% Nothing has to be recorded to know which one that was. Column 4 of every
% point is the datenum it was clicked -- the original GUI has always written
% it -- so the two lists already carry the placement order between them. Points
% are appended, so the newest on each side is its last row, and comparing those
% two stamps says which side was touched last.
%
% One point, not the pair: that is how they go down, so ctrl+z once undoes a
% misclick without throwing away its partner, and twice walks back a full pair.

sl = gui_data.curr_slice;
h  = gui_data.histology_control_points{sl};
a  = gui_data.atlas_control_points{sl};

if isempty(h) && isempty(a)
    disp('Nothing to undo on this slice.');
    return
end

if stamp_of_last(h) >= stamp_of_last(a)
    h(end,:) = [];
    gui_data.histology_control_points{sl} = h;
    sidestr = 'histology';
else
    a(end,:) = [];
    gui_data.atlas_control_points{sl} = a;
    sidestr = 'atlas';
end

fprintf('Undid the last %s point (%d histology / %d atlas left).\n', ...
    sidestr, size(h,1), size(a,1));

% Rows shift under a removal, so whatever was selected no longer means what
% it did.
gui_data.sel_side = '';
gui_data.sel_idx  = 0;

% Taking a point back still counts as having touched the slice
gui_data.provisional(sl) = false;
end


function t = stamp_of_last(pts)
% When the newest point in a list was placed. An empty list loses every
% comparison; a list from before column 4 existed still beats an empty one,
% but never a real timestamp.

if isempty(pts)
    t = -inf;
elseif size(pts,2) < 4
    t = 0;
else
    t = pts(end,4);
end
end


function set_sel_ring(hplot, pts, idx, isheld)
% Place one highlight ring, thick for the point under the mouse and thin for
% its partner. Nothing is drawn if that index does not exist on this side.

if idx < 1 || idx > size(pts,1)
    return
end

if isheld
    lw = 2;   ms = 13;
else
    lw = 1;   ms = 15;
end

set(hplot, 'XData', pts(idx,3), 'YData', pts(idx,2), ...
    'LineWidth', lw, 'MarkerSize', ms);
end


function update_window_title(gui_fig)
% Keep the active modes in the window title bar.
%
% Both modes change what a click does, and carry-forward quietly writes points
% onto slices you walk through, so leaving them visible only in the console is
% asking for a session spent in a mode nobody meant to be in.

gui_data = guidata(gui_fig);

modes = {};
if gui_data.edit_mode
    modes{end+1} = 'EDIT';
end
if gui_data.carry_forward
    modes{end+1} = 'CARRY FORWARD';
end

% The controls window is easy to close and easy to lose behind this one, so
% the way back to it rides along in the title bar.
reminder = '        [ h ] controls';

if isempty(modes)
    set(gui_fig, 'Name', [gui_data.mouse_name reminder]);
else
    set(gui_fig, 'Name', [gui_data.mouse_name '     >>  ' ...
        strjoin(modes, '  +  ') '  <<' reminder]);
end
end


function [h, a, ok, msg] = auto_refine_points(gui_data, slice, plane, h, a)
% Hand the two panels and the atlas landmarks to the matcher, take back
% refined coordinates for both sides.
%
% Point columns here are [slice/plane y x t]; the matcher speaks [x y]. Rows
% are pairs, so the two lists must have the same count, and any pair whose
% transfer lands outside the histology panel is dropped from BOTH lists so the
% pairing stays intact. The plane index and timestamps are kept as they were.
%
% The panels are built exactly as the GUI draws them -- adapthisteq on the
% atlas, the autofluorescence channel of the sample -- because that is what
% the method was measured on.

ok  = false;
msg = '';

if isempty(a) || size(h,1) ~= size(a,1)
    msg = 'histology and atlas point counts differ';
    return
end
if ~exist('landmark_refine', 'file')
    msg = 'landmark_refine.m is not on the path';
    return
end

plane    = min(max(round(plane), 1), size(gui_data.tv, 1));
atlas_im = adapthisteq(squeeze(gui_data.tv(plane, :, :)));
hist_im  = squeeze(gui_data.volume(slice, :, :, gui_data.refine_chan));
labels   = squeeze(gui_data.av(plane, :, :));

out = landmark_refine('refine', atlas_im, hist_im, a(:, [3 2]), labels);
if ~out.ok
    msg = out.message;
    return
end

[H, W] = size(hist_im);
inside = out.hist_pts(:,1) >= 1 & out.hist_pts(:,1) <= W & ...
         out.hist_pts(:,2) >= 1 & out.hist_pts(:,2) <= H;
if nnz(inside) < 3
    msg = 'the transferred points fell outside the image';
    return
end

h = h(inside, :);
a = a(inside, :);
a(:, [3 2]) = out.atlas_pts(inside, :);
h(:, [3 2]) = out.hist_pts(inside, :);
ok = true;
end


function draw_point_numbers(gui_fig)
% A small index next to every marker on both panels, so the pairing is always
% visible: the same number is the same correspondence. Orange while the slice
% is provisional, white once it is committed, like the markers themselves.

gui_data = guidata(gui_fig);
delete(gui_data.hist_labels(isvalid(gui_data.hist_labels)));
delete(gui_data.atlas_labels(isvalid(gui_data.atlas_labels)));

h = gui_data.histology_control_points{gui_data.curr_slice};
a = gui_data.atlas_control_points{gui_data.curr_slice};
if gui_data.provisional(gui_data.curr_slice)
    col = [0.95 0.55 0.10];
else
    col = [1 1 1];
end
gui_data.hist_labels  = label_points(gui_data.histology_ax, h, col);
gui_data.atlas_labels = label_points(gui_data.atlas_ax,     a, col);
guidata(gui_fig, gui_data);
end


function t = label_points(ax, pts, col)
% Text objects cannot be given PickableParts none through the plot call, so
% they are created one by one and told not to catch the mouse.
t = gobjects(size(pts,1), 1);
for k = 1:size(pts,1)
    t(k) = text(ax, pts(k,3) + 6, pts(k,2) - 5, sprintf('%d', k), ...
        'Color', col, 'FontSize', 8, 'FontWeight', 'bold', ...
        'PickableParts', 'none', 'Clipping', 'on');
end
end


function src = source_slice(gui_data, sl)
% Where r and t take their points from: the slice before if it has any,
% else the slice after, else 0.
src = 0;
if sl > 1 && ~isempty(gui_data.atlas_control_points{sl-1})
    src = sl - 1;
elseif sl < gui_data.Nslices && ~isempty(gui_data.atlas_control_points{sl+1})
    src = sl + 1;
end
end
