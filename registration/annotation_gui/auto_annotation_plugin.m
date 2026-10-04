function value = auto_annotation_plugin(event, gui_fig, gui_data, value, info)
%AUTO_ANNOTATION_PLUGIN  The automatic annotation's part of the control-point GUI.
%   LightSuite's control-point GUI (matchControlPointsInSlices) calls this when
%   it is passed in as opts.plugin, which run_register_to_atlas's 'annotate'
%   mode does when the automatic annotation is installed (auto_annotate
%   'check'). It adds to the GUI:
%
%     a       fix the atlas plane on screen as this slice's anchor (again on
%             the same plane: remove it)
%     j       jump to the next suggested anchor slice
%     k / K   accept this proposed slice and go to the next / accept every
%             orange slice, after a confirmation
%     u / U   re-propose this slice at the plane on screen / every orange
%             slice at the plane the anchors and the accepted slices give it,
%             through auto_annotate
%     ?       after a point's number: the proposal was least sure of it
%
%   and with them: the anchors in the GUI's plane prediction; the proposal
%   (auto_proposal_controlpoints.mat, from run_register_to_atlas's
%   'autoannotate') loaded orange onto every slice still empty; the stage of
%   the annotation in the window title and the slice title; and on every save
%   the files the automatic annotation reads next to the annotation:
%   plane_anchors.mat, auto_atlas_planes.mat and annotation_provenance.mat.
%   None of them ends in tform.mat, which registerSlicesToAtlas globs for. A
%   brain with no anchors and no proposal opens as it would without the
%   plugin, apart from the suggested anchors in the titles. Its settings (how
%   many anchor slices to suggest, the * mark's outlier rule) are
%   annotation_settings'.
%
%   The hook, as the GUI's help describes it:
%
%     value = AUTO_ANNOTATION_PLUGIN(event, gui_fig, gui_data, value, info)
%
%   gui_fig is the GUI's window and gui_data its state at that moment; value
%   is what the GUI goes on with, returned changed or not; info.gui holds the
%   GUI's update_slice, update_window_title and step_slice. The events:
%
%     'open'    the window is built, nothing drawn yet: loads the anchors and
%               the proposal, sets the GUI's outlier_rule, controls_extra and
%               unsaved_hint, stores it all (guidata) and draws slice 1 when
%               there is something to review
%     'key'     a key the GUI does not use (info.key): a, j, k, u above
%     'planes'  value, the planes known from points: the anchors added
%     'title'   value, the slice title's two lines: the review stage, the
%               anchor or suggested anchor of a slice without points
%     'labels'  value, the marks on the point numbers: the ? flags
%     'window'  value, the modes in the window title: the next step
%     'edit'    value is gui_data: the ? flags follow d (info.row removed
%               from info.slice), and go with c and t (info.row empty)
%     'save'    info.histology and info.atlas as saved: writes the three files
%
%   Any other event returns value unchanged.
%
%   See also matchControlPointsInSlices, auto_annotate, annotation_settings.

switch event
    case 'open'
        open_review(gui_fig, gui_data, info.gui);
    case 'key'
        on_key(gui_fig, gui_data, info.key, info.gui);
    case 'planes'
        value = add_anchor_planes(gui_data, value);
    case 'title'
        value = slice_title(gui_data, value);
    case 'labels'
        value = flag_points(gui_data, value);
    case 'window'
        value = window_modes(gui_data, value);
    case 'edit'
        value = follow_edit(value, info);
    case 'save'
        save_auto_files(gui_fig, info.histology, info.atlas, info.gui.update_slice);
end
end

% ===== Local functions: events =====

function open_review(gui_fig, gui_data, gui)
% Event 'open': the anchors and the proposal, loaded before anything is
% drawn, then slice 1 drawn as the review needs it.

annot = annotation_settings();

% what the GUI draws itself, with this project's settings and keys: the * mark's
% rule, the keys in the controls window, and how to accept an orange slice in the
% note printed when a save leaves one out
gui_data.outlier_rule   = annot.outlier_rule;
gui_data.controls_extra = controls_lines(annot);
gui_data.unsaved_hint   = 'k accepts one as it stands (K all of them), as does grabbing any of its points.';

% per slice, one flag per row: the proposed points the automatic annotation was
% least sure of, drawn with a ? after their number
gui_data.uncertain    = cell(gui_data.Nslices, 1);

% the anchors: a slice's atlas plane fixed by hand without points (a). They join
% the slices with points in the plane interpolation, so four of them spread along
% the stack put every slice within about a section of its plane; saved to
% plane_anchors.mat with the plane the interpolation gives every slice, which the
% automatic step uses as they are. Nothing changes for a brain that never used it.
gui_data.plane_anchors     = nan(gui_data.Nslices, 1);
gui_data.suggested_anchors = unique(round(linspace(1, gui_data.Nslices, ...
    annot.suggested_anchors)));
anchor_fn = fullfile(gui_data.save_path, 'plane_anchors.mat');
if exist(anchor_fn, 'file')
    A = load(anchor_fn, 'anchor_slices', 'anchor_planes');
    gui_data.plane_anchors(A.anchor_slices(:)) = A.anchor_planes(:);
    fprintf('Loaded %d plane anchor(s) from %s\n', numel(A.anchor_slices), anchor_fn);
end

% the proposal (auto_proposal_controlpoints.mat, from 'autoannotate') fills every
% slice still empty in the saved annotation, provisional: orange, and saved only
% when accepted (k) or touched, its least confident points marked ?
gui_data.proposal     = [];
gui_data.proposal_low = {};

% the plane each slice's current proposal was made at (u and U update it), so
% "plane changed" means scrolled since the last proposal, not since the first
gui_data.proposed_plane = nan(gui_data.Nslices, 1);
proposal_fn = fullfile(gui_data.save_path, 'auto_proposal_controlpoints.mat');
% the proposal, onto the slices still empty
gui_data = load_proposal(proposal_fn, gui_data);

guidata(gui_fig, gui_data);

% the GUI draws slice 1's overlay only at the first key press, so a proposal's
% first slice would look accepted: with anchors or a proposal in play draw it now;
% on a brain with no points the title bar says how to start; others open as usual
if any(~isnan(gui_data.plane_anchors)) || any(gui_data.provisional)
    gui.update_slice(gui_fig);
elseif all(cellfun(@isempty, gui_data.histology_control_points))
    gui.update_window_title(gui_fig);
end
end

function gui_data = load_proposal(proposal_fn, gui_data)
% The proposal loaded onto every slice still empty, provisional, with its ? flags;
% an accepted slice keeps the flags of its unchanged points.

if exist(proposal_fn, 'file')
    P = load(proposal_fn);
    info_fn = fullfile(gui_data.save_path, 'auto_proposal_info.mat');
    low = cell(gui_data.Nslices, 1);
    if exist(info_fn, 'file')
        I = load(info_fn, 'low_confidence');
        low(1:numel(I.low_confidence)) = I.low_confidence(:);
    end

    % a slice re-proposed in an earlier session (u, U) has its latest proposal in
    % annotation_provenance.mat: that is what was accepted, so flags and
    % provenance refer to it, not to the 'autoannotate' file
    prov_fn = fullfile(gui_data.save_path, 'annotation_provenance.mat');
    if exist(prov_fn, 'file')
        L = load(prov_fn);
        if isfield(L, 'latest_histology')
            for k = 1:min(numel(L.latest_histology), numel(P.histology_control_points))
                if ~isempty(L.latest_histology{k})
                    P.histology_control_points{k} = L.latest_histology{k};
                    P.atlas_control_points{k}     = L.latest_atlas{k};
                    low{k} = L.latest_low{k};
                end
            end
        end
    end

    % fill the empty slices; an accepted slice keeps the flags of its unchanged points
    n_filled = 0;
    for k = 1:min(gui_data.Nslices, numel(P.histology_control_points))
        if isempty(gui_data.histology_control_points{k}) ...
                && isempty(gui_data.atlas_control_points{k}) ...
                && ~isempty(P.histology_control_points{k})
            gui_data.histology_control_points{k} = P.histology_control_points{k};
            gui_data.atlas_control_points{k}     = P.atlas_control_points{k};
            gui_data.provisional(k) = true;
            gui_data.proposed_plane(k) = P.atlas_control_points{k}(1, 1);
            if k <= numel(low)
                gui_data.uncertain{k} = logical(low{k}(:));
            end
            n_filled = n_filled + 1;
        elseif ~isempty(gui_data.histology_control_points{k}) && k <= numel(low) && ...
                ~isempty(P.histology_control_points{k})
            % an accepted slice from an earlier session: its unchanged proposal
            % points keep their ?, a point moved since cannot be recognised
            h  = gui_data.histology_control_points{k};
            ph = P.histology_control_points{k};
            lk = logical(low{k}(:));
            f  = false(size(h, 1), 1);
            for r = 1:size(h, 1)
                hit = find(all(abs(ph(:, 1:3) - h(r, 1:3)) < 1e-6, 2), 1);
                if ~isempty(hit) && hit <= numel(lk)
                    f(r) = lk(hit);
                end
            end
            gui_data.uncertain{k} = f;
        end
    end

    % the latest proposal per slice, which u and U update
    gui_data.proposal     = P;
    gui_data.proposal_low = low;
    fprintf(['Loaded the automatic proposal onto %d slice(s), provisional (orange).\n' ...
             '  k accepts a slice as it is, u re-proposes it at the plane on screen,\n' ...
             '  and only accepted or touched slices are saved.\n'], n_filled);
end
end

function keys = controls_lines(annot)
% The automatic annotation's keys, on top of the right column of the GUI's
% controls window.

keys = { ...
    'AUTOMATIC ANNOTATION  (two GUI sessions)', ...
    '', ...
    sprintf('Session 1, new brain: fix %d atlas planes', annot.suggested_anchors), ...
    '  j              go to next suggested slice', ...
    '  wheel          find its atlas plane', ...
    '  a              fix that plane (again: undo)', ...
    '  s              save, close the window', ...
    '  then run run_register_to_atlas with run_mode ''autoannotate''', ...
    '', ...
    'Session 2: review, every slice starts orange', ...
    '  k              looks right: accept, next', ...
    '  K              accept ALL orange slices (asks)', ...
    '  e + drag       a point is off: move it', ...
    '  wheel, then u  wrong plane: re-propose here', ...
    '  a, then U      anchor wrong: fix it, re-propose', ...
    '                 all orange slices (k first the', ...
    '                 ones to keep)', ...
    '  n?             least sure: look there first', ...
    '  s              save; orange slices are NOT', ...
    '                 saved', ...
    ''};
end

function on_key(gui_fig, gui_data, eventdata, gui)
% Event 'key': the automatic annotation's keys. The GUI's functions are bound
% to their own names, so this reads as it did inside the GUI.

update_slice        = gui.update_slice;
update_window_title = gui.update_window_title;
step_slice          = gui.step_slice;

switch eventdata.Key

    % a: fix this slice's atlas plane, as shown, as an anchor for the automatic
    % annotation, or remove it when pressed again on the same plane; a slice with
    % points of its own already has its plane, so anchors are for the others
    case 'a'
        sl    = gui_data.curr_slice;
        plane = round(gui_data.atlas_slice);
        if ~isempty(gui_data.atlas_control_points{sl}) && ~gui_data.provisional(sl)
            disp('This slice has points, which already fix its plane.');
        elseif gui_data.plane_anchors(sl) == plane
            gui_data.plane_anchors(sl) = nan;
            fprintf('Anchor removed from slice %d.\n', sl);
        else
            gui_data.plane_anchors(sl) = plane;
            fprintf('Slice %d anchored at atlas plane %d (%d anchor(s) set).\n', ...
                sl, plane, nnz(~isnan(gui_data.plane_anchors)));
        end
        guidata(gui_fig, gui_data);
        update_window_title(gui_fig);
        update_slice(gui_fig);

    % j: jump to the next suggested anchor slice (first, last, two between)
    case 'j'
        sug  = gui_data.suggested_anchors;
        nxt  = sug(find(sug > gui_data.curr_slice, 1));
        if isempty(nxt)
            nxt = sug(1);
        end
        gui_data.curr_slice = nxt;
        gui_data.sel_side = '';
        gui_data.sel_idx  = 0;
        guidata(gui_fig, gui_data);
        update_slice(gui_fig);

    % k: keep this slice's proposed points as they are, then move on
    case 'k'
        accept_on_key(gui_data, eventdata, gui_fig, update_slice, step_slice, ...
            update_window_title);

    % u: ask the automatic annotation for this slice at the plane on screen, after
    % scrolling a proposed slice to a better plane or on an empty one; never
    % overwrites points placed or accepted by hand
    case 'u'
        propose_on_key(gui_data, eventdata, gui_fig, update_slice);

end
end

function accept_on_key(gui_data, eventdata, gui_fig, update_slice, step_slice, ...
    update_window_title)
% Keys k and K: accept this orange slice and move on, or every orange slice
% after a confirmation.

sl = gui_data.curr_slice;
if any(strcmp(eventdata.Modifier, 'shift'))
    % K: accept every orange slice as it stands, after a confirmation, for
    % a proposal already checked by eye or trusted as it is
    prov = find(gui_data.provisional(:)' & ...
                ~cellfun(@isempty, gui_data.histology_control_points(:)'));
    if isempty(prov)
        disp('No orange slice left to accept.');
        return
    end
    answer = questdlg(sprintf('Accept all %d orange slice(s) as proposed?', ...
        numel(prov)), 'Accept all', 'Accept all', 'Cancel', 'Cancel');
    if strcmp(answer, 'Accept all')
        gui_data.provisional(prov) = false;
        fprintf('Accepted %d slice(s) as proposed: %s. Save with s.\n', ...
            numel(prov), mat2str(prov));
        guidata(gui_fig, gui_data);
        update_slice(gui_fig);
    end
    return
end
if gui_data.provisional(sl) && ~isempty(gui_data.histology_control_points{sl})
    gui_data.provisional(sl) = false;
    fprintf('Slice %d accepted (%d point(s)); %d provisional slice(s) left.\n', ...
        sl, size(gui_data.histology_control_points{sl}, 1), ...
        nnz(gui_data.provisional));
    gui_data = step_slice(gui_data, +1);
    guidata(gui_fig, gui_data);
    update_window_title(gui_fig);
    update_slice(gui_fig);
else
    disp('Nothing provisional to accept on this slice.');
end
end

function propose_on_key(gui_data, eventdata, gui_fig, update_slice)
% Keys u and U: re-propose this slice at the plane on screen, or every orange
% slice at the plane the anchors and the accepted slices give it.

sl = gui_data.curr_slice;
if any(strcmp(eventdata.Modifier, 'shift'))
    % U: re-propose every orange slice, each at the plane the anchors and
    % the accepted slices give it, so a corrected anchor takes effect without
    % leaving the GUI; accepted slices are never touched, so an orange slice
    % worth keeping is accepted (k) first
    prov = find(gui_data.provisional(:))';
    if isempty(prov)
        disp('No orange slice left to re-propose.');
        return
    end

    % the known planes: the accepted slices', else the anchors'
    known = known_planes(gui_data);
    if nnz(~isnan(known)) < 2
        disp('Need at least two anchors or accepted slices to place the planes.');
        return
    end
    planes = predict_planes(gui_data, known);
    fprintf('Re-proposing %d orange slice(s) from %d anchor(s) / accepted slice(s)...\n', ...
        numel(prov), nnz(~isnan(known)));

    % the engine reads the atlas from auto_atlas_planes.mat
    if ~exist(fullfile(gui_data.save_path, 'auto_atlas_planes.mat'), 'file')
        write_atlas_planes(gui_data);
    end
    try
        out = auto_annotate('sections', gui_data.save_path, prov, planes(prov));
    catch err
        out = struct('ok', false, 'message', err.message);
    end
    if ~out.ok
        fprintf('Proposal failed: %s\n', out.message);
        return
    end

    % the new points, still orange, as each slice's proposal
    gui_data = store_proposals(prov, out, gui_data, planes);
    gui_data.sel_side = '';
    gui_data.sel_idx  = 0;
    fprintf('Done: %d slice(s) re-proposed, still orange.\n', numel(prov));
    guidata(gui_fig, gui_data);
    update_slice(gui_fig);
    return
end
% u: this slice, at the plane on screen
propose_slice(gui_data, sl, gui_fig, update_slice);
end

function known = known_planes(gui_data)
% The planes known per slice: an accepted slice's, else its anchor's, else NaN.

% the known planes: the accepted slices', else the anchors'
known = nan(gui_data.Nslices, 1);
for k = 1:gui_data.Nslices
    if ~gui_data.provisional(k) && ~isempty(gui_data.atlas_control_points{k})
        known(k) = gui_data.atlas_control_points{k}(1, 1);
    end
end
use_anchor = isnan(known) & ~isnan(gui_data.plane_anchors);
known(use_anchor) = gui_data.plane_anchors(use_anchor);
end

function gui_data = store_proposals(prov, out, gui_data, planes)
% The new points of every re-proposed slice, still orange, as its proposal.

% the new points, still orange, as each slice's proposal
for j = 1:numel(prov)
    k = prov(j);
    n = size(out.atlas{j}, 1);
    gui_data.histology_control_points{k} = [repmat(k, n, 1), out.hist{j}, ...
        zeros(n, 1)];
    gui_data.atlas_control_points{k} = [repmat(planes(k), n, 1), out.atlas{j}, ...
        zeros(n, 1)];
    gui_data.uncertain{k} = logical(out.low{j}(:));
    gui_data.proposed_plane(k) = planes(k);
    gui_data = remember_proposal(gui_data, k);
end
end

function propose_slice(gui_data, sl, gui_fig, update_slice)
% Key u: this slice proposed at the plane on screen, unless it has points placed
% or accepted by hand.

if ~isempty(gui_data.histology_control_points{sl}) && ~gui_data.provisional(sl)
    disp('This slice has hand-placed or accepted points. Press c to clear them first.');
else
    plane = round(gui_data.atlas_slice);
    fprintf('Proposing slice %d at atlas plane %d...\n', sl, plane);

    % the engine reads the atlas from auto_atlas_planes.mat
    if ~exist(fullfile(gui_data.save_path, 'auto_atlas_planes.mat'), 'file')
        write_atlas_planes(gui_data);
    end
    try
        out = auto_annotate('section', gui_data.save_path, sl, plane);
    catch err
        out = struct('ok', false, 'message', err.message);
    end
    if out.ok && ~isempty(out.atlas)
        n = size(out.atlas, 1);
        gui_data.histology_control_points{sl} = [repmat(sl, n, 1), out.hist, ...
            zeros(n, 1)];
        gui_data.atlas_control_points{sl} = [repmat(plane, n, 1), out.atlas, ...
            zeros(n, 1)];
        gui_data.provisional(sl) = true;
        gui_data.uncertain{sl}   = out.low(:);
        gui_data.proposed_plane(sl) = plane;
        gui_data = remember_proposal(gui_data, sl);
        gui_data.sel_side = '';
        gui_data.sel_idx  = 0;
        fprintf('Proposed %d point(s), %d marked ?. k to accept.\n', n, nnz(out.low));
        guidata(gui_fig, gui_data);
        update_slice(gui_fig);
    elseif out.ok
        disp('No landmarks found on this plane.');
    else
        fprintf('Proposal failed: %s\n', out.message);
    end
end
end

function planes = add_anchor_planes(gui_data, planes)
% Event 'planes': an anchor counts like a slice with points, which wins where it
% has them; with no anchors this is the GUI's own rule.

anchored_only         = isnan(planes) & ~isnan(gui_data.plane_anchors);
planes(anchored_only) = gui_data.plane_anchors(anchored_only);
end

function title_lines = slice_title(gui_data, title_lines)
% Event 'title': the slice's stage in the review on the first line; on the
% second, the anchor of a slice without points, or a request for one.

% the stage in words, since after k the view moves on and colour is easy to misread
if ~isempty(gui_data.proposal) ...
        && ~isempty(gui_data.histology_control_points{gui_data.curr_slice})
    if gui_data.provisional(gui_data.curr_slice)
        title_lines{1} = sprintf('%s  PROPOSED (orange, not saved: k accepts)', title_lines{1});
    else
        title_lines{1} = sprintf('%s  ACCEPTED', title_lines{1});
    end
end

anchorpts = gui_data.atlas_control_points{gui_data.curr_slice};
sl_now    = gui_data.curr_slice;
if ~isempty(anchorpts)
    % a proposed slice scrolled away from its plane: the points came along, but
    % they were made for the other plane
    if gui_data.provisional(sl_now) && ~isnan(gui_data.proposed_plane(sl_now)) && ...
            anchorpts(1,1) ~= gui_data.proposed_plane(sl_now)
        title_lines{2} = sprintf('%s -- plane changed, press u to re-propose here', title_lines{2});
    end
elseif ~isnan(gui_data.plane_anchors(sl_now))
    title_lines{2} = sprintf('PLANE ANCHOR at atlas %2.2f h-slice widths (a to move it here / remove)', ...
        gui_data.plane_anchors(sl_now)/gui_data.slicewidth);
elseif ismember(sl_now, gui_data.suggested_anchors) && isfield(gui_data, 'atlasindsuse')
    title_lines{2} = sprintf('SUGGESTED ANCHOR: scroll to the right plane, press a (now %2.2f)', ...
        gui_data.atlasindsuse(sl_now)/gui_data.slicewidth);
end
end

function mark = flag_points(gui_data, mark)
% Event 'labels': a ? after the number of each point the proposal was least
% sure of, red while provisional and yellow once accepted.

% least sure: the automatic annotation's confidence in the lowest quarter. The
% flag stays after the slice is accepted or edited, since reviewing a point should
% not make its flag vanish. Flags are per row: d removes its row's flag
% (follow_edit), and a point added or undone at the end is matched by padding or
% truncating here.
sl    = gui_data.curr_slice;
flags = logical(gui_data.uncertain{sl}(:));
unc   = false(numel(mark.rows), 1);
m     = min(numel(flags), numel(unc));
unc(1:m) = flags(1:m);
if gui_data.provisional(sl)
    qcol = [1 0.35 0.35];
else
    qcol = [1 0.85 0.1];
end
mark.rows  = unc;
mark.text  = '?';
mark.color = qcol;
end

function modes = window_modes(gui_data, modes)
% Event 'window': the next step of the automatic annotation, once it is in play
% (anchors, a proposal) or on a brain with no points, so a brain annotated by
% hand looks as it did.

n_anch  = nnz(~isnan(gui_data.plane_anchors));
n_sugg  = nnz(~isnan(gui_data.plane_anchors(gui_data.suggested_anchors)));
n_rev   = nnz(gui_data.provisional);
no_pts  = all(cellfun(@isempty, gui_data.histology_control_points));
if ~isempty(gui_data.proposal) && n_rev > 0
    modes{end+1} = sprintf('STEP 3 REVIEW: %d slice(s) left, k accepts, u re-proposes', n_rev);
elseif ~isempty(gui_data.proposal)
    modes{end+1} = 'REVIEW DONE: s to save, then run_register_to_atlas ''register''';
elseif n_anch > 0 && n_sugg < numel(gui_data.suggested_anchors)
    modes{end+1} = sprintf('STEP 1 ANCHORS: %d of %d suggested set (j, scroll, a)', ...
        n_sugg, numel(gui_data.suggested_anchors));
elseif n_anch > 0
    modes{end+1} = sprintf('STEP 1 DONE (%d anchors): s to save, then run_register_to_atlas ''autoannotate''', n_anch);
elseif no_pts
    modes{end+1} = 'NO POINTS YET: click to annotate, or j + a to set anchors for the automatic one';
end
end

function gui_data = follow_edit(gui_data, info)
% Event 'edit': the ? flags are per row too, so the row d removed goes from
% them, and c or t (every point of the slice replaced) clears them.

if isempty(info.row)
    gui_data.uncertain{info.slice} = [];
elseif numel(gui_data.uncertain{info.slice}) >= info.row
    gui_data.uncertain{info.slice}(info.row) = [];
end
end

function save_auto_files(gui_fig, hpts, apts, update_slice)
% Event 'save': write what the automatic annotation needs next to the
% annotation (update_slice is the GUI's).

% the files, none ending in tform.mat, which registerSlicesToAtlas globs for:
%   plane_anchors.mat          the anchors, and the plane the interpolation gives
%                              every slice, on which the automatic step proposes,
%                              so the GUI and it cannot disagree
%   auto_atlas_planes.mat      the warped atlas as this window draws it
%   annotation_provenance.mat  only when a proposal was loaded: per slice, how
%                              many saved pairs are the proposal's, unchanged
gui_data = guidata(gui_fig);
anchor_fn = fullfile(gui_data.save_path, 'plane_anchors.mat');
has_anchor = ~isnan(gui_data.plane_anchors);
if any(has_anchor) || exist(anchor_fn, 'file')
    % the interpolation is refreshed on every redraw: draw first, so the planes
    % written below include the latest anchor
    update_slice(gui_fig);
    gui_data = guidata(gui_fig);
    anchor_slices = find(has_anchor);
    anchor_planes = gui_data.plane_anchors(has_anchor);
    planes        = round(gui_data.atlasindsuse(:));
    save(anchor_fn, 'anchor_slices', 'anchor_planes', 'planes');
    fprintf('Saved %d plane anchor(s) to %s\n', numel(anchor_slices), anchor_fn);
    write_atlas_planes(gui_data);
end

if ~isempty(gui_data.proposal)
    % per slice, the saved pairs and those that are a proposal pair, unchanged
    P = gui_data.proposal;
    n_points = zeros(gui_data.Nslices, 1);
    n_auto_unchanged = zeros(gui_data.Nslices, 1);
    for k = 1:gui_data.Nslices
        h = hpts{k};
        a = apts{k};
        n_points(k) = size(h, 1);
        if isempty(h) || k > numel(P.histology_control_points) ...
                || isempty(P.histology_control_points{k})
            continue
        end
        ph = P.histology_control_points{k};
        pa = P.atlas_control_points{k};
        for r = 1:size(h, 1)
            same = all(abs(ph(:, 1:3) - h(r, 1:3)) < 1e-6, 2) ...
                & all(abs(pa(:, 1:3) - a(r, 1:3)) < 1e-6, 2);
            n_auto_unchanged(k) = n_auto_unchanged(k) + any(same);
        end
    end
    model_version = '';
    info_fn = fullfile(gui_data.save_path, 'auto_proposal_info.mat');
    if exist(info_fn, 'file')
        I = load(info_fn, 'model_version');
        model_version = strtrim(char(I.model_version));
    end
    saved = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));

    % the latest proposal of every slice (after any u or U), so the next session
    % compares and restores the ? flags against what was shown
    latest_histology = P.histology_control_points;
    latest_atlas     = P.atlas_control_points;
    latest_low       = gui_data.proposal_low;
    save(fullfile(gui_data.save_path, 'annotation_provenance.mat'), ...
        'n_points', 'n_auto_unchanged', 'model_version', 'saved', ...
        'latest_histology', 'latest_atlas', 'latest_low');
    fprintf('Provenance: %d of %d saved pair(s) are the automatic proposal, unchanged.\n', ...
        sum(n_auto_unchanged), sum(n_points));
end
end

% ===== Local functions: helpers of the keys and the save =====

function gui_data = remember_proposal(gui_data, k)
% Keep slice k's newest proposal (from u or U) as its proposal: provenance and
% the ? flags of a later session refer to what was shown and accepted.

if isempty(gui_data.proposal)
    gui_data.proposal = struct( ...
        'histology_control_points', {cell(gui_data.Nslices, 1)}, ...
        'atlas_control_points', {cell(gui_data.Nslices, 1)});
    gui_data.proposal_low = cell(gui_data.Nslices, 1);
end
gui_data.proposal.histology_control_points{k} = gui_data.histology_control_points{k};
gui_data.proposal.atlas_control_points{k}     = gui_data.atlas_control_points{k};
gui_data.proposal_low{k} = uint8(gui_data.uncertain{k}(:));
end

function planes = predict_planes(gui_data, known)
% A plane for every slice from the slices whose plane is known (NaN elsewhere),
% by update_slice's rule; used by U, where only anchors and accepted slices count.

% linear between the known slices, and beyond the outermost ones at the step the
% section spacing gives
kn   = find(~isnan(known));
allsl = (1:gui_data.Nslices)';
step = (gui_data.atlasinds(end) - gui_data.atlasinds(1)) / max(1, gui_data.Nslices - 1);
planes = nan(gui_data.Nslices, 1);
inside = allsl >= kn(1) & allsl <= kn(end);
planes(inside) = interp1(kn, known(kn), allsl(inside), 'linear');
before = allsl < kn(1);
planes(before) = known(kn(1))   + (allsl(before) - kn(1))   * step;
after  = allsl > kn(end);
planes(after)  = known(kn(end)) + (allsl(after)  - kn(end)) * step;
planes = min(max(round(planes), 1), size(gui_data.tv, 1));
end

function write_atlas_planes(gui_data)
% Save the warped atlas this window draws, for the automatic annotation, when
% it is missing or older than the rigid alignment or the cutting angle.

% rewritten only when stale, since it is tens of MB
fn  = fullfile(gui_data.save_path, 'auto_atlas_planes.mat');
src = {fullfile(gui_data.save_path, 'regopts.mat'), ...
    fullfile(gui_data.save_path, 'cutting_angle_data.mat')};
stale = ~exist(fn, 'file');
if ~stale
    d_out = dir(fn);
    for k = 1:numel(src)
        d_src = dir(src{k});
        if ~isempty(d_src) && d_src.datenum > d_out.datenum
            stale = true;
        end
    end
end
if stale
    % planes x H x W, uint8
    tv = gui_data.tv;
    save(fn, 'tv', '-v7');
    fprintf('Saved the atlas as drawn here to %s\n', fn);
end
end
