function remap_control_points(mousename, old_decisions_file, do_apply)
%REMAP_CONTROL_POINTS  Carry saved control points across a slice reorder.
%   REMAP_CONTROL_POINTS(mouse, old_decisions) prints what would move (a dry
%   run); REMAP_CONTROL_POINTS(mouse, old_decisions, true) writes the new file.
%
%   Control points live in atlas2histology_tform.mat as one cell per slice,
%   indexed by position in the ordered volume, not by the piece of tissue they
%   were placed on. A reorder in run_order_slices therefore leaves every point
%   on whatever slice now holds its old position, and nothing complains.
%
%   This reads the ordering decisions as they were when the points were placed
%   and as they are now, matches positions through the original slice index
%   both refer to, and rewrites the cell arrays in the new order. The atlas
%   plane is column 1, so an anchored slice stays anchored where it was.
%
%   Points on a slice whose flip state changed are dropped, not moved: a flip
%   mirrors the image, so the coordinates no longer land on the same tissue
%   and the slice has to be annotated again. Slices dropped from the volume,
%   or new to it, come out empty.
%
%   old_decisions is the backup of the decisions file taken before curating in
%   run_order_slices; without it there is nothing to match against, which is
%   why it has to be kept. On writing, the previous file is kept beside it as
%   atlas2histology_tform_prereorder_<date>.mat.

if nargin < 3
    do_apply = false;
end

%% Locate the mouse

cohort   = get_cohort('names', {mousename});
procpath = fullfile(cohort(1).base_dir, 'lightsuite');

newfile   = fullfile(procpath, 'volume_for_ordering_processing_decisions.txt');
tformfile = fullfile(procpath, 'atlas2histology_tform.mat');

if ~exist(old_decisions_file, 'file')
    error('Old decisions file not found:\n  %s', old_decisions_file);
end
if ~exist(newfile, 'file')
    error('Current decisions file not found:\n  %s', newfile);
end
if ~exist(tformfile, 'file')
    error('No saved control points to remap:\n  %s', tformfile);
end

%% Read both orderings

Told = readtable(old_decisions_file);
Tnew = readtable(newfile);

[seqold, flipold] = ordered_sequence(Told);
[seqnew, flipnew] = ordered_sequence(Tnew);

% the saved points, histology and atlas, in the old order
S = load(tformfile);
hist_old_pts = S.histology_control_points;
atlas_old_pts = S.atlas_control_points;

fprintf('%s\n', repmat('=', 1, 72));
fprintf('%s\n', mousename);
fprintf('  old order: %d slices in the volume\n', numel(seqold));
fprintf('  new order: %d slices in the volume\n', numel(seqnew));
fprintf('  saved points: %d slices\n', numel(hist_old_pts));

if numel(hist_old_pts) ~= numel(seqold)
    warning(['The saved points cover %d slices but the old decisions file ' ...
             'describes %d. The backup may not be the one that was in force ' ...
             'when these points were placed.'], numel(hist_old_pts), numel(seqold));
end

%% Walk the points across

hist_new_pts = repmat({zeros(0,4)}, numel(seqnew), 1);
atlas_new_pts = repmat({zeros(0,4)}, numel(seqnew), 1);

nmoved = 0;
nsame = 0;
ndropped = 0;
nlost = 0;

fprintf('\n%-6s %-9s %-9s %6s   %s\n', 'new', 'original', 'was at', 'Npts', ...
    'what happens');
fprintf('%s\n', repmat('-', 1, 72));

% each slice of the new order: q its new position, p its old one
for q = 1:numel(seqnew)

    orig = seqnew(q);
    p    = find(seqold == orig, 1);

    if isempty(p) || p > numel(hist_old_pts)
        fprintf('%-6d %-9d %-9s %6s   new to the volume, starts empty\n', ...
            q, orig, '-', '-');
        continue
    end

    npts = size(hist_old_pts{p}, 1);

    % a flipped slice loses its points
    if flipnew(orig) ~= flipold(orig)
        if npts > 0
            ndropped = ndropped + 1;
            fprintf('%-6d %-9d %-9d %6d   FLIP CHANGED, points dropped, re-annotate\n', ...
                q, orig, p, npts);
        end
        continue
    end

    hist_new_pts{q} = hist_old_pts{p};
    atlas_new_pts{q} = atlas_old_pts{p};

    if npts == 0
        continue
    end
    if p == q
        nsame = nsame + 1;
    else
        nmoved = nmoved + 1;
        fprintf('%-6d %-9d %-9d %6d   moved\n', q, orig, p, npts);
    end
end

% anything annotated in the old volume that has no place in the new one
for p = 1:min(numel(hist_old_pts), numel(seqold))
    if ~isempty(hist_old_pts{p}) && ~ismember(seqold(p), seqnew)
        nlost = nlost + 1;
        fprintf('%-6s %-9d %-9d %6d   slice dropped from the volume, points lost\n', ...
            '-', seqold(p), p, size(hist_old_pts{p},1));
    end
end

fprintf('%s\n', repmat('-', 1, 72));
fprintf('%d annotated slices stay put, %d move, %d lose their points to a flip, %d to a removal\n', ...
    nsame, nmoved, ndropped, nlost);

%% Write, or say what would have been written

if ~do_apply
    fprintf('\nDRY RUN. Nothing written. Re-run with do_apply = true to commit.\n');
    return
end

backup = fullfile(procpath, sprintf('atlas2histology_tform_prereorder_%s.mat', ...
    datestr(now, 'yyyymmdd_HHMMSS'))); %#ok<TNOW1,DATST>
copyfile(tformfile, backup);

% the loaded struct goes back whole, not the two arrays alone: some GUI saves put
% an atlas2histology_tform in the file too, which a save of two variables would drop
S.histology_control_points = hist_new_pts;
S.atlas_control_points     = atlas_new_pts;
save(tformfile, '-struct', 'S');

fprintf('\nWritten: %s\n', tformfile);
fprintf('Previous version kept at: %s\n', backup);

end

% ===== Local functions =====

function [seq, flipstate] = ordered_sequence(T)
% The original indices of the slices that stay in the volume, in their order
% there, and the flip state of every original slice.

% reorder first, then drop, as alignSliceVolume does, so the positions here are
% the ones the GUI counts in
order     = T.NewOrderOriginalIndex(:);
toremove  = T.FlipState == -1;
seq       = order(~toremove(order));
flipstate = T.FlipState(:);

end
