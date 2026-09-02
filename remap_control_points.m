function remap_control_points(mousename, old_decisions_file, do_apply)
%REMAP_CONTROL_POINTS Carry saved control points across a slice reorder.
%
% Control points live in atlas2histology_tform.mat as one cell per slice,
% indexed by POSITION in the ordered volume rather than by the piece of tissue
% they were placed on. Reordering in P1bis therefore leaves every point behind
% on whatever slice now occupies its old position, and nothing complains.
%
% This walks them across. It reads the ordering decisions as they were when
% the points went down and as they are now, matches positions through the
% original slice index that both refer to, and rewrites the cell arrays in the
% new order. The atlas plane rides along in column 1, so an anchored slice
% stays anchored where it was.
%
% Points on a slice whose FLIP state changed are dropped rather than moved: a
% flip mirrors the image, so the coordinates no longer land on the same tissue
% and that slice has to be annotated again. Slices that were dropped from the
% volume, or that are new to it, come out empty.
%
%   remap_control_points(mouse, old_decisions)          % dry run, prints only
%   remap_control_points(mouse, old_decisions, true)    % write the new file
%
% The old decisions file is the backup taken BEFORE curating in P1bis. Without
% it there is nothing to match against, which is why it has to be kept.

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

S = load(tformfile);
hold_pts = S.histology_control_points;
aold_pts = S.atlas_control_points;

fprintf('%s\n', repmat('=', 1, 72));
fprintf('%s\n', mousename);
fprintf('  old order: %d slices in the volume\n', numel(seqold));
fprintf('  new order: %d slices in the volume\n', numel(seqnew));
fprintf('  saved points: %d slices\n', numel(hold_pts));

if numel(hold_pts) ~= numel(seqold)
    warning(['The saved points cover %d slices but the old decisions file ' ...
             'describes %d. The backup may not be the one that was in force ' ...
             'when these points were placed.'], numel(hold_pts), numel(seqold));
end

%% Walk the points across

hnew_pts = repmat({zeros(0,4)}, numel(seqnew), 1);
anew_pts = repmat({zeros(0,4)}, numel(seqnew), 1);

nmoved = 0; nsame = 0; ndropped = 0; nlost = 0;

fprintf('\n%-6s %-9s %-9s %6s   %s\n', 'new', 'original', 'was at', 'Npts', 'what happens');
fprintf('%s\n', repmat('-', 1, 72));

for q = 1:numel(seqnew)

    orig = seqnew(q);
    p    = find(seqold == orig, 1);

    if isempty(p) || p > numel(hold_pts)
        fprintf('%-6d %-9d %-9s %6s   new to the volume, starts empty\n', q, orig, '-', '-');
        continue
    end

    npts = size(hold_pts{p}, 1);

    if flipnew(orig) ~= flipold(orig)
        if npts > 0
            ndropped = ndropped + 1;
            fprintf('%-6d %-9d %-9d %6d   FLIP CHANGED, points dropped, re-annotate\n', ...
                q, orig, p, npts);
        end
        continue
    end

    hnew_pts{q} = hold_pts{p};
    anew_pts{q} = aold_pts{p};

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

% Anything annotated in the old volume that has no home in the new one
for p = 1:min(numel(hold_pts), numel(seqold))
    if ~isempty(hold_pts{p}) && ~ismember(seqold(p), seqnew)
        nlost = nlost + 1;
        fprintf('%-6s %-9d %-9d %6d   slice dropped from the volume, points lost\n', ...
            '-', seqold(p), p, size(hold_pts{p},1));
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

histology_control_points = hnew_pts;
atlas_control_points     = anew_pts;
save(tformfile, 'atlas_control_points', 'histology_control_points');

fprintf('\nWritten: %s\n', tformfile);
fprintf('Previous version kept at: %s\n', backup);

end


function [seq, flipstate] = ordered_sequence(T)
% The original slice indices that survive into the volume, in the order they
% appear there. Reorder first, then drop -- exactly what alignSliceVolume does,
% so the positions here are the positions the GUI counts in.

order     = T.NewOrderOriginalIndex(:);
toremove  = T.FlipState == -1;
seq       = order(~toremove(order));
flipstate = T.FlipState(:);

end
