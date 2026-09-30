clear all
close all
clc

% /// Pipeline script #4bis: carry the SEP channel into registered space ///
%
% P4 registers five channels -- DAPI, NANO, AUTO, DIFF, MASK. The microscope
% recorded a sixth, on the green filter it names EGFP, and that one stops at
% volume_centered. It is not a second label: every mouse in this project is a
% SEP-GluA1 knock-in, so the green channel is the tagged receptor itself. Ex
% vivo it reports the whole GluA1 pool -- the pH sensitivity that makes SEP
% surface-specific in a living cell is gone in fixed, permeabilised tissue --
% while the nanobody stain reports the receptors that sat on the membrane.
% Normalising nano by it therefore asks a well-posed question, surface per unit
% receptor expressed, where nano per unit autofluorescence only asks surface
% per unit tissue. Sami expects to trust it more than the autofluorescence, and
% the only honest way to choose is to carry both to the end on every brain and
% look. So this script adds SEP as an EXTRA registered volume and overwrites
% nothing:
%
%   in   <mouse>\lightsuite\volume_centered\chan04_EGFP.tiff   raw, per slice
%                                                              (EGFP names the
%                                                              filter, SEP the
%                                                              molecule)
%   out  <mouse>\lightsuite\volume_registered_sep\
%          chan01_DAPI.tiff   the check channel, see below
%          chan02_SEP.tiff    what the analysis reads
%
% volume_aligned and volume_registered are opened read-only and stay exactly as
% they are, so every number produced so far keeps standing and the choice of
% reference channel becomes a switch in the analysis rather than a fork of the
% data. Nothing in LightSuite is modified either: the registration is not
% redone, it is RE-APPLIED.
%
% Why it can be re-applied at all. Registration happens in two stages:
%
%   centered --( per-slice 2D transform, tformslices )--> aligned
%   aligned  --( per-slice elastix B-spline + affine, then one 3D rigid )--> registered
%
% The second stage is saved in full (transform_params.mat and the elastix_*
% folders), so it costs a transformix pass and no refitting at all. The first
% stage is not saved: alignSliceVolume computes tformslices from the point
% clouds and keeps only the volume it produced. Refitting it would be the wrong
% answer -- the fit subsamples the point clouds at random, so a rerun is never
% bit-for-bit the same, and a reference channel sitting a fraction of a pixel
% off the channel it normalises is exactly the artefact one must not introduce.
% Instead each slice's transform is RECOVERED from the pair of images that fit
% already produced: the centered DAPI slice and the aligned DAPI slice are the
% same picture before and after tformslices(islice). Phase correlation gives a
% first guess on a coarse grid, an affine refinement per level brings it down
% to full resolution, and the result is scored against the very image it had to
% reproduce: on MG897 every slice came back at r = 1.00000 over tissue with no
% residual shift, so the recovered transform is not an approximation of the
% original one, it is the original one. The SEP channel is then warped with it
% and goes through the saved elastix transforms beside the DAPI.
%
% That DAPI is also why the output carries two channels. It rides all the way
% to registered space and is compared there with the DAPI P4 registered months
% ago: if those two agree voxel for voxel, the SEP volume beside it sits in
% the same space as NANO and AUTO, and nano/SEP is a ratio between channels
% that line up. The check is printed per mouse and drawn in
% data\comparisons_v2\processing_diagnostics\sep_channel\<mouse>.png.
% (A second, duller reason: with a single channel in the folder
% loadLargeSliceVolume squeezes the channel dimension away and
% generateRegisteredSliceVolume then reads slices as channels.)
%
% Cost: no refitting, but the per-slice recovery is an image registration and
% the re-application is a transformix pass per slice per channel, so roughly
% 30-45 min per brain. Run it detached and leave it.

%% User-defined parameters

% Where the project lives (worked out from where this file sits, see get_paths).
paths = get_paths();

% Cohort selection. Set mice_to_process to {} for every mouse in the groups.
% A mouse without transform_params.mat has never been registered and is
% skipped with a note rather than an error, so all three groups can be given at
% once and the script picks out the brains that are ready.
groups_to_process = {'young', 'naive', 'rws'};
mice_to_process   = {};

% Nothing here needs the GPU or a shared file, so the cohort splits cleanly
% across a few MATLAB sessions run side by side. P4BIS_MICE, a comma-separated
% list in the environment, wins over the list above when it is set -- an
% environment variable rather than a workspace one because of the clear all at
% the top of this file.
%
%   $env:P4BIS_MICE = 'MG903_SepGluA_P20,MG913_SepGluA_P20'
%   matlab -batch "cd('D:\sep_histology\code'); P4bis_add_sep_channel"
env_mice = getenv('P4BIS_MICE');
if ~isempty(env_mice)
    mice_to_process = strtrim(strsplit(env_mice, ','));
end

% Redo a mouse that already has volume_registered_sep.
overwrite = false;

% The per-slice transform is recovered in a cascade, coarsest grid first, one
% affine refinement per level (same world extent each time, fewer pixels). The
% coarse levels cost seconds and bring the fit to within a pixel; the
% full-resolution level, started from there, needs a single pyramid pass and
% lands on the transform exactly. Dropping the 1 would leave a systematic
% half-pixel offset behind -- measured, not assumed, see min_slice_corr.
recovery_levels = [4 2 1];

% A slice whose recovered transform reproduces the aligned DAPI below this
% correlation is reported by number. With the cascade above every slice of
% MG897 came back at r = 1.00000 and a residual shift of 0 px, so anything
% below 0.99 means something is wrong with that slice, not with the method.
min_slice_corr = 0.99;

%% Add paths

% No atlas is needed here and none is put on the path: every bit of geometry
% comes from transform_params.mat and from the aligned volume already on disk,
% which is what makes this step age-agnostic and safe to run over young and
% adult brains in one go.
addpath(genpath(paths.lightsuite))
addpath(genpath(paths.yaml))
addpath(genpath(paths.elastix))

%% Resolve cohort

get_cohort('verify');
if isempty(mice_to_process)
    cohort = get_cohort('groups', groups_to_process);
else
    cohort = get_cohort('names', mice_to_process);
end

diag_dir = fullfile(paths.data, 'comparisons_v2', 'processing_diagnostics', 'sep_channel');
makeNewDir(diag_dir);

fprintf('P4bis: %d mouse/mice selected.\n', numel(cohort));

%% Loop over mice

for mouse_idx = 1:numel(cohort)

    mouse_name = cohort(mouse_idx).name;
    mouse_type = cohort(mouse_idx).group;
    fprintf('\n=== %s (%s) ===\n', mouse_name, mouse_type);

    mouse_dir      = fullfile(paths.data, mouse_type, mouse_name, 'lightsuite');
    centered_dir   = fullfile(mouse_dir, 'volume_centered');
    aligned_dir    = fullfile(mouse_dir, 'volume_aligned');
    registered_dir = fullfile(mouse_dir, 'volume_registered');
    out_aligned    = fullfile(mouse_dir, 'volume_aligned_sep');
    out_registered = fullfile(mouse_dir, 'volume_registered_sep');
    work_dir       = fullfile(mouse_dir, 'sep_work');

    %% What this mouse has

    tp_name = fullfile(mouse_dir, 'transform_params.mat');
    if ~exist(tp_name, 'file')
        fprintf('  not registered yet (no transform_params.mat) -- skipped.\n');
        continue
    end
    if exist(out_registered, 'dir') && ~overwrite
        fprintf('  volume_registered_sep already there -- skipped (set overwrite = true to redo).\n');
        continue
    end

    files      = [dir(fullfile(centered_dir, '*.tif')); dir(fullfile(centered_dir, '*.tiff'))];
    chan_names = lower({files.name});
    idx_dapi   = find(contains(chan_names, 'dapi'), 1);
    idx_egfp   = find(contains(chan_names, 'egfp'), 1);
    if isempty(idx_dapi) || isempty(idx_egfp)
        error('P4bis: %s has no DAPI and/or green (EGFP) channel in\n  %s', mouse_name, centered_dir);
    end

    aligned_dapi = fullfile(aligned_dir, 'chan01_DAPI.tiff');
    info_aligned = imfinfo(aligned_dapi);

    S         = load(fullfile(mouse_dir, 'sliceinfo.mat'));
    sliceinfo = S.sliceinfo;

    %% Raw slices, in the order the alignment saw them

    % The same three lines alignSliceVolume runs before it aligns anything: the
    % ordering file flips some slices, reorders all of them and drops the ones
    % marked bad. Read it here rather than trusting the slice counts to match.
    fprintf('Loading raw DAPI and SEP... '); tic;
    rawvol    = loadLargeSliceVolume(centered_dir, [idx_dapi idx_egfp]);
    orderfile = fullfile(mouse_dir, 'volume_for_ordering_processing_decisions.txt');
    if exist(orderfile, 'file')
        tabledecisions = readtable(orderfile);
        sliceorder     = tabledecisions.NewOrderOriginalIndex;
        flipsdo        = tabledecisions.FlipState == 1;
        toremove       = tabledecisions.FlipState == -1;
    else
        sliceorder = (1:size(rawvol, 4))';
        flipsdo    = false(size(sliceorder));
        toremove   = false(size(sliceorder));
    end
    rawvol(:, :, :, flipsdo) = flip(rawvol(:, :, :, flipsdo), 2);
    rawvol                   = rawvol(:, :, :, sliceorder);
    rawvol(:, :, :, toremove(sliceorder)) = [];
    Nslices = size(rawvol, 4);
    fprintf('Done! Took %2.2f s\n', toc);

    if Nslices ~= numel(info_aligned)
        error(['P4bis: %s has %d slices after reordering but volume_aligned has %d.\n' ...
               'The ordering file and the aligned volume disagree, so no slice-to-slice\n' ...
               'correspondence can be assumed. Not touching this mouse.'], ...
               mouse_name, Nslices, numel(info_aligned));
    end

    %% Recover the per-slice transform into aligned space

    % Both frames in the world coordinates getRigidlyAlignedVolume used, so a
    % transform estimated on a coarse grid can be handed down to the next one:
    % the world extent is the same at every level, only the sampling differs.
    pxsamp   = double(sliceinfo.px_process) / double(sliceinfo.px_register);
    Rsample  = imref2d(size(rawvol, [1 2]), pxsamp, pxsamp);
    Raligned = imref2d([info_aligned(1).Height info_aligned(1).Width], pxsamp, pxsamp);

    [optimizer, metric]         = imregconfig('monomodal');
    optimizer.MaximumIterations = 300;

    tforms(Nslices, 1) = affinetform2d;
    corr_slice         = nan(Nslices, 1);
    shift_slice        = nan(Nslices, 1);

    fprintf('Recovering per-slice transforms from the DAPI pair...\n');
    rectic = tic; msg = [];
    for islice = 1:Nslices

        fixed  = single(imread(aligned_dapi, 'Index', islice, 'Info', info_aligned));
        moving = single(rawvol(:, :, 1, islice));

        best = [];
        for kdown = recovery_levels

            if kdown == 1
                fx = fixed; Rf = Raligned;
                mv = moving; Rm = Rsample;
            else
                fx = imresize(fixed,  1/kdown); Rf = imref2d(size(fx), Raligned.XWorldLimits, Raligned.YWorldLimits);
                mv = imresize(moving, 1/kdown); Rm = imref2d(size(mv), Rsample.XWorldLimits,  Rsample.YWorldLimits);
            end

            % The first level needs a starting point and phase correlation
            % provides one without a guess. imregcorr is asked for the
            % transform between the two pixel grids and moved into world
            % coordinates here, because its own composition with these
            % references fails ('Invalid transformation matrix') as soon as a
            % pixel is not one world unit wide, which here it never is.
            if isempty(best)
                tc   = imregcorr(mv, fx, 'similarity');
                best = affinetform2d(world_from_intrinsic(Rf) * tc.A / world_from_intrinsic(Rm));
            end

            % Each level keeps its refinement only if it scores better than
            % what it started from, so a level that wanders costs nothing.
            try
                cand = imregtform(mv, Rm, fx, Rf, 'affine', optimizer, metric, ...
                    'InitialTransformation', best, 'PyramidLevels', 1);
                if image_match(mv, Rm, fx, Rf, cand) > image_match(mv, Rm, fx, Rf, best)
                    best = affinetform2d(cand.A);
                end
            catch
                % a level imregtform cannot handle leaves the fit in hand standing
            end
        end

        tforms(islice) = best;
        [corr_slice(islice), shift_slice(islice)] = recovery_score(moving, Rsample, fixed, Raligned, best);

        fprintf(repmat('\b', 1, numel(msg)));
        msg = sprintf('Slice %d/%d done. Took %2.2f s. Median r %1.5f.\n', islice, Nslices, ...
            toc(rectic), median(corr_slice(1:islice)));
        fprintf(msg);
    end

    bad = find(corr_slice < min_slice_corr)';
    if isempty(bad)
        fprintf('  every slice recovered at r >= %.2f (worst %1.5f, largest residual shift %.2f px).\n', ...
            min_slice_corr, min(corr_slice), max(shift_slice));
    else
        fprintf(2, '  slice(s) %s recovered below r = %.2f (%s). See the diagnostic figure.\n', ...
            mat2str(bad), min_slice_corr, strjoin(compose('%1.3f', corr_slice(bad)'), ', '));
    end

    %% Warp both channels into aligned space

    fprintf('Warping DAPI and SEP into aligned space... '); tic;
    backvalues = recompute_backvalues(rawvol);
    alvol      = zeros([Raligned.ImageSize 2 Nslices], 'uint16');
    for islice = 1:Nslices
        for ichan = 1:2
            alvol(:, :, ichan, islice) = imwarp(rawvol(:, :, ichan, islice), Rsample, tforms(islice), ...
                'linear', 'OutputView', Raligned, 'FillValues', double(backvalues(ichan, islice)));
        end
    end
    saveLargeSliceVolume(alvol, {'DAPI', 'SEP'}, out_aligned);
    fprintf('Done! Took %2.2f s\n', toc);

    %% Re-apply the saved registration

    transformparams = load(tp_name);
    % The elastix folder was written as an absolute path, drive letter
    % included, so a brain registered on the annotation machine points at a
    % letter that may not exist here. The folder itself travels with the mouse.
    if ~exist(transformparams.tformbspline_samp20um_to_atlas_20um, 'dir')
        transformparams.tformbspline_samp20um_to_atlas_20um = fullfile(mouse_dir, 'elastix_reverse');
    end

    % generateRegisteredSliceVolume writes to <procpath>\volume_registered and
    % that path is hardcoded. LightSuite is not to be edited (the adults were
    % registered with it as it stands), so it is pointed at a working folder
    % and the result is moved into place afterwards. sliceinfo is otherwise
    % passed exactly as P4's 'register' branch passes it, stale px_atlas and
    % all, because reproducing that run is the whole point: the SEP volume has
    % to land on the same grid as volume_registered, not on a better one.
    makeNewDir(work_dir);
    si             = sliceinfo;
    si.channames   = {'DAPI', 'SEP'};
    si.procpath    = work_dir;
    si.slicevolfin = out_aligned;
    generateRegisteredSliceVolume(si, transformparams);

    if exist(out_registered, 'dir')
        rmdir(out_registered, 's');
    end
    movefile(fullfile(work_dir, 'volume_registered'), out_registered);
    rmdir(work_dir, 's');

    %% Check it landed where NANO and AUTO are

    % Every twentieth plane of the registered DAPI, this run against P4's.
    fprintf('Checking the registered DAPI against the one P4 wrote... '); tic;
    ref_dapi = fullfile(registered_dir, 'chan01_DAPI.tiff');
    new_dapi = fullfile(out_registered, 'chan01_DAPI.tiff');
    info_ref = imfinfo(ref_dapi);
    info_new = imfinfo(new_dapi);
    if numel(info_ref) ~= numel(info_new) || info_ref(1).Height ~= info_new(1).Height ...
            || info_ref(1).Width ~= info_new(1).Width
        error(['P4bis: %s registered SEP came out %dx%dx%d against volume_registered''s %dx%dx%d.\n' ...
               'The two are not on the same grid and must not be divided into one another.'], ...
               mouse_name, info_new(1).Height, info_new(1).Width, numel(info_new), ...
               info_ref(1).Height, info_ref(1).Width, numel(info_ref));
    end
    planes   = 1:20:numel(info_ref);
    corr_reg = nan(numel(planes), 1);
    for k = 1:numel(planes)
        a = single(imread(ref_dapi, 'Index', planes(k), 'Info', info_ref));
        b = single(imread(new_dapi, 'Index', planes(k), 'Info', info_new));
        m = a > 0 | b > 0;
        if nnz(m) > 100
            corr_reg(k) = corr(a(m), b(m));
        end
    end
    fprintf('Done! Took %2.2f s\n', toc);
    fprintf('  registered DAPI agreement: median r = %1.4f over %d planes (worst %1.4f).\n', ...
        median(corr_reg, 'omitnan'), nnz(~isnan(corr_reg)), min(corr_reg));

    %% Diagnostic figure

    fig = figure('Color', 'w', 'Position', [80 80 1500 850], 'Visible', 'off');
    tl  = tiledlayout(fig, 2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile(tl);
    plot(1:Nslices, corr_slice, 'o-', 'Color', [0.15 0.15 0.15], ...
        'MarkerFaceColor', [0.15 0.15 0.15], 'MarkerSize', 4);
    hold on; yline(min_slice_corr, '--', 'Color', [0.8 0.2 0.1]);
    xlabel('slice (aligned order)'); ylabel('r vs the aligned DAPI, tissue only');
    title(sprintf('per-slice transform recovery (median %1.5f, worst shift %.2f px)', ...
        median(corr_slice), max(shift_slice)));
    ylim([min(0.9, min(corr_slice) - 0.01) 1.001]); box off

    nexttile(tl);
    plot(planes - 1, corr_reg, 'o-', 'Color', [0.15 0.15 0.15], ...
        'MarkerFaceColor', [0.15 0.15 0.15], 'MarkerSize', 4);
    xlabel('registered plane (ML index)'); ylabel('r vs volume\_registered');
    title(sprintf('registered DAPI, this run vs P4 (median %1.4f)', median(corr_reg, 'omitnan')));
    ylim([min(0.9, min(corr_reg) - 0.01) 1.001]); box off

    [~, iworst] = min(corr_slice);
    fixed = single(imread(aligned_dapi, 'Index', iworst, 'Info', info_aligned));
    chk   = imwarp(single(rawvol(:, :, 1, iworst)), Rsample, tforms(iworst), 'linear', ...
                   'OutputView', Raligned, 'FillValues', 0);
    lim   = [0 quantile(fixed(fixed > 0), 0.999)];
    nexttile(tl); imagesc(fixed, lim); axis image off; colormap(gca, gray);
    title(sprintf('aligned DAPI, slice %d (the worst one)', iworst));
    nexttile(tl); imagesc(chk, lim); axis image off; colormap(gca, gray);
    title(sprintf('recovered warp of the raw DAPI (r = %1.4f)', corr_slice(iworst)));
    nexttile(tl); imagesc(abs(chk - fixed), [0 diff(lim) / 4]); axis image off; colormap(gca, hot);
    title('difference');

    mid      = round(numel(info_new) / 2);
    sep_mid = single(imread(fullfile(out_registered, 'chan02_SEP.tiff'), 'Index', mid));
    nexttile(tl); imagesc(sep_mid, [0 quantile(sep_mid(sep_mid > 0), 0.999)]);
    axis image off; colormap(gca, hot);
    title(sprintf('registered SEP, ML plane %d', mid));

    title(tl, sprintf('%s -- SEP carried into registered space with the saved transforms', ...
        strrep(mouse_name, '_', ' ')), 'FontWeight', 'bold');
    exportgraphics(fig, fullfile(diag_dir, sprintf('%s.png', mouse_name)), 'Resolution', 150);
    close(fig);

    % The numbers behind the figure, for anyone auditing this without MATLAB.
    fid = fopen(fullfile(diag_dir, sprintf('%s.txt', mouse_name)), 'w');
    fprintf(fid, 'mouse            %s (%s)\n', mouse_name, mouse_type);
    fprintf(fid, 'written          %s\n', datetime('now', 'Format', 'yyyy-MM-dd HH:mm'));
    fprintf(fid, 'slices           %d\n', Nslices);
    fprintf(fid, 'slice recovery   median r %1.5f, min %1.5f, below %.2f: %s\n', ...
        median(corr_slice), min(corr_slice), min_slice_corr, mat2str(bad));
    fprintf(fid, 'residual shift   median %.3f px, max %.3f px (%.2f um at %g um/px)\n', ...
        median(shift_slice), max(shift_slice), max(shift_slice) * double(sliceinfo.px_process), ...
        double(sliceinfo.px_process));
    fprintf(fid, 'registered DAPI  median r %1.5f, min %1.5f over %d planes\n', ...
        median(corr_reg, 'omitnan'), min(corr_reg), nnz(~isnan(corr_reg)));
    fprintf(fid, 'output           %s\n', out_registered);
    fclose(fid);

    fprintf('  wrote %s\n', out_registered);

    clear rawvol alvol tforms

end

fprintf('\nP4bis: done.\n');

%% Local functions

function M = world_from_intrinsic(R)
% The plain definition of an imref2d's intrinsic-to-world map: pixel centres
% sit at 1..N, the first one half a pixel inside the world limit.
M = [R.PixelExtentInWorldX 0 R.XWorldLimits(1) - 0.5 * R.PixelExtentInWorldX; ...
     0 R.PixelExtentInWorldY R.YWorldLimits(1) - 0.5 * R.PixelExtentInWorldY; ...
     0 0 1];
end

function r = image_match(moving, Rmoving, fixed, Rfixed, tform)
% Correlation between the fixed image and the moving one warped onto it. Used
% to choose between two candidate transforms at one level of the cascade, where
% only the comparison matters.
warped = imwarp(moving, Rmoving, tform, 'linear', 'OutputView', Rfixed, 'FillValues', 0);
m      = fixed > 0;
if nnz(m) < 100
    r = NaN;
    return
end
r = corr(double(warped(m)), double(fixed(m)));
end

function [r, shiftpx] = recovery_score(moving, Rmoving, fixed, Rfixed, tform)
% How faithfully the recovered transform reproduces the slice the original
% alignment saved. Two numbers, because they fail differently: a correlation
% over tissue only -- the background is a constant fill and would flatter it --
% and the residual translation between the two images, which is what a
% misplaced channel would actually look like, in pixels of the aligned grid.
warped = imwarp(moving, Rmoving, tform, 'linear', 'OutputView', Rfixed, 'FillValues', NaN);
bg     = quantile(fixed(fixed > 0), 0.01);
ok     = (fixed > 2 * bg) & isfinite(warped);
if nnz(ok) < 100
    r = NaN; shiftpx = NaN;
    return
end
r = corr(double(warped(ok)), double(fixed(ok)));

a = fixed;  a(~isfinite(warped)) = 0;
b = warped; b(~isfinite(b))      = 0;
sh      = imregcorr(b, a, 'translation');
shiftpx = hypot(sh.Translation(1), sh.Translation(2));
end
