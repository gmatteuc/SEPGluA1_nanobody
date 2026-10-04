function add_sep_channel(run_settings)
%ADD_SEP_CHANNEL  Carry the SEP channel of the selected mice into registered space.
%   ADD_SEP_CHANNEL(run_settings) does the work of run_add_sep_channel,
%   which sets the fields of run_settings (paths, groups_to_process,
%   mice_to_process, overwrite, recovery_levels, min_slice_corr) and says
%   what each one does. P4BIS_MICE is read there, before the call.
%
%   No atlas is put on the path: the geometry comes from transform_params.mat
%   and from the aligned volume on disk, so the step does not depend on the
%   mouse's age and runs over young and adult brains in one go. The toolboxes
%   are on the path from sep_setup_paths.

% settings of run_add_sep_channel, under the names the code below uses
paths = run_settings.paths;
groups_to_process = run_settings.groups_to_process;
mice_to_process = run_settings.mice_to_process;
overwrite = run_settings.overwrite;
recovery_levels = run_settings.recovery_levels;
min_slice_corr = run_settings.min_slice_corr;

%% Resolve cohort

% check that the registry still lists the adults in their legacy order
get_cohort('verify');

% the mice named, or else every mouse of the groups
if isempty(mice_to_process)
    cohort = get_cohort('groups', groups_to_process);
else
    cohort = get_cohort('names', mice_to_process);
end

% the folder of the diagnostic figures
diag_dir = fullfile(paths.data, 'comparisons_v2', 'processing_diagnostics', ...
    'sep_channel');
makeNewDir(diag_dir);

fprintf('run_add_sep_channel: %d mouse/mice selected.\n', numel(cohort));

%% Carry the SEP channel of each mouse

for mouse_idx = 1:numel(cohort)

    mouse_name = cohort(mouse_idx).name;
    mouse_type = cohort(mouse_idx).group;
    fprintf('\n=== %s (%s) ===\n', mouse_name, mouse_type);

    % the folders read and written
    mouse_dir = fullfile(paths.data, mouse_type, mouse_name, 'lightsuite');
    centered_dir = fullfile(mouse_dir, 'volume_centered');
    aligned_dir = fullfile(mouse_dir, 'volume_aligned');
    registered_dir = fullfile(mouse_dir, 'volume_registered');
    out_aligned = fullfile(mouse_dir, 'volume_aligned_sep');
    out_registered = fullfile(mouse_dir, 'volume_registered_sep');
    work_dir = fullfile(mouse_dir, 'sep_work');

    %% What this mouse has

    % skip a mouse not registered yet, or done already
    tp_name = fullfile(mouse_dir, 'transform_params.mat');
    if ~exist(tp_name, 'file')
        fprintf('  not registered yet (no transform_params.mat) -- skipped.\n');
        continue
    end
    if exist(out_registered, 'dir') && ~overwrite
        fprintf('  volume_registered_sep already there -- skipped (set overwrite = true to redo).\n');
        continue
    end

    % carry its SEP channel into registered space
    carry_sep_channel(centered_dir, mouse_name, aligned_dir, mouse_dir, ...
        recovery_levels, min_slice_corr, out_aligned, tp_name, work_dir, ...
        out_registered, registered_dir, diag_dir, mouse_type);

end

fprintf('\nrun_add_sep_channel: done.\n');

end

% ===== Local functions =====

function carry_sep_channel(centered_dir, mouse_name, aligned_dir, mouse_dir, ...
    recovery_levels, min_slice_corr, out_aligned, tp_name, work_dir, out_registered, ...
    registered_dir, diag_dir, mouse_type)
% One mouse's SEP channel carried into registered space, checked against its
% registered DAPI, with the diagnostic figure and its numbers.

% the DAPI and green channels among the centered volumes
files = [dir(fullfile(centered_dir, '*.tif')); ...
    dir(fullfile(centered_dir, '*.tiff'))];
chan_names = lower({files.name});
idx_dapi = find(contains(chan_names, 'dapi'), 1);
idx_egfp = find(contains(chan_names, 'egfp'), 1);
if isempty(idx_dapi) || isempty(idx_egfp)
    error('run_add_sep_channel: %s has no DAPI and/or green (EGFP) channel in\n  %s', ...
        mouse_name, centered_dir);
end

aligned_dapi = fullfile(aligned_dir, 'chan01_DAPI.tiff');
info_aligned = imfinfo(aligned_dapi);

S = load(fullfile(mouse_dir, 'sliceinfo.mat'));
sliceinfo = S.sliceinfo;

%% Raw slices, in the order the alignment saw them

[rawvol, Nslices] = load_raw_slices(centered_dir, idx_dapi, idx_egfp, mouse_dir, ...
    info_aligned, mouse_name);

%% Recover the per-slice transform into aligned space

[Rsample, Raligned, tforms, corr_slice, shift_slice, bad] = recover_slice_transforms( ...
    sliceinfo, rawvol, info_aligned, Nslices, aligned_dapi, recovery_levels, ...
    min_slice_corr);

%% Warp both channels into aligned space

warp_to_aligned(rawvol, Raligned, Nslices, Rsample, tforms, out_aligned);

%% Re-apply the saved registration

reapply_registration(tp_name, mouse_dir, work_dir, sliceinfo, out_aligned, ...
    out_registered);

%% Check it landed where NANO and AUTO are

[info_new, planes, corr_reg] = check_registered_dapi(registered_dir, out_registered, ...
    mouse_name);

%% Diagnostic figure

plot_sep_diagnostics(Nslices, corr_slice, min_slice_corr, shift_slice, planes, ...
    corr_reg, aligned_dapi, info_aligned, rawvol, Rsample, tforms, Raligned, info_new, ...
    out_registered, mouse_name, diag_dir);

write_sep_numbers(diag_dir, mouse_name, mouse_type, Nslices, corr_slice, ...
    min_slice_corr, bad, shift_slice, sliceinfo, corr_reg, out_registered);

fprintf('  wrote %s\n', out_registered);

% free memory before the next mouse
clear rawvol alvol tforms

end

function [rawvol, Nslices] = load_raw_slices(centered_dir, idx_dapi, idx_egfp, ...
    mouse_dir, info_aligned, mouse_name)
% The raw DAPI and SEP slices, flipped, reordered and with the bad ones dropped
% as the ordering file says; stops if their number differs from the aligned volume.

% the three steps alignSliceVolume takes before it aligns: the ordering file
% flips some slices, reorders all of them and drops the ones marked bad; read
% here rather than trusting the slice counts to match
fprintf('Loading raw DAPI and SEP... ');
tic;
rawvol = loadLargeSliceVolume(centered_dir, [idx_dapi idx_egfp]);
orderfile = fullfile(mouse_dir, 'volume_for_ordering_processing_decisions.txt');
if exist(orderfile, 'file')
    tabledecisions = readtable(orderfile);
    sliceorder = tabledecisions.NewOrderOriginalIndex;
    flipsdo = tabledecisions.FlipState == 1;
    toremove = tabledecisions.FlipState == -1;
else
    sliceorder = (1:size(rawvol, 4))';
    flipsdo = false(size(sliceorder));
    toremove = false(size(sliceorder));
end
rawvol(:, :, :, flipsdo) = flip(rawvol(:, :, :, flipsdo), 2);
rawvol = rawvol(:, :, :, sliceorder);
rawvol(:, :, :, toremove(sliceorder)) = [];
Nslices = size(rawvol, 4);
fprintf('Done! Took %2.2f s\n', toc);

if Nslices ~= numel(info_aligned)
    error(['run_add_sep_channel: %s has %d slices after reordering but volume_aligned has %d.\n' ...
           'The ordering file and the aligned volume disagree, so no slice-to-slice\n' ...
           'correspondence can be assumed. Not touching this mouse.'], ...
           mouse_name, Nslices, numel(info_aligned));
end

end

function [Rsample, Raligned, tforms, corr_slice, shift_slice, bad] = ...
    recover_slice_transforms(sliceinfo, rawvol, info_aligned, Nslices, aligned_dapi, ...
    recovery_levels, min_slice_corr)
% Each slice's transform into aligned space, recovered from the DAPI pair level
% by level, with its correlation and residual shift; the slices below min_slice_corr.

% both frames in the world coordinates getRigidlyAlignedVolume used: the world
% extent is the same at every level, so a transform fitted on a coarse grid
% can be handed down to the next one
pxsamp = double(sliceinfo.px_process) / double(sliceinfo.px_register);
Rsample = imref2d(size(rawvol, [1 2]), pxsamp, pxsamp);
Raligned = imref2d([info_aligned(1).Height info_aligned(1).Width], pxsamp, pxsamp);

[optimizer, metric] = imregconfig('monomodal');
optimizer.MaximumIterations = 300;

tforms(Nslices, 1) = affinetform2d;
corr_slice = nan(Nslices, 1);
shift_slice = nan(Nslices, 1);

fprintf('Recovering per-slice transforms from the DAPI pair...\n');
rectic = tic;
msg = [];
for islice = 1:Nslices

    fixed = single(imread(aligned_dapi, 'Index', islice, 'Info', info_aligned));
    moving = single(rawvol(:, :, 1, islice));

    % the transform, refined over the levels
    best = fit_slice_transform(recovery_levels, fixed, Raligned, moving, Rsample, ...
        optimizer, metric);

    tforms(islice) = best;
    [corr_slice(islice), shift_slice(islice)] = recovery_score(moving, Rsample, ...
        fixed, Raligned, best);

    % progress, on one line rewritten for each slice
    fprintf(repmat('\b', 1, numel(msg)));
    msg = sprintf('Slice %d/%d done. Took %2.2f s. Median r %1.5f.\n', islice, Nslices, ...
        toc(rectic), median(corr_slice(1:islice)));
    fprintf(msg);
end

% name the slices recovered below min_slice_corr
bad = find(corr_slice < min_slice_corr)';
if isempty(bad)
    fprintf('  every slice recovered at r >= %.2f (worst %1.5f, largest residual shift %.2f px).\n', ...
        min_slice_corr, min(corr_slice), max(shift_slice));
else
    fprintf(2, '  slice(s) %s recovered below r = %.2f (%s). See the diagnostic figure.\n', ...
        mat2str(bad), min_slice_corr, ...
        strjoin(compose('%1.3f', corr_slice(bad)'), ', '));
end

end

function best = fit_slice_transform(recovery_levels, fixed, Raligned, moving, Rsample, ...
    optimizer, metric)
% One slice's transform, refined level by level from a phase-correlation start;
% a level keeps its refinement only if it matches better.

best = [];
for kdown = recovery_levels

    % both images at this level's sampling, over the same world extent
    if kdown == 1
        fx = fixed;
        Rf = Raligned;
        mv = moving;
        Rm = Rsample;
    else
        fx = imresize(fixed,  1/kdown);
        Rf = imref2d(size(fx), Raligned.XWorldLimits, Raligned.YWorldLimits);
        mv = imresize(moving, 1/kdown);
        Rm = imref2d(size(mv), Rsample.XWorldLimits,  Rsample.YWorldLimits);
    end

    % the first level's starting point, from phase correlation, which needs
    % no guess. imregcorr is asked for the transform between the two pixel
    % grids, moved into world coordinates here: its own composition with
    % these references fails ('Invalid transformation matrix') whenever a
    % pixel is not one world unit wide, which here it never is.
    if isempty(best)
        tc = imregcorr(mv, fx, 'similarity');
        best = affinetform2d(world_from_intrinsic(Rf) * tc.A / ...
            world_from_intrinsic(Rm));
    end

    % a level keeps its refinement only if it scores better than what it
    % started from, so a level that wanders costs nothing
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

end

function warp_to_aligned(rawvol, Raligned, Nslices, Rsample, tforms, out_aligned)
% Both channels warped into aligned space with the recovered transforms, and saved.

fprintf('Warping DAPI and SEP into aligned space... ');
tic;
backvalues = recompute_backvalues(rawvol);
alvol = zeros([Raligned.ImageSize 2 Nslices], 'uint16');
for islice = 1:Nslices
    for ichan = 1:2
        alvol(:, :, ichan, islice) = imwarp(rawvol(:, :, ichan, islice), Rsample, ...
            tforms(islice), 'linear', 'OutputView', Raligned, ...
            'FillValues', double(backvalues(ichan, islice)));
    end
end
saveLargeSliceVolume(alvol, {'DAPI', 'SEP'}, out_aligned);
fprintf('Done! Took %2.2f s\n', toc);

end

function reapply_registration(tp_name, mouse_dir, work_dir, sliceinfo, out_aligned, ...
    out_registered)
% The saved registration re-applied to the aligned SEP volume, in a working folder
% moved into place.

% the elastix folder was saved as an absolute path, drive letter included,
% which may not exist here; the folder itself travels with the mouse
transformparams = load(tp_name);
if ~exist(transformparams.tformbspline_samp20um_to_atlas_20um, 'dir')
    transformparams.tformbspline_samp20um_to_atlas_20um = fullfile(mouse_dir, ...
        'elastix_reverse');
end

% generateRegisteredSliceVolume writes to <procpath>\volume_registered, a fixed
% path, and LightSuite stays as the adults were registered with it, so it
% writes into a working folder and the result is moved into place. sliceinfo is
% otherwise passed as the 'register' mode of run_register_to_atlas passes it,
% with the registered grid in px_atlas: the SEP volume must land on
% volume_registered's grid.
makeNewDir(work_dir);
si = sliceinfo;
si.px_atlas = registered_grid_um();
si.channames = {'DAPI', 'SEP'};
si.procpath = work_dir;
si.slicevolfin = out_aligned;
generateRegisteredSliceVolume(si, transformparams);

if exist(out_registered, 'dir')
    rmdir(out_registered, 's');
end
movefile(fullfile(work_dir, 'volume_registered'), out_registered);
rmdir(work_dir, 's');

end

function [info_new, planes, corr_reg] = check_registered_dapi(registered_dir, ...
    out_registered, mouse_name)
% The registered DAPI of this run against run_register_to_atlas's: the same grid,
% and the correlation of every twentieth plane.

% the registered DAPI of this run against run_register_to_atlas's: same grid,
% then every twentieth plane correlated over the voxels either one covers
fprintf('Checking the registered DAPI against the one run_register_to_atlas wrote... ');
tic;
ref_dapi = fullfile(registered_dir, 'chan01_DAPI.tiff');
new_dapi = fullfile(out_registered, 'chan01_DAPI.tiff');
info_ref = imfinfo(ref_dapi);
info_new = imfinfo(new_dapi);
if numel(info_ref) ~= numel(info_new) || info_ref(1).Height ~= info_new(1).Height ...
        || info_ref(1).Width ~= info_new(1).Width
    error(['run_add_sep_channel: %s registered SEP came out %dx%dx%d against volume_registered''s %dx%dx%d.\n' ...
           'The two are not on the same grid and must not be divided into one another.'], ...
           mouse_name, info_new(1).Height, info_new(1).Width, numel(info_new), ...
           info_ref(1).Height, info_ref(1).Width, numel(info_ref));
end
planes = 1:20:numel(info_ref);
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

end

function plot_sep_diagnostics(Nslices, corr_slice, min_slice_corr, shift_slice, ...
    planes, corr_reg, aligned_dapi, info_aligned, rawvol, Rsample, tforms, Raligned, ...
    info_new, out_registered, mouse_name, diag_dir)
% The diagnostic figure of the recovery and the registration, saved in diag_dir.

fig = figure('Color', 'w', 'Position', [80 80 1500 850], 'Visible', 'off');
tl = tiledlayout(fig, 2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

% the per-slice recovery
nexttile(tl);
plot(1:Nslices, corr_slice, 'o-', 'Color', [0.15 0.15 0.15], ...
    'MarkerFaceColor', [0.15 0.15 0.15], 'MarkerSize', 4);
hold on;
yline(min_slice_corr, '--', 'Color', [0.8 0.2 0.1]);
xlabel('slice (aligned order)');
ylabel('r vs the aligned DAPI, tissue only');
title(sprintf('per-slice transform recovery (median %1.5f, worst shift %.2f px)', ...
    median(corr_slice), max(shift_slice)));
ylim([min(0.9, min(corr_slice) - 0.01) 1.001]);
box off

% the registered DAPI against run_register_to_atlas's, plane by plane
nexttile(tl);
plot(planes - 1, corr_reg, 'o-', 'Color', [0.15 0.15 0.15], ...
    'MarkerFaceColor', [0.15 0.15 0.15], 'MarkerSize', 4);
xlabel('registered plane (ML index)');
ylabel('r vs volume\_registered');
title(sprintf('registered DAPI, this run vs run_register_to_atlas (median %1.4f)', ...
    median(corr_reg, 'omitnan')), 'Interpreter', 'none');
ylim([min(0.9, min(corr_reg) - 0.01) 1.001]);
box off

% the worst slice: its aligned DAPI, the recovered warp of the raw one, and the
% difference
[~, iworst] = min(corr_slice);
fixed = single(imread(aligned_dapi, 'Index', iworst, 'Info', info_aligned));
chk = imwarp(single(rawvol(:, :, 1, iworst)), Rsample, tforms(iworst), 'linear', ...
    'OutputView', Raligned, 'FillValues', 0);
lim = [0 quantile(fixed(fixed > 0), 0.999)];
nexttile(tl);
imagesc(fixed, lim);
axis image off;
colormap(gca, sep_palette('anatomy'));
title(sprintf('aligned DAPI, slice %d (the worst one)', iworst));
nexttile(tl);
imagesc(chk, lim);
axis image off;
colormap(gca, sep_palette('anatomy'));
title(sprintf('recovered warp of the raw DAPI (r = %1.4f)', corr_slice(iworst)));
nexttile(tl);
imagesc(abs(chk - fixed), [0 diff(lim) / 4]);
axis image off;
colormap(gca, sep_palette('intensity'));
title('difference');

% the registered SEP, middle plane
mid = round(numel(info_new) / 2);
sep_mid = single(imread(fullfile(out_registered, 'chan02_SEP.tiff'), 'Index', mid));
nexttile(tl);
imagesc(sep_mid, [0 quantile(sep_mid(sep_mid > 0), 0.999)]);
axis image off;
colormap(gca, sep_palette('intensity'));
title(sprintf('registered SEP, ML plane %d', mid));

title(tl, sprintf('%s -- SEP carried into registered space with the saved transforms', ...
    strrep(mouse_name, '_', ' ')), 'FontWeight', 'bold');
exportgraphics(fig, fullfile(diag_dir, sprintf('%s.png', mouse_name)), ...
    'Resolution', 150);
close(fig);

end

function write_sep_numbers(diag_dir, mouse_name, mouse_type, Nslices, corr_slice, ...
    min_slice_corr, bad, shift_slice, sliceinfo, corr_reg, out_registered)
% The numbers behind the figure, in a text file in diag_dir.

% the numbers behind the figure, for an audit without MATLAB
fid = fopen(fullfile(diag_dir, sprintf('%s.txt', mouse_name)), 'w');
fprintf(fid, 'mouse            %s (%s)\n', mouse_name, mouse_type);
fprintf(fid, 'written          %s\n', datetime('now', 'Format', 'yyyy-MM-dd HH:mm'));
fprintf(fid, 'slices           %d\n', Nslices);
fprintf(fid, 'slice recovery   median r %1.5f, min %1.5f, below %.2f: %s\n', ...
    median(corr_slice), min(corr_slice), min_slice_corr, mat2str(bad));
fprintf(fid, 'residual shift   median %.3f px, max %.3f px (%.2f um at %g um/px)\n', ...
    median(shift_slice), max(shift_slice), ...
    max(shift_slice) * double(sliceinfo.px_process), double(sliceinfo.px_process));
fprintf(fid, 'registered DAPI  median r %1.5f, min %1.5f over %d planes\n', ...
    median(corr_reg, 'omitnan'), min(corr_reg), nnz(~isnan(corr_reg)));
fprintf(fid, 'output           %s\n', out_registered);
fclose(fid);

end

function M = world_from_intrinsic(R)
% The intrinsic-to-world map of an imref2d, by its definition: pixel centres
% sit at 1..N, the first one half a pixel inside the world limit.

M = [R.PixelExtentInWorldX 0 R.XWorldLimits(1) - 0.5 * R.PixelExtentInWorldX; ...
     0 R.PixelExtentInWorldY R.YWorldLimits(1) - 0.5 * R.PixelExtentInWorldY; ...
     0 0 1];
end

function r = image_match(moving, Rmoving, fixed, Rfixed, tform)
% Correlation of the fixed image with the moving one warped onto it, to choose
% between two candidate transforms at one level, where only the order matters.

warped = imwarp(moving, Rmoving, tform, 'linear', 'OutputView', Rfixed, 'FillValues', 0);
m = fixed > 0;
if nnz(m) < 100
    r = NaN;
    return
end
r = corr(double(warped(m)), double(fixed(m)));
end

function [r, shiftpx] = recovery_score(moving, Rmoving, fixed, Rfixed, tform)
% How well the recovered transform reproduces the slice the alignment saved: the
% correlation over tissue, and the residual shift in pixels of the aligned grid.

% two numbers, since they fail differently: the correlation over tissue only
% (the constant background fill would flatter it), and the residual translation,
% which is what a misplaced channel would look like
warped = imwarp(moving, Rmoving, tform, 'linear', 'OutputView', Rfixed, ...
    'FillValues', NaN);
bg = quantile(fixed(fixed > 0), 0.01);
ok = (fixed > 2 * bg) & isfinite(warped);
if nnz(ok) < 100
    r = NaN;
    shiftpx = NaN;
    return
end
r = corr(double(warped(ok)), double(fixed(ok)));

% the residual shift, with the area outside the warp set to zero in both images
a = fixed;
a(~isfinite(warped)) = 0;
b = warped;
b(~isfinite(b)) = 0;
sh = imregcorr(b, a, 'translation');
shiftpx = hypot(sh.Translation(1), sh.Translation(2));
end
