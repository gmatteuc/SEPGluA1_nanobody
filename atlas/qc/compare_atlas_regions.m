clear all
close all
clc

% /// QC helper: do regions sit differently in the P20 atlas than the adult one? ///
% This is the measurement that decides whether registering the young brains to
% an age-matched atlas is worth its costs. Gross brain size is not the criterion
% -- the analysis reports per-region means, so what matters is whether region
% BOUNDARIES sit differently enough at P20 to bias those means.
%
% For every region present in both atlases it computes:
%   (1) volume fraction, region voxels over brain voxels. A fraction rather than
%       a volume, so the 11% AP scale difference between the two templates
%       cancels and only genuine differences in proportion survive
%   (2) centroid in normalised brain coordinates, again scale-free, reported as
%       a displacement in adult-equivalent mm
%   (3) cortical ribbon thickness by AP level, and ventricular volume fraction,
%       the two places a P20-versus-adult difference should be largest
%
% How to read the result: if regional proportions agree to a few percent, the
% age-matched atlas buys little and the simpler single-atlas option is
% defensible. If cortex or ventricles are systematically off, it earns its keep.

%% User-defined parameters

% Where the project lives. Derived from the location of the code rather than
% written out, so the tree can be moved or copied to another drive as is.
paths = get_paths();

% The two atlases, as get_atlas keys
atlas_key_adult = 'ccf';
atlas_key_young = 'demba_p20';

% Ignore regions too small for a stable centroid or fraction
min_voxels_adult = 5000;

% Regions called out by name in the report and on the figure
regions_of_interest = {'Isocortex', 'ventricular systems', 'Hippocampal formation', ...
                       'Striatum', 'Thalamus', 'Cerebellum', 'fiber tracts', ...
                       'Hypothalamus', 'Midbrain', 'Olfactory areas'};

% The cortical areas this project is actually about, so the overlap is also
% shown at the level the analysis reports rather than only for whole divisions.
% These are sibling areas, never a parent together with its own child, so the
% columns do not double-count. Short labels keep the matrix readable.
cortical_areas = {'Primary visual area', ...
                  'Anterolateral visual area', ...
                  'Rostrolateral visual area', ...
                  'Anteromedial visual area', ...
                  'posteromedial visual area', ...
                  'Lateral visual area', ...
                  'Postrhinal area', ...
                  'Primary somatosensory area, barrel field', ...
                  'Primary somatosensory area, lower limb', ...
                  'Primary somatosensory area, upper limb', ...
                  'Primary somatosensory area, mouth', ...
                  'Supplemental somatosensory area'};
cortical_labels = {'VISp', 'VISal (AL)', 'VISrl (RL)', 'VISam (AM)', 'VISpm (PM)', ...
                   'VISl (LM)', 'VISpor (POR)', 'SSp-bfd', 'SSp-ll', 'SSp-ul', ...
                   'SSp-m', 'SSs'};

% Cortical thickness is measured on the Isocortex ribbon, at this many evenly
% spaced coronal levels through the crop
n_thickness_levels = 40;

% Resolution the two atlases are resampled onto to compute label overlap. 20 um
% is the resolution the registration itself runs at, and keeps the volumes to a
% few hundred MB.
overlap_res_um = 20;

% How many coronal levels the template-contrast comparison samples, and how
% much the image is smoothed before measuring structure. The smoothing is what
% separates anatomical contrast from the noise of a template built from few
% brains -- sigma is in voxels at overlap_res_um.
n_contrast_levels = 60;
contrast_smooth_sigma = 2;

% Where to cut the side-by-side contrast panel, and the shared colour scale it
% is shown on, in multiples of each template's own mean brain brightness.
contrast_patch_level = 0.50;
contrast_patch_clim  = [0 2.2];

% What the shaded band on the centroid panel means, spelled out rather than
% left as jargon. These are the numbers Carey 2025 report for DeMBA itself:
% the transformations that built the atlas land 0.103-0.186 mm from expert
% consensus, and an average human expert lands 0.145-0.220 mm from it. A
% displacement inside that band is indistinguishable from how accurately the
% atlas was built in the first place.
accuracy_note = 'DeMBA''s own build accuracy, 0.10-0.19 mm (expert rater: 0.14-0.22)';

% Output
out_dir = fullfile(paths.data, 'young', 'registration_qc');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
save_figure = true;

% Color palette
adult_color = [0.35 0.35 0.35];   % gray, adult atlas
young_color = [0.95 0.55 0.10];   % strong orange, young atlas
roi_color   = [0.85 0.50 0.00];   % darker orange, called-out regions
unity_color = [0.45 0.45 0.45];   % gray, identity line
quiet_color = [0.78 0.78 0.78];   % light gray, regions that agree within 10%

%% Load and crop both atlases

atlas_adult = get_atlas(atlas_key_adult);
atlas_young = get_atlas(atlas_key_young);

% The CSVs that name the regions live with the adult atlas, and both annotation
% volumes are in the same parcellation_index space, so one set serves both.
csv_dir = atlas_adult.dir;

av_adult = niftiread(fullfile(atlas_adult.dir, atlas_adult.annotation_file));
av_young = niftiread(fullfile(atlas_young.dir, atlas_young.annotation_file));

lim_a = atlas_adult.default_aplims;
lim_y = atlas_young.default_aplims;
av_adult = av_adult(lim_a(1):lim_a(2), :, :);
av_young = av_young(lim_y(1):lim_y(2), :, :);

brain_adult = av_adult > 0;
brain_young = av_young > 0;

fprintf('adult: %s at %g um, %d brain voxels\n', ...
    mat2str(size(av_adult)), atlas_adult.res_um, nnz(brain_adult));
fprintf('young: %s at %g um, %d brain voxels\n', ...
    mat2str(size(av_young)), atlas_young.res_um, nnz(brain_young));

%% Per-region volume fraction and normalised centroid

labels = intersect(unique(av_adult(:)), unique(av_young(:)));
labels = labels(labels ~= 0);
fprintf('\nregions present in both atlases: %d\n', numel(labels));

if numel(labels) < 100
    error(['Only %d shared labels. The two annotation volumes are probably in ' ...
           'different ID spaces -- see registration_qc\\ATLAS_PARAMETERS.md.'], ...
           numel(labels));
end

[frac_adult, cent_adult, count_adult] = region_stats(av_adult, labels, brain_adult);
[frac_young, cent_young, count_young] = region_stats(av_young, labels, brain_young);

keep = count_adult >= min_voxels_adult & count_young > 0;
fprintf('regions big enough to compare: %d\n', nnz(keep));

labels_k = labels(keep);
fa = frac_adult(keep);
fy = frac_young(keep);
ca = cent_adult(keep, :);
cy = cent_young(keep, :);

% Ratio of proportions. Log2 so over- and under-representation are symmetric.
ratio = fy ./ fa;
log2_ratio = log2(ratio);

% Centroid displacement, in normalised brain units, converted to millimetres of
% the adult brain so the number means something physical.
adult_span_mm = brain_span_mm(brain_adult, atlas_adult.res_um);
disp_norm = cy - ca;
disp_mm = disp_norm .* adult_span_mm;
disp_total_mm = sqrt(sum(disp_mm.^2, 2));

fprintf('\n--- Regional proportions (volume fraction, young / adult) ---\n');
fprintf('  median ratio            %.3f\n', median(ratio));
fprintf('  within +/- 5%%           %.0f%% of regions\n', 100 * mean(abs(log2_ratio) < log2(1.05)));
fprintf('  within +/- 10%%          %.0f%% of regions\n', 100 * mean(abs(log2_ratio) < log2(1.10)));
fprintf('  within +/- 25%%          %.0f%% of regions\n', 100 * mean(abs(log2_ratio) < log2(1.25)));
fprintf('  median |log2 ratio|     %.3f  (= %.1f%% typical difference)\n', ...
    median(abs(log2_ratio)), 100 * (2^median(abs(log2_ratio)) - 1));

fprintf('\n--- Centroid displacement (adult-equivalent mm) ---\n');
fprintf('  median %.3f mm, 90th pct %.3f mm, max %.3f mm\n', ...
    median(disp_total_mm), prctile(disp_total_mm, 90), max(disp_total_mm));
fprintf('  by axis, median |shift|: AP %.3f  DV %.3f  ML %.3f mm\n', ...
    median(abs(disp_mm(:,1))), median(abs(disp_mm(:,2))), median(abs(disp_mm(:,3))));
fprintf(['  for scale, DeMBA''s own transformations validate at 0.103-0.186 mm mean\n' ...
         '  landmark error, so a shift of this size is near the atlas noise floor.\n']);

% Is the disagreement biology, or just label-transfer error on small regions?
% If it were real developmental difference it should not care how big a region
% is; if it is transfer noise, small regions suffer far more.
size_adult = count_adult(keep);
rho = corr(log10(double(size_adult)), abs(log2_ratio), 'Type', 'Spearman');
fprintf('\n--- Is the difference size-dependent? ---\n');
fprintf('  Spearman |log2 ratio| vs region size: rho = %+.3f\n', rho);

[~, order] = sort(size_adult, 'descend');
for n_top = [50 100 200]
    sel = order(1:min(n_top, numel(order)));
    fprintf(['  largest %3d regions: median difference %4.1f%%, ' ...
             '%3.0f%% within 10%%, median shift %.3f mm\n'], ...
        n_top, 100 * (2^median(abs(log2_ratio(sel))) - 1), ...
        100 * mean(abs(log2_ratio(sel)) < log2(1.10)), ...
        median(disp_total_mm(sel)));
end

%% The named regions

fprintf('\n--- Regions of interest ---\n');
fprintf('%-26s %10s %10s %8s %10s\n', 'region', 'adult %', 'young %', 'ratio', 'shift mm');
fprintf('%s\n', repmat('-', 1, 68));

roi_frac_a = nan(1, numel(regions_of_interest));
roi_frac_y = nan(1, numel(regions_of_interest));

for i = 1:numel(regions_of_interest)

    mask_a = get_allen_region_mask(csv_dir, av_adult, regions_of_interest(i), brain_adult);
    mask_y = get_allen_region_mask(csv_dir, av_young, regions_of_interest(i), brain_young);

    f_a = nnz(mask_a) / nnz(brain_adult);
    f_y = nnz(mask_y) / nnz(brain_young);
    roi_frac_a(i) = f_a;
    roi_frac_y(i) = f_y;

    c_a = normalised_centroid(mask_a, brain_adult);
    c_y = normalised_centroid(mask_y, brain_young);
    shift = norm((c_y - c_a) .* adult_span_mm);

    fprintf('%-26s %9.3f%% %9.3f%% %8.3f %10.3f\n', ...
        regions_of_interest{i}, 100 * f_a, 100 * f_y, f_y / f_a, shift);

end

%% Label overlap, once both atlases are on a common grid

% Volume fractions say whether a region is the right SIZE. They cannot say
% whether it is in the right PLACE. This does: it puts both atlases on one grid
% and asks, for each adult region, what the young atlas calls those same voxels.
%
% The common grid comes from the two crops, which were measured to span the same
% anatomy. Resampling both to one shape therefore applies the uniform AP stretch
% and the 2x in-plane difference at once. That is an assumption -- it takes the
% crop correspondence and a linear map between them as given -- but both were
% verified independently (area profile and 678-region regression, r = 0.998).

grid_common = round(size(av_adult) * atlas_adult.res_um / overlap_res_um);
fprintf('\n--- Label overlap on a common %g um grid %s ---\n', ...
    overlap_res_um, mat2str(grid_common));

av_a_rs = imresize3(av_adult, grid_common, 'Method', 'nearest');
av_y_rs = imresize3(av_young, grid_common, 'Method', 'nearest');

% Leaf-level self-agreement, in one pass: of the voxels the adult atlas assigns
% to a region, what fraction does the young atlas assign to the same region?
both_labelled = av_a_rs > 0 & av_y_rs > 0;
same_label = both_labelled & (av_a_rs == av_y_rs);

lut = zeros(double(max(labels)) + 1, 1);
lut(double(labels) + 1) = 1:numel(labels);
idx_all  = lut(double(av_a_rs(both_labelled)) + 1);
idx_same = lut(double(av_a_rs(same_label)) + 1);
tot_vox  = accumarray(idx_all(idx_all > 0),  1, [numel(labels) 1]);
hit_vox  = accumarray(idx_same(idx_same > 0), 1, [numel(labels) 1]);

self_cov = hit_vox ./ tot_vox;
big = tot_vox >= 200;
fprintf('  leaf regions scored: %d\n', nnz(big));
fprintf('  self-overlap: median %.3f, 25th pct %.3f, 75th pct %.3f\n', ...
    median(self_cov(big)), prctile(self_cov(big), 25), prctile(self_cov(big), 75));
fprintf('  regions above 0.5 self-overlap: %.0f%%\n', 100 * mean(self_cov(big) > 0.5));

% And the same thing between named regions, small enough to read as a matrix.
% Done twice: once for the major divisions, once for the cortical areas the
% analysis actually reports on.
overlap_mat = region_overlap_matrix(csv_dir, av_a_rs, av_y_rs, regions_of_interest);
n_roi_ov = numel(regions_of_interest);

fprintf('\n  major divisions, diagonal (adult region found again in the young atlas):\n');
for i = 1:n_roi_ov
    fprintf('    %-24s %.3f\n', regions_of_interest{i}, overlap_mat(i, i));
end

cortex_mat = region_overlap_matrix(csv_dir, av_a_rs, av_y_rs, cortical_areas);
n_ctx_ov = numel(cortical_areas);

fprintf('\n  cortical areas, diagonal:\n');
for i = 1:n_ctx_ov
    fprintf('    %-14s %.3f   (leaks most to %s)\n', cortical_labels{i}, ...
        cortex_mat(i, i), biggest_leak(cortex_mat(i, :), i, cortical_labels));
end
fprintf('  cortical diagonal: median %.3f, range %.3f-%.3f\n', ...
    median(diag(cortex_mat)), min(diag(cortex_mat)), max(diag(cortex_mat)));

%% Cortical ribbon thickness by AP level

fprintf('\n--- Cortical ribbon thickness ---\n');

ctx_adult = get_allen_region_mask(csv_dir, av_adult, {'Isocortex'}, brain_adult);
ctx_young = get_allen_region_mask(csv_dir, av_young, {'Isocortex'}, brain_young);

[thk_adult, lev_adult] = ribbon_thickness(ctx_adult, atlas_adult.res_um, n_thickness_levels);
[thk_young, lev_young] = ribbon_thickness(ctx_young, atlas_young.res_um, n_thickness_levels);

valid = ~isnan(thk_adult) & ~isnan(thk_young);
fprintf('  adult mean %.3f mm, young mean %.3f mm, ratio %.3f\n', ...
    mean(thk_adult(valid)), mean(thk_young(valid)), ...
    mean(thk_young(valid)) / mean(thk_adult(valid)));
fprintf('  measured as twice the 95th percentile of the distance transform\n');
fprintf('  inside the ribbon -- a proxy, but applied identically to both.\n');

%% Template contrast

% The Allen template averages 1675 brains, the DeMBA P21 anchor 12, and P20 is
% interpolated from it. That difference is visible by eye; this measures it.
% Local gradient magnitude divided by the template's own mean brain intensity,
% so the comparison is about structure and not about brightness scaling. Both
% templates are brought to the same resolution first -- otherwise a finer grid
% would score lower gradients per voxel for free.

tv_adult = niftiread(fullfile(atlas_adult.dir, atlas_adult.template_file));
tv_young = niftiread(fullfile(atlas_young.dir, atlas_young.template_file));
tv_adult = tv_adult(lim_a(1):lim_a(2), :, :);
tv_young = tv_young(lim_y(1):lim_y(2), :, :);

[grad_adult, cv_adult, noise_adult] = template_contrast(tv_adult, brain_adult, ...
    atlas_adult.res_um, overlap_res_um, n_contrast_levels, contrast_smooth_sigma);
[grad_young, cv_young, noise_young] = template_contrast(tv_young, brain_young, ...
    atlas_young.res_um, overlap_res_um, n_contrast_levels, contrast_smooth_sigma);

ia_crop_mm = size(av_adult, 1) * atlas_adult.res_um / 1000;
iy_crop_mm = size(av_young, 1) * atlas_young.res_um / 1000;
info_len_mm = ia_crop_mm;

fprintf('\n--- Template contrast ---\n');
fprintf('  global CV (std/mean inside the brain)    : adult %.3f, young %.3f, ratio %.2f\n', ...
    cv_adult, cv_young, cv_adult / cv_young);
fprintf('  local smoothed gradient (median / mean)  : adult %.4f, young %.4f, ratio %.2f\n', ...
    median(grad_adult), median(grad_young), median(grad_adult) / median(grad_young));
fprintf('  high-frequency noise / mean intensity   : adult %.4f, young %.4f, ratio %.2f\n', ...
    noise_adult, noise_young, noise_young / noise_adult);
fprintf('  brains averaged: adult 1675, DeMBA P21 anchor 12 (P20 interpolated)\n');

%% Figure

% Two rows: the three measurements across the top, and the pair of coronal
% images along the bottom, where the wide row suits their aspect ratio.
fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'pixels', ...
             'Position', [50 50 1620 900]);
tl = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% Panel 1: how far regions move between the two atlases, against the accuracy
% the atlas itself was built to.
nexttile(tl);
h_hist = histogram(disp_total_mm, 50, 'FaceColor', young_color, 'EdgeColor', 'none');
hold on
yl = ylim;
h_band = fill([0.103 0.186 0.186 0.103], [0 0 yl(2) yl(2)], unity_color, ...
              'FaceAlpha', 0.20, 'EdgeColor', 'none');
xline(median(disp_total_mm), 'Color', adult_color, 'LineWidth', 1.8);
ylim(yl)
xlabel('region centroid shift between the atlases (mm)')
ylabel('regions')
legend([h_hist h_band], {'regions', 'DeMBA intrinsic error'}, ...
       'Box', 'off', 'Location', 'northeast')
title(sprintf('Regions move %.2f mm = %.1f%% of the %.0f mm brain', ...
      median(disp_total_mm), 100 * median(disp_total_mm) / info_len_mm, info_len_mm), ...
      'FontSize', 10)
box off

% Panel 2: cortical thickness along AP
nexttile(tl);
plot(lev_adult, thk_adult, 'Color', adult_color, 'LineWidth', 2); hold on
plot(lev_young, thk_young, 'Color', young_color, 'LineWidth', 2);
xlabel('position along the crop (0 = anterior)')
ylabel('cortical ribbon thickness (mm)')
legend({'adult CCF', 'DeMBA P20'}, 'Box', 'off', 'Location', 'south')
title(sprintf('Cortical thickness matches (%.2f vs %.2f mm, ratio %.2f)', ...
      mean(thk_adult(valid)), mean(thk_young(valid)), ...
      mean(thk_young(valid)) / mean(thk_adult(valid))), 'FontSize', 10)
box off

% Panel 3: the contrast difference as a distribution. Gradient is measured on a
% lightly smoothed image, because the raw gradient is inflated by noise and the
% young template averages 12 brains against the adult's 1675 -- unsmoothed, it
% reads as the MORE contrasted of the two, which is backwards.
nexttile(tl);
edges = linspace(0, prctile([grad_adult; grad_young], 99), 60);
histogram(grad_adult, edges, 'Normalization', 'probability', ...
          'FaceColor', adult_color, 'EdgeColor', 'none', 'FaceAlpha', 0.65); hold on
histogram(grad_young, edges, 'Normalization', 'probability', ...
          'FaceColor', young_color, 'EdgeColor', 'none', 'FaceAlpha', 0.65);
xlabel('local smoothed gradient (fraction of mean brightness)')
ylabel('fraction of voxels')
legend({sprintf('adult CCF  (1675 brains, global CV %.2f)', cv_adult), ...
        sprintf('DeMBA P20  (12 brains, global CV %.2f)', cv_young)}, ...
        'Box', 'off', 'Location', 'northeast')
title({sprintf('Adult has %.1fx the local smoothed gradient', ...
               median(grad_adult) / median(grad_young)), ...
       '(how sharply brightness changes point to point)'}, 'FontSize', 10)
box off

% Bottom row: what that difference actually looks like. Both templates at the
% same matched level, each divided by its own mean brain brightness so overall
% exposure cancels, then shown on one shared colour scale. What is left is how
% much the image varies from place to place -- which is what contrast means here.
nexttile(tl, [1 3]);
patch_pair = contrast_patch(tv_adult, brain_adult, atlas_adult.res_um, ...
                            tv_young, brain_young, atlas_young.res_um, ...
                            overlap_res_um, contrast_patch_level);
imagesc(patch_pair, contrast_patch_clim);
colormap(gca, gray)
axis image off
hold on
xline(size(patch_pair, 2) / 2, 'Color', 'w', 'LineWidth', 1.5);
text(0.25, 0.05, sprintf('adult CCF (1675 brains)  global CV %.2f', cv_adult), ...
     'Units', 'normalized', 'Color', 'w', 'FontSize', 9, ...
     'HorizontalAlignment', 'center');
text(0.75, 0.05, sprintf('DeMBA P20 (12 brains)  global CV %.2f', cv_young), ...
     'Units', 'normalized', 'Color', 'w', 'FontSize', 9, ...
     'HorizontalAlignment', 'center');
cb = colorbar;
cb.Label.String = 'brightness / mean brain brightness';
title({sprintf('Adult has %.1fx the global CV', cv_adult / cv_young), ...
       '(how much light and dark differ across the brain)'}, 'FontSize', 10)

title(tl, 'Adult CCF vs DeMBA P20', 'Interpreter', 'none', 'FontWeight', 'bold');

if save_figure
    png_name = fullfile(out_dir, 'atlas_region_comparison.png');
    fig_name = fullfile(out_dir, 'atlas_region_comparison.fig');
    exportgraphics(fig, png_name, 'Resolution', 150);
    savefig(fig, fig_name);
    fprintf('\nsaved:\n  %s\n  %s\n', png_name, fig_name);
end

%% Second figure: the overlap matrices, kept for our own reference

% Deliberately separate. These measure how far apart the two atlases are with
% only a linear alignment between them, which is the worst case for the
% pipeline rather than what it would actually incur -- every brain is
% registered nonlinearly, and the two label sets are the same ontology one
% transform apart, so a perfect registration makes them agree exactly. Useful
% for ranking which regions are sensitive, misleading as a headline.

fig2 = figure('Visible', 'off', 'Color', 'w', 'Units', 'pixels', ...
              'Position', [50 50 1250 560]);
tl2 = tiledlayout(fig2, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile(tl2);
draw_overlap_matrix(overlap_mat, regions_of_interest, regions_of_interest, 7.5, 7);
title(sprintf('Major divisions: diagonal median %.2f', median(diag(overlap_mat))), ...
      'FontSize', 10)

nexttile(tl2);
draw_overlap_matrix(cortex_mat, cortical_labels, cortical_labels, 7, 6);
title(sprintf('Visual and somatosensory areas: diagonal median %.2f', ...
      median(diag(cortex_mat))), 'FontSize', 10)

title(tl2, ['Label overlap, linear alignment only -- WORST CASE, not the ' ...
            'pipeline''s error'], 'Interpreter', 'none', 'FontWeight', 'bold');

if save_figure
    png2 = fullfile(out_dir, 'atlas_label_overlap_worstcase.png');
    exportgraphics(fig2, png2, 'Resolution', 150);
    fprintf('  %s\n', png2);
end

%% Local function: volume fraction and normalised centroid for every label

function [frac, cent, count] = region_stats(av, labels, brain)

    n = numel(labels);
    frac  = zeros(n, 1);
    cent  = zeros(n, 3);
    count = zeros(n, 1);

    n_brain = nnz(brain);
    [bb_min, bb_size] = brain_box(brain);

    % One pass per AP plane, accumulating against a compact label index so
    % nothing the size of the largest label value is ever allocated.
    lut = zeros(double(max(labels)) + 1, 1);
    lut(double(labels) + 1) = 1:n;

    sum_ap = zeros(n, 1);
    sum_dv = zeros(n, 1);
    sum_ml = zeros(n, 1);

    [n_ap, n_dv, n_ml] = size(av);
    [dv_grid, ml_grid] = ndgrid(1:n_dv, 1:n_ml);

    for i = 1:n_ap
        plane = double(squeeze(av(i, :, :)));
        hit = plane > 0;
        if ~any(hit(:))
            continue
        end
        vals = plane(hit);
        idx = lut(vals + 1);
        good = idx > 0;
        if ~any(good)
            continue
        end
        idx = idx(good);
        dvv = dv_grid(hit); dvv = dvv(good);
        mlv = ml_grid(hit); mlv = mlv(good);

        count  = count  + accumarray(idx, 1,   [n 1]);
        sum_ap = sum_ap + accumarray(idx, i,   [n 1]);
        sum_dv = sum_dv + accumarray(idx, dvv, [n 1]);
        sum_ml = sum_ml + accumarray(idx, mlv, [n 1]);
    end

    frac = count / n_brain;

    cent(:, 1) = (sum_ap ./ count - bb_min(1)) / bb_size(1);
    cent(:, 2) = (sum_dv ./ count - bb_min(2)) / bb_size(2);
    cent(:, 3) = (sum_ml ./ count - bb_min(3)) / bb_size(3);

end

%% Local function: bounding box of the brain, so centroids are scale-free

function [bb_min, bb_size] = brain_box(brain)
    ap = find(squeeze(any(any(brain, 2), 3)));
    dv = find(squeeze(any(any(brain, 1), 3)));
    ml = find(squeeze(any(any(brain, 1), 2)));
    bb_min  = [ap(1) dv(1) ml(1)];
    bb_size = [ap(end) - ap(1) + 1, dv(end) - dv(1) + 1, ml(end) - ml(1) + 1];
end

%% Local function: physical size of the brain bounding box

function span_mm = brain_span_mm(brain, res_um)
    [~, bb_size] = brain_box(brain);
    span_mm = bb_size * res_um / 1000;
end

%% Local function: normalised centroid of one mask

function c = normalised_centroid(mask, brain)
    [bb_min, bb_size] = brain_box(brain);
    idx = find(mask);
    [i1, i2, i3] = ind2sub(size(mask), idx);
    c = ([mean(i1) mean(i2) mean(i3)] - bb_min) ./ bb_size;
end

%% Local function: ribbon thickness of a mask, by AP level

function [thk, levels] = ribbon_thickness(mask, res_um, n_levels)

    n_ap = size(mask, 1);
    idx = round(linspace(1, n_ap, n_levels));
    thk = nan(1, n_levels);

    for k = 1:n_levels
        plane = squeeze(mask(idx(k), :, :));
        if nnz(plane) < 50
            continue
        end
        % Distance to the nearest non-cortex voxel. In a ribbon this peaks at
        % half the thickness along the ribbon's midline, so twice a high
        % percentile of it estimates the thickness.
        d = bwdist(~plane);
        thk(k) = 2 * prctile(d(plane), 95) * res_um / 1000;
    end

    levels = (idx - 1) / (n_ap - 1);

end

%% Local function: overlap matrix between two labellings, for a set of regions

function m = region_overlap_matrix(csv_dir, av_a, av_y, region_names)
    n = numel(region_names);
    mask_a = cell(1, n);
    mask_y = cell(1, n);
    for i = 1:n
        mask_a{i} = get_allen_region_mask(csv_dir, av_a, region_names(i), av_a > 0);
        mask_y{i} = get_allen_region_mask(csv_dir, av_y, region_names(i), av_y > 0);
    end
    m = zeros(n);
    for i = 1:n
        denom = nnz(mask_a{i});
        if denom == 0
            continue
        end
        for j = 1:n
            m(i, j) = nnz(mask_a{i} & mask_y{j}) / denom;
        end
    end
end

%% Local function: which other region does a row leak into most?

function name = biggest_leak(row, self_idx, labels)
    row(self_idx) = -inf;
    [v, j] = max(row);
    if v <= 0.005
        name = 'nothing else';
    else
        name = sprintf('%s %.2f', labels{j}, v);
    end
end

%% Local function: draw one overlap matrix

function draw_overlap_matrix(m, row_labels, col_labels, tick_font, cell_font)
    imagesc(m, [0 1]);
    colormap(gca, hot);
    axis square
    n = size(m, 1);
    set(gca, 'XTick', 1:n, 'XTickLabel', col_labels, ...
             'YTick', 1:n, 'YTickLabel', row_labels, ...
             'TickLabelInterpreter', 'none', 'FontSize', tick_font);
    xtickangle(45)

    % Print each cell worth reading, in a colour that survives the hot colormap
    for i = 1:n
        for j = 1:n
            v = m(i, j);
            if v < 0.02
                continue
            end
            if v > 0.55
                txt_col = [0 0 0];
            else
                txt_col = [1 1 1];
            end
            text(j, i, sprintf('%.2f', v), 'HorizontalAlignment', 'center', ...
                 'FontSize', cell_font, 'Color', txt_col);
        end
    end

    cb = colorbar;
    cb.Label.String = 'fraction of the adult region';
    xlabel('called this in DeMBA P20')
    ylabel('adult CCF region')
end

%% Local function: structural contrast and noise of a template

function [grad_norm, global_cv, noise_level] = template_contrast(tv, brain, ...
                                        res_um, target_res_um, n_levels, smooth_sigma)
    % Raw local gradient is the wrong measure here: the young template averages
    % 12 brains against the adult's 1675, so it is noisier, and noise RAISES
    % local gradient. Measuring it directly says the young template has more
    % contrast, which is backwards. So the gradient is taken on a lightly
    % smoothed image (structure), and the part removed by that smoothing is
    % reported separately (noise).
    scale = res_um / target_res_um;

    % Global CV: coefficient of variation of brightness over the whole brain,
    % std/mean. Dimensionless, so bit depth, exposure and the histogram matching
    % the templates were given all cancel. Blind to spatial scale -- a smooth
    % front-to-back brightness ramp would score the same as fine detail.
    inside = single(tv(brain));
    mean_int = mean(inside);
    global_cv = std(inside) / mean_int;

    idx = round(linspace(1, size(tv, 1), n_levels));
    grad_norm  = [];
    noise_vals = [];

    for k = 1:n_levels
        img = single(squeeze(tv(idx(k), :, :)));
        msk = squeeze(brain(idx(k), :, :));
        if nnz(msk) < 500
            continue
        end
        if scale ~= 1
            img = imresize(img, scale, 'bilinear');
            msk = imresize(msk, scale, 'nearest');
        end
        % Erode so the brain edge itself does not dominate the gradient
        msk = imerode(msk, strel('disk', 3));
        if nnz(msk) < 200
            continue
        end
        % Local smoothed gradient: Sobel gradient magnitude of the lightly
        % smoothed image, divided by mean brain brightness. Blind to global
        % range, sensitive to how sharp the edges are. Unsmoothed it inverts,
        % because noise dominates and the 12-brain average carries more of it.
        smoothed = imgaussfilt(img, smooth_sigma);
        g = imgradient(smoothed, 'sobel');
        grad_norm  = [grad_norm;  double(g(msk)) / double(mean_int)];       %#ok<AGROW>
        noise_vals = [noise_vals; double(img(msk) - smoothed(msk))];        %#ok<AGROW>
    end

    noise_level = std(noise_vals) / double(mean_int);
end

%% Local function: matched coronal patch from each template, common scaling

function pair = contrast_patch(tv_a, brain_a, res_a, tv_b, brain_b, res_b, ...
                               target_res_um, level_frac)
    img_a = one_patch(tv_a, brain_a, res_a, target_res_um, level_frac);
    img_b = one_patch(tv_b, brain_b, res_b, target_res_um, level_frac);

    h = max(size(img_a, 1), size(img_b, 1));
    w = max(size(img_a, 2), size(img_b, 2));
    pair = [pad_centre(img_a, h, w), pad_centre(img_b, h, w)];
end

function img = one_patch(tv, brain, res_um, target_res_um, level_frac)
    i = max(1, min(size(tv, 1), round(1 + level_frac * (size(tv, 1) - 1))));
    img = single(squeeze(tv(i, :, :)));
    msk = squeeze(brain(i, :, :));

    scale = res_um / target_res_um;
    if scale ~= 1
        img = imresize(img, scale, 'bilinear');
        msk = imresize(msk, scale, 'nearest');
    end

    % Divide by this template's own mean brain brightness, so the two are
    % compared on variation rather than on exposure.
    img = img / mean(single(tv(brain)));
    img(~msk) = 0;

    rows = find(any(msk, 2));
    cols = find(any(msk, 1));
    img = img(rows(1):rows(end), cols(1):cols(end));
end

function out = pad_centre(img, h, w)
    out = zeros(h, w, 'like', img);
    [ih, iw] = size(img);
    r0 = floor((h - ih) / 2) + 1;
    c0 = floor((w - iw) / 2) + 1;
    out(r0:r0 + ih - 1, c0:c0 + iw - 1) = img;
end
