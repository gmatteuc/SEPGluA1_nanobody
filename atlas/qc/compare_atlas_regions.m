%% compare_atlas_regions
% ===== Measure how differently regions sit in the P20 atlas and the adult one =====
%
% Atlas QC, the measurement behind registering the young brains to an
% age-matched atlas, run by hand.
%
% Gross brain size is not the criterion: the analysis reports per-region
% means, so what matters is whether region boundaries sit differently enough
% at P20 to bias those means. For every region present in both atlases it
% measures:
%   1. the volume fraction, region voxels over brain voxels: a fraction rather
%      than a volume, so the 11% AP scale difference between the two templates
%      cancels and only differences in proportion remain
%   2. the centroid in normalised brain coordinates, again scale-free, reported
%      as a displacement in adult-equivalent mm
%   3. the overlap of the labels, once both atlases are resampled to one
%      grid from their crops, which were measured to span the same anatomy:
%      that applies the uniform AP stretch and the 2x in-plane difference at
%      once, taking the crop correspondence and a linear map as given (both
%      checked on their own: area profile, and the 678-region regression,
%      r = 0.998). Volume fractions say whether a region has the right size,
%      this whether it is in the right place. Between named regions it is the
%      second figure, kept apart on purpose: a linear alignment alone is the
%      worst case, not what the pipeline incurs (every brain is registered
%      nonlinearly, and the two label sets are one ontology one transform
%      apart), so it ranks the sensitive regions but misleads as a headline
%   4. the thickness of the cortical ribbon by AP level
%   5. the contrast of the two templates: the Allen template averages 1675
%      brains, the DeMBA P21 anchor 12, and P20 is interpolated from it. The
%      local gradient is taken on a lightly smoothed image, over the
%      template's own mean brain brightness, with both templates at one
%      resolution: noise raises the raw gradient, so the noisier template
%      would score as the more contrasted one, which is backwards, and a finer
%      grid would score lower gradients per voxel for free. What the smoothing
%      removes is reported apart, as noise.
% If regional proportions agree to a few percent, the age-matched atlas buys
% little and the simpler single-atlas option is defensible; if cortex or
% ventricles are systematically off, it earns its keep. Prints every number,
% and saves atlas_region_comparison.png and .fig, and
% atlas_label_overlap_worstcase.png, in young\registration_qc under the data
% root.
%
% Setup: the adult CCF against DeMBA P20. Run sep_setup_paths first, once
% per MATLAB session.

clear; clc; close all;

%% Settings

% project folders, worked out from where the code sits, so the tree can be moved
% or copied to another drive as is (SEP_DATA_ROOT points the data elsewhere)
paths = get_paths();

% the two atlases, as get_atlas keys
atlas_key_adult = 'ccf';
atlas_key_young = 'demba_p20';

% smallest adult region, in voxels, with a stable centroid and fraction
min_voxels_adult = 5000;

% regions called out by name in the report and on the figure
regions_of_interest = {'Isocortex', 'ventricular systems', 'Hippocampal formation', ...
                       'Striatum', 'Thalamus', 'Cerebellum', 'fiber tracts', ...
                       'Hypothalamus', 'Midbrain', 'Olfactory areas'};

% the cortical areas this project is about, so the overlap is also shown at the
% level the analysis reports; siblings only, never a parent with its own child,
% so the columns do not double-count; short labels keep the matrix readable
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

% evenly spaced coronal levels through the crop at which the thickness of the
% isocortex ribbon is measured
n_thickness_levels = 40;

% resolution, in um, both atlases are resampled to for the label overlap: the
% registration's own, which keeps the volumes to a few hundred MB
overlap_res_um = 20;

% coronal levels sampled for the template contrast, and the smoothing before
% measuring it (sigma in voxels at overlap_res_um): the smoothing separates
% anatomical contrast from the noise of a template built from few brains
n_contrast_levels = 60;
contrast_smooth_sigma = 2;

% where to cut the side-by-side contrast panel (fraction of the crop from the
% anterior end), and its shared colour scale, in multiples of each template's
% own mean brain brightness
contrast_patch_level = 0.50;
contrast_patch_clim  = [0 2.2];

% what the shaded band on the centroid panel means, from Carey 2025 for DeMBA
% itself: the transformations that built the atlas land 0.103-0.186 mm from
% expert consensus, and an average human expert 0.145-0.220 mm from it, so a
% displacement inside the band is as small as the atlas's own build accuracy
accuracy_note = 'DeMBA''s own build accuracy, 0.10-0.19 mm (expert rater: 0.14-0.22)';

% output folder, and whether to save the figures
out_dir = fullfile(paths.data, 'young', 'registration_qc');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end
save_figure = true;

% colours: the adult atlas grey, the young atlas orange, the called-out regions
% darker orange, the identity line grey, regions that agree within 10% light
% grey
adult_color = [0.35 0.35 0.35];
young_color = [0.95 0.55 0.10];
roi_color   = [0.85 0.50 0.00];
unity_color = [0.45 0.45 0.45];
quiet_color = [0.78 0.78 0.78];

%% Load and crop both atlases

atlas_adult = get_atlas(atlas_key_adult);
atlas_young = get_atlas(atlas_key_young);

% the CSVs that name the regions live with the adult atlas, and both annotation
% volumes are in parcellation_index space, so one set serves both
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

%% Volume fraction and normalised centroid of each region

% the labels in both atlases; few shared means two different id spaces
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

% ratio of proportions; log2, so over- and under-representation are symmetric
ratio = fy ./ fa;
log2_ratio = log2(ratio);

% centroid displacement in normalised brain units, converted to millimetres of
% the adult brain so the number means something physical
adult_span_mm = brain_span_mm(brain_adult, atlas_adult.res_um);
disp_norm = cy - ca;
disp_mm = disp_norm .* adult_span_mm;
disp_total_mm = sqrt(sum(disp_mm.^2, 2));

fprintf('\n--- Regional proportions (volume fraction, young / adult) ---\n');
fprintf('  median ratio            %.3f\n', median(ratio));
fprintf('  within +/- 5%%           %.0f%% of regions\n', ...
    100 * mean(abs(log2_ratio) < log2(1.05)));
fprintf('  within +/- 10%%          %.0f%% of regions\n', ...
    100 * mean(abs(log2_ratio) < log2(1.10)));
fprintf('  within +/- 25%%          %.0f%% of regions\n', ...
    100 * mean(abs(log2_ratio) < log2(1.25)));
fprintf('  median |log2 ratio|     %.3f  (= %.1f%% typical difference)\n', ...
    median(abs(log2_ratio)), 100 * (2^median(abs(log2_ratio)) - 1));

fprintf('\n--- Centroid displacement (adult-equivalent mm) ---\n');
fprintf('  median %.3f mm, 90th pct %.3f mm, max %.3f mm\n', ...
    median(disp_total_mm), prctile(disp_total_mm, 90), max(disp_total_mm));
fprintf('  by axis, median |shift|: AP %.3f  DV %.3f  ML %.3f mm\n', ...
    median(abs(disp_mm(:,1))), median(abs(disp_mm(:,2))), median(abs(disp_mm(:,3))));
fprintf(['  for scale, DeMBA''s own transformations validate at 0.103-0.186 mm mean\n' ...
         '  landmark error, so a shift of this size is near the atlas noise floor.\n']);

% biology or label-transfer error? A developmental difference should not
% depend on a region's size; transfer noise hits the small regions far more
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

    % the region's fraction of each brain
    mask_a = get_allen_region_mask(csv_dir, av_adult, regions_of_interest(i), brain_adult);
    mask_y = get_allen_region_mask(csv_dir, av_young, regions_of_interest(i), brain_young);

    f_a = nnz(mask_a) / nnz(brain_adult);
    f_y = nnz(mask_y) / nnz(brain_young);
    roi_frac_a(i) = f_a;
    roi_frac_y(i) = f_y;

    % and how far its centroid moves, in adult mm
    c_a = normalised_centroid(mask_a, brain_adult);
    c_y = normalised_centroid(mask_y, brain_young);
    shift = norm((c_y - c_a) .* adult_span_mm);

    fprintf('%-26s %9.3f%% %9.3f%% %8.3f %10.3f\n', ...
        regions_of_interest{i}, 100 * f_a, 100 * f_y, f_y / f_a, shift);

end

%% Label overlap on a common grid

% both atlases on one grid from their crops (see the header): for each adult
% region, what the young atlas calls the same voxels
grid_common = round(size(av_adult) * atlas_adult.res_um / overlap_res_um);
fprintf('\n--- Label overlap on a common %g um grid %s ---\n', ...
    overlap_res_um, mat2str(grid_common));

av_a_rs = imresize3(av_adult, grid_common, 'Method', 'nearest');
av_y_rs = imresize3(av_young, grid_common, 'Method', 'nearest');

% leaf-level self-agreement, in one pass: of the voxels the adult atlas assigns
% to a region, the fraction the young atlas assigns to the same region
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

% the same between named regions, small enough to read as a matrix: the major
% divisions, then the cortical areas the analysis reports on
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

[thk_adult, lev_adult] = ribbon_thickness(ctx_adult, atlas_adult.res_um, ...
    n_thickness_levels);
[thk_young, lev_young] = ribbon_thickness(ctx_young, atlas_young.res_um, ...
    n_thickness_levels);

valid = ~isnan(thk_adult) & ~isnan(thk_young);
fprintf('  adult mean %.3f mm, young mean %.3f mm, ratio %.3f\n', ...
    mean(thk_adult(valid)), mean(thk_young(valid)), ...
    mean(thk_young(valid)) / mean(thk_adult(valid)));
fprintf('  measured as twice the 95th percentile of the distance transform\n');
fprintf('  inside the ribbon -- a proxy, but applied identically to both.\n');

%% Template contrast

% the local gradient of the smoothed template over its mean brain brightness,
% both templates at one resolution (see the header)
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

% two rows: the three measurements across the top, the pair of coronal images
% along the bottom, where the wide row suits their aspect ratio
fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'pixels', ...
             'Position', [50 50 1620 900]);
tl = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% panel 1: how far regions move between the two atlases, against the accuracy
% the atlas itself was built to
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

% panel 2: cortical thickness along AP
nexttile(tl);
plot(lev_adult, thk_adult, 'Color', adult_color, 'LineWidth', 2);
hold on
plot(lev_young, thk_young, 'Color', young_color, 'LineWidth', 2);
xlabel('position along the crop (0 = anterior)')
ylabel('cortical ribbon thickness (mm)')
legend({'adult CCF', 'DeMBA P20'}, 'Box', 'off', 'Location', 'south')
title(sprintf('Cortical thickness matches (%.2f vs %.2f mm, ratio %.2f)', ...
      mean(thk_adult(valid)), mean(thk_young(valid)), ...
      mean(thk_young(valid)) / mean(thk_adult(valid))), 'FontSize', 10)
box off

% panel 3: the contrast difference as a distribution, the gradient taken on a
% lightly smoothed image (see template_contrast)
nexttile(tl);
edges = linspace(0, prctile([grad_adult; grad_young], 99), 60);
histogram(grad_adult, edges, 'Normalization', 'probability', ...
          'FaceColor', adult_color, 'EdgeColor', 'none', 'FaceAlpha', 0.65);
hold on
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

% bottom row: both templates at the same level, each over its own mean brain
% brightness so exposure cancels, on one colour scale; what is left is how much
% the image varies from place to place, which is what contrast means here
nexttile(tl, [1 3]);
patch_pair = contrast_patch(tv_adult, brain_adult, atlas_adult.res_um, ...
                            tv_young, brain_young, atlas_young.res_um, ...
                            overlap_res_um, contrast_patch_level);
imagesc(patch_pair, contrast_patch_clim);
colormap(gca, sep_palette('anatomy'))
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

%% Second figure: the overlap matrices

% the worst case of a linear alignment alone, kept apart from the first figure
% (see the header)
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

% ===== Local functions =====

function [frac, cent, count] = region_stats(av, labels, brain)
% Volume fraction of the brain, normalised centroid (AP, DV, ML) and voxel
% count of every label.

n = numel(labels);
frac  = zeros(n, 1);
cent  = zeros(n, 3);
count = zeros(n, 1);

n_brain = nnz(brain);
[bb_min, bb_size] = brain_box(brain);

% one pass per AP plane, accumulating against a compact label index so nothing
% the size of the largest label value is ever allocated
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
    dvv = dv_grid(hit);
    dvv = dvv(good);
    mlv = ml_grid(hit);
    mlv = mlv(good);

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

function [bb_min, bb_size] = brain_box(brain)
% Bounding box of the brain (first plane and size along AP, DV, ML), so
% centroids are scale-free.

ap = find(squeeze(any(any(brain, 2), 3)));
dv = find(squeeze(any(any(brain, 1), 3)));
ml = find(squeeze(any(any(brain, 1), 2)));
bb_min  = [ap(1) dv(1) ml(1)];
bb_size = [ap(end) - ap(1) + 1, dv(end) - dv(1) + 1, ml(end) - ml(1) + 1];

end

function span_mm = brain_span_mm(brain, res_um)
% Size of the brain's bounding box in mm (AP, DV, ML).

[~, bb_size] = brain_box(brain);
span_mm = bb_size * res_um / 1000;

end

function c = normalised_centroid(mask, brain)
% Centroid of one mask in the brain's bounding box, 0 to 1 along each axis.

[bb_min, bb_size] = brain_box(brain);
idx = find(mask);
[i1, i2, i3] = ind2sub(size(mask), idx);
c = ([mean(i1) mean(i2) mean(i3)] - bb_min) ./ bb_size;

end

function [thk, levels] = ribbon_thickness(mask, res_um, n_levels)
% Thickness of a ribbon-shaped mask in mm at n_levels coronal levels (NaN where
% the plane holds fewer than 50 voxels), and the levels from 0 to 1 along AP.

n_ap = size(mask, 1);
idx = round(linspace(1, n_ap, n_levels));
thk = nan(1, n_levels);

for k = 1:n_levels
    plane = squeeze(mask(idx(k), :, :));
    if nnz(plane) < 50
        continue
    end

    % the distance to the nearest voxel outside peaks at half the thickness on
    % the ribbon's midline, so twice a high percentile of it estimates the thickness
    d = bwdist(~plane);
    thk(k) = 2 * prctile(d(plane), 95) * res_um / 1000;
end

levels = (idx - 1) / (n_ap - 1);

end

function m = region_overlap_matrix(csv_dir, av_a, av_y, region_names)
% Fraction of each named region of av_a (rows) that av_y calls each named
% region (columns).

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

function name = biggest_leak(row, self_idx, labels)
% The other region a row of the overlap matrix leaks into most, with the
% fraction, or 'nothing else' below 0.005.

row(self_idx) = -inf;
[v, j] = max(row);
if v <= 0.005
    name = 'nothing else';
else
    name = sprintf('%s %.2f', labels{j}, v);
end

end

function draw_overlap_matrix(m, row_labels, col_labels, tick_font, cell_font)
% Draw one overlap matrix in hot, with each value worth reading printed in its
% cell.

imagesc(m, [0 1]);
colormap(gca, sep_palette('intensity'));
axis square
n = size(m, 1);
set(gca, 'XTick', 1:n, 'XTickLabel', col_labels, ...
         'YTick', 1:n, 'YTickLabel', row_labels, ...
         'TickLabelInterpreter', 'none', 'FontSize', tick_font);
xtickangle(45)

% each cell worth reading, in a colour that shows on hot
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

function [grad_norm, global_cv, noise_level] = template_contrast(tv, brain, ...
                                        res_um, target_res_um, n_levels, smooth_sigma)
% Local gradient of the lightly smoothed template over its mean brain brightness,
% its global CV, and its high-frequency noise (the method is in the header).

scale = res_um / target_res_um;

% global CV, std/mean over the whole brain: bit depth, exposure and histogram
% matching cancel, but a smooth front-to-back ramp scores as fine detail does
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

    % erode, so the brain's edge does not dominate the gradient
    msk = imerode(msk, strel('disk', 3));
    if nnz(msk) < 200
        continue
    end

    % Sobel gradient magnitude of the lightly smoothed image, over the mean
    % brain brightness: blind to the global range, sensitive to sharp edges
    smoothed = imgaussfilt(img, smooth_sigma);
    g = imgradient(smoothed, 'sobel');
    grad_norm  = [grad_norm;  double(g(msk)) / double(mean_int)];       %#ok<AGROW>
    noise_vals = [noise_vals; double(img(msk) - smoothed(msk))];        %#ok<AGROW>
end

noise_level = std(noise_vals) / double(mean_int);

end

function pair = contrast_patch(tv_a, brain_a, res_a, tv_b, brain_b, res_b, ...
                               target_res_um, level_frac)
% The coronal plane at level_frac of each template, side by side on one
% canvas, each over its own mean brain brightness.

img_a = one_patch(tv_a, brain_a, res_a, target_res_um, level_frac);
img_b = one_patch(tv_b, brain_b, res_b, target_res_um, level_frac);

h = max(size(img_a, 1), size(img_b, 1));
w = max(size(img_a, 2), size(img_b, 2));
pair = [pad_centre(img_a, h, w), pad_centre(img_b, h, w)];

end

function img = one_patch(tv, brain, res_um, target_res_um, level_frac)
% One template's coronal plane at level_frac, at target_res_um, over its mean
% brain brightness, cropped to the brain.

i = max(1, min(size(tv, 1), round(1 + level_frac * (size(tv, 1) - 1))));
img = single(squeeze(tv(i, :, :)));
msk = squeeze(brain(i, :, :));

scale = res_um / target_res_um;
if scale ~= 1
    img = imresize(img, scale, 'bilinear');
    msk = imresize(msk, scale, 'nearest');
end

% over this template's own mean brain brightness, so the two are compared on
% variation rather than on exposure
img = img / mean(single(tv(brain)));
img(~msk) = 0;

rows = find(any(msk, 2));
cols = find(any(msk, 1));
img = img(rows(1):rows(end), cols(1):cols(end));

end

function out = pad_centre(img, h, w)
% Centre img on a zero canvas of h x w.

out = zeros(h, w, 'like', img);
[ih, iw] = size(img);
r0 = floor((h - ih) / 2) + 1;
c0 = floor((w - iw) / 2) + 1;
out(r0:r0 + ih - 1, c0:c0 + iw - 1) = img;

end
