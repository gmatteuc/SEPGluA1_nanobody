%% check_demba_to_allen
% ===== Check that the DeMBA P20 to adult CCF transform reconciles the two =====
%
% Atlas QC, run by hand after the atlas or its crop changes.
%
% The young cohort registers to the age-matched P20 template while the adults
% stay on the adult CCF, so anything voxelwise across the two needs a way to
% carry P20 data into adult space. CCF Translator provides it, built from the
% same deformations Carey used to make DeMBA.
%
% The check is unusually strong: DeMBA's labels are the CCF labels, warped
% onto the P20 template, so warping them back should recover the adult
% annotation, and the overlap with the real CCF annotation should come out
% near 1. Anything much lower means the transform is not doing what it
% claims. Both states are computed here rather than quoted, so the figure
% stands on its own:
%   before  the two atlases aligned only by their AP crops, what a
%           registration would face with no deformation at all, the worst case
%   after   the DeMBA annotation carried into adult space by CCF Translator
% Prints the voxel agreement and the recovery of each region, and saves
% demba_to_allen_transform_qc.png and .fig in young\registration_qc under the
% data root.
%
% Setup: DeMBA P20 against the adult CCF, compared at 20 um. The transformed
% annotation comes from tmp/demba_to_allen.py, which is not in the
% repository; re-run it if the atlas or its crop changes. Run sep_setup_paths
% first, once per MATLAB session.

clear; clc; close all;

%% Settings

% project folders, worked out from where the code sits, so the tree can be moved
% or copied to another drive as is (SEP_DATA_ROOT points the data elsewhere)
paths = get_paths();

% the transformed annotation, already in adult space at 20 um
transformed_file = fullfile(paths.data, 'atlas_demba_p20', ...
                            'annotation_in_allen_space_20um.nii.gz');

% common working resolution, in um
work_res_um = 20;

% regions shown in the matrices: major divisions, then the visual and
% somatosensory areas with their short labels
regions_of_interest = {'Isocortex', 'ventricular systems', 'Hippocampal formation', ...
                       'Striatum', 'Thalamus', 'Cerebellum', 'fiber tracts', ...
                       'Hypothalamus', 'Midbrain', 'Olfactory areas'};

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

% smallest region, in voxels, with a stable number
min_voxels = 200;

% output folder, and whether to save the figure
out_dir = fullfile(paths.data, 'young', 'registration_qc');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end
save_figure = true;

% colours: before the transform grey, after it orange
before_color = [0.55 0.55 0.55];
after_color  = [0.95 0.55 0.10];

%% Load

atlas_adult = get_atlas('ccf');
atlas_young = get_atlas('demba_p20');
csv_dir = atlas_adult.dir;

if ~exist(transformed_file, 'file')
    error(['Transformed annotation not found:\n  %s\n' ...
           'Run tmp/demba_to_allen.py first (it needs the venv at ' ...
           'code\\tools\\venv_atlas).'], transformed_file);
end

fprintf('loading volumes...\n');
av_adult_full = niftiread(fullfile(atlas_adult.dir, atlas_adult.annotation_file));
av_young_full = niftiread(fullfile(atlas_young.dir, atlas_young.annotation_file));
av_trans      = niftiread(transformed_file);

% the adult atlas at the working resolution, subsampled rather than averaged:
% these are labels, and the mean of two region ids is not a region
step = work_res_um / atlas_adult.res_um;
av_adult_20 = av_adult_full(1:step:end, 1:step:end, 1:step:end);

fprintf('  adult %s at %g um -> %s at %g um\n', mat2str(size(av_adult_full)), ...
    atlas_adult.res_um, mat2str(size(av_adult_20)), work_res_um);
fprintf('  transformed DeMBA %s\n', mat2str(size(av_trans)));

if ~isequal(size(av_trans), size(av_adult_20))
    error(['Transformed volume is %s but the adult atlas at %g um is %s.\n' ...
           'They must share a grid for this comparison to mean anything.'], ...
           mat2str(size(av_trans)), work_res_um, mat2str(size(av_adult_20)));
end

%% After: both cropped to the adult analysis range

crop_after = round(atlas_adult.default_aplims / step);
adult_after = av_adult_20(crop_after(1):crop_after(2), :, :);
trans_after = av_trans(crop_after(1):crop_after(2), :, :);
fprintf('\nadult crop [%d %d] at 10 um -> [%d %d] at %g um, %s\n', ...
    atlas_adult.default_aplims, crop_after, work_res_um, mat2str(size(adult_after)));

%% Before: the two atlases aligned only by their crops

% each atlas cropped to its own AP span, then resampled to a common shape: the
% comparison with no deformation between them, the worst case, not what the
% pipeline incurs
crop_y = atlas_young.default_aplims;
young_cropped = av_young_full(crop_y(1):crop_y(2), :, :);
before_shape  = size(adult_after);
trans_before  = imresize3(young_cropped, before_shape, 'Method', 'nearest');
adult_before  = adult_after;
fprintf('before-state: DeMBA crop [%d %d] resampled %s -> %s\n', ...
    crop_y, mat2str(size(young_cropped)), mat2str(before_shape));

%% Agreement in both states

[agree_before, rec_before] = agreement(adult_before, trans_before, min_voxels);
[agree_after,  rec_after ] = agreement(adult_after,  trans_after,  min_voxels);

fprintf('\n--- voxel agreement where both label something ---\n');
fprintf('  before (crop alignment only) : %.4f\n', agree_before);
fprintf('  after  (CCF Translator)      : %.4f\n', agree_after);

fprintf('\n--- per-region recovery ---\n');
fprintf('  %-8s %8s %8s %10s %10s\n', 'state', 'median', 'n', '>0.5', '>0.8');
fprintf('  %-8s %8.4f %8d %9.0f%% %9.0f%%\n', 'before', median(rec_before), ...
    numel(rec_before), 100*mean(rec_before > 0.5), 100*mean(rec_before > 0.8));
fprintf('  %-8s %8.4f %8d %9.0f%% %9.0f%%\n', 'after', median(rec_after), ...
    numel(rec_after), 100*mean(rec_after > 0.5), 100*mean(rec_after > 0.8));

%% Overlap matrices

fprintf('\ncomputing overlap matrices...\n');
mat_div_before = region_overlap_matrix(csv_dir, adult_before, trans_before, ...
    regions_of_interest);
mat_div_after  = region_overlap_matrix(csv_dir, adult_after,  trans_after,  ...
    regions_of_interest);
mat_ctx_before = region_overlap_matrix(csv_dir, adult_before, trans_before, ...
    cortical_areas);
mat_ctx_after  = region_overlap_matrix(csv_dir, adult_after,  trans_after,  ...
    cortical_areas);

fprintf('\n--- diagonal, major divisions ---\n');
fprintf('  %-24s %8s %8s\n', 'region', 'before', 'after');
for i = 1:numel(regions_of_interest)
    fprintf('  %-24s %8.3f %8.3f\n', regions_of_interest{i}, ...
        mat_div_before(i, i), mat_div_after(i, i));
end

fprintf('\n--- diagonal, cortical areas ---\n');
fprintf('  %-14s %8s %8s\n', 'area', 'before', 'after');
for i = 1:numel(cortical_areas)
    fprintf('  %-14s %8.3f %8.3f\n', cortical_labels{i}, ...
        mat_ctx_before(i, i), mat_ctx_after(i, i));
end

%% Figure

fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'pixels', ...
             'Position', [50 50 1680 940]);
tl = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% the major divisions, before and after
nexttile(tl);
draw_overlap_matrix(mat_div_before, regions_of_interest, regions_of_interest, 7.5, 7, ...
    'called this in the transformed DeMBA');
title({'Major divisions BEFORE', sprintf('crop alignment only, diagonal median %.2f', ...
       median(diag(mat_div_before)))}, 'FontSize', 10)

nexttile(tl);
draw_overlap_matrix(mat_div_after, regions_of_interest, regions_of_interest, 7.5, 7, ...
    'called this in the transformed DeMBA');
title({'Major divisions AFTER', sprintf('CCF Translator, diagonal median %.2f', ...
       median(diag(mat_div_after)))}, 'FontSize', 10)

% the recovery of every region, before and after
nexttile(tl);
edges = linspace(0, 1, 40);
histogram(rec_before, edges, 'Normalization', 'probability', ...
          'FaceColor', before_color, 'EdgeColor', 'none', 'FaceAlpha', 0.75);
hold on
histogram(rec_after, edges, 'Normalization', 'probability', ...
          'FaceColor', after_color, 'EdgeColor', 'none', 'FaceAlpha', 0.75);
xlabel('fraction of the adult region recovered')
ylabel('fraction of regions')
legend({sprintf('before  (median %.2f)', median(rec_before)), ...
        sprintf('after   (median %.2f)', median(rec_after))}, ...
        'Box', 'off', 'Location', 'northwest')
title({sprintf('All %d regions', numel(rec_after)), ...
       sprintf('%.0f%% above 0.8 after, %.0f%% before', ...
       100*mean(rec_after > 0.8), 100*mean(rec_before > 0.8))}, 'FontSize', 10)
box off

% the visual and somatosensory areas, before and after
nexttile(tl);
draw_overlap_matrix(mat_ctx_before, cortical_labels, cortical_labels, 7, 6, ...
    'called this in the transformed DeMBA');
title({'Visual + somatosensory BEFORE', sprintf('diagonal median %.2f', ...
       median(diag(mat_ctx_before)))}, 'FontSize', 10)

nexttile(tl);
draw_overlap_matrix(mat_ctx_after, cortical_labels, cortical_labels, 7, 6, ...
    'called this in the transformed DeMBA');
title({'Visual + somatosensory AFTER', sprintf('diagonal median %.2f', ...
       median(diag(mat_ctx_after)))}, 'FontSize', 10)

% their diagonals as bars
nexttile(tl);
d_before = diag(mat_ctx_before);
d_after  = diag(mat_ctx_after);
b = barh([d_before d_after]);
b(1).FaceColor = before_color;
b(1).EdgeColor = 'none';
b(2).FaceColor = after_color;
b(2).EdgeColor = 'none';
set(gca, 'YTick', 1:numel(cortical_labels), 'YTickLabel', cortical_labels, ...
         'TickLabelInterpreter', 'none', 'FontSize', 8, 'YDir', 'reverse');
xlim([0 1]);
xlabel('fraction recovered')
legend({'before', 'after'}, 'Box', 'off', 'Location', 'southeast')
title('The areas this project reports on', 'FontSize', 10)
box off

title(tl, sprintf(['DeMBA P20 carried into adult CCF space by CCF Translator   |   ' ...
                   'voxel agreement %.3f'], agree_after), ...
      'Interpreter', 'none', 'FontWeight', 'bold');

if save_figure
    png_name = fullfile(out_dir, 'demba_to_allen_transform_qc.png');
    fig_name = fullfile(out_dir, 'demba_to_allen_transform_qc.fig');
    exportgraphics(fig, png_name, 'Resolution', 150);
    savefig(fig, fig_name);
    fprintf('\nsaved:\n  %s\n  %s\n', png_name, fig_name);
end

% ===== Local functions =====

function [agree, recovery] = agreement(ref, test, min_voxels)
% Voxel agreement where both label something, and the fraction of each region
% of ref (min_voxels or more) that test labels the same.

both = ref > 0 & test > 0;
agree = nnz(both & ref == test) / nnz(both);

labels = intersect(unique(ref(:)), unique(test(:)));
labels = labels(labels ~= 0);
recovery = [];
for k = 1:numel(labels)
    m = ref == labels(k);
    n = nnz(m);
    if n < min_voxels
        continue
    end
    recovery(end+1, 1) = nnz(m & test == labels(k)) / n; %#ok<AGROW>
end

end
