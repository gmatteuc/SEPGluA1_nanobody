clear all
close all
clc

% /// Pipeline script #9: batch comparison of LR-sum (nano or auto) vs Allen ISH ///
% Reads gene_targets.csv, loops over all genes, and for each:
%   - Downloads/caches the ISH grid from the Allen API
%   - Repairs bad ISH sections at native 200um resolution
%   - Upsamples to 10um, computes L+R, z-scores
%   - Computes 3 correlation metrics (region Spearman, region Pearson, voxel Pearson)
%   - Generates diagnostic video, scatter, paired bar charts
%   - Saves gene_result_<symbol>.mat (light) + ish_lr_sum_<symbol>.mat (heavy)
% After the loop, produces cross-gene summary figures.
% Designed to run unattended overnight on ~50 genes.

%% User-defined parameters

% Where the project lives. Derived from the location of the code rather than
% written out, so the tree can be moved or copied to another drive as is.
paths = get_paths();

% Channel to compare against ISH: 'nano' (default, surface GluA1) or 'auto'
% (autofluorescence control). Must match the channel used in the P8 cache.
channel = 'nano';

% Which P8 cache to load (must exist — run P8 first). The channel is appended
% so this points to the matching nano or auto cohort cache.
p8_merged_tag    = ['merged_naive_rws_' channel];
p8_smooth_suffix = '_nosmooth';

% Gene target list (CSV with columns: symbol, experiment_id, category, plane, description)
gene_list_path = fullfile(paths.data, 'gene_targets.csv');

% ISH grid orientation (same for all Allen coronal experiments)
ish_raw_permute = [1, 2, 3];
ish_raw_flipdim = [false, false, false];
ish_variant     = 'energy';

% Output toggles
generate_comparison_video = true;
generate_scatter_plot     = true;
generate_paired_barchart  = true;
save_diagnostic_slices    = true;

% Z-scoring
zscore_method = 'standard';

% Video labels (on diff panel only)
show_slice_labels       = true;
label_fontsize          = 7;
label_color             = [0.85 0.85 0.85];
label_smooth_halfwindow = 5;
label_min_px_per_slice  = 150;

% Region analysis
macro_divi_list = {'Isocortex','OLF','HPF','CTXsp','STR','PAL','TH','HY','MB'};
roi_erode_radius    = 3;
dist_weight_power   = 4;  % exponent for distance weighting (must match P8)
analysis_slice_range = [100, 700];

% Video clim for difference panel
video_clim_diff = [-3, 3];

% Paths
base_root     = paths.data;
p8_out_dir    = fullfile(base_root, 'comparisons', p8_merged_tag);
ish_cache_dir = fullfile(base_root, 'atlas_ish');
if ~exist(ish_cache_dir, 'dir'), mkdir(ish_cache_dir); end
summary_dir = fullfile(base_root, 'comparisons', ...
    [p8_merged_tag '_vs_ish_summary' p8_smooth_suffix]);
if ~exist(summary_dir, 'dir'), mkdir(summary_dir); end

%% Read gene target list

gene_tbl = readtable(gene_list_path, 'TextType', 'string', 'Delimiter', ',');
n_genes = height(gene_tbl);
fprintf('Loaded %d gene targets from: %s\n', n_genes, gene_list_path);

%% Load Allen atlas (once)

allenDir = fullfile(base_root, 'atlas');
addpath(allenDir);
fprintf('Loading Allen atlas...\n');
AllenVol = niftiread(fullfile(allenDir, 'annotation_10.nii.gz'));
allen_full_size = size(AllenVol);
limits = [180 1079];
AllenCrop  = AllenVol(limits(1):limits(2), :, :);
brainMask  = AllenCrop > 0;
half_atlas = AllenCrop(:, :, 1:end);
clear AllenVol

%% Load P8 nano cache (once)

cache_path = fullfile(p8_out_dir, ...
    ['mean_lr_sum_' p8_merged_tag p8_smooth_suffix '.mat']);
if ~exist(cache_path, 'file')
    error('P8 cache not found: %s\nRun P8 first.', cache_path);
end
fprintf('Loading P8 cache...\n');
S = load(cache_path);
mean_lr_sum      = S.mean_lr_sum;
brainMask_merged = S.brainMask_merged;
half_width       = S.half_width;

%% Build analysis brain mask (slice-range restricted, once)

brainMask_analysis = brainMask_merged;
brainMask_analysis(1:analysis_slice_range(1)-1, :, :) = false;
brainMask_analysis(analysis_slice_range(2)+1:end, :, :) = false;

%% Build shared ROI masks (once — reused across all genes)

fprintf('Building STRU-leaf ROI list from ontology...\n');
[roi_list, roi_acronyms, roi_macro] = ...
    build_stru_leaves_by_macro_local(allenDir, macro_divi_list);
n_rois = numel(roi_list);
fprintf('  %d STRU-leaf ROIs.\n', n_rois);

macro_groups = struct();
for r = 1:n_rois
    mkey = matlab.lang.makeValidName(roi_macro{r});
    if ~isfield(macro_groups, mkey), macro_groups.(mkey) = []; end
    macro_groups.(mkey)(end+1) = r; %#ok<SAGROW>
end

% Build ROI masks + distance weights for region aggregation (cache to disk)
% Two methods computed in parallel:
%   (A) eroded masks: binary erosion by roi_erode_radius voxels
%   (B) distance-weighted: continuous weights from bwdist, center=1, border=0
% Both are cached and both produce per-region means for every gene, so we
% can compare them systematically after the batch run.

[~, ~, n_width] = size(AllenCrop);
half_w = floor(n_width / 2);
masks_cache_path = fullfile(p8_out_dir, ...
    sprintf('roi_masks_r%d_pw%d_slices%d-%d.mat', ...
    roi_erode_radius, dist_weight_power, analysis_slice_range(1), analysis_slice_range(2)));

if exist(masks_cache_path, 'file')
    fprintf('  Loading cached ROI masks: %s\n', masks_cache_path);
    M = load(masks_cache_path);
    roi_indices     = M.roi_indices;      % cell of uint32 linear index vectors
    roi_weights_dw  = M.roi_weights_dw;   % cell of single weight vectors (dist^power)
    roi_indices_ero = M.roi_indices_ero;   % cell of uint32 linear index vectors (eroded)
    roi_px_raw      = M.roi_px_raw;
    roi_px_ero      = M.roi_px_ero;
else
    fprintf('  Building ROI masks + distance weights for %d regions...\n', n_rois);
    if roi_erode_radius > 0
        se_erode = strel('sphere', roi_erode_radius);
    end
    roi_indices     = cell(n_rois, 1);
    roi_weights_dw  = cell(n_rois, 1);
    roi_indices_ero = cell(n_rois, 1);
    roi_px_raw      = zeros(n_rois, 1);
    roi_px_ero      = zeros(n_rois, 1);
    for r = 1:n_rois
        try
            m_raw = get_allen_region_mask(allenDir, AllenCrop, roi_list(r), brainMask, '');
            if size(m_raw, 3) > half_w, m_raw = m_raw(:, :, 1:half_w); end
            idx_raw = uint32(find(m_raw));
            roi_indices{r} = idx_raw;
            roi_px_raw(r) = numel(idx_raw);
            % Erosion indices
            if roi_erode_radius > 0
                m_ero = imerode(m_raw, se_erode);
                if nnz(m_ero) == 0, m_ero = m_raw; end
            else
                m_ero = m_raw;
            end
            roi_indices_ero{r} = uint32(find(m_ero));
            roi_px_ero(r) = numel(roi_indices_ero{r});
            % Distance weights (stored only at ROI voxel positions)
            dist = single(bwdist(~m_raw));
            mx = max(dist(m_raw));
            if mx > 0
                w = (dist(m_raw) / mx) .^ dist_weight_power;
            else
                w = ones(numel(idx_raw), 1, 'single');
            end
            roi_weights_dw{r} = w;
        catch
            roi_indices{r}     = uint32([]);
            roi_weights_dw{r}  = single([]);
            roi_indices_ero{r} = uint32([]);
            warning('Could not build mask for: %s', roi_list{r});
        end
        if mod(r, 50) == 0, fprintf('    %d/%d masks done\n', r, n_rois); end
    end
    save(masks_cache_path, 'roi_indices', 'roi_weights_dw', 'roi_indices_ero', ...
        'roi_px_raw', 'roi_px_ero', ...
        'roi_list', 'roi_acronyms', 'roi_macro', '-v7.3');
    fprintf('  ROI masks cached: %s\n', masks_cache_path);
end

% Precompute nano-side region means with BOTH methods (once)
data_vol_masked = mean_lr_sum;
data_vol_masked(~brainMask_merged) = NaN;
data_vol_masked(1:analysis_slice_range(1)-1, :, :) = NaN;
data_vol_masked(analysis_slice_range(2)+1:end, :, :) = NaN;

data_z_global = zscore_in_mask(mean_lr_sum, brainMask_analysis);
data_z_global(~brainMask_merged) = NaN;
data_z_masked = data_z_global;
data_z_masked(1:analysis_slice_range(1)-1, :, :) = NaN;
data_z_masked(analysis_slice_range(2)+1:end, :, :) = NaN;

% Method A: eroded means (using sparse indices)
roi_data_mean_ero = nan(n_rois, 1);
roi_data_z_ero    = nan(n_rois, 1);
for r = 1:n_rois
    if roi_px_ero(r) == 0, continue; end
    vn = data_vol_masked(roi_indices_ero{r}); valid = ~isnan(vn);
    if any(valid), roi_data_mean_ero(r) = mean(vn(valid)); end
    vnz = data_z_masked(roi_indices_ero{r}); valid = ~isnan(vnz);
    if any(valid), roi_data_z_ero(r) = mean(vnz(valid)); end
end

% Method B: distance-weighted means (using sparse indices + weights)
roi_data_mean_dw = nan(n_rois, 1);
roi_data_z_dw    = nan(n_rois, 1);
for r = 1:n_rois
    if roi_px_raw(r) == 0, continue; end
    w = double(roi_weights_dw{r});
    vn = data_vol_masked(roi_indices{r}); valid = ~isnan(vn);
    ws = sum(w(valid));
    if ws > 0, roi_data_mean_dw(r) = sum(vn(valid) .* w(valid)) / ws; end
    vnz = data_z_masked(roi_indices{r}); valid = ~isnan(vnz);
    ws = sum(w(valid));
    if ws > 0, roi_data_z_dw(r) = sum(vnz(valid) .* w(valid)) / ws; end
end

%% Precompute label centroids (once — for video overlay)

label_centroids = [];
label_acronyms  = {};
if generate_comparison_video && show_slice_labels
    centroid_cache = fullfile(p8_out_dir, ...
        ['label_centroids_' p8_merged_tag p8_smooth_suffix '.mat']);
    if exist(centroid_cache, 'file')
        fprintf('Loading cached label centroids...\n');
        LC = load(centroid_cache);
        label_centroids = LC.label_centroids;
        label_acronyms  = LC.label_acronyms;
    else
        fprintf('Computing label centroids for %d regions...\n', n_rois);
        atlas_half_w = floor(size(AllenCrop, 3) / 2);
        n_sl = size(AllenCrop, 1);
        label_acronyms  = roi_acronyms;
        label_centroids = nan(n_sl, n_rois, 2);
        for r = 1:n_rois
            try
                m_temp = get_allen_region_mask(allenDir, AllenCrop, roi_list(r), brainMask, '');
                if size(m_temp, 3) > atlas_half_w, m_temp = m_temp(:, :, 1:atlas_half_w); end
                for z = 1:n_sl
                    sm = squeeze(m_temp(z, :, :));
                    [ri, ci] = find(sm);
                    if numel(ri) >= label_min_px_per_slice
                        label_centroids(z, r, 1) = mean(ri);
                        label_centroids(z, r, 2) = mean(ci);
                    end
                end
            catch, end
        end
        win = 2 * label_smooth_halfwindow + 1;
        for r = 1:n_rois
            for d = 1:2
                v = label_centroids(:, r, d);
                pr = ~isnan(v);
                s = movmean(v, win, 'omitnan');
                s(~pr) = NaN;
                label_centroids(:, r, d) = s;
            end
        end
        save(centroid_cache, 'label_centroids', 'label_acronyms', '-v7.3');
    end
end

% Nano voxel vector for voxel-level Pearson (precompute once)
data_vox = mean_lr_sum(brainMask_analysis);

%% Custom macro color map (shared across scatter plots)

macro_color_map = containers.Map( ...
    {'CTXsp','HPF','HY','Isocortex','MB','OLF','PAL','STR','TH'}, ...
    {[0.20 0.40 0.80], [0.85 0.33 0.10], [0.93 0.69 0.13], ...
     [0.55 0.25 0.70], [0.80 0.15 0.15], [0.30 0.75 0.93], ...
     [0.64 0.08 0.18], [0.10 0.30 0.70], [0.90 0.50 0.10]});

%% ====================== GENE LOOP ======================================

gene_results = struct();  % collector for summary

for gi = 1:n_genes
    gene_sym  = char(gene_tbl.symbol(gi));
    gene_id   = gene_tbl.experiment_id(gi);
    gene_cat  = char(gene_tbl.category(gi));
    gene_desc = char(gene_tbl.description(gi));

    fprintf('\n========== [%d/%d] %s (exp %d, %s) ==========\n', ...
        gi, n_genes, gene_sym, gene_id, gene_cat);

    out_dir = fullfile(base_root, 'comparisons', ...
        [p8_merged_tag '_vs_ish_' lower(gene_sym) p8_smooth_suffix]);
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end

    result_path  = fullfile(out_dir, ['gene_result_' gene_sym '.mat']);
    ish_vol_path = fullfile(out_dir, ['ish_lr_sum_' gene_sym '.mat']);

    % --- Check if already fully processed ---
    if exist(result_path, 'file') && exist(ish_vol_path, 'file')
        fprintf('  Already processed (gene_result + ish_lr_sum exist). Skipping.\n');
        R = load(result_path);
        gene_results(gi).symbol       = gene_sym;
        gene_results(gi).category     = gene_cat;
        gene_results(gi).r_spearman   = R.r_spearman;
        gene_results(gi).r_pearson    = R.r_pearson;
        gene_results(gi).r_voxel      = R.r_voxel_pearson;
        gene_results(gi).n_bad        = R.n_bad_sections;
        if isfield(R, 'r_spearman_ero')
            gene_results(gi).r_spearman_ero = R.r_spearman_ero;
            gene_results(gi).r_spearman_dw  = R.r_spearman_dw;
        else
            gene_results(gi).r_spearman_ero = R.r_spearman;
            gene_results(gi).r_spearman_dw  = R.r_spearman;
        end
        continue;
    end

    % --- Download / load ISH at native resolution ---
    try
        ish_vol_native = load_allen_ish_grid_native( ...
            ish_cache_dir, gene_id, ish_variant, ...
            ish_raw_permute, ish_raw_flipdim);
    catch ME
        warning('SKIPPING %s: %s', gene_sym, ME.message);
        gene_results(gi).symbol     = gene_sym;
        gene_results(gi).category   = gene_cat;
        gene_results(gi).r_spearman = NaN;
        gene_results(gi).r_pearson  = NaN;
        gene_results(gi).r_voxel    = NaN;
        gene_results(gi).n_bad      = NaN;
        continue;
    end

    % --- Repair bad ISH sections at native resolution ---
    n_native_ap = size(ish_vol_native, 1);
    ish_native_median = nan(n_native_ap, 1);
    for z = 1:n_native_ap
        sl = squeeze(ish_vol_native(z, :, :));
        vals = sl(sl > 0);
        if ~isempty(vals), ish_native_median(z) = median(vals); end
    end
    global_median = nanmedian(ish_native_median);
    bad_thresh = 0.3 * global_median;
    bad_slices = find(ish_native_median < bad_thresh | isnan(ish_native_median));
    good_slices = setdiff(1:n_native_ap, bad_slices);
    ish_native_median_orig = ish_native_median;

    if ~isempty(bad_slices) && numel(good_slices) >= 2
        fprintf('  Repairing %d bad ISH sections...\n', numel(bad_slices));
        for bi = 1:numel(bad_slices)
            bz = bad_slices(bi);
            below = good_slices(good_slices < bz);
            above = good_slices(good_slices > bz);
            if ~isempty(below) && ~isempty(above)
                z_lo = below(end); z_hi = above(1);
                alpha = (bz - z_lo) / (z_hi - z_lo);
                ish_vol_native(bz,:,:) = (1-alpha)*ish_vol_native(z_lo,:,:) + alpha*ish_vol_native(z_hi,:,:);
            elseif ~isempty(below)
                ish_vol_native(bz,:,:) = ish_vol_native(below(end),:,:);
            elseif ~isempty(above)
                ish_vol_native(bz,:,:) = ish_vol_native(above(1),:,:);
            end
        end
    end
    ish_native_median_rep = nan(n_native_ap, 1);
    for z = 1:n_native_ap
        sl = squeeze(ish_vol_native(z,:,:));
        vals = sl(sl > 0);
        if ~isempty(vals), ish_native_median_rep(z) = median(vals); end
    end
    n_bad_sections = numel(bad_slices);

    % --- Diagnostic: native profile ---
    if save_diagnostic_slices
        fh_nat = figure('Color','w','Position',[100 100 1200 450],'visible','off');
        hold on;
        bar(1:n_native_ap, ish_native_median_orig, 'FaceColor',[0.85 0.33 0.10],'EdgeColor','none','DisplayName','Original');
        if ~isempty(bad_slices)
            bar(bad_slices, ish_native_median_orig(bad_slices),'FaceColor',[0.8 0 0],'EdgeColor','k','LineWidth',1.5,'DisplayName',sprintf('Bad (%d)',numel(bad_slices)));
        end
        plot(1:n_native_ap, ish_native_median_rep,'b-o','LineWidth',2,'MarkerSize',4,'MarkerFaceColor','b','DisplayName','After repair');
        yline(bad_thresh,'k--',sprintf('threshold = %.1f',bad_thresh),'LineWidth',1.5,'FontSize',10);
        xlabel('Native ISH slice index (200 um)'); ylabel('Median expression energy');
        title(sprintf('%s — %d bad slices repaired out of %d', gene_sym, n_bad_sections, n_native_ap));
        legend('Location','northeast'); grid on; box on; set(gca,'FontSize',11);
        fb = find(ish_native_median_orig > 0, 1, 'first');
        xlim([max(1,fb-1), n_native_ap+1]); ylim([0 20]);
        exportgraphics(fh_nat, fullfile(out_dir,['Diagnostic_ISH_NativeProfile_' gene_sym '.png']),'Resolution',200);
        close(fh_nat);
    end

    % --- Upsample to 10um ---
    fprintf('  Upsampling ISH to 10um...\n');
    ish_vol_full = imresize3(single(ish_vol_native), allen_full_size, 'linear');
    ish_vol_full(ish_vol_full < 0) = 0;
    clear ish_vol_native
    ish_vol_10um = ish_vol_full(limits(1):limits(2), :, :);
    clear ish_vol_full
    ish_vol_10um(~brainMask) = 0;

    % --- Diagnostic: orientation ---
    if save_diagnostic_slices
        mid = round(size(AllenCrop, 1) / 2);
        fh = figure('Color','w','Position',[100 100 1400 450],'visible','off');
        subplot(1,3,1); imagesc(squeeze(AllenCrop(mid,:,:))); axis image off;
        title(sprintf('Atlas slice %d',mid)); colormap(gca,gray);
        subplot(1,3,2); imagesc(squeeze(ish_vol_10um(mid,:,:))); axis image off;
        title(sprintf('%s ISH slice %d',gene_sym,mid)); colormap(gca,hot);
        subplot(1,3,3);
        atl = squeeze(AllenCrop(mid,:,:));
        bd = (gradient(single(atl))~=0) & (atl>0); [br,bc] = find(bd);
        imagesc(squeeze(ish_vol_10um(mid,:,:))); axis image off;
        title([gene_sym ' + atlas']); colormap(gca,hot);
        hold on; plot(bc,br,'.','Color',[0.5 0.5 0.5],'MarkerSize',0.5);
        sgtitle(['P9 ISH orientation — ' gene_sym]);
        exportgraphics(fh, fullfile(out_dir,['Diagnostic_ISH_Orientation_' gene_sym '.png']),'Resolution',200);
        close(fh);
    end

    % --- Compute L+R ---
    [~, ish_lr_sum] = compute_lr_stats(single(ish_vol_10um));
    clear ish_vol_10um
    if ~isequal(size(ish_lr_sum), size(mean_lr_sum))
        warning('SKIPPING %s: shape mismatch after L+R.', gene_sym);
        gene_results(gi).symbol = gene_sym; gene_results(gi).category = gene_cat;
        gene_results(gi).r_spearman = NaN; gene_results(gi).r_pearson = NaN;
        gene_results(gi).r_voxel = NaN; gene_results(gi).n_bad = n_bad_sections;
        continue;
    end

    % --- Z-score ISH ---
    if strcmp(zscore_method, 'robust')
        ish_z = robust_zscore_in_mask(ish_lr_sum, brainMask_analysis);
    else
        ish_z = zscore_in_mask(ish_lr_sum, brainMask_analysis);
    end
    ish_z(~brainMask_merged) = NaN;
    diff_z = data_z_global - ish_z;

    % --- Voxel-level Pearson ---
    ish_vox = ish_lr_sum(brainMask_analysis);
    valid_vox = ~isnan(data_vox) & ~isnan(ish_vox) & ish_vox > 0;
    if sum(valid_vox) > 100
        r_voxel_pearson = corr(data_vox(valid_vox), ish_vox(valid_vox), 'type', 'Pearson');
    else
        r_voxel_pearson = NaN;
    end

    % --- Per-region ISH means with BOTH methods ---
    ish_vol_masked = ish_lr_sum;
    ish_vol_masked(~brainMask_merged) = NaN;
    ish_vol_masked(1:analysis_slice_range(1)-1,:,:) = NaN;
    ish_vol_masked(analysis_slice_range(2)+1:end,:,:) = NaN;
    ish_z_masked = ish_z;
    ish_z_masked(1:analysis_slice_range(1)-1,:,:) = NaN;
    ish_z_masked(analysis_slice_range(2)+1:end,:,:) = NaN;

    % Method A: eroded means (sparse indices)
    roi_ish_mean_ero = nan(n_rois, 1);
    roi_ish_z_ero    = nan(n_rois, 1);
    for r = 1:n_rois
        if roi_px_ero(r) == 0, continue; end
        vi = ish_vol_masked(roi_indices_ero{r}); valid = ~isnan(vi);
        if any(valid), roi_ish_mean_ero(r) = mean(vi(valid)); end
        viz = ish_z_masked(roi_indices_ero{r}); valid = ~isnan(viz);
        if any(valid), roi_ish_z_ero(r) = mean(viz(valid)); end
    end

    % Method B: distance-weighted means (sparse indices + weights)
    roi_ish_mean_dw = nan(n_rois, 1);
    roi_ish_z_dw    = nan(n_rois, 1);
    for r = 1:n_rois
        if roi_px_raw(r) == 0, continue; end
        w = double(roi_weights_dw{r});
        vi = ish_vol_masked(roi_indices{r}); valid = ~isnan(vi);
        ws = sum(w(valid));
        if ws > 0, roi_ish_mean_dw(r) = sum(vi(valid) .* w(valid)) / ws; end
        viz = ish_z_masked(roi_indices{r}); valid = ~isnan(viz);
        ws = sum(w(valid));
        if ws > 0, roi_ish_z_dw(r) = sum(viz(valid) .* w(valid)) / ws; end
    end

    % Build region tables for both methods
    region_table_ero = table(roi_list(:), roi_acronyms(:), roi_macro(:), ...
        roi_data_mean_ero(:), roi_ish_mean_ero(:), roi_data_z_ero(:), roi_ish_z_ero(:), ...
        roi_px_ero(:), ...
        'VariableNames', {'region','acronym','macro','data_mean','ish_mean','data_zmean','ish_zmean','n_voxels_eroded'});
    region_table_dw = table(roi_list(:), roi_acronyms(:), roi_macro(:), ...
        roi_data_mean_dw(:), roi_ish_mean_dw(:), roi_data_z_dw(:), roi_ish_z_dw(:), ...
        roi_px_raw(:), ...
        'VariableNames', {'region','acronym','macro','data_mean','ish_mean','data_zmean','ish_zmean','n_voxels'});

    % Default region_table for plots = distance-weighted (but both are saved)
    region_table = region_table_dw;

    % --- Region-level correlations for BOTH methods ---
    valid_ero = ~isnan(region_table_ero.data_zmean) & ~isnan(region_table_ero.ish_zmean);
    valid_dw  = ~isnan(region_table_dw.data_zmean)  & ~isnan(region_table_dw.ish_zmean);

    if sum(valid_ero) >= 3
        r_pearson_ero  = corr(region_table_ero.data_zmean(valid_ero), region_table_ero.ish_zmean(valid_ero), 'type','Pearson');
        r_spearman_ero = corr(region_table_ero.data_zmean(valid_ero), region_table_ero.ish_zmean(valid_ero), 'type','Spearman');
    else
        r_pearson_ero = NaN; r_spearman_ero = NaN;
    end
    if sum(valid_dw) >= 3
        r_pearson_dw  = corr(region_table_dw.data_zmean(valid_dw), region_table_dw.ish_zmean(valid_dw), 'type','Pearson');
        r_spearman_dw = corr(region_table_dw.data_zmean(valid_dw), region_table_dw.ish_zmean(valid_dw), 'type','Spearman');
    else
        r_pearson_dw = NaN; r_spearman_dw = NaN;
    end

    % Use distance-weighted as primary for summary
    r_spearman = r_spearman_dw;
    r_pearson  = r_pearson_dw;

    fprintf('  Eroded:   Spearman=%.3f, Pearson=%.3f\n', r_spearman_ero, r_pearson_ero);
    fprintf('  DistWt:   Spearman=%.3f, Pearson=%.3f\n', r_spearman_dw, r_pearson_dw);
    fprintf('  VoxelPearson=%.3f\n', r_voxel_pearson);

    % --- Save gene_result (light — includes BOTH methods) ---
    save(result_path, ...
        'gene_sym', 'gene_id', 'gene_cat', 'gene_desc', ...
        'r_spearman_ero', 'r_pearson_ero', 'r_spearman_dw', 'r_pearson_dw', ...
        'r_spearman', 'r_pearson', 'r_voxel_pearson', ...
        'n_bad_sections', 'ish_native_median_orig', 'ish_native_median_rep', ...
        'region_table_ero', 'region_table_dw', ...
        'analysis_slice_range', 'roi_erode_radius', 'zscore_method');

    % --- Save ISH volume (heavy, separate file) ---
    fprintf('  Saving ISH LR-sum volume...\n');
    save(ish_vol_path, 'ish_lr_sum', '-v7.3');

    % --- Save region table CSV ---
    writetable(region_table, fullfile(out_dir, ['Region_NanoVsISH_Table_' gene_sym '.csv']));

    % --- Collect for summary ---
    gene_results(gi).symbol         = gene_sym;
    gene_results(gi).category       = gene_cat;
    gene_results(gi).r_spearman     = r_spearman;
    gene_results(gi).r_pearson      = r_pearson;
    gene_results(gi).r_voxel        = r_voxel_pearson;
    gene_results(gi).r_spearman_ero = r_spearman_ero;
    gene_results(gi).r_spearman_dw  = r_spearman_dw;
    gene_results(gi).n_bad          = n_bad_sections;

    % --- Comparison video ---
    if generate_comparison_video
        fprintf('  Writing comparison video...\n');
        data_bv = mean_lr_sum(brainMask_merged & ~isnan(mean_lr_sum));
        ish_bv  = ish_lr_sum(brainMask_merged & ~isnan(ish_lr_sum));
        data_pct = [prctile(data_bv,5), prctile(data_bv,90)];
        ish_pct  = [prctile(ish_bv,2),  prctile(ish_bv,98)];
        write_3panel_comparison_video( ...
            mean_lr_sum, ish_lr_sum, diff_z, half_atlas, brainMask_merged, ...
            out_dir, ['cmp_' channel '_vs_ish_' gene_sym '.mp4'], ...
            data_pct, ish_pct, video_clim_diff, ...
            ['Merged ' channel ' LR-sum - ' strrep(p8_merged_tag,'_',' ')], ...
            [gene_sym ' ISH LR-sum - exp ' num2str(gene_id)], ...
            ['Delta z (' channel ' - ' gene_sym ')'], ...
            label_centroids, label_acronyms, label_fontsize, label_color);
    end

    % --- Scatter plot (uses distance-weighted table) ---
    if generate_scatter_plot
        nz = region_table_dw.data_zmean(valid_dw);
        iz = region_table_dw.ish_zmean(valid_dw);
        ml = region_table_dw.macro(valid_dw);
        al = region_table_dw.acronym(valid_dw);

        fig_sc = figure('Color','w','Units','Normalized','Position',[0 0 0.9 0.9],'visible','off');
        hold on;
        mu = unique(ml);
        lims = [min([nz;iz])-0.5, max([nz;iz])+0.5];
        xs_f = linspace(lims(1),lims(2),200)';
        for mi = 1:numel(mu)
            sel = strcmp(ml, mu{mi});
            if isKey(macro_color_map, mu{mi}), mc = macro_color_map(mu{mi}); else, mc = [0.5 0.5 0.5]; end
            scatter(iz(sel), nz(sel), 35, mc, 'filled', 'MarkerFaceAlpha', 0.85, 'HandleVisibility','off');
            if sum(sel) >= 3
                np = numel(find(sel)); p = polyfit(iz(sel),nz(sel),1);
                yf = polyval(p,iz(sel)); SSr = sum((nz(sel)-yf).^2); MSE = SSr/(np-2);
                xm = mean(iz(sel)); Sxx = sum((iz(sel)-xm).^2);
                yff = polyval(p,xs_f); se = sqrt(MSE*(1/np+(xs_f-xm).^2/Sxx));
                tc = tinv(0.975,np-2);
                fill([xs_f;flipud(xs_f)],[yff+tc*se;flipud(yff-tc*se)],mc,'FaceAlpha',0.06,'EdgeColor','none','HandleVisibility','off');
                plot(xs_f,yff,'-','Color',mc,'LineWidth',1.2,'HandleVisibility','off');
                le = sprintf('%s (s=%.2f)',mu{mi},p(1));
            else
                le = mu{mi};
            end
            plot(NaN,NaN,'s','MarkerFaceColor',mc,'MarkerEdgeColor','none','MarkerSize',8,'DisplayName',le);
        end
        plot(lims,lims,'k--','LineWidth',1,'HandleVisibility','off');
        xlim(lims); ylim(lims);
        xlabel(['ISH ' gene_sym ' z-score (per region)']);
        ylabel('Nano LR-sum z-score (per region)');
        title(sprintf(['Nano vs ' gene_sym ' — rS=%.3f  rP=%.3f  rVox=%.3f'], r_spearman, r_pearson, r_voxel_pearson));
        legend('Location','northwest','Interpreter','none','FontSize',8);
        grid on; axis square;
        for k = 1:numel(nz)
            text(iz(k),nz(k),['  ' al{k}],'FontSize',5,'Interpreter','none','Color',[0.3 0.3 0.3]);
        end
        exportgraphics(fig_sc, fullfile(out_dir,['Scatter_NanoVsISH_' gene_sym '.png']),'Resolution',300);
        saveas(fig_sc, fullfile(out_dir,['Scatter_NanoVsISH_' gene_sym '.fig']));
        close(fig_sc);
    end

    % --- Paired bar charts (sorted by nano + sorted by ISH) ---
    if generate_paired_barchart
        for sort_mode = 1:2
            sub = region_table_dw(valid_dw, :);
            if sort_mode == 1
                [~,ord] = sort(sub.data_zmean,'descend'); stag = channel;
            else
                [~,ord] = sort(sub.ish_zmean,'descend');  stag = 'ISH';
            end
            sub = sub(ord,:); ns = height(sub);
            fig_pb = figure('Color','w','Units','Normalized','Position',[0 0 0.6 1],'visible','off');
            bh = barh([sub.data_zmean, sub.ish_zmean],'grouped');
            bh(1).FaceColor = [0.85 0.33 0.10]; bh(1).EdgeColor = 'none';
            bh(2).FaceColor = [0.20 0.40 0.80]; bh(2).EdgeColor = 'none';
            hold on;
            sw = min(15, max(3, round(ns/10)));
            plot(movmean(sub.data_zmean,sw), 1:ns, '-', 'Color',[0.85 0.33 0.10],'LineWidth',2,'HandleVisibility','off');
            plot(movmean(sub.ish_zmean,sw),  1:ns, '-', 'Color',[0.20 0.40 0.80],'LineWidth',2,'HandleVisibility','off');
            yticks(1:ns); yticklabels(sub.acronym);
            set(gca,'YDir','reverse','FontSize',3.5,'TickLabelInterpreter','none');
            xlabel(sprintf('z-score (%s)',zscore_method),'FontSize',11);
            title([channel ' vs ' gene_sym ' ISH - sorted by ' stag],'FontSize',12);
            legend({channel,[gene_sym ' ISH']},'Location','southeast','FontSize',10);
            grid on; box on; ylim([0 ns+1]);
            sfx = ['PairedBar_' channel 'VsISH_sort' stag '_' gene_sym];
            exportgraphics(fig_pb, fullfile(out_dir,[sfx '.png']),'Resolution',300);
            saveas(fig_pb, fullfile(out_dir,[sfx '.fig']));
            close(fig_pb);
        end
    end

    clear ish_lr_sum ish_z diff_z ish_vol_masked ish_z_masked
    fprintf('  Done with %s.\n', gene_sym);
end

%% ====================== CROSS-GENE SUMMARY =============================

fprintf('\n========== Cross-gene summary ==========\n');

% Collect into a table (robust to missing fields)
sym_list     = cell(n_genes, 1);
cat_list     = cell(n_genes, 1);
rs_dw_list   = nan(n_genes, 1);
rs_ero_list  = nan(n_genes, 1);
rp_dw_list   = nan(n_genes, 1);
rp_ero_list  = nan(n_genes, 1);
rv_list      = nan(n_genes, 1);
nb_list      = nan(n_genes, 1);
for gi = 1:n_genes
    if isfield(gene_results, 'symbol') && gi <= numel(gene_results) && ~isempty(gene_results(gi).symbol)
        sym_list{gi} = gene_results(gi).symbol;
        cat_list{gi} = gene_results(gi).category;
        rv_list(gi)  = gene_results(gi).r_voxel;
        nb_list(gi)  = gene_results(gi).n_bad;
        if isfield(gene_results, 'r_spearman_ero') && ~isempty(gene_results(gi).r_spearman_ero)
            rs_ero_list(gi) = gene_results(gi).r_spearman_ero;
            rs_dw_list(gi)  = gene_results(gi).r_spearman_dw;
        else
            rs_ero_list(gi) = gene_results(gi).r_spearman;
            rs_dw_list(gi)  = gene_results(gi).r_spearman;
        end
        % Load per-gene result to get both Pearson values
        gene_sym_gi = gene_results(gi).symbol;
        gene_dir = fullfile(base_root, 'comparisons', ...
            [p8_merged_tag '_vs_ish_' lower(gene_sym_gi) p8_smooth_suffix]);
        rp_file = fullfile(gene_dir, ['gene_result_' gene_sym_gi '.mat']);
        if exist(rp_file, 'file')
            Rp = load(rp_file, 'r_pearson_ero', 'r_pearson_dw');
            if isfield(Rp, 'r_pearson_ero')
                rp_ero_list(gi) = Rp.r_pearson_ero;
                rp_dw_list(gi)  = Rp.r_pearson_dw;
            else
                rp_dw_list(gi)  = gene_results(gi).r_pearson;
                rp_ero_list(gi) = gene_results(gi).r_pearson;
            end
        else
            rp_dw_list(gi)  = gene_results(gi).r_pearson;
            rp_ero_list(gi) = gene_results(gi).r_pearson;
        end
    else
        sym_list{gi} = char(gene_tbl.symbol(gi));
        cat_list{gi} = char(gene_tbl.category(gi));
    end
end

summary_tbl = table(sym_list, cat_list, ...
    rs_dw_list, rs_ero_list, rp_dw_list, rp_ero_list, rv_list, nb_list, ...
    'VariableNames', {'symbol','category', ...
        'r_spearman_dw','r_spearman_ero','r_pearson_dw','r_pearson_ero','r_voxel_pearson','n_bad_sections'});
summary_tbl = sortrows(summary_tbl, 'r_spearman_dw', 'descend', 'MissingPlacement','last');
writetable(summary_tbl, fullfile(summary_dir, 'gene_panel_summary.csv'));
save(fullfile(summary_dir, 'gene_panel_summary.mat'), 'summary_tbl', 'gene_results');

% --- Summary bar chart: Spearman by gene, colored by category ---
valid_s = ~isnan(summary_tbl.r_spearman_dw);
st = summary_tbl(valid_s, :);
n_st = height(st);

cat_colors = containers.Map( ...
    {'AMPAR_core','auxiliary','scaffold','trafficking','excitatory','plasticity',...
     'control_inhib','control_glia','control_struct','control_vasc'}, ...
    {[0.85 0.33 0.10],  ...  % AMPAR_core    — orange (matches nano)
     [0.93 0.69 0.13],  ...  % auxiliary      — golden yellow
     [0.20 0.40 0.80],  ...  % scaffold       — medium blue
     [0.30 0.60 0.85],  ...  % trafficking    — steel blue
     [0.64 0.08 0.18],  ...  % excitatory     — dark red
     [0.80 0.15 0.15],  ...  % plasticity     — bright red
     [0.50 0.50 0.50],  ...  % control_inhib  — medium gray
     [0.70 0.70 0.70],  ...  % control_glia   — light gray
     [0.35 0.35 0.35],  ...  % control_struct — dark gray
     [0.60 0.60 0.60]});     % control_vasc   — gray

% (Removed: single-metric spearman_dw bar chart — redundant with 4-metric chart)

% --- Eroded vs distance-weighted comparison ---
fig_cmp = figure('Visible','off','Color','w','Units','Normalized','Position',[0 0 0.9 0.9]);
bar_ero_dw = [st.r_spearman_ero, st.r_spearman_dw];
bh_cmp = barh(bar_ero_dw, 'grouped');
bh_cmp(1).FaceColor = [0.75 0.68 0.50]; bh_cmp(1).DisplayName = 'Eroded (r=3)';
bh_cmp(2).FaceColor = [0.85 0.33 0.10]; bh_cmp(2).DisplayName = 'Distance-weighted';
yticks(1:n_st); yticklabels(st.symbol);
set(gca,'YDir','reverse','FontSize',7,'TickLabelInterpreter','none');
xlabel('Spearman r (region-level)','FontSize',12);
title('Eroded vs Distance-weighted — ranked by dist-weighted Spearman','FontSize',13);
legend('Location','southeast','FontSize',10);
grid on; box on; xline(0,'k-','LineWidth',0.5);
exportgraphics(fig_cmp, fullfile(summary_dir, 'correlation_eroded_vs_distweight.png'), 'Resolution', 300);
saveas(fig_cmp, fullfile(summary_dir, 'correlation_eroded_vs_distweight.fig'));
close(fig_cmp);

% --- All 4 metrics side by side ---
fig_4m = figure('Visible','off','Color','w','Units','Normalized','Position',[0 0 0.9 0.9]);
bar_data_5 = [st.r_spearman_dw, st.r_spearman_ero, st.r_pearson_dw, st.r_pearson_ero, st.r_voxel_pearson];
bh4 = barh(bar_data_5, 'grouped');
bh4(1).FaceColor = [0.85 0.33 0.10]; bh4(1).DisplayName = 'Spearman (dist-wt)';
bh4(2).FaceColor = [0.75 0.68 0.50]; bh4(2).DisplayName = 'Spearman (eroded)';
bh4(3).FaceColor = [0.20 0.40 0.80]; bh4(3).DisplayName = 'Pearson (dist-wt)';
bh4(4).FaceColor = [0.30 0.60 0.85]; bh4(4).DisplayName = 'Pearson (eroded)';
bh4(5).FaceColor = [0.64 0.08 0.18]; bh4(5).DisplayName = 'Pearson (voxel)';
yticks(1:n_st); yticklabels(st.symbol);
set(gca,'YDir','reverse','FontSize',7,'TickLabelInterpreter','none');
xlabel(['Correlation with ' channel ' LR-sum'],'FontSize',12);
title('Cross-gene: 5 correlation metrics — ranked by dist-weighted Spearman','FontSize',13);
legend('Location','southeast','FontSize',10);
grid on; box on; xline(0,'k-','LineWidth',0.5);
exportgraphics(fig_4m, fullfile(summary_dir, 'correlation_barchart_4metrics.png'), 'Resolution', 300);
saveas(fig_4m, fullfile(summary_dir, 'correlation_barchart_4metrics.fig'));
close(fig_4m);

% --- Category violin plots (3 metrics) ---

% Split auxiliary into forebrain vs cerebellar/other subgroups
split_auxiliary = true;  % toggle: false = single "Auxiliary" category
aux_forebrain_genes = {'Cacng8','Cacng3','Cnih2','Cnih3','Grm5'};

if split_auxiliary
    % Remap 'auxiliary' entries in st to sub-categories
    for k = 1:n_st
        if strcmp(st.category{k}, 'auxiliary')
            if ismember(st.symbol{k}, aux_forebrain_genes)
                st.category{k} = 'aux_forebrain';
            else
                st.category{k} = 'aux_other';
            end
        end
    end
    cats_ordered = {'AMPAR_core','aux_forebrain','aux_other','scaffold','trafficking','excitatory','plasticity','control_inhib','control_glia','control_struct'};
    cats_display = {'AMPAR core','Aux forebrain','Aux other','Scaffold','Trafficking','Excitatory','Plasticity','Ctrl inhib','Ctrl glia','Ctrl struct'};
    cat_col_list = { ...
        [0.85 0.33 0.10], ...  % AMPAR_core    — orange
        [0.93 0.69 0.13], ...  % aux_forebrain — golden yellow
        [0.75 0.68 0.50], ...  % aux_other     — grayish yellow
        [0.20 0.40 0.80], ...  % scaffold       — medium blue
        [0.30 0.60 0.85], ...  % trafficking    — steel blue
        [0.64 0.08 0.18], ...  % excitatory     — dark red
        [0.80 0.15 0.15], ...  % plasticity     — bright red
        [0.50 0.50 0.50], ...  % control_inhib  — medium gray
        [0.70 0.70 0.70], ...  % control_glia   — light gray
        [0.35 0.35 0.35]};     % control_struct — dark gray
else
    cats_ordered = {'AMPAR_core','auxiliary','scaffold','trafficking','excitatory','plasticity','control_inhib','control_glia','control_struct'};
    cats_display = {'AMPAR core','Auxiliary','Scaffold','Trafficking','Excitatory','Plasticity','Ctrl inhib','Ctrl glia','Ctrl struct'};
    cat_col_list = { ...
        [0.85 0.33 0.10], ...  % AMPAR_core    — orange
        [0.93 0.69 0.13], ...  % auxiliary      — golden yellow
        [0.20 0.40 0.80], ...  % scaffold       — medium blue
        [0.30 0.60 0.85], ...  % trafficking    — steel blue
        [0.64 0.08 0.18], ...  % excitatory     — dark red
        [0.80 0.15 0.15], ...  % plasticity     — bright red
        [0.50 0.50 0.50], ...  % control_inhib  — medium gray
        [0.70 0.70 0.70], ...  % control_glia   — light gray
        [0.35 0.35 0.35]};     % control_struct — dark gray
end
n_cats = numel(cats_ordered);

metrics = {'r_spearman_dw', 'r_spearman_ero', 'r_voxel_pearson'};
metric_titles = {'Spearman (distance-weighted)', 'Spearman (eroded)', 'Pearson (voxel-level)'};
metric_fnames = {'violin_spearman_dw', 'violin_spearman_ero', 'violin_voxel_pearson'};

for mi = 1:numel(metrics)
    metric_vals = st.(metrics{mi});

    % Build cell array of distributions per category, then sort by median descending
    distrs_unsorted = cell(1, n_cats);
    cat_medians = nan(1, n_cats);
    for ci = 1:n_cats
        sel = strcmp(st.category, cats_ordered{ci});
        vals_ci = metric_vals(sel);
        vals_ci = vals_ci(~isnan(vals_ci));
        distrs_unsorted{ci} = vals_ci;
        if ~isempty(vals_ci), cat_medians(ci) = median(vals_ci); end
    end
    [~, sort_order] = sort(cat_medians, 'descend', 'MissingPlacement', 'last');
    distrs = distrs_unsorted(sort_order);
    colors = cat_col_list(sort_order);
    cats_display_sorted = cats_display(sort_order);
    cats_ordered_sorted = cats_ordered(sort_order);

    fig_vio = figure('Visible', 'off', 'Color', 'w', 'Position', [50 50 1100 500]);
    ax = axes(fig_vio);

    inputadata.inputdistrs = distrs;
    inputpars.n_distribs       = n_cats;
    inputpars.dirstrcenters    = 1:n_cats;
    inputpars.boxplotwidth     = 0.15;
    inputpars.boxplotlinewidth = 1.5;
    inputpars.densityplotwidth = 0.35;
    inputpars.xlimtouse        = [0 n_cats + 1];
    inputpars.yimtouse         = [-0.4, 1.0];
    inputpars.scatterjitter    = 0.15;
    inputpars.scatteralpha     = 0.7;
    inputpars.scattersize      = 25;
    inputpars.xtickslabelvector = cats_display_sorted;
    inputpars.distrcolors      = colors;
    inputpars.distralpha       = 0.5;
    inputpars.xlabelstring     = '';
    inputpars.ylabelstring     = 'Correlation r';
    inputpars.titlestring      = metric_titles{mi};
    inputpars.boolscatteron    = true;
    inputpars.ks_bandwidth     = 0.12;
    inputpars.inputaxh         = ax;

    plot_violinplot(inputadata, inputpars);
    yline(0, 'k--', 'LineWidth', 0.8);

    % Add gene name labels next to each dot
    for ci = 1:n_cats
        sel = strcmp(st.category, cats_ordered_sorted{ci});
        sel_idx = find(sel);
        vals_ci = metric_vals(sel);
        syms_ci = st.symbol(sel);
        for ki = 1:numel(vals_ci)
            if ~isnan(vals_ci(ki))
                text(ax, ci + 0.22, vals_ci(ki), syms_ci{ki}, ...
                    'FontSize', 4.5, 'Interpreter', 'none', 'Color', [0.3 0.3 0.3]);
            end
        end
    end

    ylim(ax, [-0.4, 1.0]);
    set(ax, 'FontSize', 10);
    xtickangle(ax, 25);

    % One-way ANOVA: effect of category on correlation
    cat_idx_anova = zeros(numel(metric_vals), 1);
    for ci = 1:n_cats
        cat_idx_anova(strcmp(st.category, cats_ordered_sorted{ci})) = ci;
    end
    keep_anova = cat_idx_anova > 0 & ~isnan(metric_vals);
    [p_anova, tbl_anova] = anova1(metric_vals(keep_anova), ...
        cat_idx_anova(keep_anova), 'off');
    F_anova = tbl_anova{2, 5};  % F-statistic from the ANOVA table
    title(ax, sprintf('%s (ANOVA p = %.2e, F = %.1f)', ...
        metric_titles{mi}, p_anova, F_anova), 'FontSize', 12);

    ylabel(ax, 'Correlation r', 'FontSize', 12);
    grid(ax, 'on'); box(ax, 'on');

    exportgraphics(fig_vio, fullfile(summary_dir, [metric_fnames{mi} '.png']), ...
        'Resolution', 300);
    saveas(fig_vio, fullfile(summary_dir, [metric_fnames{mi} '.fig']));
    close(fig_vio);
end

fprintf('\nP9 batch complete. Summary in: %s\n', summary_dir);
fprintf('Per-gene results in: %s_vs_ish_<gene>%s/\n', ...
    fullfile(base_root,'comparisons',p8_merged_tag), p8_smooth_suffix);


%% Local functions ======================================================

function ish_vol = load_allen_ish_grid_native(cache_dir, exp_id, variant, ...
    ish_raw_permute, ish_raw_flipdim)
    local_zip = fullfile(cache_dir, sprintf('%d_grid.zip', exp_id));
    local_mhd = fullfile(cache_dir, sprintf('%d_%s.mhd', exp_id, variant));
    local_raw = fullfile(cache_dir, sprintf('%d_%s.raw', exp_id, variant));
    if ~exist(local_mhd, 'file') || ~exist(local_raw, 'file')
        if ~exist(local_zip, 'file')
            url = sprintf('http://api.brain-map.org/grid_data/download/%d', exp_id);
            fprintf('  Downloading ISH grid (exp %d)...\n', exp_id);
            try
                websave(local_zip, url);
            catch ME
                error('Failed to download exp %d: %s', exp_id, ME.message);
            end
        end
        extracted = unzip(local_zip, cache_dir);
        for k = 1:numel(extracted)
            [~,nm,ext] = fileparts(extracted{k});
            if strcmpi(nm,variant) && any(strcmpi(ext,{'.mhd','.raw'}))
                movefile(extracted{k}, fullfile(cache_dir, sprintf('%d_%s%s',exp_id,variant,ext)));
            end
        end
    end
    [dims, dtype] = parse_mhd_header(local_mhd);
    fid = fopen(local_raw, 'r', 'ieee-le');
    raw_data = fread(fid, prod(dims), ['*' dtype]);
    fclose(fid);
    ish_vol = permute(reshape(raw_data, dims(:)'), ish_raw_permute);
    for d = 1:3
        if ish_raw_flipdim(d), ish_vol = flip(ish_vol, d); end
    end
end

function [dims, dtype] = parse_mhd_header(mhd_path)
    fid = fopen(mhd_path, 'r');
    c = onCleanup(@() fclose(fid));
    dims = []; dtype = '';
    while ~feof(fid)
        line = fgetl(fid);
        if ~ischar(line), break; end
        if startsWith(line,'DimSize')
            dims = sscanf(regexprep(line,'.*=\s*',''), '%d');
        elseif startsWith(line,'ElementType')
            et = strtrim(regexprep(line,'.*=\s*',''));
            switch et
                case 'MET_FLOAT',  dtype='single';
                case 'MET_DOUBLE', dtype='double';
                case 'MET_UCHAR',  dtype='uint8';
                case 'MET_SHORT',  dtype='int16';
                case 'MET_USHORT', dtype='uint16';
                otherwise, error('Unsupported: %s',et);
            end
        end
    end
end

function vol_z = zscore_in_mask(vol, mask)
    vals = vol(mask); vals = vals(~isnan(vals));
    mu = mean(vals); sg = std(vals);
    if sg==0||isnan(sg), vol_z=nan(size(vol),'like',vol); return; end
    vol_z = (vol-mu)/sg; vol_z(~mask) = NaN;
end

function vol_z = robust_zscore_in_mask(vol, mask)
    vals = vol(mask); vals = vals(~isnan(vals));
    med = median(vals); mad_val = median(abs(vals-med));
    sg = mad_val * 1.4826;
    if sg==0||isnan(sg), vol_z=nan(size(vol),'like',vol); return; end
    vol_z = (vol-med)/sg; vol_z(~mask) = NaN;
end

function write_3panel_comparison_video(vol_a, vol_b, vol_diff, atlas_vol, ...
    brain_mask, save_dir, video_filename, ...
    clim_a, clim_b, clim_diff, ...
    title_a, title_b, title_diff, ...
    label_centroids, label_acronyms, label_fontsize, label_color)
    if ~exist(save_dir,'dir'), mkdir(save_dir); end
    vidObj = VideoWriter(fullfile(save_dir,video_filename),'MPEG-4');
    vidObj.FrameRate = 15; vidObj.Quality = 95; open(vidObj);
    n_slices = size(vol_a,1);
    n_rois_lbl = numel(label_acronyms);
    have_labels = ~isempty(label_centroids) && n_rois_lbl > 0;
    try, diff_cmap = get_color2color_colormap([0 0 1],[1 0 0]); catch, diff_cmap = jet; end
    for j = 1:n_slices
        if sum(sum(brain_mask(j,:,:)))==0, continue; end
        atlasim = single(squeeze(atlas_vol(j,:,:)));
        bnd = gradient(atlasim)~=0 & atlasim>1;
        [row,col] = ind2sub(size(atlasim), find(bnd));
        fh = figure('visible','off','units','normalized','outerposition',[0 0 1 1],'Color','k');
        set(fh,'InvertHardcopy','off');
        panels = {vol_a, vol_b, vol_diff};
        clims = {clim_a, clim_b, clim_diff};
        titles = {title_a, title_b, title_diff};
        cmaps = {hot, hot, diff_cmap};
        for p = 1:3
            subplot(1,3,p);
            d2 = squeeze(panels{p}(j,:,:));
            himg = imagesc(d2); clim(clims{p});
            set(himg,'AlphaData',squeeze(brain_mask(j,:,1:size(d2,2))));
            colormap(gca,cmaps{p});
            ax=gca; ax.Color='k'; axis equal off; hold on;
            line(col,row,'Marker','.','LineStyle','none','Color',[0.66 0.66 0.66],'MarkerSize',0.5);
            xlim([0 size(d2,2)]);
            title(titles{p},'Color','w','FontSize',11,'Interpreter','none');
            cb=colorbar; cb.Color='w'; cb.Label.Color='w'; cb.Label.FontSize=9;
            if have_labels && p==3
                for r=1:n_rois_lbl
                    cy=label_centroids(j,r,1); cx=label_centroids(j,r,2);
                    if isnan(cy)||isnan(cx), continue; end
                    text(cx,cy,label_acronyms{r},'Color',label_color,'FontSize',label_fontsize,...
                        'FontWeight','bold','HorizontalAlignment','center','VerticalAlignment','middle',...
                        'Interpreter','none','Clipping','on');
                end
            end
        end
        sgtitle(['Slice # ' num2str(j)],'Color','w','FontSize',14);
        writeVideo(vidObj, getframe(fh)); close(fh);
        if mod(j,100)==0, fprintf('    Frame %d/%d\n',j,n_slices); end
    end
    close(vidObj);
end

function [roi_names, roi_acronyms, roi_macro] = build_stru_leaves_by_macro_local(allenDir, macro_divi_list)
    parc_tbl = readtable(fullfile(allenDir,'parcellation_term.csv'),'TextType','string');
    mem_tbl  = readtable(fullfile(allenDir,'parcellation_term_set_membership.csv'),'TextType','string');
    ids = cellstr(parc_tbl.identifier);
    names = cellstr(parc_tbl.name);
    acronyms = cellstr(parc_tbl.acronym);
    parents = cellstr(parc_tbl.parent_identifier);
    mem_labels = cellstr(mem_tbl.parcellation_term_label);
    mem_sets = cellstr(mem_tbl.parcellation_term_set_label);
    set_suffix = regexprep(mem_sets,'^AllenCCF-Ontology-2017-','');
    label_to_set = containers.Map(mem_labels, set_suffix);
    children_map = containers.Map('KeyType','char','ValueType','any');
    for k=1:numel(ids)
        p=parents{k}; if isempty(p), continue; end
        if isKey(children_map,p), children_map(p)=[children_map(p),k]; else, children_map(p)=k; end
    end
    row_set = repmat({'NONE'},numel(ids),1);
    term_labels = cellstr(parc_tbl.label);
    for k=1:numel(ids)
        if isKey(label_to_set,term_labels{k}), row_set{k}=label_to_set(term_labels{k}); end
    end
    acronym_to_idx = containers.Map(acronyms, 1:numel(acronyms));
    roi_names={}; roi_acronyms={}; roi_macro={};
    for mi=1:numel(macro_divi_list)
        mac=macro_divi_list{mi};
        if ~isKey(acronym_to_idx,mac), continue; end
        root_k=acronym_to_idx(mac);
        stack=root_k; all_desc=[];
        while ~isempty(stack)
            k=stack(end); stack(end)=[];
            if isKey(children_map,ids{k})
                ch=children_map(ids{k}); all_desc=[all_desc,ch]; stack=[stack,ch]; %#ok<AGROW>
            end
        end
        all_desc=unique(all_desc);
        is_stru=cellfun(@(s)strcmp(s,'STRU'),row_set(all_desc));
        stru_desc=all_desc(is_stru);
        is_leaf=false(size(stru_desc));
        for li=1:numel(stru_desc)
            kk=stru_desc(li); sub_stack=[];
            if isKey(children_map,ids{kk}), sub_stack=children_map(ids{kk}); end
            hsc=false;
            while ~isempty(sub_stack)
                ss=sub_stack(end); sub_stack(end)=[];
                if strcmp(row_set{ss},'STRU'), hsc=true; break; end
                if isKey(children_map,ids{ss}), sub_stack=[sub_stack,children_map(ids{ss})]; end %#ok<AGROW>
            end
            is_leaf(li)=~hsc;
        end
        leaves=stru_desc(is_leaf);
        for li=1:numel(leaves)
            kk=leaves(li);
            roi_names{end+1}=names{kk}; roi_acronyms{end+1}=acronyms{kk}; roi_macro{end+1}=mac; %#ok<AGROW>
        end
    end
    roi_names=roi_names(:); roi_acronyms=roi_acronyms(:); roi_macro=roi_macro(:);
end
