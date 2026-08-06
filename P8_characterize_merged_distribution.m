clear all
close all
clc

% /// Pipeline script #8: characterize pooled distribution of a signal channel ///
% Channel is configurable via the `channel` parameter ('nano' = surface GluA1,
% 'auto' = autofluorescence control). Pools mice across groups (default:
% naive + rws) into one merged cohort, computes the hemispheric SUM of the
% chosen channel, and produces:
%   (1) a merged-group LR-sum video (+ reliability t-score panel)
%   (2) per-mouse LR-sum videos
%   (3) a regional bar chart ranked by MEAN sum intensity per Allen region
% behavior mice are excluded (noisy data).

%% User-defined parameters

% Groups to pool (reference group listed first — other groups get linearly
% aligned to the reference via per-slice profile polyfit, as in P7bis).
groups_to_merge = {'naive', 'rws'};

% Smoothing (match P7bis)
apply_smoothing = false;
smooth_sigma    = 5.0;

% Output toggles
generate_avg_video        = true;
generate_individual_videos = false;
generate_region_barchart   = true;
generate_zscore_video      = true;   % raw + z-scored threshold videos

% Enrichment thresholds
raw_threshold    = 1.5;   % threshold in original LR-sum intensity units
zscore_threshold = 0.0;   % threshold in z-score units (SD above brain mean)
show_threshold_raw    = true;   % highlight enriched regions in raw barplot
show_threshold_zscore = true;   % highlight enriched regions in z-scored barplot
show_threshold_raw_video    = true;   % highlight above-threshold contour in raw video
show_threshold_zscore_video = true;   % highlight above-threshold contour in z-scored video

% Per-mouse + SEM (heavy block, opt-in). When true: loops mice through
% lr_sum_merged using the cached sparse ROI masks, computes per-region mean
% per mouse, renders additional *_withSEM barplots (mean +/- SEM errorbars,
% optional individual dots). The cheap cohort-mean path stays untouched.
compute_per_mouse_sem = true;
show_per_mouse_dots   = true;    % overlay individual mouse points on SEM bars

% Channel to analyze: 'nano' (default, surface GluA1) or 'auto' (autofluorescence
% control). Loads <channel>_4d_normalized.mat from each cohort folder, and rolls
% the channel into merged_tag so all output filenames carry the suffix and
% nano/auto runs do not collide.
channel = 'nano'; % 'auto'

% Region mask caching (set true to force re-generation)
force_recompute_masks = false;  % set back to false after first run with new z-score code

% Region aggregation method: 'eroded' or 'distweight' or 'both'
% 'eroded'     = binary erosion by roi_erode_radius (original method)
% 'distweight' = distance-weighted mean (center=1, border=0)
% 'both'       = compute both, plot distance-weighted, save both to table
region_agg_method = 'distweight';
dist_weight_power = 4;  % exponent for distance weighting (1=linear, 2=quadratic, 4=quartic)

% Region-acronym overlay on videos (grand average + per-mouse)
show_slice_labels       = true;   % single kill-switch (both videos)
label_fontsize          = 7;
label_color             = [0.85 0.85 0.85];
label_smooth_halfwindow = 5;      % +/- N slices moving average on centroids
label_min_px_per_slice  = 150;    % min region pixels per slice to render label

% Ontology-based ROI list for bar charts + video labels
% Set of DIVI-level macros whose STRU-leaf descendants become the ROI list.
% Comment out entries to reduce the panel count.
macro_divi_list = { ...
    'Isocortex', ...   % Isocortex
    'OLF', ...         % Olfactory areas
    'HPF', ...         % Hippocampal formation
    'CTXsp', ...       % Cortical subplate (claustrum, amygdala, endopiriform, ...)
    'STR', ...         % Striatum
    'PAL', ...         % Pallidum (incl. septum, BST, SI, ...)
    'TH', ...          % Thalamus
    'HY', ...          % Hypothalamus
    'MB'};             % Midbrain

% Per-ROI 3D erosion (voxels) applied before computing the mean for bar
% charts. Prevents border contamination between adjacent areas. Centroids
% for video labels use the UN-eroded mask so placement stays natural.
roi_erode_radius = 3;

% Slice range for quantitative analyses (1-indexed coronal pseudoslices).
% Excludes noisy anterior/posterior extremes from region means and bar
% charts. Videos still render the full range.
analysis_slice_range = [100, 700];

% Paths
base_root = 'D:\sep_histology\data\';
% merged_tag includes the channel so nano and auto outputs go to separate folders
merged_tag = ['merged_' strjoin(groups_to_merge, '_') '_' channel];
out_dir    = fullfile(base_root, 'comparisons', merged_tag);
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

% Suffix appended to every output filename, reflecting smoothing state
if apply_smoothing
    smooth_suffix = sprintf('_smooth%g', smooth_sigma);
else
    smooth_suffix = '_nosmooth';
end

%% Allen atlas setup (identical to P7bis lines 66-77)

allenDir = 'D:\sep_histology\data\atlas';
addpath(allenDir);
AllenFile = fullfile(allenDir, 'annotation_10.nii.gz');
AllenVol = niftiread(AllenFile);
limits = [180 1079];
AllenCrop = AllenVol(limits(1):limits(2), :, :);
brainMask = AllenCrop > 0;
half_atlas = AllenCrop(:, :, 1:end);
clear AllenVol

%% Load normalized data for each group

G = numel(groups_to_merge);
raw_vols              = cell(G, 1);
raw_masks             = cell(G, 1);
med_profiles          = cell(G, 1);
mouse_names_per_group = cell(G, 1);

% Channel-aware filenames + dynamic field name
norm_var_name    = [channel '_4d_normalized'];                  % e.g. 'nano_4d_normalized'
norm_filename    = [channel '_4d_normalized.mat'];
bkgmask_filename = [channel '_4d_normalized_bkgmask.mat'];

for gi = 1:G
    gname = groups_to_merge{gi};
    gdir  = fullfile(base_root, gname);
    fprintf('Loading group %s (channel=%s) ...\n', gname, channel);

    S_vol  = load(fullfile(gdir, norm_filename), norm_var_name, 'current_mice');
    S_mask = load(fullfile(gdir, bkgmask_filename), 'recomputed_bkg_mask_4d');

    raw_vols{gi}  = S_vol.(norm_var_name);
    raw_masks{gi} = S_mask.recomputed_bkg_mask_4d;

    if isfield(S_vol, 'current_mice') && ~isempty(S_vol.current_mice)
        mouse_names_per_group{gi} = S_vol.current_mice;
    else
        % Fallback: synthetic names if P6bis didn't save current_mice
        n_m = size(raw_vols{gi}, 4);
        mouse_names_per_group{gi} = arrayfun( ...
            @(k) sprintf('%s_m%d', gname, k), 1:n_m, 'UniformOutput', false);
    end
    fprintf('  Loaded %d mice.\n', size(raw_vols{gi}, 4));
end

%% Per-slice background-subtracted mean profiles (used for cross-group alignment)

for gi = 1:G
    vol  = raw_vols{gi};
    mask = raw_masks{gi};
    n_slices = size(vol, 1);
    n_mice   = size(vol, 4);
    med_prof = nan(n_slices, n_mice);

    fprintf('Computing per-slice profile for %s ...\n', groups_to_merge{gi});
    for m = 1:n_mice
        for z = 1:n_slices
            img = squeeze(vol(z, :, :, m));
            bgm = squeeze(mask(z, :, :, m));
            med_prof(z, m) = nanmean(img(~bgm));
        end
    end
    med_profiles{gi} = med_prof;
end

%% NaN-tolerant 3D Gaussian smoothing (per-mouse, match P7bis)

if apply_smoothing
    fprintf('Applying NaN-robust 3D Gaussian smoothing (sigma = %.1f)...\n', smooth_sigma);
    nan_smooth_3d = @(v, sig) imgaussfilt3(fillmissing(v, 'constant', 0), sig) ./ ...
        imgaussfilt3(double(~isnan(v)), sig);

    for gi = 1:G
        vol  = raw_vols{gi};
        mask = raw_masks{gi};
        fprintf('  Group %s (%d mice)\n', groups_to_merge{gi}, size(vol, 4));
        for m = 1:size(vol, 4)
            tic
            v = vol(:, :, :, m);
            b = mask(:, :, :, m);
            is_valid_tissue = brainMask & ~b;
            v(~is_valid_tissue) = NaN;
            vs = nan_smooth_3d(v, smooth_sigma);
            vs(~is_valid_tissue) = 0;
            vs(isnan(vs)) = 0;
            vol(:, :, :, m) = vs;
            toc
        end
        raw_vols{gi} = vol;
    end
    fprintf('  Smoothing complete.\n');
end

%% Linearly align each non-reference group to the reference group

ref_gi = 1;  % first group listed is the reference
interest_region = 200:700;

mean_ref_profile = nanmean(med_profiles{ref_gi}, 2);
y_target_full    = mean_ref_profile(interest_region);

norm_profiles    = cell(G, 1);
norm_profiles{ref_gi} = med_profiles{ref_gi};

for gi = 1:G
    if gi == ref_gi, continue; end
    mp = nanmean(med_profiles{gi}, 2);
    x_source = mp(interest_region);
    valid = ~isnan(x_source) & ~isnan(y_target_full);
    p = polyfit(x_source(valid), y_target_full(valid), 1);
    slope = p(1); intercept = p(2);
    fprintf('Aligning %s -> %s: slope=%.4f, intercept=%.4f\n', ...
        groups_to_merge{gi}, groups_to_merge{ref_gi}, slope, intercept);
    raw_vols{gi}     = raw_vols{gi} .* slope + intercept;
    norm_profiles{gi} = med_profiles{gi} .* slope + intercept;
end

%% Compute a common median factor (unit scaling, match P7bis lines 206-211)

interest_region_bis = 300:500;
fact_list = zeros(G, 1);
for gi = 1:G
    fact_list(gi) = nanmean(nanmean(norm_profiles{gi}(interest_region_bis, :), 1));
end
common_fact = mean(fact_list);
fprintf('Common normalization factor: %.4f\n', common_fact);

%% Merge into a single 4D cohort

data_4d_merged    = cat(4, raw_vols{:}) ./ common_fact;
raw_mask_4d_merged = cat(4, raw_masks{:});
merged_mouse_names = [mouse_names_per_group{:}];
fprintf('Merged cohort: %d mice total\n', size(data_4d_merged, 4));

clear raw_vols raw_masks

%% Compute LR sum (discard diff)

[~, lr_sum_merged] = compute_lr_stats(data_4d_merged);

% Per-voxel summary over cohort
mean_lr_sum = nanmean(abs(lr_sum_merged), 4);

% Reliability: t = mean / SEM
sem_lr_sum = nanstd(abs(lr_sum_merged), [], 4) ./ ...
    sqrt(sum(~isnan(lr_sum_merged), 4));
sem_lr_sum(sem_lr_sum == 0) = NaN;
t_lr_sum = mean_lr_sum ./ sem_lr_sum;

%% Build hemisphere masks (match P7bis lines 281-295)

n_half = size(lr_sum_merged, 3);
n_full = size(raw_mask_4d_merged, 3);

bg_L = raw_mask_4d_merged(:, :, 1:n_half, :);
bg_R = raw_mask_4d_merged(:, :, (n_full - n_half + 1):end, :);
mask_bg_merged = logical(bg_L | flip(bg_R, 3));  % 4D per-mouse hemisphere bg

brainMask_cropped = brainMask(:, :, 1:n_half);
tissue_3d_merged  = any(~mask_bg_merged, 4);
brainMask_merged  = brainMask_cropped & tissue_3d_merged;  % 3D cohort hemi mask

%% Whole-volume z-scoring
% Compute z-score of the mean LR-sum across all brain voxels within the
% analysis slice range (consistent with P9). This gives a principled scale
% where zscore_threshold (default 2) = "N SD above brain-wide mean".

brainMask_analysis = brainMask_merged;
brainMask_analysis(1:analysis_slice_range(1)-1, :, :) = false;
brainMask_analysis(analysis_slice_range(2)+1:end, :, :) = false;

zscore_lr_sum = zscore_in_mask(mean_lr_sum, brainMask_analysis);
fprintf('Z-scored LR-sum: mean=%.4f, std=%.4f (computed over %d brain voxels)\n', ...
    mean(mean_lr_sum(brainMask_analysis), 'omitnan'), ...
    std(mean_lr_sum(brainMask_analysis), 'omitnan'), ...
    nnz(brainMask_analysis));
fprintf('  Voxels above z=%g threshold: %d (%.1f%%)\n', ...
    zscore_threshold, ...
    nnz(zscore_lr_sum > zscore_threshold & brainMask_analysis), ...
    100 * nnz(zscore_lr_sum > zscore_threshold & brainMask_analysis) / nnz(brainMask_analysis));

%% Cache merged nano outputs for downstream scripts (P9 ISH comparison, ...)
% Persist the merged half-width mean LR-sum, the reliability t-score, and
% the cohort brain mask so scripts downstream can load them in <1 s instead
% of rerunning the full P8 pipeline. Saved in out_dir with the merged tag
% + smooth suffix in the filename so nosmooth and smoothed runs coexist.

half_width = n_half;  % alias — downstream scripts expect this name

data_cache_path = fullfile(out_dir, ...
    ['mean_lr_sum_' merged_tag smooth_suffix '.mat']);
fprintf('Caching merged %s outputs to: %s\n', channel, data_cache_path);
save(data_cache_path, ...
    'mean_lr_sum', ...       % 3D, half-width (row, col, half_w): per-voxel mean L+R
    't_lr_sum', ...          % 3D, half-width: reliability t = mean / SEM
    'zscore_lr_sum', ...     % 3D, half-width: whole-brain z-scored mean LR-sum
    'brainMask_merged', ...  % 3D, half-width: cohort tissue mask
    'brainMask_analysis', ...% 3D, half-width: analysis-slice-restricted brain mask
    'zscore_threshold', ...  % scalar: enrichment threshold in z-score units
    'half_width', ...        % scalar: index bound for left hemisphere
    'n_half', ...            % alias for half_width used in some downstream code
    'merged_tag', ...        % string: e.g. 'merged_naive_rws'
    'smooth_suffix', ...     % string: '_smooth5' or '_nosmooth'
    'groups_to_merge', ...   % cell: the groups pooled into this cache
    'apply_smoothing', ...   % bool: whether spatial smoothing was applied
    'smooth_sigma', ...      % scalar: smoothing sigma (if applied)
    '-v7.3');

%% Build ontology-based ROI list (shared by bar charts + video labels)
% Walks the Allen ontology and collects, for each DIVI in macro_divi_list,
% all its STRU "leaves" (STRU nodes whose own descendants are all SUBS, i.e.
% cortical layers / sub-layers). Leaf-only selection avoids double-counting
% umbrella STRU nodes like SSp against their STRU children (SSp-bfd, ...).

roi_list      = {};
roi_acronyms  = {};
roi_macro     = {};  % macro acronym per ROI
macro_groups  = struct();  % macro_groups.(MACRO) = [row indices into roi_list]

need_roi_list = generate_region_barchart || (generate_avg_video && show_slice_labels) ...
    || (generate_individual_videos && show_slice_labels);

if need_roi_list
    fprintf('Building STRU-leaf ROI list from ontology for %d macros...\n', numel(macro_divi_list));
    [roi_list, roi_acronyms, roi_macro] = build_stru_leaves_by_macro(allenDir, macro_divi_list);
    fprintf('  Collected %d STRU-leaf ROIs across %d macros.\n', numel(roi_list), numel(macro_divi_list));
    for r = 1:numel(roi_list)
        mkey = matlab.lang.makeValidName(roi_macro{r});
        if ~isfield(macro_groups, mkey)
            macro_groups.(mkey) = [];
        end
        macro_groups.(mkey)(end+1) = r; %#ok<SAGROW>
    end
end

%% (Optional) Precompute region-acronym label positions for the videos
% Builds a per-slice smoothed centroid table for each ROI using the shared
% ontology-built roi_list. Skipped entirely if show_slice_labels = false.

label_centroids = [];
label_acronyms  = {};
centroid_cache_path = fullfile(out_dir, ...
    ['label_centroids_' merged_tag smooth_suffix '.mat']);

if (generate_avg_video || generate_individual_videos || generate_zscore_video) && show_slice_labels && ~isempty(roi_list)
    if exist(centroid_cache_path, 'file') && ~force_recompute_masks
        fprintf('Loading cached label centroids: %s\n', centroid_cache_path);
        S_lc = load(centroid_cache_path);
        label_centroids = S_lc.label_centroids;
        label_acronyms  = S_lc.label_acronyms;
        fprintf('  %d regions loaded from cache.\n', numel(label_acronyms));
    else
        fprintf('Precomputing region label positions...\n');

        label_acronyms = roi_acronyms;

        [~, ~, atlas_n_width] = size(AllenCrop);
        atlas_half_width = floor(atlas_n_width / 2);

        n_slices_lbl   = size(AllenCrop, 1);
        n_rois_lbl     = numel(roi_list);
        label_centroids = nan(n_slices_lbl, n_rois_lbl, 2);

        fprintf('  Computing per-slice centroids for %d regions...\n', n_rois_lbl);
        for r = 1:n_rois_lbl
            try
                m_temp = get_allen_region_mask(allenDir, AllenCrop, roi_list(r), brainMask, '');
                if size(m_temp, 3) > atlas_half_width
                    m_temp = m_temp(:, :, 1:atlas_half_width);
                end
                for z = 1:n_slices_lbl
                    slice_mask = squeeze(m_temp(z, :, :));
                    [ri, ci] = find(slice_mask);
                    if numel(ri) >= label_min_px_per_slice
                        label_centroids(z, r, 1) = mean(ri);
                        label_centroids(z, r, 2) = mean(ci);
                    end
                end
            catch
                warning('Could not build label mask for: %s', roi_list{r});
            end
        end

        fprintf('  Smoothing centroid trajectories (halfwindow = %d)...\n', label_smooth_halfwindow);
        win = 2 * label_smooth_halfwindow + 1;
        for r = 1:n_rois_lbl
            for d = 1:2
                v = label_centroids(:, r, d);
                present = ~isnan(v);
                s = movmean(v, win, 'omitnan');
                s(~present) = NaN;
                label_centroids(:, r, d) = s;
            end
        end

        save(centroid_cache_path, 'label_centroids', 'label_acronyms', '-v7.3');
        fprintf('  Centroid cache saved: %s\n', centroid_cache_path);
    end
end

%% Merged-group LR-sum video
% We reuse write_lr_video (unmodified). Its left panel takes clim_values and
% its right panel takes 2*clim_values. We send mean_lr_sum to the left and
% t_lr_sum (reliability) to the right — two useful views in one pass.
% If show_slice_labels is true, we dispatch to a local labeled variant
% (see bottom of file) that adds region-acronym text overlays.

if generate_avg_video
    fprintf('Writing merged-group sum video...\n');
    merged_tag_display = strrep(merged_tag, '_', ' ');
    if show_slice_labels
        write_lr_video_labeled( ...
            mean_lr_sum, ...
            t_lr_sum, ...
            half_atlas, ...
            brainMask_merged, ...
            out_dir, ...
            [['lr_sum_' channel '_'] merged_tag smooth_suffix '.mp4'], ...
            [0, 5], ...
            merged_tag_display, ...
            ['Mean LR sum (' channel ')'], ...
            'Reliability t = mean/SEM', ...
            label_centroids, ...
            label_acronyms, ...
            label_fontsize, ...
            label_color);
    else
        write_lr_video( ...
            mean_lr_sum, ...
            t_lr_sum, ...
            half_atlas, ...
            brainMask_merged, ...
            out_dir, ...
            [['lr_sum_' channel '_'] merged_tag smooth_suffix '.mp4'], ...
            [0, 5], ...
            merged_tag_display, ...
            ['Mean LR sum (' channel ')'], ...
            'Reliability t = mean/SEM');
    end
end

%% Per-mouse LR-sum videos
% One video per mouse, single full-screen panel showing the per-voxel
% LR sum (L + R), styled like the merged-group video. Rendering is done
% by a local function at the bottom of this file so we don't modify any
% shared helpers.

if generate_individual_videos
    fprintf('Writing per-mouse sum videos...\n');
    clim_sum_indiv = [0, 10];
    n_mice_merged  = size(lr_sum_merged, 4);
    for mi = 1:n_mice_merged
        mouse_name = merged_mouse_names{mi};
        mouse_tag  = matlab.lang.makeValidName(mouse_name);
        vid_name   = ['Individual_LR_Sum_' mouse_tag smooth_suffix '.mp4'];
        per_mouse_brain = ~squeeze(mask_bg_merged(:, :, :, mi)) & brainMask_cropped;
        if show_slice_labels
            write_single_mouse_sum_video_labeled( ...
                squeeze(abs(lr_sum_merged(:, :, :, mi))), ...
                per_mouse_brain, ...
                half_atlas, ...
                out_dir, ...
                vid_name, ...
                clim_sum_indiv, ...
                mouse_name, ...
                label_centroids, label_acronyms, label_fontsize, label_color);
        else
            write_single_mouse_sum_video( ...
                squeeze(abs(lr_sum_merged(:, :, :, mi))), ...
                per_mouse_brain, ...
                half_atlas, ...
                out_dir, ...
                vid_name, ...
                clim_sum_indiv, ...
                mouse_name);
        end
    end
end

%% Regional bar charts: per-macro panels + cross-DIVI summary
% Uses the ontology-built roi_list (STRU leaves per DIVI macro) and
% per-region 3D erosion to avoid contamination from adjacent-area borders.

if generate_region_barchart && ~isempty(roi_list)
    fprintf('Building regional mean-sum bar charts...\n');

    [~, ~, n_width] = size(AllenCrop);
    half_width = floor(n_width / 2);
    atlas_left = AllenCrop(:, :, 1:half_width);
    n_rois = numel(roi_list);

    region_cache_path_p8 = fullfile(out_dir, ...
        ['Region_MeanSum_Table_' merged_tag smooth_suffix '.mat']);

    if exist(region_cache_path_p8, 'file') && ~force_recompute_masks
        fprintf('  Loading cached region table: %s\n', region_cache_path_p8);
        S_rc = load(region_cache_path_p8);
        region_rank = S_rc.region_rank;
        [~, reorder] = ismember(roi_list, region_rank.region);
        roi_mean   = nan(n_rois, 1);
        roi_px_raw = zeros(n_rois, 1);
        roi_px_ero = zeros(n_rois, 1);
        valid_ro = reorder > 0;
        roi_mean(valid_ro)   = region_rank.mean_lr_sum(reorder(valid_ro));
        roi_px_raw(valid_ro) = region_rank.n_voxels_raw(reorder(valid_ro));
        if ismember('n_voxels_eroded', region_rank.Properties.VariableNames)
            roi_px_ero(valid_ro) = region_rank.n_voxels_eroded(reorder(valid_ro));
        else
            roi_px_ero = roi_px_raw;
        end
        fprintf('  %d regions loaded from cache.\n', sum(valid_ro));
    else
        roi_mean_ero = nan(n_rois, 1);
        roi_mean_dw  = nan(n_rois, 1);
        roi_px_raw   = zeros(n_rois, 1);
        roi_px_ero   = zeros(n_rois, 1);

        if roi_erode_radius > 0
            se_erode = strel('sphere', roi_erode_radius);
        end

        mean_vol_masked = mean_lr_sum;
        mean_vol_masked(~brainMask_merged) = NaN;
        mean_vol_masked(1:analysis_slice_range(1)-1, :, :) = NaN;
        mean_vol_masked(analysis_slice_range(2)+1:end, :, :) = NaN;

        fprintf('  Computing region means for %d ROIs (method: %s, slices %d-%d)...\n', ...
            n_rois, region_agg_method, analysis_slice_range(1), analysis_slice_range(2));
        for r = 1:n_rois
            try
                m_raw = get_allen_region_mask(allenDir, AllenCrop, roi_list(r), brainMask, '');
                if size(m_raw, 3) > half_width
                    m_raw = m_raw(:, :, 1:half_width);
                end
                roi_px_raw(r) = nnz(m_raw);

                % Eroded mean
                if roi_erode_radius > 0
                    m_ero = imerode(m_raw, se_erode);
                    if nnz(m_ero) == 0, m_ero = m_raw; end
                else
                    m_ero = m_raw;
                end
                roi_px_ero(r) = nnz(m_ero);
                vals_ero = mean_vol_masked(m_ero);
                vals_ero = vals_ero(~isnan(vals_ero));
                if ~isempty(vals_ero), roi_mean_ero(r) = mean(vals_ero); end

                % Distance-weighted mean (with power exponent)
                dist = single(bwdist(~m_raw));
                dist(~m_raw) = 0;
                mx = max(dist(:));
                if mx > 0, dist = (dist / mx) .^ dist_weight_power; end
                mv = m_raw & ~isnan(mean_vol_masked);
                wsum = sum(dist(mv));
                if wsum > 0
                    roi_mean_dw(r) = sum(mean_vol_masked(mv) .* double(dist(mv))) / wsum;
                end
            catch
                warning('Could not compute mean for region: %s', roi_list{r});
            end
        end

        % Select which method to use for plots
        if strcmp(region_agg_method, 'distweight')
            roi_mean = roi_mean_dw;
        else
            roi_mean = roi_mean_ero;
        end

        % Save table with both methods
        region_rank = table(roi_list(:), roi_acronyms(:), roi_macro(:), ...
            roi_mean(:), roi_mean_ero(:), roi_mean_dw(:), ...
            roi_px_raw(:), roi_px_ero(:), ...
            'VariableNames', {'region', 'acronym', 'macro', 'mean_lr_sum', ...
                              'mean_lr_sum_ero', 'mean_lr_sum_dw', ...
                              'n_voxels_raw', 'n_voxels_eroded'});
        region_rank = sortrows(region_rank, 'mean_lr_sum', 'descend', 'MissingPlacement', 'last');
        writetable(region_rank, fullfile(out_dir, ...
            ['Region_MeanSum_Table_' merged_tag smooth_suffix '.csv']));
        save(region_cache_path_p8, 'region_rank');
        fprintf('  Region table cached (both methods).\n');
    end

    % ---- Shared: z-scored region means + per-region reliability (t-score) ----
    % Z-scored region means (linear transform → arithmetic on roi_mean)
    brain_vals = mean_lr_sum(brainMask_analysis);
    brain_vals = brain_vals(~isnan(brain_vals));
    brain_mu    = mean(brain_vals);
    brain_sigma = std(brain_vals);
    roi_mean_z  = (roi_mean - brain_mu) / brain_sigma;

    % Per-region mean reliability (t-score). Use sparse mask cache if
    % available (fast); otherwise compute on the fly (slow).
    roi_mean_t = nan(n_rois, 1);
    t_vol_masked = t_lr_sum;
    t_vol_masked(~brainMask_analysis) = NaN;

    masks_cache_path = fullfile(out_dir, ...
        sprintf('roi_masks_r%d_pw%d_slices%d-%d.mat', ...
        roi_erode_radius, dist_weight_power, ...
        analysis_slice_range(1), analysis_slice_range(2)));

    if exist(masks_cache_path, 'file')
        fprintf('  Loading sparse ROI masks for t-score computation...\n');
        M = load(masks_cache_path);
        for r = 1:n_rois
            if strcmp(region_agg_method, 'distweight')
                idx_m = M.roi_indices{r};
                w   = M.roi_weights_dw{r};
            else
                idx_m = M.roi_indices_ero{r};
                w   = ones(numel(M.roi_indices_ero{r}), 1, 'single');
            end
            if isempty(idx_m), continue; end
            vals = t_vol_masked(idx_m);
            valid = ~isnan(vals);
            if ~any(valid), continue; end
            roi_mean_t(r) = sum(vals(valid) .* double(w(valid))) / sum(double(w(valid)));
        end
        clear M
    else
        fprintf('  Computing per-region t-score means (no sparse cache found)...\n');
        for r = 1:n_rois
            try
                m_raw = get_allen_region_mask(allenDir, AllenCrop, roi_list(r), brainMask, '');
                if size(m_raw, 3) > half_width
                    m_raw = m_raw(:, :, 1:half_width);
                end
                if strcmp(region_agg_method, 'distweight')
                    dist = single(bwdist(~m_raw));
                    dist(~m_raw) = 0;
                    mx = max(dist(:));
                    if mx > 0, dist = (dist / mx) .^ dist_weight_power; end
                    mv = m_raw & ~isnan(t_vol_masked);
                    wsum = sum(dist(mv));
                    if wsum > 0
                        roi_mean_t(r) = sum(t_vol_masked(mv) .* double(dist(mv))) / wsum;
                    end
                else
                    if roi_erode_radius > 0
                        m_ero = imerode(m_raw, strel('sphere', roi_erode_radius));
                        if nnz(m_ero) == 0, m_ero = m_raw; end
                    else
                        m_ero = m_raw;
                    end
                    vals = t_vol_masked(m_ero);
                    vals = vals(~isnan(vals));
                    if ~isempty(vals), roi_mean_t(r) = mean(vals); end
                end
            catch
            end
        end
    end
    fprintf('  Reliability range: t = [%.1f, %.1f]\n', ...
        min(roi_mean_t, [], 'omitnan'), max(roi_mean_t, [], 'omitnan'));

    % Shared colormap for reliability: grayscale (darker = more reliable)
    c_map_rel = flipud(gray(256));
    c_map_rel = c_map_rel(1:200, :);  % avoid pure white at the low end
    t_max_clamp = 20;  % clamp range for colormap

    % ---- Multi-panel figure: one subplot per DIVI, bars = STRU leaves ----
    macros_present = macro_divi_list(ismember(cellfun( ...
        @matlab.lang.makeValidName, macro_divi_list, 'UniformOutput', false), ...
        fieldnames(macro_groups)));
    n_macros_plot = numel(macros_present);
    if n_macros_plot > 0
        n_cols_grid = 3;
        n_rows_grid = ceil(n_macros_plot / n_cols_grid);
        fig_by_macro = figure( ...
            'Visible', 'off', ...
            'Name', ['Region_MeanSum_BarByMacro_' merged_tag smooth_suffix], ...
            'Color', 'w', 'Units', 'Normalized', ...
            'Position', [0 0 1 1]);

        for pi = 1:n_macros_plot
            macro_ac  = macros_present{pi};
            macro_key = matlab.lang.makeValidName(macro_ac);
            idx = macro_groups.(macro_key);
            sub_vals = roi_mean(idx);
            sub_acro = roi_acronyms(idx);
            sub_t    = roi_mean_t(idx);

            keep = ~isnan(sub_vals);
            sub_vals = sub_vals(keep);
            sub_acro = sub_acro(keep);
            sub_t    = sub_t(keep);
            [sub_vals, order] = sort(sub_vals, 'ascend');
            sub_acro = sub_acro(order);
            sub_t    = sub_t(order);

            subplot(n_rows_grid, n_cols_grid, pi);
            if isempty(sub_vals)
                title([macro_ac ' (no data)'], 'Interpreter', 'none');
                axis off;
                continue;
            end

            b = barh(sub_vals);
            b.FaceColor = 'flat';
            % Color by reliability t-score (grayscale): darker = more reliable
            t_clamp = max(0, min(sub_t, t_max_clamp));
            c_idx = round((t_clamp / t_max_clamp) * (size(c_map_rel, 1) - 1)) + 1;
            c_idx(isnan(c_idx)) = 1;
            for k = 1:numel(sub_vals)
                b.CData(k, :) = c_map_rel(c_idx(k), :);
            end

            lbl_fs = max(4, min(7, 120 / numel(sub_acro)));
            yticks(1:numel(sub_acro));
            set(gca, 'FontSize', lbl_fs);
            hold on;
            if show_threshold_raw
                above = sub_vals >= raw_threshold;
                highlight_enriched_labels(gca, sub_acro, above, lbl_fs);
                xline(raw_threshold, '--', sprintf('thr = %g', raw_threshold), ...
                    'Color', [0.85 0.5 0.0], 'LineWidth', 1.5, ...
                    'LabelOrientation', 'horizontal', 'FontSize', 7);
            else
                yticklabels(sub_acro);
            end

            xlabel(['Mean LR-sum ' channel]);
            title([macro_ac ' (n = ' num2str(numel(sub_vals)) ')'], 'Interpreter', 'none');
            grid on;
            ylim([0 numel(sub_acro) + 1]);
        end
        sgtitle(['Mean ' channel ' LR-sum per STRU leaf, by DIVI macro - ' ...
            strrep(merged_tag, '_', ' ')], 'FontSize', 13, 'FontWeight', 'bold');

        % Shared colorbar for reliability
        cb_ax_raw = axes('Position', [0.93 0.08 0.015 0.35], 'Visible', 'off');
        colormap(cb_ax_raw, c_map_rel);
        cb_raw = colorbar(cb_ax_raw, 'Location', 'eastoutside');
        cb_raw.Label.String = 'Reliability (t = mean/SEM)';
        cb_raw.Label.FontSize = 9;
        clim(cb_ax_raw, [0 t_max_clamp]);

        set(fig_by_macro, 'InvertHardcopy', 'off');
        saveas(fig_by_macro, fullfile(out_dir, ...
            ['Region_MeanSum_BarByMacro_' merged_tag smooth_suffix '.fig']));
        exportgraphics(fig_by_macro, fullfile(out_dir, ...
            ['Region_MeanSum_BarByMacro_' merged_tag smooth_suffix '.png']), ...
            'Resolution', 300);
    end

    % ---- Cross-DIVI summary: one bar per macro (raw + z-scored) ----
    % For each macro, pool the eroded voxels of all its STRU leaves and compute
    % a single weighted mean. Two figures: raw intensity and z-scored, with
    % the same reliability coloring and threshold highlighting as ByMacro.
    macro_names = fieldnames(macro_groups);
    n_macros = numel(macro_names);
    macro_mean      = nan(n_macros, 1);
    macro_mean_z    = nan(n_macros, 1);
    macro_mean_rel  = nan(n_macros, 1);
    macro_display   = cell(n_macros, 1);
    for mi = 1:n_macros
        mkey = macro_names{mi};
        idx  = macro_groups.(mkey);
        macro_display{mi} = mkey;
        for aa = 1:numel(macro_divi_list)
            if strcmp(matlab.lang.makeValidName(macro_divi_list{aa}), mkey)
                macro_display{mi} = macro_divi_list{aa};
                break;
            end
        end
        vals   = roi_mean(idx);
        vals_z = roi_mean_z(idx);
        vals_t = roi_mean_t(idx);
        wts    = roi_px_ero(idx);
        keep_k = ~isnan(vals) & wts > 0;
        if any(keep_k)
            macro_mean(mi)     = sum(vals(keep_k)   .* wts(keep_k)) / sum(wts(keep_k));
            macro_mean_z(mi)   = sum(vals_z(keep_k) .* wts(keep_k)) / sum(wts(keep_k));
            keep_t = keep_k & ~isnan(vals_t);
            if any(keep_t)
                macro_mean_rel(mi) = sum(vals_t(keep_t) .* wts(keep_t)) / sum(wts(keep_t));
            end
        end
    end

    % Helper: plot one cross-DIVI bar chart
    cross_divi_configs = struct( ...
        'vals',   {macro_mean,           macro_mean_z}, ...
        'thr',    {raw_threshold,        zscore_threshold}, ...
        'flag',   {show_threshold_raw,   show_threshold_zscore}, ...
        'xlabel', {['Mean LR-sum ' channel],   'Z-scored mean LR-sum'}, ...
        'title',  {['DIVI-level mean ' channel ' sum - ' strrep(merged_tag,'_',' ')], ...
                   ['DIVI-level z-scored ' channel ' sum - ' strrep(merged_tag,'_',' ')]}, ...
        'suffix', {['Region_MeanSum_BarAcrossDivi_' merged_tag smooth_suffix], ...
                   ['Region_Zscore_BarAcrossDivi_' merged_tag smooth_suffix]});

    for ci = 1:2
        cfg = cross_divi_configs(ci);
        [s_vals, s_ord] = sort(cfg.vals, 'ascend');
        s_labels = macro_display(s_ord);
        s_rel    = macro_mean_rel(s_ord);
        keep_m = ~isnan(s_vals);
        s_vals   = s_vals(keep_m);
        s_labels = s_labels(keep_m);
        s_rel    = s_rel(keep_m);

        fig_cd = figure('Visible', 'off', 'Color', 'w', 'Units', 'Normalized', 'Position', [0 0 0.45 0.6]);
        if ~isempty(s_vals)
            b = barh(s_vals);
            b.FaceColor = 'flat';
            % Color by reliability (grayscale)
            t_c = max(0, min(s_rel, t_max_clamp));
            c_i = round((t_c / t_max_clamp) * (size(c_map_rel,1)-1)) + 1;
            c_i(isnan(c_i)) = 1;
            for k = 1:numel(s_vals)
                b.CData(k,:) = c_map_rel(c_i(k),:);
            end
            yticks(1:numel(s_labels));
            set(gca, 'FontSize', 9);
            hold on;
            if cfg.flag
                above = s_vals >= cfg.thr;
                highlight_enriched_labels(gca, s_labels, above, 9);
                xline(cfg.thr, '--', sprintf('thr = %g', cfg.thr), ...
                    'Color', [0.85 0.5 0.0], 'LineWidth', 1.5, ...
                    'LabelOrientation', 'horizontal', 'FontSize', 8);
            else
                yticklabels(s_labels);
            end
            xlabel(cfg.xlabel);
            title(cfg.title, 'Interpreter', 'none');
            grid on;
            ylim([0 numel(s_labels) + 1]);
        end
        set(fig_cd, 'InvertHardcopy', 'off');
        saveas(fig_cd, fullfile(out_dir, [cfg.suffix '.fig']));
        exportgraphics(fig_cd, fullfile(out_dir, [cfg.suffix '.png']), 'Resolution', 300);
    end

    fprintf('Bar charts saved to: %s\n', out_dir);

    % ---- Z-scored bar chart: threshold line + reliability coloring ----
    % roi_mean_z and roi_mean_t already computed above.
    if n_macros_plot > 0
        fig_by_macro_z = figure( ...
            'Visible', 'off', ...
            'Name', ['Region_Zscore_BarByMacro_' merged_tag smooth_suffix], ...
            'Color', 'w', 'Units', 'Normalized', ...
            'Position', [0 0 1 1]);

        % Same grayscale reliability colormap as the raw barplot

        for pi = 1:n_macros_plot
            macro_ac  = macros_present{pi};
            macro_key = matlab.lang.makeValidName(macro_ac);
            idx = macro_groups.(macro_key);
            sub_vals = roi_mean_z(idx);
            sub_acro = roi_acronyms(idx);
            sub_t    = roi_mean_t(idx);

            keep = ~isnan(sub_vals);
            sub_vals = sub_vals(keep);
            sub_acro = sub_acro(keep);
            sub_t    = sub_t(keep);
            [sub_vals, order] = sort(sub_vals, 'ascend');
            sub_acro = sub_acro(order);
            sub_t    = sub_t(order);

            ax = subplot(n_rows_grid, n_cols_grid, pi);
            if isempty(sub_vals)
                title([macro_ac ' (no data)'], 'Interpreter', 'none');
                axis off;
                continue;
            end

            b = barh(sub_vals);
            b.FaceColor = 'flat';
            % Color by reliability t-score (grayscale): darker = more reliable
            t_clamp = max(0, min(sub_t, t_max_clamp));
            c_idx = round((t_clamp / t_max_clamp) * (size(c_map_rel, 1) - 1)) + 1;
            c_idx(isnan(c_idx)) = 1;
            for k = 1:numel(sub_vals)
                b.CData(k, :) = c_map_rel(c_idx(k), :);
            end

            lbl_fs = max(4, min(7, 120 / numel(sub_acro)));
            yticks(1:numel(sub_acro));
            set(gca, 'FontSize', lbl_fs);
            hold on;
            if show_threshold_zscore
                above = sub_vals >= zscore_threshold;
                highlight_enriched_labels(gca, sub_acro, above, lbl_fs);
                xline(zscore_threshold, '--', sprintf('z = %g', zscore_threshold), ...
                    'Color', [0.85 0.5 0.0], 'LineWidth', 1.5, ...
                    'LabelOrientation', 'horizontal', 'FontSize', 7);
            else
                yticklabels(sub_acro);
            end

            xlabel('Z-scored mean LR-sum');
            title([macro_ac ' (n = ' num2str(numel(sub_vals)) ')'], 'Interpreter', 'none');
            grid on;
            ylim([0 numel(sub_acro) + 1]);
        end
        sgtitle(['Z-scored ' channel ' LR-sum per STRU leaf, by DIVI macro - ' ...
            strrep(merged_tag, '_', ' ')], 'FontSize', 13, 'FontWeight', 'bold');

        % Shared colorbar for reliability
        cb_ax_z = axes('Position', [0.93 0.08 0.015 0.35], 'Visible', 'off');
        colormap(cb_ax_z, c_map_rel);
        cb_z = colorbar(cb_ax_z, 'Location', 'eastoutside');
        cb_z.Label.String = 'Reliability (t = mean/SEM)';
        cb_z.Label.FontSize = 9;
        clim(cb_ax_z, [0 t_max_clamp]);

        set(fig_by_macro_z, 'InvertHardcopy', 'off');
        saveas(fig_by_macro_z, fullfile(out_dir, ...
            ['Region_Zscore_BarByMacro_' merged_tag smooth_suffix '.fig']));
        exportgraphics(fig_by_macro_z, fullfile(out_dir, ...
            ['Region_Zscore_BarByMacro_' merged_tag smooth_suffix '.png']), ...
            'Resolution', 300);
        fprintf('  Z-scored bar chart saved.\n');
    end
end

%% Per-mouse + SEM barplots (additive, gated by compute_per_mouse_sem)
% Loops the per-mouse 4D LR-sum already in memory (lr_sum_merged), uses the
% cached sparse ROI masks built by P9, and produces _withSEM variants of
% all four barplots (BarByMacro raw/z, BarAcrossDivi raw/z) with mean +/-
% SEM errorbars and optional per-mouse dots. The small per-mouse-per-region
% matrix is cached so figure tweaks don't re-loop the 4D.

if compute_per_mouse_sem && generate_region_barchart && ~isempty(roi_list)
    fprintf('Computing per-mouse per-region means (compute_per_mouse_sem = true)...\n');

    pm_cache_path = fullfile(out_dir, ...
        ['per_mouse_region_means_' merged_tag smooth_suffix '.mat']);
    n_mice_pm = size(lr_sum_merged, 4);

    recompute_pm = true;
    if exist(pm_cache_path, 'file') && ~force_recompute_masks
        fprintf('  Loading cached per-mouse means: %s\n', pm_cache_path);
        S_pm = load(pm_cache_path);
        if isfield(S_pm, 'per_mouse_mean_dw') && ...
                size(S_pm.per_mouse_mean_dw, 1) == n_mice_pm && ...
                size(S_pm.per_mouse_mean_dw, 2) == n_rois
            per_mouse_mean_dw  = S_pm.per_mouse_mean_dw;
            per_mouse_mean_ero = S_pm.per_mouse_mean_ero;
            per_mouse_mice     = S_pm.per_mouse_mice;
            recompute_pm = false;
        else
            warning(['Cached per-mouse matrix is missing or shape-mismatched ' ...
                '(expected %d mice x %d regions). Recomputing.'], ...
                n_mice_pm, n_rois);
        end
        clear S_pm
    end

    if recompute_pm
        if ~exist(masks_cache_path, 'file')
            error(['Sparse ROI mask cache not found at %s. ' ...
                   'Run P9 once to build it (or set force_recompute_masks=true ' ...
                   'so P8 builds the t-score path).'], masks_cache_path);
        end
        fprintf('  Loading sparse ROI masks: %s\n', masks_cache_path);
        Mpm = load(masks_cache_path);

        per_mouse_mean_dw  = nan(n_mice_pm, n_rois);
        per_mouse_mean_ero = nan(n_mice_pm, n_rois);

        fprintf('  Looping %d mice x %d regions...\n', n_mice_pm, n_rois);
        for mi = 1:n_mice_pm
            tic_pm = tic;
            vol_i = abs(lr_sum_merged(:, :, :, mi));
            vol_i(~brainMask_merged) = NaN;
            vol_i(1:analysis_slice_range(1)-1, :, :) = NaN;
            vol_i(analysis_slice_range(2)+1:end, :, :) = NaN;

            for r = 1:n_rois
                % Distance-weighted mean
                idx_dw = Mpm.roi_indices{r};
                if ~isempty(idx_dw)
                    vals = vol_i(idx_dw);
                    w_dw = double(Mpm.roi_weights_dw{r});
                    valid = ~isnan(vals);
                    if any(valid)
                        wsum = sum(w_dw(valid));
                        if wsum > 0
                            per_mouse_mean_dw(mi, r) = ...
                                sum(double(vals(valid)) .* w_dw(valid)) / wsum;
                        end
                    end
                end
                % Eroded uniform mean
                idx_e = Mpm.roi_indices_ero{r};
                if ~isempty(idx_e)
                    vals = vol_i(idx_e);
                    valid = ~isnan(vals);
                    if any(valid)
                        per_mouse_mean_ero(mi, r) = mean(double(vals(valid)));
                    end
                end
            end
            fprintf('    mouse %d/%d (%.1fs)\n', mi, n_mice_pm, toc(tic_pm));
        end
        per_mouse_mice = merged_mouse_names;
        save(pm_cache_path, 'per_mouse_mean_dw', 'per_mouse_mean_ero', ...
            'per_mouse_mice', 'roi_list', 'roi_acronyms', 'roi_macro', ...
            'analysis_slice_range', 'dist_weight_power', 'roi_erode_radius', ...
            '-v7.3');
        fprintf('  Per-mouse means cached: %s\n', pm_cache_path);
        clear Mpm
    end

    % Pick which method to plot (matches region_agg_method choice elsewhere)
    if strcmp(region_agg_method, 'distweight')
        per_mouse_mean = per_mouse_mean_dw;
    else
        per_mouse_mean = per_mouse_mean_ero;
    end

    % Per-region cohort mean (from per-mouse) and SEM
    pm_n    = sum(~isnan(per_mouse_mean), 1)';
    pm_mean = mean(per_mouse_mean, 1, 'omitnan')';
    pm_sem  = std(per_mouse_mean, 0, 1, 'omitnan')' ./ sqrt(max(pm_n, 1));
    pm_sem(pm_n < 2) = NaN;

    % Z-scored per-mouse: apply the same linear transform the cohort uses
    % so cohort-mean bar heights stay comparable.
    per_mouse_mean_z = (per_mouse_mean - brain_mu) / brain_sigma;
    pm_mean_z = mean(per_mouse_mean_z, 1, 'omitnan')';
    pm_sem_z  = std(per_mouse_mean_z, 0, 1, 'omitnan')' ./ sqrt(max(pm_n, 1));
    pm_sem_z(pm_n < 2) = NaN;

    % Sanity check vs the existing cohort-mean bar heights
    valid_sc = ~isnan(roi_mean) & ~isnan(pm_mean);
    if any(valid_sc)
        rel_err = abs(pm_mean(valid_sc) - roi_mean(valid_sc)) ./ ...
            max(abs(roi_mean(valid_sc)), eps);
        fprintf('  Sanity (per-mouse mean vs cohort bar): median rel diff = %.4f, max = %.4f\n', ...
            median(rel_err), max(rel_err));
        if max(rel_err) > 0.10
            warning(['Per-mouse mean differs from cohort bar by up to %.1f%%. ' ...
                'Likely due to NaN voxels differing per mouse vs cohort. ' ...
                'Inspect before publishing figures.'], 100 * max(rel_err));
        end
    end

    % ---------- Per-macro By-Macro figures with SEM (raw + z-scored) ----------
    pm_configs = struct( ...
        'vals',   {pm_mean,             pm_mean_z}, ...
        'pm_mat', {per_mouse_mean,      per_mouse_mean_z}, ...
        'sem',    {pm_sem,              pm_sem_z}, ...
        'thr',    {raw_threshold,       zscore_threshold}, ...
        'flag',   {show_threshold_raw,  show_threshold_zscore}, ...
        'xlabel', {['Mean LR-sum ' channel ' (mean +/- SEM)'], 'Z-scored mean LR-sum (mean +/- SEM)'}, ...
        'sgttl',  {['Mean ' channel ' LR-sum per STRU leaf (mean +/- SEM, n=' num2str(n_mice_pm) ' mice) - ' strrep(merged_tag,'_',' ')], ...
                   ['Z-scored ' channel ' LR-sum per STRU leaf (mean +/- SEM, n=' num2str(n_mice_pm) ' mice) - ' strrep(merged_tag,'_',' ')]}, ...
        'fname',  {['Region_MeanSum_BarByMacro_' merged_tag smooth_suffix '_withSEM'], ...
                   ['Region_Zscore_BarByMacro_' merged_tag smooth_suffix '_withSEM']});

    for ci_pm = 1:numel(pm_configs)
        cfg = pm_configs(ci_pm);
        if n_macros_plot == 0, continue; end

        fig_pm = figure('Visible', 'off', 'Name', cfg.fname, 'Color', 'w', ...
            'Units', 'Normalized', 'Position', [0 0 1 1]);

        for pi = 1:n_macros_plot
            macro_ac  = macros_present{pi};
            macro_key = matlab.lang.makeValidName(macro_ac);
            idx = macro_groups.(macro_key);

            sub_vals = cfg.vals(idx);
            sub_sem  = cfg.sem(idx);
            sub_acro = roi_acronyms(idx);
            sub_t    = roi_mean_t(idx);
            sub_pm   = cfg.pm_mat(:, idx);

            keep = ~isnan(sub_vals);
            sub_vals = sub_vals(keep);
            sub_sem  = sub_sem(keep);
            sub_acro = sub_acro(keep);
            sub_t    = sub_t(keep);
            sub_pm   = sub_pm(:, keep);
            [sub_vals, order] = sort(sub_vals, 'ascend');
            sub_sem  = sub_sem(order);
            sub_acro = sub_acro(order);
            sub_t    = sub_t(order);
            sub_pm   = sub_pm(:, order);

            subplot(n_rows_grid, n_cols_grid, pi);
            if isempty(sub_vals)
                title([macro_ac ' (no data)'], 'Interpreter', 'none');
                axis off;
                continue;
            end

            b = barh(sub_vals);
            b.FaceColor = 'flat';
            t_clamp = max(0, min(sub_t, t_max_clamp));
            c_idx = round((t_clamp / t_max_clamp) * (size(c_map_rel, 1) - 1)) + 1;
            c_idx(isnan(c_idx)) = 1;
            for k = 1:numel(sub_vals)
                b.CData(k, :) = c_map_rel(c_idx(k), :);
            end

            lbl_fs = max(4, min(7, 120 / numel(sub_acro)));
            yticks(1:numel(sub_acro));
            set(gca, 'FontSize', lbl_fs);
            hold on;

            % Horizontal SEM errorbars
            errorbar(sub_vals, 1:numel(sub_vals), sub_sem, 'horizontal', ...
                'Color', 'k', 'LineStyle', 'none', 'LineWidth', 0.5, 'CapSize', 3);

            % Per-mouse dots overlay
            if show_per_mouse_dots
                n_pm_loc = size(sub_pm, 1);
                for k = 1:numel(sub_vals)
                    yk = k + (rand(1, n_pm_loc) - 0.5) * 0.3;
                    scatter(sub_pm(:, k), yk, 8, [0.3 0.3 0.3], 'filled', ...
                        'MarkerFaceAlpha', 0.5);
                end
            end

            if cfg.flag
                above = sub_vals >= cfg.thr;
                highlight_enriched_labels(gca, sub_acro, above, lbl_fs);
                xline(cfg.thr, '--', sprintf('thr = %g', cfg.thr), ...
                    'Color', [0.85 0.5 0.0], 'LineWidth', 1.5, ...
                    'LabelOrientation', 'horizontal', 'FontSize', 7);
            else
                yticklabels(sub_acro);
            end

            xlabel(cfg.xlabel);
            title([macro_ac ' (n = ' num2str(numel(sub_vals)) ')'], 'Interpreter', 'none');
            grid on;
            ylim([0 numel(sub_acro) + 1]);
        end
        sgtitle(cfg.sgttl, 'FontSize', 13, 'FontWeight', 'bold');

        % Shared reliability colorbar
        cb_ax = axes('Position', [0.93 0.08 0.015 0.35], 'Visible', 'off');
        colormap(cb_ax, c_map_rel);
        cb = colorbar(cb_ax, 'Location', 'eastoutside');
        cb.Label.String = 'Reliability (t = mean/SEM)';
        cb.Label.FontSize = 9;
        clim(cb_ax, [0 t_max_clamp]);

        set(fig_pm, 'InvertHardcopy', 'off');
        saveas(fig_pm, fullfile(out_dir, [cfg.fname '.fig']));
        exportgraphics(fig_pm, fullfile(out_dir, [cfg.fname '.png']), 'Resolution', 300);
        fprintf('  %s saved.\n', cfg.fname);
    end

    % ---------- Cross-DIVI summary with SEM (raw + z-scored) ----------
    % For each mouse and each macro, pool the eroded-pixel-weighted mean
    % across STRU leaves of that macro -> per-mouse per-macro values,
    % then mean +/- SEM across mice.
    macro_pm_mean   = nan(n_mice_pm, n_macros);
    macro_pm_mean_z = nan(n_mice_pm, n_macros);
    for mi = 1:n_macros
        mkey = macro_names{mi};
        idx  = macro_groups.(mkey);
        wts  = roi_px_ero(idx);  wts = wts(:);   % force column to match vals_k
        for k_mouse = 1:n_mice_pm
            vals_k   = per_mouse_mean(k_mouse, idx).';    % .' makes column to match wts
            vals_k_z = per_mouse_mean_z(k_mouse, idx).';
            keep_kk  = ~isnan(vals_k) & wts > 0;
            if any(keep_kk)
                macro_pm_mean(k_mouse, mi)   = sum(vals_k(keep_kk)   .* wts(keep_kk)) / sum(wts(keep_kk));
                macro_pm_mean_z(k_mouse, mi) = sum(vals_k_z(keep_kk) .* wts(keep_kk)) / sum(wts(keep_kk));
            end
        end
    end
    macro_n_pm    = sum(~isnan(macro_pm_mean), 1);
    macro_mu_pm   = mean(macro_pm_mean,   1, 'omitnan');
    macro_sem_pm  = std(macro_pm_mean,  0, 1, 'omitnan') ./ sqrt(max(macro_n_pm, 1));
    macro_sem_pm(macro_n_pm < 2) = NaN;
    macro_mu_pm_z  = mean(macro_pm_mean_z,   1, 'omitnan');
    macro_sem_pm_z = std(macro_pm_mean_z, 0, 1, 'omitnan') ./ sqrt(max(macro_n_pm, 1));
    macro_sem_pm_z(macro_n_pm < 2) = NaN;

    cd_pm_configs = struct( ...
        'mu',     {macro_mu_pm(:),       macro_mu_pm_z(:)}, ...
        'sem',    {macro_sem_pm(:),      macro_sem_pm_z(:)}, ...
        'pm_mat', {macro_pm_mean,        macro_pm_mean_z}, ...
        'thr',    {raw_threshold,        zscore_threshold}, ...
        'flag',   {show_threshold_raw,   show_threshold_zscore}, ...
        'xlabel', {['Mean LR-sum ' channel ' (mean +/- SEM)'], 'Z-scored mean LR-sum (mean +/- SEM)'}, ...
        'title',  {['DIVI-level mean ' channel ' sum (mean +/- SEM, n=' num2str(n_mice_pm) ' mice) - ' strrep(merged_tag,'_',' ')], ...
                   ['DIVI-level z-scored ' channel ' sum (mean +/- SEM, n=' num2str(n_mice_pm) ' mice) - ' strrep(merged_tag,'_',' ')]}, ...
        'suffix', {['Region_MeanSum_BarAcrossDivi_' merged_tag smooth_suffix '_withSEM'], ...
                   ['Region_Zscore_BarAcrossDivi_' merged_tag smooth_suffix '_withSEM']});

    for ci_cd = 1:numel(cd_pm_configs)
        cfg = cd_pm_configs(ci_cd);
        [s_vals, s_ord] = sort(cfg.mu, 'ascend');
        s_sem    = cfg.sem(s_ord);
        s_labels = macro_display(s_ord);
        s_rel    = macro_mean_rel(s_ord);
        s_pm     = cfg.pm_mat(:, s_ord);
        keep_m = ~isnan(s_vals);
        s_vals   = s_vals(keep_m);
        s_sem    = s_sem(keep_m);
        s_labels = s_labels(keep_m);
        s_rel    = s_rel(keep_m);
        s_pm     = s_pm(:, keep_m);

        fig_cd = figure('Visible', 'off', 'Color', 'w', 'Units', 'Normalized', 'Position', [0 0 0.45 0.6]);
        if ~isempty(s_vals)
            b = barh(s_vals);
            b.FaceColor = 'flat';
            t_c = max(0, min(s_rel, t_max_clamp));
            c_i = round((t_c / t_max_clamp) * (size(c_map_rel, 1) - 1)) + 1;
            c_i(isnan(c_i)) = 1;
            for k = 1:numel(s_vals)
                b.CData(k, :) = c_map_rel(c_i(k), :);
            end
            yticks(1:numel(s_labels));
            set(gca, 'FontSize', 9);
            hold on;

            errorbar(s_vals, 1:numel(s_vals), s_sem, 'horizontal', ...
                'Color', 'k', 'LineStyle', 'none', 'LineWidth', 0.7, 'CapSize', 5);

            if show_per_mouse_dots
                n_pm_loc = size(s_pm, 1);
                for k = 1:numel(s_vals)
                    yk = k + (rand(1, n_pm_loc) - 0.5) * 0.3;
                    scatter(s_pm(:, k), yk, 14, [0.3 0.3 0.3], 'filled', ...
                        'MarkerFaceAlpha', 0.5);
                end
            end

            if cfg.flag
                above = s_vals >= cfg.thr;
                highlight_enriched_labels(gca, s_labels, above, 9);
                xline(cfg.thr, '--', sprintf('thr = %g', cfg.thr), ...
                    'Color', [0.85 0.5 0.0], 'LineWidth', 1.5, ...
                    'LabelOrientation', 'horizontal', 'FontSize', 8);
            else
                yticklabels(s_labels);
            end
            xlabel(cfg.xlabel);
            title(cfg.title, 'Interpreter', 'none');
            grid on;
            ylim([0 numel(s_labels) + 1]);
        end
        set(fig_cd, 'InvertHardcopy', 'off');
        saveas(fig_cd, fullfile(out_dir, [cfg.suffix '.fig']));
        exportgraphics(fig_cd, fullfile(out_dir, [cfg.suffix '.png']), 'Resolution', 300);
        fprintf('  %s saved.\n', cfg.suffix);
    end
end

%% Threshold-highlighted volume videos
% Single-panel videos with percentile clim, cyan threshold contour, and
% above-threshold region labels highlighted in cyan. Placed after barplot
% analysis so the enriched region list is available.

if generate_zscore_video
    merged_tag_display = strrep(merged_tag, '_', ' ');

    % Build per-region enrichment masks (if barplot analysis ran)
    if exist('roi_mean', 'var') && exist('roi_mean_z', 'var')
        enriched_raw = roi_mean   >= raw_threshold;
        enriched_z   = roi_mean_z >= zscore_threshold;
        % Map to label_acronyms order (they share roi_acronyms)
        label_enriched_raw = false(numel(label_acronyms), 1);
        label_enriched_z   = false(numel(label_acronyms), 1);
        [~, loc] = ismember(label_acronyms, roi_acronyms);
        valid = loc > 0;
        label_enriched_raw(valid) = enriched_raw(loc(valid));
        label_enriched_z(valid)   = enriched_z(loc(valid));
    else
        label_enriched_raw = false(numel(label_acronyms), 1);
        label_enriched_z   = false(numel(label_acronyms), 1);
    end

    % Adaptive clim for raw volume (P9 style)
    data_bv = mean_lr_sum(brainMask_merged & ~isnan(mean_lr_sum));
    raw_clim = [prctile(data_bv, 5), prctile(data_bv, 90)];
    fprintf('Raw video clim (p5-p90): [%.2f, %.2f]\n', raw_clim(1), raw_clim(2));

    % Adaptive clim for z-scored volume
    z_bv = zscore_lr_sum(brainMask_analysis & ~isnan(zscore_lr_sum));
    z_clim = [prctile(z_bv, 5), prctile(z_bv, 90)];
    fprintf('Z-score video clim (p5-p90): [%.2f, %.2f]\n', z_clim(1), z_clim(2));

    % Raw intensity threshold video
    fprintf('Writing raw intensity threshold video...\n');
    write_volume_threshold_video( ...
        mean_lr_sum, mean_lr_sum, half_atlas, brainMask_merged, ...
        raw_threshold, show_threshold_raw_video, raw_clim, ...
        out_dir, ['lr_sum_threshold_' merged_tag smooth_suffix '.mp4'], ...
        merged_tag_display, ['Mean LR-sum (' channel ')'], sprintf('thr >= %g', raw_threshold), ...
        label_centroids, label_acronyms, label_fontsize, label_color, label_enriched_raw);

    % Z-scored threshold video
    fprintf('Writing z-scored threshold video...\n');
    write_volume_threshold_video( ...
        zscore_lr_sum, zscore_lr_sum, half_atlas, brainMask_merged, ...
        zscore_threshold, show_threshold_zscore_video, z_clim, ...
        out_dir, ['lr_sum_zscore_' merged_tag smooth_suffix '.mp4'], ...
        merged_tag_display, ['Z-scored ' channel], sprintf('z >= %g', zscore_threshold), ...
        label_centroids, label_acronyms, label_fontsize, label_color, label_enriched_z);
end

fprintf('P8 done. Outputs in: %s\n', out_dir);

%% Diagnostic: visualize distance weights on a specific slice

diag_slice = 620;  % the slice where HP bleedthrough is visible
diag_power = [1, 2, 4];  % linear, quadratic, quartic weighting

fprintf('Generating distance-weight diagnostic for slice %d...\n', diag_slice);

% Build a composite weight map for this slice: each ROI's distance
% transform painted into a single 2D image. Also build the nano + weighted
% nano side by side.
[~, ~, nw] = size(AllenCrop);
hw = floor(nw / 2);

% Get atlas boundaries for overlay
atl_sl = single(squeeze(AllenCrop(diag_slice, :, 1:hw)));
bnd = gradient(atl_sl) ~= 0 & atl_sl > 0;
[br, bc] = find(bnd);

% Get nano slice
data_sl = squeeze(mean_lr_sum(diag_slice, :, :));

for pi = 1:numel(diag_power)
    pw = diag_power(pi);

    % Composite distance weight map for this slice. Uses the atlas
    % annotation directly: each unique region ID in the slice gets its
    % own distance transform. No get_allen_region_mask calls needed.
    atlas_sl = squeeze(AllenCrop(diag_slice, :, 1:hw));
    region_ids = unique(atlas_sl(atlas_sl > 0));
    weight_map = zeros(size(data_sl), 'single');
    for ri = 1:numel(region_ids)
        sl_mask = atlas_sl == region_ids(ri);
        if nnz(sl_mask) < 5, continue; end
        dist_sl = single(bwdist(~sl_mask));
        dist_sl(~sl_mask) = 0;
        mx = max(dist_sl(:));
        if mx > 0, dist_sl = (dist_sl / mx) .^ pw; end
        weight_map = max(weight_map, dist_sl);
    end

    % Weighted nano
    data_weighted = data_sl .* weight_map;

    fh_diag = figure('Visible', 'off', 'Color', 'k', 'Units', 'Normalized', 'Position', [0 0 1 0.9]);

    % Panel 1: raw nano
    subplot(1, 4, 1);
    imagesc(data_sl); clim([0 5]); colormap(gca, hot);
    set(findobj(gca,'Type','image'), 'AlphaData', atl_sl > 0);
    axis image off; hold on;
    plot(bc, br, '.', 'Color', [0.66 0.66 0.66], 'MarkerSize', 0.5);
    title(['Raw ' channel], 'Color', 'w', 'FontSize', 12);
    cb = colorbar; cb.Color = 'w';
    set(gca, 'Color', 'k');

    % Panel 2: distance weights
    subplot(1, 4, 2);
    imagesc(weight_map); clim([0 1]); colormap(gca, hot);
    set(findobj(gca,'Type','image'), 'AlphaData', atl_sl > 0);
    axis image off; hold on;
    plot(bc, br, '.', 'Color', [0.66 0.66 0.66], 'MarkerSize', 0.5);
    title(sprintf('Distance weights (power=%d)', pw), 'Color', 'w', 'FontSize', 12);
    cb = colorbar; cb.Color = 'w';
    set(gca, 'Color', 'k');

    % Panel 3: weighted nano
    subplot(1, 4, 3);
    imagesc(data_weighted); clim([0 5]); colormap(gca, hot);
    set(findobj(gca,'Type','image'), 'AlphaData', atl_sl > 0);
    axis image off; hold on;
    plot(bc, br, '.', 'Color', [0.66 0.66 0.66], 'MarkerSize', 0.5);
    title(sprintf('Nano x weights^%d', pw), 'Color', 'w', 'FontSize', 12);
    cb = colorbar; cb.Color = 'w';
    set(gca, 'Color', 'k');

    % Panel 4: difference (raw - weighted) shows what gets suppressed
    subplot(1, 4, 4);
    diff_diag = data_sl - data_weighted;
    imagesc(diff_diag); clim([0 2]); colormap(gca, hot);
    set(findobj(gca,'Type','image'), 'AlphaData', atl_sl > 0);
    axis image off; hold on;
    plot(bc, br, '.', 'Color', [0.66 0.66 0.66], 'MarkerSize', 0.5);
    title('Suppressed signal (raw - weighted)', 'Color', 'w', 'FontSize', 12);
    cb = colorbar; cb.Color = 'w';
    set(gca, 'Color', 'k');

    sgtitle(sprintf('Slice %d — Distance weighting diagnostic (power = %d)', ...
        diag_slice, pw), 'Color', 'w', 'FontSize', 14);
    set(fh_diag, 'InvertHardcopy', 'off');
    exportgraphics(fh_diag, fullfile(out_dir, ...
        sprintf('Diagnostic_DistWeight_slice%d_pow%d.png', diag_slice, pw)), ...
        'Resolution', 200, 'BackgroundColor', 'k');
    close(fh_diag);
end
fprintf('  Distance weight diagnostics saved.\n');

%% Local functions

function write_single_mouse_sum_video(lr_sum_vol, brain_mask, atlas_vol, ...
    save_dir, video_filename, clim_values, mouse_name)
% Single-panel full-screen per-mouse LR-sum video. Styled to match
% write_lr_video (hot colormap, atlas boundaries overlay, half-width xlim)
% but renders only one axes so each mouse gets the full frame.

    if ~exist(save_dir, 'dir'), mkdir(save_dir); end
    full_video_path = fullfile(save_dir, video_filename);
    vidObj = VideoWriter(full_video_path, 'MPEG-4');
    vidObj.FrameRate = 15;
    vidObj.Quality   = 95;
    open(vidObj);

    n_slices = size(lr_sum_vol, 1);
    mouse_name_disp = strrep(mouse_name, '_', ' ');

    fprintf('Writing video: %s\n', video_filename);

    for j = 1:n_slices
        if sum(sum(brain_mask(j, :, :))) == 0
            continue;
        end

        atlasim = squeeze(atlas_vol(j, :, :));
        atlasim = single(atlasim);
        av_warp_boundaries = gradient(atlasim) ~= 0 & (atlasim > 1);
        [row, col] = ind2sub(size(atlasim), find(av_warp_boundaries));

        fh = figure('visible', 'off', 'units', 'normalized', ...
            'outerposition', [0 0 1 1], 'Color', 'k');
        set(fh, 'InvertHardcopy', 'off');

        h1 = imagesc(squeeze(lr_sum_vol(j, :, :)));
        clim(clim_values);
        set(h1, 'AlphaData', squeeze(brain_mask(j, :, 1:size(lr_sum_vol, 3))));
        colormap(gca, hot);

        ax1 = gca; ax1.Color = 'k';
        axis equal; axis off;
        hold on;
        line(col, row, 'Marker', '.', 'LineStyle', 'none', ...
            'Color', [0.66 0.66 0.66], 'MarkerSize', 0.5);
        xlim([0, size(lr_sum_vol, 3)]);

        cb1 = colorbar;
        cb1.Label.String   = ['LR sum (' channel ')'];
        cb1.Label.FontSize = 10;
        cb1.Color          = 'w';
        cb1.Label.Color    = 'w';

        sgtitle([mouse_name_disp [' - LR sum (' channel ') - Slice # '] num2str(j)], ...
            'Color', 'w', 'FontSize', 14);

        frame = getframe(fh);
        writeVideo(vidObj, frame);
        close(fh);

        if mod(j, 50) == 0
            fprintf('  Frame %d written...\n', j);
        end
    end

    close(vidObj);
    fprintf('Video saved: %s\n', full_video_path);
end


function write_lr_video_labeled(lr_diff_vol, lr_sum_vol, atlas_vol, brain_mask, ...
    save_dir, video_filename, clim_values, group_name, ...
    label_string_diff, label_string_sum, ...
    label_centroids, label_acronyms, label_fontsize, label_color)
% Two-panel LR video (same layout as write_lr_video) with region-acronym
% text overlays. label_centroids is [n_slices x n_rois x 2] where
% (:,:,1) = row and (:,:,2) = col; NaN entries are skipped per slice.

    if ~exist(save_dir, 'dir'), mkdir(save_dir); end
    full_video_path = fullfile(save_dir, video_filename);
    vidObj = VideoWriter(full_video_path, 'MPEG-4');
    vidObj.FrameRate = 15;
    vidObj.Quality   = 95;
    open(vidObj);
    n_slices = size(lr_diff_vol, 1);
    n_rois_lbl = numel(label_acronyms);
    have_labels = ~isempty(label_centroids) && n_rois_lbl > 0;

    fprintf('Writing labeled video: %s\n', video_filename);

    for j = 1:n_slices
        if sum(sum(brain_mask(j, :, :))) == 0
            continue;
        end

        atlasim = squeeze(atlas_vol(j, :, :));
        atlasim = single(atlasim);
        av_warp_boundaries = gradient(atlasim) ~= 0 & (atlasim > 1);
        [row, col] = ind2sub(size(atlasim), find(av_warp_boundaries));

        fh = figure('visible', 'off', 'units', 'normalized', ...
            'outerposition', [0 0 1 1], 'Color', 'k');
        set(fh, 'InvertHardcopy', 'off');

        % ---- Panel 1: diff (left) ----
        subplot(1, 2, 1);
        h1 = imagesc(squeeze(lr_diff_vol(j, :, :)));
        clim(clim_values);
        set(h1, 'AlphaData', squeeze(brain_mask(j, :, 1:size(lr_diff_vol, 3))));
        if abs(clim_values(1)) == abs(clim_values(2))
            try
                colormap(gca, get_color2color_colormap([0, 0, 1], [1, 0, 0]));
            catch
                colormap(gca, jet);
            end
        else
            colormap(gca, hot);
        end
        ax1 = gca; ax1.Color = 'k';
        axis equal; axis off;
        hold on;
        line(col, row, 'Marker', '.', 'LineStyle', 'none', ...
            'Color', [0.66 0.66 0.66], 'MarkerSize', 0.5);
        xlim([0, size(lr_diff_vol, 3)]);
        title([group_name ' - ' label_string_diff], 'Color', 'w', 'FontSize', 12);
        cb1 = colorbar;
        cb1.Label.String   = label_string_diff;
        cb1.Label.FontSize = 10;
        cb1.Color          = 'w';
        cb1.Label.Color    = 'w';

        if have_labels
            for r = 1:n_rois_lbl
                cy = label_centroids(j, r, 1);
                cx = label_centroids(j, r, 2);
                if isnan(cy) || isnan(cx), continue; end
                text(cx, cy, label_acronyms{r}, ...
                    'Color', label_color, ...
                    'FontSize', label_fontsize, ...
                    'FontWeight', 'bold', ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'middle', ...
                    'Interpreter', 'none', ...
                    'Clipping', 'on');
            end
        end

        % ---- Panel 2: sum (right) ----
        subplot(1, 2, 2);
        h2 = imagesc(squeeze(lr_sum_vol(j, :, :)));
        clim(2 * clim_values);
        set(h2, 'AlphaData', squeeze(brain_mask(j, :, 1:size(lr_sum_vol, 3))));
        if abs(clim_values(1)) == abs(clim_values(2))
            try
                colormap(gca, get_color2color_colormap([0, 0, 1], [1, 0, 0]));
            catch
                colormap(gca, jet);
            end
        else
            colormap(gca, hot);
        end
        ax2 = gca; ax2.Color = 'k';
        axis equal; axis off;
        hold on;
        line(col, row, 'Marker', '.', 'LineStyle', 'none', ...
            'Color', [0.66 0.66 0.66], 'MarkerSize', 0.5);
        xlim([0, size(lr_sum_vol, 3)]);
        title([group_name ' - ' label_string_sum], 'Color', 'w', 'FontSize', 12);
        cb2 = colorbar;
        cb2.Label.String   = label_string_sum;
        cb2.Label.FontSize = 10;
        cb2.Color          = 'w';
        cb2.Label.Color    = 'w';

        if have_labels
            for r = 1:n_rois_lbl
                cy = label_centroids(j, r, 1);
                cx = label_centroids(j, r, 2);
                if isnan(cy) || isnan(cx), continue; end
                text(cx, cy, label_acronyms{r}, ...
                    'Color', label_color, ...
                    'FontSize', label_fontsize, ...
                    'FontWeight', 'bold', ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'middle', ...
                    'Interpreter', 'none', ...
                    'Clipping', 'on');
            end
        end

        sgtitle(['Slice # ' num2str(j)], 'Color', 'w', 'FontSize', 14);

        frame = getframe(fh);
        writeVideo(vidObj, frame);
        close(fh);

        if mod(j, 50) == 0
            fprintf('  Frame %d written...\n', j);
        end
    end

    close(vidObj);
    fprintf('Video saved: %s\n', full_video_path);
end


function write_single_mouse_sum_video_labeled(lr_sum_vol, brain_mask, atlas_vol, ...
    save_dir, video_filename, clim_values, mouse_name, ...
    label_centroids, label_acronyms, label_fontsize, label_color)
% Single-panel per-mouse LR-sum video with region-acronym text overlay.
% Same rendering as write_single_mouse_sum_video plus a text pass per slice.

    if ~exist(save_dir, 'dir'), mkdir(save_dir); end
    full_video_path = fullfile(save_dir, video_filename);
    vidObj = VideoWriter(full_video_path, 'MPEG-4');
    vidObj.FrameRate = 15;
    vidObj.Quality   = 95;
    open(vidObj);

    n_slices = size(lr_sum_vol, 1);
    mouse_name_disp = strrep(mouse_name, '_', ' ');
    n_rois_lbl = numel(label_acronyms);
    have_labels = ~isempty(label_centroids) && n_rois_lbl > 0;

    fprintf('Writing labeled video: %s\n', video_filename);

    for j = 1:n_slices
        if sum(sum(brain_mask(j, :, :))) == 0
            continue;
        end

        atlasim = squeeze(atlas_vol(j, :, :));
        atlasim = single(atlasim);
        av_warp_boundaries = gradient(atlasim) ~= 0 & (atlasim > 1);
        [row, col] = ind2sub(size(atlasim), find(av_warp_boundaries));

        fh = figure('visible', 'off', 'units', 'normalized', ...
            'outerposition', [0 0 1 1], 'Color', 'k');
        set(fh, 'InvertHardcopy', 'off');

        h1 = imagesc(squeeze(lr_sum_vol(j, :, :)));
        clim(clim_values);
        set(h1, 'AlphaData', squeeze(brain_mask(j, :, 1:size(lr_sum_vol, 3))));
        colormap(gca, hot);

        ax1 = gca; ax1.Color = 'k';
        axis equal; axis off;
        hold on;
        line(col, row, 'Marker', '.', 'LineStyle', 'none', ...
            'Color', [0.66 0.66 0.66], 'MarkerSize', 0.5);
        xlim([0, size(lr_sum_vol, 3)]);

        if have_labels
            for r = 1:n_rois_lbl
                cy = label_centroids(j, r, 1);
                cx = label_centroids(j, r, 2);
                if isnan(cy) || isnan(cx), continue; end
                text(cx, cy, label_acronyms{r}, ...
                    'Color', label_color, ...
                    'FontSize', label_fontsize, ...
                    'FontWeight', 'bold', ...
                    'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'middle', ...
                    'Interpreter', 'none', ...
                    'Clipping', 'on');
            end
        end

        cb1 = colorbar;
        cb1.Label.String   = ['LR sum (' channel ')'];
        cb1.Label.FontSize = 10;
        cb1.Color          = 'w';
        cb1.Label.Color    = 'w';

        sgtitle([mouse_name_disp [' - LR sum (' channel ') - Slice # '] num2str(j)], ...
            'Color', 'w', 'FontSize', 14);

        frame = getframe(fh);
        writeVideo(vidObj, frame);
        close(fh);

        if mod(j, 50) == 0
            fprintf('  Frame %d written...\n', j);
        end
    end

    close(vidObj);
    fprintf('Video saved: %s\n', full_video_path);
end


function [roi_names, roi_acronyms, roi_macro] = build_stru_leaves_by_macro(allenDir, macro_divi_list)
% Walk the Allen CCF ontology and collect, for each DIVI-level macro in
% macro_divi_list, all its STRU-leaf descendants. A STRU leaf is a STRU
% node none of whose descendants (anywhere in the tree) are also STRU —
% i.e. its children are either SUBS (cortical layers, sub-layers) or
% nothing.
%
% Returns three parallel cell arrays, one entry per collected leaf:
%   roi_names    - full region name (e.g. "Primary somatosensory area, barrel field")
%   roi_acronyms - short acronym   (e.g. "SSp-bfd")
%   roi_macro    - parent DIVI acronym (e.g. "Isocortex")

    parc_tbl = readtable(fullfile(allenDir, 'parcellation_term.csv'), 'TextType', 'string');
    mem_tbl  = readtable(fullfile(allenDir, 'parcellation_term_set_membership.csv'), 'TextType', 'string');

    % Maps: identifier (MBA:xxx) -> row info
    ids       = cellstr(parc_tbl.identifier);
    names     = cellstr(parc_tbl.name);
    acronyms  = cellstr(parc_tbl.acronym);
    parents   = cellstr(parc_tbl.parent_identifier);
    id_to_idx = containers.Map(ids, 1:numel(ids));

    % Level membership (label -> set acronym). Labels in membership table
    % are "AllenCCF-Ontology-2017-<num>" and match parcellation_term.label.
    mem_labels = cellstr(mem_tbl.parcellation_term_label);
    mem_sets   = cellstr(mem_tbl.parcellation_term_set_label);
    % Strip "AllenCCF-Ontology-2017-" prefix to get the set suffix (DIVI/STRU/SUBS)
    set_suffix = regexprep(mem_sets, '^AllenCCF-Ontology-2017-', '');
    label_to_set = containers.Map(mem_labels, set_suffix);

    % Walk children: build parent -> children list (by identifier MBA:xxx)
    children_map = containers.Map('KeyType', 'char', 'ValueType', 'any');
    for k = 1:numel(ids)
        p = parents{k};
        if isempty(p), continue; end
        if isKey(children_map, p)
            children_map(p) = [children_map(p), k];
        else
            children_map(p) = k;
        end
    end

    % Helper: lookup set level for a given row index, via parcellation_term.label
    row_set = repmat({'NONE'}, numel(ids), 1);
    term_labels = cellstr(parc_tbl.label);
    for k = 1:numel(ids)
        if isKey(label_to_set, term_labels{k})
            row_set{k} = label_to_set(term_labels{k});
        end
    end

    % Acronym -> index (used to find macro root nodes)
    acronym_to_idx = containers.Map(acronyms, 1:numel(acronyms));

    roi_names    = {};
    roi_acronyms = {};
    roi_macro    = {};

    for mi = 1:numel(macro_divi_list)
        macro_ac = macro_divi_list{mi};
        if ~isKey(acronym_to_idx, macro_ac)
            warning('Macro acronym not found in ontology: %s', macro_ac);
            continue;
        end
        root_k = acronym_to_idx(macro_ac);
        if ~strcmp(row_set{root_k}, 'DIVI')
            warning('Macro "%s" is not a DIVI-level node (found: %s)', macro_ac, row_set{root_k});
        end

        % BFS from root through children_map to collect all descendants
        stack = root_k;
        all_desc = [];
        while ~isempty(stack)
            k = stack(end); stack(end) = [];
            id_k = ids{k};
            if isKey(children_map, id_k)
                ch = children_map(id_k);
                all_desc = [all_desc, ch]; %#ok<AGROW>
                stack = [stack, ch]; %#ok<AGROW>
            end
        end
        all_desc = unique(all_desc);

        % Keep only STRU descendants
        is_stru = cellfun(@(s) strcmp(s, 'STRU'), row_set(all_desc));
        stru_desc = all_desc(is_stru);

        % Filter to leaves: a STRU node is a leaf iff none of its descendants
        % are also STRU (cortical umbrella SSp has STRU children, so excluded;
        % SSp-bfd has only SUBS layer children, so included).
        is_leaf = false(size(stru_desc));
        for li = 1:numel(stru_desc)
            kk = stru_desc(li);
            id_k = ids{kk};
            sub_stack = [];
            if isKey(children_map, id_k)
                sub_stack = children_map(id_k);
            end
            has_stru_child = false;
            while ~isempty(sub_stack)
                ss = sub_stack(end); sub_stack(end) = [];
                if strcmp(row_set{ss}, 'STRU')
                    has_stru_child = true;
                    break;
                end
                id_ss = ids{ss};
                if isKey(children_map, id_ss)
                    sub_stack = [sub_stack, children_map(id_ss)]; %#ok<AGROW>
                end
            end
            is_leaf(li) = ~has_stru_child;
        end
        leaves = stru_desc(is_leaf);

        for li = 1:numel(leaves)
            kk = leaves(li);
            roi_names{end+1}    = names{kk};    %#ok<AGROW>
            roi_acronyms{end+1} = acronyms{kk}; %#ok<AGROW>
            roi_macro{end+1}    = macro_ac;     %#ok<AGROW>
        end
    end

    roi_names    = roi_names(:);
    roi_acronyms = roi_acronyms(:);
    roi_macro    = roi_macro(:);
end


function highlight_enriched_labels(ax_handle, acro_list, above_mask, fontsize)
% Shared helper: highlight above-threshold region labels in orange bold.
% Uses the TickLabel ruler API (R2021b+), falls back to manual text.
    enriched_color = [0.85 0.5 0.0];  % orange
    yticklabels(ax_handle, acro_list);
    drawnow;
    try
        tl = ax_handle.YAxis.TickLabels;
        for k = 1:numel(acro_list)
            if above_mask(k)
                tl(k).FontColor  = enriched_color;
                tl(k).FontWeight = 'bold';
            end
        end
    catch
        yticklabels(ax_handle, {});
        xl = xlim(ax_handle);
        for k = 1:numel(acro_list)
            clr = [0 0 0]; fw = 'normal';
            if above_mask(k), clr = enriched_color; fw = 'bold'; end
            text(ax_handle, xl(1) - 0.02*(xl(2)-xl(1)), k, acro_list{k}, ...
                'Color', clr, 'FontSize', fontsize, 'FontWeight', fw, ...
                'HorizontalAlignment', 'right', 'Units', 'data', ...
                'Interpreter', 'none', 'Clipping', 'off');
        end
    end
end


function vol_z = zscore_in_mask(vol, mask)
% Z-score a volume using mean and std computed within mask voxels.
% Voxels outside mask are set to NaN.
    vals = vol(mask); vals = vals(~isnan(vals));
    mu = mean(vals); sg = std(vals);
    if sg==0 || isnan(sg), vol_z = nan(size(vol), 'like', vol); return; end
    vol_z = (vol - mu) / sg;
    vol_z(~mask) = NaN;
end


function write_volume_threshold_video(display_vol, threshold_vol, atlas_vol, brain_mask, ...
    threshold, show_threshold, clim_vals, ...
    save_dir, video_filename, ...
    group_name, panel_label, threshold_label, ...
    label_centroids, label_acronyms, label_fontsize, label_color, ...
    label_enriched)
% Single-panel volume video with optional threshold contour overlay.
%   display_vol    — volume to render (e.g. raw mean_lr_sum)
%   threshold_vol  — volume to threshold for contour
%   show_threshold — true = cyan contour around above-threshold voxels
%   label_enriched — logical vector (n_rois x 1): true = show label in cyan
%
% Uses percentile-based clim (computed by caller, matching P9 style).

    if ~exist(save_dir, 'dir'), mkdir(save_dir); end
    full_video_path = fullfile(save_dir, video_filename);
    vidObj = VideoWriter(full_video_path, 'MPEG-4');
    vidObj.FrameRate = 15;
    vidObj.Quality   = 95;
    open(vidObj);

    n_slices   = size(display_vol, 1);
    n_rois_lbl = numel(label_acronyms);
    have_labels = ~isempty(label_centroids) && n_rois_lbl > 0;
    if nargin < 17 || isempty(label_enriched)
        label_enriched = false(n_rois_lbl, 1);
    end
    enriched_color = [0 1 1];  % cyan for above-threshold labels

    fprintf('Writing video: %s\n', video_filename);

    for j = 1:n_slices
        bm_sl = squeeze(brain_mask(j, :, :));
        if nnz(bm_sl) == 0, continue; end

        d_sl = squeeze(display_vol(j, :, :));

        % Atlas boundaries
        atl_sl = single(squeeze(atlas_vol(j, :, 1:size(d_sl, 2))));
        bnd = gradient(atl_sl) ~= 0 & atl_sl > 0;
        [br, bc] = find(bnd);

        fh = figure('visible', 'off', 'units', 'normalized', ...
            'outerposition', [0 0 1 1], 'Color', 'k');
        set(fh, 'InvertHardcopy', 'off');

        h1 = imagesc(d_sl); clim(clim_vals); colormap(gca, hot);
        set(h1, 'AlphaData', bm_sl);
        axis image off; hold on;
        plot(bc, br, '.', 'Color', [0.5 0.5 0.5], 'MarkerSize', 0.5);
        xlim([0 size(d_sl, 2)]);

        % Threshold contour overlay
        if show_threshold
            t_sl = squeeze(threshold_vol(j, :, :));
            above_mask  = t_sl >= threshold & bm_sl;
            above_perim = bwperim(above_mask);
            [ar, ac_p]  = find(above_perim);
            if ~isempty(ar)
                plot(ac_p, ar, '.', 'Color', [0 1 1], 'MarkerSize', 1.5);
            end
        end

        % Title: all info in one line, no sgtitle
        if show_threshold
            title(sprintf('%s - %s - Slice #%d - enriched (%s)', ...
                group_name, panel_label, j, threshold_label), ...
                'Color', 'w', 'FontSize', 11);
        else
            title(sprintf('%s - %s - Slice #%d', group_name, panel_label, j), ...
                'Color', 'w', 'FontSize', 11);
        end
        cb1 = colorbar; cb1.Color = 'w';
        cb1.Label.String = panel_label; cb1.Label.FontSize = 10; cb1.Label.Color = 'w';
        set(gca, 'Color', 'k');

        % Region labels: enriched in cyan, others in default color
        if have_labels
            for r = 1:n_rois_lbl
                cy = label_centroids(j, r, 1);
                cx = label_centroids(j, r, 2);
                if isnan(cy) || isnan(cx), continue; end
                if label_enriched(r)
                    lc = enriched_color;
                else
                    lc = label_color;
                end
                text(cx, cy, label_acronyms{r}, ...
                    'Color', lc, 'FontSize', label_fontsize, ...
                    'FontWeight', 'bold', 'HorizontalAlignment', 'center', ...
                    'VerticalAlignment', 'middle', 'Interpreter', 'none', ...
                    'Clipping', 'on');
            end
        end

        frame = getframe(fh);
        writeVideo(vidObj, frame);
        close(fh);

        if mod(j, 50) == 0
            fprintf('  Frame %d written...\n', j);
        end
    end

    close(vidObj);
    fprintf('Video saved: %s\n', full_video_path);
end
