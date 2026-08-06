clear all
close all
clc

% /// Pipeline script #10: side-by-side nano vs autofluorescence comparison ///
% Loads the per-mouse-per-region matrices that P8 cached for each channel and:
%   (1) Plots grouped horizontal bars (nano + auto) per region, by DIVI macro
%   (2) Overlays per-mouse dots and connects paired (same-mouse) dots
%       with thin gray lines so the paired structure is visible
%   (3) Runs paired Wilcoxon signed-rank test (signrank) per region
%   (4) Optionally applies Bonferroni correction across regions tested
%   (5) Highlights region labels in bright orange when the test passes
%   (6) Saves a CSV with per-region p-values for inspection

%% User-defined parameters

% Paths to the two channel folders (must contain per_mouse_region_means_*.mat)
nano_dir = 'D:\sep_histology\data\comparisons\merged_naive_rws_nano';
auto_dir = 'D:\sep_histology\data\comparisons\merged_naive_rws_auto';
nano_tag = 'merged_naive_rws_nano';
auto_tag = 'merged_naive_rws_auto';
smooth_suffix = '_nosmooth';

% Output directory
out_dir = 'D:\sep_histology\data\comparisons\nano_vs_auto';
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

% Aggregation method (must match what was cached: 'distweight' or 'eroded')
agg_method = 'distweight';

% Statistics
apply_bonferroni = false;    % Bonferroni correction across all regions tested
alpha            = 0.05;     % significance threshold (per-region or family-wise)

% Display options
show_per_mouse_dots = true;
show_paired_lines   = true;
sort_regions_by     = 'nano';   % 'nano' | 'auto' | 'diff' (nano - auto)
% Macros (DIVI list) — drives panel grid, must match P8's order.
macro_divi_list = {'Isocortex','OLF','HPF','CTXsp','STR','PAL','TH','HY','MB'};

% Which figures / outputs to produce (turn any off to skip)
produce_paired_per_region = true;   % Fig 1: paired nano+auto bars per STRU leaf
produce_deltaz_per_region = true;   % Fig 2: delta-z (nano_z - auto_z) bar per STRU leaf
produce_paired_per_macro  = true;   % Fig 3: paired nano+auto bars per DIVI macro
produce_deltaz_per_macro  = true;   % Fig 4: delta-z bar per DIVI macro
produce_contrast_video    = true;   % Voxelwise z-contrast slice-by-slice video

% Color palette
nano_color      = [0.95 0.55 0.10];   % strong orange (nano bar fill)
auto_color      = [0.95 0.85 0.20];   % yellow (auto bar fill)
nano_dot_color  = [0.65 0.30 0.00];   % darker orange (per-mouse dots)
auto_dot_color  = [0.70 0.60 0.00];   % darker yellow (per-mouse dots)
sig_label_color = [0.85 0.50 0.00];   % orange used in P8 enriched-label highlight
paired_line_color = [0.6 0.6 0.6];    % gray for paired-mouse connectors

%% Load per-mouse caches

nano_pm_path = fullfile(nano_dir, ['per_mouse_region_means_' nano_tag smooth_suffix '.mat']);
auto_pm_path = fullfile(auto_dir, ['per_mouse_region_means_' auto_tag smooth_suffix '.mat']);

if ~exist(nano_pm_path, 'file')
    error('Nano per-mouse cache not found:\n  %s\nRerun P8 with channel=''nano'' and compute_per_mouse_sem=true.', nano_pm_path);
end
if ~exist(auto_pm_path, 'file')
    error('Auto per-mouse cache not found:\n  %s\nRerun P8 with channel=''auto'' and compute_per_mouse_sem=true.', auto_pm_path);
end

fprintf('Loading nano cache: %s\n', nano_pm_path);
S_nano = load(nano_pm_path);
fprintf('Loading auto cache: %s\n', auto_pm_path);
S_auto = load(auto_pm_path);

% Pick per-mouse matrix per chosen aggregation
switch agg_method
    case 'distweight'
        nano_pm = S_nano.per_mouse_mean_dw;
        auto_pm = S_auto.per_mouse_mean_dw;
    case 'eroded'
        nano_pm = S_nano.per_mouse_mean_ero;
        auto_pm = S_auto.per_mouse_mean_ero;
    otherwise
        error('Unknown agg_method: %s (use ''distweight'' or ''eroded'')', agg_method);
end

% Sanity-check mice and regions match between channels
nano_mice = S_nano.per_mouse_mice;
auto_mice = S_auto.per_mouse_mice;
if ~isequal(nano_mice, auto_mice)
    error('Cohort mice mismatch between nano and auto caches. Cannot pair.');
end
n_mice = numel(nano_mice);

if ~isequal(S_nano.roi_acronyms, S_auto.roi_acronyms)
    error('ROI acronym lists differ between nano and auto caches.');
end
roi_acronyms = S_nano.roi_acronyms;
roi_macro    = S_nano.roi_macro;
n_rois       = numel(roi_acronyms);

fprintf('  %d mice paired across channels.\n', n_mice);
fprintf('  %d regions to compare.\n', n_rois);

%% Build macro -> region-index map (preserve macro_divi_list order)

macro_groups = struct();
for r = 1:n_rois
    mkey = matlab.lang.makeValidName(roi_macro{r});
    if ~isfield(macro_groups, mkey), macro_groups.(mkey) = []; end
    macro_groups.(mkey)(end+1) = r; %#ok<SAGROW>
end

%% Paired Wilcoxon per region, optional Bonferroni

% One-sided paired Wilcoxon: tests H1 = nano > auto (enrichment only).
% A region is "significant" when nano is reliably HIGHER than auto across
% the cohort. Regions where auto > nano stay non-significant regardless.
p_values = nan(n_rois, 1);
for r = 1:n_rois
    nv = nano_pm(:, r);
    av = auto_pm(:, r);
    valid = ~isnan(nv) & ~isnan(av);
    if sum(valid) >= 3   % minimum N for signrank to be meaningful
        try
            p_values(r) = signrank(nv(valid), av(valid), 'tail', 'right');
        catch
            p_values(r) = NaN;
        end
    end
end

n_tested = sum(~isnan(p_values));
if apply_bonferroni
    p_thresh = alpha / max(n_tested, 1);
    correction_label = sprintf('Bonferroni-corrected (alpha=%g, N=%d, threshold p<%.2e)', ...
        alpha, n_tested, p_thresh);
else
    p_thresh = alpha;
    correction_label = sprintf('Uncorrected (alpha=%g)', alpha);
end
sig_mask = p_values < p_thresh;

fprintf('  %s\n', correction_label);
fprintf('  Significant regions: %d / %d\n', sum(sig_mask), n_tested);

%% Save p-value table for inspection

ptbl = table( ...
    roi_acronyms(:), roi_macro(:), p_values(:), sig_mask(:), ...
    mean(nano_pm, 1, 'omitnan')', mean(auto_pm, 1, 'omitnan')', ...
    'VariableNames', {'acronym', 'macro', 'p_value', 'significant', 'nano_mean', 'auto_mean'});
ptbl = sortrows(ptbl, 'p_value', 'MissingPlacement', 'last');
csv_name = ['NanoVsAuto_paired_signrank' smooth_suffix];
if apply_bonferroni, csv_name = [csv_name '_bonferroni']; else, csv_name = [csv_name '_uncorrected']; end
writetable(ptbl, fullfile(out_dir, [csv_name '.csv']));
fprintf('Saved p-value table: %s\n', fullfile(out_dir, [csv_name '.csv']));

%% Load cohort caches (used for: brain stats for z-scoring + 3D zscore volumes for video)

nano_cohort_path = fullfile(nano_dir, ['mean_lr_sum_' nano_tag smooth_suffix '.mat']);
auto_cohort_path = fullfile(auto_dir, ['mean_lr_sum_' auto_tag smooth_suffix '.mat']);
fprintf('Loading nano cohort cache: %s\n', nano_cohort_path);
S_nano_c = load(nano_cohort_path);
fprintf('Loading auto cohort cache: %s\n', auto_cohort_path);
S_auto_c = load(auto_cohort_path);

% Brain stats per channel: from cohort-mean LR-sum within analysis-restricted brain
brain_n_vals = S_nano_c.mean_lr_sum(S_nano_c.brainMask_analysis);
brain_n_vals = brain_n_vals(~isnan(brain_n_vals));
mu_n = mean(brain_n_vals);  sg_n = std(brain_n_vals);
brain_a_vals = S_auto_c.mean_lr_sum(S_auto_c.brainMask_analysis);
brain_a_vals = brain_a_vals(~isnan(brain_a_vals));
mu_a = mean(brain_a_vals);  sg_a = std(brain_a_vals);
fprintf('  Nano brain stats: mu=%.4f, sigma=%.4f\n', mu_n, sg_n);
fprintf('  Auto brain stats: mu=%.4f, sigma=%.4f\n', mu_a, sg_a);

%% Per-mouse z-scored region values + delta-z contrast + Wilcoxon

% z-score is per-channel against that channel's own brain distribution, so
% the contrast (nano_z - auto_z) measures how much more nano is enriched
% above its own brain mean than auto is above its own. Distinct from the
% raw nano > auto test since z-scoring rescales channels independently.
nano_pm_z   = (nano_pm - mu_n) / sg_n;
auto_pm_z   = (auto_pm - mu_a) / sg_a;
contrast_pm = nano_pm_z - auto_pm_z;

p_values_dz = nan(n_rois, 1);
for r = 1:n_rois
    cv = contrast_pm(:, r);
    valid = ~isnan(cv);
    if sum(valid) >= 3
        try
            p_values_dz(r) = signrank(cv(valid), 0, 'tail', 'right');
        catch
        end
    end
end

if apply_bonferroni
    p_thresh_dz = alpha / max(sum(~isnan(p_values_dz)), 1);
else
    p_thresh_dz = alpha;
end
sig_mask_dz = p_values_dz < p_thresh_dz;
fprintf('  Delta-z significant regions: %d / %d\n', ...
    sum(sig_mask_dz), sum(~isnan(p_values_dz)));

%% Save delta-z CSV

ptbl_dz = table( ...
    roi_acronyms(:), roi_macro(:), p_values_dz(:), sig_mask_dz(:), ...
    mean(contrast_pm, 1, 'omitnan')', mean(nano_pm_z, 1, 'omitnan')', mean(auto_pm_z, 1, 'omitnan')', ...
    'VariableNames', {'acronym', 'macro', 'p_value_deltaz', 'significant_deltaz', ...
                      'contrast_mean', 'nano_z_mean', 'auto_z_mean'});
ptbl_dz = sortrows(ptbl_dz, 'p_value_deltaz', 'MissingPlacement', 'last');
csv_name_dz = ['NanoVsAuto_deltaz_signrank' smooth_suffix];
if apply_bonferroni, csv_name_dz = [csv_name_dz '_bonferroni']; else, csv_name_dz = [csv_name_dz '_uncorrected']; end
writetable(ptbl_dz, fullfile(out_dir, [csv_name_dz '.csv']));
fprintf('Saved delta-z table: %s\n', fullfile(out_dir, [csv_name_dz '.csv']));

%% Load ROI mask cache for region pixel counts (used to weight macro pooling)

masks_filename = 'roi_masks_r3_pw4_slices100-700.mat';   % default geometry tag from P8/P9
masks_path = fullfile(nano_dir, masks_filename);
if ~exist(masks_path, 'file')
    masks_path = fullfile(auto_dir, masks_filename);
end
if ~exist(masks_path, 'file')
    error('Sparse ROI mask cache not found in either channel folder: %s', masks_filename);
end
fprintf('Loading ROI mask cache (for macro pixel weights): %s\n', masks_path);
M_masks    = load(masks_path, 'roi_px_ero');
roi_px_ero = M_masks.roi_px_ero(:);
clear M_masks

%% Per-macro per-mouse pooled values (paired and delta-z) + macro-level Wilcoxon

% Figure out which macros are present (preserves DIVI list order)
macros_present = {};
for mi = 1:numel(macro_divi_list)
    mkey = matlab.lang.makeValidName(macro_divi_list{mi});
    if isfield(macro_groups, mkey)
        macros_present{end+1} = macro_divi_list{mi}; %#ok<SAGROW>
    end
end
n_macros_present = numel(macros_present);

macro_pm_nano = nan(n_mice, n_macros_present);
macro_pm_auto = nan(n_mice, n_macros_present);
macro_pm_dz   = nan(n_mice, n_macros_present);

for mi_macro = 1:n_macros_present
    mkey = matlab.lang.makeValidName(macros_present{mi_macro});
    idx  = macro_groups.(mkey);
    wts  = roi_px_ero(idx);
    for k_mouse = 1:n_mice
        n_vals = nano_pm(k_mouse, idx).';
        a_vals = auto_pm(k_mouse, idx).';
        c_vals = contrast_pm(k_mouse, idx).';
        keep_n = ~isnan(n_vals) & wts > 0;
        keep_a = ~isnan(a_vals) & wts > 0;
        keep_c = ~isnan(c_vals) & wts > 0;
        if any(keep_n)
            macro_pm_nano(k_mouse, mi_macro) = sum(n_vals(keep_n) .* wts(keep_n)) / sum(wts(keep_n));
        end
        if any(keep_a)
            macro_pm_auto(k_mouse, mi_macro) = sum(a_vals(keep_a) .* wts(keep_a)) / sum(wts(keep_a));
        end
        if any(keep_c)
            macro_pm_dz(k_mouse, mi_macro)   = sum(c_vals(keep_c) .* wts(keep_c)) / sum(wts(keep_c));
        end
    end
end

p_values_macro_paired = nan(n_macros_present, 1);
p_values_macro_dz     = nan(n_macros_present, 1);
for mi_macro = 1:n_macros_present
    nv = macro_pm_nano(:, mi_macro);
    av = macro_pm_auto(:, mi_macro);
    cv = macro_pm_dz(:, mi_macro);
    valid_pair = ~isnan(nv) & ~isnan(av);
    valid_dz   = ~isnan(cv);
    if sum(valid_pair) >= 3
        try, p_values_macro_paired(mi_macro) = signrank(nv(valid_pair), av(valid_pair), 'tail', 'right'); catch, end
    end
    if sum(valid_dz) >= 3
        try, p_values_macro_dz(mi_macro) = signrank(cv(valid_dz), 0, 'tail', 'right'); catch, end
    end
end
if apply_bonferroni
    p_thresh_macro = alpha / max(n_macros_present, 1);
else
    p_thresh_macro = alpha;
end
sig_mask_macro_paired = p_values_macro_paired < p_thresh_macro;
sig_mask_macro_dz     = p_values_macro_dz     < p_thresh_macro;

% Macro-level CSV
ptbl_macro = table( ...
    macros_present(:), p_values_macro_paired(:), sig_mask_macro_paired(:), ...
    p_values_macro_dz(:), sig_mask_macro_dz(:), ...
    mean(macro_pm_nano, 1, 'omitnan')', mean(macro_pm_auto, 1, 'omitnan')', ...
    mean(macro_pm_dz, 1, 'omitnan')', ...
    'VariableNames', {'macro', 'p_paired_nanogtauto', 'sig_paired', ...
                      'p_deltaz', 'sig_deltaz', ...
                      'nano_mean', 'auto_mean', 'contrast_mean'});
csv_name_macro = ['NanoVsAuto_macro_signrank' smooth_suffix];
if apply_bonferroni, csv_name_macro = [csv_name_macro '_bonferroni']; else, csv_name_macro = [csv_name_macro '_uncorrected']; end
writetable(ptbl_macro, fullfile(out_dir, [csv_name_macro '.csv']));
fprintf('Saved macro-level table: %s\n', fullfile(out_dir, [csv_name_macro '.csv']));

%% Figure 1: per-region paired bars (nano + auto)

n_cols_grid = 3;
n_rows_grid = ceil(n_macros_present / n_cols_grid);

if produce_paired_per_region
fig_h = figure('Visible', 'off', 'Color', 'w', 'Units', 'Normalized', 'Position', [0 0 1 1]);

for pi = 1:n_macros_present
    macro_ac = macros_present{pi};
    mkey = matlab.lang.makeValidName(macro_ac);
    idx = macro_groups.(mkey);

    sub_acro    = roi_acronyms(idx);
    sub_nano_pm = nano_pm(:, idx);
    sub_auto_pm = auto_pm(:, idx);
    sub_sig     = sig_mask(idx);

    nano_mean_v = mean(sub_nano_pm, 1, 'omitnan');
    nano_n      = sum(~isnan(sub_nano_pm), 1);
    nano_sem_v  = std(sub_nano_pm, 0, 1, 'omitnan') ./ sqrt(max(nano_n, 1));
    auto_mean_v = mean(sub_auto_pm, 1, 'omitnan');
    auto_n      = sum(~isnan(sub_auto_pm), 1);
    auto_sem_v  = std(sub_auto_pm, 0, 1, 'omitnan') ./ sqrt(max(auto_n, 1));

    % Drop regions with no data in either channel (matches P8's keep filter
    % so empty rows don't clutter the panel with bare labels).
    keep = ~isnan(nano_mean_v) & ~isnan(auto_mean_v);
    sub_acro    = sub_acro(keep);
    sub_nano_pm = sub_nano_pm(:, keep);
    sub_auto_pm = sub_auto_pm(:, keep);
    nano_mean_v = nano_mean_v(keep);
    auto_mean_v = auto_mean_v(keep);
    nano_sem_v  = nano_sem_v(keep);
    auto_sem_v  = auto_sem_v(keep);
    sub_sig     = sub_sig(keep);

    % Sort regions
    switch sort_regions_by
        case 'nano',  [~, ord] = sort(nano_mean_v, 'ascend');
        case 'auto',  [~, ord] = sort(auto_mean_v, 'ascend');
        case 'diff',  [~, ord] = sort(nano_mean_v - auto_mean_v, 'ascend');
        otherwise, error('Unknown sort_regions_by: %s', sort_regions_by);
    end

    sub_acro    = sub_acro(ord);
    sub_nano_pm = sub_nano_pm(:, ord);
    sub_auto_pm = sub_auto_pm(:, ord);
    nano_mean_v = nano_mean_v(ord);
    auto_mean_v = auto_mean_v(ord);
    nano_sem_v  = nano_sem_v(ord);
    auto_sem_v  = auto_sem_v(ord);
    sub_sig     = sub_sig(ord);

    n_sub = numel(sub_acro);
    if n_sub == 0, continue; end

    ax = subplot(n_rows_grid, n_cols_grid, pi);

    % Grouped horizontal bars: rows are regions, columns are channels
    Y = [nano_mean_v(:), auto_mean_v(:)];
    bh = barh(Y, 'grouped');
    bh(1).FaceColor = nano_color; bh(1).EdgeColor = 'none';
    bh(2).FaceColor = auto_color; bh(2).EdgeColor = 'none';
    hold on;

    % Per-bar y centres — computed manually because bh.YEndPoints semantics
    % vary across MATLAB versions for barh (some return the value-end, not
    % the position-end). MATLAB grouped-bar convention with default BarWidth:
    %   per-sub-bar width = BarWidth / n_chans
    %   sub-bar offsets within each group = ((k - (n_chans+1)/2)) * sub_w
    % For 2 channels with BarWidth=0.8: offsets are [-0.2, +0.2].
    % Convention: column 1 (nano) plots BELOW group center, column 2 (auto)
    % ABOVE — but if your MATLAB visually reverses this, just swap the two
    % y_* assignments below; bar fill colors don't change either way.
    n_chans = 2;
    bar_w = bh(1).BarWidth;
    sub_w = bar_w / n_chans;
    sub_offsets = ((1:n_chans) - (n_chans+1)/2) * sub_w;
    group_centers = 1:n_sub;
    y_nano = group_centers + sub_offsets(1);
    y_auto = group_centers + sub_offsets(2);

    % Errorbars (SEM)
    errorbar(nano_mean_v, y_nano, nano_sem_v, 'horizontal', ...
        'Color', 'k', 'LineStyle', 'none', 'LineWidth', 0.5, 'CapSize', 3);
    errorbar(auto_mean_v, y_auto, auto_sem_v, 'horizontal', ...
        'Color', 'k', 'LineStyle', 'none', 'LineWidth', 0.5, 'CapSize', 3);

    % Per-mouse paired connector lines (thin gray)
    if show_paired_lines
        for k = 1:n_sub
            for mm = 1:n_mice
                nv = sub_nano_pm(mm, k);
                av = sub_auto_pm(mm, k);
                if isnan(nv) || isnan(av), continue; end
                plot([nv, av], [y_nano(k), y_auto(k)], '-', ...
                    'Color', paired_line_color, 'LineWidth', 0.4);
            end
        end
    end

    % Per-mouse dots
    if show_per_mouse_dots
        for k = 1:n_sub
            scatter(sub_nano_pm(:, k), repmat(y_nano(k), n_mice, 1), 4, ...
                nano_dot_color, 'filled', 'MarkerFaceAlpha', 0.5);
            scatter(sub_auto_pm(:, k), repmat(y_auto(k), n_mice, 1), 4, ...
                auto_dot_color, 'filled', 'MarkerFaceAlpha', 0.5);
        end
    end

    % Highlight ytick labels: bright orange + bold for significant tests
    yticks(1:n_sub);
    lbl_fs = max(4, min(7, 120 / n_sub));
    set(gca, 'FontSize', lbl_fs);
    highlight_significant_labels(gca, sub_acro, sub_sig, lbl_fs, sig_label_color);

    xlabel('LR-sum (mean +/- SEM)');
    title(sprintf('%s (n = %d)', macro_ac, n_sub), 'Interpreter', 'none');
    grid on;
    ylim([0 n_sub + 1]);

    if pi == 1
        legend([bh(1), bh(2)], {'Nano', 'Auto'}, ...
            'Location', 'southeast', 'FontSize', 8, 'Box', 'off');
    end
end

sgtitle({['Nano vs autofluorescence per STRU leaf (paired Wilcoxon, n=' num2str(n_mice) ' mice)'], ...
         correction_label}, 'FontSize', 13, 'FontWeight', 'bold');

% Save
fname = ['Region_NanoVsAuto_BarByMacro' smooth_suffix];
if apply_bonferroni, fname = [fname '_bonferroni']; else, fname = [fname '_uncorrected']; end
set(fig_h, 'InvertHardcopy', 'off');
saveas(fig_h, fullfile(out_dir, [fname '.fig']));
exportgraphics(fig_h, fullfile(out_dir, [fname '.png']), 'Resolution', 300);
fprintf('Saved figure: %s\n', fullfile(out_dir, [fname '.png']));
end   % end of if produce_paired_per_region

%% Figure 2: per-region delta-z bars (single bar per region)

if produce_deltaz_per_region
fig_dz = figure('Visible', 'off', 'Color', 'w', 'Units', 'Normalized', 'Position', [0 0 1 1]);

for pi = 1:n_macros_present
    macro_ac = macros_present{pi};
    mkey = matlab.lang.makeValidName(macro_ac);
    idx = macro_groups.(mkey);

    sub_acro    = roi_acronyms(idx);
    sub_dz_pm   = contrast_pm(:, idx);
    sub_sig     = sig_mask_dz(idx);

    dz_mean_v = mean(sub_dz_pm, 1, 'omitnan');
    dz_n      = sum(~isnan(sub_dz_pm), 1);
    dz_sem_v  = std(sub_dz_pm, 0, 1, 'omitnan') ./ sqrt(max(dz_n, 1));

    keep = ~isnan(dz_mean_v);
    sub_acro  = sub_acro(keep);
    sub_dz_pm = sub_dz_pm(:, keep);
    dz_mean_v = dz_mean_v(keep);
    dz_sem_v  = dz_sem_v(keep);
    sub_sig   = sub_sig(keep);

    [dz_mean_v, ord] = sort(dz_mean_v, 'ascend');
    sub_acro  = sub_acro(ord);
    sub_dz_pm = sub_dz_pm(:, ord);
    dz_sem_v  = dz_sem_v(ord);
    sub_sig   = sub_sig(ord);

    n_sub = numel(sub_acro);
    if n_sub == 0, continue; end

    subplot(n_rows_grid, n_cols_grid, pi);
    b = barh(dz_mean_v);
    b.FaceColor = nano_color;     % single-color contrast bars (orange — matches nano)
    b.EdgeColor = 'none';
    hold on;

    % Errorbars (SEM) at integer y-positions
    errorbar(dz_mean_v, 1:n_sub, dz_sem_v, 'horizontal', ...
        'Color', 'k', 'LineStyle', 'none', 'LineWidth', 0.5, 'CapSize', 3);

    % Per-mouse dots
    if show_per_mouse_dots
        for k = 1:n_sub
            scatter(sub_dz_pm(:, k), repmat(k, n_mice, 1), 4, ...
                nano_dot_color, 'filled', 'MarkerFaceAlpha', 0.5);
        end
    end

    % Zero-line (contrast = 0)
    xline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.7, 'HandleVisibility', 'off');

    yticks(1:n_sub);
    lbl_fs = max(4, min(7, 120 / n_sub));
    set(gca, 'FontSize', lbl_fs);
    highlight_significant_labels(gca, sub_acro, sub_sig, lbl_fs, sig_label_color);

    xlabel('Delta-z (nano_z - auto_z)');
    title(sprintf('%s (n = %d)', macro_ac, n_sub), 'Interpreter', 'none');
    grid on;
    ylim([0 n_sub + 1]);
end

sgtitle({['Nano vs auto delta-z (per STRU leaf, paired Wilcoxon contrast > 0, n=' num2str(n_mice) ' mice)'], ...
         correction_label}, 'FontSize', 13, 'FontWeight', 'bold');

fname_dz = ['Region_NanoVsAuto_DeltaZ_BarByMacro' smooth_suffix];
if apply_bonferroni, fname_dz = [fname_dz '_bonferroni']; else, fname_dz = [fname_dz '_uncorrected']; end
set(fig_dz, 'InvertHardcopy', 'off');
saveas(fig_dz, fullfile(out_dir, [fname_dz '.fig']));
exportgraphics(fig_dz, fullfile(out_dir, [fname_dz '.png']), 'Resolution', 300);
fprintf('Saved figure: %s\n', fullfile(out_dir, [fname_dz '.png']));
end   % end of if produce_deltaz_per_region

%% Figure 3: per-macro paired bars (DIVI-level)

if produce_paired_per_macro
fig_macro = figure('Visible', 'off', 'Color', 'w', 'Units', 'Normalized', 'Position', [0 0 0.55 0.65]);

% Sort macros by nano descending
nano_macro_mean_v = mean(macro_pm_nano, 1, 'omitnan');
auto_macro_mean_v = mean(macro_pm_auto, 1, 'omitnan');
nano_macro_n   = sum(~isnan(macro_pm_nano), 1);
auto_macro_n   = sum(~isnan(macro_pm_auto), 1);
nano_macro_sem = std(macro_pm_nano, 0, 1, 'omitnan') ./ sqrt(max(nano_macro_n, 1));
auto_macro_sem = std(macro_pm_auto, 0, 1, 'omitnan') ./ sqrt(max(auto_macro_n, 1));

[~, ord_macro] = sort(nano_macro_mean_v, 'ascend');
labels_sorted        = macros_present(ord_macro);
nano_means_sorted    = nano_macro_mean_v(ord_macro);
auto_means_sorted    = auto_macro_mean_v(ord_macro);
nano_sems_sorted     = nano_macro_sem(ord_macro);
auto_sems_sorted     = auto_macro_sem(ord_macro);
nano_pm_sorted       = macro_pm_nano(:, ord_macro);
auto_pm_sorted       = macro_pm_auto(:, ord_macro);
sig_paired_sorted    = sig_mask_macro_paired(ord_macro);

n_macros_plot = numel(labels_sorted);

Y = [nano_means_sorted(:), auto_means_sorted(:)];
bh = barh(Y, 'grouped');
bh(1).FaceColor = nano_color; bh(1).EdgeColor = 'none';
bh(2).FaceColor = auto_color; bh(2).EdgeColor = 'none';
hold on;

bar_w = bh(1).BarWidth;  sub_w = bar_w / 2;
sub_offsets = ((1:2) - (2+1)/2) * sub_w;
group_centers = 1:n_macros_plot;
y_nano_m = group_centers + sub_offsets(1);
y_auto_m = group_centers + sub_offsets(2);

errorbar(nano_means_sorted, y_nano_m, nano_sems_sorted, 'horizontal', ...
    'Color', 'k', 'LineStyle', 'none', 'LineWidth', 0.7, 'CapSize', 5);
errorbar(auto_means_sorted, y_auto_m, auto_sems_sorted, 'horizontal', ...
    'Color', 'k', 'LineStyle', 'none', 'LineWidth', 0.7, 'CapSize', 5);

if show_paired_lines
    for k = 1:n_macros_plot
        for mm = 1:n_mice
            nv = nano_pm_sorted(mm, k);
            av = auto_pm_sorted(mm, k);
            if isnan(nv) || isnan(av), continue; end
            plot([nv, av], [y_nano_m(k), y_auto_m(k)], '-', ...
                'Color', paired_line_color, 'LineWidth', 0.5);
        end
    end
end

if show_per_mouse_dots
    for k = 1:n_macros_plot
        scatter(nano_pm_sorted(:, k), repmat(y_nano_m(k), n_mice, 1), 12, ...
            nano_dot_color, 'filled', 'MarkerFaceAlpha', 0.6);
        scatter(auto_pm_sorted(:, k), repmat(y_auto_m(k), n_mice, 1), 12, ...
            auto_dot_color, 'filled', 'MarkerFaceAlpha', 0.6);
    end
end

yticks(1:n_macros_plot);
set(gca, 'FontSize', 10);
highlight_significant_labels(gca, labels_sorted, sig_paired_sorted, 10, sig_label_color);
xlabel('LR-sum (mean +/- SEM)');
title({['DIVI-level nano vs auto (paired Wilcoxon, n=' num2str(n_mice) ' mice)'], correction_label}, ...
    'Interpreter', 'none', 'FontSize', 12);
grid on;
ylim([0 n_macros_plot + 1]);
legend([bh(1), bh(2)], {'Nano', 'Auto'}, 'Location', 'southeast', 'FontSize', 9, 'Box', 'off');

fname_macro = ['Macro_NanoVsAuto_Paired' smooth_suffix];
if apply_bonferroni, fname_macro = [fname_macro '_bonferroni']; else, fname_macro = [fname_macro '_uncorrected']; end
set(fig_macro, 'InvertHardcopy', 'off');
saveas(fig_macro, fullfile(out_dir, [fname_macro '.fig']));
exportgraphics(fig_macro, fullfile(out_dir, [fname_macro '.png']), 'Resolution', 300);
fprintf('Saved figure: %s\n', fullfile(out_dir, [fname_macro '.png']));
end   % end of if produce_paired_per_macro

%% Figure 4: per-macro delta-z bars

if produce_deltaz_per_macro
fig_macro_dz = figure('Visible', 'off', 'Color', 'w', 'Units', 'Normalized', 'Position', [0 0 0.45 0.55]);

dz_macro_mean_v = mean(macro_pm_dz, 1, 'omitnan');
dz_macro_n      = sum(~isnan(macro_pm_dz), 1);
dz_macro_sem    = std(macro_pm_dz, 0, 1, 'omitnan') ./ sqrt(max(dz_macro_n, 1));

[dz_macro_sorted, ord_macro_dz] = sort(dz_macro_mean_v, 'ascend');
labels_dz_sorted   = macros_present(ord_macro_dz);
dz_macro_sem_s     = dz_macro_sem(ord_macro_dz);
dz_pm_sorted       = macro_pm_dz(:, ord_macro_dz);
sig_dz_sorted      = sig_mask_macro_dz(ord_macro_dz);

b = barh(dz_macro_sorted);
b.FaceColor = nano_color;
b.EdgeColor = 'none';
hold on;

errorbar(dz_macro_sorted, 1:numel(dz_macro_sorted), dz_macro_sem_s, 'horizontal', ...
    'Color', 'k', 'LineStyle', 'none', 'LineWidth', 0.7, 'CapSize', 5);

if show_per_mouse_dots
    for k = 1:numel(dz_macro_sorted)
        scatter(dz_pm_sorted(:, k), repmat(k, n_mice, 1), 12, ...
            nano_dot_color, 'filled', 'MarkerFaceAlpha', 0.6);
    end
end

xline(0, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.7, 'HandleVisibility', 'off');

yticks(1:numel(dz_macro_sorted));
set(gca, 'FontSize', 10);
highlight_significant_labels(gca, labels_dz_sorted, sig_dz_sorted, 10, sig_label_color);
xlabel('Delta-z (nano_z - auto_z)');
title({['DIVI-level delta-z (paired Wilcoxon contrast > 0, n=' num2str(n_mice) ' mice)'], correction_label}, ...
    'Interpreter', 'none', 'FontSize', 12);
grid on;
ylim([0 numel(dz_macro_sorted) + 1]);

fname_macro_dz = ['Macro_NanoVsAuto_DeltaZ' smooth_suffix];
if apply_bonferroni, fname_macro_dz = [fname_macro_dz '_bonferroni']; else, fname_macro_dz = [fname_macro_dz '_uncorrected']; end
set(fig_macro_dz, 'InvertHardcopy', 'off');
saveas(fig_macro_dz, fullfile(out_dir, [fname_macro_dz '.fig']));
exportgraphics(fig_macro_dz, fullfile(out_dir, [fname_macro_dz '.png']), 'Resolution', 300);
fprintf('Saved figure: %s\n', fullfile(out_dir, [fname_macro_dz '.png']));
end   % end of if produce_deltaz_per_macro

%% Voxelwise z-contrast slice-by-slice video

if produce_contrast_video
fprintf('Generating voxelwise z-contrast video...\n');

% Use the cohort z-scored volumes already cached by P8 (zscore_in_mask).
% Each is z-scored within its own channel's brain mask, so the
% subtraction is a delta-z per voxel.
z_nano = S_nano_c.zscore_lr_sum;
z_auto = S_auto_c.zscore_lr_sum;
bm_int = S_nano_c.brainMask_merged & S_auto_c.brainMask_merged;
contrast_vol = z_nano - z_auto;

% Atlas boundaries from the half-width crop
allenDir = 'D:\sep_histology\data\atlas';
AllenVol = niftiread(fullfile(allenDir, 'annotation_10.nii.gz'));
limits   = [180 1079];
AllenCrop = AllenVol(limits(1):limits(2), :, :);
half_atlas = AllenCrop(:, :, 1:size(contrast_vol, 3));
clear AllenVol

% Adaptive clim symmetric around 0 (5th-95th percentile of |contrast| in brain)
abs_vals = abs(contrast_vol(bm_int & ~isnan(contrast_vol)));
clim_max = prctile(abs_vals, 95);
clim_max = max(clim_max, 1.5);   % floor so colormap is visible even for weak contrasts
clim_max = ceil(clim_max * 2) / 2;
fprintf('  Contrast clim: [-%.1f, %.1f]\n', clim_max, clim_max);

cmap_div = redblue_cmap(256);

vid_path = fullfile(out_dir, ['Contrast_video_NanoMinusAuto_z' smooth_suffix '.mp4']);
vidObj = VideoWriter(vid_path, 'MPEG-4');
vidObj.FrameRate = 15;
vidObj.Quality   = 95;
open(vidObj);

n_slices_vid = size(contrast_vol, 1);
for j = 1:n_slices_vid
    bm_sl = squeeze(bm_int(j, :, :));
    if nnz(bm_sl) == 0, continue; end

    c_sl = squeeze(contrast_vol(j, :, :));
    atl_sl = single(squeeze(half_atlas(j, :, :)));
    bnd = gradient(atl_sl) ~= 0 & atl_sl > 0;
    [br, bc] = find(bnd);

    fh = figure('visible', 'off', 'units', 'normalized', ...
        'outerposition', [0 0 1 1], 'Color', 'k');
    set(fh, 'InvertHardcopy', 'off');

    h1 = imagesc(c_sl); clim([-clim_max, clim_max]); colormap(gca, cmap_div);
    set(h1, 'AlphaData', bm_sl);
    axis image off; hold on;
    plot(bc, br, '.', 'Color', [0.5 0.5 0.5], 'MarkerSize', 0.5);
    xlim([0 size(c_sl, 2)]);
    title(sprintf('Voxelwise delta-z (nano_z - auto_z) - Slice #%d', j), ...
        'Color', 'w', 'FontSize', 11, 'Interpreter', 'none');
    cb = colorbar; cb.Color = 'w';
    cb.Label.String = 'Delta-z';
    cb.Label.FontSize = 10; cb.Label.Color = 'w';
    set(gca, 'Color', 'k');

    frame = getframe(fh);
    writeVideo(vidObj, frame);
    close(fh);
    if mod(j, 50) == 0
        fprintf('    Frame %d / %d written\n', j, n_slices_vid);
    end
end
close(vidObj);
fprintf('Saved video: %s\n', vid_path);
end   % end of if produce_contrast_video

fprintf('P10 done. Outputs in: %s\n', out_dir);

%% Local function: highlight significant region labels in orange-bold
% Mirrors P8's highlight_enriched_labels but takes an explicit color so the
% intent (significant test result) is decoupled from the threshold-based
% enrichment used in P8.

function highlight_significant_labels(ax_handle, acro_list, sig_mask, fontsize, sig_color)
    yticklabels(ax_handle, acro_list);
    drawnow;
    try
        tl = ax_handle.YAxis.TickLabels;
        for k = 1:numel(acro_list)
            if sig_mask(k)
                tl(k).FontColor  = sig_color;
                tl(k).FontWeight = 'bold';
            end
        end
    catch
        % Fallback for older MATLAB: blank yticklabels and place text() manually
        yticklabels(ax_handle, {});
        xl = xlim(ax_handle);
        for k = 1:numel(acro_list)
            clr = [0 0 0]; fw = 'normal';
            if sig_mask(k), clr = sig_color; fw = 'bold'; end
            text(ax_handle, xl(1) - 0.02*(xl(2)-xl(1)), k, acro_list{k}, ...
                'Color', clr, 'FontSize', fontsize, 'FontWeight', fw, ...
                'HorizontalAlignment', 'right', 'Units', 'data', ...
                'Interpreter', 'none', 'Clipping', 'off');
        end
    end
end

%% Local function: diverging blue-white-red colormap for the contrast video
function cmap = redblue_cmap(n)
    if nargin < 1, n = 256; end
    half = floor(n/2);
    % Blue (low) -> white (mid) -> red (high)
    blue_to_white = [linspace(0.13, 1, half)', ...
                     linspace(0.40, 1, half)', ...
                     linspace(0.67, 1, half)'];
    white_to_red  = [linspace(1, 0.85, n - half)', ...
                     linspace(1, 0.20, n - half)', ...
                     linspace(1, 0.10, n - half)'];
    cmap = [blue_to_white; white_to_red];
end
