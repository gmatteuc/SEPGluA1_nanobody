function perm = region_permutation_test(stacks, n_ctrl, geom, perm_settings)
%REGION_PERMUTATION_TEST  Region scores of the surprise maps, with an exact permutation test.
%   perm = REGION_PERMUTATION_TEST(stacks, n_ctrl, geom, perm_settings) scores
%   each region of the surprise bars with five measures, for every split of the
%   pooled mice into two groups of the original sizes, and gives each score its
%   permutation p, uncorrected and corrected over the regions. Called by
%   group_differences (run_group_differences), by per_mouse_region_values
%   (run_per_mouse_values) for the observed split alone, and by
%   tests\test_region_permutation on synthetic stacks.
%
%   stacks         a cell, one array per map (L - R, L + R): n_cand x n_mice
%                  single, each mouse's absolute value at every candidate voxel,
%                  NaN where the mouse has none; the n_ctrl control mice first
%   geom           the candidate voxels (those where a split can give a t):
%                  grid_size (AP x DV x ML of the folded maps), cand_lin (their
%                  linear indices, ascending), cand_region (their region of the
%                  bars, 0 for none), n_regions, voxel_mm (the voxel's side)
%   perm_settings  min_mice_per_group, slab_range and p_thresh (the bars'),
%                  cluster_p, cluster_connectivity, topvol_mm3, region_quantile,
%                  n_permutations ('all' or a number), n_workers, seed, and
%                  keep_cluster_voxels (optional, false when missing: true
%                  keeps the voxels of the observed split's heaviest clusters)
%
%   For each split everything that feeds the bars is computed again, as
%   group_differences computes it for the groups as they are: the Welch t of
%   the group means at each voxel where each group has min_mice_per_group mice
%   with a value, its surprise -log10 p, here with the sign of the difference
%   (experimental minus control), and the median of the surprise over +/-
%   slab_range planes on the voxels with a t. The measures, each for the
%   positive and the negative effects apart, the region's score being the
%   larger of the two with its sign:
%     share    the fraction of the voxels with a t at p < p_thresh: the bars,
%              from the median of the unsigned surprise, so with no sign
%     sum      the summed surprise of the voxels at p < p_thresh
%     q99      the region_quantile of the surprise over its voxels with a t
%              (named after it: q95 for 0.95)
%     topvol   the mean surprise of its K most surprising voxels, K the number
%              of voxels in topvol_mm3; all of them in a region with fewer
%     cluster  the mass (summed surprise) of its heaviest cluster of voxels at
%              p < cluster_p, connected within the region (bwconncomp with
%              cluster_connectivity)
%   A region's p is the fraction of the splits whose |score| reaches the
%   observed one, the observed split included; its corrected p the fraction
%   whose largest |score| over the regions reaches it (family-wise, per map and
%   measure). With groups of equal size, a split and its mirror image (the
%   groups swapped) give opposite t everywhere, so only one of the two is
%   computed and the other's scores are its own with the signs swapped.
%
%   perm holds the splits (one row per split, the mice labelled control), the
%   measures' names and, per map, the observed scores with their p and
%   corrected p, the null of every split, the observed heaviest clusters'
%   voxels and peaks, and the observed rolling median of the unsigned surprise
%   on the candidate voxels (to check it against the bars'); with
%   keep_cluster_voxels, also the heaviest clusters' voxels, as linear indices
%   into the grid (detail.cluster_voxels, regions x positive and negative).

% settings, under the names the code below uses
min_mice_per_group = perm_settings.min_mice_per_group;
n_permutations = perm_settings.n_permutations;
n_workers = perm_settings.n_workers;
seed = perm_settings.seed;

% the voxels of the observed split's heaviest clusters, only when asked: the bars
% do not need them, and their null file would carry every cluster's voxels
keep_cluster_voxels = isfield(perm_settings, 'keep_cluster_voxels') && ...
    perm_settings.keep_cluster_voxels;

% the measures, the quantile's named after it
quantile_name = sprintf('q%g', 100 * perm_settings.region_quantile);
measure_names = {'share', 'sum', quantile_name, 'topvol', 'cluster'};

% the top volume in voxels, from the voxel size
topvol_k = round(perm_settings.topvol_mm3 / geom.voxel_mm^3);

%% Splits

n_mice = size(stacks{1}, 2);
n_exp = n_mice - n_ctrl;
if n_ctrl < min_mice_per_group || n_exp < min_mice_per_group
    error(['region_permutation_test: %d control and %d experimental mice, but a t ' ...
           'needs %d per group (min_mice_per_group).'], n_ctrl, n_exp, min_mice_per_group);
end

% every split, or a random subset of them with the observed one, and the splits
% computed (one of each mirror pair)
splits = enumerate_splits(n_mice, n_ctrl, n_permutations, seed);
n_splits = size(splits.in_ctrl, 1);
n_computed = numel(splits.computed);
fprintf(['  %d splits of %d mice into %d and %d (%d in all), %d computed, %d as ' ...
         'their mirror image.\n'], n_splits, n_mice, n_ctrl, n_exp, ...
        nchoosek(n_mice, n_ctrl), n_computed, n_splits - n_computed);

%% Candidate voxels

% the blocks of AP columns the maps are computed in, and the voxels of each
% region
geom = prepare_geometry(geom, perm_settings);
fprintf(['  %d candidate voxels in %d AP columns, %d blocks; top volume %d voxels ' ...
         '(%g mm^3).\n'], numel(geom.cand_lin), geom.n_cols, geom.n_blocks, topvol_k, ...
        perm_settings.topvol_mm3);

%% Observed split

% the groups as they are, with the clusters' voxels and peaks and the rolling
% median the bars use; timed, for the estimate of the rest
t_observed = tic;
observed = split_scores(stacks, splits.in_ctrl(splits.observed, :), geom, ...
    perm_settings, topvol_k, true, keep_cluster_voxels);
observed_s = toc(t_observed);
fprintf('  observed split: %.1f s.\n', observed_s);

%% Every other split

% the splits to compute besides the observed one
other = splits.computed(splits.computed ~= splits.observed);
n_other = numel(other);
in_ctrl_other = splits.in_ctrl(other, :);

% a thread pool: its workers share the stacks, where a process pool would copy
% them into every worker
if n_workers > 0
    pool = gcp('nocreate');
    if isempty(pool)
        pool = parpool('Threads', n_workers);
    elseif ~isa(pool, 'parallel.ThreadPool')
        error(['region_permutation_test: a process pool is open, which would copy ' ...
               'the stacks into each worker. Close it (delete(gcp)) and rerun: a ' ...
               'thread pool of %d workers is started.'], n_workers);
    end
    n_workers = pool.NumWorkers;
end
fprintf('  %d splits to compute on %d workers.\n', n_other, max(n_workers, 1));

% progress: a line every twentieth of the splits, with the time left
queue = parallel.pool.DataQueue;
afterEach(queue, @progress_report);
progress_report({'start', n_other});

null_other = cell(n_other, 1);
parfor (k = 1:n_other, n_workers)
    null_other{k} = split_scores(stacks, in_ctrl_other(k, :), geom, perm_settings, ...
        topvol_k, false, false);
    send(queue, k);
end

%% P-values

% every split's positive and negative scores, the mirror images from their
% computed split
[null_pos, null_neg] = assemble_null(splits, observed, null_other, other);

% each map's observed scores, their p and corrected p
perm = struct();
perm.splits = splits;
perm.measure_names = measure_names;
perm.topvol_k = topvol_k;
perm.settings = perm_settings;
perm.observed_s = observed_s;
perm.maps = cell(numel(stacks), 1);
for m = 1:numel(stacks)
    perm.maps{m} = map_p_values(null_pos{m}, null_neg{m}, splits.observed, ...
        measure_names, observed, m);
end
end

% ===== Local functions: splits and geometry =====

function splits = enumerate_splits(n_mice, n_ctrl, n_permutations, seed)
% Every split of the pooled mice into n_ctrl control and n_mice - n_ctrl
% experimental mice (logical rows, true for control), or a random subset with
% the observed split; which to compute, and which are the mirror image of one.

% every split, the observed one being mice 1 to n_ctrl as control
all_ctrl = nchoosek(1:n_mice, n_ctrl);
n_all = size(all_ctrl, 1);
in_ctrl = false(n_all, n_mice);
for s = 1:n_all
    in_ctrl(s, all_ctrl(s, :)) = true;
end
observed_row = find(all(in_ctrl == (1:n_mice <= n_ctrl), 2));
if numel(observed_row) ~= 1
    error('region_permutation_test: the observed split is not among the %d splits.', ...
        n_all);
end

% a random subset with a fixed seed when asked, the observed split always in it
if ischar(n_permutations) || isstring(n_permutations)
    if ~strcmp(n_permutations, 'all')
        error(['region_permutation_test: n_permutations is ''%s''; use ''all'' or ' ...
               'a number of splits.'], n_permutations);
    end
    keep = (1:n_all)';
elseif n_permutations >= n_all
    keep = (1:n_all)';
else
    rng_state = rng(seed);
    others = setdiff(1:n_all, observed_row);
    drawn = others(randperm(numel(others), n_permutations - 1));
    rng(rng_state);
    keep = sort([observed_row, drawn])';
end
in_ctrl = in_ctrl(keep, :);
observed = find(keep == observed_row);

% with groups of equal size, a split's mirror image (the groups swapped) is
% computed from the split itself: the one where mouse 1 is control
n_kept = numel(keep);
source = (1:n_kept)';
is_mirror = false(n_kept, 1);
if 2 * n_ctrl == n_mice
    for s = 1:n_kept
        if ~in_ctrl(s, 1)
            [~, partner] = ismember(~in_ctrl(s, :), in_ctrl, 'rows');
            if partner > 0
                source(s) = partner;
                is_mirror(s) = true;
            end
        end
    end
end

splits = struct();
splits.in_ctrl = in_ctrl;
splits.observed = observed;
splits.source = source;
splits.is_mirror = is_mirror;
splits.computed = find(~is_mirror);
end

function geom = prepare_geometry(geom, perm_settings)
% The candidates' AP columns (a column: one DV and ML position, every plane),
% each candidate's place in an AP x columns array, blocks of columns for the
% rolling median, and the candidates of each region.

n_ap = geom.grid_size(1);
cand_lin = double(geom.cand_lin(:));
if any(diff(cand_lin) <= 0)
    error('region_permutation_test: the candidate voxels must be in ascending order.');
end

% each candidate's column and plane; the columns with a candidate, numbered in
% order (the candidates are in ascending order, so their columns are too), and
% each candidate's place in their AP x columns array (AP runs first in the
% volume, so each column is one run of n_ap voxels)
column = floor((cand_lin - 1) / n_ap) + 1;
plane = cand_lin - (column - 1) * n_ap;
column_number = cumsum([1; diff(column) > 0]);
clear column cand_lin
geom.n_cols = column_number(end);
geom.cand_zlin = uint32(plane + n_ap * (column_number - 1));

% blocks of whole columns, about block_voxels candidates each, so a worker never
% holds a whole map at once
block_voxels = 4e6;
per_column = accumarray(column_number, 1, [geom.n_cols 1]);
column_block = ceil(cumsum(per_column) / block_voxels);
[~, geom.block_col_first] = unique(column_block, 'first');
[~, geom.block_col_last] = unique(column_block, 'last');
geom.n_blocks = numel(geom.block_col_first);
cand_before = [0; cumsum(per_column)];
geom.block_first = cand_before(geom.block_col_first) + 1;
geom.block_last = cand_before(geom.block_col_last + 1);

% the candidates of each region, in region order, and where each region's run
% starts and ends
cand_region = double(geom.cand_region(:));
in_region = find(cand_region > 0);
[sorted_region, order] = sort(cand_region(in_region));
geom.region_order = uint32(in_region(order));
per_region = accumarray(sorted_region, 1, [geom.n_regions 1]);
geom.region_last = cumsum(per_region);
geom.region_first = geom.region_last - per_region + 1;

% the thresholds, as surprise
geom.thresh_sig = -log10(perm_settings.p_thresh);
geom.thresh_cluster = -log10(perm_settings.cluster_p);
end

% ===== Local functions: one split =====

function scores = split_scores(stacks, in_ctrl, geom, perm_settings, topvol_k, ...
    keep_detail, keep_cluster_voxels)
% The five measures of every region, for each map, positive and negative apart
% (n_regions x 5 each, the share in both), for one split of the mice; with
% keep_detail, also the heaviest clusters' voxels and peaks, the number of
% voxels behind each top volume, and the rolling median of the unsigned surprise;
% with keep_cluster_voxels too, the heaviest clusters' voxels themselves.

n_maps = numel(stacks);
scores = struct();
scores.pos = cell(n_maps, 1);
scores.neg = cell(n_maps, 1);
scores.detail = cell(n_maps, 1);
for m = 1:n_maps

    % the signed surprise, rolled over planes, and the share of each region
    [rolled, share, n_with_t, rolled_unsigned] = rolled_surprise(stacks{m}, in_ctrl, ...
        geom, perm_settings, keep_detail);

    % the sum, quantile and top volume of each region, positive and negative
    [pos, neg, topvol_n] = region_measures(rolled, geom, perm_settings, topvol_k);

    % the heaviest cluster of each region, positive and negative
    [mass_pos, mass_neg, cluster_n, cluster_peak, cluster_voxels] = region_clusters( ...
        rolled, n_with_t, geom, perm_settings, keep_detail, keep_cluster_voxels);

    % the five measures in the order of measure_names; the share has no sign
    scores.pos{m} = [share, pos, mass_pos];
    scores.neg{m} = [share, neg, mass_neg];
    if keep_detail
        detail = struct();
        detail.n_with_t = n_with_t;
        detail.topvol_n = topvol_n;
        detail.cluster_n = cluster_n;
        detail.cluster_peak = cluster_peak;
        detail.rolled_unsigned = rolled_unsigned;
        if keep_cluster_voxels
            detail.cluster_voxels = cluster_voxels;
        end
        scores.detail{m} = detail;
    end
end
end

function [rolled, share, n_with_t, rolled_unsigned] = rolled_surprise(stack, in_ctrl, ...
    geom, perm_settings, keep_unsigned)
% The signed surprise of the split's t on the candidate voxels, as its median over
% +/- slab_range planes on the voxels with a t (NaN without a t), and each
% region's share of voxels with a t at p < p_thresh, from the median of the
% unsigned surprise, as the bars; block by block of AP columns.

slab_range = perm_settings.slab_range;
n_ap = geom.grid_size(1);
n_cand = numel(geom.cand_lin);
rolled = nan(n_cand, 1, 'single');
unsigned_blocks = cell(geom.n_blocks, 1);
n_with_t = zeros(geom.n_regions, 1);
n_significant = zeros(geom.n_regions, 1);
for b = 1:geom.n_blocks
    rows = geom.block_first(b):geom.block_last(b);

    % the signed surprise of the block's voxels, and those with a t
    [surprise, has_t] = signed_surprise(stack(rows, in_ctrl), stack(rows, ~in_ctrl), ...
        perm_settings.min_mice_per_group);

    % the block as planes x columns, NaN off the voxels with a t (the surprise is
    % NaN there)
    n_cols_block = geom.block_col_last(b) - geom.block_col_first(b) + 1;
    zlin = geom.cand_zlin(rows) - uint32(n_ap * (geom.block_col_first(b) - 1));
    planes = nan(n_ap, n_cols_block, 'single');
    planes(zlin) = surprise;

    % the median over the planes around each plane, the window cut at the ends of
    % the volume as in the bars (movmedian's shrinking ends), kept on the voxels
    % with a t; on the unsigned surprise it is the bars' median to the bit
    signed_median = movmedian(planes, [slab_range slab_range], 1, 'omitnan');
    block_rolled = signed_median(zlin);
    block_rolled(~has_t) = NaN;
    rolled(rows) = block_rolled;
    unsigned_median = movmedian(abs(planes), [slab_range slab_range], 1, 'omitnan');
    block_unsigned = unsigned_median(zlin);
    block_unsigned(~has_t) = NaN;
    if keep_unsigned
        unsigned_blocks{b} = block_unsigned;
    end

    % each region's voxels with a t, and those significant, as the bars count them
    region = double(geom.cand_region(rows));
    in_region = region > 0;
    n_with_t = n_with_t + accumarray(region(in_region), ...
        double(~isnan(block_unsigned(in_region))), [geom.n_regions 1]);
    n_significant = n_significant + accumarray(region(in_region), ...
        double(block_unsigned(in_region) > geom.thresh_sig), [geom.n_regions 1]);
end

% 0 / 0 is NaN, for a region without a voxel with a t
share = n_significant ./ n_with_t;

% the unsigned median of every block, in order (empty unless kept)
rolled_unsigned = vertcat(unsigned_blocks{:});
end

function [surprise, has_t] = signed_surprise(x_ctrl, x_exp, min_mice_per_group)
% The surprise -log10 p of the Welch t of experimental minus control, with the
% sign of the t, at each voxel (a row: the mice of each group as columns, NaN
% without a value), NaN where a group has fewer than min_mice_per_group mice;
% the same arithmetic as group_differences' group_welch_t and welch_surprise,
% so the observed split gives their surprise to the bit.

% each group's mice with a value, mean and SEM over them
n_ctrl = single(sum(~isnan(x_ctrl), 2));
n_exp = single(sum(~isnan(x_exp), 2));
avg_ctrl = nanmean(x_ctrl, 2); %#ok<*NANMEAN>
avg_exp = nanmean(x_exp, 2);
sem_ctrl = nanstd(x_ctrl, [], 2) ./ sqrt(double(n_ctrl)); %#ok<*NANSTD>
sem_exp = nanstd(x_exp, [], 2) ./ sqrt(double(n_exp));

% a zero SEM to NaN, so the t is NaN rather than infinite
sem_ctrl(sem_ctrl == 0) = NaN;
sem_exp(sem_exp == 0) = NaN;

% the Welch t, only where each group has enough mice
sem_diff = sqrt(sem_ctrl.^2 + sem_exp.^2);
sem_diff(isnan(sem_diff) | sem_diff == 0) = NaN;
t = (avg_exp - avg_ctrl) ./ sem_diff;
has_t = n_ctrl >= min_mice_per_group & n_exp >= min_mice_per_group;
t(~has_t) = NaN;

% Welch-Satterthwaite degrees of freedom, and the surprise of the two-sided p
var_ctrl = sem_ctrl.^2;
var_exp = sem_exp.^2;
df = (var_ctrl + var_exp).^2 ./ ...
    (var_ctrl.^2 ./ (double(n_ctrl) - 1) + var_exp.^2 ./ (double(n_exp) - 1));
df(~has_t) = NaN;
p = 2 * tcdf(-abs(t), df);
surprise = -log10(p) .* sign(t);
end

function [pos, neg, topvol_n] = region_measures(rolled, geom, perm_settings, topvol_k)
% Each region's summed surprise at p < p_thresh, quantile and top volume, for
% the positive and the negative effects apart (n_regions x 3 each, magnitudes),
% and the voxels each top volume averaged; NaN for a region without a voxel
% with a t.

n_regions = geom.n_regions;
pos = nan(n_regions, 3);
neg = nan(n_regions, 3);
topvol_n = zeros(n_regions, 1);
for r = 1:n_regions

    % the region's voxels with a t
    values = rolled(geom.region_order(geom.region_first(r):geom.region_last(r)));
    values = values(~isnan(values));
    if isempty(values)
        continue
    end

    % the positive and the negative effects, each as a magnitude, zero elsewhere
    up = max(values, 0);
    down = max(-values, 0);

    % the summed surprise of the significant voxels of each sign
    pos(r, 1) = sum(values(values > geom.thresh_sig), 'double');
    neg(r, 1) = -sum(values(values < -geom.thresh_sig), 'double');

    % the quantile of each, over all the voxels with a t
    pos(r, 2) = quantile(up, perm_settings.region_quantile);
    neg(r, 2) = quantile(down, perm_settings.region_quantile);

    % the mean of the topvol_k largest of each, or of all in a smaller region
    topvol_n(r) = min(topvol_k, numel(values));
    pos(r, 3) = mean(maxk(up, topvol_n(r)), 'double');
    neg(r, 3) = mean(maxk(down, topvol_n(r)), 'double');
end
end

function [mass_pos, mass_neg, cluster_n, cluster_peak, cluster_voxels] = ...
    region_clusters(rolled, n_with_t, geom, perm_settings, keep_detail, ...
    keep_cluster_voxels)
% Each region's heaviest cluster of voxels at p < cluster_p, connected within
% the region, positive and negative apart: its mass (summed |surprise|), and
% with keep_detail its voxels and its signed peak (n_regions x 2, positive then
% negative); 0 for a region without one, NaN without a voxel with a t. With
% keep_cluster_voxels too, the clusters' voxels as linear indices into the grid
% (a cell, n_regions x 2, empty for a region without one).

n_regions = geom.n_regions;
connectivity = perm_settings.cluster_connectivity;

% the significant voxels in a region, and their label: the region and the sign
% (2r - 1 positive, 2r negative), since a cluster stays within one of each
is_sig = (rolled > geom.thresh_cluster | rolled < -geom.thresh_cluster) & ...
    geom.cand_region(:) > 0;
sig_cand = find(is_sig);
clear is_sig
sig_lin = double(geom.cand_lin(sig_cand));
sig_value = double(rolled(sig_cand));
sig_label = 2 * double(geom.cand_region(sig_cand)) - (sig_value > 0);

% the clusters and the heaviest of each label
[cluster, n_clusters] = label_clusters(sig_lin, sig_label, geom.grid_size, ...
    connectivity);
mass = accumarray(cluster, abs(sig_value), [n_clusters 1]);
n_voxels = accumarray(cluster, 1, [n_clusters 1]);
best_cluster = heaviest_by_label(cluster, sig_label, mass, n_voxels, 2 * n_regions);
mass_by_label = zeros(2 * n_regions, 1);
mass_by_label(best_cluster > 0) = mass(best_cluster(best_cluster > 0));

% per region: positive (odd labels) and negative (even labels); NaN for a
% region without a voxel with a t
mass_pos = mass_by_label(1:2:end);
mass_neg = mass_by_label(2:2:end);
mass_pos(n_with_t == 0) = NaN;
mass_neg(n_with_t == 0) = NaN;

% the voxels and the signed peak of the heaviest clusters, for the observed split
cluster_n = [];
cluster_peak = [];
cluster_voxels = {};
if keep_detail
    peak = accumarray(cluster, abs(sig_value), [n_clusters 1], @max);
    best = reshape(best_cluster, 2, n_regions)';
    cluster_n = zeros(n_regions, 2);
    cluster_peak = zeros(n_regions, 2);
    cluster_n(best > 0) = n_voxels(best(best > 0));
    cluster_peak(best > 0) = peak(best(best > 0));
    cluster_peak(:, 2) = -cluster_peak(:, 2);

    % the voxels of each, when asked
    if keep_cluster_voxels
        cluster_voxels = cell(n_regions, 2);
        for r = 1:n_regions
            for s = 1:2
                if best(r, s) > 0
                    cluster_voxels{r, s} = sig_lin(cluster == best(r, s));
                end
            end
        end
    end
end
end

function [cluster, n_clusters] = label_clusters(sig_lin, sig_label, grid_size, ...
    connectivity)
% The cluster of each significant voxel (linear indices in the grid, ascending):
% connected components of voxels of one label. The components of all the voxels
% at once, on the whole grid; one holding two labels (two regions, or both
% signs, that touch) is split by label, each part into its own components, so a
% cluster within a region is one of the region's own voxels, the other regions
% left out. Some cluster numbers end up without a voxel.

cluster = zeros(0, 1);
n_clusters = 0;
if isempty(sig_lin)
    return
end

% the components of all the significant voxels
mask = false(grid_size);
mask(sig_lin) = true;
components = bwconncomp(mask, connectivity);
clear mask
cluster = component_of_voxel(components, numel(sig_lin));
n_clusters = components.NumObjects;

% the components holding two labels or more: the variance of their labels is not
% zero (exact in double for these small integers)
label_n = accumarray(cluster, 1, [n_clusters 1]);
label_sum = accumarray(cluster, sig_label, [n_clusters 1]);
label_sumsq = accumarray(cluster, sig_label.^2, [n_clusters 1]);
is_mixed = label_n .* label_sumsq - label_sum.^2 > 0;
in_mixed = find(is_mixed(cluster));
if isempty(in_mixed)
    return
end

% each part (one component, one label) into its own components, in its bounding box
[~, ~, part] = unique([cluster(in_mixed), sig_label(in_mixed)], 'rows');
[part_sorted, by_part] = sort(part);
part_last = find([diff(part_sorted); 1]);
part_first = [1; part_last(1:end - 1) + 1];
for g = 1:numel(part_first)
    members = in_mixed(by_part(part_first(g):part_last(g)));
    [sub_cluster, n_sub] = clusters_in_box(sig_lin(members), grid_size, connectivity);
    cluster(members) = n_clusters + sub_cluster;
    n_clusters = n_clusters + n_sub;
end
end

function best_cluster = heaviest_by_label(cluster, sig_label, mass, n_voxels, n_labels)
% The heaviest cluster of each label (0 for a label without one).

best_cluster = zeros(n_labels, 1);
if isempty(cluster)
    return
end

% each cluster's label (one per cluster), then the clusters by mass, heaviest
% first, and the first of each label
cluster_label = accumarray(cluster, sig_label, [numel(mass) 1], @max);
has_voxels = find(n_voxels > 0);
[~, by_mass] = sort(mass(has_voxels), 'descend');
heaviest_first = has_voxels(by_mass);
[labels_present, first] = unique(cluster_label(heaviest_first), 'first');
best_cluster(labels_present) = heaviest_first(first);
end

function cluster = component_of_voxel(components, n_voxels)
% The component of each voxel of a bwconncomp result, in the ascending order of
% the voxels' linear indices.

per_component = cellfun(@numel, components.PixelIdxList);
pixels = vertcat(components.PixelIdxList{:});
component = repelem((1:components.NumObjects)', per_component(:));

% a column even for one component (repelem of a scalar gives a row)
[~, order] = sort(pixels);
cluster = reshape(component(order), [], 1);
if numel(cluster) ~= n_voxels
    error('region_permutation_test: %d voxels in the components, %d expected.', ...
        numel(cluster), n_voxels);
end
end

function [sub_cluster, n_sub] = clusters_in_box(voxel_lin, grid_size, connectivity)
% The connected components of a set of voxels (linear indices in the grid,
% ascending), in their bounding box: each voxel's component and their number.

[a, d, l] = ind2sub(grid_size, voxel_lin);
low = [min(a), min(d), min(l)];
box_size = [max(a), max(d), max(l)] - low + 1;
box = false(box_size);
box_lin = sub2ind(box_size, a - low(1) + 1, d - low(2) + 1, l - low(3) + 1);
box(box_lin) = true;
components = bwconncomp(box, connectivity);
n_sub = components.NumObjects;

% the box's linear order is not the grid's: back to the voxels' order
[~, order] = sort(box_lin);
in_box_order = component_of_voxel(components, numel(voxel_lin));
sub_cluster = zeros(numel(voxel_lin), 1);
sub_cluster(order) = in_box_order;
end

% ===== Local functions: p-values and progress =====

function [null_pos, null_neg] = assemble_null(splits, observed, null_other, other)
% Every split's positive and negative scores per map (n_splits x n_regions x 5),
% a mirror image taking its computed split's scores with the signs swapped.

n_splits = size(splits.in_ctrl, 1);
n_maps = numel(observed.pos);
[n_regions, n_measures] = size(observed.pos{1});

% the computed splits' scores, by split
computed_pos = cell(n_splits, 1);
computed_neg = cell(n_splits, 1);
computed_pos{splits.observed} = observed.pos;
computed_neg{splits.observed} = observed.neg;
for k = 1:numel(other)
    computed_pos{other(k)} = null_other{k}.pos;
    computed_neg{other(k)} = null_other{k}.neg;
end

null_pos = cell(n_maps, 1);
null_neg = cell(n_maps, 1);
for m = 1:n_maps
    null_pos{m} = nan(n_splits, n_regions, n_measures);
    null_neg{m} = nan(n_splits, n_regions, n_measures);
    for s = 1:n_splits

        % a mirror image: the groups swapped, so its positive effects are its
        % source's negative ones
        source = splits.source(s);
        if splits.is_mirror(s)
            null_pos{m}(s, :, :) = computed_neg{source}{m};
            null_neg{m}(s, :, :) = computed_pos{source}{m};
        else
            null_pos{m}(s, :, :) = computed_pos{source}{m};
            null_neg{m}(s, :, :) = computed_neg{source}{m};
        end
    end
end
end

function map = map_p_values(null_pos, null_neg, observed_row, measure_names, ...
    observed, m)
% One map's observed scores (the larger of the positive and the negative, with
% its sign), their permutation p and corrected p, per region and measure, with
% the null and the observed split's details.

% the scores as magnitudes: the larger of the two signs; NaN for a region without
% a voxel with a t, which never reaches a score
[n_splits, n_regions, n_measures] = size(null_pos);
null_abs = max(null_pos, null_neg);
observed_abs = reshape(null_abs(observed_row, :, :), n_regions, n_measures);
observed_pos = reshape(null_pos(observed_row, :, :), n_regions, n_measures);
observed_neg = reshape(null_neg(observed_row, :, :), n_regions, n_measures);

% the sign: positive when the positive effect is at least as large; none for the
% share, and none for a zero score
sign_obs = ones(size(observed_abs));
sign_obs(observed_neg > observed_pos) = -1;
sign_obs(observed_abs == 0 | isnan(observed_abs)) = 0;
sign_obs(:, strcmp(measure_names, 'share')) = 0;

% the region's p: the splits whose |score| reaches the observed one (the observed
% split among them); the corrected p: those whose largest |score| over the
% regions does
null_max = reshape(max(null_abs, [], 2, 'omitnan'), n_splits, n_measures);
p_perm = nan(size(observed_abs));
p_fwer = nan(size(observed_abs));
for k = 1:n_measures
    for r = 1:n_regions
        if ~isnan(observed_abs(r, k))
            p_perm(r, k) = nnz(null_abs(:, r, k) >= observed_abs(r, k)) / n_splits;
            p_fwer(r, k) = nnz(null_max(:, k) >= observed_abs(r, k)) / n_splits;
        end
    end
end

map = struct();
map.score = sign_obs .* observed_abs;
map.score(isnan(observed_abs)) = NaN;
map.score(:, strcmp(measure_names, 'share')) = observed_abs(:, ...
    strcmp(measure_names, 'share'));
map.sign = sign_obs;
map.p_perm = p_perm;
map.p_fwer = p_fwer;
map.null_pos = null_pos;
map.null_neg = null_neg;
map.null_max = null_max;
map.detail = observed.detail{m};
end

function progress_report(message)
% A line after the first finished split, then every twentieth of them and after
% the last, with the time since the start and the time left; started by a cell
% {'start', n_total}.

persistent n_total n_done t_start
if iscell(message)
    n_total = message{2};
    n_done = 0;
    t_start = tic;
    return
end
n_done = n_done + 1;
step = max(1, round(n_total / 20));
if n_done == 1 || mod(n_done, step) == 0 || n_done == n_total
    elapsed = toc(t_start);
    remaining = elapsed / n_done * (n_total - n_done);
    fprintf('    split %d of %d done, %.1f min, about %.1f min left\n', n_done, ...
        n_total, elapsed / 60, remaining / 60);
end
end
