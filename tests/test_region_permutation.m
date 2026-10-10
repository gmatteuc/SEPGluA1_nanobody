%% test_region_permutation
% ===== Check the region measures and their permutation test on synthetic maps =====
%
% Known-answer test of region_permutation_test, the permutation test of the
% plasticity comparison's region bars (run_group_differences), on small
% synthetic stacks of 5 against 5 mice: six regions, the L - R map pure noise,
% the L + R map with (a) a strong focal bump of 400 voxels inside the largest
% region and a weak shift over the whole of a small one, then (b) pure noise
% under several seeds. Checks that every split is enumerated with the observed
% one among them, that the observed split gives the bars' surprise median to
% the bit (the arithmetic of group_differences, written out below), that its
% sum, quantile, top volume and heaviest cluster equal those written out
% region by region (each region's clusters from bwconncomp of its own voxels
% of one sign), also on (c) a map whose clusters cross a region border and
% touch the other sign, and whose heaviest cluster differs at 6, 18 and 26
% connectivity, where the voxels it keeps of each heaviest cluster on request
% (keep_cluster_voxels) must be that cluster's and the t it keeps on request
% (keep_t) the t written out, that a split's mirror image is
% its negative, that the cluster
% and top-volume measures put the large region first at a corrected p < 0.05
% in (a), and that in (b) no measure reaches a corrected p < 0.05 more often
% than chance; prints how every measure ranks the two regions of (a). Prints
% PASS or FAIL for each check, and stops with an error if any fails.
%
% Run it after any change to region_permutation_test, or to the t, the surprise
% or the rolling median of group_differences. It needs no data and writes
% nothing; it opens a pool of 4 thread workers. Run sep_setup_paths first, once
% per MATLAB session (about 2 minutes).

clear; clc; close all;

%% Settings

% the grid (AP x DV x ML) and its voxel side, which makes the top volume of
% 0.1 mm^3 about 195 voxels, half the bump of (a)
grid_size = [60 40 40];
voxel_mm = 0.08;

% the mice: 5 control, then 5 experimental
n_ctrl = 5;
n_mice = 10;

% the noise of each mouse: Gaussian, smoothed in 3D, then scaled to SD 1 around
% a mean of 10; the effects of (a) in SD units, on the experimental mice only
noise_sigma = 2;
bump_amplitude = 4;
shift_amplitude = 2;

% seeds of (b)
null_seeds = 1:6;

% the permutation test as in run_group_differences, all 252 splits
perm_settings = struct();
perm_settings.min_mice_per_group = 3;
perm_settings.slab_range = 10;
perm_settings.p_thresh = 0.01;
perm_settings.cluster_p = 0.01;
perm_settings.cluster_connectivity = 18;
perm_settings.topvol_mm3 = 0.1;
perm_settings.region_quantile = 0.99;
perm_settings.n_permutations = 'all';
perm_settings.n_workers = 4;
perm_settings.seed = 0;

n_pass = 0;
n_fail = 0;

%% Synthetic atlas

% six regions and a band of brain in none, the two top planes outside the brain
[region_vol, region_names, brain] = synthetic_atlas(grid_size);
large_region = 1;
small_region = 2;
n_regions = numel(region_names);

% the bump: 16 planes (more than half the rolling median's 21) x 5 x 5 voxels,
% inside the large region
bump = false(grid_size);
bump(23:38, 10:14, 18:22) = true;
fprintf('regions: %s; large %d voxels, small %d, bump %d\n', ...
    strjoin(region_names, ', '), nnz(region_vol == large_region), ...
    nnz(region_vol == small_region), nnz(bump));

%% (a) A focal bump and a diffuse shift

% the L - R map pure noise, the L + R map with the two effects
rng(100);
vol_diff = noisy_mice(grid_size, n_mice, noise_sigma);
vol_sum = noisy_mice(grid_size, n_mice, noise_sigma);
effect = bump_amplitude * double(bump) + shift_amplitude * double(region_vol == ...
    small_region);
vol_sum(:, :, :, n_ctrl + 1:end) = vol_sum(:, :, :, n_ctrl + 1:end) + effect;
vol_diff = missing_tissue(vol_diff, brain);
vol_sum = missing_tissue(vol_sum, brain);

% the stacks and the candidate voxels, as group_differences builds them
[stacks, geom] = candidate_stacks({vol_diff, vol_sum}, brain, region_vol, n_regions, ...
    voxel_mm, perm_settings.min_mice_per_group);

t_start = tic;
perm = region_permutation_test(stacks, n_ctrl, geom, perm_settings);
fprintf('(a) ran in %.1f s\n', toc(t_start));

%% The splits

% every split of 10 mice into 5 and 5, the observed one (mice 1 to 5 control)
% among them, half of them computed
splits = perm.splits;
[n_pass, n_fail] = check(size(splits.in_ctrl, 1) == nchoosek(n_mice, n_ctrl), ...
    sprintf('%d splits, nchoosek(10, 5) = %d', size(splits.in_ctrl, 1), ...
    nchoosek(n_mice, n_ctrl)), n_pass, n_fail);
observed_row = [true(1, n_ctrl), false(1, n_mice - n_ctrl)];
[is_there, where] = ismember(observed_row, splits.in_ctrl, 'rows');
[n_pass, n_fail] = check(is_there && where == splits.observed, ...
    'the observed split is in the enumeration, where the test takes it', n_pass, n_fail);
[n_pass, n_fail] = check(size(unique(splits.in_ctrl, 'rows'), 1) == ...
    size(splits.in_ctrl, 1) && all(sum(splits.in_ctrl, 2) == n_ctrl), ...
    'the splits are distinct, each with 5 control mice', n_pass, n_fail);
[n_pass, n_fail] = check(numel(splits.computed) == nchoosek(n_mice, n_ctrl) / 2, ...
    sprintf('%d computed, the others their mirror image', numel(splits.computed)), ...
    n_pass, n_fail);

%% The observed split gives the bars' median

% the surprise median of each map as group_differences computes it, against the
% test's own on the candidate voxels, and the share of each region
for m = 1:2
    vol_m = {vol_diff, vol_sum};
    reference = bars_rolled_surprise(vol_m{m}, n_ctrl, brain, ...
        perm_settings.min_mice_per_group, perm_settings.slab_range);
    engine = perm.maps{m}.detail.rolled_unsigned;
    [n_pass, n_fail] = check(isequaln(engine, reference(geom.cand_lin)), ...
        sprintf('map %d: the observed rolled surprise equals the bars'' on all %d voxels', ...
        m, numel(engine)), n_pass, n_fail);
    reference_share = bars_share(reference, region_vol, n_regions, ...
        perm_settings.p_thresh);
    engine_share = perm.maps{m}.score(:, strcmp(perm.measure_names, 'share'));
    [n_pass, n_fail] = check(isequaln(engine_share, reference_share), ...
        sprintf('map %d: the observed share equals the bars'' in the %d regions', m, ...
        n_regions), n_pass, n_fail);
end

%% The observed split's measures, written out region by region

% the top volume in voxels, from the voxel side
topvol_k = round(perm_settings.topvol_mm3 / voxel_mm^3);
[n_pass, n_fail] = check(perm.topvol_k == topvol_k, sprintf(['top volume of %g ' ...
    'mm^3 at %g mm voxels: %d voxels'], perm_settings.topvol_mm3, voxel_mm, ...
    perm.topvol_k), n_pass, n_fail);

% the four signed measures of the observed split against each region's own,
% computed on the whole grid from the signed median of the bars' arithmetic:
% the sum, the quantile and the top volume over its voxels with a t, and the
% heaviest cluster among the bwconncomp components of its own significant
% voxels of one sign
for m = 1:2
    vol_m = {vol_diff, vol_sum};
    [~, signed_reference] = bars_rolled_surprise(vol_m{m}, n_ctrl, brain, ...
        perm_settings.min_mice_per_group, perm_settings.slab_range);
    [n_pass, n_fail] = check_measures(perm.maps{m}, splits.observed, ...
        signed_reference, region_vol, n_regions, perm_settings, topvol_k, ...
        sprintf('(a) map %d', m), n_pass, n_fail);
end

%% Clusters cut by region and sign, and their connectivity

% a map of noise with, in the experimental mice, blocks of 16 planes in the
% large region: two positive ones touching only at a corner (one cluster at
% 26-connectivity, two at 18), two touching along an edge (one at 18, two at
% 6) and a negative one touching the edge pair by a face; and a positive block
% across the border of the large region and a ventral one. The clusters of the
% observed split must be cut at the border and between the signs, and joined
% at 18-connectivity
rng(200);
vol_cut = noisy_mice(grid_size, n_mice, noise_sigma);
effect_cut = zeros(grid_size);

% the corner pair: plane 20 against plane 21, DV and ML 9 against 10
effect_cut(5:20, 5:9, 5:9) = bump_amplitude;
effect_cut(21:36, 10:14, 10:14) = bump_amplitude;

% the edge pair, the same planes: DV 8 against 9 and ML 28 against 29; the
% negative block against the first by a face, ML 24 against 25
effect_cut(5:20, 5:8, 25:28) = bump_amplitude;
effect_cut(5:20, 9:12, 29:32) = bump_amplitude;
effect_cut(5:20, 5:8, 21:24) = -bump_amplitude;

% across the border of the large region (DV 3 to 22) and the ventral one below
effect_cut(41:56, 19:26, 5:8) = bump_amplitude;
vol_cut(:, :, :, n_ctrl + 1:end) = vol_cut(:, :, :, n_ctrl + 1:end) + effect_cut;
vol_cut = missing_tissue(vol_cut, brain);
vol_cut_noise = missing_tissue(noisy_mice(grid_size, n_mice, noise_sigma), brain);
[stacks_cut, geom_cut] = candidate_stacks({vol_cut, vol_cut_noise}, brain, ...
    region_vol, n_regions, voxel_mm, perm_settings.min_mice_per_group);
observed_only = perm_settings;
observed_only.n_permutations = 1;
observed_only.keep_cluster_voxels = true;
observed_only.keep_t = true;
perm_cut = region_permutation_test(stacks_cut, n_ctrl, geom_cut, observed_only);

% the clusters of all the significant voxels that hold two regions or both signs
[~, signed_cut, t_cut] = bars_rolled_surprise(vol_cut, n_ctrl, brain, ...
    perm_settings.min_mice_per_group, perm_settings.slab_range);
n_crossing = crossing_clusters(signed_cut, region_vol, perm_settings);
[n_pass, n_fail] = check(n_crossing > 0, sprintf(['(c) %d clusters of the ' ...
    'significant voxels hold two regions or both signs'], n_crossing), n_pass, n_fail);

% the measures against each region's own, at 18-connectivity
[n_pass, n_fail] = check_measures(perm_cut.maps{1}, perm_cut.splits.observed, ...
    signed_cut, region_vol, n_regions, perm_settings, topvol_k, '(c) map 1', ...
    n_pass, n_fail);

% the voxels kept of each heaviest cluster, which per_mouse_region_values reads
[n_pass, n_fail] = check_cluster_voxels(perm_cut.maps{1}, perm_cut.splits.observed, ...
    strcmp(perm_cut.measure_names, 'cluster'), signed_cut, region_vol, n_regions, ...
    '(c) map 1', n_pass, n_fail);

% the t kept on request, which scale_free_test draws
[n_pass, n_fail] = check(isequaln(perm_cut.maps{1}.detail.t, t_cut(geom_cut.cand_lin)), ...
    sprintf('(c) map 1: the t kept on request equals the t written out on all %d voxels', ...
    numel(geom_cut.cand_lin)), n_pass, n_fail);

% the heaviest positive cluster of the large region at 6, 18 and 26: three
% different masses, so the check above tells 18 from the other two
masses = zeros(1, 3);
connectivities = [6 18 26];
for c = 1:3
    connectivity_settings = perm_settings;
    connectivity_settings.cluster_connectivity = connectivities(c);
    reference_pos = reference_measures(signed_cut, region_vol, n_regions, ...
        connectivity_settings, topvol_k);
    masses(c) = reference_pos(large_region, 4);
end
[n_pass, n_fail] = check(numel(unique(masses)) == 3, sprintf(['(c) heaviest ' ...
    'positive cluster of the large region at 6, 18, 26-connectivity: %.0f, %.0f, ' ...
    '%.0f; the test''s %.0f'], masses, perm_cut.maps{1}.null_pos(1, large_region, ...
    strcmp(perm_cut.measure_names, 'cluster'))), n_pass, n_fail);

%% A mirror image is the split's negative

% the groups swapped: the observed split of the swapped mice is the mirror of
% the observed one, computed directly
stacks_swapped = cellfun(@(s) s(:, [n_ctrl + 1:n_mice, 1:n_ctrl]), stacks, ...
    'UniformOutput', false);
one_split = perm_settings;
one_split.n_permutations = 1;
perm_swapped = region_permutation_test(stacks_swapped, n_ctrl, geom, one_split);
mirror_row = find(ismember(splits.in_ctrl, ~observed_row, 'rows'));
same = true;
for m = 1:2
    same = same && isequaln(perm_swapped.maps{m}.null_pos(1, :, :), ...
        perm.maps{m}.null_pos(mirror_row, :, :)) && ...
        isequaln(perm_swapped.maps{m}.null_neg(1, :, :), ...
        perm.maps{m}.null_neg(mirror_row, :, :));
end
[n_pass, n_fail] = check(same, ['the mirror of the observed split, computed, equals ' ...
    'the one taken from it'], n_pass, n_fail);

%% The large region first under cluster and top volume

% the L + R map: every measure's rank of the two regions, score, p and corrected p
map_sum = perm.maps{2};
fprintf('\n(a) L + R map, large region (bump) and small region (shift):\n');
fprintf('  %-8s %22s %22s\n', 'measure', 'large: rank score pfwer', ...
    'small: rank score pfwer');
for k = 1:numel(perm.measure_names)
    rank_k = region_ranks(map_sum.score(:, k));
    fprintf('  %-8s %5d %9.3g %6.3f   %5d %9.3g %6.3f   (p %.3f, %.3f)\n', ...
        perm.measure_names{k}, rank_k(large_region), map_sum.score(large_region, k), ...
        map_sum.p_fwer(large_region, k), rank_k(small_region), ...
        map_sum.score(small_region, k), map_sum.p_fwer(small_region, k), ...
        map_sum.p_perm(large_region, k), map_sum.p_perm(small_region, k));
end
fprintf('  heaviest positive cluster of the large region: %d voxels, peak %.2f\n', ...
    map_sum.detail.cluster_n(large_region, 1), ...
    map_sum.detail.cluster_peak(large_region, 1));
for name = {'cluster', 'topvol'}
    k = find(strcmp(perm.measure_names, name{1}));
    rank_k = region_ranks(map_sum.score(:, k));
    [n_pass, n_fail] = check(rank_k(large_region) == 1 && ...
        map_sum.p_fwer(large_region, k) < 0.05 && map_sum.score(large_region, k) > 0, ...
        sprintf('%s: the large region first, positive, corrected p %.3f', name{1}, ...
        map_sum.p_fwer(large_region, k)), n_pass, n_fail);
end

% the L - R map is noise: no corrected p < 0.05 there
map_diff = perm.maps{1};
[n_pass, n_fail] = check(~any(map_diff.p_fwer(:) < 0.05), ...
    sprintf('(a) L - R map, pure noise: smallest corrected p %.3f', ...
    min(map_diff.p_fwer(:))), n_pass, n_fail);

%% (b) Pure noise

% every seed, both maps: the families (a map and a measure) with a region at a
% corrected p < 0.05
n_measures = numel(perm.measure_names);
n_rejected = zeros(1, n_measures);
n_families = 0;
for seed = null_seeds
    rng(seed);
    vol_diff = missing_tissue(noisy_mice(grid_size, n_mice, noise_sigma), brain);
    vol_sum = missing_tissue(noisy_mice(grid_size, n_mice, noise_sigma), brain);
    [stacks, geom] = candidate_stacks({vol_diff, vol_sum}, brain, region_vol, n_regions, ...
        voxel_mm, perm_settings.min_mice_per_group);
    perm_null = region_permutation_test(stacks, n_ctrl, geom, perm_settings);
    for m = 1:2
        n_rejected = n_rejected + any(perm_null.maps{m}.p_fwer < 0.05, 1);
        n_families = n_families + 1;
    end
end

% at 5% per family, 4 or more of 12 has a chance under 0.003
fprintf('\n(b) pure noise, %d families per measure (seeds %s, two maps):\n', ...
    n_families, mat2str(null_seeds));
for k = 1:n_measures
    fprintf('  %-8s %d with a corrected p < 0.05\n', perm.measure_names{k}, ...
        n_rejected(k));
end
[n_pass, n_fail] = check(all(n_rejected <= 3), sprintf(['(b) no measure above ' ...
    'chance: at most %d of %d families'], max(n_rejected), n_families), n_pass, n_fail);

%% Result

delete(gcp('nocreate'));
fprintf('\ntest_region_permutation: %d passed, %d failed\n', n_pass, n_fail);
if n_fail > 0
    error('test_region_permutation: %d checks failed (listed above).', n_fail);
end

% ===== Local functions =====

function [region_vol, region_names, brain] = synthetic_atlas(grid_size)
% Six regions on the grid: a large one (the dorsal half), a small one inside the
% ventral half, and four of the ventral half around it; a band of brain in no
% region, and the two top planes of DV outside the brain.

region_vol = zeros(grid_size, 'uint8');
region_vol(:, 3:22, :) = 1;
region_vol(1:30, 23:40, 1:20) = 3;
region_vol(31:60, 23:40, 1:20) = 4;
region_vol(1:30, 23:40, 21:38) = 5;
region_vol(31:60, 23:40, 21:38) = 6;
region_vol(21:40, 26:31, 6:11) = 2;
region_names = {'large', 'small', 'ventral 1', 'ventral 2', 'ventral 3', 'ventral 4'};
brain = true(grid_size);
brain(:, 1:2, :) = false;
end

function vol = noisy_mice(grid_size, n_mice, noise_sigma)
% Each mouse's map: Gaussian noise smoothed in 3D, scaled to SD 1, around 10.

vol = zeros([grid_size n_mice], 'single');
for i = 1:n_mice
    noise = imgaussfilt3(randn(grid_size), noise_sigma);
    vol(:, :, :, i) = single(10 + noise / std(noise(:)));
end
end

function vol = missing_tissue(vol, brain)
% NaN outside the brain, and where some mice have no tissue: two mice at either
% end of the brain (4 left in their group), three control mice at one corner
% (2 left: no t there).

vol(repmat(~brain, [1 1 1 size(vol, 4)])) = NaN;
vol(1:4, :, :, 3) = NaN;
vol(57:60, :, :, 7) = NaN;
vol(1:10, :, 1:3, 1:3) = NaN;
end

function [stacks, geom] = candidate_stacks(vols, brain, region_vol, n_regions, ...
    voxel_mm, min_mice_per_group)
% The stacks of the maps on the candidate voxels (in the brain, with at least
% twice min_mice_per_group mice with a value), and their geometry.

n_mice = size(vols{1}, 4);
n_with_value = sum(~isnan(vols{1}), 4);
candidates = brain & n_with_value >= 2 * min_mice_per_group;
geom = struct();
geom.grid_size = size(brain);
geom.cand_lin = uint32(find(candidates));
geom.cand_region = region_vol(candidates);
geom.n_regions = n_regions;
geom.voxel_mm = voxel_mm;
stacks = cell(1, numel(vols));
for m = 1:numel(vols)
    stacks{m} = zeros(numel(geom.cand_lin), n_mice, 'single');
    for i = 1:n_mice
        one_mouse = vols{m}(:, :, :, i);
        stacks{m}(:, i) = abs(one_mouse(candidates));
    end
end
end

function [rolled, rolled_signed, t] = bars_rolled_surprise(vol, n_ctrl, brain, ...
    min_mice_per_group, slab_range)
% The surprise median the bars take, written out as group_differences computes
% it (group_welch_t, welch_surprise, rolling_surprise_median) on the 4D maps;
% the same median of the surprise signed by the t, which the region measures
% other than the share take; and the t, NaN without one.

x_ctrl = abs(vol(:, :, :, 1:n_ctrl));
x_exp = abs(vol(:, :, :, n_ctrl + 1:end));
n_vox_ctrl = single(sum(~isnan(x_ctrl), 4));
n_vox_exp = single(sum(~isnan(x_exp), 4));
has_t = brain & any(~isnan(x_ctrl), 4) & any(~isnan(x_exp), 4) & ...
    n_vox_ctrl >= min_mice_per_group & n_vox_exp >= min_mice_per_group;
avg_ctrl = nanmean(x_ctrl, 4); %#ok<*NANMEAN>
avg_exp = nanmean(x_exp, 4);
sem_ctrl = nanstd(x_ctrl, [], 4) ./ sqrt(double(n_vox_ctrl)); %#ok<*NANSTD>
sem_exp = nanstd(x_exp, [], 4) ./ sqrt(double(n_vox_exp));
sem_ctrl(sem_ctrl == 0) = NaN;
sem_exp(sem_exp == 0) = NaN;
sem_diff = sqrt(sem_ctrl.^2 + sem_exp.^2);
sem_diff(isnan(sem_diff) | sem_diff == 0) = NaN;
t = (avg_exp - avg_ctrl) ./ sem_diff;
t(~has_t) = NaN;
var1 = sem_ctrl.^2;
var2 = sem_exp.^2;
df = (var1 + var2).^2 ./ ...
    (var1.^2 ./ (double(n_vox_ctrl) - 1) + var2.^2 ./ (double(n_vox_exp) - 1));
df(n_vox_ctrl < min_mice_per_group | n_vox_exp < min_mice_per_group) = NaN;
surprise = -log10(2 * tcdf(-abs(t), df));

% the median over +/- slab_range planes, on the voxels with a t, of the
% surprise and of the signed surprise
rolled = slab_rolling_median(surprise, has_t, slab_range);
rolled_signed = slab_rolling_median(surprise .* sign(t), has_t, slab_range);
end

function rolled = slab_rolling_median(surprise, has_t, slab_range)
% The median over +/- slab_range planes, on the voxels with a t, as
% rolling_surprise_median of group_differences.

rolled = zeros(size(surprise), 'single');
n_slices = size(surprise, 1);
for z = 1:n_slices
    z_start = max(1, z - slab_range);
    z_end = min(n_slices, z + slab_range);
    slab_data = surprise(z_start:z_end, :, :);
    slab_data(~has_t(z_start:z_end, :, :)) = NaN;
    slab_median = nanmedian(slab_data, 1); %#ok<NANMEDIAN>
    slab_median(~has_t(z, :, :)) = NaN;
    rolled(z, :, :) = slab_median;
end
end

function [pos, neg] = reference_measures(rolled, region_vol, n_regions, ...
    perm_settings, topvol_k)
% Each region's summed surprise at p < p_thresh, quantile, top volume and
% heaviest cluster at p < cluster_p, positive and negative apart (n_regions x 4,
% magnitudes; NaN without a voxel with a t), written out region by region on
% the whole grid: the clusters are the bwconncomp components of the region's
% own significant voxels of one sign.

thresh_sig = -log10(perm_settings.p_thresh);
thresh_cluster = -log10(perm_settings.cluster_p);
pos = nan(n_regions, 4);
neg = nan(n_regions, 4);
for r = 1:n_regions
    in_region = region_vol == r;
    values = rolled(in_region & ~isnan(rolled));
    if isempty(values)
        continue
    end
    for direction = [1 -1]

        % the effects of this sign, as magnitudes
        signed_values = direction * values;
        part = max(signed_values, 0);
        sorted_part = sort(double(part), 'descend');
        measures = zeros(1, 4);
        measures(1) = sum(double(signed_values(signed_values > thresh_sig)));
        measures(2) = quantile(part, perm_settings.region_quantile);
        measures(3) = mean(sorted_part(1:min(topvol_k, numel(sorted_part))));

        % the heaviest cluster of the region's significant voxels of this sign
        components = bwconncomp(in_region & direction * rolled > thresh_cluster, ...
            perm_settings.cluster_connectivity);
        masses = cellfun(@(voxels) sum(abs(double(rolled(voxels)))), ...
            components.PixelIdxList);
        if ~isempty(masses)
            measures(4) = max(masses);
        end
        if direction > 0
            pos(r, :) = measures;
        else
            neg(r, :) = measures;
        end
    end
end
end

function n_crossing = crossing_clusters(rolled, region_vol, perm_settings)
% The bwconncomp components of all the significant voxels in a region (at
% p < cluster_p, either sign) that hold voxels of two regions or of both signs.

is_sig = abs(rolled) > -log10(perm_settings.cluster_p) & region_vol > 0;
label = 2 * double(region_vol) - (rolled > 0);
components = bwconncomp(is_sig, perm_settings.cluster_connectivity);
n_labels = cellfun(@(voxels) numel(unique(label(voxels))), components.PixelIdxList);
n_crossing = nnz(n_labels > 1);
end

function [n_pass, n_fail] = check_measures(map, observed_row, signed_rolled, ...
    region_vol, n_regions, perm_settings, topvol_k, label, n_pass, n_fail)
% The observed split's sum, quantile, top volume and heaviest cluster of every
% region, both signs, against those written out by reference_measures.

[reference_pos, reference_neg] = reference_measures(signed_rolled, region_vol, ...
    n_regions, perm_settings, topvol_k);
engine_pos = reshape(map.null_pos(observed_row, :, 2:5), n_regions, 4);
engine_neg = reshape(map.null_neg(observed_row, :, 2:5), n_regions, 4);
largest_error = max(relative_error([engine_pos; engine_neg], ...
    [reference_pos; reference_neg]));
[n_pass, n_fail] = check(largest_error < 1e-6, sprintf(['%s: the observed sum, ' ...
    'quantile, top volume and heaviest cluster equal each region''s own, both ' ...
    'signs (largest relative difference %.1g)'], label, largest_error), n_pass, ...
    n_fail);
end

function [n_pass, n_fail] = check_cluster_voxels(map, observed_row, is_cluster, ...
    signed_rolled, region_vol, n_regions, label, n_pass, n_fail)
% The voxels the test kept of every region's heaviest cluster of each sign: as
% many as the cluster has, all in the region and of the cluster's sign, their
% summed |surprise| its mass; none for a region without a cluster.

all_match = true;
n_clusters = 0;
for r = 1:n_regions
    for s = 1:2

        % the kept voxels, and the cluster's mass and sign
        voxels = map.detail.cluster_voxels{r, s};
        if s == 1
            mass = map.null_pos(observed_row, r, is_cluster);
            direction = 1;
        else
            mass = map.null_neg(observed_row, r, is_cluster);
            direction = -1;
        end

        if isempty(voxels)
            is_match = map.detail.cluster_n(r, s) == 0;
        else
            n_clusters = n_clusters + 1;
            values = direction * double(signed_rolled(voxels));
            is_match = numel(voxels) == map.detail.cluster_n(r, s) && ...
                all(region_vol(voxels) == r) && all(values > 0) && ...
                abs(sum(values) - mass) <= 1e-6 * max(1, mass);
        end
        all_match = all_match && is_match;
    end
end
[n_pass, n_fail] = check(all_match && n_clusters > 0, sprintf(['%s: the kept ' ...
    'voxels of the %d heaviest clusters are theirs (count, region, sign, mass)'], ...
    label, n_clusters), n_pass, n_fail);
end

function err = relative_error(engine, reference)
% The difference of two arrays relative to the reference (absolute below 1);
% Inf where only one of them is NaN.

err = abs(engine - reference) ./ max(1, abs(reference));
err(isnan(engine) & isnan(reference)) = 0;
err(isnan(engine) ~= isnan(reference)) = Inf;
err = err(:);
end

function share = bars_share(rolled, region_vol, n_regions, p_thresh)
% Each region's fraction of voxels with a t at p < p_thresh, as the bars.

share = nan(n_regions, 1);
for r = 1:n_regions
    values = rolled(region_vol == r);
    share(r) = nnz(values > -log10(p_thresh)) / nnz(~isnan(values));
end
end

function ranks = region_ranks(score)
% Each region's rank by |score|, largest first; NaN last.

magnitude = abs(score);
magnitude(isnan(magnitude)) = -Inf;
[~, order] = sort(magnitude, 'descend');
ranks = zeros(size(score));
ranks(order) = 1:numel(score);
end

function [n_pass, n_fail] = check(cond, label, n_pass, n_fail)
% Print PASS or FAIL for one check, and add it to the tally.

if cond
    n_pass = n_pass + 1;
    fprintf('  PASS  %s\n', label);
else
    n_fail = n_fail + 1;
    fprintf('  FAIL  %s\n', label);
end
end
