function scale_free_test(run_settings)
%SCALE_FREE_TEST  Step 3's region test on the raw asymmetry index, which no scale changes.
%   SCALE_FREE_TEST(run_settings) does the work of run_scale_free_test, which
%   sets the fields of run_settings and says what each one does.
%
%   Reads, from the control and the experimental group's folders, the
%   collected stack <channel>_4d.mat (run_collect_by_group) and the background
%   masks <channel>_4d_normalized_bkgmask.mat (run_normalise_groups), one mouse
%   at a time; from comp_out_dir, step 3's table
%   (Region_Surprise_DiffSum_<comp_tag>.csv), the leave-one-out's clusters of
%   run_per_mouse_values (Per_Mouse_LOO_<tag>.mat), whose fold 0, all the
%   mice, is step 3's cluster in the region named in advance, the
%   selection-matched test's splits of run_per_mouse_values
%   (Per_Mouse_Selection_<tag>.mat), the same cluster's mass under every split
%   of step 3's region test, and, to check the reading against, the mice's
%   maps of run_per_mouse_values (Per_Mouse_Maps_<tag>.mat) when they are
%   there. Writes into comp_out_dir, with <tag> = <comp_tag>_smooth<sigma>:
%     Scale_Free_Regions_<tag>.csv   one row per region of the bars: each
%                                    sign's heaviest cluster (mass, voxels,
%                                    peak, planes) with its own p; the score,
%                                    the heavier of the two with its sign, its
%                                    p and corrected p
%     Scale_Free_Clusters_<tag>.csv  the region named in advance: its heaviest
%                                    cluster of each sign and step 3's cluster,
%                                    where each sits and how they overlap; for
%                                    each, the one-sided p of its sign and,
%                                    for the heavier sign, the p of either sign
%                                    (a sign without a cluster has mass 0, and
%                                    every split reaches it: p 1)
%     Scale_Free_Bars_<tag>          the bars of the cluster mass, .fig and .png
%     Scale_Free_TMap_<tag>          the t map of the index in the region named
%                                    in advance, at each cluster's planes and
%                                    from above, .fig and .png
%     Scale_Free_<tag>.mat           the cache of the test
%
%   The reading. For each mouse and voxel, the asymmetry index
%       AI = |L - R| / (L + R)
%   of its collected stack less its off-tissue level (off_tissue_level, the raw
%   stack of run_per_mouse_values), NaN outside its tissue and smoothed as step
%   3 smooths, before the folding (tissue_only, compute_lr_stats); NaN where
%   L + R is not above 0, or where either side is outside the tissue. No
%   normalisation (step 2) and no alignment of the groups (step 3): a factor on
%   a mouse or on a group divides out of the ratio, so no slope can move it.
%   Its zero is each mouse's own off-tissue level.
%
%   The statistics are step 3's region test on this one map: the candidate
%   voxels and the mice's values on them (permutation_stacks), then, under
%   every split of the mice into groups of their sizes (252 for 5 against 5,
%   126 for 5 against 4), the Welch t of experimental minus control where each
%   group has min_mice_per_group mice with a value, its surprise -log10 p
%   signed by the t, its median over +/- slab_range planes, the voxels at
%   p < cluster_p, their clusters connected (cluster_connectivity) within each
%   region's own voxels of the left hemisphere, and each region's heaviest
%   cluster of each sign (region_permutation_test). A region's score is the
%   heavier of its two clusters' masses with its sign, experimental higher
%   positive; its p, the share of the splits whose score's magnitude reaches
%   it, the observed split included; its corrected p, the share whose largest
%   over the regions does. Each sign's own p is the share of the splits whose
%   heaviest cluster of that sign is as heavy. The other measures of
%   region_permutation_test come with it and are not reported: the cluster
%   mass was fixed before the test.
%
%   The test named in advance (a_priori_region, direction_named): after RWS,
%   the mass of the region's heaviest cluster where the index is higher in the
%   experimental group, against the same in every split (one-sided: RWS
%   potentiates the stimulated barrels' synapses, Gambino et al. 2014), the p
%   of either sign beside it; after behavior, the p of either sign first, the
%   direction carried over from RWS beside it.

% settings of run_scale_free_test, under the names the code below uses
paths = run_settings.paths;
ctrl_type = run_settings.ctrl_type;
exp_type = run_settings.exp_type;
behavior_mice = run_settings.behavior_mice;
a_priori_region = run_settings.a_priori_region;
direction_named = run_settings.direction_named;
apply_smoothing = run_settings.apply_smoothing;
smooth_sigma = run_settings.smooth_sigma;
min_mice_per_group = run_settings.min_mice_per_group;
slab_range = run_settings.slab_range;
cluster_p = run_settings.cluster_p;
cluster_connectivity = run_settings.cluster_connectivity;
n_permutations = run_settings.n_permutations;
permutation_workers = run_settings.permutation_workers;
t_limit = run_settings.t_limit;
force_recompute = run_settings.force_recompute;
channel = run_settings.channel;
comp_tag = run_settings.comp_tag;
ctrl_dir = run_settings.ctrl_dir;
exp_dir = run_settings.exp_dir;
comp_out_dir = run_settings.comp_out_dir;

% the experimental groups of run_group_differences
if ~ismember(exp_type, {'rws', 'behavior'})
    error('run_scale_free_test: unknown exp_type ''%s'' (use ''rws'' or ''behavior'').', ...
        exp_type);
end

% the file names carry the smoothing, as run_per_mouse_values'
if apply_smoothing
    smooth_suffix = sprintf('_smooth%g', smooth_sigma);
else
    smooth_suffix = '_nosmooth';
end
file_tag = [comp_tag smooth_suffix];

% the region test as step 3 runs it, its threshold of the share at the clusters'
% p; the top volume and the quantile at step 3's values, though only the cluster
% is read; the observed split's clusters and t kept, for the tables and the t map
perm_settings = struct();
perm_settings.min_mice_per_group = min_mice_per_group;
perm_settings.slab_range = slab_range;
perm_settings.p_thresh = cluster_p;
perm_settings.cluster_p = cluster_p;
perm_settings.cluster_connectivity = cluster_connectivity;
perm_settings.topvol_mm3 = 0.1;
perm_settings.region_quantile = 0.99;
perm_settings.n_permutations = n_permutations;
perm_settings.n_workers = permutation_workers;
perm_settings.seed = 0;
perm_settings.keep_cluster_voxels = true;
perm_settings.keep_t = true;

% what is read besides the stacks
step3_table = fullfile(comp_out_dir, ['Region_Surprise_DiffSum_' comp_tag '.csv']);
loo_file = fullfile(comp_out_dir, ['Per_Mouse_LOO_' file_tag '.mat']);
selection_file = fullfile(comp_out_dir, ['Per_Mouse_Selection_' file_tag '.mat']);
maps_file = fullfile(comp_out_dir, ['Per_Mouse_Maps_' file_tag '.mat']);

%% Atlas and regions

% the 10 um annotation on the volumes' crop, the regions of the bars, the region
% named in advance, and the box around it (its planes and slab_range planes on
% each side), the box of run_per_mouse_values' maps
A = get_atlas_crop('ccf');
brainMask = A.brainMask;
ccf_first_index = A.aplims(1) - 1;
[T_regions, valid_pixels, region_of_voxel] = surprise_regions(A.annot, paths.atlas);
clear A
T_regions.a_priori = strcmp(T_regions.acronym, a_priori_region);
if ~any(T_regions.a_priori)
    error(['run_scale_free_test: a_priori_region is %s, not one of the bars'' ' ...
           'regions. Use their atlas acronyms: %s.'], a_priori_region, ...
           strjoin(T_regions.acronym, ', '));
end
apriori_row = find(T_regions.a_priori);
masks = region_box_masks(T_regions, valid_pixels, region_of_voxel, {a_priori_region}, ...
    a_priori_region, slab_range);

% step 3's cluster in the region named in advance, from the leave-one-out's
% fold 0, its mass and p from the selection-matched test's splits, checked
% against step 3's table
step3 = step3_cluster(loo_file, selection_file, step3_table, a_priori_region, masks);

%% The test, from its cache or computed

% every control mouse saved; for behavior, the mice behavior_mice names, as in
% run_group_differences
if strcmp(exp_type, 'behavior')
    exp_named = behavior_mice;
else
    exp_named = {};
end

% what the test depends on: the stacks it read and the settings, checked when
% the cache is read
cache_settings = struct('exp_named', {exp_named}, 'apply_smoothing', apply_smoothing, ...
    'smooth_sigma', smooth_sigma, 'min_mice_per_group', min_mice_per_group, ...
    'slab_range', slab_range, 'cluster_p', cluster_p, 'cluster_connectivity', ...
    cluster_connectivity, 'n_permutations', n_permutations, 'a_priori_region', ...
    a_priori_region, 'stack_files', stack_record({ctrl_dir, exp_dir}, channel));
cache_file = fullfile(comp_out_dir, ['Scale_Free_' file_tag '.mat']);
result = [];
if exist(cache_file, 'file') && ~force_recompute
    S_cache = load(cache_file, 'result', 'cache_settings');
    if isequal(S_cache.cache_settings, cache_settings)
        fprintf('Loading the test from %s...\n', cache_file);
        result = S_cache.result;
    else
        fprintf(['The test in %s was made from other stacks or settings; it is run ' ...
                 'again.\n'], cache_file);
    end
    clear S_cache
end

if isempty(result)

    % each mouse's asymmetry index, control mice first, in step 3's order
    t_maps = tic;
    ctrl_mice = group_ai(ctrl_type, ctrl_dir, channel, {}, brainMask, ...
        apply_smoothing, smooth_sigma);
    exp_mice = group_ai(exp_type, exp_dir, channel, exp_named, brainMask, ...
        apply_smoothing, smooth_sigma);
    maps_min = toc(t_maps) / 60;

    % the reading must be run_per_mouse_values' raw stack, mouse for mouse
    check_against_per_mouse(ctrl_mice, exp_mice, maps_file, masks);

    % the candidate voxels and every mouse's index on them, as step 3 builds its
    % stacks; the maps are not needed after that
    [stacks, geom] = permutation_stacks({{ctrl_mice.ai, exp_mice.ai}}, valid_pixels, ...
        region_of_voxel, height(T_regions), min_mice_per_group);
    ctrl_mice = rmfield(ctrl_mice, 'ai');
    exp_mice = rmfield(exp_mice, 'ai');

    % the region test under every split
    fprintf('Permutation test of the regions on the asymmetry index...\n');
    t_test = tic;
    perm = region_permutation_test(stacks, numel(ctrl_mice.names), geom, perm_settings);
    clear stacks
    test_min = toc(t_test) / 60;

    % what the tables and figures read, without the values per voxel
    result = test_result(perm, geom, masks, apriori_row, ctrl_mice, exp_mice);
    result.maps_min = maps_min;
    result.test_min = test_min;
    clear perm geom
    fprintf('Saving the test to %s...\n', cache_file);
    save(cache_file, 'result', 'cache_settings', '-v7.3');
end
clear valid_pixels region_of_voxel brainMask

%% Tables

% every region, both signs; the region named in advance, its clusters and
% step 3's
T_out = region_table(result, T_regions, ctrl_type, exp_type, ccf_first_index);
T_clusters = cluster_table(result, step3, masks, apriori_row, ctrl_type, exp_type, ...
    ccf_first_index);
writetable(T_out, fullfile(comp_out_dir, ['Scale_Free_Regions_' file_tag '.csv']));
writetable(T_clusters, fullfile(comp_out_dir, ['Scale_Free_Clusters_' file_tag '.csv']));
named = named_test(result, apriori_row, a_priori_region, direction_named, ctrl_type, ...
    exp_type);
print_summary(result, T_out, T_clusters, named, comparison_line(ctrl_type, exp_type, ...
    result));

%% Figures

% the bars of the cluster mass in step 3's style, the test named in advance
% under the title; beside the region named in advance, the p of that test (with
% its direction named, the one-sided p, not the p of either sign)
if direction_named
    a_priori_test = struct('p', nan(height(T_regions), 1), 'label', ...
        sprintf('%s higher, one-sided', exp_type));
    a_priori_test.p(apriori_row) = named.p_one;
else
    a_priori_test = [];
end
fig_bars = plot_measure_bars(result.perm, T_regions, 'cluster', exp_type, ...
    ['scale_free_' file_tag], {'asymmetry index |L - R| / (L + R), raw stack'}, ...
    named.lines, a_priori_test);
saveas(fig_bars, fullfile(comp_out_dir, ['Scale_Free_Bars_' file_tag '.fig']));
exportgraphics(fig_bars, fullfile(comp_out_dir, ['Scale_Free_Bars_' file_tag '.png']), ...
    'Resolution', 300);

% the t map of the region named in advance, at each cluster's planes
comparison = struct('ctrl_type', ctrl_type, 'exp_type', exp_type, 'region', ...
    a_priori_region, 'ccf_first_index', ccf_first_index, 't_limit', t_limit, ...
    'cluster_p', cluster_p, 'slab_range', slab_range);
plot_t_map(result, step3, masks, comparison, file_tag, comp_out_dir);
fprintf('Scale-free test saved to: %s\n', comp_out_dir);

end

% ===== Local functions: the reading =====

function records = stack_record(group_dirs, channel)
% The stacks the test reads, each with its size in bytes and its date, so a
% cache made from other stacks is not read.

records = struct('file', {}, 'bytes', {}, 'datenum', {});
for g = 1:numel(group_dirs)
    for name = {[channel '_4d.mat'], [channel '_4d_normalized_bkgmask.mat']}
        file = fullfile(group_dirs{g}, name{1});
        listing = dir(file);
        if isempty(listing)
            error('run_scale_free_test: %s not found. Run steps 1 and 2 first.', file);
        end
        records(end + 1) = struct('file', file, 'bytes', listing.bytes, ...
            'datenum', listing.datenum); %#ok<AGROW>
    end
end
end

function G = group_ai(group, group_dir, channel, named_mice, brainMask, ...
    apply_smoothing, smooth_sigma)
% One group's mice, one at a time from the collected stack: the off-tissue level
% subtracted, NaN outside the tissue and smoothed, folded, and the asymmetry
% index of every voxel (AP x DV x ML of the left hemisphere x mice). The mice
% are those named, or every mouse of the background masks, in their order,
% which is step 3's.

% the collected stack and the background masks, read one mouse at a time
raw_var_name = [channel '_4d'];
raw_file = fullfile(group_dir, [channel '_4d.mat']);
mask_file = fullfile(group_dir, [channel '_4d_normalized_bkgmask.mat']);
M_raw = matfile(raw_file);
M_mask = matfile(mask_file);
S_mask = load(mask_file, 'current_mice');

% the mice, by their place in the masks and in the collected stack
if isempty(named_mice)
    names = S_mask.current_mice;
    mask_idx = 1:numel(names);
else
    [is_saved, mask_idx] = ismember(named_mice, S_mask.current_mice);
    if ~all(is_saved)
        error('run_scale_free_test: %s not among the mice of %s (%s).', ...
            strjoin(named_mice(~is_saved), ', '), mask_file, ...
            strjoin(S_mask.current_mice, ', '));
    end
    names = named_mice;
end
raw_idx = collected_index(names, group, M_raw, raw_var_name, raw_file);

n_mice = numel(names);
raw_size = size(M_raw, raw_var_name);
G = struct();
G.group = group;
G.names = names;
G.background = nan(n_mice, 1);
G.ai = nan([raw_size(1:2), floor(raw_size(3) / 2), n_mice], 'single');
for k = 1:n_mice
    t_mouse = tic;
    fprintf('%s, mouse %d of %d: %s\n', group, k, n_mice, names{k});
    raw = M_raw.(raw_var_name)(:, :, :, raw_idx(k));
    bg_mask = M_mask.recomputed_bkg_mask_4d(:, :, :, mask_idx(k));

    % less the off-tissue level, NaN outside the tissue and smoothed, as
    % run_per_mouse_values' raw stack; a raw 0 is a voxel no section reached
    G.background(k) = off_tissue_level(raw, bg_mask, brainMask);
    fprintf('  off-tissue level %.1f\n', G.background(k));
    raw(raw == 0) = NaN;
    raw = raw - G.background(k);
    raw = tissue_only(raw, bg_mask, brainMask, apply_smoothing, smooth_sigma);
    clear bg_mask

    % folded, and the index where L + R is above 0
    [lr_diff, lr_sum] = compute_lr_stats(raw);
    clear raw
    ai = abs(lr_diff) ./ lr_sum;
    ai(~(lr_sum > 0)) = NaN;
    fprintf('  %d voxels with an index, %d with L + R not above 0\n', nnz(~isnan(ai)), ...
        nnz(~isnan(lr_sum) & ~(lr_sum > 0)));
    G.ai(:, :, :, k) = ai;
    clear lr_diff lr_sum ai

    fprintf('  done in %.1f min.\n', toc(t_mouse) / 60);
end
end

function check_against_per_mouse(ctrl_mice, exp_mice, maps_file, masks)
% Stops unless every mouse's index in the box around the region named in
% advance is the index of run_per_mouse_values' raw stack there, to the bit,
% and its off-tissue level the same: the same mice, read and smoothed alike.
% Says so and goes on when those maps are not there.

if ~exist(maps_file, 'file')
    fprintf('  (%s not there: the index is not checked against it)\n', maps_file);
    return
end
S_maps = load(maps_file, 'ctrl_mice', 'exp_mice');
groups = {ctrl_mice, S_maps.ctrl_mice; exp_mice, S_maps.exp_mice};
for g = 1:2
    G = groups{g, 1};
    cached = groups{g, 2};
    if ~isequal(G.names(:), cached.names(:)) || ...
            ~isequal(G.background(:), cached.background(:))
        error(['run_scale_free_test: the %s mice or their off-tissue levels are not ' ...
               'those of %s.'], G.group, maps_file);
    end
    for k = 1:numel(G.names)

        % the cached box holds the region's columns and their mirror images, so
        % compute_lr_stats folds it as it folds the whole width
        [lr_diff, lr_sum] = compute_lr_stats(cached.box_raw{k});
        ai_cached = abs(lr_diff) ./ lr_sum;
        ai_cached(~(lr_sum > 0)) = NaN;
        ai_box = G.ai(masks.box_ap, masks.box_dv, masks.box_ml_left, k);
        if ~isequaln(ai_box, ai_cached)
            error(['run_scale_free_test: the index of %s in the box around the ' ...
                   'region is not that of run_per_mouse_values'' raw stack (%s).'], ...
                   G.names{k}, maps_file);
        end
    end
end
fprintf(['  every mouse''s index in the box around the region is that of ' ...
         'run_per_mouse_values'' raw stack, to the bit\n']);
end

% ===== Local functions: step 3's cluster and the test's results =====

function step3 = step3_cluster(loo_file, selection_file, step3_table, region, masks)
% Step 3's cluster in the region named in advance, as linear indices into the
% folded grid: the leave-one-out's fold 0 (all the mice) of run_per_mouse_values,
% which keeps it voxel for voxel (step 3 keeps no voxels). Its mass and its p,
% of either sign and of the experimental group higher alone, from the
% selection-matched test's splits (the splits of step 3's region test, its
% first the observed one, whose cluster must be fold 0's); its voxel count,
% mass and p of either sign checked against step 3's table.

if ~exist(loo_file, 'file') || ~exist(selection_file, 'file')
    error(['run_scale_free_test: %s or %s not found; they hold step 3''s cluster ' ...
           'voxel for voxel and its mass under every split. Run run_per_mouse_values ' ...
           'with the same comparison and smoothing.'], loo_file, selection_file);
end
S_loo = load(loo_file, 'loo_clusters', 'loo_box', 'loo_region');
if ~strcmp(S_loo.loo_region, region)
    error(['run_scale_free_test: the leave-one-out''s clusters in %s are of %s, not ' ...
           '%s.'], loo_file, S_loo.loo_region, region);
end
S_sel = load(selection_file, 'sm');
sm = S_sel.sm;

% fold 0's voxels, which must be the selection-matched test's observed cluster
% (both on the leave-one-out's box)
if ~isequal(sort(double(sm.cluster_voxels{1}(:))), ...
        sort(double(S_loo.loo_clusters{1}(:))))
    error(['run_scale_free_test: the selection-matched test''s cluster (%s) is not ' ...
           'the leave-one-out''s fold 0 (%s).'], selection_file, loo_file);
end

% fold 0's voxels, from the leave-one-out's box to the folded grid
box = S_loo.loo_box(1);
[ap, dv, ml] = ind2sub(box.box_size, double(S_loo.loo_clusters{1}));
lin = sub2ind(masks.folded_size, box.box_ap(ap)', box.box_dv(dv)', ...
    box.box_ml_left(ml)');

% its mass and p over the splits: of either sign, the larger of the two signs'
% masses, as step 3's region test takes it; and of the experimental group
% higher alone
step3 = struct();
step3.lin = uint32(sort(lin(:)));
step3.n = numel(step3.lin);
step3.mass = sm.cluster_mass(1, 1);
larger = max(sm.cluster_mass, [], 2);
step3.p = mean(larger >= larger(1));
step3.p_one = mean(sm.cluster_mass(:, 1) >= sm.cluster_mass(1, 1));

% step 3's table: its heaviest cluster of L - R in the region, which is fold 0's
% when it is the experimental group's
T_step3 = readtable(step3_table, 'Delimiter', ',');
row = strcmp(T_step3.acronym, region);
table_n = T_step3.lr_diff_cluster_n(row);
table_mass = T_step3.lr_diff_cluster_score(row);
table_p = T_step3.lr_diff_cluster_p_perm(row);
if table_mass <= 0 || step3.n ~= table_n || ...
        abs(step3.mass - table_mass) > 1e-6 * table_mass || abs(step3.p - table_p) > 1e-9
    error(['run_scale_free_test: the cluster of %s, %d voxels, mass %.2f, p %.4g ' ...
           '(%s, %s); step 3''s table gives %d voxels, a score of %.2f, p %.4g.'], ...
           region, step3.n, step3.mass, step3.p, loo_file, selection_file, table_n, ...
           table_mass, table_p);
end
fprintf(['  step 3''s cluster in %s: %d voxels, mass %.2f, p %.3g (either sign, as ' ...
         'its table), %.3g (one-sided)\n'], region, step3.n, step3.mass, step3.p, ...
         step3.p_one);
end

function result = test_result(perm, geom, masks, apriori_row, ctrl_mice, exp_mice)
% What the tables and figures read from the test: the mice and their off-tissue
% levels; the test without its values per voxel; the planes of every region's
% heaviest cluster of each sign; the voxels of those of the region named in
% advance; and the observed t in the box around it.

result = struct();
result.ctrl_names = ctrl_mice.names;
result.exp_names = exp_mice.names;
result.ctrl_background = ctrl_mice.background;
result.exp_background = exp_mice.background;
result.n_candidates = numel(geom.cand_lin);

% the planes of each heaviest cluster (regions x sign x first and last), NaN
% for none
detail = perm.maps{1}.detail;
n_regions = size(detail.cluster_voxels, 1);
result.cluster_planes = nan(n_regions, 2, 2);
for r = 1:n_regions
    for s = 1:2
        lin = double(detail.cluster_voxels{r, s});
        if ~isempty(lin)
            planes = mod(lin - 1, geom.grid_size(1)) + 1;
            result.cluster_planes(r, s, :) = [min(planes), max(planes)];
        end
    end
end
result.apriori_clusters = {uint32(detail.cluster_voxels{apriori_row, 1}), ...
    uint32(detail.cluster_voxels{apriori_row, 2})};

% the observed t in the box around the region named in advance
[ap, dv, ml] = ind2sub(geom.grid_size, double(geom.cand_lin));
in_box = ap >= masks.box_ap(1) & ap <= masks.box_ap(end) & ...
    dv >= masks.box_dv(1) & dv <= masks.box_dv(end) & ...
    ml >= masks.box_ml_left(1) & ml <= masks.box_ml_left(end);
box_lin = sub2ind(masks.box_size, ap(in_box) - masks.box_ap(1) + 1, ...
    dv(in_box) - masks.box_dv(1) + 1, ml(in_box) - masks.box_ml_left(1) + 1);
result.t_box = nan(masks.box_size, 'single');
result.t_box(box_lin) = detail.t(in_box);

% the test without the values per voxel
perm.maps{1}.detail = rmfield(detail, {'rolled_unsigned', 't', 'cluster_voxels'});
result.perm = perm;
end

function T_out = region_table(result, T_regions, ctrl_type, exp_type, ccf_first_index)
% One row per region of the bars: each sign's heaviest cluster, its mass,
% voxels, signed peak, its own p and its planes as Allen's CCF index; the score
% (the heavier, with its sign), its p and corrected p.

perm = result.perm;
map = perm.maps{1};
k = strcmp(perm.measure_names, 'cluster');
observed = perm.splits.observed;
n_splits = size(perm.splits.in_ctrl, 1);
null_pos = map.null_pos(:, :, k);
null_neg = map.null_neg(:, :, k);

T_out = T_regions(:, {'acronym', 'name', 'group', 'n_voxels', 'a_priori'});
T_out.n_with_t = map.detail.n_with_t;
signs = {[exp_type '_higher_'], null_pos, 1; [ctrl_type '_higher_'], null_neg, 2};
for s = 1:2
    prefix = signs{s, 1};
    null_sign = signs{s, 2};
    mass = null_sign(observed, :)';
    T_out.([prefix 'mass']) = mass;
    T_out.([prefix 'voxels']) = map.detail.cluster_n(:, signs{s, 3});
    peak = map.detail.cluster_peak(:, signs{s, 3});

    % a region without a cluster of the sign has a peak of 0, not the -0 the
    % negative sign gives it
    peak(peak == 0) = 0;
    T_out.([prefix 'peak']) = peak;
    p_sign = nan(size(mass));
    for r = 1:numel(mass)
        if ~isnan(mass(r))
            p_sign(r) = nnz(null_sign(:, r) >= mass(r)) / n_splits;
        end
    end
    T_out.([prefix 'p']) = p_sign;

    % the planes as Allen's CCF index, counted from 0: plane p of the crop is
    % index p + ccf_first_index - 1
    T_out.([prefix 'ccf_first']) = result.cluster_planes(:, s, 1) + ccf_first_index - 1;
    T_out.([prefix 'ccf_last']) = result.cluster_planes(:, s, 2) + ccf_first_index - 1;
end
T_out.score = map.score(:, k);
T_out.p = map.p_perm(:, k);
T_out.p_corrected = map.p_fwer(:, k);
end

function T_clusters = cluster_table(result, step3, masks, apriori_row, ctrl_type, ...
    exp_type, ccf_first_index)
% The region named in advance: its heaviest cluster of each sign on the index
% and step 3's cluster on the test's maps; for each, its voxels and mass, its
% one-sided p (the share of the splits whose heaviest cluster of its sign is as
% heavy; 1 for a sign without a cluster, whose mass of 0 every split reaches),
% the p of either sign for the heavier sign (the region's p, step 3's for its
% cluster), its planes, centre and distance from the midline, and its overlap
% with step 3's cluster.

perm = result.perm;
map = perm.maps{1};
k = strcmp(perm.measure_names, 'cluster');
observed = perm.splits.observed;
n_splits = size(perm.splits.in_ctrl, 1);
null_pos = map.null_pos(:, apriori_row, k);
null_neg = map.null_neg(:, apriori_row, k);

% the p of either sign belongs to the heavier sign's cluster, the score's
score_sign = sign(map.score(apriori_row, k));
p_either = [NaN; NaN];
if score_sign > 0
    p_either(1) = map.p_perm(apriori_row, k);
elseif score_sign < 0
    p_either(2) = map.p_perm(apriori_row, k);
end

rows = {
    [exp_type ' higher, index'], result.apriori_clusters{1}, null_pos(observed), ...
        nnz(null_pos >= null_pos(observed)) / n_splits, p_either(1)
    [ctrl_type ' higher, index'], result.apriori_clusters{2}, null_neg(observed), ...
        nnz(null_neg >= null_neg(observed)) / n_splits, p_either(2)
    [exp_type ' higher, step 3 (test maps)'], step3.lin, step3.mass, step3.p_one, ...
        step3.p
    };
n_rows = size(rows, 1);
T_clusters = table();
T_clusters.cluster = rows(:, 1);
T_clusters.voxels = zeros(n_rows, 1);
T_clusters.volume_mm3 = zeros(n_rows, 1);
T_clusters.mass = cell2mat(rows(:, 3));
T_clusters.p_one_sided = cell2mat(rows(:, 4));
T_clusters.p_either_sign = cell2mat(rows(:, 5));
T_clusters.ccf_first = nan(n_rows, 1);
T_clusters.ccf_last = nan(n_rows, 1);
T_clusters.centre_ccf = nan(n_rows, 1);
T_clusters.centre_dv_mm = nan(n_rows, 1);
T_clusters.midline_mm = nan(n_rows, 1);
T_clusters.overlap_step3 = zeros(n_rows, 1);
T_clusters.share_in_step3 = nan(n_rows, 1);
T_clusters.share_of_step3 = nan(n_rows, 1);
midline_ml = masks.folded_size(3) + 0.5;
for i = 1:n_rows
    lin = double(rows{i, 2});
    T_clusters.voxels(i) = numel(lin);
    T_clusters.volume_mm3(i) = numel(lin) * 1e-6;
    if isempty(lin)
        continue
    end

    % where it sits: its planes as Allen's CCF index, its centre, its distance
    % from the midline (between the folded grid's last column and its mirror)
    [ap, dv, ml] = ind2sub(masks.folded_size, lin);
    T_clusters.ccf_first(i) = min(ap) + ccf_first_index - 1;
    T_clusters.ccf_last(i) = max(ap) + ccf_first_index - 1;
    T_clusters.centre_ccf(i) = mean(ap) + ccf_first_index - 1;
    T_clusters.centre_dv_mm(i) = mean(dv) * 0.01;
    T_clusters.midline_mm(i) = (midline_ml - mean(ml)) * 0.01;

    % its voxels in step 3's cluster
    T_clusters.overlap_step3(i) = nnz(ismember(lin, double(step3.lin)));
    T_clusters.share_in_step3(i) = T_clusters.overlap_step3(i) / numel(lin);
    T_clusters.share_of_step3(i) = T_clusters.overlap_step3(i) / numel(step3.lin);
end
end

function named = named_test(result, apriori_row, region, direction_named, ctrl_type, ...
    exp_type)
% The test named in advance: the region's heaviest cluster where the index is
% higher in the experimental group, its mass against the same in every split
% (one-sided), the p of either sign, and the lines that say so under the bars'
% title, the one named first first.

perm = result.perm;
map = perm.maps{1};
k = strcmp(perm.measure_names, 'cluster');
observed = perm.splits.observed;
n_splits = size(perm.splits.in_ctrl, 1);
null_pos = map.null_pos(:, apriori_row, k);
null_neg = map.null_neg(:, apriori_row, k);

named = struct();
named.region = region;
named.n_splits = n_splits;
named.mass_exp = null_pos(observed);
named.n_reach_exp = nnz(null_pos >= named.mass_exp);
named.p_one = named.n_reach_exp / n_splits;
named.mass_ctrl = null_neg(observed);
named.p_ctrl = nnz(null_neg >= named.mass_ctrl) / n_splits;
named.score = map.score(apriori_row, k);
named.p_either = map.p_perm(apriori_row, k);
named.p_corrected = map.p_fwer(apriori_row, k);

% the lines under the bars' title
one_sided = sprintf(['%s, %s higher (one-sided): heaviest cluster mass %.0f, p = %.3g ' ...
    '(%d of %d splits as heavy)'], region, exp_type, named.mass_exp, named.p_one, ...
    named.n_reach_exp, n_splits);
either = sprintf(['%s, either sign (step 3''s score): %+.0f, p = %.3g, corrected p = ' ...
    '%.3g; %s higher mass %.0f, p = %.3g'], region, named.score, named.p_either, ...
    named.p_corrected, ctrl_type, named.mass_ctrl, named.p_ctrl);
if direction_named
    named.lines = {['Named in advance (Gambino et al. 2014): ' one_sided], either};
else
    named.lines = {['Named in advance: ' either], ...
        ['The direction carried over from RWS: ' one_sided]};
end
end

function line = comparison_line(ctrl_type, exp_type, result)
% The comparison, its mice and their off-tissue levels, for the summary.

line = sprintf('%s (%s) against %s (%s)', exp_type, strjoin(short_names( ...
    result.exp_names), ', '), ctrl_type, strjoin(short_names(result.ctrl_names), ', '));
end

function names = short_names(names)
% The mice's names without the line's suffix.

names = strrep(names, '_Gria1', '');
end

function print_summary(result, T_out, T_clusters, named, comparison)
% The run's main numbers.

fprintf('\nScale-free test: %s, %d splits\n', comparison, named.n_splits);
fprintf('  off-tissue levels: control %s; experimental %s\n', ...
    mat2str(round(result.ctrl_background')), mat2str(round(result.exp_background')));
fprintf('  %d candidate voxels; maps %.1f min, test %.1f min\n', result.n_candidates, ...
    result.maps_min, result.test_min);
for line = named.lines
    fprintf('  %s\n', line{1});
end
for i = 1:height(T_clusters)
    fprintf(['  %s: %d voxels, mass %.1f, p %.3g one-sided, %.3g either sign, CCF %g ' ...
             'to %g, %.2f mm from the midline; %d voxels in step 3''s cluster (%.2f of ' ...
             'it, %.2f of step 3''s)\n'], T_clusters.cluster{i}, T_clusters.voxels(i), ...
             T_clusters.mass(i), T_clusters.p_one_sided(i), T_clusters.p_either_sign(i), ...
             T_clusters.ccf_first(i), T_clusters.ccf_last(i), T_clusters.midline_mm(i), ...
             T_clusters.overlap_step3(i), T_clusters.share_in_step3(i), ...
             T_clusters.share_of_step3(i));
end
passing = T_out.p_corrected < 0.05;
if any(passing)
    fprintf('  regions at a corrected p < 0.05: %s\n', strjoin(compose('%s (%+.0f, %.3g)', ...
        string(T_out.acronym(passing)), T_out.score(passing), T_out.p_corrected(passing)), ', '));
else
    fprintf('  no region at a corrected p < 0.05 (smallest %.3g)\n', min(T_out.p_corrected));
end
[~, order] = sort(abs(T_out.score), 'descend', 'MissingPlacement', 'last');
top = order(1:min(5, numel(order)));
fprintf('  the five largest scores: %s\n', strjoin(compose('%s %+.0f (p %.3g)', ...
    string(T_out.acronym(top)), T_out.score(top), T_out.p(top)), ', '));
end

% ===== Local functions: the t map =====

function plot_t_map(result, step3, masks, comparison, file_tag, comp_out_dir)
% The observed Welch t of the index in the box around the region named in
% advance: a coronal view at each cluster's planes (the mean t over them), the
% index's heaviest cluster of each sign and step 3's, within window_um of the
% clusters there; and a view from above of the whole region (the mean over DV
% of the t in it); every cluster outlined where it is in the view, the region
% in grey.

% the coronal views reach this far beyond the clusters they show: about two
% barrel columns of 300 um (Lefort et al. 2009), as run_mouse_influence's maps
window_um = 600;

% the clusters on the box's grid: the index's two, then step 3's
names = {sprintf('the index''s %s-higher cluster', comparison.exp_type), ...
    sprintf('the index''s %s-higher cluster', comparison.ctrl_type), ...
    sprintf('step 3''s %s-higher cluster (test maps)', comparison.exp_type)};
voxel_lists = {result.apriori_clusters{1}, result.apriori_clusters{2}, step3.lin};

% the index's clusters drawn thicker than step 3's, so that they show where
% they lie on its edge
line_styles = {'--', ':', '-'};
line_widths = [2, 2, 1.2];
is_cluster = cell(1, 3);
for c = 1:3
    is_cluster{c} = box_mask(voxel_lists{c}, masks);
end
is_region = false(masks.box_size);
is_region(masks.band_lin(masks.band_in_region)) = true;

% the axes: ML from the midline in mm (medial to the right), DV from the top of
% the volume in mm, AP as Allen's CCF index and in mm from the CCF's front (the
% index's 10 um), so that every view keeps its proportions
x_mm = (masks.box_ml_left - (masks.folded_size(3) + 0.5)) * 0.01;
y_dv_mm = masks.box_dv * 0.01;
y_ccf = masks.box_ap + comparison.ccf_first_index - 1;
y_ap_mm = y_ccf * 0.01;
c_map = sep_palette('difference');
limits = [-comparison.t_limit comparison.t_limit];

fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'Normalized', ...
    'Position', [0 0 1 1]);

% a coronal view at each cluster's planes
for c = 1:3
    subplot(2, 2, c);
    [planes, ~, ~] = ind2sub(masks.box_size, find(is_cluster{c}));
    if isempty(planes)
        text(0.5, 0.5, {sprintf('%s:', names{c}), sprintf(['none in %s: nowhere does ' ...
            'the surprise''s'], comparison.region), sprintf(['median over +/- %d ' ...
            'planes reach p < %g with that sign'], comparison.slab_range, ...
            comparison.cluster_p), ['(the colours of the other views are each ' ...
            'voxel''s own t)']}, ...
            'HorizontalAlignment', 'center', 'Interpreter', 'none', 'FontSize', 10);
        axis off;
        continue
    end
    planes = min(planes):max(planes);

    % the window: every cluster's footprint over these planes and window_um
    % around it, within the box
    footprints = cell(1, 3);
    for d = 1:3
        footprints{d} = squeeze(any(is_cluster{d}(planes, :, :), 1));
    end
    [rows, cols] = find(footprints{1} | footprints{2} | footprints{3});
    margin = round(window_um / 10);
    dv_win = max(1, min(rows) - margin):min(masks.box_size(2), max(rows) + margin);
    ml_win = max(1, min(cols) - margin):min(masks.box_size(3), max(cols) + margin);

    % the mean t over the planes, the region and the clusters outlined
    view_t = squeeze(mean(result.t_box(planes, dv_win, ml_win), 1, 'omitnan'));
    draw_t_view(view_t, x_mm(ml_win), y_dv_mm(dv_win), c_map, limits);
    region_view = squeeze(any(is_region(planes, dv_win, ml_win), 1));
    draw_outline(x_mm(ml_win), y_dv_mm(dv_win), region_view, ...
        sep_palette('paired_lines'), 0.8, '-');
    for d = 1:3
        draw_outline(x_mm(ml_win), y_dv_mm(dv_win), footprints{d}(dv_win, ml_win), ...
            [0 0 0], line_widths(d), line_styles{d});
    end
    xlabel('ML from the midline (mm), medial >');
    ylabel('DV (mm from the top of the volume)');
    title({sprintf('%s: %d voxels', names{c}, numel(voxel_lists{c})), ...
        sprintf('the mean t over its planes, CCF %d to %d', y_ccf(planes(1)), ...
        y_ccf(planes(end)))}, 'Interpreter', 'none', 'FontSize', 10);
end

% from above: the mean over DV of the t in the region, in each AP and ML column
subplot(2, 2, 4);
t_region = result.t_box;
t_region(~is_region) = NaN;
view_t = squeeze(mean(t_region, 2, 'omitnan'));
draw_t_view(view_t, x_mm, y_ap_mm, c_map, limits);
draw_outline(x_mm, y_ap_mm, squeeze(any(is_region, 2)), sep_palette('paired_lines'), ...
    0.8, '-');
for d = 1:3
    draw_outline(x_mm, y_ap_mm, squeeze(any(is_cluster{d}, 2)), [0 0 0], ...
        line_widths(d), line_styles{d});
end
xlabel('ML from the midline (mm), medial >');
ylabel('AP (mm from the CCF''s front: CCF index / 100)');
title(sprintf('from above: the mean over DV of the t in %s', comparison.region), ...
    'FontSize', 10);
cb = colorbar;
cb.Label.String = sprintf('Welch t of the index, %s - %s', comparison.exp_type, ...
    comparison.ctrl_type);

% the title: what is shown and what the outlines are
title_line = sprintf(['Welch t of the asymmetry index |L - R| / (L + R) on the raw ' ...
    'stacks, %s - %s, in %s - %s'], comparison.exp_type, comparison.ctrl_type, ...
    comparison.region, strrep(file_tag, '_', ' '));
outline_line = sprintf(['Outlines (each where it is in the view): dashed, the ' ...
    'index''s heaviest %s-higher cluster; dotted, its heaviest %s-higher cluster; ' ...
    'thin solid black, step 3''s cluster; grey, %s. Colour: %g to %g; grey: no t.'], ...
    comparison.exp_type, comparison.ctrl_type, comparison.region, limits);
view_line = sprintf(['Coronal views: the mean t over the cluster''s planes, within ' ...
    '%d um of the clusters; from above: the mean over DV of the t in the region.'], ...
    window_um);
sgtitle({title_line, ['\rm\fontsize{10}' outline_line], ['\rm\fontsize{10}' ...
    view_line]}, 'FontSize', 13, 'FontWeight', 'bold', 'Interpreter', 'tex');
saveas(fig, fullfile(comp_out_dir, ['Scale_Free_TMap_' file_tag '.fig']));
exportgraphics(fig, fullfile(comp_out_dir, ['Scale_Free_TMap_' file_tag '.png']), ...
    'Resolution', 300);
end

function is_in = box_mask(lin, masks)
% Voxels of the folded grid (linear indices) as a mask of the box; those outside
% it are left out.

is_in = false(masks.box_size);
if isempty(lin)
    return
end
[ap, dv, ml] = ind2sub(masks.folded_size, double(lin));
ap = ap - masks.box_ap(1) + 1;
dv = dv - masks.box_dv(1) + 1;
ml = ml - masks.box_ml_left(1) + 1;
inside = ap >= 1 & ap <= masks.box_size(1) & dv >= 1 & dv <= masks.box_size(2) & ...
    ml >= 1 & ml <= masks.box_size(3);
is_in(sub2ind(masks.box_size, ap(inside), dv(inside), ml(inside))) = true;
end

function draw_t_view(img, x, y, c_map, limits)
% One view of the t on the no-data grey, its axes in mm or planes.

imagesc(x, y, img, 'AlphaData', double(~isnan(img)));
set(gca, 'Color', sep_palette('no_data'), 'FontSize', 8, 'TickDir', 'out');
axis image;
clim(limits);
colormap(gca, c_map);
hold on;
end

function draw_outline(x, y, mask, colour, width, style)
% The edge of a mask, as a contour at one half; nothing for an empty or full mask.

if ~any(mask(:)) || all(mask(:))
    return
end
contour(x, y, double(mask), [0.5 0.5], 'LineColor', colour, 'LineWidth', width, ...
    'LineStyle', style);
end
