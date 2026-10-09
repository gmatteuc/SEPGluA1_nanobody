function per_mouse_region_values(run_settings)
%PER_MOUSE_REGION_VALUES  Each mouse's values in the regions named in advance.
%   PER_MOUSE_REGION_VALUES(run_settings) does the work of run_per_mouse_values,
%   which sets the fields of run_settings and says what each one does.
%
%   Reads, from the control and the experimental group's folders,
%   <channel>_4d_normalized.mat and <channel>_4d_normalized_bkgmask.mat
%   (run_normalise_groups), the collected stack <channel>_4d.mat
%   (run_collect_by_group) and, with auto_control, the autofluorescence stack
%   auto_4d.mat, one mouse at a time. Writes into comp_out_dir, the folder of
%   run_group_differences, with <tag> = <comp_tag>_smooth<sigma>:
%     Per_Mouse_Values_<tag>.csv   one row per mouse, every value below
%     Per_Mouse_Stats_<tag>.csv    one row per value: group means and SEMs, the
%                                  p of its test, Welch t, Hedges g
%     Per_Mouse_Values_<tag>       the figure, .fig and .png: the selection-
%                                  matched values first, then the regions
%     Per_Mouse_LOO_<tag>          the leave-one-out's figure, .fig and .png;
%                                  .csv, one row per fold; .mat, the folds'
%                                  clusters, voxel by voxel
%     Per_Mouse_Maps_<tag>.mat     the cache of each mouse's maps (see below)
%     Per_Mouse_Auto_<tag>.mat     the cache of each mouse's autofluorescence
%     Per_Mouse_Selection_<tag>.mat  the selection-matched test's splits
%     Per_Mouse_LOO_Relabelled_<tag>.mat  the leave-one-out redone under every
%                                  split, with loo_relabel (see below)
%
%   Each mouse is taken as run_group_differences takes it: NaN outside its
%   tissue and smoothed (tissue_only), folded onto the left hemisphere
%   (compute_lr_stats), the experimental mice put on the control group's scale
%   by the line between the groups' plane profiles (plane_tissue_means,
%   align_exp_to_ctrl). The regions are those of the bars (surprise_regions).
%   In each region named in advance, over the mouse's own voxels with a value:
%     ai           mean |L - R| over mean (L + R), on the test's maps
%     ai_raw       the same on the collected stack, the mouse's off-tissue
%                  level subtracted first and the same smoothing, so it
%                  depends on neither the normalisation nor the alignment
%     sum_rel      mean (L + R) in the region over its mean over the isocortex
%                  (every isocortical area of the bars), test maps
%     sum_rel_raw  the same on the raw stack
%   No signed value: which hemisphere is left is not certain for every brain,
%   so the stimulated side is not known mouse by mouse, and L - R would carry
%   an arbitrary sign; the test takes |L - R| and L + R for the same reason.
%   The off-tissue level is the median of the raw stack over the mouse's
%   background voxels (run_normalise_groups' mask) outside the atlas brain that
%   a section reached: the slide around the section, about 500. The alignment
%   is a line, so a region's means on the aligned maps follow from its means
%   before it: mean |L - R| times the slope, mean (L + R) times the slope plus
%   twice the intercept, both over the common factor, which cancels in every
%   ratio. Its intercept enters the experimental mice's L + R, and so their ai
%   and sum_rel: the raw values have no such term. The zero of the test's maps
%   is the normalisation's, not the absence of signal: a mouse's off-tissue
%   level lands between about -2,400 and +500 on the control group's scale
%   (before the common factor), so the raw ai is the one to read as a ratio;
%   an ai over a mean L + R not above 0 is NaN.
%
%   Selection-matched test (selection_matched; the main per-mouse test): with
%   the true groups, the heaviest cluster of cluster_region where |L - R| is
%   higher in the experimental group (a positive t), found with all the mice by
%   region_permutation_test as run_group_differences finds it (Welch t, the
%   fewest mice per group, surprise, the rolling median, p < cluster_p,
%   connectivity), and every mouse's ai (sm_ai), ai_raw (sm_ai_raw) and, with
%   auto_control, the ai of its autofluorescence (sm_ai_auto) in it. The
%   statistic is the experimental group's mean minus the control group's. Under
%   every split of the mice into groups of the observed sizes the cluster is
%   found again with the split's groups (positive: the split's experimental
%   group higher; negative: its control group higher), every mouse read in it,
%   the same statistic taken; the mice keep the alignment of their true groups,
%   as the region test keeps it. One-sided p: the share of the splits whose
%   statistic in the positive cluster reaches the observed one, the observed
%   split included. Two-sided: per split the larger of that statistic and of
%   control minus experimental in the negative cluster, as the region test's
%   score takes the larger of its two signs. A split without a cluster, or with
%   a group without a mouse with a value, gives 0; the p with those splits left
%   out is given too. The selection is redone in every split, so the null
%   carries it as the observed statistic does; shuffling the observed values
%   would not, and its p (in the table only) is anti-conservative. The first
%   split must give run_group_differences' cluster and the cluster mass over
%   the splits its p, both checked against its table when it is there. It has
%   that test's selection and splits; what it takes from each split's cluster
%   differs: the mice's mean AI, which does not grow with the cluster's
%   extent, where that test takes the cluster's mass, which does.
%
%   Autofluorescence (auto_control): the auto stack, collected from the same
%   registered volumes as the nano stack (an auto_4d.mat without the files it
%   was read from is the older one of another registration, and is refused),
%   checked against the nano stack (the voxels a section reached must agree),
%   less its own off-tissue level over the same voxels, NaN outside the nano
%   channel's tissue and smoothed alike, as the raw stack.
%
%   Leave-one-out: each mouse in turn is left out, and the comparison is done
%   again on the others as the test does it (the alignment refitted on their
%   profiles), here by region_permutation_test on the observed split alone,
%   within cluster_region. The left-out mouse's ai (on that fold's maps,
%   loo_ai) and ai_raw (loo_ai_raw) are read in the fold's heaviest cluster
%   where |L - R| is higher in the experimental group, so no mouse is read in a
%   cluster its own data helped define: out of sample, each mouse read in a
%   cluster of four mice against five, stricter about circularity than the
%   selection-matched test by design (its p need not be the higher). Fold 0 leaves
%   no mouse out: its cluster must be that of run_group_differences, which is
%   checked against its table when it is there. A fold with no cluster gives
%   NaN. Not redone without the mouse: run_normalise_groups, which fits every
%   mouse of a group onto the group's median cortex, so without one mouse the
%   others' lines would change too. To first order they change by one line
%   common to the group, which leaves the fold's t of |L - R| as it is (the
%   fold's alignment is refitted, and a scale common to both groups does not
%   change a t); what is left is second order. With loo_relabel the
%   leave-one-out is redone under every split, each fold's cluster found again
%   with the split's groups (the fold keeping its alignment), and its p is the
%   share of the splits with a value in both groups whose difference of the
%   group means reaches the observed one.
%
%   Statistics, per value, over the mice with a value: the difference of the
%   group means (experimental minus control) and its p, one-sided for the
%   experimental group higher and two-sided; the Welch t for reference; Hedges'
%   g. For a region's values the p is the exact permutation of the values over
%   every split of the mice (252 for 5 and 5), the observed one included; for
%   the cluster values, that of their test above, the permutation of the values,
%   the Welch p and Hedges' g being kept in the table under names that say they
%   ignore the selection (anti-conservative, inflated). With direction_named, the
%   one-sided p is the test named before any number and comes first (RWS
%   potentiates the stimulated barrels' synapses and brings AMPA receptors to
%   their surface, Gambino et al. 2014); without it, the two-sided p comes
%   first.
%
%   The maps are cached: each mouse's region sums and the box around
%   cluster_region that the clusters are read in, so the statistics and the
%   figures can be redone without smoothing every mouse again;
%   force_recompute_mice redoes them. The cache is read only if each group's
%   normalised stack still holds the mice and the lines of run_normalise_groups
%   it was made from. The autofluorescence, the selection-matched test and the
%   redone leave-one-out are cached beside it, each read only when made from
%   these maps and settings.

% settings of run_per_mouse_values, under the names the code below uses
paths = run_settings.paths;
ctrl_type = run_settings.ctrl_type;
exp_type = run_settings.exp_type;
behavior_mice = run_settings.behavior_mice;
regions = run_settings.regions;
cluster_region = run_settings.cluster_region;
apply_smoothing = run_settings.apply_smoothing;
smooth_sigma = run_settings.smooth_sigma;
min_mice_per_group = run_settings.min_mice_per_group;
slab_range = run_settings.slab_range;
cluster_p = run_settings.cluster_p;
cluster_connectivity = run_settings.cluster_connectivity;
force_recompute_mice = run_settings.force_recompute_mice;
direction_named = run_settings.direction_named;
selection_matched = run_settings.selection_matched;
auto_control = run_settings.auto_control;
loo_relabel = run_settings.loo_relabel;
channel = run_settings.channel;
comp_tag = run_settings.comp_tag;
ctrl_dir = run_settings.ctrl_dir;
exp_dir = run_settings.exp_dir;
comp_out_dir = run_settings.comp_out_dir;

% the experimental groups of run_group_differences
if ~ismember(exp_type, {'rws', 'behavior'})
    error('run_per_mouse_values: unknown exp_type ''%s'' (use ''rws'' or ''behavior'').', ...
        exp_type);
end

% the clusters' region must be one of the regions named
if ~ismember(cluster_region, regions)
    error(['run_per_mouse_values: cluster_region is %s, not one of the regions named ' ...
           'in regions (%s).'], cluster_region, strjoin(regions, ', '));
end

% the autofluorescence is read in the selection-matched test's clusters of the
% nano channel
if auto_control && (~selection_matched || ~strcmp(channel, 'nano'))
    error(['run_per_mouse_values: auto_control reads the autofluorescence in the ' ...
           'clusters of the selection-matched test of the nano channel; set ' ...
           'selection_matched = true and channel = ''nano'', or auto_control = false.']);
end

% the file names carry the smoothing
if apply_smoothing
    smooth_suffix = sprintf('_smooth%g', smooth_sigma);
else
    smooth_suffix = '_nosmooth';
end
file_tag = [comp_tag smooth_suffix];

% the clusters' test: the observed split alone, on no pool; the top volume and
% the quantile at run_group_differences' values, though only the cluster is read
perm_settings = struct();
perm_settings.min_mice_per_group = min_mice_per_group;
perm_settings.slab_range = slab_range;
perm_settings.p_thresh = cluster_p;
perm_settings.cluster_p = cluster_p;
perm_settings.cluster_connectivity = cluster_connectivity;
perm_settings.topvol_mm3 = 0.1;
perm_settings.region_quantile = 0.99;
perm_settings.n_permutations = 1;
perm_settings.n_workers = 0;
perm_settings.seed = 0;
perm_settings.keep_cluster_voxels = true;

% the table of run_group_differences, which the clusters are checked against
step3_table = fullfile(comp_out_dir, ['Region_Surprise_DiffSum_' comp_tag '.csv']);

%% Atlas and regions

% the 10 um annotation on the volumes' crop, the regions of the bars, and the
% voxels of each region named in advance, of the isocortex, and of the box the
% clusters are read in
A = get_atlas_crop('ccf');
brainMask = A.brainMask;
[T_regions, valid_pixels, region_of_voxel] = surprise_regions(A.annot, paths.atlas);
clear A
masks = region_masks(T_regions, valid_pixels, region_of_voxel, regions, cluster_region, ...
    slab_range);
clear valid_pixels region_of_voxel

%% Each mouse of both groups

% every control mouse saved; for behavior, the mice behavior_mice names
if strcmp(exp_type, 'behavior')
    exp_named = behavior_mice;
else
    exp_named = {};
end

% the settings the maps depend on, checked when the cache is read; the clusters'
% region keeps the field name of the caches made before 9 October 2026, when
% only the leave-one-out read its clusters
cache_settings = struct('regions', {regions}, 'loo_region', cluster_region, ...
    'slab_range', slab_range, 'apply_smoothing', apply_smoothing, ...
    'smooth_sigma', smooth_sigma, 'exp_named', {exp_named});

% from the cache, if the groups' normalised stacks are still those it was made
% from, or one mouse at a time from the stacks
cache_file = fullfile(comp_out_dir, ['Per_Mouse_Maps_' file_tag '.mat']);
if exist(cache_file, 'file') && ~force_recompute_mice
    fprintf('Loading the mice''s maps from %s...\n', cache_file);
    S_cache = load(cache_file, 'ctrl_mice', 'exp_mice', 'cache_settings');
    if ~isequal(S_cache.cache_settings, cache_settings)
        error(['run_per_mouse_values: the maps in %s were made with other regions, ' ...
               'box, smoothing or mice. Set force_recompute_mice = true.'], cache_file);
    end
    ctrl_mice = S_cache.ctrl_mice;
    exp_mice = S_cache.exp_mice;
    clear S_cache
    check_cache_against_stack(ctrl_mice, ctrl_dir, channel, cache_file);
    check_cache_against_stack(exp_mice, exp_dir, channel, cache_file);
    maps_recomputed = false;
else
    ctrl_mice = group_mice(ctrl_type, ctrl_dir, channel, {}, brainMask, masks, ...
        apply_smoothing, smooth_sigma);
    exp_mice = group_mice(exp_type, exp_dir, channel, exp_named, brainMask, masks, ...
        apply_smoothing, smooth_sigma);
    fprintf('Saving the mice''s maps to %s...\n', cache_file);
    save(cache_file, 'ctrl_mice', 'exp_mice', 'cache_settings', '-v7.3');
    maps_recomputed = true;
end

% what the clusters' tests depend on besides the maps, checked when their caches
% are read
split_settings = struct('cache_settings', cache_settings, ...
    'ctrl_norm_params', ctrl_mice.stack_norm_params, ...
    'exp_norm_params', exp_mice.stack_norm_params, ...
    'min_mice_per_group', min_mice_per_group, 'slab_range', slab_range, ...
    'cluster_p', cluster_p, 'cluster_connectivity', cluster_connectivity);

%% Autofluorescence

% each mouse's autofluorescence in the box the clusters are read in; from its
% cache when it was made from these maps and auto stacks, or one mouse at a
% time from the stacks
ctrl_auto = [];
exp_auto = [];
auto_settings = [];
if auto_control
    auto_file = fullfile(comp_out_dir, ['Per_Mouse_Auto_' file_tag '.mat']);
    auto_settings = struct('split_settings', split_settings, ...
        'ctrl_mice', {ctrl_mice.names}, 'exp_mice', {exp_mice.names}, ...
        'ctrl_auto_files', {auto_stack_record(ctrl_dir)}, ...
        'exp_auto_files', {auto_stack_record(exp_dir)});
    if exist(auto_file, 'file') && ~maps_recomputed
        S_auto = load(auto_file, 'ctrl_auto', 'exp_auto', 'auto_settings');
        if isequal(S_auto.auto_settings, auto_settings)
            fprintf('Loading the mice''s autofluorescence from %s...\n', auto_file);
            ctrl_auto = S_auto.ctrl_auto;
            exp_auto = S_auto.exp_auto;
        end
        clear S_auto
    end
    if isempty(ctrl_auto)
        ctrl_auto = group_auto(ctrl_mice, ctrl_dir, brainMask, masks, apply_smoothing, ...
            smooth_sigma);
        exp_auto = group_auto(exp_mice, exp_dir, brainMask, masks, apply_smoothing, ...
            smooth_sigma);
        fprintf('Saving the mice''s autofluorescence to %s...\n', auto_file);
        save(auto_file, 'ctrl_auto', 'exp_auto', 'auto_settings', '-v7.3');
    end
end
clear brainMask

%% Values on the test's maps and on the raw stacks

% the alignment of run_group_differences, on all the mice of both groups
[~, ~, ~, slope, intercept, common_factor] = align_exp_to_ctrl(ctrl_mice.profiles, ...
    exp_mice.profiles);
if slope <= 0
    error(['run_per_mouse_values: the alignment of %s onto %s has slope %.3f; a ' ...
           'line that turns the profile over cannot align it.'], exp_type, ctrl_type, ...
           slope);
end

% each mouse's values, one row per mouse, control mice first
T_mice = mouse_values(ctrl_mice, exp_mice, slope, intercept, common_factor, regions);
if auto_control
    T_mice.background_auto = [ctrl_auto.background; exp_auto.background];
    T_mice.auto_reached_dice = [ctrl_auto.reached_dice; exp_auto.reached_dice];
end

%% Selection-matched test

% the cluster of all the mice and every mouse read in it, and the same under
% every split of the mice; from its cache when it was made from these maps,
% autofluorescence and settings
sm = [];
sm_place = [];
if selection_matched
    sm_file = fullfile(comp_out_dir, ['Per_Mouse_Selection_' file_tag '.mat']);
    sm_settings = struct('split_settings', split_settings, 'auto_settings', ...
        auto_settings);
    if exist(sm_file, 'file') && ~maps_recomputed
        S_sm = load(sm_file, 'sm', 'sm_settings');
        if isequal(S_sm.sm_settings, sm_settings)
            fprintf('Loading the selection-matched test''s splits from %s...\n', sm_file);
            sm = S_sm.sm;
        end
        clear S_sm
    end
    if isempty(sm)
        sm = selection_splits(ctrl_mice, exp_mice, ctrl_auto, exp_auto, masks, ...
            perm_settings);
        save(sm_file, 'sm', 'sm_settings');
    end

    % its first split must give the cluster of run_group_differences, and the
    % cluster mass over its splits the p of its region test
    fprintf('Selection-matched test in %s:\n', cluster_region);
    check_against_step3(struct('n', sm.cluster_n(1, 1), 'mass', sm.cluster_mass(1, 1)), ...
        step3_table, cluster_region);
    [sm.mass_p, sm.mass_p_one] = cluster_mass_p(sm, step3_table, cluster_region);
    check_mirror_splits(sm);

    % every mouse's values in the cluster of all the mice, and where it sits; the
    % size of the cluster where the control group is higher, for the notes
    for r = 1:numel(sm.readings)
        T_mice.(['sm_' sm.readings{r}]) = squeeze(sm.values(1, :, 1, r))';
    end
    sm_place = cluster_place(sm.cluster_voxels{1}, masks);
    sm_place.n = sm.cluster_n(1, 1);
    sm_place.n_negative = sm.cluster_n(1, 2);
    sm_place.mass = sm.cluster_mass(1, 1);
    sm_place.mass_p = sm.mass_p;
    sm_place.mass_p_one = sm.mass_p_one;
end

%% Leave-one-out clusters

% each mouse's values in the cluster its group comparison finds without it
fprintf('Leave-one-out clusters in %s...\n', cluster_region);
[T_loo, loo] = leave_one_out(ctrl_mice, exp_mice, masks, perm_settings, step3_table, ...
    cluster_region);
T_mice.loo_ai = loo.ai;
T_mice.loo_ai_raw = loo.ai_raw;
T_mice.loo_coverage = loo.coverage;
T_mice.loo_cluster_n = loo.cluster_n;

% the folds' clusters, voxel by voxel, with the box they index into, to draw
% them or to see where they sit
loo_clusters = loo.cluster_voxels;
loo_box = struct('box_ap', masks.box_ap, 'box_dv', masks.box_dv, ...
    'box_ml_left', masks.box_ml_left, 'box_size', masks.box_size, ...
    'left_out', {T_loo.left_out});
loo_region = cluster_region;
save(fullfile(comp_out_dir, ['Per_Mouse_LOO_' file_tag '.mat']), 'loo_clusters', ...
    'loo_box', 'loo_region');

%% Leave-one-out under every split

% the leave-one-out redone under every split of the mice, for the full p of its
% two values; from its cache when it was made from these stacks and settings,
% and the maps were not made again in this run
relabelled = [];
if loo_relabel
    relabel_file = fullfile(comp_out_dir, ['Per_Mouse_LOO_Relabelled_' file_tag '.mat']);
    relabel_settings = split_settings;
    if exist(relabel_file, 'file') && ~maps_recomputed
        S_relabel = load(relabel_file, 'relabelled', 'relabel_settings');
        if isequal(S_relabel.relabel_settings, relabel_settings)
            fprintf('Loading the leave-one-out under every split from %s...\n', ...
                relabel_file);
            relabelled = S_relabel.relabelled;
        end
        clear S_relabel
    end
    if isempty(relabelled)
        relabelled = loo_relabelled(ctrl_mice, exp_mice, masks, perm_settings);
        save(relabel_file, 'relabelled', 'relabel_settings');
    end

    % its first split, the groups as they are, must be the leave-one-out above
    if ~isequaln(relabelled.ai(1, :)', loo.ai) || ...
            ~isequaln(relabelled.ai_raw(1, :)', loo.ai_raw)
        error(['run_per_mouse_values: the first split of %s, the groups as they ' ...
               'are, does not give the leave-one-out''s values. Set ' ...
               'force_recompute_mice = true.'], relabel_file);
    end
end

%% Statistics and figures

% the values tested and drawn, in the figures' order
values = value_list(regions, cluster_region, sm);

% each value's group test and the tables
T_stats = value_statistics(T_mice, values, ctrl_type, exp_type, sm, relabelled);
writetable(T_mice, fullfile(comp_out_dir, ['Per_Mouse_Values_' file_tag '.csv']));
writetable(T_stats, fullfile(comp_out_dir, ['Per_Mouse_Stats_' file_tag '.csv']));
writetable(T_loo, fullfile(comp_out_dir, ['Per_Mouse_LOO_' file_tag '.csv']));
print_summary(T_mice, T_stats, T_loo, values, sm_place);

% one panel per value, one dot per mouse; the selection-matched values and the
% regions in one figure, the leave-one-out in another; their text from the
% settings
comparison = struct('ctrl_type', ctrl_type, 'exp_type', exp_type, ...
    'direction_named', direction_named, 'slope', slope, 'intercept', intercept, ...
    'perm_settings', perm_settings, 'cluster_region', cluster_region, ...
    'sm_place', sm_place);
plot_per_mouse_values(T_mice, T_stats, sm, values, regions, comparison, file_tag, ...
    comp_out_dir);
plot_loo_values(T_mice, T_stats, T_loo, values, comparison, relabelled, file_tag, ...
    comp_out_dir);
fprintf('Per-mouse values saved to: %s\n', comp_out_dir);

end

% ===== Local functions: regions and mice =====

function masks = region_masks(T_regions, valid_pixels, region_of_voxel, regions, ...
    cluster_region, slab_range)
% The voxels of each region named in advance and of the isocortex, as linear
% indices into the folded grid, and the box around cluster_region the clusters
% are read in: its planes and the slab_range planes on each side, which the
% rolling median of a voxel in it reads.

% the region of every voxel of the folded grid, 0 for none
region_vol = zeros(size(valid_pixels), 'uint16');
region_vol(valid_pixels) = region_of_voxel;

% each region named in advance, by its row of the bars' table
masks = struct();
masks.folded_size = size(valid_pixels);
masks.region_lin = cell(numel(regions), 1);
for r = 1:numel(regions)
    row = find(strcmp(T_regions.acronym, regions{r}));
    if isempty(row)
        error(['run_per_mouse_values: the region %s is not one of the bars'' ' ...
               'regions. Use their atlas acronyms: %s.'], regions{r}, ...
               strjoin(T_regions.acronym, ', '));
    end
    masks.region_lin{r} = uint32(find(region_vol == row));
    fprintf('  %s: %d voxels in the left hemisphere.\n', regions{r}, ...
        numel(masks.region_lin{r}));
end

% the isocortex: every isocortical area of the bars
iso_rows = find(strcmp(T_regions.group, 'Isocortex'));
masks.iso_lin = uint32(find(ismember(region_vol, iso_rows)));
fprintf('  isocortex: %d areas, %d voxels in the left hemisphere.\n', ...
    numel(iso_rows), numel(masks.iso_lin));

% the band of the clusters: the region and slab_range planes on each side along
% AP, within the atlas
in_region = region_vol == find(strcmp(T_regions.acronym, cluster_region));
band = imdilate(in_region, true(2 * slab_range + 1, 1)) & valid_pixels;
clear region_vol

% the band's bounding box: its planes and rows, its columns of the left
% hemisphere and, after them, their mirror images in the right hemisphere, so
% compute_lr_stats pairs each column of the box with its mirror image as it
% pairs them across the whole width
[ap, dv, ml] = ind2sub(size(band), find(band));
masks.box_ap = min(ap):max(ap);
masks.box_dv = min(dv):max(dv);
masks.box_ml_left = min(ml):max(ml);
masks.box_ml = [masks.box_ml_left, flip(2 * size(band, 3) + 1 - masks.box_ml_left)];
clear ap dv ml

% the band and the region in the box's folded grid
band_box = band(masks.box_ap, masks.box_dv, masks.box_ml_left);
region_box = in_region(masks.box_ap, masks.box_dv, masks.box_ml_left);
masks.box_size = size(band_box);
masks.band_lin = uint32(find(band_box));
masks.band_in_region = region_box(masks.band_lin);
fprintf(['  box of the clusters around %s: %d x %d x %d voxels, %d in the band, %d in ' ...
         'the region.\n'], cluster_region, masks.box_size, numel(masks.band_lin), ...
        nnz(masks.band_in_region));
end

function G = group_mice(group, group_dir, channel, named_mice, brainMask, masks, ...
    apply_smoothing, smooth_sigma)
% One group's mice, one at a time from the stacks: each mouse's plane profile,
% its sums in the regions on the normalised and on the raw stack, its
% off-tissue level, and both smoothed stacks in the clusters' box. The mice
% are those named, or every mouse the normalised stack holds.

% the normalised stack, its background masks and the collected stack, read one
% mouse at a time
norm_var_name = [channel '_4d_normalized'];
raw_var_name = [channel '_4d'];
norm_file = fullfile(group_dir, [channel '_4d_normalized.mat']);
mask_file = fullfile(group_dir, [channel '_4d_normalized_bkgmask.mat']);
raw_file = fullfile(group_dir, [channel '_4d.mat']);
M_norm = matfile(norm_file);
M_mask = matfile(mask_file);
M_raw = matfile(raw_file);
S_norm = load(norm_file, 'current_mice', 'norm_params');
S_mask = load(mask_file, 'current_mice');
if ~isequal(S_mask.current_mice, S_norm.current_mice)
    error(['run_per_mouse_values: the background masks of %s are of %s, the ' ...
           'normalised stack of %s. Rerun run_normalise_groups.'], group, ...
           strjoin(S_mask.current_mice, ', '), strjoin(S_norm.current_mice, ', '));
end

% the mice, by their place in the normalised stack
if isempty(named_mice)
    names = S_norm.current_mice;
    norm_idx = 1:numel(names);
else
    [is_saved, norm_idx] = ismember(named_mice, S_norm.current_mice);
    if ~all(is_saved)
        error('run_per_mouse_values: %s not among the mice of %s (%s).', ...
            strjoin(named_mice(~is_saved), ', '), norm_file, ...
            strjoin(S_norm.current_mice, ', '));
    end
    names = named_mice;
end

% their place in the collected stack, which run_collect_by_group fills in the
% cohort table's order
raw_idx = collected_index(names, group, M_raw, raw_var_name, raw_file);

% one mouse at a time
n_mice = numel(names);
G = struct();
G.group = group;
G.names = names;

% the normalised stack's mice and their lines, which the cache is checked against
G.stack_mice = S_norm.current_mice;
G.stack_norm_params = S_norm.norm_params;
G.profiles = nan(size(M_raw, raw_var_name, 1), n_mice);
G.test = cell(n_mice, 1);
G.raw = cell(n_mice, 1);
G.background = nan(n_mice, 1);
G.background_planes = nan(size(M_raw, raw_var_name, 1), n_mice);
G.box = cell(n_mice, 1);
G.box_raw = cell(n_mice, 1);
for k = 1:n_mice
    t_mouse = tic;
    fprintf('%s, mouse %d of %d: %s\n', group, k, n_mice, names{k});

    % the mouse's normalised volume, background mask and collected volume
    vol = M_norm.(norm_var_name)(:, :, :, norm_idx(k));
    bg_mask = M_mask.recomputed_bkg_mask_4d(:, :, :, norm_idx(k));
    raw = M_raw.(raw_var_name)(:, :, :, raw_idx(k));

    % the normalised volume must be this collected one through the mouse's line
    check_raw_against_normalised(raw, vol, S_norm.norm_params(norm_idx(k), :), names{k});

    % its plane profile, before the smoothing, as run_group_differences takes it
    G.profiles(:, k) = plane_tissue_means(vol, bg_mask);

    % the normalised volume: its tissue smoothed, folded, summed in the regions,
    % and kept in the box
    vol = tissue_only(vol, bg_mask, brainMask, apply_smoothing, smooth_sigma);
    [lr_diff, lr_sum] = compute_lr_stats(vol);
    G.test{k} = region_sums(lr_diff, lr_sum, masks);
    G.box{k} = vol(masks.box_ap, masks.box_dv, masks.box_ml);
    clear vol lr_diff lr_sum

    % the raw volume less its off-tissue level, through the same steps; a raw 0 is
    % a voxel no section reached, NaN as run_normalise_groups makes it
    [G.background(k), G.background_planes(:, k)] = off_tissue_level(raw, bg_mask, ...
        brainMask);
    fprintf(['  off-tissue level %.1f (planes: 10th to 90th percentile %.1f to ' ...
             '%.1f)\n'], G.background(k), prctile(G.background_planes(:, k), 10), ...
            prctile(G.background_planes(:, k), 90));
    raw(raw == 0) = NaN;
    raw = raw - G.background(k);
    raw = tissue_only(raw, bg_mask, brainMask, apply_smoothing, smooth_sigma);
    [lr_diff, lr_sum] = compute_lr_stats(raw);
    G.raw{k} = region_sums(lr_diff, lr_sum, masks);
    G.box_raw{k} = raw(masks.box_ap, masks.box_dv, masks.box_ml);
    clear raw lr_diff lr_sum bg_mask

    fprintf('  done in %.1f min.\n', toc(t_mouse) / 60);
end
end

function raw_idx = collected_index(names, group, M_raw, raw_var_name, raw_file)
% The mice's places in a group's collected nano stack, which run_collect_by_group
% fills in the cohort table's order.

cohort = get_cohort('groups', {group});
raw_names = {cohort.name};
raw_size = size(M_raw, raw_var_name);
if numel(raw_size) < 4
    raw_size(4) = 1;
end
if raw_size(4) ~= numel(raw_names)
    error(['run_per_mouse_values: %s holds %d mice, the cohort table %d for %s. ' ...
           'Rerun run_collect_by_group.'], raw_file, raw_size(4), numel(raw_names), ...
           group);
end
[~, raw_idx] = ismember(names, raw_names);
end

function source_files = auto_stack_record(group_dir)
% The registered volumes a group's auto stack was collected from, as
% run_collect_by_group records them; an auto_4d.mat without that record is the
% older one of another registration (19 November 2025), and is refused.

auto_file = fullfile(group_dir, 'auto_4d.mat');
if ~exist(auto_file, 'file')
    error(['run_per_mouse_values: no autofluorescence stack in %s. Collect it with ' ...
           'run_collect_by_group, channels = {''auto''}, or set auto_control = ' ...
           'false.'], auto_file);
end
held = who('-file', auto_file);
if ~all(ismember({'auto_4d', 'collected_mice', 'source_files'}, held))
    error(['run_per_mouse_values: %s does not record the registered volumes it was ' ...
           'read from: it is the older auto stack, of the registration before the ' ...
           'one behind the nano stack (collect_by_group''s help). Collect it again ' ...
           'with run_collect_by_group, channels = {''auto''}.'], auto_file);
end
S = load(auto_file, 'source_files');
source_files = S.source_files;
end

function A = group_auto(G, group_dir, brainMask, masks, apply_smoothing, smooth_sigma)
% One group's autofluorescence, mouse by mouse, in the clusters' box: the auto
% stack less the mouse's off-tissue level, NaN outside the tissue of its nano
% channel and smoothed, as the raw nano stack; after checking that the voxels a
% section reached are those of the nano stack.

% the auto stack, the collected nano stack and the background masks, read one
% mouse at a time
auto_file = fullfile(group_dir, 'auto_4d.mat');
raw_file = fullfile(group_dir, 'nano_4d.mat');
mask_file = fullfile(group_dir, 'nano_4d_normalized_bkgmask.mat');
M_auto = matfile(auto_file);
M_raw = matfile(raw_file);
M_mask = matfile(mask_file);
S_auto = load(auto_file, 'collected_mice', 'source_files');
S_mask = load(mask_file, 'current_mice');

% the mice's places in each
[is_collected, auto_idx] = ismember(G.names, S_auto.collected_mice);
if ~all(is_collected)
    error('run_per_mouse_values: %s not among the mice of %s (%s).', ...
        strjoin(G.names(~is_collected), ', '), auto_file, ...
        strjoin(S_auto.collected_mice, ', '));
end
raw_idx = collected_index(G.names, G.group, M_raw, 'nano_4d', raw_file);
[~, mask_idx] = ismember(G.names, S_mask.current_mice);

n_mice = numel(G.names);
A = struct();
A.group = G.group;
A.names = G.names;
A.source_files = S_auto.source_files(auto_idx);
A.background = nan(n_mice, 1);
A.reached_dice = nan(n_mice, 1);
A.box = cell(n_mice, 1);
for k = 1:n_mice
    t_mouse = tic;
    fprintf('%s autofluorescence, mouse %d of %d: %s\n', G.group, k, n_mice, G.names{k});
    auto = M_auto.auto_4d(:, :, :, auto_idx(k));
    raw = M_raw.nano_4d(:, :, :, raw_idx(k));
    bg_mask = M_mask.recomputed_bkg_mask_4d(:, :, :, mask_idx(k));

    % one registration: the voxels a section reached (not 0) must be the same in
    % both channels; another registration, two voxels away, would move the edge
    % of every section
    reached_nano = raw > 0;
    reached_auto = auto > 0;
    A.reached_dice(k) = 2 * nnz(reached_nano & reached_auto) / ...
        (nnz(reached_nano) + nnz(reached_auto));
    fprintf('  voxels a section reached: nano %d, auto %d, in both %d (Dice %.6f)\n', ...
        nnz(reached_nano), nnz(reached_auto), nnz(reached_nano & reached_auto), ...
        A.reached_dice(k));
    if A.reached_dice(k) < 0.999
        error(['run_per_mouse_values: the autofluorescence of %s (%s) does not reach ' ...
               'the voxels of its nano stack (Dice %.4f): not the same registration.'], ...
               G.names{k}, A.source_files{k}, A.reached_dice(k));
    end
    clear raw

    % less its off-tissue level, over the voxels the nano channel's level is
    % taken on, then NaN outside the tissue and smoothed, as the raw nano stack;
    % a voxel either channel did not reach is no tissue
    [A.background(k), plane_levels] = off_tissue_level(auto, bg_mask, brainMask);
    fprintf(['  off-tissue level %.1f (planes: 10th to 90th percentile %.1f to ' ...
             '%.1f)\n'], A.background(k), prctile(plane_levels, 10), ...
            prctile(plane_levels, 90));
    auto(~(reached_nano & reached_auto)) = NaN;
    clear reached_nano reached_auto
    auto = auto - A.background(k);
    auto = tissue_only(auto, bg_mask, brainMask, apply_smoothing, smooth_sigma);
    A.box{k} = auto(masks.box_ap, masks.box_dv, masks.box_ml);
    clear auto bg_mask

    fprintf('  done in %.1f min.\n', toc(t_mouse) / 60);
end
end

function check_cache_against_stack(G, group_dir, channel, cache_file)
% Stops unless the group's normalised stack still holds the mice and the lines
% of run_normalise_groups the cached maps were made from: after step 2 is run
% again, the cache would otherwise give the old maps without a word.

norm_file = fullfile(group_dir, [channel '_4d_normalized.mat']);
S_norm = load(norm_file, 'current_mice', 'norm_params');
if ~isfield(G, 'stack_norm_params') || ~isequal(G.stack_mice, S_norm.current_mice) ...
        || ~isequal(G.stack_norm_params, S_norm.norm_params)
    error(['run_per_mouse_values: the maps of %s in %s were not made from the ' ...
           'normalised stack now in %s (its mice or their lines differ, or the ' ...
           'cache predates the check). Set force_recompute_mice = true.'], G.group, ...
           cache_file, norm_file);
end
end

function check_raw_against_normalised(raw, vol, norm_params, mouse_name)
% Stops unless the normalised volume is the collected one through the mouse's
% line, (raw - intercept) / slope with NaN where the raw is 0, to the bit, as
% run_normalise_groups computes it: so the collected stack's mouse is the one
% normalised, and the raw values are those the test's maps come from.

expected = (raw - norm_params(2)) / norm_params(1);
expected(raw == 0) = NaN;
if ~isequaln(expected, vol)
    error(['run_per_mouse_values: the normalised volume of %s is not its collected ' ...
           'volume through its line, (raw - %g) / %g. Are the collected and the ' ...
           'normalised stack from the same run of run_collect_by_group?'], ...
           mouse_name, norm_params(2), norm_params(1));
end
fprintf('  normalised = (raw - %.1f) / %.4f, to the bit\n', norm_params(2), ...
    norm_params(1));
end

function [level, plane_levels] = off_tissue_level(raw, bg_mask, brainMask)
% A mouse's off-tissue level in its collected stack: the median over its
% background voxels outside the atlas brain that a section reached, over the
% whole stack, and plane by plane (NaN for a plane with fewer than 1000 of them).

% run_normalise_groups' mask marks, plane by plane, the voxels below the knee
% between background and tissue; outside the atlas brain they are the slide
% around the section, and a raw 0 is a voxel no section reached
is_off = bg_mask & ~brainMask & raw > 0;
level = double(median(raw(is_off)));

% plane by plane, for the spread along AP; 1000 voxels give a stable median
n_planes = size(raw, 1);
plane_levels = nan(n_planes, 1);
for z = 1:n_planes
    plane_raw = raw(z, :, :);
    plane_off = is_off(z, :, :);
    if nnz(plane_off) >= 1000
        plane_levels(z) = double(median(plane_raw(plane_off)));
    end
end
end

function sums = region_sums(lr_diff, lr_sum, masks)
% A mouse's folded maps summed in each region named in advance, over its voxels
% with a value: their number, the sums of |L - R| and of L + R, and the region's
% voxels; and the number and the sum of L + R in the isocortex.

n_regions = numel(masks.region_lin);
sums = struct();
sums.n = zeros(n_regions, 1);
sums.n_region = zeros(n_regions, 1);
sums.abs_diff = zeros(n_regions, 1);
sums.sum = zeros(n_regions, 1);
for r = 1:n_regions

    % the region's voxels where the mouse has a value (L + R misses the same ones)
    d = lr_diff(masks.region_lin{r});
    s = lr_sum(masks.region_lin{r});
    has_value = ~isnan(d);
    sums.n(r) = nnz(has_value);
    sums.n_region(r) = numel(d);
    sums.abs_diff(r) = sum(abs(d(has_value)), 'double');
    sums.sum(r) = sum(s(has_value), 'double');
end

% the isocortex
s = lr_sum(masks.iso_lin);
has_value = ~isnan(s);
sums.iso_n = nnz(has_value);
sums.iso_sum = sum(s(has_value), 'double');
end

% ===== Local functions: values and clusters =====

function T_mice = mouse_values(ctrl_mice, exp_mice, slope, intercept, common_factor, ...
    regions)
% One row per mouse, control mice first: its off-tissue level, and in each
% region its coverage and its values on the test's maps and on the raw stack.

groups = {ctrl_mice, exp_mice};
rows = cell(0, 1);
for g = 1:2
    G = groups{g};

    % the line of the test: the experimental group through the alignment, the
    % control group as it is
    if g == 2
        line_slope = slope;
        line_intercept = intercept;
    else
        line_slope = 1;
        line_intercept = 0;
    end

    for k = 1:numel(G.names)
        row = struct();
        row.mouse = G.names{k};
        row.group = G.group;
        row.background = G.background(k);
        row.background_q10 = prctile(G.background_planes(:, k), 10);
        row.background_q90 = prctile(G.background_planes(:, k), 90);

        % the isocortex's mean L + R on the test's maps and on the raw stack
        t = G.test{k};
        w = G.raw{k};
        iso_sum = (line_slope * t.iso_sum / t.iso_n + 2 * line_intercept) / common_factor;
        iso_sum_raw = w.iso_sum / w.iso_n;

        for r = 1:numel(regions)
            tag = region_tag(regions{r});

            % the means on the test's maps: the line, then the common factor
            abs_diff = line_slope * t.abs_diff(r) / t.n(r) / common_factor;
            lr_sum = (line_slope * t.sum(r) / t.n(r) + 2 * line_intercept) / ...
                common_factor;

            % the means on the raw stack
            abs_diff_raw = w.abs_diff(r) / w.n(r);
            lr_sum_raw = w.sum(r) / w.n(r);

            row.(['coverage_' tag]) = t.n(r) / t.n_region(r);
            row.(['abs_diff_' tag]) = abs_diff;
            row.(['sum_' tag]) = lr_sum;
            row.(['ai_' tag]) = asymmetry_index(abs_diff, lr_sum);
            row.(['sum_rel_' tag]) = lr_sum / iso_sum;
            row.(['abs_diff_raw_' tag]) = abs_diff_raw;
            row.(['sum_raw_' tag]) = lr_sum_raw;
            row.(['ai_raw_' tag]) = asymmetry_index(abs_diff_raw, lr_sum_raw);
            row.(['sum_rel_raw_' tag]) = lr_sum_raw / iso_sum_raw;
        end
        rows{end + 1, 1} = row; %#ok<AGROW>
    end
end
T_mice = struct2table([rows{:}]', 'AsArray', true);
end

function tag = region_tag(acronym)
% A region's acronym as a column name: lower case, '-' as '_'.

tag = lower(strrep(acronym, '-', '_'));
end

function sm = selection_splits(ctrl_mice, exp_mice, ctrl_auto, exp_auto, masks, ...
    perm_settings)
% The selection-matched test's clusters and readings, under every split of the
% pooled mice into groups of the observed sizes, the groups as they are first:
% the heaviest cluster of the box's region where |L - R| is higher in the
% split's experimental group (positive) and in its control group (negative),
% found with all the mice by region_permutation_test; their voxels and mass;
% and every mouse's ai on the test's maps, its ai_raw and, with the
% autofluorescence, its ai_auto in each (NaN without a cluster). The mice keep
% the alignment of their true groups, as the test keeps it under its splits: a
% split relabels the mice, not their scales.

n_ctrl = numel(ctrl_mice.names);
n_mice = n_ctrl + numel(exp_mice.names);
names = [ctrl_mice.names(:); exp_mice.names(:)];
boxes = [ctrl_mice.box(:); exp_mice.box(:)];
boxes_raw = [ctrl_mice.box_raw(:); exp_mice.box_raw(:)];
is_exp = (1:n_mice)' > n_ctrl;

% the alignment of all the mice, as run_group_differences fits it, their
% |L - R| on the band and the candidates
[~, ~, ~, slope, intercept, common_factor] = align_exp_to_ctrl(ctrl_mice.profiles, ...
    exp_mice.profiles);
stack = band_stack(boxes, is_exp, slope, intercept, common_factor, masks);
[geom, is_cand] = band_geometry(stack, masks, perm_settings);
stack = stack(is_cand, :);

% every mouse's folded maps in the box, read in every split's clusters: on the
% test's maps, on the raw stack and on the autofluorescence
readings = {'ai', 'ai_raw'};
if ~isempty(ctrl_auto)
    readings{end + 1} = 'ai_auto';
    boxes_auto = [ctrl_auto.box(:); exp_auto.box(:)];
end
n_readings = numel(readings);
lr_diff = cell(n_mice, n_readings);
lr_sum = cell(n_mice, n_readings);
for k = 1:n_mice
    [lr_diff{k, 1}, lr_sum{k, 1}] = aligned_lr(boxes{k}, is_exp(k), slope, ...
        intercept, common_factor);
    [lr_diff{k, 2}, lr_sum{k, 2}] = compute_lr_stats(boxes_raw{k});
    if n_readings == 3
        [lr_diff{k, 3}, lr_sum{k, 3}] = compute_lr_stats(boxes_auto{k});
    end
end

% every split as a logical row, true for the mice labelled control; nchoosek's
% first is mice 1 to n_ctrl, the groups as they are
ctrl_sets = nchoosek(1:n_mice, n_ctrl);
n_splits = size(ctrl_sets, 1);
in_ctrl = false(n_splits, n_mice);
for s = 1:n_splits
    in_ctrl(s, ctrl_sets(s, :)) = true;
end

% per split, its positive and negative cluster (columns 1 and 2); the values
% as split x mouse x cluster x reading
sm = struct();
sm.names = names;
sm.in_ctrl = in_ctrl;
sm.readings = readings;
sm.slope = slope;
sm.intercept = intercept;
sm.cluster_n = zeros(n_splits, 2);
sm.cluster_mass = zeros(n_splits, 2);
sm.values = nan(n_splits, n_mice, 2, n_readings);
sm.cluster_voxels = cell(1, 2);

% the test without its progress lines, which would come once per split
quiet_settings = perm_settings;
quiet_settings.quiet = true;
fprintf('Selection-matched test: the cluster search under each of %d splits...\n', ...
    n_splits);
t_start = tic;
for s = 1:n_splits

    % the split's groups, its control mice first
    order = [find(in_ctrl(s, :)), find(~in_ctrl(s, :))];
    perm = region_permutation_test({stack(:, order)}, n_ctrl, geom, quiet_settings);
    map = perm.maps{1};
    is_cluster = strcmp(perm.measure_names, 'cluster');
    sm.cluster_mass(s, :) = [map.null_pos(perm.splits.observed, 1, is_cluster), ...
        map.null_neg(perm.splits.observed, 1, is_cluster)];
    sm.cluster_n(s, :) = map.detail.cluster_n(1, :);

    % every mouse in each of its clusters
    for c = 1:2
        cluster_lin = map.detail.cluster_voxels{1, c};
        if isempty(cluster_lin)
            continue
        end
        for k = 1:n_mice
            for r = 1:n_readings
                sm.values(s, k, c, r) = cluster_ai(lr_diff{k, r}, lr_sum{k, r}, ...
                    cluster_lin);
            end
        end
        if s == 1
            sm.cluster_voxels{c} = cluster_lin;
        end
    end
    if s == 1 || mod(s, round(n_splits / 10)) == 0 || s == n_splits
        elapsed_min = toc(t_start) / 60;
        fprintf('  split %d of %d: %.1f min, about %.1f min left\n', s, n_splits, ...
            elapsed_min, elapsed_min / s * (n_splits - s));
    end
end
end

function check_mirror_splits(sm)
% With groups of equal size, a split and its mirror image (the groups swapped)
% give opposite t everywhere, so one's positive cluster is the other's negative
% one, with the same voxels, mass and values; a warning where they differ.

n_mice = size(sm.in_ctrl, 2);
if 2 * nnz(sm.in_ctrl(1, :)) ~= n_mice
    return
end
[~, mirror] = ismember(~sm.in_ctrl, sm.in_ctrl, 'rows');
same_mass = isequal(sm.cluster_mass(:, 1), sm.cluster_mass(mirror, 2));
same_values = isequaln(sm.values(:, :, 1, :), sm.values(mirror, :, 2, :));
if same_mass && same_values
    fprintf(['  every split''s positive cluster is its mirror image''s negative ' ...
             'one: the same mass and values.\n']);
else
    warning(['run_per_mouse_values: some splits'' positive clusters differ from ' ...
             'their mirror images'' negative ones (mass the same: %d, values the ' ...
             'same: %d).'], same_mass, same_values);
end
end

function [p_mass, p_mass_one] = cluster_mass_p(sm, step3_table, cluster_region)
% The p of the region's cluster mass over the splits, as the region test of
% run_group_differences takes it (the larger of the two signs, reached by the
% splits' larger one, the observed split included), checked against its table
% when it is there: the selection-matched test's splits are those of that test.
% Beside it, for the like of the one-sided p of the mice's values, the p of the
% positive cluster's mass alone (the experimental group higher); not a test of
% run_group_differences, whose p is the first.

larger = max(sm.cluster_mass, [], 2);
p_mass = mean(larger >= larger(1));
p_mass_one = mean(sm.cluster_mass(:, 1) >= sm.cluster_mass(1, 1));
fprintf(['  cluster mass over the %d splits: p %.4f (either sign), %.4f (the positive ' ...
         'cluster alone)\n'], numel(larger), p_mass, p_mass_one);
if ~exist(step3_table, 'file')
    return
end
T_step3 = readtable(step3_table);
step3_p = T_step3.lr_diff_cluster_p_perm(strcmp(T_step3.acronym, cluster_region));
if abs(p_mass - step3_p) <= 1e-9
    fprintf('  the p of run_group_differences'' region test: %.4f, the same.\n', ...
        step3_p);
else
    warning(['run_per_mouse_values: the cluster mass of %s gives p %.4f over the ' ...
             'selection-matched test''s splits; run_group_differences gave %.4f (%s).'], ...
             cluster_region, p_mass, step3_p, step3_table);
end
end

function [T_loo, loo] = leave_one_out(ctrl_mice, exp_mice, masks, perm_settings, ...
    step3_table, cluster_region)
% Fold 0, the comparison as it is, then each mouse left out in turn, control
% mice first: the fold's alignment, its heaviest cluster where |L - R| is higher
% in the experimental group and where it sits, and the left-out mouse's ai on
% the fold's maps and ai_raw in that cluster; the clusters' voxels, one cell per
% fold, fold 0 first.

n_ctrl = numel(ctrl_mice.names);
n_exp = numel(exp_mice.names);
n_folds = n_ctrl + n_exp;
loo = struct();
loo.ai = nan(n_folds, 1);
loo.ai_raw = nan(n_folds, 1);
loo.coverage = nan(n_folds, 1);
loo.cluster_n = zeros(n_folds, 1);
loo.cluster_voxels = cell(n_folds + 1, 1);
rows = cell(n_folds + 1, 1);
full_cluster = [];
for f = 0:n_folds

    % the mice kept, and the one left out (none in fold 0)
    keep_ctrl = true(1, n_ctrl);
    keep_exp = true(1, n_exp);
    if f >= 1 && f <= n_ctrl
        keep_ctrl(f) = false;
        left_out = ctrl_mice.names{f};
        left_group = ctrl_mice.group;
        left_box = ctrl_mice.box{f};
        left_box_raw = ctrl_mice.box_raw{f};
    elseif f > n_ctrl
        keep_exp(f - n_ctrl) = false;
        left_out = exp_mice.names{f - n_ctrl};
        left_group = exp_mice.group;
        left_box = exp_mice.box{f - n_ctrl};
        left_box_raw = exp_mice.box_raw{f - n_ctrl};
    else
        left_out = '';
        left_group = '';
    end
    fprintf('Fold %d of %d, left out: %s\n', f, n_folds, left_out);

    % the fold's alignment, fitted on its own mice
    [~, ~, ~, slope, intercept, common_factor] = align_exp_to_ctrl( ...
        ctrl_mice.profiles(:, keep_ctrl), exp_mice.profiles(:, keep_exp));

    % its heaviest cluster in the region where |L - R| is higher in the
    % experimental group, as the test finds it
    [cluster_lin, cluster] = fold_cluster(ctrl_mice.box(keep_ctrl), ...
        exp_mice.box(keep_exp), slope, intercept, common_factor, masks, perm_settings);
    if f == 0
        full_cluster = cluster_lin;
        check_against_step3(cluster, step3_table, cluster_region);
    end

    % the share of its voxels in the cluster of all the mice
    overlap = NaN;
    if ~isempty(cluster_lin)
        overlap = numel(intersect(cluster_lin, full_cluster)) / numel(cluster_lin);
    end

    % the left-out mouse's values in it: on the fold's maps, the experimental
    % group through the fold's line; and on its raw stack
    ai = NaN;
    ai_raw = NaN;
    coverage = NaN;
    if f >= 1 && ~isempty(cluster_lin)
        is_exp = f > n_ctrl;
        [lr_diff, lr_sum] = aligned_lr(left_box, is_exp, slope, intercept, ...
            common_factor);
        [ai, coverage] = cluster_ai(lr_diff, lr_sum, cluster_lin);
        [lr_diff, lr_sum] = compute_lr_stats(left_box_raw);
        ai_raw = cluster_ai(lr_diff, lr_sum, cluster_lin);
    end
    if f >= 1
        loo.ai(f) = ai;
        loo.ai_raw(f) = ai_raw;
        loo.coverage(f) = coverage;
        loo.cluster_n(f) = cluster.n;
    end
    if f >= 1 && isempty(cluster_lin)
        fprintf('  no cluster with a positive t in %s without %s.\n', cluster_region, ...
            left_out);
    end
    if f >= 1 && ~isempty(cluster_lin) && isnan(ai)
        fprintf(['  %s has no AI on the fold''s maps in its cluster: its mean L + R ' ...
                 'there is not above 0.\n'], left_out);
    end
    loo.cluster_voxels{f + 1} = cluster_lin;

    % where the cluster sits, in the planes and voxels of the folded maps
    place = cluster_place(cluster_lin, masks);

    rows{f + 1} = struct('fold', f, 'left_out', left_out, 'left_out_group', ...
        left_group, 'n_ctrl', nnz(keep_ctrl), 'n_exp', nnz(keep_exp), 'slope', slope, ...
        'intercept', intercept, 'cluster_n', cluster.n, 'cluster_mass', cluster.mass, ...
        'cluster_peak', cluster.peak, 'cluster_plane_first', place.plane_first, ...
        'cluster_plane_last', place.plane_last, 'cluster_centre_plane', ...
        place.centre(1), 'cluster_centre_dv', place.centre(2), 'cluster_centre_ml', ...
        place.centre(3), 'overlap_with_all_mice', overlap, 'left_out_ai', ai, ...
        'left_out_ai_raw', ai_raw, 'left_out_coverage', coverage);
end
T_loo = struct2table([rows{:}]', 'AsArray', true);
end

function place = cluster_place(cluster_lin, masks)
% A cluster's first and last plane and its centre (plane, DV, ML), in voxels of
% the folded maps (the planes of the volumes' crop, ML in the left hemisphere);
% NaN for no cluster.

place = struct('plane_first', NaN, 'plane_last', NaN, 'centre', nan(1, 3));
if isempty(cluster_lin)
    return
end

% from the box's grid back to the folded maps'
[ap, dv, ml] = ind2sub(masks.box_size, double(cluster_lin));
ap = ap + masks.box_ap(1) - 1;
dv = dv + masks.box_dv(1) - 1;
ml = ml + masks.box_ml_left(1) - 1;
place.plane_first = min(ap);
place.plane_last = max(ap);
place.centre = [mean(ap), mean(dv), mean(ml)];
end

function [lr_diff, lr_sum] = aligned_lr(box, is_exp, slope, intercept, common_factor)
% A mouse's box on the common scale, as run_group_differences puts a whole
% volume on it, the experimental group through the line first, then folded.

if is_exp
    vol = ((box .* slope) + intercept) ./ common_factor;
else
    vol = box ./ common_factor;
end
[lr_diff, lr_sum] = compute_lr_stats(vol);
end

function [cluster_lin, cluster] = fold_cluster(boxes_ctrl, boxes_exp, slope, ...
    intercept, common_factor, masks, perm_settings)
% The heaviest cluster in the box's region where |L - R| is higher in the
% experimental group (a positive t), for the mice given, by
% region_permutation_test on the groups as they are: its voxels (linear indices
% into the box's folded grid, empty for none), voxels, mass and peak.

% every mouse's |L - R| on the band, control mice first, and the candidates
n_ctrl = numel(boxes_ctrl);
is_exp = [false(n_ctrl, 1); true(numel(boxes_exp), 1)];
stack = band_stack([boxes_ctrl(:); boxes_exp(:)], is_exp, slope, intercept, ...
    common_factor, masks);
[geom, is_cand] = band_geometry(stack, masks, perm_settings);
perm = region_permutation_test({stack(is_cand, :)}, n_ctrl, geom, perm_settings);

% the positive cluster
map = perm.maps{1};
is_cluster = strcmp(perm.measure_names, 'cluster');
cluster_lin = map.detail.cluster_voxels{1, 1};
cluster = struct();
cluster.n = map.detail.cluster_n(1, 1);
cluster.mass = map.null_pos(perm.splits.observed, 1, is_cluster);
cluster.peak = map.detail.cluster_peak(1, 1);
fprintf('  heaviest cluster with a positive t: %d voxels, mass %.2f, peak %.2f\n', ...
    cluster.n, cluster.mass, cluster.peak);
end

function stack = band_stack(boxes, is_exp, slope, intercept, common_factor, masks)
% The mice's |L - R| on the band, one column per mouse in the order given, each
% on the common scale (the experimental mice through the line), as the test's
% stacks.

stack = zeros(numel(masks.band_lin), numel(boxes), 'single');
for k = 1:numel(boxes)
    lr_diff = aligned_lr(boxes{k}, is_exp(k), slope, intercept, common_factor);
    stack(:, k) = abs(lr_diff(masks.band_lin));
end
end

function [geom, is_cand] = band_geometry(stack, masks, perm_settings)
% The candidates, as the test takes them: voxels where at least twice
% min_mice_per_group mice have a value, the fewest with which a t is possible;
% the region's voxels labelled 1, the band around it 0, which the rolling median
% reads but no cluster takes. Which mice have a value does not depend on their
% groups, so a fold's candidates serve every split of its mice.

is_cand = sum(~isnan(stack), 2) >= 2 * perm_settings.min_mice_per_group;
geom = struct();
geom.grid_size = masks.box_size;
geom.cand_lin = masks.band_lin(is_cand);
geom.cand_region = uint16(masks.band_in_region(is_cand));
geom.n_regions = 1;
geom.voxel_mm = 0.01;
end

function relabelled = loo_relabelled(ctrl_mice, exp_mice, masks, perm_settings)
% The leave-one-out redone under every split of the pooled mice into groups of
% the observed sizes, the groups as they are first: for each split and each
% mouse left out, the heaviest cluster among the other mice where |L - R| is
% higher in the split's experimental group, and the left-out mouse's ai on the
% fold's maps and ai_raw in it (NaN without a cluster). Each fold keeps the
% alignment fitted on its mice's true groups, as the test keeps the alignment
% under its splits: a split relabels the mice, not their scales.

n_ctrl = numel(ctrl_mice.names);
n_mice = n_ctrl + numel(exp_mice.names);
names = [ctrl_mice.names(:); exp_mice.names(:)];
boxes = [ctrl_mice.box(:); exp_mice.box(:)];
boxes_raw = [ctrl_mice.box_raw(:); exp_mice.box_raw(:)];
profiles = [ctrl_mice.profiles, exp_mice.profiles];
is_exp = (1:n_mice)' > n_ctrl;

% every split as a logical row, true for the mice labelled control; nchoosek's
% first is mice 1 to n_ctrl, the groups as they are
ctrl_sets = nchoosek(1:n_mice, n_ctrl);
n_splits = size(ctrl_sets, 1);
in_ctrl = false(n_splits, n_mice);
for s = 1:n_splits
    in_ctrl(s, ctrl_sets(s, :)) = true;
end

relabelled = struct();
relabelled.names = names;
relabelled.in_ctrl = in_ctrl;
relabelled.ai = nan(n_splits, n_mice);
relabelled.ai_raw = nan(n_splits, n_mice);

% the test without its progress lines, which would come once per call
quiet_settings = perm_settings;
quiet_settings.quiet = true;
fprintf('Leave-one-out redone under each of %d splits, %d folds each...\n', n_splits, ...
    n_mice);
t_start = tic;
for f = 1:n_mice

    % the fold's mice, their alignment on their true groups, their |L - R| on
    % the band and the candidates
    kept = setdiff(1:n_mice, f);
    kept_exp = is_exp(kept);
    [~, ~, ~, slope, intercept, common_factor] = align_exp_to_ctrl( ...
        profiles(:, kept(~kept_exp)), profiles(:, kept(kept_exp)));
    stack = band_stack(boxes(kept), kept_exp, slope, intercept, common_factor, masks);
    [geom, is_cand] = band_geometry(stack, masks, perm_settings);
    stack = stack(is_cand, :);

    % the left-out mouse on the fold's maps and on its raw stack
    [lr_diff, lr_sum] = aligned_lr(boxes{f}, is_exp(f), slope, intercept, ...
        common_factor);
    [lr_diff_raw, lr_sum_raw] = compute_lr_stats(boxes_raw{f});

    % each split's groups among the fold's mice, its control mice first
    for s = 1:n_splits
        split_ctrl = in_ctrl(s, kept);
        order = [find(split_ctrl), find(~split_ctrl)];
        perm = region_permutation_test({stack(:, order)}, nnz(split_ctrl), geom, ...
            quiet_settings);
        cluster_lin = perm.maps{1}.detail.cluster_voxels{1, 1};
        if ~isempty(cluster_lin)
            relabelled.ai(s, f) = cluster_ai(lr_diff, lr_sum, cluster_lin);
            relabelled.ai_raw(s, f) = cluster_ai(lr_diff_raw, lr_sum_raw, cluster_lin);
        end
    end
    elapsed_min = toc(t_start) / 60;
    fprintf('  fold %d of %d (without %s): %.1f min, about %.1f min left\n', f, n_mice, ...
        names{f}, elapsed_min, elapsed_min / f * (n_mice - f));
end
end

function [ai, coverage] = cluster_ai(lr_diff, lr_sum, cluster_lin)
% A mouse's mean |L - R| over its mean L + R in a cluster, over the cluster's
% voxels where it has a value, and the share of the cluster's voxels it has.

d = lr_diff(cluster_lin);
s = lr_sum(cluster_lin);
has_value = ~isnan(d);
coverage = nnz(has_value) / numel(cluster_lin);
ai = asymmetry_index(mean(abs(d(has_value)), 'double'), mean(s(has_value), 'double'));
end

function ai = asymmetry_index(mean_abs_diff, mean_sum)
% Mean |L - R| over mean (L + R); NaN where the mean L + R is not above 0. On
% the test's maps the zero is the normalisation's, not the absence of signal (a
% mouse's off-tissue level lands between about -2,400 and +500 on the control
% group's scale), so a dim spot of a mouse can have a sum at or below 0.

if mean_sum > 0
    ai = mean_abs_diff / mean_sum;
else
    ai = NaN;
end
end

function check_against_step3(cluster, step3_table, cluster_region)
% The cluster of all the mice must be the region's cluster that
% run_group_differences found, with the same voxels and mass, when its table is
% there; a warning where it differs.

if ~exist(step3_table, 'file')
    fprintf('  no table of run_group_differences to check the cluster against (%s).\n', ...
        step3_table);
    return
end
T_step3 = readtable(step3_table);
row = strcmp(T_step3.acronym, cluster_region);
if T_step3.lr_diff_cluster_sign(row) <= 0
    fprintf(['  run_group_differences'' heaviest %s cluster is negative; the positive ' ...
             'one is not in its table.\n'], cluster_region);
    return
end
step3_n = T_step3.lr_diff_cluster_n(row);
step3_mass = T_step3.lr_diff_cluster_score(row);
if cluster.n == step3_n && abs(cluster.mass - step3_mass) <= 1e-6 * step3_mass
    fprintf(['  all the mice give the cluster of run_group_differences: %d voxels, ' ...
             'mass %.2f.\n'], cluster.n, cluster.mass);
else
    warning(['run_per_mouse_values: all the mice give a %s cluster of %d voxels, mass ' ...
             '%.2f; run_group_differences found %d voxels, mass %.2f (%s).'], ...
             cluster_region, cluster.n, cluster.mass, step3_n, step3_mass, step3_table);
end
end

% ===== Local functions: statistics =====

function values = value_list(regions, cluster_region, sm)
% The values tested and drawn, in the figures' order: the selection-matched
% values in the cluster of all the mice; per region, the ai and the relative
% L + R on the raw stack and the test's maps; then the leave-one-out's. Each
% with its column, its kind (which test its p is), its title and its axis label.

values = struct('name', {}, 'kind', {}, 'region', {}, 'title', {}, 'label', {});
ai_label = 'mean |L - R| / mean (L + R)';
sum_label = 'L + R relative to isocortex';

% the selection-matched values, raw stack first
if ~isempty(sm)
    readings = {
        'ai_raw',  'AI, raw stack'
        'ai',      'AI, test maps'
        'ai_auto', 'autofluorescence AI'
        };
    for k = 1:size(readings, 1)
        if ismember(readings{k, 1}, sm.readings)
            values(end + 1) = struct('name', ['sm_' readings{k, 1}], 'kind', ...
                'selection', 'region', cluster_region, 'title', readings{k, 2}, ...
                'label', ai_label); %#ok<AGROW>
        end
    end
end

% the regions' values, raw stack first
kinds = {
    'ai_raw_',      'AI, raw stack',                ai_label
    'ai_',          'AI, test maps',                ai_label
    'sum_rel_raw_', 'L + R / isocortex, raw stack', sum_label
    'sum_rel_',     'L + R / isocortex, test maps', sum_label
    };
for r = 1:numel(regions)
    for k = 1:size(kinds, 1)
        values(end + 1) = struct('name', [kinds{k, 1} region_tag(regions{r})], ...
            'kind', 'region', 'region', regions{r}, 'title', kinds{k, 2}, ...
            'label', kinds{k, 3}); %#ok<AGROW>
    end
end

% the leave-one-out's, raw stack first
values(end + 1) = struct('name', 'loo_ai_raw', 'kind', 'loo', 'region', ...
    cluster_region, 'title', 'AI in the leave-one-out cluster, raw stack', 'label', ...
    ai_label);
values(end + 1) = struct('name', 'loo_ai', 'kind', 'loo', 'region', cluster_region, ...
    'title', 'AI in the leave-one-out cluster, test maps', 'label', ai_label);
end

function T_stats = value_statistics(T_mice, values, ctrl_type, exp_type, sm, relabelled)
% Per value, over the mice with a value: each group's mean and SEM, the
% difference of the means, the p of its test (one-sided for the experimental
% group higher, and two-sided), the Welch t and Hedges' g. A region's test is
% the exact permutation of its values. A cluster value's test redoes the
% cluster under every split (the selection-matched test, the leave-one-out
% redone); the permutation of its values, which keeps the cluster the true
% groups chose, the Welch p and Hedges' g ignore that choice, and are kept in
% columns named anti-conservative and inflated, the plain columns left NaN.

is_ctrl = strcmp(T_mice.group, ctrl_type);
rows = cell(numel(values), 1);
for v = 1:numel(values)
    x = T_mice.(values(v).name);
    x_ctrl = x(is_ctrl & ~isnan(x));
    x_exp = x(~is_ctrl & ~isnan(x));

    % the permutation of the values, and the Welch t and Hedges' g for reference
    [p_two_shuffled, p_one_shuffled, n_splits] = exact_permutation(x_ctrl, x_exp);
    [~, p_welch, ~, welch] = ttest2(x_exp, x_ctrl, 'Vartype', 'unequal');
    reference = struct('t', welch.tstat, 'df', welch.df, 'p', p_welch, 'g', ...
        hedges_g(x_ctrl, x_exp));

    % the value's test
    test = struct('p_one', p_one_shuffled, 'p_two', p_two_shuffled, 'n_splits', ...
        n_splits, 'n_empty', 0, 'n_empty_two', 0, 'p_one_without', NaN, ...
        'p_two_without', NaN, 'down', NaN, 'either', NaN, 'either_cluster', '', ...
        'null_median', NaN);
    anticonservative = [NaN NaN];
    inflated = [NaN NaN];
    switch values(v).kind
        case 'region'
            description = 'values permuted over the splits';
        case 'selection'
            description = ['selection-matched: the cluster search redone in every ' ...
                'split, every mouse read in its cluster; a split without a cluster ' ...
                'or a group without a value counts as 0; two-sided: the search in ' ...
                'either direction. The shuffled and Welch p and Hedges g ignore the ' ...
                'selection (anti-conservative, inflated)'];
            reading = strrep(values(v).name, 'sm_', '');
            test = selection_p(sm, reading);
        case 'loo'
            description = ['leave-one-out redone in every split; splits without a ' ...
                'value in a group left out; two-sided: |difference|, every fold''s ' ...
                'cluster on the split''s experimental-higher side. The shuffled and ' ...
                'Welch p and Hedges g ignore the selection (anti-conservative, ' ...
                'inflated)'];
            test = relabelled_test(relabelled, strrep(values(v).name, 'loo_', ''));
    end
    if ~strcmp(values(v).kind, 'region')
        anticonservative = [p_one_shuffled, p_two_shuffled];
        inflated = [reference.p, reference.g];
        reference = struct('t', NaN, 'df', NaN, 'p', NaN, 'g', NaN);
    end

    rows{v} = struct('value', values(v).name, 'region', values(v).region, ...
        'title', values(v).title, 'test', description, 'ctrl', ctrl_type, ...
        'exp', exp_type, 'n_ctrl', numel(x_ctrl), 'n_exp', numel(x_exp), ...
        'mean_ctrl', mean(x_ctrl), 'sem_ctrl', std(x_ctrl) / sqrt(numel(x_ctrl)), ...
        'mean_exp', mean(x_exp), 'sem_exp', std(x_exp) / sqrt(numel(x_exp)), ...
        'difference', mean(x_exp) - mean(x_ctrl), ...
        'difference_negative_cluster', test.down, ...
        'difference_either_direction', test.either, ...
        'either_direction_cluster', test.either_cluster, ...
        'null_median_difference', test.null_median, ...
        'p_exp_higher', test.p_one, 'p_two_sided', test.p_two, ...
        'n_splits', test.n_splits, 'n_splits_empty', test.n_empty, ...
        'n_splits_empty_two_sided', test.n_empty_two, ...
        'p_exp_higher_without_empty', test.p_one_without, ...
        'p_two_sided_without_empty', test.p_two_without, ...
        'p_shuffled_exp_higher_anticonservative', anticonservative(1), ...
        'p_shuffled_two_sided_anticonservative', anticonservative(2), ...
        'welch_t', reference.t, 'welch_df', reference.df, 'welch_p', reference.p, ...
        'hedges_g', reference.g, 'welch_p_anticonservative', inflated(1), ...
        'hedges_g_inflated', inflated(2));
end
T_stats = struct2table([rows{:}]', 'AsArray', true);
end

function [p_two, p_one, n_splits] = exact_permutation(x_ctrl, x_exp)
% The exact permutation p of the difference of the group means, experimental
% minus control: the share of all the splits of the pooled mice into groups of
% the observed sizes whose difference reaches the observed one, in either
% direction (two-sided) or upwards (one-sided), the observed split included.

pooled = [x_ctrl(:); x_exp(:)];
n_mice = numel(pooled);
n_ctrl = numel(x_ctrl);
if n_ctrl < 2 || n_mice - n_ctrl < 2
    p_two = NaN;
    p_one = NaN;
    n_splits = 0;
    return
end

% every split; nchoosek's first is mice 1 to n_ctrl as control, the groups as
% they are
ctrl_sets = nchoosek(1:n_mice, n_ctrl);
n_splits = size(ctrl_sets, 1);
differences = zeros(n_splits, 1);
for s = 1:n_splits
    in_ctrl = false(n_mice, 1);
    in_ctrl(ctrl_sets(s, :)) = true;
    differences(s) = mean(pooled(~in_ctrl)) - mean(pooled(in_ctrl));
end
observed = differences(1);
p_two = mean(abs(differences) >= abs(observed));
p_one = mean(differences >= observed);
end

function test = selection_p(sm, reading)
% The selection-matched test of one reading: per split the difference of the
% group means in its positive cluster (experimental minus control: up) and in
% its negative cluster (control minus experimental: down), each over the mice
% with a value; the one-sided p, the share of the splits whose up reaches the
% observed one, and the two-sided, whose larger of up and down reaches the
% observed one, the observed split included. A split without the cluster, or
% with a group without a value, gives 0; the p with those splits left out
% beside. Also which cluster gives the observed split's larger difference, and
% the null's median difference over the splits with a value: every split's
% cluster separates its own groups, so that median is not 0.

r = strcmp(sm.readings, reading);
n_splits = size(sm.in_ctrl, 1);
up = nan(n_splits, 1);
down = nan(n_splits, 1);
for s = 1:n_splits
    labelled_ctrl = sm.in_ctrl(s, :);
    up(s) = group_difference(sm.values(s, :, 1, r), ~labelled_ctrl, labelled_ctrl);
    down(s) = group_difference(sm.values(s, :, 2, r), labelled_ctrl, ~labelled_ctrl);
end

% the empty splits as 0; the larger of the two directions
is_empty = isnan(up);
is_empty_two = isnan(up) & isnan(down);
up_counted = up;
up_counted(is_empty) = 0;
down_counted = down;
down_counted(isnan(down)) = 0;
either = max(up_counted, down_counted);

test = struct();
test.up = up;
test.down = down(1);
test.either = either(1);
if down_counted(1) > up_counted(1)
    test.either_cluster = 'ctrl higher';
else
    test.either_cluster = 'exp higher';
end
test.null_median = median(up(~is_empty));
test.null_up = up_counted;
test.n_splits = n_splits;
test.n_empty = nnz(is_empty);
test.n_empty_two = nnz(is_empty_two);
test.p_one = mean(up_counted >= up_counted(1));
test.p_two = mean(either >= either(1));
test.p_one_without = mean(up(~is_empty) >= up(1));
test.p_two_without = mean(either(~is_empty_two) >= either(1));
if is_empty(1)
    test.p_one_without = NaN;
end
end

function d = group_difference(x, is_a, is_b)
% The mean of x over the mice of a with a value minus that over the mice of b;
% NaN when either has none.

has_value = ~isnan(x(:)');
if ~any(is_a & has_value) || ~any(is_b & has_value)
    d = NaN;
    return
end
d = mean(x(is_a & has_value)) - mean(x(is_b & has_value));
end

function test = relabelled_test(relabelled, reading)
% The test of a leave-one-out value, from the leave-one-out redone under every
% split; NaN when it was not redone.

test = struct('p_one', NaN, 'p_two', NaN, 'n_splits', 0, 'n_empty', 0, ...
    'n_empty_two', NaN, 'p_one_without', NaN, 'p_two_without', NaN, 'down', NaN, ...
    'either', NaN, 'either_cluster', '', 'null_median', NaN);
if isempty(relabelled)
    return
end
[test.p_one, test.p_two, test.n_splits] = relabelled_p(relabelled.(reading), ...
    relabelled.in_ctrl);
test.n_empty = size(relabelled.in_ctrl, 1) - test.n_splits;
test.p_one_without = test.p_one;
test.p_two_without = test.p_two;
end

function [p_one, p_two, n_splits] = relabelled_p(values, in_ctrl)
% The p of the difference of the group means (experimental minus control) over
% the splits of the redone leave-one-out, each split with its own values (one
% row per split, the groups as they are first) and over the mice with a value:
% the share of the splits with a difference that reach the groups' own,
% upwards (one-sided) or in either direction (two-sided).

n_splits_all = size(values, 1);
differences = nan(n_splits_all, 1);
for s = 1:n_splits_all
    has_value = ~isnan(values(s, :));
    differences(s) = mean(values(s, ~in_ctrl(s, :) & has_value)) - ...
        mean(values(s, in_ctrl(s, :) & has_value));
end

% a split whose group has no mouse with a value gives no difference
has_difference = ~isnan(differences);
n_splits = nnz(has_difference);
observed = differences(1);
p_one = mean(differences(has_difference) >= observed);
p_two = mean(abs(differences(has_difference)) >= abs(observed));
end

function d = swapped_difference(relabelled, reading)
% The difference of the group means (the labelled experimental group minus the
% labelled control group) of a leave-one-out value under the split that swaps
% the groups, which only groups of equal size have; NaN without it.

d = NaN;
if isempty(relabelled)
    return
end
in_ctrl = relabelled.in_ctrl;
mirror = find(all(in_ctrl == ~in_ctrl(1, :), 2), 1);
if isempty(mirror)
    return
end
d = group_difference(relabelled.(reading)(mirror, :), ~in_ctrl(mirror, :), ...
    in_ctrl(mirror, :));
end

function [r_pooled, r_within] = channel_correlation(x, y, is_ctrl)
% The Pearson r of two values over the mice with both, pooled and within the
% groups (each group's mean taken out of both first, so a difference between
% the groups does not make it).

has_both = ~isnan(x) & ~isnan(y);
x = x(has_both);
y = y(has_both);
is_ctrl = is_ctrl(has_both);
r = corrcoef(x, y);
r_pooled = r(1, 2);
for g = [true false]
    in_group = is_ctrl == g;
    x(in_group) = x(in_group) - mean(x(in_group));
    y(in_group) = y(in_group) - mean(y(in_group));
end
r = corrcoef(x, y);
r_within = r(1, 2);
end

function g = hedges_g(x_ctrl, x_exp)
% Hedges' g of experimental minus control: the difference of the means over the
% pooled SD, with the small-sample correction.

n_ctrl = numel(x_ctrl);
n_exp = numel(x_exp);
sd_pooled = sqrt(((n_ctrl - 1) * var(x_ctrl) + (n_exp - 1) * var(x_exp)) / ...
    (n_ctrl + n_exp - 2));
d = (mean(x_exp) - mean(x_ctrl)) / sd_pooled;
g = d * (1 - 3 / (4 * (n_ctrl + n_exp) - 9));
end

function print_summary(T_mice, T_stats, T_loo, values, sm_place)
% The values of every mouse, the cluster of all the mice, the leave-one-out
% folds and each value's test.

fprintf('\nPer-mouse values:\n');
columns = [{'mouse', 'group'}, {values.name}];
disp(T_mice(:, columns));
if ~isempty(sm_place)
    fprintf(['Cluster of all the mice: %d voxels, mass %.2f (its region test''s p ' ...
             '%.4f, the positive cluster alone %.4f), planes %d to %d, centre at plane ' ...
             '%.1f, DV %.1f, ML %.1f.\n'], sm_place.n, sm_place.mass, sm_place.mass_p, ...
            sm_place.mass_p_one, sm_place.plane_first, sm_place.plane_last, ...
            sm_place.centre);
end
fprintf('Leave-one-out folds:\n');
disp(T_loo);
fprintf('Group tests:\n');
disp(T_stats(:, {'value', 'n_ctrl', 'n_exp', 'mean_ctrl', 'mean_exp', 'difference', ...
    'null_median_difference', 'p_exp_higher', 'p_two_sided', 'n_splits', ...
    'n_splits_empty', 'p_exp_higher_without_empty', 'p_two_sided_without_empty', ...
    'difference_either_direction', 'either_direction_cluster', 'welch_p', 'hedges_g'}));
end

% ===== Local functions: figures =====

function plot_per_mouse_values(T_mice, T_stats, sm, values, regions, comparison, ...
    file_tag, comp_out_dir)
% The selection-matched values in the first row (one dot per mouse with its
% name, each group's mean and SEM, the p of the test), with the splits of the
% test under both its statistic and step 3's, then a row per region named in
% advance, then the notes: which p tests what.

% four columns, not the three of docs/STYLE.md, so a row is one region; the
% notes in the last row
n_cols_grid = 4;
n_rows_grid = numel(regions) + 2;
fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'Normalized', ...
    'Position', [0 0 1 1]);

% the selection-matched values, then the splits under both statistics
is_selection = strcmp({values.kind}, 'selection');
panel = 0;
for v = find(is_selection)
    panel = panel + 1;
    subplot(n_rows_grid, n_cols_grid, panel);
    draw_value(T_mice, T_stats, values(v), comparison, ...
        sprintf('%s cluster of all %d mice (step 3''s): %s', values(v).region, ...
        height(T_mice), values(v).title));
end
if ~isempty(sm)
    subplot(n_rows_grid, n_cols_grid, n_cols_grid);
    draw_selection_splits(sm, comparison);
end

% a row per region, after the first: the whole region, nothing chosen
is_region = strcmp({values.kind}, 'region');
region_values = values(is_region);
for v = 1:numel(region_values)
    subplot(n_rows_grid, n_cols_grid, n_cols_grid + v);
    draw_value(T_mice, T_stats, region_values(v), comparison, ...
        sprintf('Whole %s (no selection): %s', region_values(v).region, ...
        region_values(v).title));
end

% the notes, in two columns of the last row
[notes_left, notes_right] = notes_lines(T_mice, T_stats, regions, comparison, ...
    nnz(is_selection | is_region), file_tag);
last_row = (n_rows_grid - 1) * n_cols_grid;
subplot(n_rows_grid, n_cols_grid, last_row + [1 2]);
axis off;
text(0, 1.1, notes_left, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
    'FontSize', 8, 'Interpreter', 'none');
subplot(n_rows_grid, n_cols_grid, last_row + [3 4]);
axis off;
text(0, 1.1, notes_right, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
    'FontSize', 8, 'Interpreter', 'none');

% the title: the comparison, its alignment and the tests
title_line = sprintf('Per-mouse values - %s', strrep(file_tag, '_', ' '));
align_line = sprintf(['%s aligned onto %s by the line %.3f x %+.1f; first row: the ' ...
                      'cluster of step 3''s test, its p redoing the selection in every ' ...
                      'split; other rows: whole regions, the values permuted'], ...
                      comparison.exp_type, comparison.ctrl_type, comparison.slope, ...
                      comparison.intercept);
sgtitle({title_line, ['\rm\fontsize{11}' align_line]}, 'FontSize', 14, ...
    'FontWeight', 'bold');

% save it
saveas(fig, fullfile(comp_out_dir, ['Per_Mouse_Values_' file_tag '.fig']));
exportgraphics(fig, fullfile(comp_out_dir, ['Per_Mouse_Values_' file_tag '.png']), ...
    'Resolution', 300);
end

function plot_loo_values(T_mice, T_stats, T_loo, values, comparison, relabelled, ...
    file_tag, comp_out_dir)
% The leave-one-out's values, one panel each, with the p of the leave-one-out
% redone under every split, and its notes: a supplementary figure, each mouse
% read out of sample.

is_loo = find(strcmp({values.kind}, 'loo'));
fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'Normalized', ...
    'Position', [0 0 1 0.5]);
for k = 1:numel(is_loo)
    subplot(1, 4, k);
    draw_value(T_mice, T_stats, values(is_loo(k)), comparison, ...
        sprintf('%s: %s', values(is_loo(k)).region, values(is_loo(k)).title));
end

% the notes
subplot(1, 4, [3 4]);
axis off;
text(0, 1, loo_notes_lines(T_stats, T_loo, comparison, relabelled, file_tag), ...
    'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 8, ...
    'Interpreter', 'none');

title_line = sprintf('Per-mouse values, leave-one-out (each mouse out of sample) - %s', ...
    strrep(file_tag, '_', ' '));
sgtitle(title_line, 'FontSize', 14, 'FontWeight', 'bold');
saveas(fig, fullfile(comp_out_dir, ['Per_Mouse_LOO_' file_tag '.fig']));
exportgraphics(fig, fullfile(comp_out_dir, ['Per_Mouse_LOO_' file_tag '.png']), ...
    'Resolution', 300);
end

function draw_value(T_mice, T_stats, value, comparison, panel_title)
% One value's panel: the dots, names, means and SEMs, its title, and the p of
% its test under the title.

is_ctrl = strcmp(T_mice.group, comparison.ctrl_type);
short_names = cellfun(@(n) strtok(n, '_'), T_mice.mouse, 'UniformOutput', false);
x = T_mice.(value.name);
draw_value_panel(x(is_ctrl), x(~is_ctrl), short_names(is_ctrl), ...
    short_names(~is_ctrl), comparison.ctrl_type, comparison.exp_type);
stat = T_stats(strcmp(T_stats.value, value.name), :);
title(panel_title, 'FontSize', 10, 'Interpreter', 'none');
subtitle(stat_lines(stat, value.kind, comparison), 'FontSize', 8, 'Interpreter', 'none');
ylabel(value.label, 'FontSize', 9);
end

function draw_value_panel(x_ctrl, x_exp, names_ctrl, names_exp, ctrl_type, exp_type)
% Both groups' mice as dots with their names, and each group's mean and SEM to
% the left of its dots.

hold on;
box on;
grid on;

% the values' range with a margin on each side, so no dot or name sits on the
% frame; the gap between two names a fifteenth of it, about a line of 7 points
% in a panel of the grid
all_values = [x_ctrl(:); x_exp(:)];
all_values = all_values(~isnan(all_values));
if isempty(all_values)
    all_values = 0;
end
value_range = max(all_values) - min(all_values);
if ~(value_range > 0)
    value_range = 1;
end
y_limits = [min(all_values) - 0.12 * value_range, max(all_values) + 0.12 * value_range];
label_gap = diff(y_limits) / 15;

draw_group(1, x_ctrl, names_ctrl, sep_palette('control'), ...
    sep_palette('control_mean'), label_gap);
draw_group(2, x_exp, names_exp, sep_palette('experimental'), ...
    sep_palette('experimental_mean'), label_gap);
ylim(y_limits);
xlim([0.5 2.8]);
xticks([1 2]);
xticklabels({ctrl_type, exp_type});
set(gca, 'FontSize', 9);
end

function draw_group(x_pos, values, names, dot_colour, mean_colour, label_gap)
% One group: its dots at x_pos, their names to the right, kept label_gap apart
% and each joined to its dot by a thin line, and its mean and SEM just to the
% left.

has_value = ~isnan(values);
values = values(has_value);
names = names(has_value);
if isempty(values)
    return
end
scatter(x_pos * ones(numel(values), 1), values, 40, dot_colour, 'filled', ...
    'MarkerFaceAlpha', 0.85, 'MarkerEdgeColor', 'none');

% the mean as a short line, the SEM as a bar, left of the dots
group_mean = mean(values);
group_sem = std(values) / sqrt(numel(values));
plot(x_pos - 0.2 + [-0.1 0.1], [group_mean group_mean], '-', 'Color', mean_colour, ...
    'LineWidth', 2);
errorbar(x_pos - 0.2, group_mean, group_sem, 'Color', mean_colour, 'LineWidth', 1.2, ...
    'CapSize', 6);

% the names, moved apart where two dots are close, each joined to its dot
label_y = spread_labels(values, label_gap);
for k = 1:numel(values)
    plot(x_pos + [0.03 0.075], [values(k) label_y(k)], '-', 'Color', mean_colour, ...
        'LineWidth', 0.4);
end
text(x_pos + 0.08 * ones(numel(values), 1), label_y, names, 'FontSize', 7, ...
    'Color', mean_colour, 'Interpreter', 'none', 'VerticalAlignment', 'middle');
end

function label_y = spread_labels(y, min_gap)
% Heights for the names of dots at heights y, at least min_gap apart, each as
% close to its dot as that allows, the whole set centred on the dots.

[sorted_y, order] = sort(y(:));
placed = sorted_y;
for i = 2:numel(placed)
    placed(i) = max(placed(i), placed(i - 1) + min_gap);
end
placed = placed - (mean(placed) - mean(sorted_y));
label_y = zeros(numel(y), 1);
label_y(order) = placed;
end

function draw_selection_splits(sm, comparison)
% The splits of the selection-matched test under both statistics taken from
% their clusters: step 3's, the mass of the split's cluster where its
% experimental group is higher (x, log), and this test's, the difference of
% the group means of the raw AI in it (y); the splits without that cluster in
% a strip at the left, at the 0 they count as; the observed split marked, with
% dashed lines at its two values; how many splits reach each, and both.

test = selection_p(sm, 'ai_raw');
mass = sm.cluster_mass(:, 1);
difference = test.null_up;
has_cluster = mass > 0;
reach_mass = mass >= mass(1);
reach_difference = difference >= difference(1);

% the log axis over the clusters' masses, the strip half a decade below them
first_decade = floor(log10(min(mass(has_cluster))));
last_decade = ceil(log10(max(mass)));
strip_x = 10 ^ (first_decade - 0.5);

hold on;
box on;
grid on;
set(gca, 'XScale', 'log', 'XMinorGrid', 'off', 'XMinorTick', 'off');
scatter(mass(has_cluster), difference(has_cluster), 14, sep_palette('paired_lines'), ...
    'filled', 'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'none');
scatter(strip_x, 0, 30, sep_palette('paired_lines'), 's', 'filled');
text(strip_x * 1.4, 0, sprintf('%d', nnz(~has_cluster)), 'FontSize', 7, ...
    'VerticalAlignment', 'middle', 'Color', [0.3 0.3 0.3]);
xline(10 ^ (first_decade - 0.25), ':', 'Color', [0.5 0.5 0.5]);
xline(mass(1), '--', 'Color', sep_palette('experimental_mean'), 'LineWidth', 0.8);
yline(difference(1), '--', 'Color', sep_palette('experimental_mean'), 'LineWidth', 0.8);
scatter(mass(1), difference(1), 50, sep_palette('experimental_mean'), 'filled');

% the strip's tick, then a tick per decade, as powers of ten so they stay level
xlim(10 .^ [first_decade - 0.75, last_decade]);
powers = first_decade:last_decade;
xticks([strip_x, 10 .^ powers]);
xticklabels([{'none'}, arrayfun(@(k) sprintf('10^{%d}', k), powers, ...
    'UniformOutput', false)]);
set(gca, 'FontSize', 9, 'TickLabelInterpreter', 'tex', 'XTickLabelRotation', 0);
xlabel('mass of the split''s cluster (step 3''s score)', 'FontSize', 9);
ylabel(sprintf('mean raw AI, ''%s'' - ''%s''', comparison.exp_type, ...
    comparison.ctrl_type), 'FontSize', 9);

% the counts, the p of each statistic beside them; the observed split in red
n_splits = numel(mass);
title(sprintf('The same %d splits: cluster mass and AI difference', n_splits), ...
    'FontSize', 10);
subtitle({
    sprintf('mass >= observed: %d (p %.3f; either sign %.3f)', nnz(reach_mass), ...
        comparison.sm_place.mass_p_one, comparison.sm_place.mass_p)
    sprintf('AI difference >= observed: %d (p %.3f); both: %d', ...
        nnz(reach_difference), test.p_one, nnz(reach_mass & reach_difference))
    }, 'FontSize', 8);
end

function lines = stat_lines(stat, kind, comparison)
% The p of a value under its panel's title, the one named before any number
% first: for a region, the permutation of its values, the Welch p and Hedges'
% g beside it; for the selection-matched values, their test with the p
% without the splits that have no cluster in brackets, the null's median
% difference against the observed one, and, when the larger difference of the
% search in either direction is in the cluster where the control group is
% higher, that difference; for the leave-one-out's, the leave-one-out redone
% under every split.

switch kind
    case 'region'
        lines = {sprintf('perm p %s', p_pair(stat.p_exp_higher, stat.p_two_sided, ...
            comparison, 'two-sided')), sprintf('Welch p %.3f, Hedges g %.2f, n %d and %d', ...
            stat.welch_p, stat.hedges_g, stat.n_ctrl, stat.n_exp)};
    case 'selection'
        without = [stat.p_exp_higher_without_empty, stat.p_two_sided_without_empty];
        if ~comparison.direction_named
            without = flip(without);
        end
        lines = {
            sprintf('p %s [%.3f, %.3f]', p_pair(stat.p_exp_higher, stat.p_two_sided, ...
                comparison, 'either direction'), without)
            sprintf('null median %+.3f, observed %+.3f; n %d and %d', ...
                stat.null_median_difference, stat.difference, stat.n_ctrl, stat.n_exp)
            };
        if strcmp(stat.either_direction_cluster, 'ctrl higher')
            lines{end + 1} = sprintf('either direction: %s-higher cluster (%d voxels), %+.3f', ...
                comparison.ctrl_type, comparison.sm_place.n_negative, ...
                stat.difference_either_direction);
        end
    case 'loo'
        lines = {
            sprintf('p %s', p_pair(stat.p_exp_higher, stat.p_two_sided, comparison, ...
                'two-sided'))
            sprintf('every fold redone in each split; %d of %d with a value', ...
                stat.n_splits, stat.n_splits + stat.n_splits_empty)
            sprintf('n %d and %d', stat.n_ctrl, stat.n_exp)
            };
end
end

function text_p = p_pair(p_one, p_two, comparison, two_sided_name)
% A one-sided p (experimental group higher) and a two-sided one, under the name
% given: the one-sided first when its direction was named before any number,
% the two-sided first otherwise.

one_sided = sprintf('%.3f %s > %s', p_one, comparison.exp_type, comparison.ctrl_type);
two_sided = sprintf('%.3f %s', p_two, two_sided_name);
if comparison.direction_named
    text_p = [one_sided ', ' two_sided];
else
    text_p = [two_sided ', ' one_sided];
end
end

function [left, right] = notes_lines(T_mice, T_stats, regions, comparison, n_values, ...
    file_tag)
% The notes of the main figure, in two columns: which p tests what, then what
% the values are, the mice with little of a region and the comparison's caveat;
% the settings and the cluster from the run, not typed.

ctrl_type = comparison.ctrl_type;
exp_type = comparison.exp_type;
perm_settings = comparison.perm_settings;
place = comparison.sm_place;
short_name = @(names) cellfun(@(n) strtok(n, '_'), names, 'UniformOutput', false);
stat_of = @(name) T_stats(strcmp(T_stats.value, name), :);

left = {'Which p tests what'};
if ~isempty(place)
    sm_stats = stat_of('sm_ai_raw');
    left = [left; {
        sprintf(['First row: the heaviest cluster of %s where |L - R| is higher in %s ' ...
                 '(a positive t'], comparison.cluster_region, exp_type)
        sprintf(['  at p < %g; median over +/- %d planes; %d-connected; a t where ' ...
                 'each group has %d mice),'], perm_settings.cluster_p, ...
                 perm_settings.slab_range, perm_settings.cluster_connectivity, ...
                 perm_settings.min_mice_per_group)
        sprintf(['  found with all the mice as step 3: %d voxels, planes %d to %d. ' ...
                 'Each mouse''s AI read in it.'], place.n, place.plane_first, ...
                 place.plane_last)
        sprintf(['  p: the same search under each of the %d splits of the mice, every ' ...
                 'mouse read in that'], sm_stats.n_splits)
        '  split''s cluster; the share of the splits whose difference of the group means'
        sprintf(['  reaches the observed one (the %d without a cluster count as 0; in ' ...
                 'brackets, the p'], sm_stats.n_splits_empty)
        sprintf(['  without them). Either direction: per split the larger of that ' ...
                 'difference and of %s'], ctrl_type)
        sprintf('  minus %s in its cluster where %s is higher.', exp_type, ctrl_type)
        ['Same search, same splits as step 3''s region test; only the number taken ' ...
         'from each']
        '  split''s cluster differs. Step 3 takes its mass, which grows with its extent'
        sprintf('  (here p %s; step 3''s p: either sign);', ...
            p_pair(place.mass_p_one, place.mass_p, comparison, 'either sign'))
        '  this takes the mice''s mean AI, which does not. Every split''s cluster separates'
        '  its own groups by construction: hence the null median under each p (over the'
        '  splits with a cluster). Top right: both statistics over the same splits.'
        }];
end
left = [left; {
    sprintf(['Other rows, %s: whole regions, no selection; p: the values permuted ' ...
             'over the splits.'], strjoin(regions, ' and '))
    }];
if comparison.direction_named
    left{end + 1} = sprintf(['One-sided p (%s > %s) first: named before any number ' ...
        '(RWS potentiates the'], exp_type, ctrl_type);
    left{end + 1} = '  stimulated barrels, Gambino et al. 2014).';
else
    left{end + 1} = sprintf(['Two-sided p first: the one-sided direction (%s > %s) ' ...
        'is carried over from RWS,'], exp_type, ctrl_type);
    left{end + 1} = '  not named for this comparison before any number.';
end
loo_stats = stat_of('loo_ai_raw');
left = [left; {
    sprintf(['Leave-one-out (Per_Mouse_LOO_%s): each mouse read out of sample, in a ' ...
             'cluster'], file_tag)
    sprintf(['  found without it; stricter about circularity by design. Its p, raw ' ...
             'stack: %s'], p_pair(loo_stats.p_exp_higher, loo_stats.p_two_sided, ...
             comparison, 'two-sided'))
    sprintf(['  (README: how it compares). %d values here, their p not corrected ' ...
             'across them.'], n_values)
    }];

right = {
    'AI: mean |L - R| over mean (L + R) over the cluster''s or the region''s voxels'
    '  where the mouse has a value; L + R / isocortex: the region''s mean L + R over'
    '  the isocortex''s.'
    'Raw stack: the collected stack less the mouse''s off-tissue level (median of its'
    '  background voxels outside the atlas brain), smoothed as step 3: its zero is no'
    '  signal, so the raw AI is the ratio to read. Test maps: the maps of step 3, the'
    '  experimental group through the alignment line; their zero is the'
    '  normalisation''s, so the AI there is the index as the test sees it.'
    'No signed value: left and right are not certain for every brain, so the'
    '  stimulated side is unknown mouse by mouse (the test too takes |L - R|).'
    };

% the autofluorescence in the cluster, against the nano channel's raw AI: a
% misregistration of the surface would raise |L - R| in both, in the same mice
if any(strcmp(T_stats.value, 'sm_ai_auto'))
    auto_stats = stat_of('sm_ai_auto');
    nano_stats = stat_of('sm_ai_raw');
    [r_pooled, r_within] = channel_correlation(T_mice.sm_ai_raw, T_mice.sm_ai_auto, ...
        strcmp(T_mice.group, ctrl_type));
    right = [right; {
        'Autofluorescence (third panel): the channel registered with the nano stack, less'
        '  its off-tissue level, smoothed in the nano channel''s tissue, read in the same'
        '  clusters and tested the same way. A misregistration of the surface would'
        '  raise |L - R| in both channels, in the same mice. Here:'
        sprintf('  %s - %s %+.3f (nano, raw stack, %+.3f), p %s;', exp_type, ...
            ctrl_type, auto_stats.difference, nano_stats.difference, ...
            p_pair(auto_stats.p_exp_higher, auto_stats.p_two_sided, comparison, 'either'))
        sprintf(['  its AI against the nano raw AI over the mice r %.2f, within the ' ...
                 'groups r %.2f.'], r_pooled, r_within)
        }];
end

% the mice without a test-map AI in the cluster, and with little of a region
if any(strcmp(T_mice.Properties.VariableNames, 'sm_ai'))
    no_ai = isnan(T_mice.sm_ai);
    if any(no_ai)
        right{end + 1} = ['No AI on the test maps (mean L + R not above 0) in the ' ...
            'cluster: ' strjoin(short_name(T_mice.mouse(no_ai)), ', ')];
    end
end
for r = 1:numel(regions)
    coverage = T_mice.(['coverage_' region_tag(regions{r})]);
    is_low = coverage < 0.5;
    if any(is_low)
        right{end + 1} = sprintf('Less than half of %s: %s', regions{r}, ...
            strjoin(short_name(T_mice.mouse(is_low)), ', ')); %#ok<AGROW>
    end
end

% what the comparison itself leaves open
right = [right(:); comparison_caveat(exp_type, comparison.slope)];
end

function lines = loo_notes_lines(T_stats, T_loo, comparison, relabelled, file_tag)
% The notes of the leave-one-out's figure: what it is, how it differs from the
% selection-matched test, how its p is made, what the same procedure gives
% with the groups swapped, and the folds without a cluster or an AI; the
% settings and the numbers from the run, not typed.

ctrl_type = comparison.ctrl_type;
exp_type = comparison.exp_type;
perm_settings = comparison.perm_settings;
short_name = @(names) cellfun(@(n) strtok(n, '_'), names, 'UniformOutput', false);
loo_stats = T_stats(strcmp(T_stats.value, 'loo_ai_raw'), :);
loo_stats_test = T_stats(strcmp(T_stats.value, 'loo_ai'), :);
sm_stats = T_stats(strcmp(T_stats.value, 'sm_ai_raw'), :);

lines = {
    sprintf(['Leave-one-out: each mouse read in the heaviest cluster of %s where ' ...
             '|L - R| is higher in %s'], comparison.cluster_region, exp_type)
    sprintf(['  (a positive t at p < %g; median over +/- %d planes; %d-connected; a t ' ...
             'where each'], perm_settings.cluster_p, perm_settings.slab_range, ...
             perm_settings.cluster_connectivity)
    sprintf(['  group has %d mice) of the comparison redone without it, the ' ...
             'alignment refitted.'], perm_settings.min_mice_per_group)
    'Out of sample: no mouse is read in a cluster its own data helped choose, and each'
    '  cluster comes from one mouse fewer.'
    };
if ~isempty(sm_stats)
    lines = [lines; {
        '  Stricter about circularity by design than the selection-matched test, whose p,'
        sprintf('  raw stack, is %s', p_pair(sm_stats.p_exp_higher, ...
            sm_stats.p_two_sided, comparison, 'either direction'))
        sprintf('  (Per_Mouse_Values_%s);', file_tag)
        '  why the leave-one-out''s p can be the lower: README.'
        }];
end
if loo_stats.n_splits > 0
    lines = [lines; {
        sprintf(['p: every fold redone under each of the %d splits of the mice, its ' ...
                 'cluster found'], loo_stats.n_splits + loo_stats.n_splits_empty)
        sprintf(['  again with the split''s groups; the share of the %d splits with a ' ...
                 'value in both'], loo_stats.n_splits)
        '  groups whose difference of the group means reaches the observed one.'
        sprintf(['  Two-sided: whose |difference| reaches it, every fold''s cluster on ' ...
                 'the split''s %s-higher'], exp_type)
        '  side (not the selection-matched test''s either direction, which searches both).'
        }];

    % the same procedure with the groups swapped, when the splits include it
    swapped = [swapped_difference(relabelled, 'ai_raw'), ...
        swapped_difference(relabelled, 'ai')];
    if ~any(isnan(swapped))
        lines = [lines; {
            sprintf(['With the groups swapped (its mirror split) the same procedure ' ...
                     'reads %s-higher'], ctrl_type)
            sprintf(['  clusters: %s - %s %+.3f on the raw stack (%+.3f, test maps), ' ...
                     'against'], ctrl_type, exp_type, swapped)
            sprintf('  %s - %s %+.3f (%+.3f) with the groups as they are.', exp_type, ...
                ctrl_type, loo_stats.difference, loo_stats_test.difference)
            }];
    end
else
    lines{end + 1} = 'Not redone under the splits (loo_relabel = false): no test.';
end
if comparison.direction_named
    lines{end + 1} = sprintf('One-sided p (%s > %s) first, named before any number.', ...
        exp_type, ctrl_type);
else
    lines{end + 1} = sprintf(['Two-sided p first: the one-sided direction (%s > %s) ' ...
        'is carried over from RWS.'], exp_type, ctrl_type);
end

% the folds without a cluster, and the mice without an AI in theirs
no_cluster = T_loo.fold > 0 & T_loo.cluster_n == 0;
if any(no_cluster)
    lines{end + 1} = ['No cluster without: ' ...
        strjoin(short_name(T_loo.left_out(no_cluster)), ', ')];
else
    lines{end + 1} = 'Every fold has a cluster.';
end
no_ai = T_loo.fold > 0 & T_loo.cluster_n > 0 & isnan(T_loo.left_out_ai);
if any(no_ai)
    lines{end + 1} = ['No AI on the test maps (mean L + R not above 0) in its fold''s ' ...
        'cluster: ' strjoin(short_name(T_loo.left_out(no_ai)), ', ')];
end
lines{end + 1} = sprintf('Values and the other tests: Per_Mouse_Values_%s.', file_tag);
end

function lines = comparison_caveat(exp_type, slope)
% The caveat of a comparison, for the notes: after behavior, the alignment line
% is steep, and the broad L + R increase it shows may come from it.

switch exp_type
    case 'behavior'
        lines = {
            sprintf(['Caveat: the line aligning behavior onto naive has slope %.2f; ' ...
                     'the broad'], slope)
            '  L + R increase after behavior may come from the normalisation (open in'
            '  docs/ROADMAP.md). The raw values carry neither the line nor step 2.'
            };
    otherwise
        lines = {};
end
end
