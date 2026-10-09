function per_mouse_region_values(run_settings)
%PER_MOUSE_REGION_VALUES  Each mouse's values in the regions named in advance.
%   PER_MOUSE_REGION_VALUES(run_settings) does the work of run_per_mouse_values,
%   which sets the fields of run_settings and says what each one does.
%
%   Reads, from the control and the experimental group's folders,
%   <channel>_4d_normalized.mat and <channel>_4d_normalized_bkgmask.mat
%   (run_normalise_groups) and the collected stack <channel>_4d.mat
%   (run_collect_by_group), one mouse at a time. Writes into comp_out_dir, the
%   folder of run_group_differences, with <tag> = <comp_tag>_smooth<sigma>:
%     Per_Mouse_Values_<tag>.csv   one row per mouse, every value below
%     Per_Mouse_Stats_<tag>.csv    one row per value: group means and SEMs, the
%                                  exact permutation p, Welch t, Hedges g
%     Per_Mouse_LOO_<tag>.csv      one row per leave-one-out fold
%     Per_Mouse_LOO_<tag>.mat      the folds' clusters, voxel by voxel
%     Per_Mouse_Values_<tag>       the figure, .fig and .png
%     Per_Mouse_Maps_<tag>.mat     the cache of each mouse's maps (see below)
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
%   Leave-one-out: each mouse in turn is left out, and the comparison is done
%   again on the others as the test does it (the alignment refitted on their
%   profiles; Welch t, the fewest mice per group, surprise, the rolling median,
%   p < cluster_p, connectivity), here by region_permutation_test on the
%   observed split alone, within loo_region. The left-out mouse's ai (on that
%   fold's maps, loo_ai) and ai_raw (loo_ai_raw) are read in the fold's heaviest
%   cluster where |L - R| is higher in the experimental group (a positive t),
%   so no mouse is read in a cluster its own data helped define. Fold 0 leaves
%   no mouse out: its cluster must be that of run_group_differences, which is
%   checked against its table when it is there. A fold with no cluster gives
%   NaN. Not redone without the mouse: run_normalise_groups, which fits every
%   mouse of a group onto the group's median cortex, so without one mouse the
%   others' lines would change too. To first order they change by one line
%   common to the group, which leaves the fold's t of |L - R| as it is (the
%   fold's alignment is refitted, and a scale common to both groups does not
%   change a t); what is left is second order.
%
%   Statistics, per value, over the mice with a value: the exact permutation p
%   of the difference of the group means (experimental minus control), over
%   every split of the pooled mice into groups of the observed sizes (252 for 5
%   and 5), the observed one included, one-sided for the experimental group
%   higher and two-sided; the Welch t for reference; Hedges' g. With
%   direction_named, the one-sided p is the test named before any number and
%   comes first (RWS potentiates the stimulated barrels' synapses and brings
%   AMPA receptors to their surface, Gambino et al. 2014); without it, the
%   two-sided p comes first. For the leave-one-out values this p shuffles the
%   values, each read in its fold's cluster under the true groups. With
%   loo_relabel the leave-one-out is also redone under every split, each fold's
%   cluster found again with the split's groups (the fold keeping its
%   alignment, as the test keeps it under its splits), and the p of the
%   difference of the group means is taken over those: the full test, which
%   leads for those two values.
%
%   The maps are cached: each mouse's region sums and the box around loo_region
%   that the leave-one-out reads, so the statistics and the figure can be redone
%   without smoothing every mouse again; force_recompute_mice redoes them. The
%   cache is read only if each group's normalised stack still holds the mice
%   and the lines of run_normalise_groups it was made from. The redone
%   leave-one-out is cached beside it, Per_Mouse_LOO_Relabelled_<tag>.mat.

% settings of run_per_mouse_values, under the names the code below uses
paths = run_settings.paths;
ctrl_type = run_settings.ctrl_type;
exp_type = run_settings.exp_type;
behavior_mice = run_settings.behavior_mice;
regions = run_settings.regions;
loo_region = run_settings.loo_region;
apply_smoothing = run_settings.apply_smoothing;
smooth_sigma = run_settings.smooth_sigma;
min_mice_per_group = run_settings.min_mice_per_group;
slab_range = run_settings.slab_range;
cluster_p = run_settings.cluster_p;
cluster_connectivity = run_settings.cluster_connectivity;
force_recompute_mice = run_settings.force_recompute_mice;
direction_named = run_settings.direction_named;
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

% the leave-one-out's region must be one of the regions named
if ~ismember(loo_region, regions)
    error(['run_per_mouse_values: loo_region is %s, not one of the regions named in ' ...
           'regions (%s).'], loo_region, strjoin(regions, ', '));
end

% the file names carry the smoothing
if apply_smoothing
    smooth_suffix = sprintf('_smooth%g', smooth_sigma);
else
    smooth_suffix = '_nosmooth';
end
file_tag = [comp_tag smooth_suffix];

% the leave-one-out's test: the observed split alone, on no pool; the top volume
% and the quantile at run_group_differences' values, though only the cluster is read
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

%% Atlas and regions

% the 10 um annotation on the volumes' crop, the regions of the bars, and the
% voxels of each region named in advance, of the isocortex, and of the box the
% leave-one-out reads
A = get_atlas_crop('ccf');
brainMask = A.brainMask;
[T_regions, valid_pixels, region_of_voxel] = surprise_regions(A.annot, paths.atlas);
clear A
masks = region_masks(T_regions, valid_pixels, region_of_voxel, regions, loo_region, ...
    slab_range);
clear valid_pixels region_of_voxel

%% Each mouse of both groups

% every control mouse saved; for behavior, the mice behavior_mice names
if strcmp(exp_type, 'behavior')
    exp_named = behavior_mice;
else
    exp_named = {};
end

% the settings the maps depend on, checked when the cache is read
cache_settings = struct('regions', {regions}, 'loo_region', loo_region, ...
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

%% Leave-one-out clusters

% each mouse's values in the cluster its group comparison finds without it
fprintf('Leave-one-out clusters in %s...\n', loo_region);
step3_table = fullfile(comp_out_dir, ['Region_Surprise_DiffSum_' comp_tag '.csv']);
[T_loo, loo] = leave_one_out(ctrl_mice, exp_mice, masks, perm_settings, step3_table, ...
    loo_region);
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
save(fullfile(comp_out_dir, ['Per_Mouse_LOO_' file_tag '.mat']), 'loo_clusters', ...
    'loo_box', 'loo_region');

%% Leave-one-out under every split

% the leave-one-out redone under every split of the mice, for the full p of its
% two values; from its cache when it was made from these stacks and settings,
% and the maps were not made again in this run
relabelled = [];
if loo_relabel
    relabel_file = fullfile(comp_out_dir, ['Per_Mouse_LOO_Relabelled_' file_tag '.mat']);
    relabel_settings = struct('cache_settings', cache_settings, ...
        'ctrl_norm_params', ctrl_mice.stack_norm_params, ...
        'exp_norm_params', exp_mice.stack_norm_params, ...
        'min_mice_per_group', min_mice_per_group, 'slab_range', slab_range, ...
        'cluster_p', cluster_p, 'cluster_connectivity', cluster_connectivity);
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

%% Statistics and figure

% the values tested and drawn, in the figure's order
values = value_list(regions, loo_region);

% each value's group test, with the full p of the leave-one-out's values, and
% the tables
T_stats = value_statistics(T_mice, values, ctrl_type, exp_type);
T_stats = add_relabelled_p(T_stats, relabelled);
writetable(T_mice, fullfile(comp_out_dir, ['Per_Mouse_Values_' file_tag '.csv']));
writetable(T_stats, fullfile(comp_out_dir, ['Per_Mouse_Stats_' file_tag '.csv']));
writetable(T_loo, fullfile(comp_out_dir, ['Per_Mouse_LOO_' file_tag '.csv']));
print_summary(T_mice, T_stats, T_loo, values);

% one panel per value, one dot per mouse; the figure's text from the settings
comparison = struct('ctrl_type', ctrl_type, 'exp_type', exp_type, ...
    'direction_named', direction_named, 'slope', slope, 'intercept', intercept, ...
    'perm_settings', perm_settings);
plot_per_mouse_values(T_mice, T_stats, T_loo, values, regions, comparison, file_tag, ...
    comp_out_dir);
fprintf('Per-mouse values saved to: %s\n', comp_out_dir);

end

% ===== Local functions: regions and mice =====

function masks = region_masks(T_regions, valid_pixels, region_of_voxel, regions, ...
    loo_region, slab_range)
% The voxels of each region named in advance and of the isocortex, as linear
% indices into the folded grid, and the box around loo_region the leave-one-out
% reads: its planes and the slab_range planes on each side, which the rolling
% median of a voxel in it reads.

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

% the band of the leave-one-out: the region and slab_range planes on each side
% along AP, within the atlas
in_region = region_vol == find(strcmp(T_regions.acronym, loo_region));
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
fprintf(['  leave-one-out box around %s: %d x %d x %d voxels, %d in the band, %d in ' ...
         'the region.\n'], loo_region, masks.box_size, numel(masks.band_lin), ...
        nnz(masks.band_in_region));
end

function G = group_mice(group, group_dir, channel, named_mice, brainMask, masks, ...
    apply_smoothing, smooth_sigma)
% One group's mice, one at a time from the stacks: each mouse's plane profile,
% its sums in the regions on the normalised and on the raw stack, its
% off-tissue level, and both smoothed stacks in the leave-one-out box. The mice
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

% one mouse at a time
n_mice = numel(names);
G = struct();
G.group = group;
G.names = names;

% the normalised stack's mice and their lines, which the cache is checked against
G.stack_mice = S_norm.current_mice;
G.stack_norm_params = S_norm.norm_params;
G.profiles = nan(raw_size(1), n_mice);
G.test = cell(n_mice, 1);
G.raw = cell(n_mice, 1);
G.background = nan(n_mice, 1);
G.background_planes = nan(raw_size(1), n_mice);
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

% ===== Local functions: values and leave-one-out =====

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

function [T_loo, loo] = leave_one_out(ctrl_mice, exp_mice, masks, perm_settings, ...
    step3_table, loo_region)
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
        check_against_step3(cluster, step3_table, loo_region);
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
        fprintf('  no cluster with a positive t in %s without %s.\n', loo_region, ...
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

function check_against_step3(cluster, step3_table, loo_region)
% Fold 0 must find the region's cluster that run_group_differences found, with
% the same voxels and mass, when its table is there; a warning where it differs.

if ~exist(step3_table, 'file')
    fprintf('  no table of run_group_differences to check fold 0 against (%s).\n', ...
        step3_table);
    return
end
T_step3 = readtable(step3_table);
row = strcmp(T_step3.acronym, loo_region);
if T_step3.lr_diff_cluster_sign(row) <= 0
    fprintf(['  run_group_differences'' heaviest %s cluster is negative; fold 0''s ' ...
             'positive one is not in its table.\n'], loo_region);
    return
end
step3_n = T_step3.lr_diff_cluster_n(row);
step3_mass = T_step3.lr_diff_cluster_score(row);
if cluster.n == step3_n && abs(cluster.mass - step3_mass) <= 1e-6 * step3_mass
    fprintf(['  fold 0 gives the cluster of run_group_differences: %d voxels, mass ' ...
             '%.2f.\n'], cluster.n, cluster.mass);
else
    warning(['run_per_mouse_values: fold 0 gives a %s cluster of %d voxels, mass ' ...
             '%.2f; run_group_differences found %d voxels, mass %.2f (%s).'], ...
             loo_region, cluster.n, cluster.mass, step3_n, step3_mass, step3_table);
end
end

% ===== Local functions: statistics =====

function values = value_list(regions, loo_region)
% The values tested and drawn, in the figure's order: per region, the ai and
% the relative L + R on the test's maps and the raw stack; then the
% leave-one-out's. Each with its column, its title and its axis label.

values = struct('name', {}, 'region', {}, 'title', {}, 'label', {});
ai_label = 'mean |L - R| / mean (L + R)';
sum_label = 'L + R relative to isocortex';
kinds = {
    'ai_',          'AI, test maps',                ai_label
    'ai_raw_',      'AI, raw stack',                ai_label
    'sum_rel_',     'L + R / isocortex, test maps', sum_label
    'sum_rel_raw_', 'L + R / isocortex, raw stack', sum_label
    };
for r = 1:numel(regions)
    for k = 1:size(kinds, 1)
        values(end + 1) = struct('name', [kinds{k, 1} region_tag(regions{r})], ...
            'region', regions{r}, 'title', kinds{k, 2}, 'label', kinds{k, 3}); %#ok<AGROW>
    end
end
values(end + 1) = struct('name', 'loo_ai', 'region', loo_region, 'title', ...
    'AI in the leave-one-out cluster, test maps', 'label', ai_label);
values(end + 1) = struct('name', 'loo_ai_raw', 'region', loo_region, 'title', ...
    'AI in the leave-one-out cluster, raw stack', 'label', ai_label);
end

function T_stats = value_statistics(T_mice, values, ctrl_type, exp_type)
% Per value, over the mice with a value: each group's mean and SEM, the exact
% permutation p (one-sided for the experimental group higher, the direction
% named before the experiment, and two-sided), the Welch t and Hedges' g.

is_ctrl = strcmp(T_mice.group, ctrl_type);
rows = cell(numel(values), 1);
for v = 1:numel(values)
    x = T_mice.(values(v).name);
    x_ctrl = x(is_ctrl & ~isnan(x));
    x_exp = x(~is_ctrl & ~isnan(x));

    % the permutation p, and the Welch t for reference
    [p_two, p_one, n_splits] = exact_permutation(x_ctrl, x_exp);
    [~, p_welch, ~, welch] = ttest2(x_exp, x_ctrl, 'Vartype', 'unequal');

    rows{v} = struct('value', values(v).name, 'region', values(v).region, ...
        'title', values(v).title, 'ctrl', ctrl_type, 'exp', exp_type, ...
        'n_ctrl', numel(x_ctrl), 'n_exp', numel(x_exp), 'mean_ctrl', mean(x_ctrl), ...
        'sem_ctrl', std(x_ctrl) / sqrt(numel(x_ctrl)), 'mean_exp', mean(x_exp), ...
        'sem_exp', std(x_exp) / sqrt(numel(x_exp)), ...
        'difference', mean(x_exp) - mean(x_ctrl), 'p_perm_exp_higher', p_one, ...
        'p_perm_two_sided', p_two, 'n_splits', n_splits, 'welch_t', welch.tstat, ...
        'welch_df', welch.df, 'welch_p', p_welch, 'hedges_g', hedges_g(x_ctrl, x_exp));
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

function T_stats = add_relabelled_p(T_stats, relabelled)
% The full p of the leave-one-out's values, from the leave-one-out redone under
% every split, as columns of the statistics: NaN for the other values, and for
% every value when it was not redone.

T_stats.p_relabelled_exp_higher = nan(height(T_stats), 1);
T_stats.p_relabelled_two_sided = nan(height(T_stats), 1);
T_stats.n_splits_relabelled = zeros(height(T_stats), 1);
if isempty(relabelled)
    return
end
loo_values = {'loo_ai', 'ai'; 'loo_ai_raw', 'ai_raw'};
for v = 1:size(loo_values, 1)
    row = strcmp(T_stats.value, loo_values{v, 1});
    [p_one, p_two, n_splits] = relabelled_p(relabelled.(loo_values{v, 2}), ...
        relabelled.in_ctrl);
    T_stats.p_relabelled_exp_higher(row) = p_one;
    T_stats.p_relabelled_two_sided(row) = p_two;
    T_stats.n_splits_relabelled(row) = n_splits;
end
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

function print_summary(T_mice, T_stats, T_loo, values)
% The values of every mouse, the leave-one-out folds and each value's test.

fprintf('\nPer-mouse values:\n');
columns = [{'mouse', 'group'}, {values.name}];
disp(T_mice(:, columns));
fprintf('Leave-one-out folds:\n');
disp(T_loo);
fprintf('Group tests:\n');
disp(T_stats(:, {'value', 'n_ctrl', 'n_exp', 'mean_ctrl', 'mean_exp', ...
    'p_perm_exp_higher', 'p_perm_two_sided', 'n_splits', 'welch_p', 'hedges_g', ...
    'p_relabelled_exp_higher', 'p_relabelled_two_sided', 'n_splits_relabelled'}));
end

% ===== Local functions: figure =====

function plot_per_mouse_values(T_mice, T_stats, T_loo, values, regions, comparison, ...
    file_tag, comp_out_dir)
% One panel per value, one dot per mouse with its name, each group's mean and
% SEM, and the p of the value: a row per region named in advance, then the
% leave-one-out's values and the notes.

ctrl_type = comparison.ctrl_type;
exp_type = comparison.exp_type;

% a row per region, a column per kind of value, the leave-one-out in the last
% row: four columns, not the three of docs/STYLE.md, so a row is one region
n_kinds = 4;
n_rows_grid = numel(regions) + 1;
n_cols_grid = n_kinds;
n_splits = max(T_stats.n_splits);

fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'Normalized', ...
    'Position', [0 0 1 1]);
is_ctrl = strcmp(T_mice.group, ctrl_type);
short_names = cellfun(@(n) strtok(n, '_'), T_mice.mouse, 'UniformOutput', false);
for v = 1:numel(values)

    % the panel: the region's row, or the start of the last row for the
    % leave-one-out
    if startsWith(values(v).name, 'loo_')
        panel = (n_rows_grid - 1) * n_cols_grid + v - numel(regions) * n_kinds;
    else
        panel = v;
    end
    subplot(n_rows_grid, n_cols_grid, panel);
    panel_title = sprintf('%s: %s', values(v).region, values(v).title);

    % the dots, names, means and SEMs, and the p
    x = T_mice.(values(v).name);
    draw_value_panel(x(is_ctrl), x(~is_ctrl), short_names(is_ctrl), ...
        short_names(~is_ctrl), ctrl_type, exp_type);
    stat = T_stats(strcmp(T_stats.value, values(v).name), :);
    title(panel_title, 'FontSize', 10, 'Interpreter', 'none');
    subtitle(stat_lines(stat, comparison), 'FontSize', 8, 'Interpreter', 'none');
    ylabel(values(v).label, 'FontSize', 9);
end

% the notes, beside the leave-one-out's panels
subplot(n_rows_grid, n_cols_grid, ...
    [(n_rows_grid - 1) * n_cols_grid + 3, n_rows_grid * n_cols_grid]);
axis off;
text(0, 1.05, notes_lines(T_mice, T_stats, T_loo, regions, comparison, ...
    numel(values)), 'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 8, ...
    'Interpreter', 'none');

% the title: the comparison, its alignment and its test
title_line = sprintf('Per-mouse values in the regions named in advance - %s', ...
    strrep(file_tag, '_', ' '));
align_line = sprintf(['%s aligned onto %s by the line %.3f x %+.1f; p: exact ' ...
                      'permutation of the difference of the group means, over every ' ...
                      'split of the mice with a value (%d for all of them)'], exp_type, ...
                      ctrl_type, comparison.slope, comparison.intercept, n_splits);
sgtitle({title_line, ['\rm\fontsize{11}' align_line]}, 'FontSize', 14, ...
    'FontWeight', 'bold');

% save it
saveas(fig, fullfile(comp_out_dir, ['Per_Mouse_Values_' file_tag '.fig']));
exportgraphics(fig, fullfile(comp_out_dir, ['Per_Mouse_Values_' file_tag '.png']), ...
    'Resolution', 300);
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

function lines = stat_lines(stat, comparison)
% The p of a value under its panel's title: its exact permutation p, one-sided
% (experimental group higher) and two-sided, the one named before any number
% first; for the leave-one-out's values the p of the leave-one-out redone under
% every split leads, and the p of the values shuffled follows.

shuffled = p_pair(stat.p_perm_exp_higher, stat.p_perm_two_sided, comparison);
rest = sprintf('Welch p %.3f, Hedges g %.2f, n %d and %d', stat.welch_p, ...
    stat.hedges_g, stat.n_ctrl, stat.n_exp);
if stat.n_splits_relabelled > 0
    redone = p_pair(stat.p_relabelled_exp_higher, stat.p_relabelled_two_sided, ...
        comparison);
    lines = {sprintf('p %s, folds redone per split', redone), ...
        sprintf('values shuffled: p %s', shuffled), rest};
else
    lines = {sprintf('perm p %s', shuffled), rest};
end
end

function text_p = p_pair(p_one, p_two, comparison)
% A one-sided p (experimental group higher) and a two-sided one: the one-sided
% first when its direction was named before any number, the two-sided first
% otherwise.

one_sided = sprintf('%.3f %s > %s', p_one, comparison.exp_type, comparison.ctrl_type);
two_sided = sprintf('%.3f two-sided', p_two);
if comparison.direction_named
    text_p = [one_sided ', ' two_sided];
else
    text_p = [two_sided ', ' one_sided];
end
end

function lines = notes_lines(T_mice, T_stats, T_loo, regions, comparison, n_values)
% The notes beside the leave-one-out's panels: what each value is, why none is
% signed, how the leave-one-out and the p are made, the folds without a
% cluster, the mice with little of a region, and the comparison's caveat; the
% settings from the run, not typed.

ctrl_type = comparison.ctrl_type;
exp_type = comparison.exp_type;
perm_settings = comparison.perm_settings;
short_name = @(names) cellfun(@(n) strtok(n, '_'), names, 'UniformOutput', false);

lines = {
    'AI: mean |L - R| over mean (L + R) in the region, each mouse over its own voxels;'
    '  L + R / isocortex: the region''s mean L + R over the isocortex''s.'
    'Raw stack: the collected stack less the mouse''s off-tissue level (median of its'
    '  background voxels outside the atlas brain), smoothed as step 3: its zero is no'
    '  signal, so the raw AI is the ratio to read. Test maps: the maps of step 3, the'
    '  experimental group through the alignment line; their zero is the'
    '  normalisation''s, so the AI there is the index as the test sees it.'
    'No signed value: left and right are not certain for every brain, so the'
    '  stimulated side is unknown mouse by mouse (the test too takes |L - R|).'
    sprintf(['Leave-one-out: each mouse read in the heaviest cluster where |L - R| ' ...
             'is higher in %s'], exp_type)
    sprintf(['  (a positive t at p < %g; median over +/- %d planes; %d-connected; a t ' ...
             'where each'], perm_settings.cluster_p, perm_settings.slab_range, ...
             perm_settings.cluster_connectivity)
    sprintf('  group has %d mice) of the comparison redone without it.', ...
        perm_settings.min_mice_per_group)
    };

% which p comes first, and why
if comparison.direction_named
    lines{end + 1} = sprintf(['One-sided p (%s > %s) first: the direction named ' ...
        'before any number (RWS'], exp_type, ctrl_type);
    lines{end + 1} = '  potentiates the stimulated barrels, Gambino et al. 2014).';
else
    lines{end + 1} = sprintf(['Two-sided p first: the one-sided direction (%s > %s) ' ...
        'is carried over'], exp_type, ctrl_type);
    lines{end + 1} = '  from RWS, not named for this comparison before any number.';
end
% the leave-one-out redone under every split: how many splits, and how many of
% them leave each group a mouse with a value
n_with_value = T_stats.n_splits_relabelled(T_stats.n_splits_relabelled > 0);
if ~isempty(n_with_value)
    n_splits_all = nchoosek(height(T_mice), nnz(strcmp(T_mice.group, ctrl_type)));
    if min(n_with_value) == max(n_with_value)
        with_value = sprintf('%d', min(n_with_value));
    else
        with_value = sprintf('%d to %d', min(n_with_value), max(n_with_value));
    end
    lines{end + 1} = sprintf(['Leave-one-out p: every fold redone under each of the %d ' ...
        'splits of the mice'], n_splits_all);
    lines{end + 1} = sprintf(['  (%s with a value in both groups); "values shuffled" ' ...
        'keeps each fold''s'], with_value);
    lines{end + 1} = '  cluster as the true groups gave it.';
end
lines{end + 1} = sprintf(['%d values per comparison; their p are not corrected ' ...
    'across them.'], n_values);

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

% the mice with less than half of a region
for r = 1:numel(regions)
    coverage = T_mice.(['coverage_' region_tag(regions{r})]);
    is_low = coverage < 0.5;
    if any(is_low)
        lines{end + 1} = sprintf('Less than half of %s: %s', regions{r}, ...
            strjoin(short_name(T_mice.mouse(is_low)), ', ')); %#ok<AGROW>
    end
end

% what the comparison itself leaves open
lines = [lines(:); comparison_caveat(exp_type, comparison.slope)];
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
