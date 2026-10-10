function mouse_influence(run_settings)
%MOUSE_INFLUENCE  Each mouse's part in the barrel-field cluster of the plasticity comparison.
%   MOUSE_INFLUENCE(run_settings) does the work of run_mouse_influence, which
%   sets the fields of run_settings and says what each one does.
%
%   Reads, from comp_out_dir, what run_per_mouse_values wrote for the same
%   comparison and smoothing, <tag> = <comp_tag>_smooth<sigma>: each mouse's
%   maps in the box around cluster_region (Per_Mouse_Maps_<tag>.mat, checked
%   against the groups' normalised stacks) and its autofluorescence there
%   (Per_Mouse_Auto_<tag>.mat); each mouse's collected stack (nano_4d.mat of
%   its group), unsmoothed, in the window of the maps; and, to check against,
%   the leave-one-out's clusters (Per_Mouse_LOO_<tag>.mat), the mice's values
%   (Per_Mouse_Values_<tag>.csv) and the table of run_group_differences.
%   Writes into comp_out_dir:
%     Influence_Folds_<tag>        .csv, one row per fold (all the mice, then
%                                  each mouse left out), and the figure of the
%                                  folds and the slopes
%     Influence_Slopes_<tag>.csv   one row per slope of the alignment, the
%                                  comparison of all the mice redone at it
%     Influence_Maps_<tag>_<reading>  the mice's maps around the cluster, one
%                                  figure per reading (raw, auto, test)
%     Influence_Profiles_<tag>     each mouse's asymmetry along AP through the
%                                  cluster, with its sections: .csv, one row
%                                  per plane, and the figure
%     Influence_Consistency_<tag>  .csv, one row per reading and level, and the
%                                  figure
%     Influence_Mice_<tag>.csv     one row per mouse: its values in the cluster,
%                                  where it is above every mouse of the other
%                                  group, where its peak asymmetry sits
%     Influence_Cluster_<tag>.csv  the cluster's size and place
%     Influence_<tag>.mat          the cache of the folds and of the counts
%                                  under every split
%
%   Readings. raw: |L - R| / (L + R) on the collected stack less the mouse's
%   off-tissue level, smoothed as step 3 and folded (the raw stack of
%   run_per_mouse_values): no normalisation, no alignment; the reading to look
%   at. auto: the same on the autofluorescence. test: |L - R| on the test's
%   maps (the experimental mice through the alignment of all the mice), what
%   step 3's t compares; their L + R carries the normalisation's zero, so no
%   ratio is taken there. unscaled, for the counts only: the test's maps at
%   slope 1, each group on its own scale of step 2 (the control mice as in
%   test, each experimental mouse's |L - R| divided by the slope). In a view or
%   over the cluster, the ratio is of the means, as run_per_mouse_values reads
%   its asymmetry index. The views: a
%   coronal one, the mean over the cluster's planes; one from above, the mean
%   over the cluster's depth below the pia in each column of the window; a
%   profile along AP, the mean over the cluster's coronal footprint in each
%   plane. Each mouse's sections show in its collected stack, unsmoothed: the
%   registered volume holds a section over the planes nearest it, so the
%   change from one plane to the next is about 0 within a section and peaks
%   between two. Each mouse's coverage along AP: the share of the footprint's
%   voxels where its raw stack has a value.
%
%   Folds. Fold 0 keeps every mouse; fold k leaves mouse k out, control mice
%   first. Each is step 3's comparison of its mice: the alignment refitted on
%   their profiles, then region_permutation_test on the band around
%   cluster_region under every split of the fold's mice into groups of its
%   sizes (252 for 5 and 5, 126 for 5 and 4). Its heaviest cluster where
%   |L - R| is higher in the experimental group (a positive t), its voxels,
%   mass and share of the cluster of all the mice; two p: positive, the share of
%   the splits whose positive cluster is as heavy; either sign, step 3's (each
%   split's larger sign against the observed larger, whose sign is given).
%   Fold 0 must be step 3's cluster and p, and each fold's cluster the
%   leave-one-out's of run_per_mouse_values. As there, run_normalise_groups is
%   not redone without the mouse. The folds are then run again with the
%   alignment held at fold 0's: on the test's maps an experimental mouse's
%   |L - R| is its own times the alignment's slope (the intercept cancels in
%   L - R, the common factor in the t), so a fold's slope moves every
%   experimental mouse at once, and holding it shows what the mouse's own
%   values do.
%
%   Slopes. The comparison of all the mice redone with the alignment's slope at
%   each of sweep_slopes, its intercept and common factor at fold 0's: both
%   signs' heaviest clusters, their p, and step 3's p of either sign. The
%   splits of step 3 keep the fitted slope, so its p does not carry the
%   slope's own uncertainty.
%
%   Counts. At every voxel of cluster_region where every mouse has a value
%   (complete voxels), c, the number of experimental mice above every control
%   mouse, from 0 to the experimental group's size. On exchangeable mice a
%   voxel reaches c >= k with probability C(n_exp, k) / C(n_mice, k), 1/252
%   for five of five, which is also the mean share over every split. For each
%   level k: the voxels with c >= k against that expectation, which a group
%   that is more asymmetric everywhere exceeds; and the largest patch of voxels
%   with c >= k, connected within the region (cluster_connectivity), which asks
%   whether they lie together. Under every split of the mice the same, with the
%   split's groups; a level's p is the share of the splits that reach the
%   observed one, the observed split included. Every level is given, none
%   chosen.
%
%   Size. The cluster of all the mice in micrometres: the extent of its voxels
%   along AP, DV and ML, from the first voxel's edge to the last's; its depth
%   below the pia, the distance from a voxel's centre to the atlas brain's
%   outer edge (its unlabelled holes filled);
%   its widths along the surface and its thickness through it, the surface's
%   normal being the depth's mean gradient over the cluster; its layers in the
%   atlas; against one barrel column, about 300 um across with layer 2/3 from
%   128 to 418 um below the pia (the C2 column of the mouse, Lefort et al.
%   2009). Its planes also as Allen's CCF index, counted from 0 (plane p of
%   the 10 um annotation counted from 1 is index p - 1). Each mouse's peak
%   asymmetry over the complete voxels (raw), and for the experimental mice the
%   peak of their excess over the highest control mouse, with their distance
%   from the cluster: one voxel's maximum of a ratio, which lands where L + R is
%   low, at the region's edges, so it says only where the cluster is not. Each
%   mouse's tissue profile over the planes the alignment is fitted on, its mean
%   and its standard deviation along AP, which set the slope.

% settings of run_mouse_influence, under the names the code below uses
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
sweep_slopes = run_settings.sweep_slopes;
window_um = run_settings.window_um;
n_workers = run_settings.n_workers;
force_recompute = run_settings.force_recompute;
channel = run_settings.channel;
comp_tag = run_settings.comp_tag;
ctrl_dir = run_settings.ctrl_dir;
exp_dir = run_settings.exp_dir;
comp_out_dir = run_settings.comp_out_dir;

% the experimental groups of run_group_differences
if ~ismember(exp_type, {'rws', 'behavior'})
    error('run_mouse_influence: unknown exp_type ''%s'' (use ''rws'' or ''behavior'').', ...
        exp_type);
end

% the autofluorescence is read beside the nano channel
if ~strcmp(channel, 'nano')
    error(['run_mouse_influence: channel is ''%s''; the autofluorescence is read beside ' ...
           'the nano channel, so channel must be ''nano''.'], channel);
end

% the file names carry the smoothing, as run_per_mouse_values'
if apply_smoothing
    smooth_suffix = sprintf('_smooth%g', smooth_sigma);
else
    smooth_suffix = '_nosmooth';
end
file_tag = [comp_tag smooth_suffix];

% the clusters' test under every split of a fold's mice, keeping the clusters'
% voxels; the top volume and the quantile at run_group_differences' values,
% though only the cluster is read
perm_settings = struct();
perm_settings.min_mice_per_group = min_mice_per_group;
perm_settings.slab_range = slab_range;
perm_settings.p_thresh = cluster_p;
perm_settings.cluster_p = cluster_p;
perm_settings.cluster_connectivity = cluster_connectivity;
perm_settings.topvol_mm3 = 0.1;
perm_settings.region_quantile = 0.99;
perm_settings.n_permutations = 'all';
perm_settings.n_workers = n_workers;
perm_settings.seed = 0;
perm_settings.keep_cluster_voxels = true;
perm_settings.quiet = true;

% what is read besides the caches, to check against
step3_table = fullfile(comp_out_dir, ['Region_Surprise_DiffSum_' comp_tag '.csv']);
loo_file = fullfile(comp_out_dir, ['Per_Mouse_LOO_' file_tag '.mat']);
values_file = fullfile(comp_out_dir, ['Per_Mouse_Values_' file_tag '.csv']);

%% Atlas, the clusters' box and its anatomy

% the 10 um annotation on the volumes' crop, the regions of the bars, the box of
% run_per_mouse_values around the clusters' region, and each voxel's depth below
% the pia and layer in it
A = get_atlas_crop('ccf');
[T_regions, valid_pixels, region_of_voxel] = surprise_regions(A.annot, paths.atlas);
masks = region_box_masks(T_regions, valid_pixels, region_of_voxel, regions, ...
    cluster_region, slab_range);
clear valid_pixels region_of_voxel
anatomy = box_anatomy(A, masks, paths.atlas);
clear A

%% The mice's maps

% every control mouse saved; for behavior, the mice behavior_mice names, as in
% run_per_mouse_values
if strcmp(exp_type, 'behavior')
    exp_named = behavior_mice;
else
    exp_named = {};
end

% the maps of run_per_mouse_values, made with these settings from the groups'
% normalised stacks as they are now, and the autofluorescence made from them
cache_settings = struct('regions', {regions}, 'loo_region', cluster_region, ...
    'slab_range', slab_range, 'apply_smoothing', apply_smoothing, ...
    'smooth_sigma', smooth_sigma, 'exp_named', {exp_named});
maps_file = fullfile(comp_out_dir, ['Per_Mouse_Maps_' file_tag '.mat']);
[ctrl_mice, exp_mice] = load_maps(maps_file, cache_settings, ctrl_dir, exp_dir, channel);
split_settings = struct('cache_settings', cache_settings, ...
    'ctrl_norm_params', ctrl_mice.stack_norm_params, ...
    'exp_norm_params', exp_mice.stack_norm_params, ...
    'min_mice_per_group', min_mice_per_group, 'slab_range', slab_range, ...
    'cluster_p', cluster_p, 'cluster_connectivity', cluster_connectivity);
auto_file = fullfile(comp_out_dir, ['Per_Mouse_Auto_' file_tag '.mat']);
[ctrl_auto, exp_auto, auto_settings] = load_auto(auto_file, ctrl_mice, exp_mice, ...
    split_settings);

% the mice, control mice first, as every table and figure lists them
names = [ctrl_mice.names(:); exp_mice.names(:)];
n_ctrl = numel(ctrl_mice.names);
n_mice = numel(names);

% each mouse's tissue profile over the planes the alignment is fitted on: its
% mean and its spread along AP, which set the slope
profiles_all = [ctrl_mice.profiles, exp_mice.profiles];
fit_planes = align_exp_to_ctrl(profiles_all(:, 1:n_ctrl), profiles_all(:, n_ctrl + 1:end));
tissue = struct();
tissue.mean = mean(profiles_all(fit_planes, :), 1, 'omitnan')';
tissue.sd = std(profiles_all(fit_planes, :), 0, 1, 'omitnan')';
clear profiles_all

%% Each mouse left out

% the folds, the slopes and the counts under every split, from their cache when
% it was made from these maps, autofluorescence and settings; a part missing
% from it, or made at other slopes, is computed and added
influence_file = fullfile(comp_out_dir, ['Influence_' file_tag '.mat']);
influence_settings = struct('split_settings', split_settings, 'auto_settings', ...
    auto_settings);
cached = load_influence(influence_file, influence_settings, force_recompute);

% the folds with the alignment refitted, as step 3 would fit it on the mice kept
if isempty(cached.folds)
    cached.folds = fold_influence(ctrl_mice, exp_mice, masks, perm_settings, []);
    save_influence(influence_file, cached, influence_settings);
end
folds = cached.folds;

% fold 0 must be step 3's cluster and p, each fold's cluster the leave-one-out's
fprintf('Folds in %s:\n', cluster_region);
check_fold_zero(folds, step3_table, cluster_region);
check_against_loo(folds, loo_file);
cluster_lin = folds.cluster_voxels{1};
if isempty(cluster_lin)
    error(['run_mouse_influence: with all the mice there is no cluster in %s where ' ...
           '|L - R| is higher in %s; nothing to describe.'], cluster_region, exp_type);
end

% the folds again with the alignment of all the mice held; fold 0 must be the
% same as with it refitted
alignment = struct('slope', folds.slope(1), 'intercept', folds.intercept(1), ...
    'common_factor', folds.common_factor(1));
if isempty(cached.folds_held)
    fprintf('Each mouse left out, the alignment held at that of all the mice:\n');
    cached.folds_held = fold_influence(ctrl_mice, exp_mice, masks, perm_settings, ...
        alignment);
    save_influence(influence_file, cached, influence_settings);
end
folds_held = cached.folds_held;
check_held_fold_zero(folds, folds_held);

%% The alignment's slope

% the comparison of all the mice at each slope of sweep_slopes
if isempty(cached.sweep) || ~isequal(cached.sweep.slope(:), sweep_slopes(:))
    fprintf('All the mice at %d slopes of the alignment:\n', numel(sweep_slopes));
    cached.sweep = slope_sweep(ctrl_mice, exp_mice, masks, perm_settings, ...
        alignment, sweep_slopes);
    save_influence(influence_file, cached, influence_settings);
end
sweep = cached.sweep;

%% Each mouse's readings

% the window of the maps around the cluster of all the mice, and every mouse's
% readings: its views in the window, its value at every voxel of the region,
% its mean over the cluster
win = map_window(cluster_lin, masks, anatomy, window_um);
[views, values, cluster_values, profiles] = mouse_readings(ctrl_mice, exp_mice, ...
    ctrl_auto, exp_auto, alignment, anatomy.region_lin, cluster_lin, win);
clear ctrl_mice exp_mice ctrl_auto exp_auto
check_against_values(cluster_values, names, values_file);

% the sections along AP: each mouse's collected stack, unsmoothed, changes from
% one plane to the next only between two sections
plane_changes = [section_changes(ctrl_type, ctrl_dir, names(1:n_ctrl), win, masks), ...
    section_changes(exp_type, exp_dir, names(n_ctrl + 1:end), win, masks)];

%% Counts under every split

% for each reading, the experimental mice above every control mouse at each
% complete voxel, the voxels and the largest patch at each level under every
% split
readings = {'raw', 'test', 'unscaled', 'auto'};
if ~has_counts(cached.counts, readings)
    cached.counts = struct();
    for r = 1:numel(readings)
        fprintf('Counts of the %s reading under every split:\n', readings{r});
        cached.counts.(readings{r}) = count_splits(values.(readings{r}), n_ctrl, ...
            anatomy.region_lin, masks.box_size, cluster_connectivity, cluster_lin);
    end
    save_influence(influence_file, cached, influence_settings);
end
counts = cached.counts;
clear cached

%% Size and peaks

% the cluster in micrometres, against a barrel column; each mouse in it, where
% it is above every mouse of the other group, and its peak
cluster = cluster_size(cluster_lin, masks, anatomy);
T_mice = mouse_table(names, n_ctrl, ctrl_type, exp_type, values, cluster_values, ...
    counts, readings, tissue, anatomy, cluster_lin, win);

%% Tables and figures

% the tables
T_folds = fold_table(folds, folds_held, names, n_ctrl, ctrl_type, exp_type, win);
T_slopes = slope_table(sweep, folds, win);
T_counts = count_table(counts, readings);
T_cluster = struct2table(cluster, 'AsArray', true);
T_profiles = profile_table(profiles, plane_changes, win, names, ...
    anatomy.ccf_first_index);
writetable(T_folds, fullfile(comp_out_dir, ['Influence_Folds_' file_tag '.csv']));
writetable(T_slopes, fullfile(comp_out_dir, ['Influence_Slopes_' file_tag '.csv']));
writetable(T_profiles, fullfile(comp_out_dir, ['Influence_Profiles_' file_tag '.csv']));
writetable(T_counts, fullfile(comp_out_dir, ['Influence_Consistency_' file_tag '.csv']));
writetable(T_mice, fullfile(comp_out_dir, ['Influence_Mice_' file_tag '.csv']));
writetable(T_cluster, fullfile(comp_out_dir, ['Influence_Cluster_' file_tag '.csv']));
print_summary(T_folds, T_slopes, T_counts, T_mice, cluster);

% the figures, their text from the settings and the run
comparison = struct('ctrl_type', ctrl_type, 'exp_type', exp_type, 'cluster_region', ...
    cluster_region, 'perm_settings', perm_settings, 'n_ctrl', n_ctrl, 'n_mice', ...
    n_mice, 'names', {names}, 'ccf_first_index', anatomy.ccf_first_index, ...
    'box_first_plane', masks.box_ap(1), 'window_um', window_um);
plot_folds(T_folds, T_slopes, comparison, file_tag, comp_out_dir);
plot_maps(views, cluster_values, win, comparison, file_tag, comp_out_dir);
plot_consistency(counts, readings, T_mice, win, masks, comparison, file_tag, ...
    comp_out_dir);
plot_profiles(profiles, plane_changes, win, comparison, file_tag, comp_out_dir);
fprintf('Influence of each mouse saved to: %s\n', comp_out_dir);

end

% ===== Local functions: anatomy and loading =====

function anatomy = box_anatomy(A, masks, atlas_dir)
% The box's anatomy on its folded grid: each voxel's depth below the pia in um
% (from its centre to the atlas brain's outer edge; NaN outside the brain), its
% atlas value, the region's voxels, the layer of each atlas value, and where the
% box sits in the CCF.

% the brain without its holes: some voxels inside it carry no label, and their
% edges would count as pia, each with a shell of false depths around it
brain = imfill(A.brainMask, 'holes');

% the brain around the box, 1.2 mm beyond it on every side within the left
% hemisphere: the cortex's thickness there, so the nearest edge of the brain of
% every voxel of the region is in it
pad = 120;
around_ap = max(1, masks.box_ap(1) - pad):min(size(brain, 1), masks.box_ap(end) + pad);
around_dv = max(1, masks.box_dv(1) - pad):min(size(brain, 2), masks.box_dv(end) + pad);
around_ml = max(1, masks.box_ml_left(1) - pad):min(masks.folded_size(3), ...
    masks.box_ml_left(end) + pad);
distance = bwdist(~brain(around_ap, around_dv, around_ml));

% back to the box, in um from the edge: a voxel beside the outside is one voxel
% from the outside voxel's centre, half a voxel from the edge between them
in_ap = masks.box_ap - around_ap(1) + 1;
in_dv = masks.box_dv - around_dv(1) + 1;
in_ml = masks.box_ml_left - around_ml(1) + 1;
anatomy = struct();
anatomy.depth_um = single((distance(in_ap, in_dv, in_ml) - 0.5) * 10);
clear distance
in_brain = brain(masks.box_ap, masks.box_dv, masks.box_ml_left);
anatomy.depth_um(~in_brain) = NaN;
clear brain

% the atlas values and the region in the box
anatomy.annot = A.annot(masks.box_ap, masks.box_dv, masks.box_ml_left);
anatomy.region_lin = masks.band_lin(masks.band_in_region);
anatomy.region = false(masks.box_size);
anatomy.region(anatomy.region_lin) = true;

% the layer of each atlas value: its term at the ontology's substructure level
T_members = readtable(fullfile(atlas_dir, ...
    'parcellation_to_parcellation_term_membership.csv'));
is_layer = strcmp(T_members.parcellation_term_set_name, 'substructure');
anatomy.layer_index = T_members.parcellation_index(is_layer);
anatomy.layer_name = T_members.parcellation_term_acronym(is_layer);

% the crop's first plane as Allen's CCF index, counted from 0 (aplims counts
% the 10 um annotation's planes from 1), and the midline in the folded grid's
% columns (between the last column of the left hemisphere and the first of the
% right)
anatomy.ccf_first_index = A.aplims(1) - 1;
anatomy.midline_ml = masks.folded_size(3) + 0.5;
end

function [ctrl_mice, exp_mice] = load_maps(maps_file, cache_settings, ctrl_dir, ...
    exp_dir, channel)
% The mice's maps from the cache of run_per_mouse_values, after checking that it
% was made with these settings and from the groups' normalised stacks as they
% are now.

if ~exist(maps_file, 'file')
    error(['run_mouse_influence: no maps of run_per_mouse_values in %s. Run ' ...
           'run_per_mouse_values first, with the same comparison and smoothing.'], ...
           maps_file);
end
fprintf('Loading the mice''s maps from %s...\n', maps_file);
S_maps = load(maps_file, 'ctrl_mice', 'exp_mice', 'cache_settings');
if ~isequal(S_maps.cache_settings, cache_settings)
    error(['run_mouse_influence: the maps in %s were made with other regions, box, ' ...
           'smoothing or mice than these settings. Make the settings those of ' ...
           'run_per_mouse_values, or rerun it.'], maps_file);
end
check_maps_cache(S_maps.ctrl_mice, ctrl_dir, channel, maps_file);
check_maps_cache(S_maps.exp_mice, exp_dir, channel, maps_file);
ctrl_mice = S_maps.ctrl_mice;
exp_mice = S_maps.exp_mice;
end

function [ctrl_auto, exp_auto, auto_settings] = load_auto(auto_file, ctrl_mice, ...
    exp_mice, split_settings)
% The mice's autofluorescence from the cache of run_per_mouse_values, after
% checking that it was made from these maps and settings.

if ~exist(auto_file, 'file')
    error(['run_mouse_influence: no autofluorescence of run_per_mouse_values in %s. ' ...
           'Run run_per_mouse_values with auto_control = true first.'], auto_file);
end
fprintf('Loading the mice''s autofluorescence from %s...\n', auto_file);
S_auto = load(auto_file, 'ctrl_auto', 'exp_auto', 'auto_settings');
is_same = isequal(S_auto.auto_settings.split_settings, split_settings) && ...
    isequal(S_auto.ctrl_auto.names, ctrl_mice.names) && ...
    isequal(S_auto.exp_auto.names, exp_mice.names);
if ~is_same
    error(['run_mouse_influence: the autofluorescence in %s was made from other maps, ' ...
           'mice or settings. Rerun run_per_mouse_values.'], auto_file);
end
ctrl_auto = S_auto.ctrl_auto;
exp_auto = S_auto.exp_auto;
auto_settings = S_auto.auto_settings;
end

function cached = load_influence(influence_file, influence_settings, force_recompute)
% The parts of the cache of the folds, the slopes and the counts: those it
% holds when it was made from these maps, autofluorescence and settings, empty
% otherwise.

parts = {'folds', 'folds_held', 'sweep', 'counts'};
cached = struct();
for p = 1:numel(parts)
    cached.(parts{p}) = [];
end
if ~exist(influence_file, 'file') || force_recompute
    return
end
held = who('-file', influence_file);
S_influence = load(influence_file, 'influence_settings');
if ~isequal(S_influence.influence_settings, influence_settings)
    return
end
for p = 1:numel(parts)
    if ismember(parts{p}, held)
        S_part = load(influence_file, parts{p});
        cached.(parts{p}) = S_part.(parts{p});
    end
end
fprintf('Loading %s from %s...\n', strjoin(parts(ismember(parts, held)), ', '), ...
    influence_file);
end

function save_influence(influence_file, cached, influence_settings)
% The cache of the folds, the slopes and the counts, with the settings it was
% made with.

folds = cached.folds;
folds_held = cached.folds_held;
sweep = cached.sweep;
counts = cached.counts;
fprintf('Saving the folds, the slopes and the counts to %s...\n', influence_file);
save(influence_file, 'folds', 'folds_held', 'sweep', 'counts', 'influence_settings', ...
    '-v7.3');
end

function is_there = has_counts(counts, readings)
% Whether the cached counts hold every reading, with the voxels at each level.

is_there = isstruct(counts) && all(isfield(counts, readings));
for r = 1:numel(readings)
    is_there = is_there && isfield(counts.(readings{r}), 'n_at_level');
end
end

% ===== Local functions: folds =====

function folds = fold_influence(ctrl_mice, exp_mice, masks, perm_settings, held)
% Fold 0, the comparison as it is, then each mouse left out in turn, control
% mice first: the fold's alignment, refitted on its mice or, when held is
% given, held at held's slope, intercept and common factor; and
% one_comparison under every split of its mice.

n_ctrl = numel(ctrl_mice.names);
n_mice = n_ctrl + numel(exp_mice.names);
boxes = [ctrl_mice.box(:); exp_mice.box(:)];
profiles = [ctrl_mice.profiles, exp_mice.profiles];
is_exp = (1:n_mice)' > n_ctrl;
n_folds = n_mice + 1;

folds = struct();
folds.left_out = (0:n_mice)';
folds.n_ctrl = zeros(n_folds, 1);
folds.n_exp = zeros(n_folds, 1);
folds.slope = nan(n_folds, 1);
folds.intercept = nan(n_folds, 1);
folds.common_factor = nan(n_folds, 1);
folds.n_splits = zeros(n_folds, 1);
folds.cluster_n = zeros(n_folds, 1);
folds.cluster_mass = zeros(n_folds, 1);
folds.cluster_peak = zeros(n_folds, 1);
folds.negative_n = zeros(n_folds, 1);
folds.negative_mass = zeros(n_folds, 1);
folds.p_positive = nan(n_folds, 1);
folds.p_either = nan(n_folds, 1);
folds.either_sign = zeros(n_folds, 1);
folds.cluster_voxels = cell(n_folds, 1);
folds.null_positive = cell(n_folds, 1);
folds.null_negative = cell(n_folds, 1);
t_start = tic;
for f = 0:n_mice

    % the fold's mice, control mice first, and their alignment: refitted on
    % their groups, or held
    kept = setdiff(1:n_mice, f);
    kept_exp = is_exp(kept);
    if isempty(held)
        [~, ~, ~, slope, intercept, common_factor] = align_exp_to_ctrl( ...
            profiles(:, kept(~kept_exp)), profiles(:, kept(kept_exp)));
        alignment = struct('slope', slope, 'intercept', intercept, ...
            'common_factor', common_factor);
    else
        alignment = held;
    end

    % the test under every split of them
    result = one_comparison(boxes(kept), kept_exp, alignment, masks, perm_settings);

    row = f + 1;
    folds.n_ctrl(row) = nnz(~kept_exp);
    folds.n_exp(row) = nnz(kept_exp);
    folds.slope(row) = alignment.slope;
    folds.intercept(row) = alignment.intercept;
    folds.common_factor(row) = alignment.common_factor;
    folds.n_splits(row) = result.n_splits;
    folds.cluster_n(row) = result.cluster_n;
    folds.cluster_mass(row) = result.cluster_mass;
    folds.cluster_peak(row) = result.cluster_peak;
    folds.negative_n(row) = result.negative_n;
    folds.negative_mass(row) = result.negative_mass;
    folds.p_positive(row) = result.p_positive;
    folds.p_either(row) = result.p_either;
    folds.either_sign(row) = result.either_sign;
    folds.cluster_voxels{row} = result.cluster_voxels;
    folds.null_positive{row} = result.null_positive;
    folds.null_negative{row} = result.null_negative;

    elapsed_min = toc(t_start) / 60;
    fprintf(['  fold %d of %d: %d control and %d experimental mice, slope %.3f, ' ...
             'cluster %d voxels, mass %.1f, p %.4f (positive), %.4f (either sign); ' ...
             '%.1f min\n'], f, n_mice, folds.n_ctrl(row), folds.n_exp(row), ...
            alignment.slope, folds.cluster_n(row), folds.cluster_mass(row), ...
            folds.p_positive(row), folds.p_either(row), elapsed_min);
end
end

function sweep = slope_sweep(ctrl_mice, exp_mice, masks, perm_settings, alignment, ...
    slopes)
% The comparison of all the mice with the alignment's slope at each of slopes,
% its intercept and common factor those of alignment (neither changes a t):
% one_comparison under every split, both signs.

boxes = [ctrl_mice.box(:); exp_mice.box(:)];
n_ctrl = numel(ctrl_mice.names);
is_exp = (1:numel(boxes))' > n_ctrl;
n_slopes = numel(slopes);

sweep = struct();
sweep.slope = slopes(:);
sweep.n_splits = zeros(n_slopes, 1);
sweep.cluster_n = zeros(n_slopes, 1);
sweep.cluster_mass = zeros(n_slopes, 1);
sweep.negative_n = zeros(n_slopes, 1);
sweep.negative_mass = zeros(n_slopes, 1);
sweep.p_positive = nan(n_slopes, 1);
sweep.p_negative = nan(n_slopes, 1);
sweep.p_either = nan(n_slopes, 1);
sweep.either_sign = zeros(n_slopes, 1);
sweep.cluster_voxels = cell(n_slopes, 1);
t_start = tic;
for s = 1:n_slopes
    at_slope = alignment;
    at_slope.slope = slopes(s);
    result = one_comparison(boxes, is_exp, at_slope, masks, perm_settings);
    sweep.n_splits(s) = result.n_splits;
    sweep.cluster_n(s) = result.cluster_n;
    sweep.cluster_mass(s) = result.cluster_mass;
    sweep.negative_n(s) = result.negative_n;
    sweep.negative_mass(s) = result.negative_mass;
    sweep.p_positive(s) = result.p_positive;
    sweep.p_negative(s) = result.p_negative;
    sweep.p_either(s) = result.p_either;
    sweep.either_sign(s) = result.either_sign;
    sweep.cluster_voxels{s} = result.cluster_voxels;

    elapsed_min = toc(t_start) / 60;
    fprintf(['  slope %.2f: cluster %d voxels, mass %.1f, p %.4f; other sign mass ' ...
             '%.1f, p %.4f; either sign p %.4f (%+d); %.1f min\n'], slopes(s), ...
            result.cluster_n, result.cluster_mass, result.p_positive, ...
            result.negative_mass, result.p_negative, result.p_either, ...
            result.either_sign, elapsed_min);
end
end

function result = one_comparison(boxes, is_exp, alignment, masks, perm_settings)
% Step 3's comparison of the given mice on the given alignment:
% region_permutation_test on the band under every split of the mice into
% groups of their sizes; the heaviest cluster of each sign (voxels, mass, peak),
% the p of each sign's mass and step 3's p of either sign, and every split's
% masses.

stack = band_stack(boxes, is_exp, alignment.slope, alignment.intercept, ...
    alignment.common_factor, masks);
[geom, is_cand] = band_geometry(stack, masks, perm_settings);
stack = stack(is_cand, :);
perm = region_permutation_test({stack}, nnz(~is_exp), geom, perm_settings);
clear stack

% the cluster's measure under every split, positive and negative
map = perm.maps{1};
is_cluster = strcmp(perm.measure_names, 'cluster');
observed = perm.splits.observed;
result = struct();
result.null_positive = map.null_pos(:, 1, is_cluster);
result.null_negative = map.null_neg(:, 1, is_cluster);
result.n_splits = numel(result.null_positive);
result.cluster_n = map.detail.cluster_n(1, 1);
result.cluster_mass = result.null_positive(observed);
result.cluster_peak = map.detail.cluster_peak(1, 1);
result.cluster_voxels = map.detail.cluster_voxels{1, 1};
result.negative_n = map.detail.cluster_n(1, 2);
result.negative_mass = result.null_negative(observed);

% each sign's p, the share of the splits whose cluster of that sign is as
% heavy; and step 3's, each split's larger sign against the observed larger
result.p_positive = mean(result.null_positive >= result.cluster_mass);
result.p_negative = mean(result.null_negative >= result.negative_mass);
result.p_either = map.p_perm(1, is_cluster);
result.either_sign = map.sign(1, is_cluster);
end

function check_held_fold_zero(folds, folds_held)
% Fold 0 keeps every mouse, so holding its alignment changes nothing: the same
% cluster, mass and p; a warning where it differs.

is_same = isequal(sort(folds.cluster_voxels{1}(:)), ...
    sort(folds_held.cluster_voxels{1}(:))) && ...
    folds.cluster_mass(1) == folds_held.cluster_mass(1) && ...
    folds.p_either(1) == folds_held.p_either(1);
if is_same
    fprintf('  fold 0 with the alignment held is fold 0 refitted.\n');
else
    warning(['run_mouse_influence: fold 0 with the alignment held (%d voxels, mass ' ...
             '%.2f) differs from fold 0 refitted (%d voxels, mass %.2f).'], ...
             folds_held.cluster_n(1), folds_held.cluster_mass(1), folds.cluster_n(1), ...
             folds.cluster_mass(1));
end
end

function check_fold_zero(folds, step3_table, cluster_region)
% Fold 0, all the mice, must give the region's cluster of run_group_differences,
% its voxels, mass and p, when its table is there; a warning where it differs.

if ~exist(step3_table, 'file')
    fprintf('  no table of run_group_differences to check fold 0 against (%s).\n', ...
        step3_table);
    return
end
T_step3 = readtable(step3_table);
row = strcmp(T_step3.acronym, cluster_region);
step3_p = T_step3.lr_diff_cluster_p_perm(row);
same_p = abs(folds.p_either(1) - step3_p) <= 1e-9;
if T_step3.lr_diff_cluster_sign(row) > 0
    step3_n = T_step3.lr_diff_cluster_n(row);
    step3_mass = T_step3.lr_diff_cluster_score(row);
    same_cluster = folds.cluster_n(1) == step3_n && ...
        abs(folds.cluster_mass(1) - step3_mass) <= 1e-6 * step3_mass;
else
    step3_n = NaN;
    step3_mass = NaN;
    same_cluster = true;
end
if same_p && same_cluster
    fprintf(['  fold 0 gives the cluster of run_group_differences: %d voxels, mass ' ...
             '%.2f, p %.4f.\n'], folds.cluster_n(1), folds.cluster_mass(1), ...
            folds.p_either(1));
else
    warning(['run_mouse_influence: fold 0 gives a %s cluster of %d voxels, mass %.2f, ' ...
             'p %.4f; run_group_differences found %d voxels, mass %.2f, p %.4f (%s).'], ...
             cluster_region, folds.cluster_n(1), folds.cluster_mass(1), ...
             folds.p_either(1), step3_n, step3_mass, step3_p, step3_table);
end
end

function check_against_loo(folds, loo_file)
% Each fold's cluster must be the leave-one-out's of run_per_mouse_values, found
% there on the observed split alone: the same voxels; a warning where one differs.

if ~exist(loo_file, 'file')
    fprintf('  no leave-one-out of run_per_mouse_values to check the folds against (%s).\n', ...
        loo_file);
    return
end
S_loo = load(loo_file, 'loo_clusters');
n_folds = numel(folds.cluster_voxels);
if numel(S_loo.loo_clusters) ~= n_folds
    warning(['run_mouse_influence: %s holds %d folds, this run %d; the folds are not ' ...
             'checked.'], loo_file, numel(S_loo.loo_clusters), n_folds);
    return
end
is_same = false(n_folds, 1);
for f = 1:n_folds
    is_same(f) = isequal(sort(folds.cluster_voxels{f}(:)), sort(S_loo.loo_clusters{f}(:)));
end
if all(is_same)
    fprintf('  every fold''s cluster is the leave-one-out''s, voxel for voxel (%d folds).\n', ...
        n_folds);
else
    warning(['run_mouse_influence: the clusters of folds %s differ from the ' ...
             'leave-one-out''s in %s.'], mat2str(find(~is_same)' - 1), loo_file);
end
end

% ===== Local functions: readings =====

function win = map_window(cluster_lin, masks, anatomy, window_um)
% The window of the maps around the cluster of all the mice, on the box's folded
% grid: its bounding box and window_um on each side, within the box; the
% cluster's planes (the coronal view) and its depth below the pia (the view
% from above); the outlines of the cluster and of the region in both views;
% the axes in um from the cluster's centre.

[ap, dv, ml] = ind2sub(masks.box_size, double(cluster_lin));
margin = round(window_um / 10);
win = struct();
win.ap = max(1, min(ap) - margin):min(masks.box_size(1), max(ap) + margin);
win.dv = max(1, min(dv) - margin):min(masks.box_size(2), max(dv) + margin);
win.ml = max(1, min(ml) - margin):min(masks.box_size(3), max(ml) + margin);
win.planes = min(ap):max(ap);
win.centre = [mean(ap), mean(dv), mean(ml)];
win.box_size = masks.box_size;
win.box_first_plane = masks.box_ap(1);
win.n_cluster = numel(cluster_lin);

% the view from above: in each column of the window, the voxels of the
% cluster's depth below the pia
depth = anatomy.depth_um(cluster_lin);
win.depth_range = [min(depth), max(depth)];
depth_window = anatomy.depth_um(win.ap, :, win.ml);
win.in_band = depth_window >= win.depth_range(1) & depth_window <= win.depth_range(2);

% the outlines: the cluster, and the region (over the cluster's depth, from above)
is_cluster = false(masks.box_size);
is_cluster(cluster_lin) = true;
win.cluster_coronal = squeeze(any(is_cluster(win.planes, win.dv, win.ml), 1));
win.cluster_above = squeeze(any(is_cluster(win.ap, :, win.ml), 2));
win.region_coronal = squeeze(any(anatomy.region(win.planes, win.dv, win.ml), 1));
win.region_above = squeeze(any(anatomy.region(win.ap, :, win.ml) & win.in_band, 2));

% the axes in um from the cluster's centre: ML positive medially, DV ventrally,
% AP posteriorly
win.x_ml_um = (win.ml - win.centre(3)) * 10;
win.y_dv_um = (win.dv - win.centre(2)) * 10;
win.y_ap_um = (win.ap - win.centre(1)) * 10;
end

function [views, values, cluster_values, profiles] = mouse_readings(ctrl_mice, ...
    exp_mice, ctrl_auto, exp_auto, alignment, region_lin, cluster_lin, win)
% Every mouse's three readings (raw, auto, test), control mice first: its coronal
% view and its view from above in the window, its value at each voxel of the
% region, its mean over the cluster, and its profile along AP through the
% cluster's coronal footprint; its coverage along AP, the share of the
% footprint's voxels with a raw value; and, for the counts, its value at each
% voxel of the region on the test's maps at slope 1 (unscaled).

n_ctrl = numel(ctrl_mice.names);
n_mice = n_ctrl + numel(exp_mice.names);
boxes = [ctrl_mice.box(:); exp_mice.box(:)];
boxes_raw = [ctrl_mice.box_raw(:); exp_mice.box_raw(:)];
boxes_auto = [ctrl_auto.box(:); exp_auto.box(:)];

readings = {'raw', 'auto', 'test'};
views = struct();
values = struct();
cluster_values = struct();
profiles = struct();
for r = 1:numel(readings)
    views.(readings{r}).coronal = nan(numel(win.dv), numel(win.ml), n_mice);
    views.(readings{r}).above = nan(numel(win.ap), numel(win.ml), n_mice);
    values.(readings{r}) = nan(numel(region_lin), n_mice, 'single');
    cluster_values.(readings{r}) = nan(n_mice, 1);
    profiles.(readings{r}) = nan(numel(win.ap), n_mice);
end
values.unscaled = nan(numel(region_lin), n_mice, 'single');
profiles.coverage = nan(numel(win.ap), n_mice);
fprintf('Reading every mouse in the window and the region...\n');
for k = 1:n_mice
    is_exp = k > n_ctrl;

    % the raw stack and the autofluorescence, as the asymmetry index
    [lr_diff, lr_sum] = compute_lr_stats(boxes_raw{k});
    [views.raw.coronal(:, :, k), views.raw.above(:, :, k), values.raw(:, k), ...
        cluster_values.raw(k), profiles.raw(:, k)] = read_mouse(abs(lr_diff), lr_sum, ...
        region_lin, cluster_lin, win);
    profiles.coverage(:, k) = footprint_mean(~isnan(lr_diff(win.ap, win.dv, win.ml)), ...
        win.cluster_coronal);
    [lr_diff, lr_sum] = compute_lr_stats(boxes_auto{k});
    [views.auto.coronal(:, :, k), views.auto.above(:, :, k), values.auto(:, k), ...
        cluster_values.auto(k), profiles.auto(:, k)] = read_mouse(abs(lr_diff), ...
        lr_sum, region_lin, cluster_lin, win);

    % the test's maps, as |L - R|
    lr_diff = aligned_box_lr(boxes{k}, is_exp, alignment.slope, alignment.intercept, ...
        alignment.common_factor);
    [views.test.coronal(:, :, k), views.test.above(:, :, k), values.test(:, k), ...
        cluster_values.test(k), profiles.test(:, k)] = read_mouse(abs(lr_diff), [], ...
        region_lin, cluster_lin, win);

    % the test's maps at slope 1, each group on its own scale of step 2
    lr_diff = aligned_box_lr(boxes{k}, is_exp, 1, alignment.intercept, ...
        alignment.common_factor);
    values.unscaled(:, k) = abs(lr_diff(region_lin));
    clear lr_diff lr_sum
end
end

function [coronal, above, region_values, cluster_value, profile] = read_mouse( ...
    abs_diff, lr_sum, region_lin, cluster_lin, win)
% One reading of one mouse: the coronal view over the cluster's planes, the view
% from above over its depth and the profile along AP through its coronal
% footprint, each the ratio of the means along the view (the mean |L - R| alone
% without lr_sum); the value at the region's voxels; and over the cluster.

% |L - R|, averaged along each view and over the cluster's voxels with a value
coronal_diff = squeeze(mean(abs_diff(win.planes, win.dv, win.ml), 1, 'omitnan'));
above_diff = band_mean(abs_diff(win.ap, :, win.ml), win.in_band);
profile_diff = footprint_mean(abs_diff(win.ap, win.dv, win.ml), win.cluster_coronal);
cluster_diff = double(abs_diff(cluster_lin));
has_value = ~isnan(cluster_diff);
if isempty(lr_sum)
    coronal = coronal_diff;
    above = above_diff;
    profile = profile_diff;
    region_values = abs_diff(region_lin);
    cluster_value = mean(cluster_diff(has_value));
    return
end

% L + R the same way, and the ratios; L + R misses the same voxels as |L - R|
coronal_sum = squeeze(mean(lr_sum(win.planes, win.dv, win.ml), 1, 'omitnan'));
above_sum = band_mean(lr_sum(win.ap, :, win.ml), win.in_band);
profile_sum = footprint_mean(lr_sum(win.ap, win.dv, win.ml), win.cluster_coronal);
cluster_sum = double(lr_sum(cluster_lin));
coronal = mean_ratio(coronal_diff, coronal_sum);
above = mean_ratio(above_diff, above_sum);
profile = mean_ratio(profile_diff, profile_sum);
region_values = mean_ratio(abs_diff(region_lin), lr_sum(region_lin));
cluster_value = mean_ratio(mean(cluster_diff(has_value)), mean(cluster_sum(has_value)));
end

function m = footprint_mean(part, footprint)
% The mean of each plane of part (AP x DV x ML) over the voxels of a coronal
% footprint (DV x ML) with a value, NaN where none has one: one value per plane.

columns = reshape(part, size(part, 1), []);
m = mean(double(columns(:, footprint(:))), 2, 'omitnan');
end

function changes = section_changes(group, group_dir, names, win, masks)
% Each mouse's change from one plane of the window to the next in its collected
% stack (unsmoothed, as registered), in the left hemisphere's voxels of the
% cluster's coronal footprint: the mean absolute change over the voxels a
% section reached in both planes, over their mean. The registered volume holds
% each section over the planes nearest it, so the change is about 0 within a
% section and peaks between two: (planes - 1) x mice.

raw_file = fullfile(group_dir, 'nano_4d.mat');
M_raw = matfile(raw_file);
raw_idx = collected_index(names, group, M_raw, 'nano_4d', raw_file);

% the window in the collected stack's grid: the folded grid's columns are those
% of the left hemisphere
planes = win.ap + masks.box_ap(1) - 1;
rows = win.dv + masks.box_dv(1) - 1;
columns = win.ml + masks.box_ml_left(1) - 1;
changes = nan(numel(planes) - 1, numel(names));
fprintf('Sections of %s along AP, from %s...\n', group, raw_file);
for k = 1:numel(names)
    part = double(M_raw.nano_4d(planes, rows, columns, raw_idx(k)));
    part = reshape(part, numel(planes), []);
    part = part(:, win.cluster_coronal(:));

    % a raw 0 is a voxel no section reached
    part(part == 0) = NaN;
    step = abs(diff(part, 1, 1));
    changes(:, k) = mean(step, 2, 'omitnan') / mean(part(:), 'omitnan');
end
end

function m = band_mean(part, in_band)
% The mean of each column of part (AP x DV x ML) over its voxels in the band,
% NaN where none has a value: AP x ML.

part(~in_band) = NaN;
m = squeeze(mean(part, 2, 'omitnan'));
end

function ratio = mean_ratio(mean_diff, mean_sum)
% |L - R| over L + R; NaN where L + R is not above 0.

ratio = mean_diff ./ mean_sum;
ratio(~(mean_sum > 0)) = NaN;
end

function check_against_values(cluster_values, names, values_file)
% Each mouse's raw and autofluorescence asymmetry in the cluster must be its
% selection-matched value of run_per_mouse_values (sm_ai_raw, sm_ai_auto), read
% in the same cluster; a warning where one differs.

if ~exist(values_file, 'file')
    fprintf('  no values of run_per_mouse_values to check the cluster''s against (%s).\n', ...
        values_file);
    return
end
T_values = readtable(values_file);
[is_listed, row] = ismember(names, T_values.mouse);
if ~all(is_listed) || ~ismember('sm_ai_raw', T_values.Properties.VariableNames)
    fprintf('  %s lacks these mice or the selection-matched values; not checked.\n', ...
        values_file);
    return
end
largest_raw = max(abs(cluster_values.raw - T_values.sm_ai_raw(row)));
largest_auto = max(abs(cluster_values.auto - T_values.sm_ai_auto(row)));
if largest_raw <= 1e-12 && largest_auto <= 1e-12
    fprintf(['  every mouse''s asymmetry in the cluster is its selection-matched value ' ...
             'of run_per_mouse_values (raw and autofluorescence).\n']);
else
    warning(['run_mouse_influence: the mice''s asymmetry in the cluster differs from ' ...
             'their selection-matched values in %s by up to %.3g (raw), %.3g (auto).'], ...
             values_file, largest_raw, largest_auto);
end
end

% ===== Local functions: counts =====

function counts = count_splits(values, n_ctrl, region_lin, box_size, connectivity, ...
    cluster_lin)
% One reading's counts: at each complete voxel of the region (every mouse with a
% value), the number of the split's experimental mice above every one of its
% control mice; the voxels at each level (count >= level) and the size of the
% largest connected patch of them, under every split of the mice, the groups
% as they are first; for the observed split, the counts and each level's patch;
% and the voxels expected at each level on exchangeable mice.

n_mice = size(values, 2);
n_exp = n_mice - n_ctrl;

% the complete voxels, and each voxel's mice from the most asymmetric down
is_complete = all(~isnan(values), 2);
complete_lin = region_lin(is_complete);
values = values(is_complete, :);
[~, order] = sort(values, 2, 'descend');
order = uint8(order);
clear values

% every split as a logical row, true for the mice labelled control; nchoosek's
% first is mice 1 to n_ctrl, the groups as they are
ctrl_sets = nchoosek(1:n_mice, n_ctrl);
n_splits = size(ctrl_sets, 1);
in_ctrl = false(n_splits, n_mice);
for s = 1:n_splits
    in_ctrl(s, ctrl_sets(s, :)) = true;
end

counts = struct();
counts.n_region = numel(region_lin);
counts.complete_lin = complete_lin;
counts.cluster_complete = nnz(ismember(cluster_lin, complete_lin));
counts.in_ctrl = in_ctrl;
counts.largest = zeros(n_splits, n_exp);
counts.n_at_level = zeros(n_splits, n_exp);
counts.patch_lin = cell(1, n_exp);
t_start = tic;
for s = 1:n_splits

    % the experimental mice above every control mouse: those ranked before the
    % first control mouse
    labelled_ctrl = in_ctrl(s, :);
    ranked_ctrl = labelled_ctrl(order);
    [~, first_ctrl] = max(ranked_ctrl, [], 2);
    count = uint8(first_ctrl - 1);
    if s == 1
        counts.observed = count;
    end

    % the voxels at each level, and their largest patch
    for k = 1:n_exp
        counts.n_at_level(s, k) = nnz(count >= k);
        [counts.largest(s, k), patch_lin] = largest_patch(count >= k, complete_lin, ...
            box_size, connectivity);
        if s == 1
            counts.patch_lin{k} = patch_lin;
        end
    end
    if s == 1 || mod(s, round(n_splits / 10)) == 0 || s == n_splits
        elapsed_min = toc(t_start) / 60;
        fprintf('  split %d of %d: %.1f min, about %.1f min left\n', s, n_splits, ...
            elapsed_min, elapsed_min / s * (n_splits - s));
    end
end
counts.p = mean(counts.largest >= counts.largest(1, :), 1);
counts.p_at_level = mean(counts.n_at_level >= counts.n_at_level(1, :), 1);

% the voxels expected at each level on exchangeable mice: the top k of a
% voxel's mice are all experimental with probability C(n_exp, k) / C(n_mice, k),
% which is also the mean over every split, the check below
counts.expected_at_level = zeros(1, n_exp);
for k = 1:n_exp
    counts.expected_at_level(k) = numel(complete_lin) * nchoosek(n_exp, k) / ...
        nchoosek(n_mice, k);
end
largest_gap = max(abs(mean(counts.n_at_level, 1) - counts.expected_at_level));
if largest_gap > 1e-6 * numel(complete_lin)
    warning(['run_mouse_influence: the voxels at a level, averaged over every split, ' ...
             'differ from their expectation on exchangeable mice by up to %.1f.'], ...
             largest_gap);
end

% the observed split in the cluster: its voxels at each level, and those of
% each level's largest patch
is_cluster = ismember(complete_lin, cluster_lin);
counts.cluster_at_level = zeros(1, n_exp);
counts.patch_cluster = zeros(1, n_exp);
for k = 1:n_exp
    counts.cluster_at_level(k) = nnz(counts.observed(is_cluster) >= k);
    counts.patch_cluster(k) = numel(intersect(counts.patch_lin{k}, cluster_lin));
end
end

function [n_largest, patch_lin] = largest_patch(is_member, voxel_lin, box_size, ...
    connectivity)
% The size of the largest connected set of the voxels voxel_lin(is_member) in
% the box's grid, and its voxels; 0 and none when there is no voxel.

mask = false(box_size);
mask(voxel_lin(is_member)) = true;
components = bwconncomp(mask, connectivity);
clear mask
patch_lin = zeros(0, 1);
n_largest = 0;
if components.NumObjects == 0
    return
end
sizes = cellfun(@numel, components.PixelIdxList);
[n_largest, largest] = max(sizes);
patch_lin = components.PixelIdxList{largest};
end

function above = above_other_group(values, n_ctrl)
% For each mouse, at each voxel, whether it is above every mouse of the other
% group (NaN compares false).

n_mice = size(values, 2);
is_ctrl = (1:n_mice) <= n_ctrl;
highest_ctrl = max(values(:, is_ctrl), [], 2);
highest_exp = max(values(:, ~is_ctrl), [], 2);
above = false(size(values));
for k = 1:n_mice
    if is_ctrl(k)
        above(:, k) = values(:, k) > highest_exp;
    else
        above(:, k) = values(:, k) > highest_ctrl;
    end
end
end

% ===== Local functions: size and peaks =====

function cluster = cluster_size(cluster_lin, masks, anatomy)
% The cluster in micrometres: its volume, the extent of its voxels along AP, DV
% and ML, its depth below the pia, its widths along the surface and thickness
% through it, its layers, where it sits, and how it compares with one barrel
% column.

% a barrel column of the mouse: about 300 um across, its layer 2/3 from 128 to
% 418 um below the pia (the C2 column, Lefort et al. 2009, Neuron 61:301)
column_diameter_um = 300;
column_l23_top_um = 128;
column_l23_bottom_um = 418;

[ap, dv, ml] = ind2sub(masks.box_size, double(cluster_lin));
cluster = struct();
cluster.n_voxels = numel(cluster_lin);
cluster.volume_mm3 = numel(cluster_lin) * 1e-6;

% planes in the volumes' crop and as Allen's CCF index (from 0), and the extents
% along the axes, from the first voxel's edge to the last's
crop_ap = ap + masks.box_ap(1) - 1;
cluster.plane_first = min(crop_ap);
cluster.plane_last = max(crop_ap);
cluster.ccf_index_first = cluster.plane_first + anatomy.ccf_first_index - 1;
cluster.ccf_index_last = cluster.plane_last + anatomy.ccf_first_index - 1;
cluster.ap_extent_um = (max(ap) - min(ap) + 1) * 10;
cluster.dv_extent_um = (max(dv) - min(dv) + 1) * 10;
cluster.ml_extent_um = (max(ml) - min(ml) + 1) * 10;
crop_ml = ml + masks.box_ml_left(1) - 1;
cluster.centre_from_midline_mm = (anatomy.midline_ml - mean(crop_ml)) * 0.01;

% its depth below the pia
depth = double(anatomy.depth_um(cluster_lin));
cluster.depth_min_um = min(depth);
cluster.depth_median_um = median(depth);
cluster.depth_max_um = max(depth);

% the surface's normal: the depth's mean gradient over the cluster (gradient's
% first output runs along the second dimension, DV, its second along AP)
near_ap = max(1, min(ap) - 2):min(masks.box_size(1), max(ap) + 2);
near_dv = max(1, min(dv) - 2):min(masks.box_size(2), max(dv) + 2);
near_ml = max(1, min(ml) - 2):min(masks.box_size(3), max(ml) + 2);
[grad_dv, grad_ap, grad_ml] = gradient(double(anatomy.depth_um(near_ap, near_dv, near_ml)));
near_lin = sub2ind([numel(near_ap), numel(near_dv), numel(near_ml)], ...
    ap - near_ap(1) + 1, dv - near_dv(1) + 1, ml - near_ml(1) + 1);
normal = [mean(grad_ap(near_lin), 'omitnan'), mean(grad_dv(near_lin), 'omitnan'), ...
    mean(grad_ml(near_lin), 'omitnan')];
normal = normal / norm(normal);
cluster.normal_ap = normal(1);
cluster.normal_dv = normal(2);
cluster.normal_ml = normal(3);

% the widths along the surface, one along AP and one across it, and the
% thickness through it, from the first voxel's edge to the last's
along_ap = [1 0 0] - normal(1) * normal;
along_ap = along_ap / norm(along_ap);
across = cross(normal, along_ap);
positions = [ap, dv, ml] * 10;
cluster.width_along_ap_um = range(positions * along_ap') + 10;
cluster.width_across_um = range(positions * across') + 10;
cluster.thickness_um = range(positions * normal') + 10;

% its layers in the atlas
cluster.layers = layer_shares(anatomy.annot(cluster_lin), anatomy);

% against one barrel column: its layer 2/3 as a cylinder
column_l23_mm3 = pi * (column_diameter_um / 2)^2 * ...
    (column_l23_bottom_um - column_l23_top_um) * 1e-9;
cluster.column_diameter_um = column_diameter_um;
cluster.column_l23_um = sprintf('%d to %d', column_l23_top_um, column_l23_bottom_um);
cluster.column_l23_volume_mm3 = column_l23_mm3;
cluster.share_of_column_l23 = cluster.volume_mm3 / column_l23_mm3;
cluster.widest_over_column = max([cluster.width_along_ap_um, cluster.width_across_um]) / ...
    column_diameter_um;
cluster.column_source = 'Lefort et al. 2009, Neuron 61:301 (mouse C2 barrel column)';
end

function text_layers = layer_shares(annot_values, anatomy)
% The layers of a set of voxels, as 'name share%, ...', the largest first.

layers = layer_of(annot_values, anatomy);
[names, ~, which_layer] = unique(layers);
shares = accumarray(which_layer, 1) / numel(layers);
[shares, order] = sort(shares, 'descend');
names = names(order);
parts = cell(numel(names), 1);
for k = 1:numel(names)
    parts{k} = sprintf('%s %.0f%%', names{k}, 100 * shares(k));
end
text_layers = strjoin(parts, ', ');
end

function layers = layer_of(annot_values, anatomy)
% The atlas layer of each atlas value ('' for a value without one).

[has_layer, row] = ismember(double(annot_values(:)), anatomy.layer_index);
layers = repmat({''}, numel(annot_values), 1);
layers(has_layer) = anatomy.layer_name(row(has_layer));
end

function T_mice = mouse_table(names, n_ctrl, ctrl_type, exp_type, values, ...
    cluster_values, counts, readings, tissue, anatomy, cluster_lin, win)
% One row per mouse, control mice first: its tissue profile's mean and spread
% along AP over the alignment's planes; its values over the cluster; the share
% of the region's and of the cluster's complete voxels where it is above every
% mouse of the other group, per reading (in the cluster, chosen on these mice,
% no chance level applies); and its peak raw asymmetry and, for the
% experimental mice, the peak of its excess over the highest control mouse:
% where each sits, in um from the cluster's centre, below the pia, in which
% layer, and how far from the cluster.

n_mice = numel(names);

% where each mouse is above every mouse of the other group, per reading, over
% the complete voxels of the region and of the cluster
above = struct();
in_cluster = struct();
for r = 1:numel(readings)
    complete_lin = counts.(readings{r}).complete_lin;
    is_complete = ismember(anatomy.region_lin, complete_lin);
    above.(readings{r}) = above_other_group(values.(readings{r})(is_complete, :), n_ctrl);
    in_cluster.(readings{r}) = ismember(complete_lin, cluster_lin);
end

% the raw reading's complete voxels, for the peaks
complete_lin = counts.raw.complete_lin;
raw = values.raw(ismember(anatomy.region_lin, complete_lin), :);
highest_ctrl = max(raw(:, 1:n_ctrl), [], 2);

rows = cell(n_mice, 1);
for k = 1:n_mice
    row = struct();
    row.mouse = names{k};
    if k <= n_ctrl
        row.group = ctrl_type;
    else
        row.group = exp_type;
    end
    row.tissue_mean = tissue.mean(k);
    row.tissue_sd_ap = tissue.sd(k);
    row.cluster_ai_raw = cluster_values.raw(k);
    row.cluster_ai_auto = cluster_values.auto(k);
    row.cluster_abs_diff_test = cluster_values.test(k);
    for r = 1:numel(readings)
        is_above = above.(readings{r})(:, k);
        row.(['above_other_region_' readings{r}]) = mean(is_above);
        row.(['above_other_cluster_' readings{r}]) = mean(is_above(in_cluster.(readings{r})));
    end

    % its peak raw asymmetry
    [row.peak_ai_raw, peak] = max(raw(:, k));
    row = add_place(row, 'peak', complete_lin(peak), anatomy, cluster_lin, win);

    % its peak excess over the highest control mouse, for an experimental mouse
    if k > n_ctrl
        [row.excess_peak_raw, peak] = max(raw(:, k) - highest_ctrl);
        row = add_place(row, 'excess', complete_lin(peak), anatomy, cluster_lin, win);
    else
        row.excess_peak_raw = NaN;
        row = add_place(row, 'excess', [], anatomy, cluster_lin, win);
    end
    rows{k} = row;
end
T_mice = struct2table([rows{:}]', 'AsArray', true);
end

function row = add_place(row, prefix, voxel_lin, anatomy, cluster_lin, win)
% Where a voxel of the box sits: in um from the cluster's centre along AP, DV
% and ML, below the pia, its layer, and its distance from the nearest voxel of
% the cluster (0 inside it); NaN and '' for no voxel.

if isempty(voxel_lin)
    row.([prefix '_ap_um']) = NaN;
    row.([prefix '_dv_um']) = NaN;
    row.([prefix '_ml_um']) = NaN;
    row.([prefix '_depth_um']) = NaN;
    row.([prefix '_layer']) = '';
    row.([prefix '_from_cluster_um']) = NaN;
    return
end
box_size = size(anatomy.region);
[ap, dv, ml] = ind2sub(box_size, double(voxel_lin));
[c_ap, c_dv, c_ml] = ind2sub(box_size, double(cluster_lin));
row.([prefix '_ap_um']) = (ap - win.centre(1)) * 10;
row.([prefix '_dv_um']) = (dv - win.centre(2)) * 10;
row.([prefix '_ml_um']) = (ml - win.centre(3)) * 10;
row.([prefix '_depth_um']) = double(anatomy.depth_um(voxel_lin));
row.([prefix '_layer']) = layer_of(anatomy.annot(voxel_lin), anatomy);
row.([prefix '_layer']) = row.([prefix '_layer']){1};
row.([prefix '_from_cluster_um']) = 10 * min(sqrt((c_ap - ap).^2 + (c_dv - dv).^2 + ...
    (c_ml - ml).^2));
end

% ===== Local functions: tables and summary =====

function T_folds = fold_table(folds, folds_held, names, n_ctrl, ctrl_type, exp_type, win)
% One row per fold: the mouse left out and its group, the fold's mice and
% alignment, its cluster (voxels, volume, mass, peak, planes, centre in um from
% the cluster of all the mice's, its voxels shared with that cluster), the
% cluster of the other sign, and their p; then the same with the alignment held
% at that of all the mice (held_).

n_folds = numel(folds.left_out);
full_cluster = folds.cluster_voxels{1};
rows = cell(n_folds, 1);
for f = 1:n_folds
    row = struct();
    row.fold = folds.left_out(f);
    if folds.left_out(f) == 0
        row.left_out = '';
        row.left_out_group = '';
    else
        row.left_out = names{folds.left_out(f)};
        if folds.left_out(f) <= n_ctrl
            row.left_out_group = ctrl_type;
        else
            row.left_out_group = exp_type;
        end
    end
    row.n_ctrl = folds.n_ctrl(f);
    row.n_exp = folds.n_exp(f);
    row.slope = folds.slope(f);
    row.intercept = folds.intercept(f);
    row.n_splits = folds.n_splits(f);

    % its cluster, and what it shares with the cluster of all the mice
    cluster_lin = folds.cluster_voxels{f};
    n_shared = numel(intersect(cluster_lin, full_cluster));
    row.cluster_n = folds.cluster_n(f);
    row.cluster_volume_mm3 = folds.cluster_n(f) * 1e-6;
    row.cluster_mass = folds.cluster_mass(f);
    row.cluster_peak = folds.cluster_peak(f);
    row.shared_with_all_mice_n = n_shared;
    row.share_of_its_cluster = n_shared / max(numel(cluster_lin), 1);
    row.share_of_all_mice_cluster = n_shared / numel(full_cluster);
    [row.plane_first, row.plane_last, row.centre_ap_um, row.centre_dv_um, ...
        row.centre_ml_um] = cluster_centre(cluster_lin, win);

    % the other sign, and the p
    row.negative_n = folds.negative_n(f);
    row.negative_mass = folds.negative_mass(f);
    row.p_positive = folds.p_positive(f);
    row.n_splits_reaching_positive = round(folds.p_positive(f) * folds.n_splits(f));
    row.p_negative = mean(folds.null_negative{f} >= folds.negative_mass(f));
    row.p_either_sign = folds.p_either(f);
    row.either_sign = folds.either_sign(f);

    % the same with the alignment held
    held_lin = folds_held.cluster_voxels{f};
    row.held_slope = folds_held.slope(f);
    row.held_cluster_n = folds_held.cluster_n(f);
    row.held_cluster_mass = folds_held.cluster_mass(f);
    row.held_shared_with_all_mice_n = numel(intersect(held_lin, full_cluster));
    row.held_share_of_all_mice_cluster = row.held_shared_with_all_mice_n / ...
        numel(full_cluster);
    row.held_negative_n = folds_held.negative_n(f);
    row.held_negative_mass = folds_held.negative_mass(f);
    row.held_p_positive = folds_held.p_positive(f);
    row.held_p_negative = mean(folds_held.null_negative{f} >= folds_held.negative_mass(f));
    row.held_p_either_sign = folds_held.p_either(f);
    row.held_either_sign = folds_held.either_sign(f);
    rows{f} = row;
end
T_folds = struct2table([rows{:}]', 'AsArray', true);
end

function T_slopes = slope_table(sweep, folds, win)
% One row per slope of the alignment, the comparison of all the mice at it:
% the slopes of sweep and the fitted one (fold 0), in order; each sign's
% heaviest cluster, its voxels shared with the cluster at the fitted slope and
% its p, and step 3's p of either sign.

slope = [sweep.slope; folds.slope(1)];
is_fitted = [false(numel(sweep.slope), 1); true];
cluster_n = [sweep.cluster_n; folds.cluster_n(1)];
cluster_mass = [sweep.cluster_mass; folds.cluster_mass(1)];
cluster_voxels = [sweep.cluster_voxels; folds.cluster_voxels(1)];
negative_n = [sweep.negative_n; folds.negative_n(1)];
negative_mass = [sweep.negative_mass; folds.negative_mass(1)];
p_positive = [sweep.p_positive; folds.p_positive(1)];
p_negative = [sweep.p_negative; mean(folds.null_negative{1} >= folds.negative_mass(1))];
p_either_sign = [sweep.p_either; folds.p_either(1)];
either_sign = [sweep.either_sign; folds.either_sign(1)];
n_splits = [sweep.n_splits; folds.n_splits(1)];

full_cluster = folds.cluster_voxels{1};
shared_with_fitted_n = cellfun(@(c) numel(intersect(c, full_cluster)), cluster_voxels);
plane_first = nan(numel(slope), 1);
plane_last = nan(numel(slope), 1);
for s = 1:numel(slope)
    [plane_first(s), plane_last(s)] = cluster_centre(cluster_voxels{s}, win);
end
T_slopes = table(slope, is_fitted, n_splits, cluster_n, cluster_mass, ...
    shared_with_fitted_n, plane_first, plane_last, negative_n, negative_mass, ...
    p_positive, p_negative, p_either_sign, either_sign);
T_slopes = sortrows(T_slopes, 'slope');
end

function [plane_first, plane_last, centre_ap, centre_dv, centre_ml] = ...
    cluster_centre(cluster_lin, win)
% A cluster's first and last plane, in the planes of the volumes' crop, and its
% centre in um from the centre of the cluster of all the mice; NaN for no
% cluster.

plane_first = NaN;
plane_last = NaN;
centre_ap = NaN;
centre_dv = NaN;
centre_ml = NaN;
if isempty(cluster_lin)
    return
end
[ap, dv, ml] = ind2sub(win.box_size, double(cluster_lin));
plane_first = min(ap) + win.box_first_plane - 1;
plane_last = max(ap) + win.box_first_plane - 1;
centre_ap = (mean(ap) - win.centre(1)) * 10;
centre_dv = (mean(dv) - win.centre(2)) * 10;
centre_ml = (mean(ml) - win.centre(3)) * 10;
end

function T_counts = count_table(counts, readings)
% One row per reading and level: the complete voxels; the observed voxels at
% the level, those expected on exchangeable mice and their ratio, the null's
% median and the p; the level's largest patch, what it shares with the
% cluster, the null's median and 95th percentile, and the p.

rows = cell(0, 1);
for r = 1:numel(readings)
    C = counts.(readings{r});
    n_levels = size(C.largest, 2);
    for k = 1:n_levels
        null = C.largest(:, k);
        row = struct();
        row.reading = readings{r};
        row.level = k;
        row.n_region = C.n_region;
        row.n_complete = numel(C.complete_lin);
        row.cluster_complete = C.cluster_complete;
        row.voxels_at_level = C.n_at_level(1, k);
        row.voxels_expected = C.expected_at_level(k);
        row.voxels_over_expected = C.n_at_level(1, k) / C.expected_at_level(k);
        row.voxels_null_median = median(C.n_at_level(:, k));
        row.voxels_p = C.p_at_level(k);
        row.cluster_at_level = C.cluster_at_level(k);
        row.largest_patch = C.largest(1, k);
        row.patch_in_cluster = C.patch_cluster(k);
        row.null_median = median(null);
        row.null_q95 = prctile(null, 95);
        row.n_splits = numel(null);
        row.n_splits_reaching = nnz(null >= null(1));
        row.p = C.p(k);
        rows{end + 1, 1} = row; %#ok<AGROW>
    end
end
T_counts = struct2table([rows{:}]', 'AsArray', true);
end

function T_profiles = profile_table(profiles, plane_changes, win, names, ...
    ccf_first_index)
% One row per plane of the window: the plane in the volumes' crop and as
% Allen's CCF index (from 0), its distance from the cluster's centre, and for
% each mouse its raw, autofluorescence and test profile through the cluster's
% footprint, the share of the footprint with a raw value, and the change of its
% collected stack to the next plane (NaN for the last).

plane = (win.ap + win.box_first_plane - 1)';
T_profiles = table(plane, plane + ccf_first_index - 1, win.y_ap_um(:), ...
    'VariableNames', {'plane', 'ccf_index', 'ap_um'});
short = short_names(names);
for k = 1:numel(names)
    T_profiles.([short{k} '_raw']) = profiles.raw(:, k);
    T_profiles.([short{k} '_auto']) = profiles.auto(:, k);
    T_profiles.([short{k} '_test']) = profiles.test(:, k);
    T_profiles.([short{k} '_coverage']) = profiles.coverage(:, k);
    T_profiles.([short{k} '_section_change']) = [plane_changes(:, k); NaN];
end
end

function print_summary(T_folds, T_slopes, T_counts, T_mice, cluster)
% The folds, the slopes, the counts, the mice and the cluster's size.

fprintf('\nEach mouse left out, the alignment refitted:\n');
disp(T_folds(:, {'fold', 'left_out', 'slope', 'cluster_n', 'cluster_mass', ...
    'share_of_all_mice_cluster', 'centre_ml_um', 'negative_mass', 'p_positive', ...
    'p_either_sign', 'either_sign'}));
fprintf('Each mouse left out, the alignment held at that of all the mice:\n');
disp(T_folds(:, {'fold', 'left_out', 'held_cluster_n', 'held_cluster_mass', ...
    'held_share_of_all_mice_cluster', 'held_negative_mass', 'held_p_positive', ...
    'held_p_either_sign', 'held_either_sign'}));
fprintf('All the mice at each slope of the alignment:\n');
disp(T_slopes);
fprintf('Experimental mice above every control mouse, voxels and largest patch per level:\n');
disp(T_counts(:, {'reading', 'level', 'n_complete', 'voxels_at_level', ...
    'voxels_over_expected', 'voxels_p', 'largest_patch', 'patch_in_cluster', ...
    'null_median', 'null_q95', 'p'}));
fprintf('The mice:\n');
disp(T_mice(:, {'mouse', 'group', 'tissue_mean', 'tissue_sd_ap', 'cluster_ai_raw', ...
    'cluster_ai_auto', 'above_other_region_raw', 'above_other_cluster_raw', ...
    'above_other_region_auto', 'peak_ai_raw', 'peak_depth_um', 'peak_layer', ...
    'peak_from_cluster_um', 'excess_peak_raw', 'excess_from_cluster_um'}));
fprintf(['The cluster: %d voxels, %.4f mm^3; planes %d to %d (CCF index %d to %d); ' ...
         '%d x %d x %d um along AP, DV and ML; %.0f to %.0f um below the pia ' ...
         '(median %.0f); along the surface %.0f um (AP) by %.0f um, %.0f um thick; ' ...
         '%s; %.2f of the layer 2/3 of a barrel column, its widest %.2f of the ' ...
         'column''s diameter.\n'], cluster.n_voxels, cluster.volume_mm3, ...
        cluster.plane_first, cluster.plane_last, cluster.ccf_index_first, ...
        cluster.ccf_index_last, cluster.ap_extent_um, cluster.dv_extent_um, ...
        cluster.ml_extent_um, cluster.depth_min_um, cluster.depth_max_um, ...
        cluster.depth_median_um, cluster.width_along_ap_um, cluster.width_across_um, ...
        cluster.thickness_um, cluster.layers, cluster.share_of_column_l23, ...
        cluster.widest_over_column);
end

% ===== Local functions: figures =====

function plot_folds(T_folds, T_slopes, comparison, file_tag, comp_out_dir)
% The folds, one bar each, named by the mouse left out and coloured by its
% group: the cluster's voxels with those it shares with the cluster of all the
% mice, its mass, its p, each with the alignment refitted (bars, dots) and held
% (diamonds); all the mice at each slope of the alignment, both signs' mass and
% p; and the notes.

labels = [{'none'}; short_names(T_folds.left_out(2:end))];
[bar_colours, inner_colours] = fold_colours(T_folds, comparison);
fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'Normalized', ...
    'Position', [0 0 1 1]);

% the cluster's voxels, those shared with the cluster of all the mice inside
subplot(2, 3, 1);
draw_fold_bars(T_folds.cluster_n, T_folds.shared_with_all_mice_n, ...
    T_folds.held_cluster_n, bar_colours, inner_colours, labels);
ylabel('voxels of 10 um (1,000 voxels = 0.001 mm^3)', 'FontSize', 9);
title({'Cluster of the mice kept: its voxels', ['\rm\fontsize{9}inner bars: shared ' ...
       'with all the mice''s; diamonds: the alignment held']}, 'FontSize', 10);

% its mass
subplot(2, 3, 2);
draw_fold_bars(T_folds.cluster_mass, [], T_folds.held_cluster_mass, bar_colours, [], ...
    labels);
ylabel('cluster mass (summed surprise)', 'FontSize', 9);
title({'Its mass, step 3''s score', '\rm\fontsize{9}diamonds: the alignment held'}, ...
    'FontSize', 10);

% its p
subplot(2, 3, 3);
draw_fold_p(T_folds, bar_colours, labels);

% all the mice at each slope: the mass of each sign's heaviest cluster, and its p
subplot(2, 3, 4);
draw_slope_masses(T_slopes, T_folds, comparison);
subplot(2, 3, 5);
draw_slope_p(T_slopes, comparison);

% the notes
subplot(2, 3, 6);
axis off;
text(-0.05, 1.02, fold_notes(T_folds, comparison), 'Units', 'normalized', ...
    'VerticalAlignment', 'top', 'FontSize', 7, 'Interpreter', 'none', ...
    'FontName', 'FixedWidth');

title_line = sprintf(['Each mouse left out, and the alignment''s slope: the %s cluster ' ...
                      'where |L - R| is higher in %s - %s'], comparison.cluster_region, ...
                      comparison.exp_type, strrep(file_tag, '_', ' '));
sgtitle(title_line, 'FontSize', 14, 'FontWeight', 'bold');
saveas(fig, fullfile(comp_out_dir, ['Influence_Folds_' file_tag '.fig']));
exportgraphics(fig, fullfile(comp_out_dir, ['Influence_Folds_' file_tag '.png']), ...
    'Resolution', 300);
end

function [bar_colours, inner_colours] = fold_colours(T_folds, comparison)
% A fold's bar in the colour of the left-out mouse's group, grey for all the
% mice; the inner bars in the group's darker colour.

n_folds = height(T_folds);
bar_colours = repmat(sep_palette('paired_lines'), n_folds, 1);
inner_colours = repmat(0.6 * sep_palette('paired_lines'), n_folds, 1);
is_ctrl = strcmp(T_folds.left_out_group, comparison.ctrl_type);
is_exp = strcmp(T_folds.left_out_group, comparison.exp_type);
bar_colours(is_ctrl, :) = repmat(sep_palette('control'), nnz(is_ctrl), 1);
bar_colours(is_exp, :) = repmat(sep_palette('experimental'), nnz(is_exp), 1);
inner_colours(is_ctrl, :) = repmat(sep_palette('control_mean'), nnz(is_ctrl), 1);
inner_colours(is_exp, :) = repmat(sep_palette('experimental_mean'), nnz(is_exp), 1);
end

function draw_fold_bars(heights, inner, held, bar_colours, inner_colours, labels)
% One bar per fold, with an inner, narrower bar when given and a diamond at the
% held alignment's value, and a dotted line at the first fold's height, all the
% mice's.

hold on;
box on;
grid on;
n_folds = numel(heights);
for f = 1:n_folds
    bar(f, heights(f), 0.7, 'FaceColor', bar_colours(f, :), 'EdgeColor', 'none');
    if ~isempty(inner)
        bar(f, inner(f), 0.3, 'FaceColor', inner_colours(f, :), 'EdgeColor', 'none');
    end
end
scatter(1:n_folds, held, 40, 'd', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', ...
    'LineWidth', 1);
yline(heights(1), ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 0.8);
xlim([0.4, n_folds + 0.6]);
xticks(1:n_folds);
xticklabels(labels);
xtickangle(45);
xlabel('mouse left out', 'FontSize', 9);
set(gca, 'FontSize', 9);
end

function draw_fold_p(T_folds, bar_colours, labels)
% Each fold's p on a log axis: the positive cluster's filled, either sign's
% open (a square where its larger sign is the control group's), the positive
% cluster's with the alignment held a diamond; a dashed line at 0.05 and,
% under each of the first two, a dotted one at its smallest possible p.

hold on;
box on;
grid on;
n_folds = height(T_folds);
for f = 1:n_folds

    % the smallest p: one split for the positive cluster; for either sign, two
    % with groups of equal size, a split and its mirror being as heavy
    n_smallest_either = 1 + (T_folds.n_ctrl(f) == T_folds.n_exp(f));
    plot(f - 0.2 + [-0.12 0.12], [1 1] / T_folds.n_splits(f), ':', ...
        'Color', [0.5 0.5 0.5]);
    plot(f + [-0.12 0.12], [1 1] * n_smallest_either / T_folds.n_splits(f), ':', ...
        'Color', [0.5 0.5 0.5]);
    scatter(f - 0.2, T_folds.p_positive(f), 45, bar_colours(f, :), 'filled');
    if T_folds.either_sign(f) < 0
        marker = 's';
    else
        marker = 'o';
    end
    scatter(f, T_folds.p_either_sign(f), 45, bar_colours(f, :), marker, 'LineWidth', 1.2);
    scatter(f + 0.2, T_folds.held_p_positive(f), 40, 'd', 'MarkerFaceColor', 'w', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 1);
end
yline(0.05, '--', 'Color', [0.4 0.4 0.4]);
set(gca, 'YScale', 'log', 'FontSize', 9);
ylim([2e-3 1.2]);
yticks([0.005 0.01 0.02 0.05 0.1 0.2 0.5 1]);
xlim([0.4, n_folds + 0.6]);
xticks(1:n_folds);
xticklabels(labels);
xtickangle(45);
xlabel('mouse left out', 'FontSize', 9);
ylabel('p over every split of the mice kept', 'FontSize', 9);
title({'p of the cluster''s mass', ['\rm\fontsize{9}filled: the positive cluster; ' ...
       'open: either sign (step 3''s); diamonds: positive, the alignment held']}, ...
       'FontSize', 10);
end

function draw_slope_masses(T_slopes, T_folds, comparison)
% All the mice at each slope of the alignment: the mass of the heaviest cluster
% where |L - R| is higher in each group; the fitted slope dashed, and the
% slopes refitted without each mouse as triangles along the bottom, in the
% colour of its group.

hold on;
box on;
grid on;
plot(T_slopes.slope, T_slopes.cluster_mass, '-o', 'Color', ...
    sep_palette('experimental_mean'), 'MarkerFaceColor', ...
    sep_palette('experimental_mean'), 'LineWidth', 1.4, 'MarkerSize', 5);
plot(T_slopes.slope, T_slopes.negative_mass, '-o', 'Color', sep_palette('control_mean'), ...
    'MarkerFaceColor', sep_palette('control_mean'), 'LineWidth', 1.4, 'MarkerSize', 5);
y_top = 1.1 * max([T_slopes.cluster_mass; T_slopes.negative_mass]);
fitted = T_slopes.slope(T_slopes.is_fitted);
plot([fitted fitted], [0 y_top], '--', 'Color', [0.4 0.4 0.4]);
[bar_colours, ~] = fold_colours(T_folds, comparison);
for f = 2:height(T_folds)
    scatter(T_folds.slope(f), 0.03 * y_top, 40, bar_colours(f, :), '^', 'filled');
end
ylim([0 y_top]);
xlim([min(T_slopes.slope) - 0.05, max(T_slopes.slope) + 0.05]);
set(gca, 'FontSize', 9);
xlabel(sprintf('slope of the alignment (%s |L - R| times it; 1: each group on its own scale)', ...
    comparison.exp_type), 'FontSize', 9);
ylabel('cluster mass (summed surprise)', 'FontSize', 9);
legend({sprintf('%s higher', comparison.exp_type), sprintf('%s higher', ...
    comparison.ctrl_type)}, 'Location', 'northwest', 'FontSize', 8);
title({'All the mice at each slope of the alignment', ['\rm\fontsize{9}dashed: the ' ...
       'fitted slope; triangles: refitted without each mouse']}, 'FontSize', 10);
end

function draw_slope_p(T_slopes, comparison)
% All the mice at each slope: the p of each sign's heaviest cluster over every
% split, on a log axis, a dashed line at 0.05 and a dotted one at one split.

hold on;
box on;
grid on;
plot(T_slopes.slope, T_slopes.p_positive, '-o', 'Color', ...
    sep_palette('experimental_mean'), 'MarkerFaceColor', ...
    sep_palette('experimental_mean'), 'LineWidth', 1.4, 'MarkerSize', 5);
plot(T_slopes.slope, T_slopes.p_negative, '-o', 'Color', sep_palette('control_mean'), ...
    'MarkerFaceColor', sep_palette('control_mean'), 'LineWidth', 1.4, 'MarkerSize', 5);
fitted = T_slopes.slope(T_slopes.is_fitted);
plot([fitted fitted], [2e-3 1.2], '--', 'Color', [0.4 0.4 0.4]);
yline(0.05, '--', 'Color', [0.4 0.4 0.4]);
yline(1 / T_slopes.n_splits(1), ':', 'Color', [0.5 0.5 0.5]);
set(gca, 'YScale', 'log', 'FontSize', 9);
ylim([2e-3 1.2]);
yticks([0.005 0.01 0.02 0.05 0.1 0.2 0.5 1]);
xlim([min(T_slopes.slope) - 0.05, max(T_slopes.slope) + 0.05]);
xlabel('slope of the alignment', 'FontSize', 9);
ylabel(sprintf('p over every split (%d)', T_slopes.n_splits(1)), 'FontSize', 9);
legend({sprintf('%s higher', comparison.exp_type), sprintf('%s higher', ...
    comparison.ctrl_type)}, 'Location', 'northeast', 'FontSize', 8);
title({'Each sign''s p at each slope', ['\rm\fontsize{9}the splits keep the slope: ' ...
       'step 3''s p does not carry its uncertainty']}, 'FontSize', 10);
end

function lines = fold_notes(T_folds, comparison)
% The notes of the folds' figure, from the settings and the run.

settings = comparison.perm_settings;
split_counts = strjoin(arrayfun(@num2str, unique(T_folds.n_splits)', ...
    'UniformOutput', false), ', ');
lines = {
    sprintf(['Each fold: step 3''s comparison of the mice kept, the alignment (the line ' ...
             'from the %s'], comparison.exp_type)
    sprintf(['  profile onto the %s one) refitted on them; a t where each group has %d ' ...
             'mice'], comparison.ctrl_type, settings.min_mice_per_group)
    sprintf(['  with a value; its surprise''s median over +/- %d planes; voxels at ' ...
             'p < %g,'], settings.slab_range, settings.cluster_p)
    sprintf(['  %d-connected within %s. The cluster: the heaviest where |L - R| is ' ...
             'higher'], settings.cluster_connectivity, comparison.cluster_region)
    sprintf('  in %s.', comparison.exp_type)
    sprintf(['Held: the same with the alignment of all the mice. On the test''s maps ' ...
             'the slope'])
    sprintf(['  multiplies every %s mouse''s |L - R| (the intercept cancels in L - R, ' ...
             'the'], comparison.exp_type)
    '  common factor in the t), so a refitted slope moves all of them at once.'
    sprintf(['p positive: the share of every split of the mice kept into groups of ' ...
             'their sizes'])
    sprintf(['  (%s splits) whose positive cluster is as heavy, the observed split ' ...
             'included.'], split_counts)
    'p either sign: step 3''s, each split''s larger sign against the observed larger;'
    sprintf('  an open square where the observed larger is the %s-higher cluster.', ...
        comparison.ctrl_type)
    'Dotted: the smallest p possible, one split for the positive p; for either sign,'
    '  two with groups of equal size, a split and its mirror being as heavy.'
    'Run_normalise_groups is not redone without the mouse (README).'
    ''
    sprintf('%-15s %5s %6s %5s %6s %6s  | %13s %5s %6s', '', 'slope', 'voxels', ...
        'kept', 'p', 'either', 'held: voxels', 'kept', 'p')
    };

% the folds, in numbers
for f = 1:height(T_folds)
    if T_folds.fold(f) == 0
        name = 'all the mice';
    else
        name = ['without ' strtok(T_folds.left_out{f}, '_')];
    end
    lines{end + 1} = sprintf('%-15s %5.2f %6d %4.0f%% %6.3f %6.3f%s | %13d %4.0f%% %6.3f', ...
        name, T_folds.slope(f), T_folds.cluster_n(f), ...
        100 * T_folds.share_of_all_mice_cluster(f), T_folds.p_positive(f), ...
        T_folds.p_either_sign(f), sign_note(T_folds.either_sign(f)), ...
        T_folds.held_cluster_n(f), 100 * T_folds.held_share_of_all_mice_cluster(f), ...
        T_folds.held_p_positive(f)); %#ok<AGROW>
end
lines{end + 1} = sprintf(['kept: the share of the cluster of all the mice; * either sign''s ' ...
                          'larger is %s higher'], comparison.ctrl_type);
end

function note = sign_note(either_sign)
% A mark after a fold's either-sign p when its larger sign is the control
% group's.

note = ' ';
if either_sign < 0
    note = '*';
end
end

function plot_maps(views, cluster_values, win, comparison, file_tag, comp_out_dir)
% The mice's maps around the cluster, one figure per reading: both asymmetry
% indices (raw and auto) on one colour scale, the test's |L - R| on its own,
% each up to the 99th percentile of what its figures show.

pooled = [views.raw.coronal(:); views.raw.above(:); views.auto.coronal(:); ...
    views.auto.above(:)];
ai_limit = prctile(pooled(~isnan(pooled)), 99);
pooled = [views.test.coronal(:); views.test.above(:)];
test_limit = prctile(pooled(~isnan(pooled)), 99);

readings = {
    'raw',  'asymmetry |L - R| / (L + R), raw stack less its off-tissue level', ...
        'AI', ai_limit, 'the raw and the autofluorescence AI'
    'auto', 'asymmetry |L - R| / (L + R), autofluorescence less its off-tissue level', ...
        'AI', ai_limit, 'the raw and the autofluorescence AI'
    'test', '|L - R| on the test''s maps (common scale)', '|L - R|', test_limit, ...
        'this figure'
    };
for r = 1:size(readings, 1)
    reading = readings{r, 1};
    style = struct('description', readings{r, 2}, 'value_name', readings{r, 3}, ...
        'limit', readings{r, 4}, 'limit_of', readings{r, 5});
    draw_map_figure(views.(reading), cluster_values.(reading), win, comparison, style, ...
        sprintf('Influence_Maps_%s_%s', file_tag, reading), comp_out_dir);
end
end

function draw_map_figure(view, cluster_value, win, comparison, style, file_name, ...
    comp_out_dir)
% One reading's maps: per group a row of coronal views over the cluster's planes
% and a row of views from above over its depth, a mouse per column, control
% mice in the first two rows; the cluster outlined in blue, the region dotted.

% a column per mouse of a group, not the three of docs/STYLE.md, so a row is
% one group and one view
n_ctrl = comparison.n_ctrl;
n_mice = comparison.n_mice;
n_cols_grid = max(n_ctrl, n_mice - n_ctrl);
fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'Normalized', ...
    'Position', [0 0 1 1]);
names = short_names(comparison.names);
for k = 1:n_mice
    if k <= n_ctrl
        first_row = 0;
        column = k;
        group = comparison.ctrl_type;
        name_colour = sep_palette('control_mean');
    else
        first_row = 2;
        column = k - n_ctrl;
        group = comparison.exp_type;
        name_colour = sep_palette('experimental_mean');
    end

    % the coronal view
    subplot(4, n_cols_grid, first_row * n_cols_grid + column);
    draw_view(view.coronal(:, :, k), win.x_ml_um, win.y_dv_um, win.cluster_coronal, ...
        win.region_coronal, style.limit, sep_palette('intensity'));
    title(sprintf('%s (%s): cluster %.3f', names{k}, group, cluster_value(k)), ...
        'FontSize', 9, 'Color', name_colour, 'Interpreter', 'none');
    if column == 1
        ylabel('coronal: DV from centre (um)', 'FontSize', 8);
    end

    % the view from above
    ax = subplot(4, n_cols_grid, (first_row + 1) * n_cols_grid + column);
    draw_view(view.above(:, :, k), win.x_ml_um, win.y_ap_um, win.cluster_above, ...
        win.region_above, style.limit, sep_palette('intensity'));
    xlabel('ML from centre (um), medial >', 'FontSize', 8);
    if column == 1
        ylabel('from above: AP from centre (um)', 'FontSize', 8);
    end
end

% one colour bar for every panel
cb = colorbar(ax, 'Position', [0.93 0.3 0.01 0.4]);
cb.Label.String = style.value_name;
cb.FontSize = 9;

% the title: the reading, the views, the outlines and the scale
planes_crop = win.planes([1 end]) + win.box_first_plane - 1;
planes_ccf = planes_crop + comparison.ccf_first_index - 1;
title_line = sprintf('Each mouse''s %s around the %s cluster - %s', style.description, ...
    comparison.cluster_region, strrep(file_name, '_', ' '));
view_line = sprintf(['Rows: %s coronal, from above, %s coronal, from above. Coronal: ' ...
                     'the mean over the cluster''s planes %d to %d (CCF index %d to %d); ' ...
                     'from above: over %.0f to %.0f um below the pia, the cluster''s ' ...
                     'depth.'], comparison.ctrl_type, comparison.exp_type, planes_crop, ...
                     planes_ccf, win.depth_range);
outline_line = sprintf(['Blue: the cluster of all %d mice; dotted grey: %s; grey: no ' ...
                        'value. Titles: the mouse''s %s over the cluster. Colour: 0 to ' ...
                        '%.3g, the 99th percentile of %s.'], n_mice, ...
                        comparison.cluster_region, style.value_name, style.limit, ...
                        style.limit_of);
sgtitle({title_line, ['\rm\fontsize{10}' view_line], ['\rm\fontsize{10}' outline_line]}, ...
    'FontSize', 13, 'FontWeight', 'bold');
saveas(fig, fullfile(comp_out_dir, [file_name '.fig']));
exportgraphics(fig, fullfile(comp_out_dir, [file_name '.png']), 'Resolution', 300);
end

function draw_view(img, x_um, y_um, cluster_mask, region_mask, limit, colour_map)
% One view: the image on the no-data grey, its axes in um, the region's outline
% dotted and the cluster's in blue.

imagesc(x_um, y_um, img, 'AlphaData', double(~isnan(img)));
set(gca, 'Color', sep_palette('no_data'), 'FontSize', 7, 'TickDir', 'out');
axis image;
clim([0 limit]);
colormap(gca, colour_map);
hold on;
draw_outline(x_um, y_um, region_mask, sep_palette('paired_lines'), 0.8, ':');
draw_outline(x_um, y_um, cluster_mask, sep_palette('outline'), 1.2, '-');
end

function draw_outline(x_um, y_um, mask, colour, width, style)
% The edge of a mask, as a contour at one half; nothing for an empty or full mask.

if ~any(mask(:)) || all(mask(:))
    return
end
contour(x_um, y_um, double(mask), [0.5 0.5], 'LineColor', colour, 'LineWidth', width, ...
    'LineStyle', style);
end

function plot_consistency(counts, readings, T_mice, win, masks, comparison, file_tag, ...
    comp_out_dir)
% The counts, a row per reading: the number of experimental mice above every
% control mouse, averaged over the cluster's planes and over its depth; the
% largest patch at each level against every split, with the voxels at the level
% over their expectation; each mouse's share of the region, and of the cluster,
% where it is above every mouse of the other group.

n_exp = comparison.n_mice - comparison.n_ctrl;
names = short_names(comparison.names);
count_colours = sep_palette('counts');
descriptions = struct('raw', 'raw AI', 'test', 'test''s |L - R|', 'unscaled', ...
    'test''s |L - R| at slope 1', 'auto', 'autofluorescence AI');
fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'Normalized', ...
    'Position', [0 0 1 1]);

% four columns, not the three of docs/STYLE.md, so a row is one reading
n_rows_grid = numel(readings);
for r = 1:n_rows_grid
    C = counts.(readings{r});

    % the counts in the box, NaN off the complete voxels, and their two views
    count_box = nan(masks.box_size, 'single');
    count_box(C.complete_lin) = C.observed;
    coronal = squeeze(mean(count_box(win.planes, win.dv, win.ml), 1, 'omitnan'));
    above = band_mean(count_box(win.ap, :, win.ml), win.in_band);
    clear count_box

    ax = subplot(n_rows_grid, 4, (r - 1) * 4 + 1);
    draw_view(coronal, win.x_ml_um, win.y_dv_um, win.cluster_coronal, ...
        win.region_coronal, n_exp, count_colours);
    title(sprintf('%s: coronal, cluster''s planes', descriptions.(readings{r})), ...
        'FontSize', 9);
    ylabel('DV from centre (um)', 'FontSize', 8);
    subplot(n_rows_grid, 4, (r - 1) * 4 + 2);
    draw_view(above, win.x_ml_um, win.y_ap_um, win.cluster_above, win.region_above, ...
        n_exp, count_colours);
    title(sprintf('%s: from above, cluster''s depth', descriptions.(readings{r})), ...
        'FontSize', 9);
    ylabel('AP from centre (um)', 'FontSize', 8);
    if r == n_rows_grid
        xlabel('ML from centre (um), medial >', 'FontSize', 8);
    end

    % the largest patch at each level against every split
    subplot(n_rows_grid, 4, (r - 1) * 4 + 3);
    draw_count_null(C, comparison);
    title(sprintf('%s: largest patch at each level', descriptions.(readings{r})), ...
        'FontSize', 9);

    % each mouse's share above every mouse of the other group, over the region
    % and in the cluster
    subplot(n_rows_grid, 4, (r - 1) * 4 + 4);
    draw_mouse_shares(T_mice.(['above_other_region_' readings{r}]), ...
        T_mice.(['above_other_cluster_' readings{r}]), names, comparison);
    title(sprintf('%s: share above the other group', descriptions.(readings{r})), ...
        'FontSize', 9);
end
cb = colorbar(ax, 'Position', [0.045 0.35 0.008 0.3]);
cb.Label.String = sprintf('%s mice above every %s mouse (mean)', comparison.exp_type, ...
    comparison.ctrl_type);
cb.FontSize = 8;

% the title: what is counted, the null, the chance of a mouse's share
raw = counts.raw;
title_line = sprintf(['How many %s mice are above every %s mouse, voxel by voxel, in %s ' ...
                      '- %s'], comparison.exp_type, comparison.ctrl_type, ...
                      comparison.cluster_region, strrep(file_tag, '_', ' '));
count_line = sprintf(['Counted at the voxels where all %d mice have a value (raw: %d of ' ...
                      '%d, %d of the cluster''s %d); views: the mean count, blue the ' ...
                      'cluster of all the mice, dotted grey %s. Slope 1: each group on ' ...
                      'its own scale of step 2.'], comparison.n_mice, ...
                      numel(raw.complete_lin), raw.n_region, raw.cluster_complete, ...
                      win.n_cluster, comparison.cluster_region);
null_line = sprintf(['Third column: the largest connected patch (%d-connected) where at ' ...
                     'least k mice are above, under each of the %d splits (grey) and ' ...
                     'observed (red), its p; x: the voxels at the level over those ' ...
                     'expected on exchangeable mice, its p.'], ...
                     comparison.perm_settings.cluster_connectivity, size(raw.largest, 1));
share_line = sprintf(['Fourth column: bars, the share of the region''s voxels where a ' ...
                      'mouse is above every mouse of the other group, chance %.2f on ' ...
                      'exchangeable mice (dashed); dots, the same in the cluster, chosen ' ...
                      'on these mice, so no chance level applies.'], 1 / (n_exp + 1));
sgtitle({title_line, ['\rm\fontsize{10}' count_line], ['\rm\fontsize{10}' null_line], ...
    ['\rm\fontsize{10}' share_line]}, 'FontSize', 13, 'FontWeight', 'bold');
saveas(fig, fullfile(comp_out_dir, ['Influence_Consistency_' file_tag '.fig']));
exportgraphics(fig, fullfile(comp_out_dir, ['Influence_Consistency_' file_tag '.png']), ...
    'Resolution', 300);
end

function draw_count_null(C, comparison)
% At each level, every split's largest patch as grey dots, the observed one as
% a red dot, and above it its p and the voxels at the level over their
% expectation with that ratio's p; log axis of the voxels plus one, so a split
% without a patch shows.

hold on;
box on;
grid on;
[n_splits, n_levels] = size(C.largest);
jitter = ((1:n_splits)' / n_splits - 0.5) * 0.5;
for k = 1:n_levels
    scatter(k + jitter, C.largest(:, k) + 1, 8, sep_palette('paired_lines'), 'filled', ...
        'MarkerFaceAlpha', 0.5);
    scatter(k, C.largest(1, k) + 1, 50, sep_palette('experimental_mean'), 'filled');
    over_expected = C.n_at_level(1, k) / C.expected_at_level(k);
    text(k, 3 * (max(C.largest(:, k)) + 1), {sprintf('p %.3f', C.p(k)), ...
        sprintf('x%.2f, p %.3f', over_expected, C.p_at_level(k))}, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 6);
end
set(gca, 'YScale', 'log', 'FontSize', 8);
xlim([0.4, n_levels + 0.6]);
xticks(1:n_levels);
ylim([0.7, 300 * (max(C.largest(:)) + 1)]);
xlabel(sprintf('level: at least k %s mice above every %s mouse', comparison.exp_type, ...
    comparison.ctrl_type), 'FontSize', 8);
ylabel('largest patch (voxels + 1)', 'FontSize', 8);
end

function draw_mouse_shares(region_shares, cluster_shares, names, comparison)
% Each mouse's share of the region's complete voxels where it is above every
% mouse of the other group, a bar coloured by its group with a dashed line at
% the chance of one mouse among the other group's; and the same in the
% cluster, an open dot, without a chance level since the cluster was chosen on
% these mice.

hold on;
box on;
grid on;
n_mice = numel(region_shares);
n_ctrl = comparison.n_ctrl;
for k = 1:n_mice
    if k <= n_ctrl
        colour = sep_palette('control');
        n_other = n_mice - n_ctrl;
    else
        colour = sep_palette('experimental');
        n_other = n_ctrl;
    end
    bar(k, region_shares(k), 0.7, 'FaceColor', colour, 'EdgeColor', 'none');
    plot(k + [-0.4 0.4], [1 1] / (n_other + 1), '--', 'Color', [0.4 0.4 0.4]);
end
scatter(1:n_mice, cluster_shares, 30, 'o', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', ...
    'k', 'LineWidth', 1);
ylim([0 1]);
xlim([0.4, n_mice + 0.6]);
xticks(1:n_mice);
xticklabels(names);
xtickangle(45);
set(gca, 'FontSize', 8);
ylabel('share of the voxels', 'FontSize', 8);
end

function plot_profiles(profiles, plane_changes, win, comparison, file_tag, ...
    comp_out_dir)
% Each mouse's asymmetry along AP through the cluster's coronal footprint, a
% panel per mouse, control mice in the first row: its raw AI, its
% autofluorescence AI dashed, the highest and lowest control mouse's raw AI in
% grey, the cluster's planes shaded, the planes where part of the footprint has
% no tissue in the no-data grey (the darker, the less tissue), and at the
% bottom the change of its collected stack from plane to plane, whose peaks
% are the edges between its sections.

% a column per mouse of a group, not the three of docs/STYLE.md, so a row is
% one group
n_ctrl = comparison.n_ctrl;
n_mice = comparison.n_mice;
n_cols_grid = max(n_ctrl, n_mice - n_ctrl);
names = short_names(comparison.names);
ap_um = win.y_ap_um(:);
edge_um = (ap_um(1:end - 1) + ap_um(2:end)) / 2;
ctrl_low = min(profiles.raw(:, 1:n_ctrl), [], 2);
ctrl_high = max(profiles.raw(:, 1:n_ctrl), [], 2);
y_top = 1.1 * max([profiles.raw(:); profiles.auto(:)]);
cluster_edges_um = (win.planes([1 end]) - win.centre(1)) * 10 + [-5 5];
fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'Normalized', ...
    'Position', [0 0 1 0.7]);
for k = 1:n_mice
    if k <= n_ctrl
        panel = k;
        colour = sep_palette('control_mean');
        group = comparison.ctrl_type;
    else
        panel = n_cols_grid + k - n_ctrl;
        colour = sep_palette('experimental_mean');
        group = comparison.exp_type;
    end
    subplot(2, n_cols_grid, panel);
    hold on;
    box on;
    grid on;

    % the cluster's planes shaded, behind everything
    fill(cluster_edges_um([1 2 2 1]), [0 0 y_top y_top], [0.9 0.9 0.9], 'EdgeColor', ...
        'none');

    % the planes where part of the footprint has no tissue, opaque where none
    % has: beside a gap the asymmetry can come from the tissue's edge
    for i = find(profiles.coverage(:, k) < 1)'
        fill(ap_um(i) + [-5 5 5 -5], [0 0 y_top y_top], sep_palette('no_data'), ...
            'EdgeColor', 'none', 'FaceAlpha', 1 - profiles.coverage(i, k));
    end

    % the sections' edges at the bottom, a third of the height at the largest
    change = plane_changes(:, k) / max(plane_changes(:, k)) * y_top / 3;
    plot(edge_um, change, '-', 'Color', [0.35 0.35 0.35], 'LineWidth', 0.8);

    % the control mice's range, then the mouse
    plot(ap_um, ctrl_low, '-', 'Color', sep_palette('paired_lines'), 'LineWidth', 0.8);
    plot(ap_um, ctrl_high, '-', 'Color', sep_palette('paired_lines'), 'LineWidth', 0.8);
    plot(ap_um, profiles.auto(:, k), '--', 'Color', ...
        sep_palette('autofluorescence_dots'), 'LineWidth', 1);
    plot(ap_um, profiles.raw(:, k), '-', 'Color', colour, 'LineWidth', 1.6);
    xlim(ap_um([1 end]));
    ylim([0 y_top]);
    set(gca, 'FontSize', 8);
    title(sprintf('%s (%s)', names{k}, group), 'FontSize', 9, 'Color', colour);
    xlabel('AP from the cluster''s centre (um), posterior >', 'FontSize', 8);
    if panel == 1 || panel == n_cols_grid + 1
        ylabel('AI over the cluster''s footprint', 'FontSize', 8);
    end
end

% the title: what each line is
title_line = sprintf(['Each mouse''s asymmetry along AP through the %s cluster''s ' ...
                      'coronal footprint - %s'], comparison.cluster_region, ...
                      strrep(file_tag, '_', ' '));
lines_line = sprintf(['Thick: the raw AI (|L - R| / (L + R) per plane); dashed: the ' ...
                      'autofluorescence AI; grey: the lowest and highest %s mouse''s ' ...
                      'raw AI; shaded: the cluster''s planes; darker grey: planes where ' ...
                      'part of the footprint has no tissue, opaque where none has.'], ...
                      comparison.ctrl_type);
section_line = ['Bottom: the collected stack''s change from plane to plane in the ' ...
                'footprint, unsmoothed; it peaks between two sections, which the ' ...
                'registered volume holds each over the planes nearest it.'];
sgtitle({title_line, ['\rm\fontsize{10}' lines_line], ['\rm\fontsize{10}' section_line]}, ...
    'FontSize', 13, 'FontWeight', 'bold');
saveas(fig, fullfile(comp_out_dir, ['Influence_Profiles_' file_tag '.fig']));
exportgraphics(fig, fullfile(comp_out_dir, ['Influence_Profiles_' file_tag '.png']), ...
    'Resolution', 300);
end

function names = short_names(names)
% Mice's names without the line's suffix ('CGF027_Gria1' as 'CGF027').

names = cellfun(@(n) strtok(n, '_'), names, 'UniformOutput', false);
end
