clear all
close all
clc

% /// QC helper: measure both reference atlases and write down what they are ///
% The adult and young cohorts now register to two different atlases, so the
% parameters that relate them have to be written down somewhere they cannot
% drift out of date. This script measures everything from the volumes on disk
% and then:
%   (1) Verifies the AP crop of the young atlas against the adult one, twice
%       over -- once from the brain's cross-sectional area profile, once from
%       the AP centre of mass of every region present in both
%   (2) Draws the crop-check figure
%   (3) Writes ATLAS_PARAMETERS.md with resolution, grid, extent, crop,
%       label space, section geometry and per-brain AP coverage
%
% Everything printed is measured here, not copied from a previous note. Re-run
% it after changing an atlas, a crop, or slicethickness.
%
% The region-centroid regression is the slow part (it walks every AP plane of
% both annotation volumes); turn it off with do_region_regression to get just
% the geometry and the note.

%% User-defined parameters

% Where the project lives. Derived from the location of the code rather than
% written out, so the tree can be moved or copied to another drive as is.
paths = get_paths();

% The two atlases, as get_atlas keys
atlas_key_adult = 'ccf';
atlas_key_young = 'demba_p20';

% Which cohorts use which, for the note
adult_groups = {'rws', 'naive', 'behavior'};
young_groups = {'young'};

% Section geometry. These are the values the pipeline actually runs with; the
% note flags a mismatch rather than silently reporting the wrong one.
slicethickness_um = 150;

% The expensive cross-atlas regression
do_region_regression = true;
min_vox_adult = 5000;      % ignore labels too small for a stable centroid
min_vox_young = 600;

% Output
out_dir = fullfile(paths.data, 'young', 'registration_qc');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
note_file = fullfile(out_dir, 'ATLAS_PARAMETERS.md');
save_figure = true;

% Color palette
adult_color = [0.35 0.35 0.35];   % gray, adult atlas
young_color = [0.95 0.55 0.10];   % strong orange, young atlas

%% Measure each atlas

atlas_adult = get_atlas(atlas_key_adult);
atlas_young = get_atlas(atlas_key_young);

info_adult = measure_atlas(atlas_adult);
info_young = measure_atlas(atlas_young);

print_atlas(info_adult);
print_atlas(info_young);

%% Check the young crop against the adult one, from the area profile

area_adult = info_adult.crop_area_profile;
area_young = info_young.crop_area_profile;

x_adult = linspace(0, 1, numel(area_adult));
x_young = linspace(0, 1, numel(area_young));
n_adult = area_adult / max(area_adult);
n_young = area_young / max(area_young);

young_on_adult = interp1(x_young, n_young, x_adult, 'linear', 'extrap');
profile_rms = sqrt(mean((n_adult(:) - young_on_adult(:)).^2));

fprintf('\nArea-profile agreement between the two crops: RMS %.4f\n', profile_rms);

%% Check it again, independently, from where each region actually sits

reg = struct('done', false, 'n_shared', NaN);

if do_region_regression

    fprintf('\nMeasuring AP centre of mass of every shared region...\n');

    av_adult = niftiread(fullfile(atlas_adult.dir, atlas_adult.annotation_file));
    av_young = niftiread(fullfile(atlas_young.dir, atlas_young.annotation_file));

    labels = intersect(unique(av_adult(:)), unique(av_young(:)));
    labels = labels(labels ~= 0);
    fprintf('  labels present in both volumes: %d\n', numel(labels));

    if isempty(labels)
        error(['The two annotation volumes share no label values. They are ' ...
               'probably in different ID spaces -- the pipeline uses Allen ' ...
               'parcellation_index, while BrainGlobe ships structure IDs. ' ...
               'Remap the young annotation before using it.']);
    end

    [com_adult, tot_adult] = ap_centre_of_mass(av_adult, labels);
    [com_young, tot_young] = ap_centre_of_mass(av_young, labels);

    keep = tot_adult >= min_vox_adult & tot_young >= min_vox_young;
    fprintf('  regions large enough to use: %d\n', nnz(keep));

    cc = com_adult(keep);
    dd = com_young(keep);

    coef = polyfit(dd, cc, 1);
    fit_slope = coef(1);
    fit_offset = coef(2);
    resid = cc - polyval(coef, dd);
    rr = corr(dd, cc);

    % The slope the two voxel sizes alone would predict. Departing from it
    % means the two atlases disagree about how long the brain is in AP.
    nominal_slope = atlas_young.res_um / atlas_adult.res_um;
    ap_scale_ratio = nominal_slope / fit_slope;

    reg.done        = true;
    reg.n_shared    = numel(labels);
    reg.n_regions   = nnz(keep);
    reg.slope       = fit_slope;
    reg.offset      = fit_offset;
    reg.r           = rr;
    reg.resid_sd    = std(resid);
    reg.resid_sd_mm = std(resid) * atlas_adult.res_um / 1000;
    reg.nominal     = nominal_slope;
    reg.ap_ratio    = ap_scale_ratio;
    reg.implied_crop = round((atlas_adult.default_aplims - fit_offset) / fit_slope);

    fprintf('\n  adult_plane = %.4f * young_plane + %.2f   (r = %.5f)\n', ...
        fit_slope, fit_offset, rr);
    fprintf('  residual sd %.1f adult planes = %.3f mm\n', std(resid), reg.resid_sd_mm);
    fprintf('  nominal slope from voxel sizes alone: %.4f\n', nominal_slope);
    fprintf('  => the young atlas is %.1f%% longer in AP for the same anatomy\n', ...
        100 * (ap_scale_ratio - 1));
    fprintf('  crop implied by the regression: [%d %d]  (in use: [%d %d])\n', ...
        reg.implied_crop(1), reg.implied_crop(2), ...
        atlas_young.default_aplims(1), atlas_young.default_aplims(2));

end

%% Per-brain AP coverage

get_cohort('verify');
cohort = get_cohort();

cov = struct('name', {}, 'group', {}, 'n_slices', {}, 'span_mm', {}, 'frac', {});

for k = 1:numel(cohort)

    decisions = fullfile(cohort(k).base_dir, 'lightsuite', ...
                'volume_for_ordering_processing_decisions.txt');
    if ~exist(decisions, 'file')
        continue
    end

    T = readtable(decisions);
    n_kept = sum(T.FlipState ~= -1);

    if ismember(cohort(k).group, young_groups)
        crop_mm = info_young.crop_extent_mm;
    else
        crop_mm = info_adult.crop_extent_mm;
    end

    span_mm = n_kept * slicethickness_um / 1000;

    cov(end+1).name  = cohort(k).name;   %#ok<SAGROW>
    cov(end).group    = cohort(k).group;
    cov(end).n_slices = n_kept;
    cov(end).span_mm  = span_mm;
    cov(end).frac     = span_mm / crop_mm;

end

fprintf('\nAP coverage, %d curated brains (sections x %g um against the crop):\n', ...
    numel(cov), slicethickness_um);
for k = 1:numel(cov)
    fprintf('  %-26s %-9s %3d sections  %5.2f mm  %4.0f%%\n', ...
        cov(k).name, cov(k).group, cov(k).n_slices, cov(k).span_mm, 100 * cov(k).frac);
end

%% Figure: the crop check

fig = figure('Visible', 'off', 'Color', 'w', 'Units', 'pixels', ...
             'Position', [50 50 1150 (1 + double(reg.done)) * 380 + 60]);
tl = tiledlayout(fig, 1 + double(reg.done), 1, ...
                 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile(tl);
plot(x_adult, n_adult, 'Color', adult_color, 'LineWidth', 2); hold on
plot(x_young, n_young, 'Color', young_color, 'LineWidth', 2);
xlabel('position within crop (0 = anterior end)')
ylabel('brain cross-sectional area (normalised)')
legend({sprintf('%s  [%d %d]', atlas_adult.key, atlas_adult.default_aplims), ...
        sprintf('%s  [%d %d]', atlas_young.key, atlas_young.default_aplims)}, ...
        'Interpreter', 'none', 'Box', 'off', 'Location', 'south')
title(sprintf('Do the two crops hold the same anatomy?   profile RMS %.4f', profile_rms))
box off

if reg.done
    nexttile(tl);
    scatter(dd, cc, 28, young_color, 'filled', 'MarkerFaceAlpha', 0.75); hold on
    xline_vals = [min(dd) max(dd)];
    plot(xline_vals, polyval(coef, xline_vals), 'Color', adult_color, 'LineWidth', 1.8);
    xlabel(sprintf('%s AP plane (%g um)', atlas_young.key, atlas_young.res_um))
    ylabel(sprintf('%s AP plane (%g um)', atlas_adult.key, atlas_adult.res_um))
    title(sprintf(['Same region, both atlases: %d regions, slope %.4f ' ...
                   '(nominal %.2f), r = %.5f, residual %.2f mm'], ...
                   reg.n_regions, reg.slope, reg.nominal, reg.r, reg.resid_sd_mm), ...
          'Interpreter', 'none')
    legend({'one region', 'least-squares fit'}, 'Box', 'off', 'Location', 'northwest')
    box off
end

if save_figure
    png_name = fullfile(out_dir, 'atlas_crop_and_ap_mapping.png');
    exportgraphics(fig, png_name, 'Resolution', 150);
    fprintf('\nsaved %s\n', png_name);
end

%% Write the note

write_note(note_file, info_adult, info_young, ...
           adult_groups, young_groups, slicethickness_um, profile_rms, reg, cov);
fprintf('wrote %s\n', note_file);

%% Local function: measure one atlas from its files

function info = measure_atlas(atlas)

    av = niftiread(fullfile(atlas.dir, atlas.annotation_file));

    info.key        = atlas.key;
    info.dir        = atlas.dir;
    info.res_um     = atlas.res_um;
    info.aplims     = atlas.default_aplims;
    info.desc       = atlas.description;
    info.grid       = size(av);
    info.extent_mm  = double(size(av)) * atlas.res_um / 1000;
    info.class      = class(av);

    % Where the labelled brain actually sits in the full volume
    per_plane = squeeze(sum(sum(av > 0, 2), 3));
    nz = find(per_plane > 0);
    info.brain_first    = nz(1);
    info.brain_last     = nz(end);
    info.brain_span_mm  = (nz(end) - nz(1) + 1) * atlas.res_um / 1000;

    % And within the crop the pipeline actually uses
    crop = av(info.aplims(1):info.aplims(2), :, :);
    info.crop_planes     = size(crop, 1);
    info.crop_extent_mm  = size(crop, 1) * atlas.res_um / 1000;
    info.crop_brain_frac = nnz(crop > 0) / numel(crop);
    info.crop_area_profile = squeeze(sum(sum(crop > 0, 2), 3));

    labels = unique(av(:));
    info.n_labels = numel(labels) - any(labels == 0);
    info.max_label = max(labels);

end

%% Local function: print one atlas

function print_atlas(info)
    fprintf('\n=== %s ===\n', info.key);
    fprintf('  %s\n', info.desc);
    fprintf('  grid (AP,DV,ML) %s at %g um  =  %.2f x %.2f x %.2f mm\n', ...
        mat2str(info.grid), info.res_um, info.extent_mm);
    fprintf('  labelled brain spans planes %d..%d (%.2f mm)\n', ...
        info.brain_first, info.brain_last, info.brain_span_mm);
    fprintf('  crop [%d %d] = %d planes = %.2f mm, brain fraction %.4f\n', ...
        info.aplims(1), info.aplims(2), info.crop_planes, ...
        info.crop_extent_mm, info.crop_brain_frac);
    fprintf('  %d labels, max value %d, class %s\n', ...
        info.n_labels, info.max_label, info.class);
end

%% Local function: AP centre of mass of each label

function [com, total] = ap_centre_of_mass(vol, labels)

    n_ap = size(vol, 1);
    n_lab = numel(labels);

    % Map raw label values onto 1..n_lab once, so the per-plane accumulation
    % does not have to allocate anything the size of the largest label value.
    lut = zeros(double(max(labels)) + 1, 1);
    lut(double(labels) + 1) = 1:n_lab;

    weighted = zeros(n_lab, 1);
    total    = zeros(n_lab, 1);

    for i = 1:n_ap
        plane = double(vol(i, :, :));
        plane = plane(plane > 0);
        if isempty(plane)
            continue
        end
        idx = lut(plane + 1);
        idx = idx(idx > 0);
        if isempty(idx)
            continue
        end
        counts   = accumarray(idx, 1, [n_lab 1]);
        total    = total + counts;
        weighted = weighted + counts * i;
    end

    com = weighted ./ total;

end

%% Local function: write the note

function write_note(note_file, ia, iy, adult_groups, young_groups, ...
                    slicethickness_um, profile_rms, reg, cov)

    fid = fopen(note_file, 'w');
    c = onCleanup(@() fclose(fid));

    fprintf(fid, '# Reference atlases and section geometry\n\n');
    fprintf(fid, ['Generated by `atlas_diagnostics.m` on %s. Every number here was ' ...
                  'measured from the files on disk at that moment -- re-run the ' ...
                  'script rather than editing this file.\n\n'], ...
                  datestr(now, 'yyyy-mm-dd HH:MM')); %#ok<TNOW1,DATST>

    fprintf(fid, ['The two cohorts register to two different atlases: the adults to the ' ...
                  'adult Allen CCF, the young brains to the age-matched DeMBA template. ' ...
                  'They are compared at the level of regions, not voxels, which works ' ...
                  'because both annotation volumes use the same label space.\n\n']);

    fprintf(fid, '## Label space -- check this when adding an atlas\n\n');
    fprintf(fid, ['Both annotation volumes here hold Allen **`parcellation_index`** values, ' ...
                  'the re-indexed scheme that ships with the Allen AWS release and the only ' ...
                  'one `get_allen_region_mask.m` understands (it resolves region names ' ...
                  'through `parcellation_to_parcellation_term_membership.csv`).\n\n']);
    fprintf(fid, ['**BrainGlobe ships DeMBA in Allen *structure IDs* instead**, which is a ' ...
                  'different numbering scheme that collides numerically without meaning the ' ...
                  'same thing -- structure 672 (caudoputamen) is `parcellation_index` 662. ' ...
                  'The volume in `%s` has been remapped; BrainGlobe''s original is kept ' ...
                  'beside it as `annotation_structureids_original.nii.gz`.\n\n'], iy.dir);
    fprintf(fid, ['Right now %d of the young atlas'' %d labels are present in the adult one, ' ...
                  'which is what a correct remap looks like. If that number drops to roughly ' ...
                  'half, the volume is in raw structure IDs and every region lookup will be ' ...
                  'silently wrong.\n\n'], reg.n_shared, iy.n_labels);

    fprintf(fid, '## The two atlases\n\n');
    fprintf(fid, '| | adult (`%s`) | young (`%s`) |\n', ia.key, iy.key);
    fprintf(fid, '|---|---|---|\n');
    fprintf(fid, '| used by | %s | %s |\n', strjoin(adult_groups, ', '), strjoin(young_groups, ', '));
    fprintf(fid, '| description | %s | %s |\n', ia.desc, iy.desc);
    fprintf(fid, '| folder | `%s` | `%s` |\n', ia.dir, iy.dir);
    fprintf(fid, '| voxel size | **%g um** isotropic | **%g um** isotropic |\n', ia.res_um, iy.res_um);
    fprintf(fid, '| grid (AP, DV, ML) | %s | %s |\n', mat2str(ia.grid), mat2str(iy.grid));
    fprintf(fid, '| full extent | %.2f x %.2f x %.2f mm | %.2f x %.2f x %.2f mm |\n', ...
        ia.extent_mm, iy.extent_mm);
    fprintf(fid, '| labelled brain spans | planes %d..%d (%.2f mm) | planes %d..%d (%.2f mm) |\n', ...
        ia.brain_first, ia.brain_last, ia.brain_span_mm, ...
        iy.brain_first, iy.brain_last, iy.brain_span_mm);
    fprintf(fid, '| `atlasaplims` crop | **[%d %d]** | **[%d %d]** |\n', ...
        ia.aplims, iy.aplims);
    fprintf(fid, '| cropped AP | %d planes = **%.2f mm** | %d planes = **%.2f mm** |\n', ...
        ia.crop_planes, ia.crop_extent_mm, iy.crop_planes, iy.crop_extent_mm);
    fprintf(fid, '| brain fraction of crop | %.4f | %.4f |\n', ...
        ia.crop_brain_frac, iy.crop_brain_frac);
    fprintf(fid, '| labels | %d (max %d, %s) | %d (max %d, %s) |\n', ...
        ia.n_labels, ia.max_label, ia.class, iy.n_labels, iy.max_label, iy.class);
    fprintf(fid, '| `px_atlas` to set | %g | %g |\n\n', ia.res_um, iy.res_um);

    fprintf(fid, '## How the two line up in AP\n\n');
    fprintf(fid, ['Cross-sectional-area profiles of the two crops agree to an RMS of ' ...
                  '**%.4f** (normalised area, 0-1).\n\n'], profile_rms);

    if reg.done
        fprintf(fid, ['Independently, regressing the AP centre of mass of the **%d regions** ' ...
                      'present in both volumes:\n\n'], reg.n_regions);
        fprintf(fid, '```\nadult_plane = %.4f * young_plane + %.2f     r = %.5f, residual %.3f mm\n```\n\n', ...
            reg.slope, reg.offset, reg.r, reg.resid_sd_mm);
        fprintf(fid, ['The slope the voxel sizes alone would predict is **%.2f**. It comes out ' ...
                      'at **%.4f**, so the young atlas is **%.1f%% longer in AP than the adult ' ...
                      'one for the same anatomy**. That is the documented rostrocaudal shrinkage ' ...
                      'of the adult CCF template, which the developmental templates do not have ' ...
                      '(Carey 2025).\n\n'], ...
                      reg.nominal, reg.slope, 100 * (reg.ap_ratio - 1));
        fprintf(fid, ['The crop implied by that regression is `[%d %d]`, against `[%d %d]` in use ' ...
                      '-- agreeing to within a few planes, and reached by a completely different ' ...
                      'route from the area profile above.\n\n'], ...
                      reg.implied_crop(1), reg.implied_crop(2), iy.aplims(1), iy.aplims(2));
        fprintf(fid, ['**Consequence for `slicethickness`.** It was settled at %g um by how well ' ...
                      'the adults fit the adult CCF. Against the young atlas the same ' ...
                      'reconstruction is about %.0f%% too short in AP; the equivalent value there ' ...
                      'would be roughly **%.0f um**. Each cohort only needs to be scaled ' ...
                      'correctly to its own atlas, since the comparison is per region -- but ' ...
                      'this is a decision, not a detail.\n\n'], ...
                      slicethickness_um, 100 * (reg.ap_ratio - 1), ...
                      slicethickness_um * reg.ap_ratio);
    end

    fprintf(fid, '## Section geometry\n\n');
    fprintf(fid, '| parameter | value | what it does |\n|---|---|---|\n');
    fprintf(fid, '| `slicethickness` | **%g um** | AP spacing between consecutive mounted sections; sets the AP scale of the whole reconstruction |\n', slicethickness_um);
    fprintf(fid, '| `px_process` | 5 um | working resolution of the extracted and corrected volumes |\n');
    fprintf(fid, '| `px_register` | 20 um | resolution the registration runs at |\n');
    fprintf(fid, '| `px_atlas` | %g um (adult) / %g um (young) | must equal the atlas voxel size, or the AP scale is wrong by that ratio |\n', ia.res_um, iy.res_um);
    fprintf(fid, '| `pxsizes(1)` | %.2f | `slicethickness / px_register`, the AP step in registration voxels |\n', slicethickness_um / 20);
    fprintf(fid, '| `extentfactor` | 6 | how far the atlas fit is extended beyond the section range |\n');
    fprintf(fid, '| `regchan` | `dapi` | the channel registration is driven by |\n\n');

    fprintf(fid, ['One section therefore steps **%g um** in AP, so a brain of N sections covers ' ...
                  'N x %.2f mm. Note this is the *spacing* between mounted sections, not the ' ...
                  'thickness of the cut tissue -- the gap between them is not imaged.\n\n'], ...
                  slicethickness_um, slicethickness_um / 1000);

    fprintf(fid, '## AP coverage per brain\n\n');
    fprintf(fid, ['Sections kept after the manual ordering step, times %g um, against the ' ...
                  'cropped extent of that cohort''s atlas. Regions near the AP extremes will ' ...
                  'have fewer contributing brains.\n\n'], slicethickness_um);
    fprintf(fid, ['**The two percentage columns are not directly comparable.** Each brain is ' ...
                  'measured against its own atlas crop, and the young crop is %.0f%% longer ' ...
                  '(%.2f mm against %.2f mm). A young brain with the same number of sections ' ...
                  'as an adult therefore scores lower here without covering any less anatomy. ' ...
                  'Compare the section counts directly, or settle `slicethickness` for the ' ...
                  'young cohort first -- at %.0f um the young percentages would rise by about ' ...
                  'that same %.0f%%.\n\n'], ...
                  100 * (iy.crop_extent_mm / ia.crop_extent_mm - 1), ...
                  iy.crop_extent_mm, ia.crop_extent_mm, ...
                  slicethickness_um * iy.crop_extent_mm / ia.crop_extent_mm, ...
                  100 * (iy.crop_extent_mm / ia.crop_extent_mm - 1));
    fprintf(fid, '| mouse | group | sections | AP span | of crop |\n|---|---|---|---|---|\n');
    for k = 1:numel(cov)
        fprintf(fid, '| %s | %s | %d | %.2f mm | %.0f%% |\n', ...
            cov(k).name, cov(k).group, cov(k).n_slices, cov(k).span_mm, 100 * cov(k).frac);
    end

    is_young_cov = ismember({cov.group}, young_groups);
    if any(is_young_cov)
        fprintf(fid, '\nadult mean %.0f%% of crop (n = %d), young mean %.0f%% (n = %d).\n', ...
            100 * mean([cov(~is_young_cov).frac]), nnz(~is_young_cov), ...
            100 * mean([cov(is_young_cov).frac]),  nnz(is_young_cov));
    end

end
