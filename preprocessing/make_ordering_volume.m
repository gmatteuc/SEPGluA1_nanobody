%% make_ordering_volume
% ===== Rebuild the slice-ordering composite with chosen colours =====
%
% A tool, run by hand. volume_for_ordering.tiff is the colour composite that
% SliceOrderEditor shows. generateSliceVolume writes it with the first three
% channels as red, green and blue, in the order they are stored; the
% registration channel comes first, so DAPI is red and the autofluorescence
% blue, against the usual DAPI in blue. Only the ordering window reads the
% composite (every later step opens the channels by file name: chan01_DAPI,
% chan02_Cy5, chan03_Cy3), so its colours can change without changing any result.
%
% Rebuilds the composite as generateSliceVolume does (the same resize,
% background normalisation and 99th percentile scaling), with the colours chosen
% below, from the channels already in lightsuite\volume_centered\, so
% run_extract_and_center need not run again. The slice order is untouched: page
% N is the same section as before, so a curation already done stays valid, and
% the run stops if the page count differs. A rerun of run_extract_and_center
% writes the original colours back; run this again after it.
%
% Setup: the whole young cohort, nano in green. Run sep_setup_paths first, once
% per MATLAB session.

clear; clc; close all;

%% Settings

% mice to rebuild, from the cohort registry (get_cohort): the groups ('rws',
% 'naive', 'behavior', 'young'), or the mice named, which take precedence ({} =
% the groups)
groups_to_process = {'young'};
mice_to_process = {};

% the channel shown in each colour, by role: 'dapi', 'nano', 'auto', 'egfp' or
% 'none'; the autofluorescence in red, the SEP-GluA1 signal in green, the nuclei
% in blue, as usual
red_channel = 'auto';
green_channel = 'nano';
blue_channel = 'dapi';

%% Resolve cohort

% check that the registry still lists the adults in their legacy order
get_cohort('verify');

% the mice named, or else every mouse of the groups
if isempty(mice_to_process)
    cohort = get_cohort('groups', groups_to_process);
else
    cohort = get_cohort('names', mice_to_process);
end
fprintf('Rebuilding ordering volume for %d mouse/mice.\n', numel(cohort));
fprintf('Mapping: R = %s, G = %s, B = %s\n\n', red_channel, green_channel, blue_channel);

wanted_roles = {red_channel, green_channel, blue_channel};

%% Rebuild each mouse's composite

for mouse_idx = 1:numel(cohort)

    % the mouse, its centred channels and the composite
    mousename = cohort(mouse_idx).name;
    procpath = fullfile(cohort(mouse_idx).base_dir, 'lightsuite');
    vc_dir = fullfile(procpath, 'volume_centered');
    out_file = fullfile(procpath, 'volume_for_ordering.tiff');

    fprintf('=== %s ===\n', mousename);

    % its extraction record; a mouse without one is skipped
    sliceinfo_file = fullfile(procpath, 'sliceinfo.mat');
    if ~exist(sliceinfo_file, 'file')
        warning('No sliceinfo.mat for %s, skipping.', mousename);
        continue
    end
    S = load(sliceinfo_file);
    sliceinfo = S.sliceinfo;

    % say so when the composite is being curated
    decisions = fullfile(procpath, 'volume_for_ordering_processing_decisions.txt');
    if exist(decisions, 'file')
        fprintf('  note: this mouse already has a decisions file. Slice order is\n');
        fprintf('        unchanged, so your curation stays valid; only colours change.\n');
    end

    % the size of the composite, by generateSliceVolume's formula
    scale_hw = ceil(sliceinfo.size_proc * sliceinfo.px_process / sliceinfo.px_register);
    n_slices = sliceinfo.Nslices;
    scalesize = [scale_hw n_slices];

    volproc = zeros([scale_hw 3 n_slices], 'uint8');

    for slot = 1:3
        role = lower(wanted_roles{slot});
        if strcmp(role, 'none')
            continue
        end

        % the channel's centred stack
        ch_idx = channel_index_for_role(role, sliceinfo.channames);
        ch_file = fullfile(vc_dir, sprintf('chan%02d_%s.tiff', ch_idx, ...
            sliceinfo.channames{ch_idx}));
        if ~exist(ch_file, 'file')
            error('Channel file not found: %s', ch_file);
        end

        % read the whole stack
        info = imfinfo(ch_file);
        stack = zeros(info(1).Height, info(1).Width, numel(info), 'single');
        for z = 1:numel(info)
            stack(:, :, z) = single(imread(ch_file, z));
        end

        % resize, normalise to the background and scale to the 99th percentile,
        % as generateSliceVolume does
        stack = imresize3(stack, scalesize);

        backproc = median(single(sliceinfo.backvalues(ch_idx, :)));
        if backproc > 0
            stack = (stack - backproc) ./ backproc;
        end

        maxval = quantile(stack, 0.99, 'all');
        volproc(:, :, slot, :) = uint8(255 * stack ./ maxval);

        fprintf('  %s <- %s (%s)\n', upper(colour_name(slot)), role, ...
            sliceinfo.channames{ch_idx});
    end

    % the page count must match, or a decisions file would point at the wrong
    % sections
    if exist(out_file, 'file')
        old_pages = numel(imfinfo(out_file));
        if old_pages ~= n_slices
            error(['Refusing to overwrite: existing volume has %d pages but %d ' ...
                   'were rebuilt. Slice indices would no longer line up.'], ...
                   old_pages, n_slices);
        end
        delete(out_file);
    end

    % write it as generateSliceVolume does
    options.compress = 'lzw';
    options.message = false;
    options.color = true;
    options.big = false;
    saveastiff(volproc, out_file, options);

    fprintf('  wrote %s (%d slices)\n\n', out_file, n_slices);

end

fprintf('Done.\n');

% ===== Local functions =====

function idx = channel_index_for_role(role, channames)
% The index in channames of the channel with a role; roles name what the channel
% is, not its dye, so nobody has to remember that the nano is Cy5.

switch role
    case 'dapi'
        dye = 'DAPI';
    case 'nano'
        dye = 'Cy5';
    case 'auto'
        dye = 'Cy3';
    case 'egfp'
        dye = 'EGFP';
    otherwise
        error('Unknown channel role "%s" (use dapi, nano, auto, egfp or none).', role);
end

idx = find(strcmpi(channames, dye), 1);
if isempty(idx)
    error('Channel "%s" (role %s) not found in this mouse. Available: %s', ...
        dye, role, strjoin(cellstr(channames), ', '));
end
end

function name = colour_name(slot)
% The colour of a slot, for the printout.

names = {'red', 'green', 'blue'};
name = names{slot};
end
