close all
clear all
clc

% /// Rebuilds volume_for_ordering.tiff with a chosen channel-to-colour mapping ///
%
% The ordering volume is the RGB composite that SliceOrderEditor displays. It is
% written by generateSliceVolume, which simply drops the first three channels
% into R, G, B in whatever order they happen to be stored. Because the
% registration channel is promoted to position 1, that puts DAPI in red and
% autofluorescence in blue, which reads backwards against the usual convention
% of DAPI in blue.
%
% Nothing in the analysis depends on this: every downstream script opens the
% channels by filename (chan01_DAPI, chan02_Cy5, chan03_Cy3). The composite is
% used only by the ordering GUI. So the colours are free to change, and doing so
% cannot affect any result.
%
% This script re-does exactly what generateSliceVolume does to build the
% composite (same resize, same background normalisation, same 99th percentile
% scaling), only with the colour assignment under your control. It reads the
% already-extracted volume_centered channels, so P1 does not need rerunning.
%
% IMPORTANT: slice order is untouched. Page N of the rebuilt volume is the same
% section as page N of the old one, so any curation you have already done stays
% valid.
%
% CAVEAT: rerunning P1 calls generateSliceVolume again, which will overwrite the
% composite with the original mapping. Rerun this script afterwards if that
% happens.

%% User-defined parameters

% Where the project lives. Derived from the location of the code rather than
% written out, so the tree can be moved or copied to another drive as is.
paths = get_paths();

% Cohort selection (mice come from the shared registry get_cohort.m).
% Set mice_to_process to {} to rebuild every mouse in groups_to_process.
groups_to_process = {'young'};
mice_to_process   = {};

% Which channel goes into which display colour.
% Use role names: 'dapi' | 'nano' | 'auto' | 'egfp' | 'none'
red_channel   = 'auto';    % autofluorescence
green_channel = 'nano';    % SEP-GluA1, the signal of interest
blue_channel  = 'dapi';    % nuclei, conventionally blue

%% Add paths

lightsuiteDir = paths.lightsuite;
addpath(genpath(lightsuiteDir))

%% Resolve cohort

get_cohort('verify');
if isempty(mice_to_process)
    cohort = get_cohort('groups', groups_to_process);
else
    cohort = get_cohort('names', mice_to_process);
end
fprintf('Rebuilding ordering volume for %d mouse/mice.\n', numel(cohort));
fprintf('Mapping: R = %s, G = %s, B = %s\n\n', red_channel, green_channel, blue_channel);

wanted_roles = {red_channel, green_channel, blue_channel};

%% Loop over mice

for mouse_idx = 1:numel(cohort)

    mousename = cohort(mouse_idx).name;
    procpath  = fullfile(cohort(mouse_idx).base_dir, 'lightsuite');
    vc_dir    = fullfile(procpath, 'volume_centered');
    out_file  = fullfile(procpath, 'volume_for_ordering.tiff');

    fprintf('=== %s ===\n', mousename);

    sliceinfo_file = fullfile(procpath, 'sliceinfo.mat');
    if ~exist(sliceinfo_file, 'file')
        warning('No sliceinfo.mat for %s, skipping.', mousename);
        continue
    end
    S = load(sliceinfo_file);
    sliceinfo = S.sliceinfo;

    % Warn rather than silently changing a file being curated
    decisions = fullfile(procpath, 'volume_for_ordering_processing_decisions.txt');
    if exist(decisions, 'file')
        fprintf('  note: this mouse already has a decisions file. Slice order is\n');
        fprintf('        unchanged, so your curation stays valid; only colours change.\n');
    end

    % Size of the composite, same formula generateSliceVolume uses
    scale_hw  = ceil(sliceinfo.size_proc * sliceinfo.px_process / sliceinfo.px_register);
    n_slices  = sliceinfo.Nslices;
    scalesize = [scale_hw n_slices];

    volproc = zeros([scale_hw 3 n_slices], 'uint8');

    for slot = 1:3
        role = lower(wanted_roles{slot});
        if strcmp(role, 'none')
            continue
        end

        ch_idx = channel_index_for_role(role, sliceinfo.channames);
        ch_file = fullfile(vc_dir, sprintf('chan%02d_%s.tiff', ch_idx, sliceinfo.channames{ch_idx}));
        if ~exist(ch_file, 'file')
            error('Channel file not found: %s', ch_file);
        end

        % Read the whole stack for this channel
        info = imfinfo(ch_file);
        stack = zeros(info(1).Height, info(1).Width, numel(info), 'single');
        for z = 1:numel(info)
            stack(:,:,z) = single(imread(ch_file, z));
        end

        % Same processing generateSliceVolume applies to the composite
        stack = imresize3(stack, scalesize);

        backproc = median(single(sliceinfo.backvalues(ch_idx, :)));
        if backproc > 0
            stack = (stack - backproc) ./ backproc;
        end

        maxval = quantile(stack, 0.99, 'all');
        volproc(:, :, slot, :) = uint8(255 * stack ./ maxval);

        fprintf('  %s <- %s (%s)\n', upper(colour_name(slot)), role, sliceinfo.channames{ch_idx});
    end

    % Page count must match, otherwise an existing decisions file would point
    % at the wrong sections
    if exist(out_file, 'file')
        old_pages = numel(imfinfo(out_file));
        if old_pages ~= n_slices
            error(['Refusing to overwrite: existing volume has %d pages but %d ' ...
                   'were rebuilt. Slice indices would no longer line up.'], old_pages, n_slices);
        end
        delete(out_file);
    end

    options.compress = 'lzw';
    options.message  = false;
    options.color    = true;
    options.big      = false;
    saveastiff(volproc, out_file, options);

    fprintf('  wrote %s (%d slices)\n\n', out_file, n_slices);

end

fprintf('Done.\n');

%% Local function: map a role name to its index in sliceinfo.channames

function idx = channel_index_for_role(role, channames)
% Roles are named after what the channel is, not the dye it was acquired with,
% so the calling code does not have to remember that nano is Cy5 and autofluo
% is Cy3.

switch role
    case 'dapi', dye = 'DAPI';
    case 'nano', dye = 'Cy5';
    case 'auto', dye = 'Cy3';
    case 'egfp', dye = 'EGFP';
    otherwise
        error('Unknown channel role "%s" (use dapi, nano, auto, egfp or none).', role);
end

idx = find(strcmpi(channames, dye), 1);
if isempty(idx)
    error('Channel "%s" (role %s) not found in this mouse. Available: %s', ...
        dye, role, strjoin(cellstr(channames), ', '));
end
end

%% Local function: colour name for a slot, just for the printout

function name = colour_name(slot)
names = {'red', 'green', 'blue'};
name = names{slot};
end
