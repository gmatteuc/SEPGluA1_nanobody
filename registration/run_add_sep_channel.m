%% run_add_sep_channel
% ===== Carry the SEP channel into registered space =====
%
% Registration, step 2 of 2:
%   1. run_register_to_atlas   register the sections to the atlas
%   2. run_add_sep_channel     carry the SEP channel into registered space  <- this script
%
% run_register_to_atlas registers five channels (DAPI, NANO, AUTO, DIFF, MASK).
% The microscope recorded a sixth, on the green filter it names EGFP, and that one
% stops at volume_centered. Every mouse is a SEP-GluA1 knock-in, so the green
% channel was meant as the tagged receptor itself, the whole GluA1 pool (fixed,
% permeabilised tissue loses the pH sensitivity that makes SEP surface-specific in
% a living cell), and nano/SEP as surface receptor per receptor expressed.
% Measured directly (run_sep_channel_check), the green channel tracks the
% autofluorescence instead (rho 0.79 across the ten adults), so nano/SEP is not a
% surface fraction. The SEP channel is added as an extra registered volume, and
% nothing is overwritten:
%
%   in   <mouse>\lightsuite\volume_centered\chan04_EGFP.tiff   raw, per slice
%                                                              (EGFP names the
%                                                              filter, SEP the
%                                                              molecule)
%   out  <mouse>\lightsuite\volume_registered_sep\
%          chan01_DAPI.tiff   the check channel, see below
%          chan02_SEP.tiff    what the analysis reads
%
% volume_aligned and volume_registered are only read, so every earlier result
% stands, and the reference channel is a switch in the analysis rather than a
% fork of the data. LightSuite is not changed either: the registration is not
% redone, it is re-applied. The registration has two stages:
%
%   centered --( per-slice 2D transform, tformslices )--> aligned
%   aligned  --( per-slice elastix B-spline + affine, then one 3D rigid )--> registered
%
% The second is saved in full (transform_params.mat and the elastix_* folders),
% so it costs a transformix pass and no refitting. The first is not saved:
% alignSliceVolume fits tformslices to point clouds and keeps only the volume it
% made. Refitting would be wrong: the fit subsamples the point clouds at random,
% so a rerun is never the same, and a reference channel a fraction of a pixel off
% the channel it normalises is the artefact to avoid. Instead each slice's
% transform is recovered from the two images that fit produced, the centered and
% the aligned DAPI slice: phase correlation gives a first guess on a coarse grid,
% an affine refinement per level brings it to full resolution, and the result is
% scored against the image it had to reproduce. On MG897 every slice came back at
% r = 1.00000 over tissue with no residual shift, so the recovered transform is
% the original one. The SEP channel is warped with it and goes through the saved
% elastix transforms beside the DAPI.
%
% That DAPI is the check: it reaches registered space and is compared there with
% the DAPI run_register_to_atlas registered. If the two agree voxel for voxel, the
% SEP volume beside it is in the same space as NANO and AUTO. The check is printed
% per mouse and drawn in
% data\comparisons_v2\processing_diagnostics\sep_channel\<mouse>.png. The DAPI is
% also needed because with a single channel in the folder loadLargeSliceVolume
% squeezes the channel dimension away, and generateRegisteredSliceVolume then
% reads slices as channels.
%
% About 30-45 min per brain, an image registration per slice and a transformix
% pass per slice and channel: run it detached.
%
% Setup: every registered mouse of the young, naive and rws groups. Run
% sep_setup_paths first, once per MATLAB session; the code is in
% pipeline\add_sep_channel.m.

clear; clc; close all;

%% Settings

% project folders, worked out from where the code sits (see get_paths)
paths = get_paths();

% groups to process ('rws', 'naive', 'behavior', 'young'); a mouse without
% transform_params.mat has never been registered and is skipped with a note, so
% several groups can be given at once and the brains that are ready are picked out
groups_to_process = {'young', 'naive', 'rws'};

% mice to process, by name ({} = every mouse of groups_to_process)
mice_to_process = {};

% P4BIS_MICE, a comma-separated list in the environment, replaces the list above
% when it is set, so the cohort can be split across MATLAB sessions run side by
% side (nothing here needs the GPU or a shared file); an environment variable,
% since the clear at the top of this file empties the workspace:
%   $env:P4BIS_MICE = 'MG903_SepGluA_P20,MG913_SepGluA_P20'
%   matlab -batch "cd('D:\sep_histology\code'); sep_setup_paths; run_add_sep_channel"
env_mice = getenv('P4BIS_MICE');
if ~isempty(env_mice)
    mice_to_process = strtrim(strsplit(env_mice, ','));
end

% redo a mouse that already has volume_registered_sep
overwrite = false;

% levels of the per-slice recovery, coarsest grid first, one affine refinement per
% level (the same world extent, fewer pixels). The coarse levels take seconds and
% bring the fit within a pixel; the full-resolution level, started from there,
% needs one pyramid pass and lands on the transform exactly. Without the 1 a
% systematic half-pixel offset stays (measured, see min_slice_corr).
recovery_levels = [4 2 1];

% correlation with the aligned DAPI below which a slice's recovered transform is
% reported by number. With the levels above, every slice of MG897 came back at
% r = 1.00000 and a residual shift of 0 px, so below 0.99 something is wrong with
% the slice, not with the method.
min_slice_corr = 0.99;

%% Run

% pass the settings to the code, under the same names
run_settings = struct();
run_settings.paths = paths;
run_settings.groups_to_process = groups_to_process;
run_settings.mice_to_process = mice_to_process;
run_settings.overwrite = overwrite;
run_settings.recovery_levels = recovery_levels;
run_settings.min_slice_corr = min_slice_corr;
add_sep_channel(run_settings);
