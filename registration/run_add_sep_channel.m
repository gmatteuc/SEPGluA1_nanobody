clear all
close all
clc

% /// Registration, step 2 of 2: carry the SEP channel into registered space ///
%
% run_register_to_atlas registers five channels -- DAPI, NANO, AUTO, DIFF, MASK.
% The microscope recorded a sixth, on the green filter it names EGFP, and that
% one stops at volume_centered. It is not a second label: every mouse in this
% project is a SEP-GluA1 knock-in, so the green channel is the tagged receptor
% itself. Ex vivo it reports the whole GluA1 pool -- the pH sensitivity that
% makes SEP surface-specific in a living cell is gone in fixed, permeabilised
% tissue -- while the nanobody stain reports the receptors that sat on the
% membrane. Normalising nano by it therefore asks a well-posed question, surface
% per unit receptor expressed, where nano per unit autofluorescence only asks
% surface per unit tissue. Sami expects to trust it more than the
% autofluorescence, and the only honest way to choose is to carry both to the
% end on every brain and look. So this script adds SEP as an EXTRA registered
% volume and overwrites nothing:
%
%   in   <mouse>\lightsuite\volume_centered\chan04_EGFP.tiff   raw, per slice
%                                                              (EGFP names the
%                                                              filter, SEP the
%                                                              molecule)
%   out  <mouse>\lightsuite\volume_registered_sep\
%          chan01_DAPI.tiff   the check channel, see below
%          chan02_SEP.tiff    what the analysis reads
%
% volume_aligned and volume_registered are opened read-only and stay exactly as
% they are, so every number produced so far keeps standing and the choice of
% reference channel becomes a switch in the analysis rather than a fork of the
% data. Nothing in LightSuite is modified either: the registration is not
% redone, it is RE-APPLIED.
%
% Why it can be re-applied at all. Registration happens in two stages:
%
%   centered --( per-slice 2D transform, tformslices )--> aligned
%   aligned  --( per-slice elastix B-spline + affine, then one 3D rigid )--> registered
%
% The second stage is saved in full (transform_params.mat and the elastix_*
% folders), so it costs a transformix pass and no refitting at all. The first
% stage is not saved: alignSliceVolume computes tformslices from the point
% clouds and keeps only the volume it produced. Refitting it would be the wrong
% answer -- the fit subsamples the point clouds at random, so a rerun is never
% bit-for-bit the same, and a reference channel sitting a fraction of a pixel
% off the channel it normalises is exactly the artefact one must not introduce.
% Instead each slice's transform is RECOVERED from the pair of images that fit
% already produced: the centered DAPI slice and the aligned DAPI slice are the
% same picture before and after tformslices(islice). Phase correlation gives a
% first guess on a coarse grid, an affine refinement per level brings it down
% to full resolution, and the result is scored against the very image it had to
% reproduce: on MG897 every slice came back at r = 1.00000 over tissue with no
% residual shift, so the recovered transform is not an approximation of the
% original one, it is the original one. The SEP channel is then warped with it
% and goes through the saved elastix transforms beside the DAPI.
%
% That DAPI is also why the output carries two channels. It rides all the way to
% registered space and is compared there with the DAPI that
% run_register_to_atlas registered months ago: if those two agree voxel for
% voxel, the SEP volume beside it sits in the same space as NANO and AUTO, and
% nano/SEP is a ratio between channels that line up. The check is printed per
% mouse and drawn in
% data\comparisons_v2\processing_diagnostics\sep_channel\<mouse>.png. (A second,
% duller reason: with a single channel in the folder loadLargeSliceVolume
% squeezes the channel dimension away and generateRegisteredSliceVolume then
% reads slices as channels.)
%
% Cost: no refitting, but the per-slice recovery is an image registration and
% the re-application is a transformix pass per slice per channel, so roughly
% 30-45 min per brain. Run it detached and leave it.
%
% Run sep_setup_paths first, once per MATLAB session. The settings are below,
% the code is in pipeline\add_sep_channel.m.

%% User-defined parameters

% Where the project lives (worked out from where this file sits, see get_paths).
paths = get_paths();

% Cohort selection. Set mice_to_process to {} for every mouse in the groups.
% A mouse without transform_params.mat has never been registered and is
% skipped with a note rather than an error, so all three groups can be given at
% once and the script picks out the brains that are ready.
groups_to_process = {'young', 'naive', 'rws'};
mice_to_process   = {};

% Nothing here needs the GPU or a shared file, so the cohort splits cleanly
% across a few MATLAB sessions run side by side. P4BIS_MICE, a comma-separated
% list in the environment, wins over the list above when it is set -- an
% environment variable rather than a workspace one because of the clear all at
% the top of this file.
%
%   $env:P4BIS_MICE = 'MG903_SepGluA_P20,MG913_SepGluA_P20'
%   matlab -batch "cd('D:\sep_histology\code'); sep_setup_paths; run_add_sep_channel"
env_mice = getenv('P4BIS_MICE');
if ~isempty(env_mice)
    mice_to_process = strtrim(strsplit(env_mice, ','));
end

% Redo a mouse that already has volume_registered_sep.
overwrite = false;

% The per-slice transform is recovered in a cascade, coarsest grid first, one
% affine refinement per level (same world extent each time, fewer pixels). The
% coarse levels cost seconds and bring the fit to within a pixel; the
% full-resolution level, started from there, needs a single pyramid pass and
% lands on the transform exactly. Dropping the 1 would leave a systematic
% half-pixel offset behind -- measured, not assumed, see min_slice_corr.
recovery_levels = [4 2 1];

% A slice whose recovered transform reproduces the aligned DAPI below this
% correlation is reported by number. With the cascade above every slice of
% MG897 came back at r = 1.00000 and a residual shift of 0 px, so anything
% below 0.99 means something is wrong with that slice, not with the method.
min_slice_corr = 0.99;

%% Run

% The settings above go to the code under the same names
run_settings = struct();
run_settings.paths = paths;
run_settings.groups_to_process = groups_to_process;
run_settings.mice_to_process = mice_to_process;
run_settings.overwrite = overwrite;
run_settings.recovery_levels = recovery_levels;
run_settings.min_slice_corr = min_slice_corr;
add_sep_channel(run_settings);
