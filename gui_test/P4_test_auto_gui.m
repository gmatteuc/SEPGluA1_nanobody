close all
clear all
clc

% /// Test launcher for the control-point GUI with automatic landmark help ///
%
% Opens matchControlPointsInSlices_auto -- a COPY of the production GUI with
% three additions, nothing else changed:
%
%   p   carry-forward now REFINES the carried points: the atlas side is nudged
%       onto the nearest structure boundary of the new plane, the histology
%       side is re-derived through the image match field (landmark_refine).
%       Falls back to the plain copy if the matcher is unavailable or fails.
%   r   refine the points already on the current slice, in place
%   a   toggle the automatic refinement on/off (on by default)
%
%   every point carries its index number, on both panels, so the pairing is
%   always visible.
%
% The production GUI, the registration code and the pipeline scripts are not
% touched. This writes atlas2histology_tform.mat exactly where the production
% GUI does, so back that file up before testing on a mouse you care about --
% or point it at a mouse whose annotation you are happy to redo.
%
% Expect ~5-7 s per refined carry on CPU: Python start-up plus four LoFTR
% passes. The window title says AUTO-REFINE while it is on.

%% User-defined parameters

% This lives one folder below the code root, so put the root on the path
% before asking it for anything -- then it runs from wherever it is opened.
addpath(fileparts(fileparts(mfilename('fullpath'))));
paths = get_paths();
mouse_to_test = 'MG903_SepGluA_P20';
atlas_key     = 'demba_p20';

%% Paths

addpath(paths.code);
addpath(fullfile(paths.code, 'gui_test'));
addpath(genpath(paths.lightsuite));
addpath(genpath(paths.yaml));
atlas = get_atlas(atlas_key);        % puts the right atlas on the path

%% Open

cohort    = get_cohort('names', {mouse_to_test});
mouse_dir = fullfile(cohort(1).base_dir, 'lightsuite');

regopts_name = fullfile(mouse_dir, 'regopts.mat');
if ~exist(regopts_name, 'file')
    error('no regopts.mat for %s -- run P4 with run_mode = ''align'' first', mouse_to_test);
end
opts = load(regopts_name);

tform_name = fullfile(mouse_dir, 'atlas2histology_tform.mat');
if exist(tform_name, 'file')
    fprintf('NOTE: control points exist and will be overwritten on save:\n  %s\n', tform_name);
end

fprintf('opening the TEST GUI on %s against %s\n', mouse_to_test, atlas_key);
matchControlPointsInSlices_auto(opts);
