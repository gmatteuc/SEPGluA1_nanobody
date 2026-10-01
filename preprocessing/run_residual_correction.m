clear all
close all
clc

% /// Preprocessing, step 4 of 6: scale the autofluorescence onto the nano channel ///
% For each selected mouse, slice by slice, the autofluorescence channel (Cy3)
% is fitted to the nano channel (Cy5) on reference pixels, and the scaled
% autofluorescence is subtracted from the nano:
%   (1) loads the centered nano and autofluorescence volumes
%       (lightsuite\volume_centered\chan02_Cy5.tiff and chan03_Cy3.tiff)
%   (2) per slice, picks reference pixels and the background on the nano
%       image (select_reference_pixels) and fits nano against
%       autofluorescence on the reference pixels (robustfit, bisquare): a
%       slope and an intercept
%   (3) applies the fit two ways, 'slicewise' (each slice's own) and
%       'global' (the mean over the slices): the scaled autofluorescence,
%       and the nano minus it, with negative values set to 0
%   (4) saves, per way, in lightsuite\correction_output\:
%       corrected_volume_<type>.mat (read by run_register_to_atlas),
%       scaled_auto_volume_<type>.mat (read by run_annotate_artifacts and
%       run_register_to_atlas), and a video of the relative difference,
%       (nano - scaled auto) / scaled auto
% The figures of each slice's fit go to correction_output\diagnostic_plots\,
% with the figures of the reference-pixel selection when savePlotBkg is on
% (doPlotBkg draws them; with it off the run stops, as select_reference_pixels
% then returns no figure). A video of the nano / autofluorescence ratio is
% written when saveRatioMap is on.
%
% Run sep_setup_paths first, once per MATLAB session. The settings are below,
% the code is in pipeline\residual_correction.m.

%% User-defined parameters

% Where the project lives. Derived from the location of the code rather than
% written out, so the tree can be moved or copied to another drive as is.
paths = get_paths();

% Cohort selection (mice come from the shared registry get_cohort.m).
% Set mice_to_process to {} to process every mouse in groups_to_process.
groups_to_process = {'young'};                  % 'rws' | 'naive' | 'behavior' | 'young'
mice_to_process   = {'MG909_SepGluA_P20', 'MG910_SepGluA_P20', 'MG911_SepGluA_P16', ...
                     'MG912_SepGluA_P20', 'MG913_SepGluA_P20', 'MG914_SepGluA_P28'};                         % {} = all mice in groups_to_process

% Reference atlas
atlas_key = 'ccf';

doPlotBkg = true;
savePlotBkg = true;
saveRatioMap = false;

%% Run

% The settings above go to the code under the same names
run_settings = struct();
run_settings.paths = paths;
run_settings.groups_to_process = groups_to_process;
run_settings.mice_to_process = mice_to_process;
run_settings.atlas_key = atlas_key;
run_settings.doPlotBkg = doPlotBkg;
run_settings.savePlotBkg = savePlotBkg;
run_settings.saveRatioMap = saveRatioMap;
residual_correction(run_settings);
