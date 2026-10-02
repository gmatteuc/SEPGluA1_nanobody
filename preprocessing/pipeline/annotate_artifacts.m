function annotate_artifacts(run_settings)
%ANNOTATE_ARTIFACTS  Open the artifact annotation window on each selected mouse.
%   ANNOTATE_ARTIFACTS(run_settings) does the work of
%   run_annotate_artifacts, which sets the fields of run_settings (paths,
%   groups_to_process, mice_to_process, correction_type) and says what each
%   one does.
%
%   ArtifactAnnotator holds MATLAB (uiwait) until its window is closed, then
%   the next mouse opens, as in the script. The window keeps what it needs
%   in its own data and saves the masks itself, so nothing it does depends
%   on this function's variables.
%
%   Each mouse's scaled_auto_volume_<correction_type>.mat is loaded with
%   load() and no output, so the variables passed to the window
%   (scaledautoVol, nanoVol, bg_mask_vol, slice_data) are assigned nowhere
%   in the code: a function finds them only because no function on the
%   path has their names. The file also holds correction_type, which
%   replaces the setting (with the same value) when it is loaded, as it did
%   in the script.

% settings of run_annotate_artifacts, under the names the code below uses
paths = run_settings.paths;
groups_to_process = run_settings.groups_to_process;
mice_to_process = run_settings.mice_to_process;
correction_type = run_settings.correction_type;

%% Resolve cohort

% check that the registry still lists the adults in their legacy order
get_cohort('verify');

% the mice named, or else every mouse of the groups
if isempty(mice_to_process)
    cohort = get_cohort('groups', groups_to_process);
else
    cohort = get_cohort('names', mice_to_process);
end
fprintf('run_annotate_artifacts: %d mouse/mice selected.\n', numel(cohort));

%% Annotate each mouse

for mouse_idx = 1:numel(cohort)

    % the mouse's name and group
    mouse_name = cohort(mouse_idx).name;
    mouse_type = cohort(mouse_idx).group;

    % its correction folder, where run_residual_correction writes and
    % run_register_to_atlas reads
    base_dir = fullfile(paths.data, mouse_type);
    output_dir = fullfile(base_dir, mouse_name, 'lightsuite', 'correction_output');
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    % load what run_residual_correction saved: scaledautoVol, nanoVol, bg_mask_vol,
    % slice_data, average_slope, average_intercept and correction_type
    matfile_name = fullfile(output_dir, sprintf('scaled_auto_volume_%s.mat', ...
        correction_type));
    if ~exist(matfile_name, 'file')
        fprintf('File not found for %s: %s\n', mouse_name, matfile_name);
        continue;
    end
    load(matfile_name);
    [H, W, Z] = size(nanoVol);

    % open the annotation window, which holds MATLAB until it is closed
    ArtifactAnnotator(scaledautoVol, nanoVol, bg_mask_vol, slice_data, mouse_name, ...
        output_dir, correction_type);

    % free the mouse's volumes
    clear scaledautoVol nanoVol bg_mask_vol slice_data average_slope average_intercept

end

end
