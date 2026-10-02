function extract_and_center(run_settings)
%EXTRACT_AND_CENTER  Extract and centre the sections of the selected mice.
%   EXTRACT_AND_CENTER(run_settings) does the work of
%   run_extract_and_center, which sets the fields of run_settings
%   (groups_to_process, mice_to_process, atlas_key) and says what each one
%   does.

% settings of run_extract_and_center, under the names the code below uses
groups_to_process = run_settings.groups_to_process;
mice_to_process = run_settings.mice_to_process;
atlas_key = run_settings.atlas_key;

%% Add paths

% the toolboxes, the reader of the raw .czi files included (BioformatsImage and
% its bfmatlab), are on the path from sep_setup_paths; only the atlas is added here
atlas = get_atlas(atlas_key);
atlas_dir = atlas.dir;
addpath(atlas_dir)

%% Resolve cohort

% check that the registry still lists the adults in their legacy order
get_cohort('verify');

% the mice named, or else every mouse of the groups
if isempty(mice_to_process)
    cohort = get_cohort('groups', groups_to_process);
else
    cohort = get_cohort('names', mice_to_process);
end
fprintf('run_extract_and_center: %d mouse/mice selected.\n', numel(cohort));

%% Extract each mouse

% a mouse that fails is reported and the run goes on, so one bad brain does not
% cost an unattended batch; the failures are listed again at the end
failed_mice = {};

for mouse_idx = 1:numel(cohort)

    % the mouse and its folder
    mousename = cohort(mouse_idx).name;
    mouse_dir = cohort(mouse_idx).base_dir;
    fprintf('\n=== %s (group %s) ===\n%s\n', mousename, cohort(mouse_idx).group, mouse_dir);

    try

        if ~exist(mouse_dir, 'dir')
            error('Mouse dir not found: %s\nCopy the raw .czi from the lab share first.', ...
                mouse_dir);
        end

        % read the extraction settings: from the mouse folder, else from its
        % lightsuite folder, else LightSuite's defaults, which the adults ran on
        settings_path = resolve_settings_path(mouse_dir);
        sliceinfo = parseSettingsFile(settings_path);
        fprintf('  settings: %s\n', settings_path);
        fprintf('  slicethickness=%g px_process=%g px_register=%g px_atlas=%g regchan=%s\n', ...
            sliceinfo.slicethickness, sliceinfo.px_process, sliceinfo.px_register, ...
            sliceinfo.px_atlas, sliceinfo.regchan);

        % list the .czi files and find the sections in them (LightSuite)
        sliceinfo.mousename = mousename;
        filelistcheck = dir(fullfile(mouse_dir, '*.czi'));
        if isempty(filelistcheck)
            error('No .czi found in %s\nCopy the raw files from the lab share first.', ...
                mouse_dir);
        end
        filepaths = fullfile({filelistcheck(:).folder}', {filelistcheck(:).name}');
        sliceinfo.filepaths = filepaths;
        fprintf('  %d .czi file(s)\n', numel(filepaths));
        sliceinfo = getSliceInfo(sliceinfo);

        % write every channel, centred, with the ordering composite and sliceinfo.mat
        slicevol = generateSliceVolume(sliceinfo, sliceinfo.regchan); %#ok<NASGU>

        % write volume_ordered.tiff from the decisions file of run_order_slices
        % when there is one, otherwise in the extracted order
        generateReordedVolume(sliceinfo);

        fprintf('  done: %s\n', mousename);

    catch err
        fprintf('  FAILED (%s): %s\n', mousename, err.message);
        failed_mice{end+1} = mousename; %#ok<AGROW>
    end

end

%% Report

if isempty(failed_mice)
    fprintf('\nP1 finished: all %d mouse/mice processed.\n', numel(cohort));
else
    fprintf('\nP1 finished with %d failure(s): %s\n', ...
        numel(failed_mice), strjoin(failed_mice, ', '));
end

end

% ===== Local functions =====

function settings_path = resolve_settings_path(mouse_dir)
% The first local_settings.txt found, in the mouse folder, then in its lightsuite
% folder; when there is none, the path in the mouse folder, so that
% parseSettingsFile warns and applies its defaults (slicethickness 150,
% px_process 5, px_register 20, px_atlas 10, regchan 'dapi'), the values the
% adults were processed with.

candidates = { ...
    fullfile(mouse_dir, 'local_settings.txt'), ...
    fullfile(mouse_dir, 'lightsuite', 'local_settings.txt')};

settings_path = candidates{1};
for k = 1:numel(candidates)
    if exist(candidates{k}, 'file')
        settings_path = candidates{k};
        return
    end
end
end
