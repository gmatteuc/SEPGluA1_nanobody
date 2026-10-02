function order_slices(run_settings)
%ORDER_SLICES  Curate the slice order of one mouse, or apply the curated order.
%   ORDER_SLICES(run_settings) does the work of run_order_slices, which sets
%   the fields of run_settings (mice_to_process, run_mode) and says what
%   each one, and each run mode, does.
%
%   In 'edit' mode SliceOrderEditor holds MATLAB (uiwait) until its window
%   is closed, so this function returns only then, as the script did. The
%   window keeps what it needs in its own data and writes the decisions
%   file itself, so nothing it does depends on this function's variables.
%
%   In 'apply' mode generateReordedVolume (LightSuite) takes its folder from
%   sliceinfo.mat, not from the cohort: it reads the decisions file and
%   volume_for_ordering.tiff, and deletes and rewrites volume_ordered.tiff,
%   at the absolute paths run_extract_and_center stored there (procpath,
%   volorder). In the production data tree that is this mouse's folder. In
%   a copy of it (a check tree, a copy on another drive) it is still the
%   original folder, which 'apply' then overwrites, and the 'rebuilt:'
%   line below names a file it did not write. Never run 'apply' on a
%   copied mouse folder until this is fixed (it is on the bug list).

% settings of run_order_slices, under the names the code below uses
mice_to_process = run_settings.mice_to_process;
run_mode = run_settings.run_mode;

%% Resolve cohort

% check that the registry still lists the adults in their legacy order
get_cohort('verify');
cohort = get_cohort('names', mice_to_process);

%% Run the selected mode

for mouse_idx = 1:numel(cohort)

    % the mouse, its ordering composite and its decisions file
    mousename = cohort(mouse_idx).name;
    procpath = fullfile(cohort(mouse_idx).base_dir, 'lightsuite');
    volorder = fullfile(procpath, 'volume_for_ordering.tiff');
    decisions = fullfile(procpath, 'volume_for_ordering_processing_decisions.txt');

    fprintf('\n=== %s (group %s) ===\n', mousename, cohort(mouse_idx).group);

    if ~exist(volorder, 'file')
        error('Ordering volume not found:\n  %s\nRun run_extract_and_center for this mouse first.', ...
            volorder);
    end

    switch lower(run_mode)

        case 'edit'

            % one window at a time, which saves the decisions file itself
            if numel(cohort) > 1
                error('run_mode ''edit'' opens one GUI at a time; select a single mouse.');
            end
            if exist(decisions, 'file')
                fprintf('  NOTE: a decisions file already exists and will be overwritten on save:\n    %s\n', ...
                    decisions);
            end
            fprintf('  opening SliceOrderEditor on:\n    %s\n', volorder);
            fprintf('  curate, then SAVE and CLOSE the GUI, set run_mode = ''apply'' and re-run.\n');
            SliceOrderEditor(volorder);

        case 'apply'

            % the decisions file and sliceinfo.mat, both needed
            if ~exist(decisions, 'file')
                error(['No decisions file found:\n  %s\n' ...
                       'Run this script with run_mode = ''edit'' first.'], decisions);
            end
            sliceinfo_name = fullfile(procpath, 'sliceinfo.mat');
            if ~exist(sliceinfo_name, 'file')
                error('sliceinfo.mat not found:\n  %s\nRun run_extract_and_center for this mouse first.', ...
                    sliceinfo_name);
            end
            S = load(sliceinfo_name);
            sliceinfo = S.sliceinfo;

            % report the decisions (slices flipped, slices not in their place)
            T = readtable(decisions);
            fprintf('  decisions: %d slices, %d flipped, %d reordered\n', ...
                height(T), sum(T.FlipState == 1), ...
                sum(T.NewOrderOriginalIndex(:)' ~= 1:height(T)));

            % rebuild volume_ordered.tiff, in the folder sliceinfo.mat names (see above)
            generateReordedVolume(sliceinfo);
            fprintf('  rebuilt: %s\n', fullfile(procpath, 'volume_ordered.tiff'));

        otherwise
            error('Unknown run_mode: %s (use ''edit'' or ''apply'').', run_mode);
    end

end

end
