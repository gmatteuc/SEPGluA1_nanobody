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
%   In 'apply' mode generateReordedVolume (LightSuite) reads the decisions
%   file and volume_for_ordering.tiff, and deletes and rewrites
%   volume_ordered.tiff, at the paths it finds in sliceinfo (procpath,
%   volorder). sliceinfo.mat stores them as absolute paths, written by
%   run_extract_and_center; in a copy of the data tree (a check tree, a copy
%   on another drive) they still name the original folder. Both are set to
%   this mouse's folder in the cohort before the call, so 'apply' reads and
%   writes only there, wherever the tree sits.

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

            % rebuild volume_ordered.tiff from the decisions file
            apply_decisions(decisions, procpath);

        otherwise
            error('Unknown run_mode: %s (use ''edit'' or ''apply'').', run_mode);
    end

end

end

% ===== Local functions =====

function apply_decisions(decisions, procpath)
% Mode 'apply': volume_ordered.tiff rebuilt from the decisions file, in the
% mouse's folder in the cohort (procpath).

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

% report the decisions (slices flipped, slices not in their place, slices
% discarded, FlipState -1)
T = readtable(decisions);
fprintf('  decisions: %d slices, %d flipped, %d reordered, %d discarded\n', ...
    height(T), sum(T.FlipState == 1), ...
    sum(T.NewOrderOriginalIndex(:)' ~= 1:height(T)), sum(T.FlipState == -1));

% the mouse's folder in the cohort, in place of the absolute paths stored in
% sliceinfo.mat, which in a copied tree name the original folder (see above)
if ~strcmpi(sliceinfo.procpath, procpath)
    fprintf('  sliceinfo.mat names %s; this mouse''s folder is used instead.\n', ...
        sliceinfo.procpath);
end
sliceinfo.procpath = procpath;
sliceinfo.volorder = fullfile(procpath, 'volume_for_ordering.tiff');

% rebuild volume_ordered.tiff
generateReordedVolume(sliceinfo);
fprintf('  rebuilt: %s\n', fullfile(procpath, 'volume_ordered.tiff'));

end
