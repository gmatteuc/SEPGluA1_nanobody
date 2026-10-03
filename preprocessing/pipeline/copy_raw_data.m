function copy_raw_data(run_settings)
%COPY_RAW_DATA  Copy the raw .czi files of the selected mice from the lab share.
%   COPY_RAW_DATA(run_settings) does the work of run_copy_raw_data, which
%   sets the fields of run_settings (groups_to_process, mice_to_process,
%   share_root, do_copy) and says what each one does.

% settings of run_copy_raw_data, under the names the code below uses
groups_to_process = run_settings.groups_to_process;
mice_to_process = run_settings.mice_to_process;
share_root = run_settings.share_root;
do_copy = run_settings.do_copy;

%% Resolve cohort

% check that the registry still lists the adults in their legacy order
get_cohort('verify');

% the mice named, or else every mouse of the groups
if isempty(mice_to_process)
    cohort = get_cohort('groups', groups_to_process);
else
    cohort = get_cohort('names', mice_to_process);
end
fprintf('run_copy_raw_data: %d mouse/mice selected.\n', numel(cohort));

%% Copy each mouse

% mice copied and checked, and mice with a problem
n_ok = 0;
n_bad = 0;

for mouse_idx = 1:numel(cohort)

    % the mouse's folder on the share, and its local folder
    mousename = cohort(mouse_idx).name;
    src_dir = fullfile(share_root, mousename, cohort(mouse_idx).share_subdir);
    dst_dir = cohort(mouse_idx).base_dir;

    fprintf('\n=== %s (group %s) ===\n', mousename, cohort(mouse_idx).group);
    fprintf('  src: %s\n', src_dir);
    fprintf('  dst: %s\n', dst_dir);

    % never write anywhere on the share
    assert_local_destination(dst_dir, share_root);

    % the .czi files on the share; a mouse without any is skipped
    src_files = dir(fullfile(src_dir, '*.czi'));
    if isempty(src_files)
        warning('No .czi found under %s -- skipping.', src_dir);
        n_bad = n_bad + 1;
        continue
    end
    src_bytes = sum([src_files.bytes]);
    fprintf('  %d .czi, %.2f GB\n', numel(src_files), src_bytes/1024^3);

    % a dry run stops here
    if ~do_copy
        fprintf('  (dry run, nothing copied)\n');
        continue
    end

    if ~exist(dst_dir, 'dir')
        mkdir(dst_dir);
    end

    % copy the .czi files
    status = robocopy_czi(src_dir, dst_dir);

    % robocopy exits with 0 to 7 on success, 8 or more on failure
    if status >= 8
        warning('robocopy reported failure (exit %d) for %s.', status, mousename);
        n_bad = n_bad + 1;
        continue
    end

    % check that every file on the share arrived with the same size
    ok = verify_copy(src_files, dst_dir);

    if ok
        fprintf('  VERIFIED: %d/%d files, %.2f GB\n', numel(src_files), ...
            numel(src_files), src_bytes/1024^3);
        n_ok = n_ok + 1;
    else
        n_bad = n_bad + 1;
    end

end

%% Report

fprintf('\n%s\n', repmat('=', [1 60]));
fprintf('run_copy_raw_data done: %d mouse/mice verified, %d with problems.\n', ...
    n_ok, n_bad);
fprintf('%s\n', repmat('=', [1 60]));

end

% ===== Local functions =====

function status = robocopy_czi(src_dir, dst_dir)
% The .czi files copied from src_dir into dst_dir by robocopy; its exit status.

% /Z restartable, /R:3 three retries, /W:10 ten seconds between them, /NP no
% progress, /NDL no folder list, /NJH no job header; never /MIR or /MOV
cmd = sprintf('robocopy "%s" "%s" *.czi /Z /R:3 /W:10 /NP /NDL /NJH', ...
    src_dir, dst_dir);
t0 = tic;
[status, out] = system(cmd);
fprintf('%s', out);
fprintf('  robocopy exit %d, %.1f min\n', status, toc(t0)/60);

end

function ok = verify_copy(src_files, dst_dir)
% Whether every file on the share arrived in dst_dir with its size; prints the others.

% check that every file on the share arrived with the same size
ok = true;
for k = 1:numel(src_files)
    dst_file = fullfile(dst_dir, src_files(k).name);
    if ~exist(dst_file, 'file')
        fprintf('  MISSING  %s\n', src_files(k).name);
        ok = false;
    else
        dinfo = dir(dst_file);
        if dinfo.bytes ~= src_files(k).bytes
            fprintf('  SIZE MISMATCH  %s: src %d vs dst %d\n', ...
                src_files(k).name, src_files(k).bytes, dinfo.bytes);
            ok = false;
        end
    end
end

end

function assert_local_destination(dst_dir, share_root)
% Stop if the destination is on the drive of the raw-data share, which is read
% only.

share_drive = upper(extractBefore([share_root ':'], ':'));
dst_drive = upper(extractBefore([dst_dir ':'], ':'));

if strcmp(dst_drive, share_drive)
    error(['Refusing to write to the raw-data share.\n' ...
           '  destination: %s\n  share root : %s\n' ...
           'Raw acquisition data is read-only.'], dst_dir, share_root);
end
end
