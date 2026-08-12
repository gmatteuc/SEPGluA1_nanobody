close all
clear all
clc

% /// Pipeline script #0: copy raw .czi from the lab share to local storage ///
% For each selected mouse:
%   (1) Resolves the source dir on the share, honouring the per-mouse
%       share_subdir in the registry (some brains keep their .czi under
%       Anatomy\Axioscan, others at the mouse root)
%   (2) Copies *.czi into <base_root>\<group>\<name>\ so every downstream
%       script sees ONE layout, with the .czi at the mouse root
%   (3) Verifies each file arrived with a byte-identical size
%
% WHY THIS EXISTS: the raw data on the share is READ-ONLY, and getSliceInfo
% creates its 'lightsuite' working folder NEXT TO the .czi it is given. Point
% P1 at the share and it would write there. Copying first is mandatory, not
% stylistic.
%
% Transfer uses robocopy (restartable, resumes rather than restarts). No /MIR
% and no /MOV are ever passed, so the source cannot be modified.

%% User-defined parameters

% Cohort selection (mice come from the shared registry get_cohort.m).
% Set mice_to_process to {} to copy every mouse in groups_to_process.
groups_to_process = {'young'};
mice_to_process   = {'MG909_SepGluA_P20', 'MG910_SepGluA_P20', 'MG911_SepGluA_P16', ...
                     'MG912_SepGluA_P20', 'MG913_SepGluA_P20', 'MG914_SepGluA_P28'};

% Root of the raw data on the lab share (READ-ONLY - never written to)
share_root = 'S:\ElboustaniLab\#SHARE\Data';

% Set false for a dry run that reports what would be copied
do_copy = true;

%% Resolve cohort

get_cohort('verify');
if isempty(mice_to_process)
    cohort = get_cohort('groups', groups_to_process);
else
    cohort = get_cohort('names', mice_to_process);
end
fprintf('P0: %d mouse/mice selected.\n', numel(cohort));

%% Loop over mice

n_ok = 0;
n_bad = 0;

for mouse_idx = 1:numel(cohort)

    mousename = cohort(mouse_idx).name;
    src_dir   = fullfile(share_root, mousename, cohort(mouse_idx).share_subdir);
    dst_dir   = cohort(mouse_idx).base_dir;

    fprintf('\n=== %s (group %s) ===\n', mousename, cohort(mouse_idx).group);
    fprintf('  src: %s\n', src_dir);
    fprintf('  dst: %s\n', dst_dir);

    % Guard: never write anywhere on the share
    assert_local_destination(dst_dir, share_root);

    src_files = dir(fullfile(src_dir, '*.czi'));
    if isempty(src_files)
        warning('No .czi found under %s -- skipping.', src_dir);
        n_bad = n_bad + 1;
        continue
    end
    src_bytes = sum([src_files.bytes]);
    fprintf('  %d .czi, %.2f GB\n', numel(src_files), src_bytes/1024^3);

    if ~do_copy
        fprintf('  (dry run, nothing copied)\n');
        continue
    end

    if ~exist(dst_dir, 'dir')
        mkdir(dst_dir);
    end

    % /Z restartable, /R:3 retries, /W:10 wait, /NP quiet progress,
    % /NDL no dir list, /NJH no job header. Deliberately NO /MIR and NO /MOV.
    cmd = sprintf('robocopy "%s" "%s" *.czi /Z /R:3 /W:10 /NP /NDL /NJH', src_dir, dst_dir);
    t0 = tic;
    [status, out] = system(cmd);
    fprintf('%s', out);
    fprintf('  robocopy exit %d, %.1f min\n', status, toc(t0)/60);

    % robocopy uses 0-7 for success, 8+ for failure
    if status >= 8
        warning('robocopy reported failure (exit %d) for %s.', status, mousename);
        n_bad = n_bad + 1;
        continue
    end

    % Verify every source file arrived at the same size
    ok = true;
    for k = 1:numel(src_files)
        d = fullfile(dst_dir, src_files(k).name);
        if ~exist(d, 'file')
            fprintf('  MISSING  %s\n', src_files(k).name);
            ok = false;
        else
            dinfo = dir(d);
            if dinfo.bytes ~= src_files(k).bytes
                fprintf('  SIZE MISMATCH  %s: src %d vs dst %d\n', ...
                    src_files(k).name, src_files(k).bytes, dinfo.bytes);
                ok = false;
            end
        end
    end

    if ok
        fprintf('  VERIFIED: %d/%d files, %.2f GB\n', numel(src_files), numel(src_files), src_bytes/1024^3);
        n_ok = n_ok + 1;
    else
        n_bad = n_bad + 1;
    end

end

fprintf('\n%s\n', repmat('=', [1 60]));
fprintf('P0 done: %d mouse/mice verified, %d with problems.\n', n_ok, n_bad);
fprintf('%s\n', repmat('=', [1 60]));

%% Local function: refuse to write anywhere on the read-only share

function assert_local_destination(dst_dir, share_root)
% The raw data share must never be written to. Fail loudly rather than risk it.

share_drive = upper(extractBefore([share_root ':'], ':'));
dst_drive   = upper(extractBefore([dst_dir ':'], ':'));

if strcmp(dst_drive, share_drive)
    error(['Refusing to write to the raw-data share.\n' ...
           '  destination: %s\n  share root : %s\n' ...
           'Raw acquisition data is read-only.'], dst_dir, share_root);
end
end
