function out = landmark_refine(mode, atlas_im, hist_im, atlas_xy, labels, opts)
%LANDMARK_REFINE Automatic landmark help for the control-point GUI.
%
%   out = landmark_refine('refine',  atlas_im, hist_im, atlas_xy, labels)
%   out = landmark_refine('suggest', atlas_im, hist_im, [],       labels, opts)
%
% 'refine' is the adjusted carry-over: atlas_xy are the previous slice's
% atlas landmarks (Nx2, [x y] in pixels of atlas_im). They are nudged onto
% the nearest structure boundary of this plane, and the histology positions
% are re-derived through the image match field. 'suggest' proposes points
% from scratch; opts.n (default 10) and opts.mirror (default true).
%
% atlas_im and hist_im are the two panels exactly as the GUI displays them,
% uint8 HxW, same size. labels is the annotation plane for atlas_im, or []
% to skip boundary snapping. All coordinates are [x y], i.e. [column row];
% the GUI stores [plane y x t], so convert at the call site.
%
% Returns a struct with ok, message, atlas_pts, hist_pts and per-point
% diagnostics (see landmark_refine/README.md). Never throws on a matching
% failure: check out.ok and fall back to a plain copy.
%
% The Python side lives in landmark_refine/ next to this file and runs in its
% own virtual environment; run setup_landmark_refine.ps1 once per machine.
% One .mat file each way and a single system() call -- measured to cost
% ~1.5 s over the ~2.5 s of matching, so nothing fancier was worth having.

if nargin < 5, labels = []; end
if nargin < 6, opts = struct(); end

out = struct('ok', false, 'message', '', 'atlas_pts', zeros(0,2), 'hist_pts', zeros(0,2));

here   = fileparts(mfilename('fullpath'));
pydir  = fullfile(here, 'landmark_refine');
py     = landmark_refine_python(pydir);
if isempty(py)
    out.message = ['no Python interpreter found: run setup_landmark_refine.ps1 or set ' ...
                   'LANDMARK_REFINE_PYTHON'];
    return
end

% -------------------------------------------------------------- request
req.mode      = char(mode);
req.atlas     = uint8(atlas_im);
req.hist      = uint8(hist_im);
req.atlas_pts = double(reshape(atlas_xy, [], 2));
if ~isempty(labels), req.labels = double(labels); else, req.labels = []; end
if isfield(opts, 'snap_r'), req.snap_r = double(opts.snap_r); end
if isfield(opts, 'n'),      req.n      = double(opts.n);      end
if isfield(opts, 'mirror'), req.mirror = double(opts.mirror); end

tmp   = tempname;                                   % unique per call
reqf  = [tmp '_req.mat'];
respf = [tmp '_resp.mat'];
save(reqf, '-struct', 'req', '-v7');
cleaner = onCleanup(@() delete_quiet({reqf, respf}));

% ----------------------------------------------------------------- call
cmd = sprintf('"%s" "%s" "%s" "%s"', py, fullfile(pydir, 'cli.py'), reqf, respf);
[status, log] = system(cmd);

if ~exist(respf, 'file')
    out.message = sprintf('python produced no response (status %d):\n%s', status, strtrim(log));
    return
end
r = load(respf);

% ---------------------------------------------------------------- unpack
out.ok      = logical(r.ok);
out.message = strtrim(char(r.message));
out.atlas_pts = double(r.atlas_pts);
out.hist_pts  = double(r.hist_pts);
for f = {'moved', 'n_local', 'local_rms', 'score', 'n_matches', 'n_inliers', 'det', ...
         'mean_conf', 'seconds'}
    if isfield(r, f{1}), out.(f{1}) = double(r.(f{1})); end
end
if isfield(out, 'moved'), out.moved = logical(out.moved); end
if ~out.ok && isempty(out.message)
    out.message = sprintf('python failed (status %d):\n%s', status, strtrim(log));
end
end


function py = landmark_refine_python(pydir)
% The interpreter, in order of preference: an explicit override, the venv the
% setup script creates, then the sandbox venv the experiments were run in.
cands = { getenv('LANDMARK_REFINE_PYTHON'), ...
          fullfile(pydir, '.venv', 'Scripts', 'python.exe'), ...
          'D:\sep_histology\sandbox_landmark_matching\.venv\Scripts\python.exe' };
py = '';
for k = 1:numel(cands)
    if ~isempty(cands{k}) && exist(cands{k}, 'file')
        py = cands{k};
        return
    end
end
end


function delete_quiet(files)
for k = 1:numel(files)
    if exist(files{k}, 'file'), delete(files{k}); end
end
end
