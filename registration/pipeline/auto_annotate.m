function out = auto_annotate(mode, procpath, slice, plane)
%AUTO_ANNOTATE  Automatic control points for one brain, from Python.
%   out = AUTO_ANNOTATE('propose', procpath) runs the whole brain: it reads the
%   sections, the GUI's warped atlas (auto_atlas_planes.mat) and the anchor
%   planes (plane_anchors.mat) from the lightsuite folder, and writes
%   auto_proposal_controlpoints.mat and auto_proposal_info.mat there. Never
%   atlas2histology_tform.mat: the GUI loads the proposal for review, and only
%   what is accepted there is saved as the annotation.
%
%   out = AUTO_ANNOTATE('section', procpath, slice, plane) proposes one section
%   at the given atlas plane (both 1-based), for the GUI when a proposed
%   section's plane is changed by hand. out.atlas and out.hist are n x 2
%   [y x], out.low marks the least confident points.
%
%   out = AUTO_ANNOTATE('sections', procpath, slices, planes) proposes several
%   slices at once, each at its plane (vectors, 1-based), for the GUI's U.
%   out.atlas, out.hist and out.low are cells, one per slice.
%
%   out = AUTO_ANNOTATE('check') runs nothing: out.ok says whether the engine
%   is installed, by the test every other mode starts with (its interpreter is
%   found). run_register_to_atlas's 'annotate' mode gives the GUI the
%   automatic annotation's keys (auto_annotation_plugin) only then.
%
%   out.ok is false, with out.message saying why, on any failure: a missing
%   interpreter, a missing file, a Python error. Callers check it.
%
%   The Python side lives in registration/auto_annotation/ (get_paths'
%   auto_annotation) and runs in its own virtual environment; run
%   registration/auto_annotation/setup.ps1 once per machine, or point
%   AUTO_ANNOTATION_PYTHON at an interpreter that has torch.

out = struct('ok', false, 'message', '');

% the engine's folder and its interpreter
p     = get_paths();
pydir = p.auto_annotation;
py    = auto_annotation_python(pydir);
if isempty(py)
    out.message = ['no Python interpreter found: run registration\auto_annotation\setup.ps1 ' ...
                   'or set AUTO_ANNOTATION_PYTHON'];
    return
end
cli = fullfile(pydir, 'cli.py');

switch mode
    case 'check'
        % the interpreter was found above; nothing to run
        out.ok = true;

    case 'propose'
        % the proposal file must exist afterwards, whatever the exit status
        cmd = sprintf('"%s" "%s" propose "%s"', py, cli, procpath);
        [status, log] = system(cmd, '-echo');
        if status ~= 0 ...
                || ~exist(fullfile(procpath, 'auto_proposal_controlpoints.mat'), 'file')
            out.message = sprintf('python failed (status %d):\n%s', status, strtrim(log));
            return
        end
        out.ok = true;

    case 'section'
        % the answer comes back in a temporary file, deleted on return
        respf   = [tempname '_resp.mat'];
        cleaner = onCleanup(@() delete_quiet(respf));
        cmd = sprintf('"%s" "%s" section "%s" %d %d "%s"', py, cli, procpath, ...
                      round(slice), round(plane), respf);
        [status, log] = system(cmd);
        if ~exist(respf, 'file')
            out.message = sprintf('python produced no response (status %d):\n%s', status, strtrim(log));
            return
        end
        r = load(respf);
        out.ok      = logical(r.ok);
        out.message = strtrim(char(r.message));
        if out.ok
            out.atlas = double(r.atlas);
            out.hist  = double(r.hist);
            out.low   = logical(r.low(:));
        end

    case 'sections'
        % here slice and plane are vectors: every orange slice, at its plane; the
        % request and the answer go through temporary files, deleted on return
        reqf    = [tempname '_req.mat'];
        respf   = [tempname '_resp.mat'];
        cleaner = onCleanup(@() delete_quiet({reqf, respf}));

        % saved by name, as cli.py reads them
        slices = round(slice(:));
        planes = round(plane(:));
        save(reqf, 'slices', 'planes', '-v7');
        cmd = sprintf('"%s" "%s" sections "%s" "%s" "%s"', py, cli, procpath, ...
            reqf, respf);
        [status, log] = system(cmd);
        if ~exist(respf, 'file')
            out.message = sprintf('python produced no response (status %d):\n%s', status, strtrim(log));
            return
        end
        r = load(respf);
        out.ok      = logical(r.ok);
        out.message = strtrim(char(r.message));
        if out.ok
            out.atlas = r.atlas;
            out.hist = r.hist;
            out.low = r.low;
        end

    otherwise
        out.message = sprintf('unknown mode ''%s'' (check | propose | section | sections)', mode);
end
end

% ===== Local functions =====

function py = auto_annotation_python(pydir)
% The interpreter: AUTO_ANNOTATION_PYTHON when set, else the environment
% setup.ps1 makes in the engine's folder, which moves with the code.

cands = { getenv('AUTO_ANNOTATION_PYTHON'), ...
          fullfile(pydir, '.venv', 'Scripts', 'python.exe') };
py = '';
for k = 1:numel(cands)
    if ~isempty(cands{k}) && exist(cands{k}, 'file')
        py = cands{k};
        return
    end
end
end

function delete_quiet(f)
% Delete the files that exist of f, a file name or a cell of them.

if ischar(f)
    f = {f};
end
for k = 1:numel(f)
    if exist(f{k}, 'file')
        delete(f{k});
    end
end
end
