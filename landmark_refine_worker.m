function [alive, folder] = landmark_refine_worker(action)
%LANDMARK_REFINE_WORKER Start, stop or query the persistent matcher process.
%
%   landmark_refine_worker('start')    launch it if not already running
%   landmark_refine_worker('stop')     ask it to exit
%   alive = landmark_refine_worker('status')
%
% The worker (landmark_refine/serve.py) imports torch and loads the model
% once, then answers requests dropped into a folder under tempdir. With it
% running, landmark_refine() takes ~0.4 s on the GPU instead of several
% seconds of Python start-up per call. Without it, landmark_refine() still
% works, through cli.py, just slower -- so the GUI never depends on it.
%
% It keeps running after the GUI closes; that is deliberate, so reopening is
% instant. It idles at no CPU. 'stop' it when you are done for the day, or
% just leave it.

if nargin < 1, action = 'status'; end

folder = fullfile(tempdir, 'landmark_refine_worker');
here   = fileparts(mfilename('fullpath'));
pydir  = fullfile(here, 'landmark_refine');

alive = is_alive(folder);

switch lower(action)

    case 'status'
        if nargout == 0
            if alive, fprintf('landmark_refine worker: alive (%s)\n', folder);
            else,     fprintf('landmark_refine worker: not running\n'); end
        end

    case 'start'
        if alive
            if nargout == 0, disp('landmark_refine worker already running.'); end
            return
        end
        py = find_python(pydir);
        if isempty(py)
            warning('landmark_refine_worker:nopython', ...
                'no Python interpreter found: run setup_landmark_refine.ps1');
            return
        end
        if ~exist(folder, 'dir'), mkdir(folder); end
        if exist(fullfile(folder, 'stop'), 'file'), delete(fullfile(folder, 'stop')); end
        logf = fullfile(folder, 'worker.log');
        % pythonw, not python, and start without /B: a console-less process
        % launched through start owns itself, so it survives the shell that
        % system() opened and closed. With /B the child shared that shell's
        % console and, from the MATLAB desktop, died with it two seconds in --
        % pid written, no heartbeat, empty log. The worker writes its own log.
        pyw = strrep(py, 'python.exe', 'pythonw.exe');
        if ~exist(pyw, 'file'), pyw = py; end
        cmd = sprintf('start "" "%s" "%s" "%s"', pyw, fullfile(pydir, 'serve.py'), folder);
        system(cmd);
        fprintf('starting the landmark_refine worker (model load, a few seconds)');
        for k = 1:60                          % up to ~30 s
            pause(0.5); fprintf('.');
            alive = is_alive(folder);
            if alive, break, end
        end
        fprintf('\n');
        if alive
            disp('landmark_refine worker ready.');
        else
            warning('landmark_refine_worker:timeout', ...
                'worker did not report in; see %s', logf);
            if exist(logf, 'file')
                txt = fileread(logf);
                fprintf('--- last lines of the worker log ---\n%s\n', ...
                    strtrim(txt(max(1, end-1500):end)));
            end
        end

    case 'stop'
        if ~alive
            if nargout == 0, disp('landmark_refine worker is not running.'); end
            return
        end
        fclose(fopen(fullfile(folder, 'stop'), 'w'));
        for k = 1:20
            pause(0.25);
            if ~is_alive(folder), break, end
        end
        alive = is_alive(folder);
        if ~alive, disp('landmark_refine worker stopped.'); end

    otherwise
        error('landmark_refine_worker: unknown action ''%s''', action);
end
end


function alive = is_alive(folder)
% A heartbeat younger than a few seconds means a process is looping on it
hb = fullfile(folder, 'heartbeat');
alive = false;
if exist(hb, 'file')
    d = dir(hb);
    alive = (now - d.datenum) * 86400 < 4;    %#ok<TNOW1>
end
end


function py = find_python(pydir)
% An explicit override, else the venv setup_landmark_refine.ps1 creates next
% to the package -- relative to the code folder, so it moves with it.
cands = { getenv('LANDMARK_REFINE_PYTHON'), ...
          fullfile(pydir, '.venv', 'Scripts', 'python.exe') };
py = '';
for k = 1:numel(cands)
    if ~isempty(cands{k}) && exist(cands{k}, 'file'), py = cands{k}; return, end
end
end
