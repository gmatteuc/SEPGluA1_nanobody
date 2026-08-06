function results = run_transformix(transformix_exe, transform_params, output_dir, varargin)
%RUN_TRANSFORMIX Wrapper to execute the transformix command-line tool.
%
%   results = run_transformix(transformix_exe, transform_params, output_dir, ...)
%
%   Applies a pre-computed transformation (defined by transform_params)
%   using the specified transformix executable, placing outputs in output_dir.
%   Can operate on an input image, a set of points, or both simultaneously.
%
%   Required Inputs:
%       transformix_exe  - Full path to the transformix executable.
%       transform_params - Path(s) to the Transform Parameter file(s) (.txt).
%                          Provide a single string for one file, or a cell
%                          array of strings {'param1.txt', 'param0.txt'}
%                          for chained transforms (order matters: final -> initial).
%       output_dir       - Path to the directory where results will be saved.
%                          Will be created if it doesn't exist.
%
%   Optional Name-Value Pair Inputs:
%       'InputImage'     - Path to the input image file (.mhd, .tif, etc.)
%                          to be transformed. Default: ''.
%       'InputPoints'    - Points to be transformed. Can be:
%                          * An Nx2 or Nx3 matrix of coordinates.
%                          * A path to a text file formatted for transformix's -def.
%                          Default: [].
%       'PointType'      - Type of points if providing a matrix for 'InputPoints'.
%                          Either 'Point' (world coords) or 'Index'.
%                          Default: 'Point'.
%       'Verbose'        - Set to true to display command and status messages.
%                          Default: false.
%       'OtherTransformixArgs' - Cell array of additional arguments to pass
%                                directly to transformix (e.g., {'-threads', '4'}).
%                                Default: {}.
%       'LoadOutputs'    - Set to true to load the transformed image and points
%                          into the results struct. Default: false.
%
%   Outputs (returned in 'results' struct):
%       Status              - Exit status of the transformix command (0 = success).
%       CommandOutput       - Standard output/error text from transformix.
%       CommandString       - The exact command executed.
%       LogFile             - Path to the 'transformix.log' file.
%       TransformedImageFile- Path to the output image file (e.g., 'result.mhd').
%                             (Populated if 'InputImage' was provided).
%       TransformedPointsFile- Path to the 'outputpoints.txt' file.
%                              (Populated if 'InputPoints' was provided).
%       TransformedImage    - Loaded transformed image data (if 'LoadOutputs' is true).
%       TransformedPoints   - Loaded transformed points matrix (if 'LoadOutputs' is true).
%
%   Example Usage:
%   % --- Define inputs ---
%   exe_path = 'C:\path\to\transformix.exe'; % Or '/usr/bin/transformix'
%   tp_file = 'elastix_output/TransformParameters.0.txt';
%   out_dir = 'transformix_output';
%   image_in = 'moving_image.mhd';
%   points_in = [10 20 5; 15 30 8]; % Nx3 matrix
%
%   % --- Run transformix on both image and points ---
%   results = run_transformix(exe_path, tp_file, out_dir, ...
%                             'InputImage', image_in, ...
%                             'InputPoints', points_in, ...
%                             'PointType', 'Point', ...
%                             'Verbose', true, ...
%                             'LoadOutputs', true);
%
%   % --- Check results ---
%   if results.Status == 0
%       disp('Transformix completed successfully.');
%       disp(['Transformed image saved to: ', results.TransformedImageFile]);
%       disp('Transformed points:');
%       disp(results.TransformedPoints);
%   else
%       disp('Transformix failed.');
%       disp(results.CommandOutput);
%   end
%
% Requires helper functions: writePointsFile, readTransformedPointsFile, mhd_read (or similar)

    % --- Input Parsing ---
    p = inputParser;
    p.KeepUnmatched = true; % Allow other args for future use if needed

    addRequired(p, 'transformix_exe', @(x) ischar(x) || isstring(x));
    addRequired(p, 'transform_params', @(x) ischar(x) || isstring(x) || iscellstr(x));
    addRequired(p, 'output_dir', @(x) ischar(x) || isstring(x));

    addParameter(p, 'InputImage', '', @(x) ischar(x) || isstring(x));
    addParameter(p, 'InputPoints', [], @(x) isnumeric(x) || ischar(x) || isstring(x));
    addParameter(p, 'PointType', 'Point', @(x) ismember(lower(x), {'point', 'index'}));
    addParameter(p, 'Verbose', false, @islogical);
    addParameter(p, 'OtherTransformixArgs', {}, @iscell);
    addParameter(p, 'LoadOutputs', false, @islogical);

    parse(p, transformix_exe, transform_params, output_dir, varargin{:});

    % --- Validate Inputs ---
    if ~isfile(p.Results.transformix_exe)
        error('run_transformix:exeNotFound', 'Transformix executable not found at: %s', p.Results.transformix_exe);
    end

    if ischar(p.Results.transform_params) || isstring(p.Results.transform_params)
        transform_param_files = {char(p.Results.transform_params)};
    else % Cell array
        transform_param_files = cellstr(p.Results.transform_params);
    end
    for i = 1:length(transform_param_files)
        if ~isfile(transform_param_files{i})
            error('run_transformix:paramFileNotFound', 'Transform parameter file not found: %s', transform_param_files{i});
        end
    end

    process_image = ~isempty(p.Results.InputImage);
    process_points = ~isempty(p.Results.InputPoints);

    if ~process_image && ~process_points
        error('run_transformix:noInput', 'At least one of InputImage or InputPoints must be provided.');
    end

    if process_image && ~isfile(p.Results.InputImage)
        error('run_transformix:imageNotFound', 'Input image file not found: %s', p.Results.InputImage);
    end

    points_are_matrix = isnumeric(p.Results.InputPoints);
    input_points_file = ''; % Path to the points file transformix will use
    temp_points_file = '';  % Keep track if we created a temporary file

    if process_points
        if points_are_matrix
            if size(p.Results.InputPoints, 2) < 2 || size(p.Results.InputPoints, 2) > 3
                error('run_transformix:invalidPointsMatrix', 'InputPoints matrix must be Nx2 or Nx3.');
            end
            % Create a temporary file for the points matrix
            temp_points_file = fullfile(tempdir, ['transformix_points_' datestr(now,'yymmddHHMMSSFFF') '.txt']);
            try
                writePointsFile(temp_points_file, p.Results.InputPoints, p.Results.PointType);
                input_points_file = temp_points_file;
                 if p.Results.Verbose
                    fprintf('Created temporary points file: %s\n', input_points_file);
                end
            catch ME
                error('run_transformix:writePointsFailed', 'Failed to write temporary points file: %s\n%s', temp_points_file, ME.message);
            end
        else % Points provided as a file path
            input_points_file = char(p.Results.InputPoints);
            if ~isfile(input_points_file)
                error('run_transformix:pointsFileNotFound', 'Input points file not found: %s', input_points_file);
            end
        end
    end

    % --- Prepare Output Directory ---
    output_dir_abs = char(p.Results.output_dir);
    if ~isfolder(output_dir_abs)
        [status, msg, msgID] = mkdir(output_dir_abs);
        if ~status
            error('run_transformix:mkdirFailed', 'Failed to create output directory: %s\n%s (%s)', output_dir_abs, msg, msgID);
        end
        if p.Results.Verbose
             fprintf('Created output directory: %s\n', output_dir_abs);
        end
    end

    % --- Build Command String ---
    % Use full paths with quotes to handle spaces
    cmd = sprintf('"%s"', p.Results.transformix_exe);
    cmd = sprintf('%s -out "%s"', cmd, output_dir_abs);
    cmd = sprintf('%s -tp "%s"', cmd, transform_param_files{1}); % Only the first TP file is passed directly

    % Handle chained transforms: copy TPs to output dir and link them
    if length(transform_param_files) > 1
        copied_tp_files = cell(size(transform_param_files));
        for i = 1:length(transform_param_files)
            [~, fname, ext] = fileparts(transform_param_files{i});
            dest_path = fullfile(output_dir_abs, [fname, ext]);
            try
                copyfile(transform_param_files{i}, dest_path, 'f'); % Overwrite if exists
                copied_tp_files{i} = dest_path;
            catch ME
                 error('run_transformix:copyTpFailed', 'Failed to copy transform parameter file to output directory: %s\n%s', transform_param_files{i}, ME.message);
            end
        end
        % Modify copied files to chain correctly (points to files in output dir)
        for i = 1:(length(copied_tp_files) - 1)
             try
                changeParameterInElastixFile(copied_tp_files{i}, 'InitialTransformParametersFileName', copied_tp_files{i+1});
             catch ME
                 warning('run_transformix:linkTpFailed', 'Failed to automatically link transform parameter file %s.\n Ensure InitialTransformParametersFileName is set correctly or points to "NoInitialTransform".\n%s', copied_tp_files{i}, ME.message);
             end
        end
        % Update command to use the copied (and potentially modified) first TP file
        cmd = sprintf('%s -out "%s" -tp "%s"', sprintf('"%s"', p.Results.transformix_exe), output_dir_abs, copied_tp_files{1});
    end


    if process_image
        cmd = sprintf('%s -in "%s"', cmd, p.Results.InputImage);
    end

    if process_points
        cmd = sprintf('%s -def "%s"', cmd, input_points_file);
    end

    % Add any other arguments
    for i = 1:length(p.Results.OtherTransformixArgs)
        cmd = sprintf('%s %s', cmd, p.Results.OtherTransformixArgs{i});
    end

    % --- Execute Transformix ---
    if p.Results.Verbose
        fprintf('Executing Transformix...\nCommand:\n%s\n', cmd);
    end

    tic;
    [status, cmdout] = system(cmd);
    elapsedTime = toc;

    if p.Results.Verbose
        fprintf('Transformix finished in %.2f seconds.\n', elapsedTime);
        fprintf('Status: %d\n', status);
        disp('--- Transformix Output ---');
        disp(cmdout);
        disp('--------------------------');
    end

    % --- Process Outputs ---
    results.Status = status;
    results.CommandOutput = cmdout;
    results.CommandString = cmd;
    results.LogFile = fullfile(output_dir_abs, 'transformix.log');
    results.TransformedImageFile = '';
    results.TransformedPointsFile = '';
    results.TransformedImage = [];
    results.TransformedPoints = [];


    if status == 0 % Success
        if process_image
            % Find result image (could be .mhd, .tif, etc.)
            d_img = dir(fullfile(output_dir_abs, 'result.*'));
            % Remove .raw file if present
            raw_idx = cellfun(@(x) endsWith(x, '.raw', 'IgnoreCase', true), {d_img.name});
            d_img(raw_idx) = [];
            if ~isempty(d_img)
                results.TransformedImageFile = fullfile(output_dir_abs, d_img(1).name);
                if p.Results.LoadOutputs
                     try
                        if endsWith(results.TransformedImageFile, '.mhd', 'IgnoreCase', true)
                            results.TransformedImage = mhd_read(results.TransformedImageFile);
                        elseif endsWith(results.TransformedImageFile, {'.tif', '.tiff'}, 'IgnoreCase', true)
                             % Assuming load3Dtiff or similar exists
                            results.TransformedImage = load3Dtiff(results.TransformedImageFile);
                        else
                             warning('run_transformix:unknownImageFormat', 'Loading for image format %s not implemented.', results.TransformedImageFile);
                        end
                     catch ME
                         warning('run_transformix:loadImageFailed', 'Failed to load transformed image: %s\n%s', results.TransformedImageFile, ME.message);
                     end
                end
            elseif p.Results.Verbose
                warning('run_transformix:noResultImage', 'Transformix succeeded but no result image file found in %s.', output_dir_abs);
            end
        end

        if process_points
            results.TransformedPointsFile = fullfile(output_dir_abs, 'outputpoints.txt');
            if isfile(results.TransformedPointsFile)
                if p.Results.LoadOutputs
                    try
                        results.TransformedPoints = readTransformedPointsFile(results.TransformedPointsFile);
                    catch ME
                         warning('run_transformix:loadPointsFailed', 'Failed to load transformed points file: %s\n%s', results.TransformedPointsFile, ME.message);
                    end
                end
            elseif p.Results.Verbose
                 warning('run_transformix:noResultPoints', 'Transformix succeeded but outputpoints.txt not found in %s.', output_dir_abs);
            end
        end
    else % Failure
        warning('run_transformix:executionFailed', 'Transformix command failed with status %d. Check CommandOutput and LogFile for details.', status);
    end

    % --- Cleanup ---
    if ~isempty(temp_points_file) && isfile(temp_points_file)
        delete(temp_points_file);
        if p.Results.Verbose
            fprintf('Deleted temporary points file: %s\n', temp_points_file);
        end
    end

end % function run_transformix