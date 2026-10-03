%% explore_czi_G
% ===== Print the metadata of .czi files =====
%
% A tool, run by hand. For each .czi file named below, prints the number of
% series and of scenes (with each scene's name, centre and size), and for each
% series its dimensions, voxel sizes, channel names, the stage position of each
% tile (or, for a stitched image, the reader's optimal tile size) and the global
% metadata keys about scenes, tiles, positions, overlap and grid; it also shows
% each series' first plane. The files are read with Bio-Formats (bfGetReader), on
% the path from sep_setup_paths (third_party\BioformatsImage).
%
% Setup: three files of MG705 on the lab share. Run sep_setup_paths first, once
% per MATLAB session.

%% Settings

% the .czi files to look into: their folder and their names; they are only read,
% so the folder can be on the lab share (raw data, read only) or a mouse folder of
% the local copy (data\<group>\<mouse>\, see run_copy_raw_data)
folderPath = 'S:\ElboustaniLab\#SHARE\Data\MG705_Gria1\Anatomy\Axioscan\20250706\';
fileNames = {'MG705_SEP_nAB_2WD_1.czi', 'MG705_SEP_nAB_2WD_2.czi', ...
    'MG705_SEP_nAB_2WD_3.czi'};

%% Explore the files

for f = 1:length(fileNames)
    filePath = fullfile(folderPath, fileNames{f});
    fprintf('Exploring file: %s\n', fileNames{f});

    % open the reader and its metadata; the global metadata (a Java hash map) is
    % the file's, the same for every series
    reader = bfGetReader(filePath);
    omeMeta = reader.getMetadataStore();
    globalMeta = reader.getGlobalMetadata();

    % number of series
    numSeries = reader.getSeriesCount();
    fprintf('Number of series: %d\n', numSeries);

    % number of scenes, and each one's name, centre and size, from the original
    % metadata: the adult and the young files alike name them
    % 'Information|Image|SizeS' and 'Information|Image|S|Scene|<field> #<scene>'
    % (a missing key reads as empty, so as NaN)
    sizeS = str2double(char(globalMeta.get('Information|Image|SizeS')));
    if isnan(sizeS)
        fprintf('Scene information (SizeS) not found in metadata.\n');
    else
        fprintf('Number of scenes (SizeS): %d\n', sizeS);

        % the scene number in a key is zero-padded to the digits of the count
        % ('#01' to '#17' for 17 scenes)
        nDigits = numel(num2str(sizeS));
        for scene = 1:sizeS
            sceneKey = sprintf('Information|Image|S|Scene|%%s #%0*d', nDigits, scene);
            sceneName = char(globalMeta.get(sprintf(sceneKey, 'Name')));
            centre = scene_pair(globalMeta, sprintf(sceneKey, 'CenterPosition'));
            extent = scene_pair(globalMeta, sprintf(sceneKey, 'ContourSize'));
            fprintf('  Scene %d (%s): centre %.2f, %.2f µm, size %.0f x %.0f µm\n', ...
                scene, sceneName, centre(1), centre(2), extent(1), extent(2));
        end
    end

    for s = 1:numSeries

        % Bio-Formats counts series from 0
        reader.setSeries(s-1);

        % dimensions; the plane count is Z x C x T
        sizeX = reader.getSizeX();
        sizeY = reader.getSizeY();
        sizeZ = reader.getSizeZ();
        sizeC = reader.getSizeC();
        sizeT = reader.getSizeT();
        numImages = reader.getImageCount();

        fprintf('\nSeries %d:\n', s);
        fprintf('  Dimensions: X=%d, Y=%d, Z=%d, C=%d, T=%d, Total Planes=%d\n', ...
                sizeX, sizeY, sizeZ, sizeC, sizeT, numImages);

        % voxel sizes, if recorded
        try
            voxelX = omeMeta.getPixelsPhysicalSizeX(s-1).value( ...
                ome.units.UNITS.MICROMETER).doubleValue();
            voxelY = omeMeta.getPixelsPhysicalSizeY(s-1).value( ...
                ome.units.UNITS.MICROMETER).doubleValue();
            voxelZ = omeMeta.getPixelsPhysicalSizeZ(s-1).value( ...
                ome.units.UNITS.MICROMETER).doubleValue();
            fprintf('  Voxel sizes (µm): X=%.4f, Y=%.4f, Z=%.4f\n', voxelX, voxelY, ...
                voxelZ);
        catch
            fprintf('  Voxel sizes not available.\n');
        end

        % channel names
        fprintf('  Channels:\n');
        for c = 1:sizeC
            channelName = char(omeMeta.getChannelName(s-1, c-1));
            if isempty(channelName)
                channelName = 'Unnamed';
            end
            fprintf('    Channel %d: %s\n', c, channelName);
        end

        % more planes than channels: tiles or positions, so print each one's stage
        % position; otherwise a stitched image, read by regions
        if numImages > sizeC
            fprintf('  Multiple planes detected - likely individual tiles:\n');
            for p = 1:numImages

                % the plane's index, from 0
                iPlane = p - 1;

                % its stage position
                try
                    posX = omeMeta.getPlanePositionX(s-1, iPlane).value( ...
                        ome.units.UNITS.MICROMETER).doubleValue();
                    posY = omeMeta.getPlanePositionY(s-1, iPlane).value( ...
                        ome.units.UNITS.MICROMETER).doubleValue();
                    posZ = omeMeta.getPlanePositionZ(s-1, iPlane).value( ...
                        ome.units.UNITS.MICROMETER).doubleValue();
                    fprintf('    Plane/Tile %d: Position X=%.2f µm, Y=%.2f µm, Z=%.2f µm\n', ...
                        p, posX, posY, posZ);
                catch
                    fprintf('    Plane/Tile %d: Position not available.\n', p);
                end

                % its time from the start, if recorded (nothing printed otherwise)
                try
                    deltaT = omeMeta.getPlaneDeltaT(s-1, iPlane).doubleValue();
                    fprintf('      DeltaT: %.2f s\n', deltaT);
                catch
                end
            end
        else
            fprintf('  Single plane per channel - likely stitched image. Use region reading for sub-tiles.\n');

            % the tile size the reader reads best, for reading by regions
            optTileW = reader.getOptimalTileWidth();
            optTileH = reader.getOptimalTileHeight();
            fprintf('  Optimal tile size for reading: Width=%d, Height=%d\n', optTileW, ...
                optTileH);
        end

        % show the series' first plane
        reader.setSeries(s-1);
        img = bfGetPlane(reader, 1);
        figure;
        imshow(img, []);
        title(sprintf('Series %d, Plane 1', s));

        % the global metadata keys about tiling, with their values
        fprintf('Analyzing relevant metadata keys (containing Scene, Tile, Position, Overlap, Grid):\n');

        if ~isempty(globalMeta)

            % collect the values of every key about scenes, tiles, positions,
            % overlap or grid, in a containers.Map: a key such as
            % 'Global Information|Image|SizeS #1' is not a valid struct field name
            metaKeys = cell(globalMeta.keySet().toArray());
            relevantKeys = {};
            keyValues = containers.Map('KeyType', 'char', 'ValueType', 'any');
            for k = 1:length(metaKeys)
                key = char(metaKeys{k});
                if contains(lower(key), {'scene', 'tile', 'position', 'overlap', 'grid'})
                    value = char(globalMeta.get(key));
                    relevantKeys{end+1} = key; %#ok<SAGROW>
                    if isKey(keyValues, key)
                        keyValues(key) = [keyValues(key), {value}];
                    else
                        keyValues(key) = {value};
                    end
                end
            end

            % print each key once, with how often it appears and its distinct values
            uniqueKeys = unique(relevantKeys);
            fprintf('Found %d unique relevant keys:\n', length(uniqueKeys));
            for i = 1:length(uniqueKeys)
                key = uniqueKeys{i};
                values = keyValues(key);
                uniqueValues = unique(values);

                fprintf('  Key: %s\n', key);
                fprintf('    Number of occurrences: %d\n', length(values));
                fprintf('    Unique values (%d):\n', length(uniqueValues));
                for v = 1:length(uniqueValues)
                    fprintf('      %s\n', uniqueValues{v});
                end
                fprintf('\n');
            end
        else
            fprintf('No global metadata available or inaccessible.\n');
        end

    end

    % close the reader
    reader.close();

    fprintf('\n--------------------------------------------------\n');
end

% off: examples to adapt by hand, not run
%
% read and save one tile or plane (bfGetPlane counts planes from 1):
% reader = bfGetReader(filePath);
% reader.setSeries(0);  % First series
% iPlane = 1;  % First plane/tile (1-based for bfGetPlane)
% img = bfGetPlane(reader, iPlane);  % Reads the image data for that plane
% figure; imshow(img, []); title('Sample Tile');
% % Save as TIFF
% imwrite(img, 'sample_tile.tif');
% reader.close();
%
% read a region of a large stitched image (planes == C); openBytes(iPlane, x, y,
% w, h) returns a byte array, to convert to a matrix, e.g. img =
% typecast(javaArrayToMatlab(reader.openBytes(0, 0, 0, 512, 512)), 'uint16'):
% reader = bfGetReader(filePath);
% reader.setSeries(0);
% javaMethod('openBytes', reader, 0, 0, 0, 512, 512);
% reader.close();

% ===== Local functions =====

function pair = scene_pair(globalMeta, key)
% The two numbers of a scene field written "x,y" in the global metadata; NaN NaN
% when the key is missing.

pair = [NaN NaN];
value = char(globalMeta.get(key));
if ~isempty(value)
    pair = str2double(strsplit(value, ','));
end
end
