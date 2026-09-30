% Script to explore and recover tiles from CZI files using Bio-Formats in MATLAB
% Assumes Bio-Formats toolbox is installed and added to MATLAB's Java path.
% Download from: https://www.openmicroscopy.org/bio-formats/downloads/
% Add to path: javaaddpath('/path/to/bioformats_package.jar');

% Define the folder path and select a file (change as needed)
folderPath = 'S:\ElboustaniLab\#SHARE\Data\MG705_Gria1\Anatomy\Axioscan\20250706\';
fileNames = {'MG705_SEP_nAB_2WD_1.czi', 'MG705_SEP_nAB_2WD_2.czi', 'MG705_SEP_nAB_2WD_3.czi'};

% Loop over each file
for f = 1:length(fileNames)
    filePath = fullfile(folderPath, fileNames{f});
    fprintf('Exploring file: %s\n', fileNames{f});
    
    % Initialize reader and metadata
    reader = bfGetReader(filePath);
    omeMeta = reader.getMetadataStore();
    
    % Get number of series
    numSeries = reader.getSeriesCount();
    fprintf('Number of series: %d\n', numSeries);

    % Extract number of scenes if available (from original metadata)
    try
        sizeS = str2double(globalMeta.get('Global Information|Image|SizeS #1'));
        fprintf('Number of scenes (SizeS): %d\n', sizeS);
        % Get scene positions
        for scene = 1:sizeS
            xKey = sprintf('Global Information|Image|S|Scene|Position|X #%d', scene);
            yKey = sprintf('Global Information|Image|S|Scene|Position|Y #%d', scene);
            zKey = sprintf('Global Information|Image|S|Scene|Position|Z #%d', scene);
            posX = str2double(globalMeta.get(xKey));
            posY = str2double(globalMeta.get(yKey));
            posZ = str2double(globalMeta.get(zKey));
            fprintf('  Scene %d Position: X=%.2f, Y=%.2f, Z=%.2f µm\n', scene, posX, posY, posZ);
        end
    catch
        fprintf('Scene information (SizeS) not found in metadata.\n');
    end
    
    % Loop over each series
    for s = 1:numSeries
        reader.setSeries(s-1);  % 0-based index
        
        % Get dimensions
        sizeX = reader.getSizeX();
        sizeY = reader.getSizeY();
        sizeZ = reader.getSizeZ();
        sizeC = reader.getSizeC();
        sizeT = reader.getSizeT();
        numImages = reader.getImageCount();  % Total planes = Z * C * T
        
        fprintf('\nSeries %d:\n', s);
        fprintf('  Dimensions: X=%d, Y=%d, Z=%d, C=%d, T=%d, Total Planes=%d\n', ...
                sizeX, sizeY, sizeZ, sizeC, sizeT, numImages);
        
        % Get voxel sizes if available
        try
            voxelX = omeMeta.getPixelsPhysicalSizeX(s-1).value(ome.units.UNITS.MICROMETER).doubleValue();
            voxelY = omeMeta.getPixelsPhysicalSizeY(s-1).value(ome.units.UNITS.MICROMETER).doubleValue();
            voxelZ = omeMeta.getPixelsPhysicalSizeZ(s-1).value(ome.units.UNITS.MICROMETER).doubleValue();
            fprintf('  Voxel sizes (µm): X=%.4f, Y=%.4f, Z=%.4f\n', voxelX, voxelY, voxelZ);
        catch
            fprintf('  Voxel sizes not available.\n');
        end
        
        % Get channel names
        fprintf('  Channels:\n');
        for c = 1:sizeC
            channelName = char(omeMeta.getChannelName(s-1, c-1));
            if isempty(channelName)
                channelName = 'Unnamed';
            end
            fprintf('    Channel %d: %s\n', c, channelName);
        end
        
        % If multiple planes, assume they are tiles/scenes; recover positions and parameters
        if numImages > sizeC  % More planes than channels suggests tiles or multi-position
            fprintf('  Multiple planes detected - likely individual tiles:\n');
            for p = 1:numImages
                % Calculate effective plane index (0-based)
                iPlane = p - 1;
                
                % Get position (stage position for the tile)
                try
                    posX = omeMeta.getPlanePositionX(s-1, iPlane).value(ome.units.UNITS.MICROMETER).doubleValue();
                    posY = omeMeta.getPlanePositionY(s-1, iPlane).value(ome.units.UNITS.MICROMETER).doubleValue();
                    posZ = omeMeta.getPlanePositionZ(s-1, iPlane).value(ome.units.UNITS.MICROMETER).doubleValue();
                    fprintf('    Plane/Tile %d: Position X=%.2f µm, Y=%.2f µm, Z=%.2f µm\n', p, posX, posY, posZ);
                catch
                    fprintf('    Plane/Tile %d: Position not available.\n', p);
                end
                
                % Optional: DeltaT or other parameters
                try
                    deltaT = omeMeta.getPlaneDeltaT(s-1, iPlane).doubleValue();
                    fprintf('      DeltaT: %.2f s\n', deltaT);
                catch
                    % Ignore if not available
                end
            end
        else
            fprintf('  Single plane per channel - likely stitched image. Use region reading for sub-tiles.\n');
            % For stitched/large images, get optimal tile size for region reading
            optTileW = reader.getOptimalTileWidth();
            optTileH = reader.getOptimalTileHeight();
            fprintf('  Optimal tile size for reading: Width=%d, Height=%d\n', optTileW, optTileH);
        end

        reader.setSeries(s-1);  % 0-based index
        img = bfGetPlane(reader, 1); 
        figure; 
        imshow(img, []); 
        title(sprintf('Series %d, Plane 1', s));

% Enhanced metadata analysis for tiling information
fprintf('Analyzing relevant metadata keys (containing Scene, Tile, Position, Overlap, Grid):\n');
globalMeta = reader.getGlobalMetadata(); % Returns a Java HashMap or similar

if ~isempty(globalMeta)
    % Convert global metadata keys to a cell array
    metaKeys = cell(globalMeta.keySet().toArray());
    relevantKeys = {};
    keyValues = struct(); % Store key-value pairs
    
    % Collect all relevant key-value pairs
    for k = 1:length(metaKeys)
        key = char(metaKeys{k});
        if contains(lower(key), {'scene', 'tile', 'position', 'overlap', 'grid'})
            value = char(globalMeta.get(key));
            relevantKeys{end+1} = key; %#ok<SAGROW>
            keyValues.(key) = [keyValues.(key), {value}]; % Append value to the key's array
        end
    end
    
    % Remove duplicates from relevantKeys and get unique keys
    uniqueKeys = unique(relevantKeys);
    fprintf('Found %d unique relevant keys:\n', length(uniqueKeys));
    for i = 1:length(uniqueKeys)
        key = uniqueKeys{i};
        values = keyValues.(key); % Get all values for this key
        uniqueValues = unique(values); % Get unique values
        
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
    
    % Close the reader
    reader.close();
    
    fprintf('\n--------------------------------------------------\n');
end

% Example to extract and save a specific tile/plane (uncomment and adjust)
% reader = bfGetReader(filePath);
% reader.setSeries(0);  % First series
% iPlane = 1;  % First plane/tile (1-based for bfGetPlane)
% img = bfGetPlane(reader, iPlane);  % Reads the image data for that plane
% figure; imshow(img, []); title('Sample Tile');
% % Save as TIFF
% imwrite(img, 'sample_tile.tif');
% reader.close();

% For large stitched images (if planes == C), example to read a region/tile:
% reader = bfGetReader(filePath);
% reader.setSeries(0);
% javaMethod('openBytes', reader, 0, 0, 0, 512, 512);  % openBytes(iPlane, x, y, w, h) - returns byte array, convert to matrix
% % To convert byte[] to image: img = typecast(javaArrayToMatlab(reader.openBytes(0, 0, 0, 512, 512)), 'uint16'); adjust type
% reader.close();