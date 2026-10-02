function vol = loadVolume(inArg, ~)
%LOADVOLUME  Read a volume from a multi-page TIFF, or from a folder of TIFFs.
%   vol = LOADVOLUME(inArg) reads every page of the TIFF file inArg, or every
%   .tif file of the folder inArg in natural order (natsortfiles), one plane
%   each, into an H x W x n single array. inArg may also be a cell holding
%   the path; a second input (the callers pass 1) is ignored.

% the callers pass the path in a cell
if iscell(inArg)
    inArg = inArg{1};
end

if isfolder(inArg)

    % a folder: one plane per .tif file, in natural order
    files = dir(fullfile(inArg, '*.tif'));
    files = natsortfiles({files.name});
    Z = numel(files);
    info = imfinfo(fullfile(inArg, files{1}));
    vol = zeros(info.Height, info.Width, Z, 'single');
    for z = 1:Z
        vol(:, :, z) = single(imread(fullfile(inArg, files{z})));
    end
else

    % a multi-page TIFF: one plane per page
    info = imfinfo(inArg);
    Z = numel(info);
    vol = zeros(info(1).Height, info(1).Width, Z, 'single');
    for z = 1:Z
        vol(:, :, z) = single(imread(inArg, z));
    end
end

end
