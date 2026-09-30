function vol = loadVolume(inArg, ~)

if iscell(inArg), inArg = inArg{1}; end

if isfolder(inArg)
    files = dir(fullfile(inArg,'*.tif'));
    files = natsortfiles({files.name});
    Z = numel(files);
    info = imfinfo(fullfile(inArg,files{1}));
    vol = zeros(info.Height, info.Width, Z, 'single');
    for z = 1:Z
        vol(:,:,z) = single( imread(fullfile(inArg,files{z})) );
    end
else
    info = imfinfo(inArg);
    Z = numel(info);
    vol = zeros(info(1).Height, info(1).Width, Z, 'single');
    for z = 1:Z
        vol(:,:,z) = single( imread(inArg, z) );
    end
end

end
