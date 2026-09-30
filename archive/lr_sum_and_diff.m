function lr_sum_and_diff(inTiff, outSumTiff, outDiffTiff)

    info = imfinfo(inTiff);
    N = numel(info);
    if N == 0, error('No pages found in %s', inTiff); end

    % Read first slice to infer sizes/channels and LR axis
    I1 = single(imread(inTiff, 1));      
    sz = size(I1);
    if numel(sz) == 2, sz(3) = 1; end
    rows = sz(1); cols = sz(2); spp = sz(3);

    % Determine LR axis (NB: as the dimension equal to 800, to be fixed if changing resolution)
    dimsEqual800 = find([rows, cols] == 800);
    if isempty(dimsEqual800)
        error(['Auto-detect failed: neither rows (%d) nor cols (%d) equals 800. ', ...
               'Ensure your LR dimension is 800.'], rows, cols);
    elseif numel(dimsEqual800) > 1
        xDim = 2;
    else
        xDim = dimsEqual800;
    end

    % Prepare outputs
    if exist(outSumTiff,'file'),  delete(outSumTiff);  end
    if exist(outDiffTiff,'file'), delete(outDiffTiff); end
    ts = Tiff(outSumTiff,'w8');
    td = Tiff(outDiffTiff,'w8');

    % Pet tags from a sample frame
    function set_tags(t, frame)
        r = size(frame,1); c = size(frame,2);
        if ndims(frame) >= 3, sp = size(frame,3); else, sp = 1; end %#ok<ISMAT>
        tag.ImageLength          = r;
        tag.ImageWidth           = c;
        tag.Compression          = Tiff.Compression.LZW;
        tag.Photometric          = Tiff.Photometric.MinIsBlack;
        tag.PlanarConfiguration  = Tiff.PlanarConfiguration.Chunky;
        tag.SamplesPerPixel      = sp;
        if sp == 1
            tag.BitsPerSample    = 32;
            tag.SampleFormat     = Tiff.SampleFormat.IEEEFP;
        else
            tag.BitsPerSample    = repmat(32,1,sp);
            tag.SampleFormat     = repmat(Tiff.SampleFormat.IEEEFP,1,sp);
        end
        bytesPerPixel = 4 * sp;
        rps = max(16, floor(2^16 / max(1, c * bytesPerPixel)));
        tag.RowsPerStrip = min(r, rps);
        t.setTag(tag);
    end

    try
        for z = 1:N
            I = single(imread(inTiff, z));   
            F = flip(I, xDim);
            valid = isfinite(I) & isfinite(F);
            % |L+R|
            Sframe = nan(size(I), 'single');
            Sframe(valid) = abs(I(valid) + F(valid));
            % |L-R|
            Dframe = nan(size(I), 'single');
            Dframe(valid) = abs(I(valid) - F(valid));
            % Midline handling for odd LR size: set to 0
            curSz = size(I);
            if mod(curSz(xDim), 2) == 1
                mid = (curSz(xDim)+1)/2;
                idx = repmat({':'}, 1, ndims(I)); idx{xDim} = mid;
                Sframe(idx{:}) = 0;
                Dframe(idx{:}) = 0;
            end
            if z > 1
                ts.writeDirectory();
                td.writeDirectory();
            end
            set_tags(ts, Sframe); ts.write(Sframe);
            set_tags(td, Dframe); td.write(Dframe);
        end
    catch ME
        try ts.close(); end
        try td.close(); end 
        rethrow(ME);
    end

    ts.close();
    td.close();
end

