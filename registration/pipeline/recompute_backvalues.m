function backvalues = recompute_backvalues(input_vol)
%RECOMPUTE_BACKVALUES  Background value of each channel and slice, as LightSuite has it.
%   backvalues = RECOMPUTE_BACKVALUES(input_vol) takes a slice volume
%   (H x W x channel x slice) and returns a channel x slice uint16 array: the
%   1st percentile of each image's nonzero pixels, 0 for an image with none.
%   LightSuite keeps these as sliceinfo.backvalues and fills with them when it
%   warps a slice; they are recomputed for channels LightSuite did not extract
%   (the bridged ones of register_to_atlas, the SEP of add_sep_channel).

% the 1st percentile of the nonzero pixels of each channel and slice
backvalues = zeros(size(input_vol,3), size(input_vol,4), 'uint16');
for c = 1:size(input_vol,3)
    for s = 1:size(input_vol,4)
        img = input_vol(:, :, c, s);
        non_zero_vals = double(img(img > 0));
        if ~isempty(non_zero_vals)
            q = quantile(non_zero_vals, 0.01, 'all');
            backvalues(c, s) = uint16(q);
        else
            backvalues(c, s) = 0;
        end
    end
end

end
