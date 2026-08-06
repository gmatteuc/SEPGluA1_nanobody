function backvalues = recompute_backvalues(input_vol)

% Compute backvalues as in Lightsuite
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

