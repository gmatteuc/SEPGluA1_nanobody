function set_lr_colormap(clim_values)
%SET_LR_COLORMAP  Colormap of one panel of a left-right video, on the current axes.
%   SET_LR_COLORMAP(clim_values) gives the current axes the blue-red difference
%   map when the colour limits clim_values are symmetric (jet if that map
%   cannot be made), and hot otherwise.
%
%   Run by write_lr_video, write_lr_video_surpmask and
%   write_lr_video_surpmask_rolling.

if abs(clim_values(1)) == abs(clim_values(2))
    try
        colormap(gca, sep_palette('difference'));
    catch
        colormap(gca, jet);
    end
else
    colormap(gca, sep_palette('intensity'));
end

end
