function c = sep_palette(name)
%SEP_PALETTE  The project's colours and colormaps, one place for every figure.
%   c = SEP_PALETTE(name) returns the colormap (n x 3) or the colour (1 x 3,
%   RGB from 0 to 1) that the figures use for name. docs/STYLE.md (Figures)
%   says what each is for.
%
%   Colormaps:
%     'intensity'    hot: nano, autofluorescence, a reading
%     'anatomy'      gray: the atlas, anatomy, raw images
%     'difference'   blue through white to red, 256 levels, for a difference
%                    or a t map drawn with symmetric limits, so zero is white
%     'bars'         200 greys from light (0.78) to black, for bars coloured
%                    by their value, darker for more; never white
%   hot and gray take the length of the current figure's colormap, as they do
%   when called directly.
%
%   Colours:
%     'nano', 'autofluorescence'             the two channels, orange and yellow
%     'nano_dots', 'autofluorescence_dots'   their per-mouse dots
%     'paired_lines'                         lines joining paired mice
%     'young', 'naive', 'rws'                the groups
%     'control', 'experimental'              the plasticity comparison's two
%                                            groups (MATLAB's default blue and
%                                            orange)
%     'control_mean', 'experimental_mean'    their group means, darker

switch name

    % colormaps
    case 'intensity'
        c = hot;
    case 'anatomy'
        c = gray;
    case 'difference'
        c = get_color2color_colormap([0 0 1], [1 0 0]);
    case 'bars'

        % the darkest 200 of 256 grey levels, lightest first
        levels = gray(256);
        c = flipud(levels(1:200, :));

    % the two channels
    case 'nano'
        c = [0.95 0.55 0.10];
    case 'autofluorescence'
        c = [0.95 0.85 0.20];
    case 'nano_dots'
        c = [0.65 0.30 0.00];
    case 'autofluorescence_dots'
        c = [0.70 0.60 0.00];
    case 'paired_lines'
        c = [0.6 0.6 0.6];

    % the groups (#c0392b, #555555, #9a9a9a)
    case 'young'
        c = [192 57 43] / 255;
    case 'naive'
        c = [85 85 85] / 255;
    case 'rws'
        c = [154 154 154] / 255;

    % the plasticity comparison
    case 'control'
        c = [0 0.45 0.74];
    case 'experimental'
        c = [0.85 0.33 0.10];
    case 'control_mean'
        c = [0 0.2 0.5];
    case 'experimental_mean'
        c = [0.64 0.08 0.18];

    otherwise
        error('sep_palette: unknown name ''%s''; help sep_palette lists the names', ...
            name);
end

end
