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
%     'counts'       magma, 256 levels, for counts of brains (n maps)
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
%     'control_higher', 'experimental_higher'
%                                            a cluster where that group's
%                                            |L - R| is the higher, dark grey
%                                            #4d4d4d and red #c0392b
%     'no_data'                              flat grey #bfbfbf where a map has no
%                                            value, which no data colormap gives
%     'outline'                              blue #3a6db5, a cluster's outline
%                                            over an intensity map
%     'mid_grey', 'note_grey'                a mark and a note drawn on a figure,
%                                            #9a9a9a and 0.4 (plotting.py's
%                                            MID_GREY and NOTE_GREY)
%     'shading_grey'                         #eeeeee, a band shaded behind the
%                                            data

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
    case 'counts'

        % matplotlib's magma at nine even steps, interpolated to 256 levels
        anchors = [
            0.0015 0.0005 0.0139
            0.1131 0.0655 0.2768
            0.3167 0.0717 0.4854
            0.5128 0.1482 0.5076
            0.7164 0.2150 0.4753
            0.9043 0.3196 0.3881
            0.9867 0.5356 0.3822
            0.9969 0.7696 0.5349
            0.9871 0.9914 0.7495
            ];
        c = interp1(linspace(0, 1, 9), anchors, linspace(0, 1, 256));

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

    % a cluster where one group is the higher, apart from the groups' blue and
    % orange (#4d4d4d, #c0392b)
    case 'control_higher'
        c = [77 77 77] / 255;
    case 'experimental_higher'
        c = [192 57 43] / 255;

    % maps: no value, and an outline over them (#bfbfbf, #3a6db5)
    case 'no_data'
        c = [191 191 191] / 255;
    case 'outline'
        c = [58 109 181] / 255;

    % greys for what a figure marks, notes or shades (#9a9a9a, 0.4 and #eeeeee)
    case 'mid_grey'
        c = [154 154 154] / 255;
    case 'note_grey'
        c = [0.4 0.4 0.4];
    case 'shading_grey'
        c = [238 238 238] / 255;

    otherwise
        error('sep_palette: unknown name ''%s''; help sep_palette lists the names', ...
            name);
end

end
