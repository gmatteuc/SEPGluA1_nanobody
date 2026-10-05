function inputaxh = plot_violinplot(inputadata, inputpars)
%PLOT_VIOLINPLOT  Draw a violin with a box plot for each of several distributions.
%   inputaxh = PLOT_VIOLINPLOT(inputadata, inputpars) draws one violin per
%   distribution of inputadata.inputdistrs (a cell array; an empty cell leaves
%   its place empty) on the axes inputpars.inputaxh, and returns those axes.
%   Each violin is a kernel density estimate (ksdensity, 500 points) mirrored
%   about its centre, scaled to the same half-width; over it, a box from the
%   25th to the 75th percentile, the median as a thick line in the violin's
%   colour, and whiskers to the most extreme values within 1.5 interquartile
%   ranges of the box.
%
%   Fields of inputpars:
%     n_distribs          number of distributions
%     dirstrcenters       x position of each violin
%     boxplotwidth        box width, in x units
%     boxplotlinewidth    line width of box and whiskers (the median is 4 times it)
%     densityplotwidth    half-width of each violin at its widest, in x units
%     xlimtouse           x limits
%     yimtouse            y limits
%     scatterjitter       width of the uniform x jitter of the dots
%     scatteralpha        dot transparency
%     scattersize         dot size
%     xtickslabelvector   tick label of each violin
%     distrcolors         colour of each violin and its median (cell of RGB triplets)
%     distralpha          violin transparency
%     xlabelstring        x label
%     ylabelstring        y label
%     titlestring         title
%     boolscatteron       also draw every value as a black dot
%     ks_bandwidth        kernel bandwidth, in data units
%     inputaxh            the axes to draw on
%
%   Used by run_compare_with_allen_ish. Giulio, 2021.

% unpack the inputs
inputdistrs = inputadata.inputdistrs;
n_distribs = inputpars.n_distribs;
dirstrcenters = inputpars.dirstrcenters;
boxplotwidth = inputpars.boxplotwidth;
boxplotlinewidth = inputpars.boxplotlinewidth;
densityplotwidth = inputpars.densityplotwidth;
xlimtouse = inputpars.xlimtouse;
yimtouse = inputpars.yimtouse;
scatterjitter = inputpars.scatterjitter;
scatteralpha = inputpars.scatteralpha;
scattersize = inputpars.scattersize;
xtickslabelvector = inputpars.xtickslabelvector;
distrcolors = inputpars.distrcolors;
distralpha = inputpars.distralpha;
xlabelstring = inputpars.xlabelstring;
ylabelstring = inputpars.ylabelstring;
titlestring = inputpars.titlestring;
boolscatteron = inputpars.boolscatteron;
ks_bandwidth = inputpars.ks_bandwidth;
inputaxh = inputpars.inputaxh;

% kernel density of each distribution
ks_y = cell(1, n_distribs);
ks_x = cell(1, n_distribs);
ks_n_bins = 500;

% off: the bandwidth is an input (inputpars.ks_bandwidth)
% ks_bandwidth=0.05;

faces = cell(1, n_distribs);
vertsA = cell(1, n_distribs);
vertsB = cell(1, n_distribs);
for distribs_idx = 1:n_distribs
    if not(isempty(inputdistrs{distribs_idx}))

        % estimate the density
        [ks_y{distribs_idx}, ks_x{distribs_idx}] = ksdensity(inputdistrs{distribs_idx}, ...
            'NumPoints', ks_n_bins, 'bandwidth', ks_bandwidth);

        % rescale its height to the violin's half-width
        ks_y{distribs_idx} = (ks_y{distribs_idx} ./ max(ks_y{distribs_idx})) ...
            * densityplotwidth;

        % faces joining each pair of neighbouring density points to the same
        % points at zero height
        qqq = (1:ks_n_bins - 1)';
        faces{distribs_idx} = [qqq, qqq + 1, qqq + ks_n_bins + 1, qqq + ks_n_bins];
    end
end

% patch vertices from the density, x and y swapped so the violins stand upright
for distribs_idx = 1:n_distribs
    if not(isempty(inputdistrs{distribs_idx}))

        % left side of the violin
        vertsA{distribs_idx} = fliplr([ks_x{distribs_idx}', ...
            -ks_y{distribs_idx}' + dirstrcenters(distribs_idx); ...
            ks_x{distribs_idx}', ones(ks_n_bins, 1) * dirstrcenters(distribs_idx)]);

        % right side of the violin
        vertsB{distribs_idx} = fliplr([ks_x{distribs_idx}', ...
            +ks_y{distribs_idx}' + dirstrcenters(distribs_idx); ...
            ks_x{distribs_idx}', ones(ks_n_bins, 1) * dirstrcenters(distribs_idx)]);
    end
end

% box plot statistics of each distribution
quartiles = cell(1, n_distribs);
iqr = cell(1, n_distribs);
Xs = cell(1, n_distribs);
whiskers = cell(1, n_distribs);
Y = cell(1, n_distribs);
box_pos = cell(1, n_distribs);
for distribs_idx = 1:n_distribs
    if not(isempty(inputdistrs{distribs_idx}))

        % 25th and 75th percentiles and median, and the most extreme values
        % within 1.5 interquartile ranges of the box
        quartiles{distribs_idx} = quantile(inputdistrs{distribs_idx}, [0.25 0.75 0.5]);
        iqr{distribs_idx} = quartiles{distribs_idx}(2) - quartiles{distribs_idx}(1);
        Xs{distribs_idx} = sort(inputdistrs{distribs_idx});
        temp_w1 = min(Xs{distribs_idx}(Xs{distribs_idx} > ...
            (quartiles{distribs_idx}(1) - (1.5 * iqr{distribs_idx}))));
        temp_w2 = max(Xs{distribs_idx}(Xs{distribs_idx} < ...
            (quartiles{distribs_idx}(2) + (1.5 * iqr{distribs_idx}))));

        % a whisker with no value in its range ends at the box
        if not(isempty(temp_w1))
            whiskers{distribs_idx}(1) = temp_w1;
        else
            whiskers{distribs_idx}(1) = quartiles{distribs_idx}(1);
        end
        if not(isempty(temp_w2))
            whiskers{distribs_idx}(2) = temp_w2;
        else
            whiskers{distribs_idx}(2) = quartiles{distribs_idx}(2);
        end

        % Y: 25th and 75th percentiles, median, lower and upper whisker
        Y{distribs_idx} = [quartiles{distribs_idx}, whiskers{distribs_idx}];

        % box position: left, bottom, width, height
        box_pos{distribs_idx} = [ ...
            dirstrcenters(distribs_idx) - (boxplotwidth * 0.5), ...
            Y{distribs_idx}(1), ...
            boxplotwidth, ...
            Y{distribs_idx}(2) - Y{distribs_idx}(1) ...
            ];
    end
end

% draw on the given axes
currax = inputaxh;
hold(currax, 'on')
for distribs_idx = 1:n_distribs
    if not(isempty(inputdistrs{distribs_idx}))

        % the two halves of the violin
        patch(currax, 'Faces', faces{distribs_idx}, 'Vertices', vertsA{distribs_idx}, ...
            'FaceVertexCData', distrcolors{distribs_idx}, 'FaceColor', 'flat', ...
            'EdgeColor', 'none', 'FaceAlpha', distralpha);
        patch(currax, 'Faces', faces{distribs_idx}, 'Vertices', vertsB{distribs_idx}, ...
            'FaceVertexCData', distrcolors{distribs_idx}, 'FaceColor', 'flat', ...
            'EdgeColor', 'none', 'FaceAlpha', distralpha);
    end
end
if boolscatteron
    for distribs_idx = 1:n_distribs
        if not(isempty(inputdistrs{distribs_idx}))

            % every value as a dot, jittered in x
            scatter(currax, dirstrcenters(distribs_idx) + (scatterjitter * ...
                rand(size(inputdistrs{distribs_idx})) - scatterjitter / 2), ...
                inputdistrs{distribs_idx}, ...
                'MarkerFaceColor', [0, 0, 0], 'MarkerEdgeColor', [0, 0, 0], ...
                'MarkerFaceAlpha', scatteralpha, 'MarkerEdgeAlpha', 0, ...
                'SizeData', scattersize);
        end
    end
end
for distribs_idx = 1:n_distribs
    if not(isempty(inputdistrs{distribs_idx}))

        % box
        hrect = rectangle(currax, 'Position', box_pos{distribs_idx});
        set(hrect, 'EdgeColor', [0, 0, 0])
        set(hrect, 'LineWidth', boxplotlinewidth);

        % median line
        line(currax, ...
            [box_pos{distribs_idx}(1), box_pos{distribs_idx}(1) + boxplotwidth], ...
            [Y{distribs_idx}(3) Y{distribs_idx}(3)], 'col', distrcolors{distribs_idx}, ...
            'LineWidth', 4 * boxplotlinewidth);

        % whiskers
        line(currax, [dirstrcenters(distribs_idx), dirstrcenters(distribs_idx)], ...
            [Y{distribs_idx}(2) Y{distribs_idx}(5)], 'col', [0, 0, 0], ...
            'LineWidth', boxplotlinewidth);
        line(currax, [dirstrcenters(distribs_idx), dirstrcenters(distribs_idx)], ...
            [Y{distribs_idx}(1) Y{distribs_idx}(4)], 'col', [0, 0, 0], ...
            'LineWidth', boxplotlinewidth);
    end
end

% axis limits, labels and ticks
xlim(currax, xlimtouse)
ylim(currax, yimtouse)
ylabel(ylabelstring)
xlabel(xlabelstring)
title(titlestring)
xticks(dirstrcenters);
xticklabelstouse = xtickslabelvector;
xticklabels(xticklabelstouse);
set(currax, 'fontsize', 12)
hold(currax, 'off')

end
