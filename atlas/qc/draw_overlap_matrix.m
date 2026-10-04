function draw_overlap_matrix(m, row_labels, col_labels, tick_font, cell_font, x_label)
%DRAW_OVERLAP_MATRIX  Draw an overlap matrix in hot, with its values in the cells.
%   DRAW_OVERLAP_MATRIX(m, row_labels, col_labels, tick_font, cell_font,
%   x_label) draws m (from REGION_OVERLAP_MATRIX, rows the adult CCF regions)
%   on the current axes with limits 0 to 1: the rows labelled row_labels and
%   the columns col_labels, in tick labels of size tick_font; each value of
%   0.02 or more printed in its cell (size cell_font), black on the bright
%   cells and white on the dark ones; a colorbar; x_label says what the
%   columns are.
%
%   Run by the atlas QC scripts check_demba_to_allen and compare_atlas_regions.

imagesc(m, [0 1]);
colormap(gca, sep_palette('intensity'));
axis square
n = size(m, 1);
set(gca, 'XTick', 1:n, 'XTickLabel', col_labels, ...
         'YTick', 1:n, 'YTickLabel', row_labels, ...
         'TickLabelInterpreter', 'none', 'FontSize', tick_font);
xtickangle(45)

% each cell worth reading, in a colour that shows on hot
for i = 1:n
    for j = 1:n
        v = m(i, j);
        if v < 0.02
            continue
        end
        if v > 0.55
            txt_col = [0 0 0];
        else
            txt_col = [1 1 1];
        end
        text(j, i, sprintf('%.2f', v), 'HorizontalAlignment', 'center', ...
             'FontSize', cell_font, 'Color', txt_col);
    end
end

cb = colorbar;
cb.Label.String = 'fraction of the adult region';
xlabel(x_label)
ylabel('adult CCF region')

end
