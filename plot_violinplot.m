function inputaxh = plot_violinplot(inputadata,inputpars)
% Giulio 2021

% reassign inputs
inputdistrs=inputadata.inputdistrs;
n_distribs=inputpars.n_distribs;
dirstrcenters=inputpars.dirstrcenters;
boxplotwidth=inputpars.boxplotwidth;
boxplotlinewidth=inputpars.boxplotlinewidth;
densityplotwidth=inputpars.densityplotwidth;
xlimtouse=inputpars.xlimtouse;
yimtouse=inputpars.yimtouse;
scatterjitter=inputpars.scatterjitter;
scatteralpha=inputpars.scatteralpha;
scattersize=inputpars.scattersize;
xtickslabelvector=inputpars.xtickslabelvector;
distrcolors=inputpars.distrcolors;
distralpha=inputpars.distralpha;
xlabelstring=inputpars.xlabelstring;
ylabelstring=inputpars.ylabelstring;
titlestring=inputpars.titlestring;
boolscatteron=inputpars.boolscatteron;
ks_bandwidth=inputpars.ks_bandwidth;
inputaxh=inputpars.inputaxh;


% initialize density estimation variables
ks_y=cell(1,n_distribs);
ks_x=cell(1,n_distribs);
ks_n_bins=500;
% ks_bandwidth=0.05;
faces=cell(1,n_distribs);
vertsA=cell(1,n_distribs);
vertsB=cell(1,n_distribs);
for distribs_idx=1:n_distribs
    if not(isempty(inputdistrs{distribs_idx}))
        % estimate density using 'ksdensity' function
        [ks_y{distribs_idx}, ks_x{distribs_idx}] = ksdensity(inputdistrs{distribs_idx}, 'NumPoints', ks_n_bins, 'bandwidth', ks_bandwidth);
        % rescale height of estimated density
        ks_y{distribs_idx}=(ks_y{distribs_idx}./max(ks_y{distribs_idx}))*densityplotwidth;
        % define the faces to connect each adjacent estimated density points and the corresponding points at zero height
        qqq = (1:ks_n_bins - 1)';
        faces{distribs_idx} = [qqq, qqq + 1, qqq + ks_n_bins + 1, qqq + ks_n_bins];
    end
end
% calculate patch vertices from kernel density and plip x and y
for distribs_idx = 1:n_distribs
    if not(isempty(inputdistrs{distribs_idx}))
        % left side of the violin
        vertsA{distribs_idx} = fliplr([ks_x{distribs_idx}', -ks_y{distribs_idx}' + dirstrcenters(distribs_idx); ks_x{distribs_idx}', ones(ks_n_bins, 1) * dirstrcenters(distribs_idx)]);
        % right side of the violin
        vertsB{distribs_idx} = fliplr([ks_x{distribs_idx}', +ks_y{distribs_idx}' + dirstrcenters(distribs_idx); ks_x{distribs_idx}', ones(ks_n_bins, 1) * dirstrcenters(distribs_idx)]);
    end
end
% initialize boxplot summary statistics variables
quartiles=cell(1,n_distribs);
iqr=cell(1,n_distribs);
Xs=cell(1,n_distribs);
whiskers=cell(1,n_distribs);
Y=cell(1,n_distribs);
box_pos=cell(1,n_distribs);
for distribs_idx=1:n_distribs
    if not(isempty(inputdistrs{distribs_idx}))
        % compute boxplot summary statistics
        quartiles{distribs_idx} = quantile(inputdistrs{distribs_idx}, [0.25 0.75 0.5]);
        iqr{distribs_idx} = quartiles{distribs_idx}(2) - quartiles{distribs_idx}(1);
        Xs{distribs_idx} = sort(inputdistrs{distribs_idx});
        temp_w1=min( Xs{distribs_idx}(Xs{distribs_idx} > (quartiles{distribs_idx}(1) - (1.5*iqr{distribs_idx}))) );
        temp_w2=max( Xs{distribs_idx}(Xs{distribs_idx} < (quartiles{distribs_idx}(2) + (1.5*iqr{distribs_idx} ))) );
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
        Y{distribs_idx} = [quartiles{distribs_idx}, whiskers{distribs_idx}];
        % fill in box position vector
        box_pos{distribs_idx} = [...
            dirstrcenters(distribs_idx)-(boxplotwidth * 0.5),...
            Y{distribs_idx}(1),...
            boxplotwidth,...
            Y{distribs_idx}(2)-Y{distribs_idx}(1)...
            ];
    end
end
% get current axis handle
currax=inputaxh;
hold(currax,'on')
for distribs_idx = 1:n_distribs
    if not(isempty(inputdistrs{distribs_idx}))
        % draw violin patches
        patch(currax,'Faces', faces{distribs_idx}, 'Vertices', vertsA{distribs_idx},...
            'FaceVertexCData', distrcolors{distribs_idx}, 'FaceColor', 'flat', 'EdgeColor', 'none', 'FaceAlpha', distralpha);
        patch(currax,'Faces', faces{distribs_idx}, 'Vertices', vertsB{distribs_idx},...
            'FaceVertexCData', distrcolors{distribs_idx}, 'FaceColor', 'flat', 'EdgeColor', 'none', 'FaceAlpha', distralpha);
    end
end
if boolscatteron
    for distribs_idx = 1:n_distribs
        if not(isempty(inputdistrs{distribs_idx}))
            % draw scatter dots
            scatter(currax,dirstrcenters(distribs_idx)+(scatterjitter*rand(size(inputdistrs{distribs_idx}))-scatterjitter/2),inputdistrs{distribs_idx},...
                'MarkerFaceColor', [0,0,0], 'MarkerEdgeColor', [0,0,0],'MarkerFaceAlpha', scatteralpha, 'MarkerEdgeAlpha', 0, 'SizeData', scattersize);
        end
    end
end
for distribs_idx = 1:n_distribs
    if not(isempty(inputdistrs{distribs_idx}))
        % draw box
        hrect = rectangle(currax,'Position', box_pos{distribs_idx});
        set(hrect, 'EdgeColor', [0,0,0])
        set(hrect, 'LineWidth', boxplotlinewidth);
        % draw median line
        line(currax,[box_pos{distribs_idx}(1),box_pos{distribs_idx}(1)+boxplotwidth],[Y{distribs_idx}(3) Y{distribs_idx}(3)],'col',distrcolors{distribs_idx},'LineWidth',4*boxplotlinewidth);
        % draw whiskers
        line(currax,[dirstrcenters(distribs_idx),dirstrcenters(distribs_idx)],[Y{distribs_idx}(2) Y{distribs_idx}(5)],'col',[0,0,0],'LineWidth',boxplotlinewidth);
        line(currax,[dirstrcenters(distribs_idx),dirstrcenters(distribs_idx)],[Y{distribs_idx}(1) Y{distribs_idx}(4)],'col',[0,0,0],'LineWidth',boxplotlinewidth);
    end
end
% embellish axes and add labels
xlim(currax,xlimtouse)
ylim(currax,yimtouse)
ylabel(ylabelstring)
xlabel(xlabelstring)
title(titlestring)
xticks(dirstrcenters);
xticklabelstouse=xtickslabelvector;
xticklabels(xticklabelstouse);
set(currax,'fontsize',12)
hold(currax,'off')

end

