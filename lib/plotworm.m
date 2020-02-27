function plotworm(obj)

    set(gcf, 'Position', [2068 1 1200 800]);

    % Plotting samples at original position.
    ax1 = subplot('Position', [0.05, 0.125, 0.9, 0.8]);
    [X,Y] = scatterX(obj.contractionSequence(:, :, 1), 'colorAssignment', obj.clusterAssignments(obj.iteration, :),...
        'dimensionalityReductionMethod', 'tsne', ...
        'labels', obj.options.labelfn(obj.clusterAssignments(obj.iteration, :), obj.channels));
        colormap(ax1, distinguishable_colors(length(unique(obj.clusterAssignments(obj.iteration, :)))));

    % Draw convex hull around the groups
    groups = groupby([X,Y], obj.clusterAssignments(obj.iteration, :));
    hull(ax1, groups, distinguishable_colors(max(obj.clusterAssignments(obj.iteration, :))));
    
    plotcolorbar(obj);

    % Plotting Header Line
    subplot('Position', [0.05, 0.925, 0.935, 0.125], 'Visible', 'off')
    toWrite = ['Iteration ' num2str(obj.iteration) ...
                ', \sigma = ' num2str(obj.currentSigma) ...
                ', #Clusters = ' num2str(length(unique(obj.clusterAssignments(obj.iteration, :)))) ...
                ', #Samples = ' num2str(size(obj.dataContracted, 1))];
    text(0, 0.5, toWrite, 'FontSize', 20);
    savefig(obj);
end

function groups = groupby(M, assigments)
    groups = {};
    for group = 1:max(assigments)
        groups{group} = M(find(assigments == group), :);
    end
end

function plotcolorbar(obj)
    bar = subplot('Position', [0.05, 0.05, 0.9, 0.05]);
    branches = max(obj.clusterAssignments(end,:));
    cbranch = 1*ones(length(obj.clusterAssignments(end, obj.clusterAssignments(end,:) == 1)),1);  

    for cluster = 2:branches
       group = obj.clusterAssignments(end, obj.clusterAssignments(end,:) == cluster);
       cbranch=[cbranch; cluster*ones(length(group),1)];
    end
    
    imagesc(bar, cbranch');
    colormap(bar, distinguishable_colors(max(obj.clusterAssignments(obj.iteration, :))));
    set(bar,'xtick', []);
    set(bar,'ytick', []);
end

function hull(fig, groups, colors)
    hold on

    for i = 1:size(groups, 2)
        group = groups{i};
        if ( size(group, 1) > 2)
            k = convhull(group(:, 1), group(:, 2));
            plot(fig, group(k,1), group(k,2),'-', 'Color', colors(i,:));
        end
    end
    
    hold off;
end
