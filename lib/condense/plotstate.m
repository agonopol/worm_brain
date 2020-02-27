function plotstate(obj)

    set(gcf, 'Position', [2068 1 1200 800]);
    if (obj.options.phateEmbedding)
        % Plotting samples at original position.
        ax1 = subplot('Position', [0.05, 0.125, 0.425, 0.8]);
        [~, npca] = size(obj.contractionSequence(:, :, 1));
        embedding = phate(obj.contractionSequence(:, :, 1), 'npca', npca, 'mds_method', 'cmds');
        scatterX(embedding, 'colorAssignment', obj.clusterAssignments(obj.iteration, :));
        colormap(ax1, distinguishable_colors(length(unique(obj.clusterAssignments(obj.iteration, :)))));
        % Plotting samples at contracted position.
        ax2 = subplot('Position', [0.525, 0.125, 0.425, 0.8]);
        scatterX(embedding);
        colormap(ax2, 'gray');
        hold on;
        [centroids, sizes] = stats(embedding, obj.clusterAssignments(obj.iteration, :));
        scatter(centroids(:, 1), centroids(:, 2), sizes, distinguishable_colors(max(obj.clusterAssignments(obj.iteration, :))), 'o', 'filled');
        hold off;
    else
        % Plotting samples at original position.
        ax1 = subplot('Position', [0.05, 0.125, 0.425, 0.8]);
        scatterX(obj.contractionSequence(:, :, 1), 'colorAssignment', obj.clusterAssignments(obj.iteration, :),...
            'dimensionalityReductionMethod', 'tsne', ...
            'sizeAssignment', obj.options.sizefn(obj.clusterAssignments(obj.iteration, :), obj.channels), ...
            'labels', obj.options.labelfn(obj.clusterAssignments(obj.iteration, :), obj.channels));
            colormap(ax1, distinguishable_colors(length(unique(obj.clusterAssignments(obj.iteration, :)))));
        % Plotting samples at contracted position.
        ax2 = subplot('Position', [0.525, 0.125, 0.425, 0.8]);
        sizeAssignment = sqrt(cellfun(@size, obj.sampleIndices, repmat({2}, 1, length(obj.sampleIndices))));
        scatterX(obj.dataContracted, ...
            'colorAssignment', 1:max(obj.clusterAssignments(obj.iteration, :)), ...
            'sizeAssignment', sizeAssignment');
            colormap(ax2, distinguishable_colors(max(obj.clusterAssignments(obj.iteration, :))));
    end
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
    % Plotting Header Line
    subplot('Position', [0.05, 0.925, 0.935, 0.125], 'Visible', 'off')
    toWrite = ['Iteration ' num2str(obj.iteration) ...
                ', \sigma = ' num2str(obj.currentSigma) ...
                ', #Clusters = ' num2str(length(unique(obj.clusterAssignments(obj.iteration, :)))) ...
                ', #Samples = ' num2str(size(obj.dataContracted, 1))];
    text(0, 0.5, toWrite, 'FontSize', 20);
    savefig(obj);
end


