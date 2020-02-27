function condensed = condense(adj, neurons, options)

    markov = rownorm(adj);
    [V,D] = eig(markov);
    d = min(size(V, 1), 50);
    matrix = bsxfun(@times,V(:,1:d)',diag(D(1:d, :)))';
    

    options.sizefn = @(clusters, labels) centrality(adj, clusters)';
    options.labelfn = @(clusters, labels) label(labels, centrality(adj, clusters)');
    options.plotfn = @plotworm;

    condenser = ContractionClustering(matrix, cellstr(neurons), options);
    condensed = condenser.contract();
    
    % Write out cluster assigments
    
    f = fopen(strcat(options.destination, 'cluster_assigments.csv'), 'w');
    fprintf(f, join(string(neurons'), ','));
    fprintf(f, '\n');

    for i = 1:min(size(condensed.clusterAssignments))
        fprintf(f, join(string(condensed.clusterAssignments(i,:)), ','));
        fprintf(f, '\n');
    end
    fclose(f);

    % Write out cluster centrality analysis
    
    f = fopen(strcat(options.destination, 'cluster_centrality.csv'), 'w');
    fprintf(f, join(string(neurons'), ','));
    fprintf(f, '\n');
    
    
    for i = 1:min(size(condensed.clusterAssignments))
        fprintf(f, join(string(centrality(adj, condensed.clusterAssignments(i,:))), ','));
        fprintf(f, '\n');
    end        
    fclose(f);
    
    quality(adj, 'modularity', ...
            'maxk', 50, ...
            'output', strcat(condensed.options.destination, 'modularity-comparison.png'), ...
            'condensation', condensed.clusterAssignments, ... 
            'k-means', kmeanscompare(adj, condensed.clusterAssignments), ...
            'agglomerative', hierarchicalcompare(adj, condensed.clusterAssignments));
    

end