# worm  
The algorithm is the following:  

Read the worm csv file by:  
  + Read csv file line by line  
  + Create an zero matrix for every neuron in union(set(neuron 1) and set(neuron 2) column values)  
  + For every EMSection:  
    + Add weight value if exists for neuron 1 and neuron 2 pair  

Given the weighted adjacency matrix create a special embedding of size neuron x neuron:  
  + Create a spacial embedding by taking the eigen decomposition and taking the top 50 eigen vectors of the markov matrix:  
    + Create a marov matrix by row norm decomposition  
      + //markov = rownorm(adj);  
      + Take the eigenvalue decomposition of the markov matrix  
        + [V,D] = eig(markov);  
        + // Create an embedding by creating taking the top 50 eigen vectors corresponding to highest eigen values, creating a  
        + d = min(size(V, 1), 50);  
        + matrix = bsxfun(@times,V(:,1:d)',diag(D(1:d, :)))';  
        
While the number of clusters is > 1  
  + calculate affinity matrix:  
    + For every neuron x neuron pair, create a diffusion using given sigma and calculate the affinities to other neurons  
      + // obj = obj.calcAffinities();  
    + merge clusters:  
      + For every neuron x neuron pair that is less then a given epsilon, collapse the pairs into one cluster  
      + // obj = obj.performContractionMove();  
      + // obj = obj.mergeEpsilonClusters();  
      + // obj = obj.assignClusters();  
    + adjust sigma:  
      + using nuclear norm stabilization:  
      + obj = obj.controlSigma();  
      + if within the last 10 iterations the sum of eigenvalues changed by 5% then increase sigma by 1.1  

To run:  

The main driver is brain.m, this sets up which CSV file to read in the line:  

  + files = dir('data/*.csv'); //change the 'data/*.csv' to a specific 'path to csv file.csv'  

line:  
  + condensed = condense(adj, neurons, options);  

Generates the condensation clustering and its history:  

The history of the merging process is kept in the condensed.contractionSequence(:, :, 1...steps), this is the spacial coordinates of the eigen vector embedding of the input data, this is what you would use if you wanted to plot it with phate  

condensed.clusterAssignments(:, 1...steps) is the cluster assigments after each step, you can use those to check the quality scores and centrality based  

quality comparison is done inside the lib/condense.m file  

given the original adj matrix:  

quality(adj, 'modularity', ...  
            'maxk', 50, ...  
            'output', strcat(condensed.options.destination, 'modularity-comparison.png'), ...  
            'condensation', condensed.clusterAssignments, ...   
            'k-means', kmeanscompare(adj, condensed.clusterAssignments), ...  
            'agglomerative', hierarchicalcompare(adj, condensed.clusterAssignments));  

This compare the quality given the cluster assigments over the number of clusters/steps  
