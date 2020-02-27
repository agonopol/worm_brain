close all;
clear;
clc;

addpath('./lib');
loaddeps();

options = Options();
options.clusterAssignmentMethod = 'none';
options.epsilonClusterIdentificationMethod = 'constantEpsilon';
options.frequencyMergingEpsilonClusters = 'always'; %always,uponMetastability%
options.controlSigmaMethod = 'nuclearNormStabilization'; %nuclearNormStabilization,movementStabilization
options.numDiffusionSteps = 3;
options.phateEmbedding = false;

files = dir('data/*.csv');

jsh = arrayfun(@(x) sprintf('JSH%03d', x), 1:282, 'UniformOutput', false);
n2u = [arrayfun(@(x) sprintf('N2U_%03d', x), 2:182, 'UniformOutput', false) ...
       arrayfun(@(x) sprintf('N2U_VC_%03d', x), 1:34, 'UniformOutput', false) ];


for file = files'
    
    path = fullfile(file.folder, file.name);
    adj = readworm(path);
    
    if not(isempty(strfind(path, 'N2U')))
        rows = find(arrayfun(@(x1) any(strcmp(x1, n2u)), adj.EMSection));
    else
        rows = find(arrayfun(@(x1) any(strcmp(x1, jsh)), adj.EMSection));
    end
    
    [neurons, adj] = weightedadj(adj(rows, {'Neuron1','Neuron2', 'EMSection', 'Weight'}));
    
    [~, name, ~] = fileparts(path);
    options.destination = fullfile(pwd(), 'results', name, '//');
    [dest, ~, ~] = fileparts(options.destination);
    mkdir_if_not_exists(dest);
    
    condense(adj, neurons, options);
    
    close all force;
    close all hidden;
      
end
