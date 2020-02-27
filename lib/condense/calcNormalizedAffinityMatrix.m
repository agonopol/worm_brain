function [normalizedAffinityMatrix, Q] = calcNormalizedAffinityMatrix(distanceMatrix, zeroIndicatorMatrix, varargin)
% CALCNORMALIZEDAFFINITYMATRIX Returns the row stochastic affinity matrix
%   given the distance matrix.
%
%   Author: Tobias Welp, modifying a version by David v. Dijk
%   Year  : 2016 based on original work from 2015

    % default arguments
    sigma = 1;
    exponent = 2;
    mode = 'rstochastic';

    % parsing of variable argument list
    for i=1:length(varargin)-1
        if (strcmp(varargin{i}, 'sigma'))
            sigma = varargin{i+1};
        end
        if (strcmp(varargin{i}, 'exponent'))
            exponent = varargin{i+1};
        end
        if (strcmp(varargin{i}, 'weights'))
            weights = varargin{i+1};
        end
        if (strcmp(varargin{i}, 'mode'))
            mode = varargin{i+1};
            if ~ismember(mode, {'rstochastic', 'cstochastic', 'fstochastic'})
                warning(['calcNormalAffinityMatrix: Invalid choice for ' ...
                         'mode. Using default value "rstochastic"']);
                mode = 'rstochastic';
            end
        end
    end

    normalizedAffinityMatrix = calcAffinityMatrix(distanceMatrix, ...
                                                  zeroIndicatorMatrix, ...
                                                  'sigma', sigma, ...
                                                  'exponent', exponent);

    n = size(normalizedAffinityMatrix, 1);

    if (~exist('weights', 'var'))
        weights = ones(size(distanceMatrix, 1), 1);
    end

    weightedNormalizedAffinityMatrix = normalizedAffinityMatrix * sparse(1:n , 1:n, weights);
    switch (mode)
        case 'rstochastic'
            rowSums = sum(weightedNormalizedAffinityMatrix, 2);
            degreeMatrix = sparse(1:n, 1:n, rowSums);
            Q = degreeMatrix;
            normalizedAffinityMatrix = (degreeMatrix^-1) * normalizedAffinityMatrix;
        case 'cstochastic'
            columnSums = sum(weightedNormalizedAffinityMatrix, 1);
            degreeMatrix = sparse(1:n, 1:n, columnSums);
            Q = degreeMatrix;
            normalizedAffinityMatrix = normalizedAffinityMatrix * (degreeMatrix^-1);
        case 'fstochastic'
            rowSums = sum(weightedNormalizedAffinityMatrix, 2);
            columnSums = sum(weightedNormalizedAffinityMatrix, 1);
            sqrtDegreesRows = sparse(1:n, 1:n, sqrt(rowSums));
            sqrtDegreesColumns = sparse(1:n, 1:n, sqrt(columnSums));
            Q = sqrtDegreesRows;
            normalizedAffinityMatrix =   (sqrtDegreesRows^-1) ...
                                       * normalizedAffinityMatrix ...
                                       * (sqrtDegreesColumns^-1);
end
