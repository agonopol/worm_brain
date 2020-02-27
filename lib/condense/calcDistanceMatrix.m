function [distanceMatrix, indicatorMatrix] = calcDistanceMatrix(samples, varargin)
% CALCDISTANCEMATRIX Returns distance matrix between the samples
%   specified by the argument 'samples'. Argument 'samples' should
%   be n-by-m where n is the number of samples and m the number of
%   dimensions of the space the samples are taken from. The return
%   value 'distanceMatrix' returns the calculated distances. The
%   return value 'indicatorMatrix' contains the values where the
%   distance between points is zero. This is required to differentiate
%   between real zero distance and the implicit infinity distance
%   that is represent as zero in the sparse matrix if k is selected
%   finite.
%
%   Author: Tobias Welp, modifying a version by Ingo Buerk
%   Year  : 2016 based on original work from 2011/2012

    % default arguments
    k_knn = size(samples, 1);
    distfun = 'euclidean';
    type_k_knn = 'normal';
    verbosityLevel = 1;
    indentationLevel = 0;
    mode = 'normal';
    lengthPartitions = 0.2;

    % parsing of variable argument list
    for i=1:length(varargin)-1
        if (strcmp(varargin{i}, 'k_knn'))
            k_knn = varargin{i+1};
            if (k_knn == 'full')
                k_knn = size(samples, 1);
            end
        end
        if (strcmp(varargin{i}, 'mode'))
            mode = varargin{i+1};
        end
        if (strcmp(varargin{i}, 'lengthPartitions'))
            lengthPartitions = varargin{i+1};
        end
        if (strcmp(varargin{i}, 'n_pca'))
            n_pca = varargin{i+1};
        end
        if (strcmp(varargin{i}, 'distfun'))
            distfun = varargin{i+1};
        end
        if (strcmp(varargin{i}, 'type_k_knn'))
            type_k_knn = varargin{i+1};
            if ~ismember(type_k_knn, {'normal', 'mutual'})
                warning(['calcDistanceMatrix: Invalid choice for argument' ...
                         'type_k_knn. Using default value "normal"']);
                type_k_knn = 'normal';
            end
        end
        if (strcmp(varargin{i}, 'verbosityLevel'))
            verbosityLevel = varargin{i+1};
        end
        if (strcmp(varargin{i}, 'indentationLevel'))
            indentationLevel = varargin{i+1};
        end
    end

    if exist('n_pca', 'var')
        if (n_pca < size(samples, 2))
            if (verbosityLevel > 0)
                disp([indent(indentationLevel) 'CalcDistanceMatrix: PCA to ' num2str(n_pca) ' dimensions.'])
            end
            samples = real(pcaMaaten(samples, n_pca));
            % Remove potentially present NaN columns
            indicesToRemove = any(isnan(samples));
            samples = samples(:, ~indicesToRemove);
        end
    end

    n = size(samples, 1);
    indicatorMatrix = sparse(1:n, 1:n, ones(1, n));

    switch (mode)
        case 'normal'
            numDistancesConsidered = min(n, k_knn);

            % Preallocate memory
            indi = zeros(1, numDistancesConsidered * n);
            indj = zeros(1, numDistancesConsidered * n);
            inds = zeros(1, numDistancesConsidered * n);

            for ii = 1:n
                % Compute i-th column of distance matrix
                dist = pdist2(samples(ii, :), samples, distfun);

                % Sort row by distance
                [s, O] = sort(dist, 'ascend');

                % Save indices and value of the k
                indi(1, (ii-1)*numDistancesConsidered+1:ii*numDistancesConsidered) = ii;
                indj(1, (ii-1)*numDistancesConsidered+1:ii*numDistancesConsidered) = O(1:numDistancesConsidered);
                inds(1, (ii-1)*numDistancesConsidered+1:ii*numDistancesConsidered) = s(1:numDistancesConsidered);

                id = find(s == 0);
                if ~isempty(id)
                    indicatorMatrix(ii, O(id)) = 1;
                end
            end
        case 'partitioning'
            %% populate partions
            d = size(samples, 2);
            cell2Partition = {};
            partition2Cells = {};
            for j=1:d
                cell2PartitionD = containers.Map('KeyType', 'double', 'ValueType', 'double');
                partition2CellsD = containers.Map('KeyType', 'double', 'ValueType', 'any');
                for i=1:n
                    partition = floor(samples(i, j)/lengthPartitions);
                    cell2PartitionD(i) = partition;
                    if (partition2CellsD.isKey(partition))
                        partition2CellsD(partition) = [partition2CellsD(partition) i];
                    else
                        partition2CellsD(partition) = [i];
                    end
                end
                cell2Partition{j} = cell2PartitionD;
                partition2Cells{j} = partition2CellsD;
            end
            %for j = 1:d
            %    for i = 1:n
            %        disp(['Dim ' num2str(j) ': Partition of cell ' num2str(i) ': ' num2str(cell2Partition{j}(i))])
            %    end
            %end
            %% calc distances to cells in same and neighbor partitions
            indi = [];
            indj = [];
            inds = [];
            for i=1:n
                candidates = [];
                for j=1:d
                    partition = cell2Partition{j}(i);
                    candidatesD = partition2Cells{j}(partition);
                    if (partition2Cells{j}.isKey(partition-1))
                        candidatesD = [candidatesD partition2Cells{j}(partition-1)];
                    end
                    if (partition2Cells{j}.isKey(partition+1))
                        candidatesD = [candidatesD partition2Cells{j}(partition+1)];
                    end
                    if (j==1)
                        candidates = candidatesD;
                    else
                        candidates = intersect(candidates, candidatesD);
                    end
                end
                numDistancesConsidered = min(length(candidates), k_knn);
                dist = pdist2(samples(i, :), samples(candidates, :), distfun);
                [s, O] = sort(dist, 'ascend');
                indi = [indi repmat(i, 1, numDistancesConsidered)];
                indj = [indj candidates(O(1:numDistancesConsidered))];
                inds = [inds s(1:numDistancesConsidered)];

                id = find(s == 0);
                if ~isempty(id)
                    indicatorMatrix(i, candidates(O(id))) = 1;
                end

            end
    end
    % Create sparse matrix
    distanceMatrix = sparse(indi, indj, inds, n, n);

    %Construct either normal or mutual graph
    switch (type_k_knn)
        case 'normal'
            distanceMatrix = max(distanceMatrix, distanceMatrix');
        case 'mutual'
            distanceMatrix = min(distanceMatrix, distanceMatrix');
    end
end
