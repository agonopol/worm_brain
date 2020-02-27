function [assigments] = hierarchicalcompare(matrix, partitions, varargin)
    maxk = max(max(partitions));
    for i=1:2:length(varargin)-1
        if (strcmp(varargin{i}, 'maxk'))
            maxk = varargin{i+1};
        end
    end
    assigments = arrayfun(@(k) clusterdata(matrix, 'maxclust', k), ...
        1:maxk, ...
        'UniformOutput', false);
    assigments = cell2mat(assigments)';
end
