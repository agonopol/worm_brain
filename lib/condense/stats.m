function [centroids, sizes] = stats(samples, assigments)
    centroids = [];
    sizes = [];
    for cluster = 1:max(assigments)
        index = find(assigments == cluster);
        data = samples(index,:);
        centroids = [centroids; mean(data', 2)'];
        sizes = [sizes, length(index)];
    end
end