function l = label(labels, scores)
    l = arrayfun(@(i) sprintf('%s: %0.2f', ...
        labels{i}, scores(i)), 1:size(scores,1), 'UniformOutput', false);    
end