function [neurons, matrix] = weightedadj(data)
    neurons = union(unique(data.Neuron1), unique(data.Neuron2));
    neurons = cellfun(@flip, sort(cellfun(@flip, cellstr(neurons), 'UniformOutput', false)), 'UniformOutput', false);
    neurons = categorical(neurons);
    sections = unique(data.EMSection);
    matrix = zeros(length(neurons), length(neurons));
    for i = 1:length(sections)
        section = data(data.EMSection == sections(i), {'Neuron1', 'Neuron2', 'Weight'});
        for n = 1:height(section)
            n1 = find(neurons == section(n, :).Neuron1);
            n2 = find(neurons == section(n, :).Neuron2);
            weight = section(n, :).Weight;
            if isnan(weight)
                matrix(n1,n2) = matrix(n1,n2) + 1;
                matrix(n2,n1) = matrix(n2,n1) + 1;
            else
                matrix(n1,n2) = matrix(n1,n2) + weight;
                matrix(n2,n1) = matrix(n2,n1) + weight;
            end
        end
    end
end

