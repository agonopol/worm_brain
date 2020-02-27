function [q, scores] = affmodularity(data, assigment)    

    A = CalculateAffinity(rownorm(data));
    [q, scores] = modularity(A, assigment);
end

function [affinity] = CalculateAffinity(data)
    % set the parameters
    sigma = 1;
    D = squareform(pdist(data));

    for i=1:size(data,1)    
        for j=1:size(data,1)
            dist = D(i,j);
            affinity(i,j) = exp(-dist/(2*sigma^2));
        end
    end
end

function N = rownorm(M)
    rows = size(M, 1);
    N = spdiags (sum (M,2), 0, rows, rows) \ M ;
end