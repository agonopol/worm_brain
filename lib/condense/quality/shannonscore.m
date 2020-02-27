function [ Q ] = shannonscore( data, assigments )
     
     Q = zeros(size(data, 1), 1);
     for i = 1:max(assigments)
         C = data(assigments == i, :);
         if ( size(C, 1) == 1 ) 
              Q(assigments == i) = 0.0;
         else
            A = CalculateAffinity(C);
            Q(assigments == i) = shannon(eig(A) / sum(eig(A)));
         end

     end
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