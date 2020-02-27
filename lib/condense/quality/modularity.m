function [q, scores] = modularity(adj, assigment)

    deg = sum(adj);
    m = sum(deg)/2;
    assigment = assigment(:);

    A = modmat(adj);
    Q = zeros(1, max(unique(assigment)));
    scores = zeros(size(assigment));

    for v = 1:max(unique(assigment))
        s = ones(size(assigment));
        s(assigment == v) = 1;
        s(assigment ~= v) = - 1;
        Q(v) = ( s' * A * s / (4*m) );
        scores( assigment == v ) = Q(v);
    end

    % q = assigment' * modmat(adj) * assigment / (4*m);
    q = sum(Q);
    scores = scores';
end

function [b] = modmat(adj)
    n = length(adj);
    deg = sum(adj,1);
    m = sum(deg)/2;
    b = zeros(size(adj));
    for k1 = 1:n
        for k2 = 1:n
            b(k1, k2) = adj(k1, k2) - deg(k1) * deg(k2) / (2 * m);
        end
    end
end