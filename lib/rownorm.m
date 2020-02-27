function N = rownorm(M)
    rows = size(M, 1);
    N = spdiags (sum (M,2), 0, rows, rows) \ M ;
end