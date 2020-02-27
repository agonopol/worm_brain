function product = weightedMultiply(lHS, rHS, weights)
    product = lHS * diag(weights) * rHS;
end
