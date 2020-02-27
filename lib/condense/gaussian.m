function [ M ] = gaussian( M, sigma, k )
%SIMGAUSSIAN Calculates Gaussian similarity on matrix
%   simGaussian(M, sigma) returns a matrix of the same size as
%   the distance matrix M, which contains similarity values
%   that are computed by using a Gaussian similarity function
%   with parameter sigma.
%
%   Author: Ingo Buerk
%   Year  : 2011/2012
%   Bachelor Thesis

    if nargin < 3
        k = 2;
    end
    if nargin < 2
        sigma = 1;
    end

    M = exp(-M.^k ./ (2*sigma^2));

end