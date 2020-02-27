function [zscored] = zscorep(data, p)
    zscored = zscore(data);
    z = -sqrt(2) * erfcinv(p*2);
    zscored ( zscored >= z ) = z;
    zscored ( zscored <= -z ) = -z;
end