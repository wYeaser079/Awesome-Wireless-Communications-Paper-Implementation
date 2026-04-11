function eta_mk = equal_power_allocation(gamma_mk, N)
%EQUAL_POWER_ALLOCATION  Baseline "no power control" (Fig. 2 of Ngo et
%  al. 2018).  Each AP uses its full transmit budget 1/N and splits it
%  equally among all users:
%
%      eta_{m,k} = 1 / ( N * sum_{k'} gamma_{m,k'} )     for all k
%
%  so that sum_k eta_{m,k} gamma_{m,k} = 1/N exactly (constraint (8)
%  holds with equality).

    denom = N * sum(gamma_mk, 2);              % M x 1
    eta_mk = bsxfun(@rdivide, ones(size(gamma_mk)), denom);
end
