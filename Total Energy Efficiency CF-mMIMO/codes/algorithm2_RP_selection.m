function [gamma_hat, A_k, num_APs_per_user] = algorithm2_RP_selection( ...
                                                eta_mk, beta_mk, gamma_mk, rho_pct)
%ALGORITHM2_RP_SELECTION  Received-power-based AP selection of
%  Ngo-Tran-Duong-Matthaiou-Larsson IEEE TGCN 2018, Section V-A.
%
%  For every user k, rank the APs by the fraction of the desired-signal
%  power they contribute, p(m,k) = gamma_{m,k}*eta_{m,k}/sum_m' ...,
%  and keep the smallest set whose cumulative p(.,k) reaches rho%.
%
%  INPUTS:
%     eta_mk   - M-by-K optimal power control from a PRIOR run of
%                Algorithm 1 on the un-selected problem.
%     beta_mk  - M-by-K large-scale fading (only used to preserve shape)
%     gamma_mk - M-by-K MMSE variance
%     rho_pct  - percentile threshold (e.g. 95)
%
%  OUTPUTS:
%     gamma_hat         - M-by-K modified variance with zeros on removed
%                         (m,k) pairs
%     A_k               - 1xK cell of serving AP index lists
%     num_APs_per_user  - K-by-1 vector

    [M, K] = size(beta_mk);
    gamma_hat = zeros(M, K);
    A_k       = cell(1, K);
    num_APs_per_user = zeros(K, 1);

    p_mk = gamma_mk .* eta_mk;                            % p-numerator
    colsum = sum(p_mk, 1);                                % 1 x K
    colsum(colsum == 0) = 1;                              % avoid /0
    p_mk = bsxfun(@rdivide, p_mk, colsum);

    target = rho_pct / 100;
    for k = 1:K
        [ps, idx] = sort(p_mk(:, k), 'descend');
        csum = cumsum(ps);
        Mk = find(csum >= target, 1, 'first');
        if isempty(Mk), Mk = M; end
        chosen = idx(1:Mk);
        A_k{k} = chosen;
        num_APs_per_user(k) = Mk;
        gamma_hat(chosen, k) = gamma_mk(chosen, k);
    end
end
