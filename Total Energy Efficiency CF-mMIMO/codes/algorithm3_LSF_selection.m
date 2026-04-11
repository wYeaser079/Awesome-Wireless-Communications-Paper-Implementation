function [gamma_hat, A_k, num_APs_per_user] = algorithm3_LSF_selection( ...
                                                beta_mk, gamma_mk, rho_pct)
%ALGORITHM3_LSF_SELECTION  Largest-Large-Scale-Fading based AP selection
%  of Ngo-Tran-Duong-Matthaiou-Larsson IEEE TGCN 2018, Section V-B.
%
%  For every user k, rank the APs by their large-scale fading
%  beta_{m,k} in descending order and keep the smallest set whose
%  cumulative sum reaches rho% of the total Σ_m beta_{m,k}.
%
%  INPUTS:
%     beta_mk  - M-by-K large-scale fading (linear)
%     gamma_mk - M-by-K MMSE variance
%     rho_pct  - percentile threshold (e.g. 95)
%
%  OUTPUTS:
%     gamma_hat         - M-by-K modified variance: gamma_hat(m,k) =
%                         gamma(m,k) if AP m serves user k, 0 otherwise.
%                         Plug this directly into Algorithm 1 in place
%                         of gamma.
%     A_k               - 1xK cell, A_k{k} lists the AP indices serving k
%     num_APs_per_user  - K-by-1 vector of |A_k|

    [M, K] = size(beta_mk);
    gamma_hat = zeros(M, K);
    A_k       = cell(1, K);
    num_APs_per_user = zeros(K, 1);

    target = rho_pct / 100;
    for k = 1:K
        bk = beta_mk(:, k);
        [bs, idx] = sort(bk, 'descend');
        csum = cumsum(bs) / sum(bk);
        Mk = find(csum >= target, 1, 'first');
        if isempty(Mk), Mk = M; end
        chosen = idx(1:Mk);
        A_k{k} = chosen;
        num_APs_per_user(k) = Mk;
        gamma_hat(chosen, k) = gamma_mk(chosen, k);
    end
end
