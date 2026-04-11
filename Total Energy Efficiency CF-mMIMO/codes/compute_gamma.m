function gamma_mk = compute_gamma(beta_mk, pilotseq, tau_p, rho_p)
%COMPUTE_GAMMA  Evaluate the MMSE per-component variance gamma_{m,k} of
%  the channel estimate, given the large-scale fading beta_{m,k} and the
%  pilot assignment.  Implements Eq. (5) of Ngo et al. 2018:
%
%      gamma_{m,k} = tau_p * rho_p * beta_{m,k}^2 /
%                     ( tau_p * rho_p * sum_{k'} beta_{m,k'} |phi_{k'}^H phi_k|^2 + 1 )
%
%  INPUTS:
%     beta_mk  - M-by-K large-scale fading matrix (linear)
%     pilotseq - tau_p-by-K pilot matrix (unit-norm columns)
%     tau_p    - pilot length (samples)
%     rho_p    - normalised UL pilot SNR
%
%  OUTPUT:
%     gamma_mk - M-by-K MMSE per-component variance

    [M, K] = size(beta_mk);

    % |phi_{k'}^H phi_k|^2 for all pairs  (K x K, real-valued)
    inner = abs(pilotseq' * pilotseq).^2;       % K x K

    % For each m,k: sum_{k'} beta_{m,k'} |phi_{k'}^H phi_k|^2
    % = (beta_mk * inner)(m,k)     (M-by-K)
    den = beta_mk * inner;

    gamma_mk = tau_p * rho_p * (beta_mk.^2) ./ (tau_p * rho_p * den + 1);
    % Numerical safety
    gamma_mk = max(gamma_mk, 0);

    if any(~isfinite(gamma_mk(:)))
        error('compute_gamma:NaN', 'Non-finite gamma entries produced.');
    end

    % Force consistency (M,K) output dimensions
    if ~isequal(size(gamma_mk), [M, K])
        gamma_mk = reshape(gamma_mk, M, K);
    end
end
