function R_ZF = compute_SE_ZF(beta, alpha, pilot_index, p, eta)
%COMPUTE_SE_ZF ZF spectral efficiency via closed-form (Theorem 1, Eq. 16-20).
%   R_ZF = compute_SE_ZF(beta, alpha, pilot_index, p)
%   R_ZF = compute_SE_ZF(beta, alpha, pilot_index, p, eta)
%
%   Eq. 16:  R_k^ZF ~ log2(1 + rho_u * Omega(k) * eta_k
%                              / [rho_u * SUM_i eta_i * c_i + 1])
%
%   Eq. 17:  c_i = max_m (beta_mi - alpha_mi)
%   Eq. 18:  Omega(k) = (N * S1^2 - S2) / S1
%                where S1 = SUM_{m in M_k} alpha_mk
%                      S2 = SUM_{m in M_k} alpha_mk^2
%   Eq. 19:  M_k = {1,...,M} \ {m_j* : j != k}
%   Eq. 20:  m_j* = argmax_m beta_mj
%
%   Inputs:
%       beta        - M x K large-scale fading coefficients
%       alpha       - M x K channel estimate variances (Eq. 8)
%       pilot_index - K x 1 pilot assignment vector
%       p           - Parameter struct from params()
%       eta         - K x 1 power control coefficients (optional, default = 1)
%
%   Output:
%       R_ZF        - K x 1 per-user spectral efficiency [bits/s/Hz]

    [M, K] = size(beta);

    % Default: equal power control
    if nargin < 5
        eta = ones(K, 1);
    end

    % --- Eq. 17: c_i = max_m (beta_mi - alpha_mi) for each user i ---
    c_err = max(beta - alpha, [], 1).';   % K x 1

    % --- Eq. 20: m_j* = argmax_m beta_mj for each user j ---
    [~, dominant_ap] = max(beta, [], 1);  % 1 x K row vector

    % --- Denominator of Eq. 16 (shared across all users) ---
    denom = p.rho_u * sum(eta .* c_err) + 1;

    R_ZF = zeros(K, 1);

    for k = 1:K
        % --- Eq. 19: M_k = {1,...,M} \ {m_j* : j != k} ---
        excluded = dominant_ap([1:k-1, k+1:K]);
        Mk = setdiff(1:M, excluded);

        % --- Eq. 18: Omega(k) ---
        alpha_Mk = alpha(Mk, k);            % alpha_mk for m in M_k
        S1 = sum(alpha_Mk);                  % SUM_{m in M_k} alpha_mk
        S2 = sum(alpha_Mk.^2);               % SUM_{m in M_k} alpha_mk^2
        Omega_k = (p.N * S1^2 - S2) / S1;

        % --- Eq. 16: SINR ---
        sinr_k = p.rho_u * Omega_k * eta(k) / denom;

        % --- Spectral efficiency with pre-log factor ---
        R_ZF(k) = p.prelog * log2(1 + sinr_k);
    end

end
