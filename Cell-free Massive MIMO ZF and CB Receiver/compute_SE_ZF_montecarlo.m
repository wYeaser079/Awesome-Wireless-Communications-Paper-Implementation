function R_ZF_mc = compute_SE_ZF_montecarlo(beta, alpha, c_coeff, pilot_index, p, eta)
%COMPUTE_SE_ZF_MONTECARLO ZF spectral efficiency via Monte Carlo (Eq. 15).
%   R_ZF_mc = compute_SE_ZF_montecarlo(beta, alpha, c_coeff, pilot_index, p)
%   R_ZF_mc = compute_SE_ZF_montecarlo(beta, alpha, c_coeff, pilot_index, p, eta)
%
%   For each channel realization:
%     1. Generate G, G_hat via generate_channel          (Eq. 1-4)
%     2. Build ZF detector A = G_hat*(G_hat^H*G_hat)^-1  (Section 2.2)
%     3. Compute exact SINR_k from Eq. 15
%     4. R_k = log2(1 + SINR_k)
%   Then average over realizations and apply pre-log.
%
%   Eq. 15:
%     SINR_k^ZF = rho_u * eta_k
%       / SUM_m ||a_mk||^2 * [rho_u * SUM_i eta_i*(beta_mi - alpha_mi) + 1]
%
%   Purpose: Validates Theorem 1 (Eq. 16) — reproduces Figure 1 of the paper.
%
%   Inputs:
%       beta        - M x K large-scale fading coefficients
%       alpha       - M x K channel estimate variances (Eq. 8)
%       c_coeff     - M x K MMSE filter coefficients   (Eq. 5)
%       pilot_index - K x 1 pilot assignment vector
%       p           - Parameter struct from params()
%       eta         - K x 1 power control coefficients (optional, default = 1)
%
%   Output:
%       R_ZF_mc     - K x 1 per-user spectral efficiency [bits/s/Hz]

    [M, K] = size(beta);
    N = p.N;

    if nargin < 6
        eta = ones(K, 1);
    end

    % Pre-compute per-AP weighted estimation error (does not change across realizations)
    %   err_m(m) = SUM_i eta_i * (beta_mi - alpha_mi)
    err_m = (beta - alpha) * eta;   % M x 1

    R_sum = zeros(K, 1);

    for r = 1:p.num_channel_real
        % --- Generate one channel realization (Eq. 1-4) ---
        [~, G_hat] = generate_channel(beta, alpha, c_coeff, pilot_index, p);

        % --- Build ZF detector (Section 2.2) ---
        %   A = G_hat * (G_hat^H * G_hat)^{-1}   (MN x K)
        A = G_hat / (G_hat' * G_hat);

        % --- Compute SINR for each user (Eq. 15) ---
        for k = 1:K
            a_k = A(:, k);                          % MN x 1

            % Reshape into N x M: each column is the N-element sub-vector a_mk
            a_k_mat = reshape(a_k, N, M);           % N x M

            % Per-AP norms: ||a_mk||^2 for m = 1,...,M
            a_mk_sq = sum(abs(a_k_mat).^2, 1).';   % M x 1

            % Denominator of Eq. 15:
            %   SUM_m ||a_mk||^2 * [rho_u * SUM_i eta_i*(beta_mi - alpha_mi) + 1]
            denom = sum(a_mk_sq .* (p.rho_u * err_m + 1));

            % SINR (Eq. 15)
            sinr_k = p.rho_u * eta(k) / denom;

            % Accumulate log-rate for averaging
            R_sum(k) = R_sum(k) + log2(1 + sinr_k);
        end
    end

    % Average over realizations and apply pre-log factor
    R_ZF_mc = p.prelog * R_sum / p.num_channel_real;

end
