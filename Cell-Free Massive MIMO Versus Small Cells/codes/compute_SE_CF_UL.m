function [R, throughput] = compute_SE_CF_UL(beta, gamma, pilot_index, p, eta)
%COMPUTE_SE_CF_UL Cell-Free Massive MIMO uplink achievable rate (Theorem 2, Eq. 27).
%   [R, throughput] = compute_SE_CF_UL(beta, gamma, pilot_index, p, eta)
%
%   Eq. 27:
%                        rho_u * eta_k * (SUM_m gamma_mk)^2
%   SINR_k = ---------------------------------------------------------------
%             rho_u * SUM_{k'!=k, same pilot} eta_k' * (SUM_m gamma_mk * beta_mk'/beta_mk)^2
%           + rho_u * SUM_{k'=1}^K eta_k' * SUM_m gamma_mk * beta_mk'
%           + SUM_m gamma_mk
%
%   Inputs:
%       beta        - M x K large-scale fading coefficients
%       gamma       - M x K channel estimate variances (from compute_gamma)
%       pilot_index - K x 1 pilot assignment vector
%       p           - Parameter struct from params()
%       eta         - K x 1 UL power control coefficients (optional)
%                     Default: eta_k = 1 for all k (full power)
%
%   Outputs:
%       R           - K x 1 spectral efficiency [bits/s/Hz]
%       throughput  - K x 1 net throughput [bits/s] (Eq. 56)

    [M, K] = size(beta);

    % Default: full power (no power control)
    if nargin < 5 || isempty(eta)
        eta = ones(K, 1);
    end

    R = zeros(K, 1);

    for k = 1:K
        % --- Numerator: rho_u * eta_k * (SUM_m gamma_mk)^2 ---
        signal = p.rho_u_cf * eta(k) * (sum(gamma(:,k)))^2;

        % --- Denominator Term 1: Pilot contamination (coherent interference) ---
        %   SUM over k' != k with same pilot
        pilot_contam = 0;
        for kp = 1:K
            if kp == k
                continue;
            end
            if pilot_index(kp) ~= pilot_index(k)
                continue;
            end
            % eta_k' * (SUM_m gamma_mk * beta_mk'/beta_mk)^2
            term = sum(gamma(:,k) .* beta(:,kp) ./ beta(:,k));
            pilot_contam = pilot_contam + eta(kp) * term^2;
        end
        pilot_contam = p.rho_u_cf * pilot_contam;

        % --- Denominator Term 2: Non-coherent interference ---
        %   SUM_{k'=1}^K eta_k' * SUM_m gamma_mk * beta_mk'
        noncoh = 0;
        for kp = 1:K
            noncoh = noncoh + eta(kp) * sum(gamma(:,k) .* beta(:,kp));
        end
        noncoh = p.rho_u_cf * noncoh;

        % --- Denominator Term 3: Noise ---
        %   SUM_m gamma_mk
        noise = sum(gamma(:,k));

        % --- SINR (Eq. 27) ---
        sinr_k = signal / (pilot_contam + noncoh + noise);

        R(k) = log2(1 + sinr_k);
    end

    % --- Net throughput (Eq. 56) ---
    throughput = p.B * p.prelog_cf * R;

end
