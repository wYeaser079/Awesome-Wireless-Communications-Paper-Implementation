function [R, throughput] = compute_SE_CF_DL(beta, gamma, pilot_index, p, eta)
%COMPUTE_SE_CF_DL Cell-Free Massive MIMO downlink achievable rate (Theorem 1, Eq. 24).
%   [R, throughput] = compute_SE_CF_DL(beta, gamma, pilot_index, p, eta)
%
%   Eq. 24:
%                        rho_d * (SUM_m sqrt(eta_mk) * gamma_mk)^2
%   SINR_k = ---------------------------------------------------------------
%             rho_d * SUM_{k'!=k, same pilot} (SUM_m sqrt(eta_mk') * gamma_mk' * beta_mk/beta_mk')^2
%           + rho_d * SUM_{k'=1}^K SUM_m eta_mk' * gamma_mk' * beta_mk
%           + 1
%
%   Inputs:
%       beta        - M x K large-scale fading coefficients
%       gamma       - M x K channel estimate variances (from compute_gamma)
%       pilot_index - K x 1 pilot assignment vector
%       p           - Parameter struct from params()
%       eta         - M x K DL power control coefficients (optional)
%                     Default: equal power, eta_mk = 1 / sum_k' gamma_mk'
%
%   Outputs:
%       R           - K x 1 spectral efficiency [bits/s/Hz]
%       throughput  - K x 1 net throughput [bits/s] (Eq. 56)

    [M, K] = size(beta);

    % Default: equal power allocation (no power control)
    if nargin < 5 || isempty(eta)
        eta = zeros(M, K);
        for m = 1:M
            eta(m, :) = 1 / sum(gamma(m, :));
        end
    end

    R = zeros(K, 1);

    for k = 1:K
        % --- Numerator: rho_d * (SUM_m sqrt(eta_mk) * gamma_mk)^2 ---
        signal = p.rho_d_cf * (sum(sqrt(eta(:,k)) .* gamma(:,k)))^2;

        % --- Denominator Term 1: Pilot contamination (coherent interference) ---
        %   SUM over k' != k with same pilot
        pilot_contam = 0;
        for kp = 1:K
            if kp == k
                continue;
            end
            if pilot_index(kp) ~= pilot_index(k)
                continue;  % Orthogonal pilots: no contribution
            end
            % SUM_m sqrt(eta_mk') * gamma_mk' * beta_mk / beta_mk'
            term = sum(sqrt(eta(:,kp)) .* gamma(:,kp) .* beta(:,k) ./ beta(:,kp));
            pilot_contam = pilot_contam + term^2;
        end
        pilot_contam = p.rho_d_cf * pilot_contam;

        % --- Denominator Term 2: Non-coherent interference ---
        %   SUM_{k'=1}^K SUM_m eta_mk' * gamma_mk' * beta_mk
        noncoh = 0;
        for kp = 1:K
            noncoh = noncoh + sum(eta(:,kp) .* gamma(:,kp) .* beta(:,k));
        end
        noncoh = p.rho_d_cf * noncoh;

        % --- SINR (Eq. 24) ---
        sinr_k = signal / (pilot_contam + noncoh + 1);

        R(k) = log2(1 + sinr_k);
    end

    % --- Net throughput (Eq. 56) ---
    throughput = p.B * p.prelog_cf * R;

end
