function eta = power_control_GP_CB(beta, alpha, pilot_index, p, R_bar)
%POWER_CONTROL_GP_CB Algorithm 4: GP for CB single-user rate max (Eq. 44-45).
%   eta = power_control_GP_CB(beta, alpha, pilot_index, p)
%   eta = power_control_GP_CB(beta, alpha, pilot_index, p, R_bar)
%
%   Eq. 44:  maximize    R_1^CB
%            subject to  R_k^CB >= R_bar_k,  k = 2,...,K   (QoS for others)
%                        0 <= eta_k <= 1                    (power)
%
%   CB SINR (Lemma 1, Eq. 33):
%     SINR_k = eta_k * S_k / D_k
%     S_k    = ||l_kk||^2 = SUM_m alpha_mk^2               (signal)
%     D_k    = SUM_j eta_j * w_{jk} + noise_k               (user-specific denom)
%
%   In the log domain (x_k = ln eta_k):
%     ln(SINR_k) = ln(S_k) + x_k - ln(D_k)
%     ln(D_k) = log_sum_exp([ln(w_{jk}) + x_j ; ln(noise_k)])  (convex)
%
%   So max ln(SINR_1) is a concave maximization — directly solvable by CVX.
%
%   Requires: CVX (http://cvxr.com/cvx/)
%
%   Inputs:
%       beta, alpha, pilot_index, p - Standard parameters
%       R_bar  - K x 1 QoS targets [bits/s/Hz] (default = 1 for all)
%
%   Output:
%       eta    - K x 1 optimized power coefficients

    [M, K] = size(beta);

    if nargin < 5
        R_bar = ones(K, 1);
    end

    % --- Precompute CB parameters (Lemma 1, Eq. 33-35) ---

    % Signal: S_k = ||l_kk||^2 (Eq. 34 with i = k)
    S = sum(alpha.^2, 1).';                    % K x 1

    % Interference weight matrix w(i,k)
    w = zeros(K, K);
    for k = 1:K
        for i = 1:K
            % Eq. 35: non-coherent term (always present)
            noncoh = sum(beta(:, i) .* alpha(:, k)) / p.N;

            % Eq. 34: coherent term (only co-pilot, i != k)
            coh = 0;
            if i ~= k && pilot_index(i) == pilot_index(k)
                l_ik = alpha(:, k) .* beta(:, i) ./ beta(:, k);
                coh = sum(l_ik.^2);
            end

            w(i, k) = coh + noncoh;
        end
    end

    % Noise: noise_k = ||l_kk||_1 / (N * rho_u)
    noise = sum(alpha, 1).' / (p.N * p.rho_u);    % K x 1

    % QoS: SINR_k >= gamma_k
    gamma = 2.^(R_bar / p.prelog) - 1;

    % --- Solve via CVX (log-domain GP, Eq. 45) ---
    cvx_begin quiet
        variable x(K, 1)

        % User 1's denominator:
        %   ln(D_1) = log_sum_exp([ln(w_{j1}) + x_j ; ln(noise_1)])
        lse1 = log_sum_exp([log(w(:,1)) + x; log(noise(1))]);

        % Maximize ln(SINR_1) = ln(S_1) + x_1 - ln(D_1)
        maximize(log(S(1)) + x(1) - lse1)

        subject to
            % QoS for users 2,...,K
            for k = 2:K
                lse_k = log_sum_exp([log(w(:,k)) + x; log(noise(k))]);
                log(S(k)) + x(k) - lse_k >= log(gamma(k));
            end

            % Power: eta_k = exp(x_k) <= 1
            x <= 0;
    cvx_end

    % Extract solution
    if ~isempty(strfind(cvx_status, 'Solved'))
        eta = exp(x);
    else
        warning('GP_CB: CVX failed (status: %s). Returning equal power.', cvx_status);
        eta = ones(K, 1);
    end

end
