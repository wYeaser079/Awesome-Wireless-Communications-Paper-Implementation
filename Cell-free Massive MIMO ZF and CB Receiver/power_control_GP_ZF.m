function eta = power_control_GP_ZF(beta, alpha, pilot_index, p, R_bar)
%POWER_CONTROL_GP_ZF Algorithm 2: GP for ZF single-user rate max (Eq. 30-31).
%   eta = power_control_GP_ZF(beta, alpha, pilot_index, p)
%   eta = power_control_GP_ZF(beta, alpha, pilot_index, p, R_bar)
%
%   Eq. 30:  maximize    R_1^ZF
%            subject to  R_k^ZF >= R_bar_k,  k = 2,...,K   (QoS for others)
%                        0 <= eta_k <= 1                    (power)
%
%   ZF SINR (Eq. 16):
%     SINR_k = a_k * eta_k / D,   D = rho_u * SUM_i eta_i * c_i + 1
%
%   The common denominator D makes this a Geometric Program (GP).
%   In the log domain (x_k = ln eta_k):
%     ln(SINR_k) = ln(a_k) + x_k - ln(D)
%     ln(D) = log_sum_exp([ln(rho_u*c_i) + x_i ; 0])  (convex)
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
        R_bar = ones(K, 1);           % Default: 1 bit/s/Hz QoS
    end

    % --- Precompute ZF parameters (Eq. 16-20) ---

    % Eq. 17: c_i = max_m (beta_mi - alpha_mi)
    c_err = max(beta - alpha, [], 1).';       % K x 1

    % Eq. 20: dominant AP per user
    [~, dominant_ap] = max(beta, [], 1);      % 1 x K

    % Eq. 18: a_k = rho_u * Omega_k
    a = zeros(K, 1);
    for k = 1:K
        excluded = dominant_ap([1:k-1, k+1:K]);
        Mk = setdiff(1:M, excluded);          % Eq. 19
        alpha_Mk = alpha(Mk, k);
        S1 = sum(alpha_Mk);
        S2 = sum(alpha_Mk.^2);
        Omega_k = (p.N * S1^2 - S2) / S1;    % Eq. 18
        a(k) = p.rho_u * Omega_k;
    end

    % QoS: SINR_k >= gamma_k
    gamma = 2.^(R_bar / p.prelog) - 1;

    % --- Solve via CVX (log-domain GP, Eq. 31) ---
    % x_k = ln(eta_k) => power constraint: x_k <= 0

    % Constants for log_sum_exp of common denominator:
    %   ln(D) = log_sum_exp([ ln(rho_u*c_1)+x_1, ..., ln(rho_u*c_K)+x_K, 0 ])
    lse_const = log(p.rho_u * c_err);         % K x 1

    cvx_begin quiet
        variable x(K, 1)

        % Common denominator: ln(D) as log_sum_exp
        lse = log_sum_exp([lse_const + x; 0]);

        % Maximize ln(SINR_1) = ln(a_1) + x_1 - ln(D)
        maximize(log(a(1)) + x(1) - lse)

        subject to
            % QoS for users 2,...,K:
            %   ln(a_k) + x_k - ln(D) >= ln(gamma_k)
            for k = 2:K
                log(a(k)) + x(k) - lse >= log(gamma(k));
            end

            % Power: eta_k = exp(x_k) <= 1
            x <= 0;
    cvx_end

    % Extract solution
    if ~isempty(strfind(cvx_status, 'Solved'))
        eta = exp(x);
    else
        warning('GP_ZF: CVX failed (status: %s). Returning equal power.', cvx_status);
        eta = ones(K, 1);
    end

end
