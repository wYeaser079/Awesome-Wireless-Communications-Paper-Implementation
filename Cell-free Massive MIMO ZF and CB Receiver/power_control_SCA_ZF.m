function eta = power_control_SCA_ZF(beta, alpha, pilot_index, p, R_bar)
%POWER_CONTROL_SCA_ZF Algorithm 1: SCA for ZF total rate maximization (Eq. 21-29).
%   eta = power_control_SCA_ZF(beta, alpha, pilot_index, p)
%   eta = power_control_SCA_ZF(beta, alpha, pilot_index, p, R_bar)
%
%   Eq. 21:  maximize    SUM_k R_k^ZF({eta_k})
%            subject to  R_k^ZF >= R_bar_k,  for all k   (QoS constraint)
%                        0 <= eta_k <= 1                  (power constraint)
%
%   where R_k^ZF = prelog * log2(1 + a_k * eta_k / D)  (Theorem 1, Eq. 16)
%         a_k = rho_u * Omega_k                          (Eq. 18)
%         D   = rho_u * SUM_i eta_i * c_i + 1            (common denominator)
%         c_i = max_m (beta_mi - alpha_mi)                (Eq. 17)
%
%   DC programming approach:
%     SUM_k R_k = SUM_k log(a_k*eta_k + D) - K*log(D)
%               = f(eta) - g(eta)        (concave - concave)
%
%   At each SCA iteration n, linearize g around eta^(n):
%     g(eta) ~ g(eta^(n)) + nabla_g^T * (eta - eta^(n))
%     nabla_g_j = K * rho_u * c_j / D^(n)
%
%   The resulting subproblem is concave maximization (solvable by CVX).
%   Convergence: monotonically increasing, bounded above (~4 iterations).
%
%   Requires: CVX (http://cvxr.com/cvx/)
%
%   Inputs:
%       beta        - M x K large-scale fading coefficients
%       alpha       - M x K channel estimate variances (Eq. 8)
%       pilot_index - K x 1 pilot assignment vector
%       p           - Parameter struct from params()
%       R_bar       - K x 1 QoS targets [bits/s/Hz] (default = 0, no QoS)
%
%   Output:
%       eta         - K x 1 optimized power control coefficients

    [M, K] = size(beta);

    if nargin < 5
        R_bar = zeros(K, 1);
    end

    % --- Precompute ZF closed-form parameters (Eq. 16-20) ---

    % Eq. 17: c_i = max_m (beta_mi - alpha_mi) for each user i
    c_err = max(beta - alpha, [], 1).';       % K x 1

    % Eq. 20: m_j* = argmax_m beta_mj (dominant AP for user j)
    [~, dominant_ap] = max(beta, [], 1);      % 1 x K

    % Eq. 18: a_k = rho_u * Omega_k for each user k
    a = zeros(K, 1);
    for k = 1:K
        excluded = dominant_ap([1:k-1, k+1:K]);
        Mk = setdiff(1:M, excluded);          % Eq. 19: M_k
        alpha_Mk = alpha(Mk, k);
        S1 = sum(alpha_Mk);
        S2 = sum(alpha_Mk.^2);
        Omega_k = (p.N * S1^2 - S2) / S1;    % Eq. 18
        a(k) = p.rho_u * Omega_k;
    end

    % QoS thresholds: SINR_k >= gamma_k
    gamma = 2.^(R_bar / p.prelog) - 1;

    % --- SCA iterations (Algorithm 1) ---
    eta      = ones(K, 1);       % Initial point: equal power control
    max_iter = 20;
    epsilon  = 1e-4;

    for iter = 1:max_iter
        % Current denominator D^(n) = rho_u * SUM_i eta_i * c_i + 1
        D_n = p.rho_u * (c_err' * eta) + 1;

        % Gradient of g(eta) = K*ln(D) at eta^(n):
        %   nabla_j = K * rho_u * c_j / D_n
        d = K * p.rho_u * c_err / D_n;         % K x 1

        % --- Solve convex subproblem via CVX ---
        cvx_begin quiet
            variable eta_new(K, 1)

            % Objective: f(eta) - nabla_g' * eta  (concave + affine = concave)
            obj = -d' * eta_new;
            for k = 1:K
                % f_k = log(a_k*eta_k + D) where D = rho_u*c_err'*eta + 1
                obj = obj + log(a(k) * eta_new(k) ...
                            + p.rho_u * (c_err' * eta_new) + 1);
            end

            maximize(obj)

            subject to
                % QoS: a_k*eta_k >= gamma_k * D  (linear in eta)
                %   => a_k*eta_k - gamma_k*rho_u*(c_err'*eta) >= gamma_k
                for k = 1:K
                    if gamma(k) > 0
                        a(k) * eta_new(k) ...
                            - gamma(k) * p.rho_u * (c_err' * eta_new) >= gamma(k);
                    end
                end

                % Power constraints
                eta_new >= 0;
                eta_new <= 1;
        cvx_end

        % Check CVX status
        if isempty(strfind(cvx_status, 'Solved'))
            warning('SCA_ZF: CVX failed at iteration %d (status: %s).', ...
                    iter, cvx_status);
            break;
        end

        % Check convergence
        if norm(eta_new - eta) < epsilon
            eta = eta_new;
            break;
        end
        eta = eta_new;
    end

end
