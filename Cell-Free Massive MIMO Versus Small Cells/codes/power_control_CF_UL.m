function eta = power_control_CF_UL(beta, gamma, pilot_index, p)
%POWER_CONTROL_CF_UL Max-min UL power control for Cell-Free (Eq. 36-37).
%   eta = power_control_CF_UL(beta, gamma, pilot_index, p)
%
%   Solves the quasi-linear max-min problem (Eq. 36):
%     max_{eta_k}  min_k  R_u,k^cf
%     s.t.  0 <= eta_k <= 1
%
%   For fixed target SINR t, the feasibility constraints are linear in eta_k.
%   Uses bisection on t, with linprog at each step.
%
%   Inputs:
%       beta        - M x K large-scale fading coefficients
%       gamma       - M x K channel estimate variances
%       pilot_index - K x 1 pilot assignment vector
%       p           - Parameter struct from params()
%
%   Output:
%       eta         - K x 1 UL power control coefficients

    [M, K] = size(beta);

    % --- Precompute interference terms ---
    % Signal strength for each user: S_k = (SUM_m gamma_mk)^2
    S = zeros(K, 1);
    for k = 1:K
        S(k) = (sum(gamma(:,k)))^2;
    end

    % Pilot contamination: PC(k,k') = (SUM_m gamma_mk * beta_mk'/beta_mk)^2 * xi(k,k')
    PC = zeros(K, K);
    for k = 1:K
        for kp = 1:K
            if kp ~= k && pilot_index(k) == pilot_index(kp)
                PC(k, kp) = (sum(gamma(:,k) .* beta(:,kp) ./ beta(:,k)))^2;
            end
        end
    end

    % Non-coherent interference: NC(k,k') = SUM_m gamma_mk * beta_mk'
    NC = zeros(K, K);
    for k = 1:K
        for kp = 1:K
            NC(k, kp) = sum(gamma(:,k) .* beta(:,kp));
        end
    end

    % Noise term: N_k = SUM_m gamma_mk
    N = zeros(K, 1);
    for k = 1:K
        N(k) = sum(gamma(:,k));
    end

    % --- Bisection ---
    t_min = 0;

    % Upper bound: compute with equal power
    R_eq = compute_SE_CF_UL(beta, gamma, pilot_index, p, ones(K,1));
    t_max = 2^(max(R_eq)) - 1;
    t_max = max(t_max, 1);

    eta_best = ones(K, 1);

    for iter = 1:p.bisection_maxiter
        t = (t_min + t_max) / 2;

        % Build LP: for each k, SINR_k >= t
        % rho_u * eta_k * S_k >= t * [rho_u * SUM_{k'!=k} eta_k'*(PC_kk'+NC_kk')
        %                              + rho_u * eta_k * NC_kk + N_k]
        %
        % Rearranging: A * eta >= b
        A = zeros(K, K);
        b = zeros(K, 1);

        for k = 1:K
            for kp = 1:K
                if kp == k
                    A(k, k) = p.rho_u_cf * S(k) - t * p.rho_u_cf * NC(k, k);
                else
                    A(k, kp) = -t * p.rho_u_cf * (PC(k, kp) + NC(k, kp));
                end
            end
            b(k) = t * N(k);
        end

        % LP feasibility: find eta s.t. A*eta >= b, 0 <= eta <= 1
        % linprog minimizes f'*eta s.t. A_ineq*eta <= b_ineq, lb <= eta <= ub
        % Convert A*eta >= b to -A*eta <= -b
        options = optimoptions('linprog', 'Display', 'off');
        f = zeros(K, 1);  % Feasibility check: any objective
        [eta_trial, ~, exitflag] = linprog(f, -A, -b, [], [], ...
                                           zeros(K,1), ones(K,1), options);

        if exitflag == 1
            t_min = t;
            eta_best = eta_trial;
        else
            t_max = t;
        end

        if t_max - t_min < p.bisection_tol
            break;
        end
    end

    eta = eta_best;

end
