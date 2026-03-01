function eta = power_control_CF_DL(beta, gamma, pilot_index, p)
%POWER_CONTROL_CF_DL Max-min DL power control for Cell-Free (Algorithm 2).
%   eta = power_control_CF_DL(beta, gamma, pilot_index, p)
%
%   Solves the max-min fairness problem (Eq. 32):
%     max_{eta_mk}  min_k  R_d,k^cf
%     s.t.  SUM_k eta_mk * gamma_mk <= 1,  for all m
%           eta_mk >= 0
%
%   Uses bisection on the target SINR. At each bisection step, solves
%   an SOCP feasibility problem (Eq. 35) using CVX if available.
%   Falls back to heuristic SINR balancing if CVX is not available.
%
%   Inputs:
%       beta        - M x K large-scale fading coefficients
%       gamma       - M x K channel estimate variances
%       pilot_index - K x 1 pilot assignment vector
%       p           - Parameter struct from params()
%
%   Output:
%       eta         - M x K DL power control coefficients

    [M, K] = size(beta);

    % Check if CVX is available
    has_cvx = exist('cvx_begin', 'file') > 0;

    if has_cvx
        eta = dl_power_control_cvx(beta, gamma, pilot_index, p);
    else
        eta = dl_power_control_heuristic(beta, gamma, pilot_index, p);
    end
end


function eta = dl_power_control_cvx(beta, gamma, pilot_index, p)
%DL_POWER_CONTROL_CVX Optimal DL max-min power control via bisection + SOCP.
%   Implements Algorithm 2 from the paper.

    [M, K] = size(beta);

    % Precompute pilot overlap indicator
    xi = zeros(K, K);   % xi(k',k) = |phi_k'^H phi_k|^2
    for k = 1:K
        for kp = 1:K
            if pilot_index(k) == pilot_index(kp)
                xi(kp, k) = 1;
            end
        end
    end

    % Bisection bounds
    t_min = 0;
    t_max = 1000;  % Upper bound on achievable SINR

    % Reduce t_max: compute SINR with equal power to get reasonable bound
    eta_eq = zeros(M, K);
    for m = 1:M
        eta_eq(m,:) = 1 / sum(gamma(m,:));
    end
    [R_eq, ~] = compute_SE_CF_DL(beta, gamma, pilot_index, p, eta_eq);
    t_max = 2^(max(R_eq)) - 1;
    t_max = max(t_max, 1);

    eta_best = eta_eq;

    for iter = 1:p.bisection_maxiter
        t = (t_min + t_max) / 2;

        % Solve SOCP feasibility problem (Eq. 35)
        [feasible, eta_trial] = check_feasibility_socp(beta, gamma, xi, p, t, M, K);

        if feasible
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


function [feasible, eta_out] = check_feasibility_socp(beta, gamma, xi, p, t, M, K)
%CHECK_FEASIBILITY_SOCP Solve the SOCP feasibility problem (Eq. 35).

    feasible = false;
    eta_out = [];

    try
        cvx_begin quiet
            variable s(M, K) nonneg        % s_mk = sqrt(eta_mk)
            variable rho_s(K, K) nonneg     % slack rho_{k'k}
            variable theta(M) nonneg        % slack theta_m

            % For each user k: ||v_k|| <= (1/sqrt(t)) * SUM_m gamma_mk * s_mk
            for k = 1:K
                % Build v_k components
                % v_k1: pilot contamination terms
                vk1_entries = [];
                for kp = 1:K
                    if kp ~= k && xi(kp, k) > 0
                        vk1_entries = [vk1_entries; rho_s(kp, k)];  %#ok
                    end
                end

                % v_k2: non-coherent + noise terms
                vk2 = sqrt(beta(:, k)) .* theta;

                % Full v_k vector
                if isempty(vk1_entries)
                    v_k = [vk2; 1/sqrt(p.rho_d_cf)];
                else
                    v_k = [vk1_entries; vk2; 1/sqrt(p.rho_d_cf)];
                end

                % SOC constraint
                norm(v_k) <= (1/sqrt(t)) * sum(gamma(:,k) .* s(:,k));
            end

            % Per-AP power constraint: SUM_k gamma_mk * s_mk^2 <= theta_m^2
            for m = 1:M
                norm(sqrt(gamma(m,:)') .* s(m,:)') <= theta(m);
                theta(m) <= 1;
            end

            % Slack variable constraints
            for kp = 1:K
                for k = 1:K
                    if kp ~= k && xi(kp, k) > 0
                        sum(gamma(:,kp) .* (beta(:,k) ./ beta(:,kp)) .* s(:,kp)) <= rho_s(kp, k);
                    end
                end
            end

        cvx_end

        if strfind(cvx_status, 'Solved') %#ok
            feasible = true;
            eta_out = s.^2;
        end
    catch
        feasible = false;
    end
end


function eta = dl_power_control_heuristic(beta, gamma, pilot_index, p)
%DL_POWER_CONTROL_HEURISTIC Heuristic max-min DL power control.
%   Iterative SINR balancing without CVX.
%   Uses coordinate descent: adjust each user's power to equalize SINRs.

    [M, K] = size(beta);

    % Start with equal power
    eta = zeros(M, K);
    for m = 1:M
        eta(m,:) = 1 / sum(gamma(m,:));
    end

    for iter = 1:20
        % Compute current rates
        [R, ~] = compute_SE_CF_DL(beta, gamma, pilot_index, p, eta);
        target_rate = min(R);

        if target_rate <= 0
            break;
        end

        % Scale power coefficients: increase power for weak users
        for k = 1:K
            if R(k) > 0
                scale = target_rate / R(k);
                % Smooth scaling to avoid oscillation
                scale = sqrt(scale);
                eta(:,k) = eta(:,k) * scale;
            end
        end

        % Re-normalize to satisfy per-AP power constraints
        for m = 1:M
            total = sum(eta(m,:) .* gamma(m,:));
            if total > 1
                eta(m,:) = eta(m,:) / total;
            end
        end
    end

end
