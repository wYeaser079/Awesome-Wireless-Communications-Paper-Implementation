function [alpha_d, alpha_u] = power_control_SC(beta, pilot_index, ap_selected, p)
%POWER_CONTROL_SC Max-min power control for small-cell system (Eq. 45-46, 50).
%   [alpha_d, alpha_u] = power_control_SC(beta, pilot_index, ap_selected, p)
%
%   Downlink (Eq. 45-46): max min_k mu_bar_mk,  s.t. 0 <= alpha_d,k <= 1
%   Uplink (Eq. 50):      max min_k omega_bar_mk, s.t. 0 <= alpha_u,k <= 1
%
%   Both are quasi-linear and solved using bisection with LP feasibility checks.
%
%   Inputs:
%       beta        - M x K large-scale fading coefficients
%       pilot_index - K x 1 pilot assignment vector
%       ap_selected - K x 1 AP selection vector
%       p           - Parameter struct from params()
%
%   Outputs:
%       alpha_d     - K x 1 DL power control coefficients
%       alpha_u     - K x 1 UL power control coefficients

    K = size(beta, 2);

    % --- Downlink power control ---
    alpha_d = sc_power_control_bisection(beta, pilot_index, ap_selected, p, 'DL');

    % --- Uplink power control ---
    alpha_u = sc_power_control_bisection(beta, pilot_index, ap_selected, p, 'UL');

end


function alpha = sc_power_control_bisection(beta, pilot_index, ap_selected, p, direction)
%SC_POWER_CONTROL_BISECTION Bisection-based max-min power control for small cells.

    K = size(beta, 2);

    % Precompute per-user channel estimation quality and interference terms
    mu_or_omega = zeros(K, 1);
    beta_self   = zeros(K, 1);   % beta_{mk,k}
    beta_cross  = zeros(K, K);   % beta_{mk,k'} for k' != k

    for k = 1:K
        mk = ap_selected(k);
        beta_self(k) = beta(mk, k);

        same_pilot = (pilot_index == pilot_index(k));

        if strcmp(direction, 'DL')
            denom_est = p.tau_d_sc * p.rho_p_sc * sum(beta(mk, same_pilot)) + 1;
            mu_or_omega(k) = p.tau_d_sc * p.rho_p_sc * beta(mk, k)^2 / denom_est;
        else
            denom_est = p.tau_u_sc * p.rho_p_sc * sum(beta(mk, same_pilot)) + 1;
            mu_or_omega(k) = p.tau_u_sc * p.rho_p_sc * beta(mk, k)^2 / denom_est;
        end

        for kp = 1:K
            beta_cross(k, kp) = beta(mk, kp);
        end
    end

    if strcmp(direction, 'DL')
        rho = p.rho_d_sc;
    else
        rho = p.rho_u_sc;
    end

    % Bisection on target mu_bar (or omega_bar)
    t_min = 0;
    t_max = 100;    % Generous upper bound
    alpha_best = ones(K, 1);

    for iter = 1:p.bisection_maxiter
        t = (t_min + t_max) / 2;

        % For each k, mu_bar_k >= t means:
        %   rho * alpha_k * mu_k >= t * [rho * alpha_k * (beta_self_k - mu_k)
        %                                + rho * SUM_{k'!=k} alpha_k' * beta_cross(k,k') + 1]
        %
        % Rearranging: A * alpha >= b
        A = zeros(K, K);
        b = zeros(K, 1);

        for k = 1:K
            A(k, k) = rho * mu_or_omega(k) - t * rho * (beta_self(k) - mu_or_omega(k));
            for kp = 1:K
                if kp ~= k
                    A(k, kp) = -t * rho * beta_cross(k, kp);
                end
            end
            b(k) = t;
        end

        % LP feasibility check
        options = optimoptions('linprog', 'Display', 'off');
        f = zeros(K, 1);
        [alpha_trial, ~, exitflag] = linprog(f, -A, -b, [], [], ...
                                             zeros(K,1), ones(K,1), options);

        if exitflag == 1
            t_min = t;
            alpha_best = alpha_trial;
        else
            t_max = t;
        end

        if t_max - t_min < p.bisection_tol * 0.01
            break;
        end
    end

    alpha = alpha_best;
end
