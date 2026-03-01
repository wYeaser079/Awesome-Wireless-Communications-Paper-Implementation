function [R, throughput] = compute_SE_SC_UL(beta, pilot_index, ap_selected, p, alpha_u)
%COMPUTE_SE_SC_UL Small-cell uplink achievable rate (Eq. 47).
%   [R, throughput] = compute_SE_SC_UL(beta, pilot_index, ap_selected, p, alpha_u)
%
%   Eq. 47:  R_u,k^sc = -(log2(e)) * exp(1/omega_bar) * Ei(-1/omega_bar)
%
%   where omega_bar is the effective SINR parameter (Eq. 48) and omega_mk
%   is the UL channel estimation quality (Eq. 49).
%
%   Inputs:
%       beta        - M x K large-scale fading coefficients
%       pilot_index - K x 1 pilot assignment (for UL pilot contamination)
%       ap_selected - K x 1 vector, ap_selected(k) = index of AP serving user k
%       p           - Parameter struct from params()
%       alpha_u     - K x 1 UL power control coefficients (optional, default = 1)
%
%   Outputs:
%       R           - K x 1 spectral efficiency [bits/s/Hz]
%       throughput  - K x 1 net throughput [bits/s] (Eq. 57)

    K = size(beta, 2);

    if nargin < 5 || isempty(alpha_u)
        alpha_u = ones(K, 1);
    end

    R = zeros(K, 1);

    for k = 1:K
        mk = ap_selected(k);   % AP serving user k

        % --- Channel estimation quality omega_mk (Eq. 49) ---
        %   omega_mk = tau_u^sc * rho_u,p^sc * beta_mk^2
        %            / (tau_u^sc * rho_u,p^sc * SUM_{k': same pilot} beta_mkk' + 1)
        same_pilot = (pilot_index == pilot_index(k));

        % Note: beta(mk, k') is the path loss from AP mk to user k'
        denom_est = p.tau_u_sc * p.rho_p_sc * sum(beta(mk, same_pilot)) + 1;
        omega_mk = p.tau_u_sc * p.rho_p_sc * beta(mk, k)^2 / denom_est;

        % --- Effective SINR parameter omega_bar (Eq. 48) ---
        %   omega_bar = rho_u^sc * alpha_u,k * omega_mk
        %             / (rho_u^sc * alpha_u,k * (beta_mkk - omega_mk)
        %                + rho_u^sc * SUM_{k'!=k} alpha_u,k' * beta_mkk' + 1)
        interf = 0;
        for kp = 1:K
            if kp ~= k
                interf = interf + alpha_u(kp) * beta(mk, kp);
            end
        end

        numer = p.rho_u_sc * alpha_u(k) * omega_mk;
        denom = p.rho_u_sc * alpha_u(k) * (beta(mk, k) - omega_mk) ...
              + p.rho_u_sc * interf + 1;
        omega_bar = numer / denom;

        % --- Rate using Ei function (Eq. 47) ---
        if omega_bar > 0
            R(k) = -log2(exp(1)) * exp(1/omega_bar) * expint_ei(-1/omega_bar);
        else
            R(k) = 0;
        end
    end

    % --- Net throughput (Eq. 57) ---
    throughput = p.B * p.prelog_sc * R;

end


function y = expint_ei(x)
%EXPINT_EI Compute the exponential integral Ei(x) for x < 0.
    if x < 0
        y = -expint(-x);
    elseif x == 0
        y = -inf;
    else
        y = -expint(-x);
    end
end
