function [R, throughput] = compute_SE_SC_DL(beta, pilot_index, ap_selected, p, alpha_d)
%COMPUTE_SE_SC_DL Small-cell downlink achievable rate (Eq. 42-43).
%   [R, throughput] = compute_SE_SC_DL(beta, pilot_index, ap_selected, p, alpha_d)
%
%   Eq. 43:  R_d,k^sc = -(log2(e)) * exp(1/mu_bar) * Ei(-1/mu_bar)
%
%   where mu_bar is the effective SINR parameter (Eq. 44) and Ei is the
%   exponential integral function.
%
%   Inputs:
%       beta        - M x K large-scale fading coefficients
%       pilot_index - K x 1 pilot assignment (for DL pilot contamination)
%       ap_selected - K x 1 vector, ap_selected(k) = index of AP serving user k
%       p           - Parameter struct from params()
%       alpha_d     - K x 1 DL power control coefficients (optional, default = 1)
%
%   Outputs:
%       R           - K x 1 spectral efficiency [bits/s/Hz]
%       throughput  - K x 1 net throughput [bits/s] (Eq. 57)

    K = size(beta, 2);

    if nargin < 5 || isempty(alpha_d)
        alpha_d = ones(K, 1);
    end

    R = zeros(K, 1);

    for k = 1:K
        mk = ap_selected(k);   % AP serving user k

        % --- Channel estimation quality mu_mk (Eq. 40) ---
        %   mu_mk = tau_d^sc * rho_d,p^sc * beta_mk^2
        %         / (tau_d^sc * rho_d,p^sc * SUM_{k': same pilot} beta_mk'k + 1)
        same_pilot = (pilot_index == pilot_index(k));
        denom_est = p.tau_d_sc * p.rho_p_sc * sum(beta(mk, same_pilot)) + 1;
        mu_mk = p.tau_d_sc * p.rho_p_sc * beta(mk, k)^2 / denom_est;

        % --- Effective SINR parameter mu_bar (Eq. 44) ---
        %   mu_bar = rho_d^sc * alpha_d,k * mu_mk
        %          / (rho_d^sc * alpha_d,k * (beta_mk - mu_mk)
        %             + rho_d^sc * SUM_{k'!=k} alpha_d,k' * beta_mk'k + 1)
        interf = 0;
        for kp = 1:K
            if kp ~= k
                interf = interf + alpha_d(kp) * beta(mk, kp);
            end
        end

        numer = p.rho_d_sc * alpha_d(k) * mu_mk;
        denom = p.rho_d_sc * alpha_d(k) * (beta(mk, k) - mu_mk) ...
              + p.rho_d_sc * interf + 1;
        mu_bar = numer / denom;

        % --- Rate using Ei function (Eq. 43) ---
        %   R_d,k^sc = -(log2(e)) * exp(1/mu_bar) * Ei(-1/mu_bar)
        if mu_bar > 0
            R(k) = -log2(exp(1)) * exp(1/mu_bar) * expint_ei(-1/mu_bar);
        else
            R(k) = 0;
        end
    end

    % --- Net throughput (Eq. 57) ---
    throughput = p.B * p.prelog_sc * R;

end


function y = expint_ei(x)
%EXPINT_EI Compute the exponential integral Ei(x) for x < 0.
%   Ei(x) = -E1(-x) for x < 0, where E1 is MATLAB's expint function.
%
%   MATLAB's expint(z) computes E1(z) = integral from z to inf of e^(-t)/t dt.
%   The relation: Ei(x) = -E1(-x) for real x < 0.

    if x < 0
        y = -expint(-x);
    elseif x == 0
        y = -inf;
    else
        % For x > 0, Ei(x) = -E1(-x) (principal value)
        % This case should not occur in our context
        y = -expint(-x);
    end
end
