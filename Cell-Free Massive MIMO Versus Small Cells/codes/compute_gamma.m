function gamma = compute_gamma(beta, pilot_index, p)
%COMPUTE_GAMMA Compute channel estimate variance gamma_mk (Eq. 8).
%   gamma = compute_gamma(beta, pilot_index, p)
%
%   gamma_mk = E{|g_hat_mk|^2} = sqrt(tau_cf * rho_p) * beta_mk * c_mk
%            = tau_cf * rho_p * beta_mk^2 / (tau_cf * rho_p * SUM_{k' in P_k} beta_mk' + 1)
%
%   where P_k is the set of users sharing user k's pilot (Eq. 4, 8).
%
%   Inputs:
%       beta        - M x K large-scale fading coefficients (linear)
%       pilot_index - K x 1 pilot assignment vector
%       p           - Parameter struct from params()
%
%   Output:
%       gamma       - M x K channel estimate variances

    [M, K] = size(beta);
    gamma = zeros(M, K);

    tp = p.tau_cf * p.rho_p_cf;  % tau_cf * rho_p

    for k = 1:K
        % Users sharing the same pilot as user k
        same_pilot = (pilot_index == pilot_index(k));

        for m = 1:M
            denom = tp * sum(beta(m, same_pilot)) + 1;
            gamma(m, k) = tp * beta(m, k)^2 / denom;
        end
    end

end
