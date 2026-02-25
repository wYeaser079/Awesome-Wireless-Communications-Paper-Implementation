function [alpha, c_coeff] = estimate_channel(beta, pilot_index, p)
%ESTIMATE_CHANNEL Compute MMSE channel estimation parameters (Eq. 5, 8).
%   [alpha, c_coeff] = estimate_channel(beta, pilot_index, p)
%
%   Computes the MMSE filter coefficient c_mk (Eq. 5) and channel estimate
%   variance alpha_mk (Eq. 8). These depend ONLY on large-scale fading
%   coefficients (beta) and pilot assignments — not on channel realizations.
%   They are computed once per network setup and reused across all Monte
%   Carlo channel realizations.
%
%   Eq. 5:  c_mk = sqrt(tau*rho_p)*beta_mk / (tau*rho_p * SUM_{i in P_k} beta_mi + 1)
%   Eq. 8:  alpha_mk = tau*rho_p*beta_mk^2  / (tau*rho_p * SUM_{i in P_k} beta_mi + 1)
%
%   where P_k = {i : pilot_index(i) == pilot_index(k)} (co-pilot users).
%
%   Inputs:
%       beta        - M x K large-scale fading coefficients
%       pilot_index - K x 1 pilot assignment vector
%       p           - Parameter struct from params()
%
%   Outputs:
%       alpha       - M x K channel estimate variances
%       c_coeff     - M x K MMSE filter coefficients

    [M, K] = size(beta);
    alpha   = zeros(M, K);
    c_coeff = zeros(M, K);

    for k = 1:K
        % Set of users sharing user k's pilot (including k itself)
        same_pilot = (pilot_index == pilot_index(k));

        for m = 1:M
            % Denominator (shared by Eq. 5 and Eq. 8):
            %   tau*rho_p * sum of beta_mi over co-pilot users + 1
            denom = p.tau * p.rho_p * sum(beta(m, same_pilot)) + 1;

            % Eq. 5: MMSE filter coefficient
            c_coeff(m, k) = sqrt(p.tau * p.rho_p) * beta(m, k) / denom;

            % Eq. 8: Channel estimate variance
            alpha(m, k) = p.tau * p.rho_p * beta(m, k)^2 / denom;
        end
    end

end
