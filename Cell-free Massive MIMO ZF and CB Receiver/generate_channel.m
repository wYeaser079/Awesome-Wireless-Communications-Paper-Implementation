function [G, G_hat] = generate_channel(beta, alpha, c_coeff, pilot_index, p)
%GENERATE_CHANNEL Generate one channel realization and MMSE estimates.
%   [G, G_hat] = generate_channel(beta, alpha, c_coeff, pilot_index, p)
%
%   For one Monte Carlo realization:
%     1. Generates true channels g_mk = sqrt(beta_mk)*h_mk        (Eq. 1)
%     2. Simulates pilot transmission at each AP                    (Eq. 2)
%     3. De-spreads pilot signal for each user                      (Eq. 3)
%     4. Computes MMSE channel estimates g_hat_mk = c_mk*Y_tilde   (Eq. 4)
%
%   Inputs:
%       beta        - M x K large-scale fading coefficients
%       alpha       - M x K channel estimate variances (from estimate_channel)
%       c_coeff     - M x K MMSE filter coefficients   (from estimate_channel)
%       pilot_index - K x 1 pilot assignment vector
%       p           - Parameter struct from params()
%
%   Outputs:
%       G           - MN x K true channel matrix
%       G_hat       - MN x K estimated channel matrix

    M = p.M;
    K = p.K;
    N = p.N;

    % --- Step 1: Generate true channels (Eq. 1) ---
    %   g_mk = sqrt(beta_mk) * h_mk,  where h_mk ~ CN(0, I_N)
    G = zeros(M*N, K);

    for m = 1:M
        for k = 1:K
            h_mk = (1/sqrt(2)) * (randn(N, 1) + 1i * randn(N, 1));
            rows = (m-1)*N + 1 : m*N;
            G(rows, k) = sqrt(beta(m, k)) * h_mk;
        end
    end

    % --- Steps 2-4: Pilot transmission + MMSE estimation ---
    G_hat = zeros(M*N, K);

    for m = 1:M
        rows = (m-1)*N + 1 : m*N;

        % Pilot noise at AP m: W_m,p is N x tau, entries ~ CN(0, 1)
        W_mp = (1/sqrt(2)) * (randn(N, p.tau) + 1i * randn(N, p.tau));

        % Pre-compute de-spread noise for each of the tau pilots:
        %   n_t = W_m,p * phi_t  (N x 1, still CN(0, I_N))
        % Using phi_t = e_t (standard basis), this is just column t of W_mp.
        % Users sharing the same pilot see the same de-spread noise.

        for k = 1:K
            % De-spread signal (Eq. 3):
            %   Y_tilde_mk = sqrt(tau*rho_p) * sum_{i in P_k} g_mi + n_t
            %
            % where P_k = {i : pilot_index(i) == pilot_index(k)}
            same_pilot = (pilot_index == pilot_index(k));
            signal = sqrt(p.tau * p.rho_p) * sum(G(rows, same_pilot), 2);
            noise  = W_mp(:, pilot_index(k));

            Y_tilde = signal + noise;

            % MMSE estimate (Eq. 4): g_hat_mk = c_mk * Y_tilde_mk
            G_hat(rows, k) = c_coeff(m, k) * Y_tilde;
        end
    end

end
