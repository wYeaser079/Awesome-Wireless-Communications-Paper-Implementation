function eta = power_control_SCA_CB(beta, alpha, pilot_index, p, R_bar)
%POWER_CONTROL_SCA_CB Algorithm 3: SCA for CB total rate maximization (Eq. 36-43).
%   eta = power_control_SCA_CB(beta, alpha, pilot_index, p)
%   eta = power_control_SCA_CB(beta, alpha, pilot_index, p, R_bar)
%
%   Eq. 36:  maximize    SUM_k R_k^CB({eta_k})
%            subject to  R_k^CB >= R_bar_k,  for all k   (QoS constraint)
%                        0 <= eta_k <= 1                  (power constraint)
%
%   where R_k^CB = prelog * log2(1 + eta_k * S_k / D_k)  (Lemma 1, Eq. 33)
%         S_k   = ||l_kk||^2 = SUM_m alpha_mk^2          (signal, Eq. 34)
%         D_k   = SUM_j eta_j * w_{jk} + noise_k         (interference + noise)
%         noise_k = SUM_m alpha_mk / (N * rho_u)
%
%   DC programming: R_k = log(eta_k*S_k + D_k) - log(D_k)
%   At each SCA iteration n, linearize -log(D_k) around eta^(n):
%     nabla_j (-log D_k) = -w_{jk} / D_k^(n)
%
%   Requires: CVX (http://cvxr.com/cvx/)
%
%   Inputs:
%       beta        - M x K large-scale fading coefficients
%       alpha       - M x K channel estimate variances (Eq. 8)
%       pilot_index - K x 1 pilot assignment vector
%       p           - Parameter struct from params()
%       R_bar       - K x 1 QoS targets [bits/s/Hz] (default = 0)
%
%   Output:
%       eta         - K x 1 optimized power control coefficients

    [M, K] = size(beta);

    if nargin < 5
        R_bar = zeros(K, 1);
    end

    % --- Precompute CB parameters (Lemma 1, Eq. 33-35) ---

    % Signal: S_k = ||l_kk||^2 = SUM_m alpha_mk^2 (Eq. 34 with i = k)
    S = sum(alpha.^2, 1).';                    % K x 1

    % Interference weight matrix w(i,k):
    %   For i != k, same pilot: w(i,k) = ||l_ik||^2 + ||theta_ik||^2 / N
    %   For i != k, diff pilot: w(i,k) = ||theta_ik||^2 / N
    %   For i  = k:             w(k,k) = ||theta_kk||^2 / N
    w = zeros(K, K);
    for k = 1:K
        for i = 1:K
            % Eq. 35: non-coherent — ||theta_ik||^2 / N = SUM_m beta_mi*alpha_mk / N
            noncoh = sum(beta(:, i) .* alpha(:, k)) / p.N;

            % Eq. 34: coherent — ||l_ik||^2 (only if co-pilot and i != k)
            coh = 0;
            if i ~= k && pilot_index(i) == pilot_index(k)
                l_ik = alpha(:, k) .* beta(:, i) ./ beta(:, k);
                coh = sum(l_ik.^2);
            end

            w(i, k) = coh + noncoh;
        end
    end

    % Noise: noise_k = ||l_kk||_1 / (N * rho_u) = SUM_m alpha_mk / (N * rho_u)
    noise = sum(alpha, 1).' / (p.N * p.rho_u);    % K x 1

    % QoS thresholds
    gamma = 2.^(R_bar / p.prelog) - 1;

    % --- SCA iterations (Algorithm 3) ---
    eta      = ones(K, 1);
    max_iter = 20;
    epsilon  = 1e-4;

    for iter = 1:max_iter
        % Per-user denominators at current point: D_k^(n) = w(:,k)'*eta + noise(k)
        D_n = w' * eta + noise;                 % K x 1

        % --- Solve convex subproblem via CVX ---
        cvx_begin quiet
            variable eta_new(K, 1)

            % Objective: SUM_k [log(eta_k*S_k + D_k) - nabla_k' * eta]
            obj = 0;
            for k = 1:K
                % Concave part: log(eta_k * S_k + w(:,k)' * eta + noise_k)
                obj = obj + log(eta_new(k) * S(k) ...
                            + w(:,k)' * eta_new + noise(k));

                % Linearized penalty: -w(:,k)' * eta / D_k^(n)
                obj = obj - w(:,k)' * eta_new / D_n(k);
            end

            maximize(obj)

            subject to
                % QoS: S_k*eta_k - gamma_k * D_k >= gamma_k * noise_k  (linear)
                for k = 1:K
                    if gamma(k) > 0
                        S(k) * eta_new(k) ...
                            - gamma(k) * (w(:,k)' * eta_new) >= gamma(k) * noise(k);
                    end
                end

                eta_new >= 0;
                eta_new <= 1;
        cvx_end

        % Check CVX status
        if isempty(strfind(cvx_status, 'Solved'))
            warning('SCA_CB: CVX failed at iteration %d (status: %s).', ...
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
