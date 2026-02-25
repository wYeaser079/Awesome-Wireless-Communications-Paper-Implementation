function R_CB = compute_SE_CB(beta, alpha, pilot_index, p, eta)
%COMPUTE_SE_CB CB spectral efficiency via closed-form (Lemma 1, Eq. 33-35).
%   R_CB = compute_SE_CB(beta, alpha, pilot_index, p)
%   R_CB = compute_SE_CB(beta, alpha, pilot_index, p, eta)
%
%   Eq. 33:  SINR_k^CB = eta_k * ||l_kk||^2
%            / [ SUM_{i!=k} eta_i*||l_ik||^2
%                + (1/N) * SUM_i eta_i*||theta_ik||^2
%                + (1/(N*rho_u)) * ||l_kk||_1 ]
%
%   Eq. 34:  l_ik = |phi_k^H phi_i| * [alpha_1k*beta_1i/beta_1k, ..., alpha_Mk*beta_Mi/beta_Mk]^T
%   Eq. 35:  theta_ik = [sqrt(beta_1i * alpha_1k), ..., sqrt(beta_Mi * alpha_Mk)]^T
%
%   Inputs:
%       beta        - M x K large-scale fading coefficients
%       alpha       - M x K channel estimate variances (Eq. 8)
%       pilot_index - K x 1 pilot assignment vector
%       p           - Parameter struct from params()
%       eta         - K x 1 power control coefficients (optional, default = 1)
%
%   Output:
%       R_CB        - K x 1 per-user spectral efficiency [bits/s/Hz]

    [M, K] = size(beta);

    % Default: equal power control
    if nargin < 5
        eta = ones(K, 1);
    end

    R_CB = zeros(K, 1);

    for k = 1:K
        % --- Desired signal term: ||l_kk||_2^2 (Eq. 34 with i = k) ---
        %   l_kk = [alpha_1k, ..., alpha_Mk]^T
        %   ||l_kk||_2^2 = SUM_m alpha_mk^2
        signal = sum(alpha(:, k).^2);

        % --- Coherent interference: SUM_{i!=k} eta_i * ||l_ik||_2^2 ---
        coh_interf = 0;
        for i = 1:K
            if i == k
                continue;
            end

            % |phi_k^H * phi_i| = 1 if same pilot, 0 if orthogonal
            if pilot_index(i) ~= pilot_index(k)
                continue;   % Orthogonal pilots: l_ik = 0, no coherent interference
            end

            % Eq. 34: l_ik (same pilot case, |phi_k^H phi_i| = 1)
            l_ik = alpha(:, k) .* beta(:, i) ./ beta(:, k);
            coh_interf = coh_interf + eta(i) * sum(l_ik.^2);
        end

        % --- Non-coherent interference: (1/N) * SUM_{i=1}^K eta_i * ||theta_ik||_2^2 ---
        noncoh_interf = 0;
        for i = 1:K
            % Eq. 35: ||theta_ik||_2^2 = SUM_m beta_mi * alpha_mk
            noncoh_interf = noncoh_interf + eta(i) * sum(beta(:, i) .* alpha(:, k));
        end
        noncoh_interf = noncoh_interf / p.N;

        % --- Noise term: (1/(N*rho_u)) * ||l_kk||_1 ---
        %   ||l_kk||_1 = SUM_m alpha_mk
        noise = sum(alpha(:, k)) / (p.N * p.rho_u);

        % --- SINR (Eq. 33) ---
        sinr_k = eta(k) * signal / (coh_interf + noncoh_interf + noise);

        % --- Spectral efficiency with pre-log factor ---
        R_CB(k) = p.prelog * log2(1 + sinr_k);
    end

end
