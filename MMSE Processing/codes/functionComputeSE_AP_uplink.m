function [SE_MR, SE_MMSE, sumSE_SIC] = functionComputeSE_AP_uplink(...
    Hhat, H, B, C, R, tau_c, tau_p, nbrOfRealizations, L, K, N, p)
%functionComputeSE_AP_uplink Compute uplink SE for cell-free mMIMO across
%all four cooperation levels with MR and MMSE/L-MMSE combining.
%
%   [SE_MR, SE_MMSE, sumSE_SIC] = functionComputeSE_AP_uplink(...)
%
%   Inputs:
%       Hhat              - N x nbrOfRealizations x L x K channel estimates
%       H                 - N x nbrOfRealizations x L x K true channels
%       B                 - N x N x L x K estimation error covariance
%       C                 - N x N x L x K estimate covariance
%       R                 - N x N x L x K spatial correlation matrices
%       tau_c             - Coherence block length
%       tau_p             - Pilot length
%       nbrOfRealizations - Number of channel realizations
%       L                 - Number of APs
%       K                 - Number of UEs
%       N                 - Number of antennas per AP
%       p                 - K x 1 transmit powers
%
%   Outputs:
%       SE_MR    - K x 4 matrix of SE values [Level4, Level3, Level2, Level1] with MR
%       SE_MMSE  - K x 4 matrix of SE values [Level4, Level3, Level2, Level1] with MMSE/L-MMSE
%       sumSE_SIC - Scalar sum SE with MMSE-SIC at Level 4
%
% Reference:
% Emil Bjornson, Luca Sanguinetti, "Making Cell-Free Massive MIMO
% Competitive With MMSE Processing and Centralized Implementation,"
% IEEE Trans. Wireless Commun., vol. 19, no. 2, pp. 1257-1271, Feb. 2020.

% Pre-log factor
prelogFactor = (tau_c - tau_p) / tau_c;

% Noise variance (normalized)
sigma2 = 1;

%% Compute Level 4 SE (Fully Centralized) - Proposition 1
% Also compute MMSE-SIC sum SE - Proposition 4

SE_L4_MR = zeros(K, 1);
SE_L4_MMSE = zeros(K, 1);
sumSE_SIC_realizations = zeros(nbrOfRealizations, 1);

% Precompute the error + noise covariance E = sum_i p_i*C_i + sigma^2*I
E = sigma2 * eye(L*N);
for i = 1:K
    Ci = zeros(L*N, L*N);
    for l = 1:L
        rows = (l-1)*N+1 : l*N;
        Ci(rows, rows) = B(:, :, l, i);
    end
    E = E + p(i) * Ci;
end

for n = 1:nbrOfRealizations

    % Stack all channel estimates and true channels into collective vectors
    Hhat_stacked = zeros(L*N, K);
    H_stacked = zeros(L*N, K);

    for k = 1:K
        for l = 1:L
            rows = (l-1)*N+1 : l*N;
            Hhat_stacked(rows, k) = Hhat(:, n, l, k);
            H_stacked(rows, k) = H(:, n, l, k);
        end
    end

    %--- Level 4: MR combining ---
    for k = 1:K
        v_MR = Hhat_stacked(:, k);

        % Signal power
        numerator = p(k) * abs(v_MR' * Hhat_stacked(:, k))^2;

        % Interference + noise (Eq. 12 denominator)
        denominator = 0;
        for i = 1:K
            if i ~= k
                denominator = denominator + p(i) * abs(v_MR' * Hhat_stacked(:, i))^2;
            end
        end
        denominator = denominator + real(v_MR' * E * v_MR);

        SE_L4_MR(k) = SE_L4_MR(k) + log2(1 + numerator / denominator);
    end

    %--- Level 4: MMSE combining (Eq. 13-14) ---
    % Matrix Sigma = sum_i p_i * hhat_i * hhat_i^H + E
    Sigma = E;
    for i = 1:K
        Sigma = Sigma + p(i) * (Hhat_stacked(:, i) * Hhat_stacked(:, i)');
    end
    Sigma_inv = inv(Sigma);

    for k = 1:K
        % SINR = p_k * hhat_k^H * (Sigma - p_k*hhat_k*hhat_k^H)^{-1} * hhat_k
        % Using matrix inversion lemma or directly from Eq. 14
        % B_k = Sigma - p_k * hhat_k * hhat_k^H
        % SINR = p_k * hhat_k^H * B_k^{-1} * hhat_k

        % More efficient: use Sigma_inv and matrix inversion lemma
        % (Sigma - p_k*hhat_k*hhat_k^H)^{-1} = Sigma_inv + ...
        hhat_k = Hhat_stacked(:, k);
        temp = Sigma_inv * hhat_k;
        SINR_MMSE = p(k) * real(hhat_k' * temp) / (1 - p(k) * real(hhat_k' * temp));

        % Ensure non-negative SINR
        SINR_MMSE = max(SINR_MMSE, 0);

        SE_L4_MMSE(k) = SE_L4_MMSE(k) + log2(1 + SINR_MMSE);
    end

    %--- MMSE-SIC sum SE (Proposition 4, Eq. 41) ---
    % sum SE = log2 det(I_K + P * Hhat^H * E^{-1} * Hhat)
    P_diag = diag(p);
    E_inv = inv(E);
    sumSE_SIC_realizations(n) = real(log2(det(eye(K) + P_diag * Hhat_stacked' * E_inv * Hhat_stacked)));
end

SE_L4_MR = prelogFactor * SE_L4_MR / nbrOfRealizations;
SE_L4_MMSE = prelogFactor * SE_L4_MMSE / nbrOfRealizations;
sumSE_SIC = prelogFactor * mean(sumSE_SIC_realizations);


%% Compute Levels 3, 2, 1 with MR and L-MMSE combining

% For Levels 3/2/1, we need local combining at each AP, then centralized processing

% Precompute local signal statistics for Level 3 LSFD
% g_ki = [v_k1^H h_i1, ..., v_kL^H h_iL]^T

% Store per-realization local estimates and statistics
% For L-MMSE combining
gMR = zeros(L, K, K, nbrOfRealizations);   % g_{ki} with MR combining for UE k
gMMSE = zeros(L, K, K, nbrOfRealizations);

% Store combining vector norms for noise term D_k
normVkl_MR = zeros(L, K, nbrOfRealizations);
normVkl_MMSE = zeros(L, K, nbrOfRealizations);

for n = 1:nbrOfRealizations
    for l = 1:L
        % Local signal covariance: Sigma_l = sum_i p_i*(hhat_il*hhat_il^H + C_il) + sigma^2*I
        Sigma_l = sigma2 * eye(N);
        for i = 1:K
            hhat_il = Hhat(:, n, l, i);
            Sigma_l = Sigma_l + p(i) * (hhat_il * hhat_il' + B(:, :, l, i));
        end
        Sigma_l_inv = inv(Sigma_l);

        for k = 1:K
            hhat_kl = Hhat(:, n, l, k);

            % MR combining: v_kl = hhat_kl
            v_MR_kl = hhat_kl;

            % L-MMSE combining: v_kl = p_k * Sigma_l^{-1} * hhat_kl (Eq. 16)
            v_LMMSE_kl = p(k) * (Sigma_l_inv * hhat_kl);

            % Store combining vector norms
            normVkl_MR(l, k, n) = norm(v_MR_kl)^2;
            normVkl_MMSE(l, k, n) = norm(v_LMMSE_kl)^2;

            % Compute g_{ki,l} = v_kl^H * h_il for all UEs i
            for i = 1:K
                h_il = H(:, n, l, i);
                gMR(l, k, i, n) = v_MR_kl' * h_il;
                gMMSE(l, k, i, n) = v_LMMSE_kl' * h_il;
            end
        end
    end
end


%% Level 3: Local processing + LSFD (Proposition 2, Corollary 2)

SE_L3_MR = zeros(K, 1);
SE_L3_MMSE = zeros(K, 1);

for k = 1:K
    % Compute statistical quantities for LSFD

    %--- MR combining ---
    % E{g_kk}: L x 1 mean effective channel
    mean_gkk_MR = mean(squeeze(gMR(:, k, k, :)), 2);

    % E{g_ki * g_ki^H}: L x L for all i
    sumGG_MR = zeros(L, L);
    for i = 1:K
        g_ki = squeeze(gMR(:, k, i, :));  % L x nbrOfRealizations
        sumGG_MR = sumGG_MR + p(i) * (g_ki * g_ki') / nbrOfRealizations;
    end

    % D_k for MR
    Dk_MR = diag(mean(squeeze(normVkl_MR(:, k, :)), 2));

    % LSFD weights (Eq. 21)
    Bk_MR = sumGG_MR + sigma2 * Dk_MR;
    ak_MR = Bk_MR \ mean_gkk_MR;

    % SINR (Eq. 22)
    num = p(k) * abs(ak_MR' * mean_gkk_MR)^2;
    den = real(ak_MR' * (sumGG_MR + sigma2 * Dk_MR - p(k) * (mean_gkk_MR * mean_gkk_MR')) * ak_MR);
    SE_L3_MR(k) = prelogFactor * log2(1 + num / max(den, 1e-20));

    %--- L-MMSE combining ---
    mean_gkk_MMSE = mean(squeeze(gMMSE(:, k, k, :)), 2);

    sumGG_MMSE = zeros(L, L);
    for i = 1:K
        g_ki = squeeze(gMMSE(:, k, i, :));
        sumGG_MMSE = sumGG_MMSE + p(i) * (g_ki * g_ki') / nbrOfRealizations;
    end

    Dk_MMSE = diag(mean(squeeze(normVkl_MMSE(:, k, :)), 2));

    Bk_MMSE = sumGG_MMSE + sigma2 * Dk_MMSE;
    ak_MMSE = Bk_MMSE \ mean_gkk_MMSE;

    num = p(k) * abs(ak_MMSE' * mean_gkk_MMSE)^2;
    den = real(ak_MMSE' * (sumGG_MMSE + sigma2 * Dk_MMSE - p(k) * (mean_gkk_MMSE * mean_gkk_MMSE')) * ak_MMSE);
    SE_L3_MMSE(k) = prelogFactor * log2(1 + num / max(den, 1e-20));
end


%% Level 2: Local processing + simple average (Corollary 3)

SE_L2_MR = zeros(K, 1);
SE_L2_MMSE = zeros(K, 1);

for k = 1:K
    % Level 2 = Level 3 with a_kl = 1 for all l (uniform weights)
    ak_uniform = ones(L, 1);

    %--- MR ---
    mean_gkk_MR = mean(squeeze(gMR(:, k, k, :)), 2);

    sumGG_MR = zeros(L, L);
    for i = 1:K
        g_ki = squeeze(gMR(:, k, i, :));
        sumGG_MR = sumGG_MR + p(i) * (g_ki * g_ki') / nbrOfRealizations;
    end
    Dk_MR = diag(mean(squeeze(normVkl_MR(:, k, :)), 2));

    num = p(k) * abs(ak_uniform' * mean_gkk_MR)^2;
    den = real(ak_uniform' * (sumGG_MR + sigma2 * Dk_MR - p(k) * (mean_gkk_MR * mean_gkk_MR')) * ak_uniform);
    SE_L2_MR(k) = prelogFactor * log2(1 + num / max(den, 1e-20));

    %--- L-MMSE ---
    mean_gkk_MMSE = mean(squeeze(gMMSE(:, k, k, :)), 2);

    sumGG_MMSE = zeros(L, L);
    for i = 1:K
        g_ki = squeeze(gMMSE(:, k, i, :));
        sumGG_MMSE = sumGG_MMSE + p(i) * (g_ki * g_ki') / nbrOfRealizations;
    end
    Dk_MMSE = diag(mean(squeeze(normVkl_MMSE(:, k, :)), 2));

    num = p(k) * abs(ak_uniform' * mean_gkk_MMSE)^2;
    den = real(ak_uniform' * (sumGG_MMSE + sigma2 * Dk_MMSE - p(k) * (mean_gkk_MMSE * mean_gkk_MMSE')) * ak_uniform);
    SE_L2_MMSE(k) = prelogFactor * log2(1 + num / max(den, 1e-20));
end


%% Level 1: Small cells (Corollary 4)
% Each UE is served by the AP that gives it the highest SE

SE_L1_MR = zeros(K, 1);
SE_L1_MMSE = zeros(K, 1);

for k = 1:K
    SE_perAP_MR = zeros(L, 1);
    SE_perAP_MMSE = zeros(L, 1);

    for l = 1:L
        % Compute SE at AP l for UE k using MR and L-MMSE
        SINR_sum_MR = 0;
        SINR_sum_MMSE = 0;

        for n = 1:nbrOfRealizations
            hhat_kl = Hhat(:, n, l, k);

            %--- MR combining ---
            v_MR = hhat_kl;
            if norm(v_MR) < 1e-15
                continue;
            end

            sig_MR = p(k) * abs(v_MR' * Hhat(:, n, l, k))^2;
            interf_MR = 0;
            for i = 1:K
                if i ~= k
                    interf_MR = interf_MR + p(i) * abs(v_MR' * Hhat(:, n, l, i))^2;
                end
            end
            noise_MR = 0;
            for i = 1:K
                noise_MR = noise_MR + p(i) * real(v_MR' * B(:, :, l, i) * v_MR);
            end
            noise_MR = noise_MR + sigma2 * (v_MR' * v_MR);

            SINR_sum_MR = SINR_sum_MR + log2(1 + real(sig_MR) / max(real(interf_MR + noise_MR), 1e-20));

            %--- L-MMSE combining ---
            Sigma_l = sigma2 * eye(N);
            for i = 1:K
                hhat_il = Hhat(:, n, l, i);
                Sigma_l = Sigma_l + p(i) * (hhat_il * hhat_il' + B(:, :, l, i));
            end

            v_MMSE = p(k) * (Sigma_l \ hhat_kl);
            if norm(v_MMSE) < 1e-15
                continue;
            end

            sig_MMSE = p(k) * abs(v_MMSE' * Hhat(:, n, l, k))^2;
            interf_MMSE = 0;
            for i = 1:K
                if i ~= k
                    interf_MMSE = interf_MMSE + p(i) * abs(v_MMSE' * Hhat(:, n, l, i))^2;
                end
            end
            noise_MMSE = 0;
            for i = 1:K
                noise_MMSE = noise_MMSE + p(i) * real(v_MMSE' * B(:, :, l, i) * v_MMSE);
            end
            noise_MMSE = noise_MMSE + sigma2 * (v_MMSE' * v_MMSE);

            SINR_sum_MMSE = SINR_sum_MMSE + log2(1 + real(sig_MMSE) / max(real(interf_MMSE + noise_MMSE), 1e-20));
        end

        SE_perAP_MR(l) = prelogFactor * SINR_sum_MR / nbrOfRealizations;
        SE_perAP_MMSE(l) = prelogFactor * SINR_sum_MMSE / nbrOfRealizations;
    end

    % Select the AP with maximum SE
    SE_L1_MR(k) = max(SE_perAP_MR);
    SE_L1_MMSE(k) = max(SE_perAP_MMSE);
end


%% Assemble output
SE_MR = [SE_L4_MR, SE_L3_MR, SE_L2_MR, SE_L1_MR];
SE_MMSE = [SE_L4_MMSE, SE_L3_MMSE, SE_L2_MMSE, SE_L1_MMSE];

end
