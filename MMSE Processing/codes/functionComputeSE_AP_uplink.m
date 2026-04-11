function [SE_MR, SE_MMSE, sumSE_SIC] = functionComputeSE_AP_uplink(...
    Hhat, H, B, C, R, tau_c, tau_p, nbrOfRealizations, L, K, N, p)
%functionComputeSE_AP_uplink Compute uplink SE for cell-free mMIMO across
%all four cooperation levels with MR and MMSE/L-MMSE combining.
%
%   [SE_MR, SE_MMSE, sumSE_SIC] = functionComputeSE_AP_uplink(...)
%
%   Vectorised implementation for speed (per-(n,l) matrix operations).
%
%   Inputs:
%       Hhat              - N x nbrOfRealizations x L x K channel estimates
%       H                 - N x nbrOfRealizations x L x K true channels
%       B                 - N x N x L x K estimation error covariance
%       C                 - N x N x L x K estimate covariance (unused)
%       R                 - N x N x L x K spatial correlation matrices (unused)
%       tau_c             - Coherence block length
%       tau_p             - Pilot length
%       nbrOfRealizations - Number of channel realizations
%       L                 - Number of APs
%       K                 - Number of UEs
%       N                 - Number of antennas per AP
%       p                 - K x 1 transmit powers
%
% Reference: Bjornson & Sanguinetti, IEEE TWC, 2020.

prelogFactor = (tau_c - tau_p) / tau_c;
sigma2 = 1;
p_col = p(:);

%% Level 4 (Proposition 1) and MMSE-SIC (Proposition 4)
SE_L4_MR = zeros(K, 1);
SE_L4_MMSE = zeros(K, 1);
sumSE_SIC_realizations = zeros(nbrOfRealizations, 1);

% E = sum_i p_i*B_i (block diagonal) + sigma^2 * I
E = sigma2 * eye(L*N);
for l = 1:L
    rows = (l-1)*N+1 : l*N;
    Bl_sum = zeros(N, N);
    for i = 1:K
        Bl_sum = Bl_sum + p(i) * B(:, :, l, i);
    end
    E(rows, rows) = E(rows, rows) + Bl_sum;
end
E_inv = inv(E);

progressStep = max(1, floor(nbrOfRealizations / 10));
for n = 1:nbrOfRealizations
    if mod(n, progressStep) == 0 || n == nbrOfRealizations
        fprintf('    [L4] %d/%d\n', n, nbrOfRealizations);
    end
    % Stack channels across APs for this realization: (L*N) x K
    Hhat_stk = reshape(permute(Hhat(:, n, :, :), [1 3 2 4]), L*N, K);
    H_stk    = reshape(permute(H(:, n, :, :),    [1 3 2 4]), L*N, K);

    % --- Level 4 MR ---
    % inner products Hhat_stk' * Hhat_stk (K x K), diagonal is signal
    G = Hhat_stk' * Hhat_stk;                      % K x K
    signalPow = p_col .* abs(diag(G)).^2;          % K x 1
    absG2 = abs(G).^2;                             % K x K
    interfPow = absG2 * p_col - p_col .* abs(diag(G)).^2;  % exclude self
    % noise: diag( Hhat_stk' * E * Hhat_stk )
    noisePow = real(sum(conj(Hhat_stk) .* (E * Hhat_stk), 1)).';
    SINR_MR = signalPow ./ (interfPow + noisePow);
    SE_L4_MR = SE_L4_MR + log2(1 + max(SINR_MR, 0));

    % --- Level 4 MMSE (Eq. 14 via matrix inversion lemma) ---
    Sigma = E + Hhat_stk * diag(p_col) * Hhat_stk';
    Sigma_inv = inv(Sigma);
    % t_k = hhat_k' * Sigma_inv * hhat_k
    tvec = real(sum(conj(Hhat_stk) .* (Sigma_inv * Hhat_stk), 1)).';
    denom = max(1 - p_col .* tvec, 1e-12);
    SINR_MMSE = p_col .* tvec ./ denom;
    SE_L4_MMSE = SE_L4_MMSE + log2(1 + max(SINR_MMSE, 0));

    % --- MMSE-SIC sum SE (Prop. 4) ---
    sumSE_SIC_realizations(n) = real(log2(det(eye(K) + diag(p_col) * (Hhat_stk' * (E_inv * Hhat_stk)))));
end

SE_L4_MR   = prelogFactor * SE_L4_MR / nbrOfRealizations;
SE_L4_MMSE = prelogFactor * SE_L4_MMSE / nbrOfRealizations;
sumSE_SIC  = prelogFactor * mean(sumSE_SIC_realizations);


%% Level 2/3 precomputation (vectorised per (l, n))
% gMR(l, k, i, n)   = v_MR_kl'   * h_il      (MR: v_kl = hhat_kl)
% gMMSE(l, k, i, n) = v_LMMSE_kl' * h_il     (v_LMMSE_kl = p_k * Sigma_l^-1 * hhat_kl)
gMR   = zeros(L, K, K, nbrOfRealizations);
gMMSE = zeros(L, K, K, nbrOfRealizations);
normVkl_MR   = zeros(L, K, nbrOfRealizations);
normVkl_MMSE = zeros(L, K, nbrOfRealizations);

% Precompute per-AP error covariance sum: Bsum_l = sigma2*I + sum_i p_i * B_il
Bsum = zeros(N, N, L);
for l = 1:L
    Bl = sigma2 * eye(N);
    for i = 1:K
        Bl = Bl + p(i) * B(:, :, l, i);
    end
    Bsum(:, :, l) = Bl;
end

for n = 1:nbrOfRealizations
    if mod(n, progressStep) == 0 || n == nbrOfRealizations
        fprintf('    [L2/L3] %d/%d\n', n, nbrOfRealizations);
    end
    for l = 1:L
        % Channel matrices at this (n,l)
        Hhat_nl = squeeze(Hhat(:, n, l, :));   % N x K
        H_nl    = squeeze(H(:, n, l, :));      % N x K
        if N == 1
            Hhat_nl = Hhat_nl(:).';            % 1 x K
            H_nl    = H_nl(:).';
        end

        % Sigma_l = Bsum_l + sum_i p_i * hhat_il * hhat_il' = Bsum_l + Hhat_nl * diag(p) * Hhat_nl'
        Sigma_l = Bsum(:, :, l) + Hhat_nl * (p_col .* Hhat_nl');
        Sigma_l_inv = inv(Sigma_l);

        % MR: V_MR = Hhat_nl
        V_MR = Hhat_nl;                        % N x K
        % LMMSE: V_LMMSE = Sigma_l_inv * Hhat_nl * diag(p)
        V_LMMSE = Sigma_l_inv * Hhat_nl .* p_col.';  % N x K

        % Norms (K x 1)
        normVkl_MR(l, :, n)   = sum(abs(V_MR).^2, 1);
        normVkl_MMSE(l, :, n) = sum(abs(V_LMMSE).^2, 1);

        % g_{ki,l} = V(:,k)' * H_nl(:,i)  ->  G_all = V' * H_nl  (K x K)
        G_MR   = V_MR'   * H_nl;
        G_LMMSE = V_LMMSE' * H_nl;

        gMR(l, :, :, n)   = reshape(G_MR,   [1 K K]);
        gMMSE(l, :, :, n) = reshape(G_LMMSE, [1 K K]);
    end
end


%% Level 3 (LSFD, Corollary 2) and Level 2 (uniform, Corollary 3)
SE_L3_MR   = zeros(K, 1);
SE_L3_MMSE = zeros(K, 1);
SE_L2_MR   = zeros(K, 1);
SE_L2_MMSE = zeros(K, 1);

for k = 1:K
    % --- MR stats ---
    gkk = squeeze(gMR(:, k, k, :));           % L x nbrOfRealizations
    mean_gkk_MR = mean(gkk, 2);
    sumGG_MR = zeros(L, L);
    for i = 1:K
        gki = squeeze(gMR(:, k, i, :));
        sumGG_MR = sumGG_MR + p(i) * (gki * gki') / nbrOfRealizations;
    end
    Dk_MR = diag(mean(squeeze(normVkl_MR(:, k, :)), 2));

    Bk_MR = sumGG_MR + sigma2 * Dk_MR;
    ak_MR = Bk_MR \ mean_gkk_MR;

    num = p(k) * abs(ak_MR' * mean_gkk_MR)^2;
    den = real(ak_MR' * (sumGG_MR + sigma2*Dk_MR - p(k)*(mean_gkk_MR*mean_gkk_MR')) * ak_MR);
    SE_L3_MR(k) = prelogFactor * log2(1 + num / max(den, 1e-20));

    ak_uniform = ones(L, 1);
    num = p(k) * abs(ak_uniform' * mean_gkk_MR)^2;
    den = real(ak_uniform' * (sumGG_MR + sigma2*Dk_MR - p(k)*(mean_gkk_MR*mean_gkk_MR')) * ak_uniform);
    SE_L2_MR(k) = prelogFactor * log2(1 + num / max(den, 1e-20));

    % --- L-MMSE stats ---
    gkk = squeeze(gMMSE(:, k, k, :));
    mean_gkk_MMSE = mean(gkk, 2);
    sumGG_MMSE = zeros(L, L);
    for i = 1:K
        gki = squeeze(gMMSE(:, k, i, :));
        sumGG_MMSE = sumGG_MMSE + p(i) * (gki * gki') / nbrOfRealizations;
    end
    Dk_MMSE = diag(mean(squeeze(normVkl_MMSE(:, k, :)), 2));

    Bk_MMSE = sumGG_MMSE + sigma2 * Dk_MMSE;
    ak_MMSE = Bk_MMSE \ mean_gkk_MMSE;

    num = p(k) * abs(ak_MMSE' * mean_gkk_MMSE)^2;
    den = real(ak_MMSE' * (sumGG_MMSE + sigma2*Dk_MMSE - p(k)*(mean_gkk_MMSE*mean_gkk_MMSE')) * ak_MMSE);
    SE_L3_MMSE(k) = prelogFactor * log2(1 + num / max(den, 1e-20));

    num = p(k) * abs(ak_uniform' * mean_gkk_MMSE)^2;
    den = real(ak_uniform' * (sumGG_MMSE + sigma2*Dk_MMSE - p(k)*(mean_gkk_MMSE*mean_gkk_MMSE')) * ak_uniform);
    SE_L2_MMSE(k) = prelogFactor * log2(1 + num / max(den, 1e-20));
end


%% Level 1 (Small cells, Corollary 4) - vectorised
% Each UE is served by the AP that gives the highest per-user SE with
% local MR or local L-MMSE combining.

SE_perAP_MR   = zeros(L, K);
SE_perAP_MMSE = zeros(L, K);

apProgressStep = max(1, floor(L / 10));
for l = 1:L
    if mod(l, apProgressStep) == 0 || l == L
        fprintf('    [L1] AP %d/%d\n', l, L);
    end
    SINRsum_MR   = zeros(K, 1);
    SINRsum_MMSE = zeros(K, 1);

    for n = 1:nbrOfRealizations
        Hhat_nl = squeeze(Hhat(:, n, l, :));
        H_nl    = squeeze(H(:, n, l, :));
        if N == 1
            Hhat_nl = Hhat_nl(:).';
            H_nl    = H_nl(:).';
        end

        Sigma_l = Bsum(:, :, l) + Hhat_nl * (p_col .* Hhat_nl');
        Sigma_l_inv = inv(Sigma_l);

        % MR combining: v_k = hhat_kl
        V_MR = Hhat_nl;
        % Inner products V_MR' * H_nl  (K x K)
        Vm_H = V_MR' * H_nl;
        absVH2_MR = abs(Vm_H).^2;
        sig_MR    = p_col .* abs(diag(Vm_H)).^2;
        total_MR  = absVH2_MR * p_col;
        interf_MR = total_MR - sig_MR;
        % noise: v_k' * Bsum_l * v_k + sigma2 already in Bsum_l? Bsum includes sigma2*I.
        % noise term here = v_k' * (sigma2*I + sum_i p_i * B_il) * v_k
        noise_MR = real(sum(conj(V_MR) .* (Bsum(:,:,l) * V_MR), 1)).';
        SINRsum_MR = SINRsum_MR + log2(1 + sig_MR ./ max(interf_MR + noise_MR, 1e-20));

        % L-MMSE combining: V = p_k * Sigma_l_inv * hhat_k
        V_LMMSE = Sigma_l_inv * Hhat_nl .* p_col.';
        Vm_H = V_LMMSE' * H_nl;
        absVH2_MMSE = abs(Vm_H).^2;
        sig_MMSE    = p_col .* abs(diag(Vm_H)).^2;
        total_MMSE  = absVH2_MMSE * p_col;
        interf_MMSE = total_MMSE - sig_MMSE;
        noise_MMSE  = real(sum(conj(V_LMMSE) .* (Bsum(:,:,l) * V_LMMSE), 1)).';
        SINRsum_MMSE = SINRsum_MMSE + log2(1 + sig_MMSE ./ max(interf_MMSE + noise_MMSE, 1e-20));
    end

    SE_perAP_MR(l, :)   = prelogFactor * (SINRsum_MR   / nbrOfRealizations).';
    SE_perAP_MMSE(l, :) = prelogFactor * (SINRsum_MMSE / nbrOfRealizations).';
end

SE_L1_MR   = max(SE_perAP_MR,   [], 1).';
SE_L1_MMSE = max(SE_perAP_MMSE, [], 1).';


%% Assemble output
SE_MR   = [SE_L4_MR,   SE_L3_MR,   SE_L2_MR,   SE_L1_MR];
SE_MMSE = [SE_L4_MMSE, SE_L3_MMSE, SE_L2_MMSE, SE_L1_MMSE];

end
