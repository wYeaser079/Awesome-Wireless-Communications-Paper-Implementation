function SE = functionComputeSE_L2_LMMSE_ngo(gainOverNoisedB, tau_c, tau_p, ...
    pilotIndex, p, K, L, nbrOfRealizations)
%functionComputeSE_L2_LMMSE_ngo Compute uplink SE for cell-free Level 2
%with L-MMSE combining using Monte Carlo simulation (N=1, three-slope model).
%
%   Used for Figure 4 to show the improvement of L-MMSE over MR at Level 2.
%
%   Inputs:
%       gainOverNoisedB   - L x K gain-over-noise in dB
%       tau_c             - Coherence block length
%       tau_p             - Pilot length
%       pilotIndex        - K x 1 pilot indices
%       p                 - K x 1 transmit powers
%       K, L              - Number of UEs and APs
%       nbrOfRealizations - Number of Monte Carlo realizations
%
%   Output:
%       SE - K x 1 per-user SE (bit/s/Hz)

sigma2 = 1;
prelogFactor = (tau_c - tau_p) / tau_c;
N = 1;  % Single antenna per AP

beta = db2pow(gainOverNoisedB);  % L x K

% Generate channels and estimates
R = zeros(N, N, L, K);
for l = 1:L
    for k = 1:K
        R(:, :, l, k) = beta(l, k);
    end
end

[Hhat, H, B, ~] = functionChannelEstimates(R, nbrOfRealizations, L, K, N, tau_p, pilotIndex, p);

% Compute Level 2 SE with L-MMSE combining using uniform weights
% g_ki = [v_k1^H h_i1, ..., v_kL^H h_iL]^T with L-MMSE combining

gMMSE = zeros(L, K, K, nbrOfRealizations);  % g_{ki,l}
normVkl = zeros(L, K, nbrOfRealizations);

for n = 1:nbrOfRealizations
    for l = 1:L
        % Local covariance
        Sigma_l = sigma2;
        for i = 1:K
            hhat_il = Hhat(:, n, l, i);
            Sigma_l = Sigma_l + p(i) * (abs(hhat_il)^2 + B(:, :, l, i));
        end

        for k = 1:K
            hhat_kl = Hhat(:, n, l, k);
            v_LMMSE = p(k) * hhat_kl / Sigma_l;
            normVkl(l, k, n) = abs(v_LMMSE)^2;

            for i = 1:K
                h_il = H(:, n, l, i);
                gMMSE(l, k, i, n) = conj(v_LMMSE) * h_il;
            end
        end
    end
end

% Level 2: uniform weights
SE = zeros(K, 1);
ak = ones(L, 1);

for k = 1:K
    mean_gkk = mean(squeeze(gMMSE(:, k, k, :)), 2);

    sumGG = zeros(L, L);
    for i = 1:K
        g_ki = squeeze(gMMSE(:, k, i, :));
        sumGG = sumGG + p(i) * (g_ki * g_ki') / nbrOfRealizations;
    end

    Dk = diag(mean(squeeze(normVkl(:, k, :)), 2));

    num = p(k) * abs(ak' * mean_gkk)^2;
    den = real(ak' * (sumGG + sigma2 * Dk - p(k) * (mean_gkk * mean_gkk')) * ak);
    SE(k) = prelogFactor * log2(1 + num / max(den, 1e-20));
end

end
