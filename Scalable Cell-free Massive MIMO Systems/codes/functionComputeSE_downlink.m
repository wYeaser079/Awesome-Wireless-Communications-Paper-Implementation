function [SE_MR, SE_LP_MMSE, SE_MR_perfect, SE_LP_MMSE_perfect, ...
    SE_P_MMSE, SE_P_MMSE_perfect] = ...
    functionComputeSE_downlink(Hhat, H, D, B, C, tau_c, tau_p, ...
    nbrOfRealizations, N, K, L, p, rho_dist, R, pilotIndex, rho_central)
%functionComputeSE_downlink Compute downlink spectral efficiency for
%different precoding schemes in Cell-Free Massive MIMO.
%
%   [SE_MR, SE_LP_MMSE, SE_MR_perfect, SE_LP_MMSE_perfect, ...
%    SE_P_MMSE, SE_P_MMSE_perfect] = ...
%       functionComputeSE_downlink(Hhat, H, D, B, C, tau_c, tau_p, ...
%       nbrOfRealizations, N, K, L, p, rho_dist, R, pilotIndex, rho_central)
%
%   Implements three precoding schemes with hardening and perfect-CSI bounds:
%     - MR: Distributed maximum ratio precoding
%     - LP-MMSE: Distributed local partial MMSE precoding
%     - P-MMSE: Centralized partial MMSE precoding
%
%   INPUTS:
%       Hhat, H, D, B, C  - Channel data (see functionChannelEstimates)
%       tau_c, tau_p       - Coherence block parameters
%       nbrOfRealizations  - Number of channel realizations
%       N, K, L            - System dimensions
%       p                  - UL transmit powers (for combining design)
%       rho_dist           - L x K distributed DL power allocation
%       R                  - N x N x L x K spatial correlation matrices
%       pilotIndex         - K x 1 pilot assignment indices
%       rho_central        - K x 1 centralized DL power allocation
%
%   OUTPUTS:
%       SE_MR, SE_LP_MMSE              - Hardening bound SE (K x 1)
%       SE_MR_perfect, SE_LP_MMSE_perfect - Perfect CSI SE (K x 1)
%       SE_P_MMSE, SE_P_MMSE_perfect   - Centralized P-MMSE SE (K x 1)
%
%   REFERENCES:
%   [1] E. Bjornson, L. Sanguinetti, "Scalable Cell-Free Massive MIMO
%       Systems," IEEE Trans. Commun., vol. 68, no. 7, pp. 4247-4261, 2020.

prelogFactor = (tau_c - tau_p) / tau_c;

if isscalar(p)
    p = p * ones(K, 1);
end
p = p(:);

%% ====================================================================
%  PASS 1: Compute combining vectors and accumulate E{||v||^2} for
%  normalization (needed for LP-MMSE and P-MMSE precoding)
%  ====================================================================

% Scaling factors for LP-MMSE: E{||v_{il}||^2} per AP per UE
scaling_LP = zeros(L, K);

% Scaling factors for P-MMSE: E{||v_i||^2} per UE (over cluster)
scaling_PM = zeros(K, 1);

for nreal = 1:nbrOfRealizations

    % LP-MMSE norm accumulation
    for l = 1:L
        servedUEs = find(D(l, :) == 1);
        if isempty(servedUEs), continue; end

        localMatrix = eye(N);
        for idx = 1:length(servedUEs)
            j = servedUEs(idx);
            hhat_jl = Hhat(:, nreal, l, j);
            localMatrix = localMatrix + p(j) * (hhat_jl * hhat_jl' + C(:, :, l, j));
        end

        for idx = 1:length(servedUEs)
            i = servedUEs(idx);
            hhat_il = Hhat(:, nreal, l, i);
            v_il = p(i) * (localMatrix \ hhat_il);
            scaling_LP(l, i) = scaling_LP(l, i) + real(v_il' * v_il) / nbrOfRealizations;
        end
    end

    % P-MMSE norm accumulation
    for i = 1:K
        servAPs = find(D(:, i) == 1);
        La = length(servAPs);
        if La == 0, continue; end

        Pi = false(K, 1);
        for l_idx = 1:La
            l = servAPs(l_idx);
            Pi(D(l, :) == 1) = true;
        end

        hhat_i = zeros(N * La, 1);
        for l_idx = 1:La
            l = servAPs(l_idx);
            hhat_i((l_idx-1)*N+1 : l_idx*N) = Hhat(:, nreal, l, i);
        end

        PMMSE_mat = eye(N * La);
        usersInPi = find(Pi);
        for idx = 1:length(usersInPi)
            j = usersInPi(idx);
            hhat_j = zeros(N * La, 1);
            Cj = zeros(N * La);
            for l_idx = 1:La
                l = servAPs(l_idx);
                hhat_j((l_idx-1)*N+1 : l_idx*N) = Hhat(:, nreal, l, j);
                Cj((l_idx-1)*N+1:l_idx*N, (l_idx-1)*N+1:l_idx*N) = C(:, :, l, j);
            end
            PMMSE_mat = PMMSE_mat + p(j) * (hhat_j * hhat_j' + Cj);
        end

        v_i = p(i) * (PMMSE_mat \ hhat_i);
        scaling_PM(i) = scaling_PM(i) + real(v_i' * v_i) / nbrOfRealizations;
    end
end

%% ====================================================================
%  PASS 2: Compute effective channel gains using properly normalized
%  precoding vectors
%  ====================================================================

gain_MR = zeros(nbrOfRealizations, K, K);
gain_LP = zeros(nbrOfRealizations, K, K);
gain_PM = zeros(nbrOfRealizations, K, K);

for nreal = 1:nbrOfRealizations

    %% MR Precoding: w_{il} = hhat_{il} * sqrt(rho_{il} / tr(B_{il}))
    for i = 1:K
        w_i = zeros(N * L, 1);
        for l = 1:L
            if D(l, i) == 1
                trB = real(trace(B(:, :, l, i)));
                if trB > 0
                    w_i((l-1)*N+1 : l*N) = Hhat(:, nreal, l, i) * sqrt(rho_dist(l, i) / trB);
                end
            end
        end
        for k = 1:K
            h_k = zeros(N * L, 1);
            for l = 1:L
                h_k((l-1)*N+1 : l*N) = H(:, nreal, l, k);
            end
            gain_MR(nreal, k, i) = h_k' * w_i;
        end
    end

    %% LP-MMSE Precoding: normalized by E{||v||^2}
    for i = 1:K
        w_i = zeros(N * L, 1);
        for l = 1:L
            if D(l, i) == 0, continue; end
            servedUEs = find(D(l, :) == 1);

            localMatrix = eye(N);
            for idx = 1:length(servedUEs)
                j = servedUEs(idx);
                hhat_jl = Hhat(:, nreal, l, j);
                localMatrix = localMatrix + p(j) * (hhat_jl * hhat_jl' + C(:, :, l, j));
            end

            hhat_il = Hhat(:, nreal, l, i);
            v_il = p(i) * (localMatrix \ hhat_il);

            % Normalize using EXPECTED norm (from Pass 1)
            if scaling_LP(l, i) > 0
                w_i((l-1)*N+1 : l*N) = v_il * sqrt(rho_dist(l, i) / scaling_LP(l, i));
            end
        end
        for k = 1:K
            h_k = zeros(N * L, 1);
            for l = 1:L
                h_k((l-1)*N+1 : l*N) = H(:, nreal, l, k);
            end
            gain_LP(nreal, k, i) = h_k' * w_i;
        end
    end

    %% P-MMSE Precoding: normalized by E{||v||^2}, centralized power
    for i = 1:K
        servAPs = find(D(:, i) == 1);
        La = length(servAPs);
        if La == 0, continue; end

        Pi = false(K, 1);
        for l_idx = 1:La
            l = servAPs(l_idx);
            Pi(D(l, :) == 1) = true;
        end

        hhat_i = zeros(N * La, 1);
        for l_idx = 1:La
            l = servAPs(l_idx);
            hhat_i((l_idx-1)*N+1 : l_idx*N) = Hhat(:, nreal, l, i);
        end

        PMMSE_mat = eye(N * La);
        usersInPi = find(Pi);
        for idx = 1:length(usersInPi)
            j = usersInPi(idx);
            hhat_j = zeros(N * La, 1);
            Cj = zeros(N * La);
            for l_idx = 1:La
                l = servAPs(l_idx);
                hhat_j((l_idx-1)*N+1 : l_idx*N) = Hhat(:, nreal, l, j);
                Cj((l_idx-1)*N+1:l_idx*N, (l_idx-1)*N+1:l_idx*N) = C(:, :, l, j);
            end
            PMMSE_mat = PMMSE_mat + p(j) * (hhat_j * hhat_j' + Cj);
        end

        v_i = p(i) * (PMMSE_mat \ hhat_i);

        % Normalize and apply centralized power
        w_i_full = zeros(N * L, 1);
        if scaling_PM(i) > 0
            for l_idx = 1:La
                l = servAPs(l_idx);
                w_i_full((l-1)*N+1 : l*N) = v_i((l_idx-1)*N+1 : l_idx*N) ...
                    * sqrt(rho_central(i) / scaling_PM(i));
            end
        end

        for k = 1:K
            h_k = zeros(N * L, 1);
            for l = 1:L
                h_k((l-1)*N+1 : l*N) = H(:, nreal, l, k);
            end
            gain_PM(nreal, k, i) = h_k' * w_i_full;
        end
    end
end

%% ====================================================================
%  Compute SE from effective channel gains
%  ====================================================================

[SE_MR, SE_MR_perfect] = computeDLSE(gain_MR, K, prelogFactor);
[SE_LP_MMSE, SE_LP_MMSE_perfect] = computeDLSE(gain_LP, K, prelogFactor);
[SE_P_MMSE, SE_P_MMSE_perfect] = computeDLSE(gain_PM, K, prelogFactor);

end


function [SE_hardening, SE_perfect] = computeDLSE(gains, K, prelogFactor)
%computeDLSE Compute DL SE using hardening bound and perfect-CSI bound.

nR = size(gains, 1);
SE_hardening = zeros(K, 1);
SE_perfect = zeros(K, 1);

% Perfect CSI (instantaneous SINR per realization)
SE_inst = zeros(nR, K);
for nreal = 1:nR
    for k = 1:K
        sigPow = abs(gains(nreal, k, k))^2;
        intPow = sum(abs(gains(nreal, k, :)).^2) - sigPow;
        SE_inst(nreal, k) = prelogFactor * log2(1 + sigPow / max(intPow + 1, 1e-16));
    end
end
SE_perfect = mean(SE_inst, 1)';

% Hardening bound (Proposition 3)
for k = 1:K
    meanDes = mean(gains(:, k, k));
    denom = 0;
    for i = 1:K
        denom = denom + mean(abs(gains(:, k, i)).^2);
    end
    denom = denom - abs(meanDes)^2 + 1;
    SE_hardening(k) = prelogFactor * log2(1 + abs(meanDes)^2 / max(denom, 1e-16));
end

end
