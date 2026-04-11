function [SE_MR, SE_MMMSE] = functionComputeSE_BS_uplink(...
    R_BS, gainOverNoisedB_BS, tau_c, tau_p, nbrOfRealizations, ...
    nbrOfBSs, M, K, Kc, p_cellular)
%functionComputeSE_BS_uplink Compute uplink SE for cellular Massive MIMO
%with MR and M-MMSE combining. Vectorised implementation.
%
%   [SE_MR, SE_MMMSE] = functionComputeSE_BS_uplink(...)
%
% Reference: Bjornson & Sanguinetti, IEEE TWC, 2020, Eqs. (31)-(36).

sigma2 = 1;
prelogFactor = (tau_c - tau_p) / tau_c;

% Assign UEs to cells based on strongest average channel
UEcellAssignment = zeros(K, 1);
for k = 1:K
    avgGain = zeros(nbrOfBSs, 1);
    for j = 1:nbrOfBSs
        avgGain(j) = gainOverNoisedB_BS(j, k, j);
    end
    [~, UEcellAssignment(k)] = max(avgGain);
end

UEsInCell = cell(nbrOfBSs, 1);
for j = 1:nbrOfBSs
    UEsInCell{j} = find(UEcellAssignment == j);
end

SE_MR    = zeros(K, 1);
SE_MMMSE = zeros(K, 1);

for j = 1:nbrOfBSs
    localUEs = UEsInCell{j};
    nbrLocalUEs = length(localUEs);
    if nbrLocalUEs == 0; continue; end

    % All UE indices (ordered cell-by-cell)
    allUEindices = [];
    for l = 1:nbrOfBSs
        allUEindices = [allUEindices; UEsInCell{l}]; %#ok<AGROW>
    end
    nbrAllUEs = length(allUEindices);

    % Generate channel realizations from every UE to BS j
    H_j = zeros(M, nbrOfRealizations, nbrAllUEs);
    for idx = 1:nbrAllUEs
        k = allUEindices(idx);
        cellOfUE = UEcellAssignment(k);
        Rjk = R_BS(:, :, j, k, cellOfUE);
        W = (randn(M, nbrOfRealizations) + 1j*randn(M, nbrOfRealizations)) / sqrt(2);
        try
            sqrtR = chol(Rjk + 1e-15*eye(M), 'lower');
        catch
            [U, D] = eig(Rjk);
            sqrtR = U * sqrt(max(D, 0));
        end
        H_j(:, :, idx) = sqrtR * W;
    end

    % MMSE channel estimation with pilot contamination
    Hhat_j = zeros(M, nbrOfRealizations, nbrAllUEs);
    C_j    = zeros(M, M, nbrAllUEs);

    for localIdx = 1:nbrLocalUEs
        k = localUEs(localIdx);
        localIdxInCell = find(UEsInCell{j} == k);

        Psi = sigma2 * eye(M);
        coPilotGlobalIdx = [];
        for l = 1:nbrOfBSs
            if localIdxInCell <= length(UEsInCell{l})
                coPilotUE = UEsInCell{l}(localIdxInCell);
                coPilotGlobalIdx = [coPilotGlobalIdx; find(allUEindices == coPilotUE)]; %#ok<AGROW>
                cellOfCoPilot = UEcellAssignment(coPilotUE);
                Psi = Psi + p_cellular * tau_p * R_BS(:, :, j, coPilotUE, cellOfCoPilot);
            end
        end
        Psi_inv = inv(Psi);

        z_pilot = zeros(M, nbrOfRealizations);
        for cpIdx = coPilotGlobalIdx'
            z_pilot = z_pilot + sqrt(p_cellular * tau_p) * H_j(:, :, cpIdx);
        end
        z_pilot = z_pilot + (randn(M, nbrOfRealizations) + 1j*randn(M, nbrOfRealizations)) / sqrt(2);

        for cpIdx = coPilotGlobalIdx'
            cpUE = allUEindices(cpIdx);
            cellOfCp = UEcellAssignment(cpUE);
            Rjcp = R_BS(:, :, j, cpUE, cellOfCp);
            A = sqrt(p_cellular * tau_p) * Rjcp * Psi_inv;
            Hhat_j(:, :, cpIdx) = A * z_pilot;
            C_j(:, :, cpIdx) = p_cellular * tau_p * Rjcp * Psi_inv * Rjcp;
        end
    end

    % Error + noise covariance (does not depend on n or k)
    E_j = sigma2 * eye(M);
    for idx = 1:nbrAllUEs
        k = allUEindices(idx);
        cellOfUE = UEcellAssignment(k);
        Bj_idx = R_BS(:, :, j, k, cellOfUE) - C_j(:, :, idx);
        E_j = E_j + p_cellular * Bj_idx;
    end

    % Indices within allUEindices that belong to cell j's local UEs
    localGlobalIdx = zeros(nbrLocalUEs, 1);
    for li = 1:nbrLocalUEs
        localGlobalIdx(li) = find(allUEindices == localUEs(li));
    end

    % Per-UE SE accumulators
    SE_MR_local    = zeros(nbrLocalUEs, 1);
    SE_MMMSE_local = zeros(nbrLocalUEs, 1);

    bsProgressStep = max(1, floor(nbrOfRealizations / 10));
    for n = 1:nbrOfRealizations
        if mod(n, bsProgressStep) == 0 || n == nbrOfRealizations
            fprintf('    [BS%d] %d/%d\n', j, n, nbrOfRealizations);
        end
        % Stack all estimates at BS j for this realization: M x nbrAllUEs
        Hh = squeeze(Hhat_j(:, n, :));
        if size(Hh, 2) == 1 && nbrAllUEs == 1
            Hh = Hh(:);  % keep as column
        end

        % Sigma is independent of k — compute once per realization
        Sigma = E_j + p_cellular * (Hh * Hh');
        Sigma_inv = inv(Sigma);

        % Inner products for MR: G = Hh' * Hh  (nbrAllUEs x nbrAllUEs)
        G = Hh' * Hh;
        absG2 = abs(G).^2;

        for li = 1:nbrLocalUEs
            kIdx = localGlobalIdx(li);
            hhat_k = Hh(:, kIdx);

            %--- MR ---
            sig = p_cellular * abs(G(kIdx, kIdx))^2;
            interf = p_cellular * (sum(absG2(kIdx, :)) - absG2(kIdx, kIdx));
            noise_term = real(hhat_k' * E_j * hhat_k);
            SE_MR_local(li) = SE_MR_local(li) + log2(1 + real(sig) / max(real(interf + noise_term), 1e-20));

            %--- M-MMSE (matrix-inversion lemma form) ---
            t = real(hhat_k' * (Sigma_inv * hhat_k));
            denom = max(1 - p_cellular * t, 1e-12);
            SINR_val = max(p_cellular * t / denom, 0);
            SE_MMMSE_local(li) = SE_MMMSE_local(li) + log2(1 + SINR_val);
        end
    end

    for li = 1:nbrLocalUEs
        k = localUEs(li);
        SE_MR(k)    = prelogFactor * SE_MR_local(li)    / nbrOfRealizations;
        SE_MMMSE(k) = prelogFactor * SE_MMMSE_local(li) / nbrOfRealizations;
    end
end

end
