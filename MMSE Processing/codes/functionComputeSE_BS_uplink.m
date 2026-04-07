function [SE_MR, SE_MMMSE] = functionComputeSE_BS_uplink(...
    R_BS, gainOverNoisedB_BS, tau_c, tau_p, nbrOfRealizations, ...
    nbrOfBSs, M, K, Kc, p_cellular)
%functionComputeSE_BS_uplink Compute uplink SE for cellular Massive MIMO
%with MR and M-MMSE combining.
%
%   [SE_MR, SE_MMMSE] = functionComputeSE_BS_uplink(...)
%
%   Inputs:
%       R_BS              - M x M x nbrOfBSs x K x nbrOfBSs correlation matrices
%       gainOverNoisedB_BS- nbrOfBSs x K x nbrOfBSs gain-over-noise
%       tau_c             - Coherence block length
%       tau_p             - Pilot length (= Kc for cellular)
%       nbrOfRealizations - Number of channel realizations
%       nbrOfBSs          - Number of BSs
%       M                 - Number of antennas per BS
%       K                 - Total number of UEs
%       Kc                - Number of UEs per cell
%       p_cellular        - Transmit power per UE (scalar, Watts)
%
%   Outputs:
%       SE_MR    - K x 1 SE with MR combining
%       SE_MMMSE - K x 1 SE with M-MMSE combining
%
% Reference:
% Emil Bjornson, Luca Sanguinetti, "Making Cell-Free Massive MIMO
% Competitive With MMSE Processing and Centralized Implementation,"
% IEEE Trans. Wireless Commun., vol. 19, no. 2, pp. 1257-1271, Feb. 2020.
% Equations (31)-(36).

sigma2 = 1;
prelogFactor = (tau_c - tau_p) / tau_c;

% Assign UEs to cells based on strongest average channel
UEcellAssignment = zeros(K, 1);
for k = 1:K
    avgGain = zeros(nbrOfBSs, 1);
    for j = 1:nbrOfBSs
        avgGain(j) = gainOverNoisedB_BS(j, k, j);  % Gain from BS j to UE k
    end
    [~, UEcellAssignment(k)] = max(avgGain);
end

% Organize UEs per cell
UEsInCell = cell(nbrOfBSs, 1);
for j = 1:nbrOfBSs
    UEsInCell{j} = find(UEcellAssignment == j);
end

% Initialize outputs
SE_MR = zeros(K, 1);
SE_MMMSE = zeros(K, 1);

% Process each BS
for j = 1:nbrOfBSs

    localUEs = UEsInCell{j};
    nbrLocalUEs = length(localUEs);

    if nbrLocalUEs == 0
        continue;
    end

    %% Generate channels and estimates for BS j
    % Channels from ALL UEs in ALL cells to BS j
    allUEindices = [];
    for l = 1:nbrOfBSs
        allUEindices = [allUEindices; UEsInCell{l}]; %#ok<AGROW>
    end

    nbrAllUEs = length(allUEindices);

    H_j = zeros(M, nbrOfRealizations, nbrAllUEs);
    Hhat_j = zeros(M, nbrOfRealizations, nbrAllUEs);
    C_j = zeros(M, M, nbrAllUEs);

    for idx = 1:nbrAllUEs
        k = allUEindices(idx);
        cellOfUE = UEcellAssignment(k);

        Rjk = R_BS(:, :, j, k, cellOfUE);

        % Generate channel realizations
        W = (randn(M, nbrOfRealizations) + 1j * randn(M, nbrOfRealizations)) / sqrt(2);
        try
            sqrtR = chol(Rjk + 1e-15*eye(M), 'lower');
        catch
            [U, D] = eig(Rjk);
            sqrtR = U * sqrt(max(D, 0));
        end
        H_j(:, :, idx) = sqrtR * W;
    end

    % MMSE channel estimation with pilot contamination
    % In cellular, pilot reuse factor 1: UE k in cell j uses same pilot as
    % UE k in cell l (for same local index within cell)
    for localIdx = 1:nbrLocalUEs
        k = localUEs(localIdx);

        % Find co-pilot UEs (same local index in other cells)
        localIdxInCell = find(UEsInCell{j} == k);

        % Compute Psi for pilot of this UE
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

        % Generate pilot observation
        z_pilot = zeros(M, nbrOfRealizations);
        for cpIdx = coPilotGlobalIdx'
            z_pilot = z_pilot + sqrt(p_cellular * tau_p) * H_j(:, :, cpIdx);
        end
        z_pilot = z_pilot + (randn(M, nbrOfRealizations) + 1j*randn(M, nbrOfRealizations))/sqrt(2);

        % Compute estimates for all co-pilot UEs at BS j
        for cpIdx = coPilotGlobalIdx'
            cpUE = allUEindices(cpIdx);
            cellOfCp = UEcellAssignment(cpUE);
            Rjcp = R_BS(:, :, j, cpUE, cellOfCp);
            A = sqrt(p_cellular * tau_p) * Rjcp * Psi_inv;
            Hhat_j(:, :, cpIdx) = A * z_pilot;
            C_j(:, :, cpIdx) = p_cellular * tau_p * Rjcp * Psi_inv * Rjcp;
        end
    end

    % Estimation error covariance
    B_j = zeros(M, M, nbrAllUEs);
    for idx = 1:nbrAllUEs
        k = allUEindices(idx);
        cellOfUE = UEcellAssignment(k);
        B_j(:, :, idx) = R_BS(:, :, j, k, cellOfUE) - C_j(:, :, idx);
    end

    %% Compute SE for UEs in cell j
    % Error + noise covariance
    E_j = sigma2 * eye(M);
    for idx = 1:nbrAllUEs
        E_j = E_j + p_cellular * B_j(:, :, idx);
    end

    for localIdx = 1:nbrLocalUEs
        k = localUEs(localIdx);
        kIdx = find(allUEindices == k);

        SINR_MR_sum = 0;
        SINR_MMMSE_sum = 0;

        for n = 1:nbrOfRealizations
            hhat_k = Hhat_j(:, n, kIdx);

            % Build M-MMSE matrix: Sigma = sum p_i * hhat_i * hhat_i^H + E
            Sigma = E_j;
            for idx = 1:nbrAllUEs
                Sigma = Sigma + p_cellular * (Hhat_j(:, n, idx) * Hhat_j(:, n, idx)');
            end
            Sigma_inv = inv(Sigma);

            %--- MR combining ---
            v_MR = hhat_k;
            if norm(v_MR) > 1e-15
                sig = p_cellular * abs(v_MR' * hhat_k)^2;
                interf = 0;
                for idx = 1:nbrAllUEs
                    if idx ~= kIdx
                        interf = interf + p_cellular * abs(v_MR' * Hhat_j(:, n, idx))^2;
                    end
                end
                noise_term = real(v_MR' * E_j * v_MR);
                SINR_MR_sum = SINR_MR_sum + log2(1 + real(sig) / max(real(interf + noise_term), 1e-20));
            end

            %--- M-MMSE combining ---
            temp = Sigma_inv * hhat_k;
            SINR_val = p_cellular * real(hhat_k' * temp) / (1 - p_cellular * real(hhat_k' * temp));
            SINR_val = max(SINR_val, 0);
            SINR_MMMSE_sum = SINR_MMMSE_sum + log2(1 + SINR_val);
        end

        SE_MR(k) = prelogFactor * SINR_MR_sum / nbrOfRealizations;
        SE_MMMSE(k) = prelogFactor * SINR_MMMSE_sum / nbrOfRealizations;
    end
end

end
