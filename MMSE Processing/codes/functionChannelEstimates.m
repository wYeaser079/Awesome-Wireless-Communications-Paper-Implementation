function [Hhat, H, B, C] = functionChannelEstimates(R, nbrOfRealizations, L, K, N, tau_p, pilotIndex, p)
%functionChannelEstimates Generate channel realizations and compute MMSE
%channel estimates for cell-free mMIMO.
%
%   [Hhat, H, B, C] = functionChannelEstimates(R, nbrOfRealizations, L, K, N, tau_p, pilotIndex, p)
%
%   Inputs:
%       R                 - N x N x L x K spatial correlation matrices
%       nbrOfRealizations - Number of channel realizations
%       L                 - Number of APs
%       K                 - Number of UEs
%       N                 - Number of antennas per AP
%       tau_p             - Number of orthogonal pilots
%       pilotIndex        - K x 1 vector of pilot indices (1 to tau_p)
%       p                 - K x 1 vector of UE transmit powers (W)
%
%   Outputs:
%       Hhat - N x nbrOfRealizations x L x K MMSE channel estimates
%       H    - N x nbrOfRealizations x L x K true channel realizations
%       B    - N x N x L x K estimation error covariance (C_kl)
%       C    - N x N x L x K estimate covariance
%
% Reference:
% Emil Bjornson, Luca Sanguinetti, "Making Cell-Free Massive MIMO
% Competitive With MMSE Processing and Centralized Implementation,"
% IEEE Trans. Wireless Commun., vol. 19, no. 2, pp. 1257-1271, Feb. 2020.
% Equations (4)-(6).

% Noise variance normalized to 1 (since R already contains gain-over-noise)
sigma2 = 1;

% Initialize output arrays
Hhat = zeros(N, nbrOfRealizations, L, K);
H = zeros(N, nbrOfRealizations, L, K);
B = zeros(N, N, L, K);   % Estimation error covariance C_kl
C = zeros(N, N, L, K);   % Estimate covariance

%% Generate channel realizations
for l = 1:L
    for k = 1:K
        Rkl = R(:, :, l, k);

        % Generate i.i.d. CN(0,I) samples
        W = (randn(N, nbrOfRealizations) + 1j * randn(N, nbrOfRealizations)) / sqrt(2);

        % Apply spatial correlation: h_kl = sqrtm(R_kl) * w
        if N == 1
            H(:, :, l, k) = sqrt(Rkl) * W;
        else
            % Use Cholesky for numerical stability
            try
                sqrtR = chol(Rkl + 1e-15 * eye(N), 'lower');
            catch
                [U, D] = eig(Rkl);
                sqrtR = U * sqrt(max(D, 0));
            end
            H(:, :, l, k) = sqrtR * W;
        end
    end
end

%% Compute MMSE channel estimates
for l = 1:L

    % Group pilots: for each pilot index t, find which UEs use it
    for t = 1:tau_p

        % Find UEs using pilot t
        usersWithPilot = find(pilotIndex == t);

        if isempty(usersWithPilot)
            continue;
        end

        % Compute Psi_tl = sum_{i in P_k} p_i * tau_p * R_il + sigma^2 * I_N (Eq. 5)
        Psi_tl = sigma2 * eye(N);
        for i = usersWithPilot'
            Psi_tl = Psi_tl + p(i) * tau_p * R(:, :, l, i);
        end

        % Invert Psi
        Psi_tl_inv = inv(Psi_tl);

        % Generate pilot observation z_tl (Eq. 3)
        z_tl = zeros(N, nbrOfRealizations);
        for i = usersWithPilot'
            z_tl = z_tl + sqrt(p(i) * tau_p) * H(:, :, l, i);
        end
        % Add noise
        z_tl = z_tl + (randn(N, nbrOfRealizations) + 1j * randn(N, nbrOfRealizations)) / sqrt(2);

        % Compute MMSE estimates for each UE using pilot t (Eq. 4)
        for k = usersWithPilot'
            Rkl = R(:, :, l, k);

            % MMSE estimate: hhat_kl = sqrt(p_k * tau_p) * R_kl * Psi^{-1} * z_tl
            A_kl = sqrt(p(k) * tau_p) * Rkl * Psi_tl_inv;
            Hhat(:, :, l, k) = A_kl * z_tl;

            % Estimate covariance: p_k * tau_p * R_kl * Psi^{-1} * R_kl
            C(:, :, l, k) = p(k) * tau_p * Rkl * Psi_tl_inv * Rkl;

            % Estimation error covariance: R_kl - C_kl (Eq. 6)
            B(:, :, l, k) = Rkl - C(:, :, l, k);
        end
    end
end

end
