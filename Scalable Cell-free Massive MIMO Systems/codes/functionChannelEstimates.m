function [Hhat, H, B, C] = functionChannelEstimates(R, nbrOfRealizations, L, K, N, tau_p, pilotIndex, p)
%functionChannelEstimates Generate channel realizations and compute MMSE
%channel estimates for Cell-Free Massive MIMO.
%
%   [Hhat, H, B, C] = functionChannelEstimates(R, nbrOfRealizations, ...
%       L, K, N, tau_p, pilotIndex, p)
%
%   Generates i.i.d. correlated Rayleigh fading channel realizations and
%   computes the MMSE channel estimates using Eq. (3)-(4) from [1].
%
%   INPUTS:
%       R                  - N x N x L x K spatial correlation matrices
%                            (normalized by noise power)
%       nbrOfRealizations  - Number of channel realizations
%       L                  - Number of APs
%       K                  - Number of UEs
%       N                  - Number of antennas per AP
%       tau_p              - Number of orthogonal pilot sequences
%       pilotIndex         - K x 1 vector of pilot assignments
%       p                  - K x 1 vector of UL transmit powers
%                            (normalized by noise power)
%
%   OUTPUTS:
%       Hhat - N x nbrOfRealizations x L x K array of MMSE channel estimates
%       H    - N x nbrOfRealizations x L x K array of true channel realizations
%       B    - N x N x L x K array of estimate correlation matrices
%       C    - N x N x L x K array of estimation error correlation matrices
%
%   REFERENCES:
%   [1] E. Bjornson, L. Sanguinetti, "Scalable Cell-Free Massive MIMO
%       Systems," IEEE Trans. Commun., vol. 68, no. 7, pp. 4247-4261, 2020.

%% Pre-allocate arrays
Hhat = zeros(N, nbrOfRealizations, L, K);
H = zeros(N, nbrOfRealizations, L, K);
B = zeros(N, N, L, K);
C = zeros(N, N, L, K);

% Ensure p is a column vector
if size(p, 2) > size(p, 1)
    p = p.';
end

%% Generate channel realizations and MMSE estimates
for l = 1:L

    % --- Compute MMSE estimation matrices per pilot ---
    % Group UEs by their pilot assignment
    % For each pilot t, compute Psi_tl = sum_{i in S_t} tau_p*p_i*R_il + I_N
    PsiInv = zeros(N, N, tau_p);
    for t = 1:tau_p
        Psi_tl = eye(N); % sigma^2 * I_N (already normalized: sigma^2 = 1)
        usersOnPilot = find(pilotIndex == t);
        for idx = 1:length(usersOnPilot)
            i = usersOnPilot(idx);
            Psi_tl = Psi_tl + tau_p * p(i) * R(:, :, l, i);
        end
        PsiInv(:, :, t) = inv(Psi_tl); % Psi_{tl}^{-1}
    end

    % --- Generate channels and estimates for each UE ---
    for k = 1:K
        % Get the correlation matrix for AP l, UE k
        Rkl = R(:, :, l, k);

        % Compute the square root of Rkl for channel generation
        % h_{kl} = sqrtm(R_{kl}) * w, where w ~ CN(0, I_N)
        if N == 1
            sqrtR = sqrt(Rkl);
        else
            sqrtR = sqrtm(Rkl);
        end

        % Generate i.i.d. complex Gaussian samples
        W = (randn(N, nbrOfRealizations) + 1j * randn(N, nbrOfRealizations)) / sqrt(2);

        % True channel realization: h_{kl} ~ CN(0, R_{kl})
        H(:, :, l, k) = sqrtR * W;

        % Compute MMSE estimate correlation matrix B_{kl}
        t_k = pilotIndex(k); % Pilot index of UE k
        B(:, :, l, k) = p(k) * tau_p * Rkl * PsiInv(:, :, t_k) * Rkl;

        % Estimation error correlation matrix C_{kl} = R_{kl} - B_{kl}
        C(:, :, l, k) = Rkl - B(:, :, l, k);
    end

    % --- Generate MMSE channel estimates ---
    % The received pilot signal at AP l for pilot t is:
    %   y^{pilot}_{tl} = sum_{i in S_t} sqrt(tau_p * p_i) * h_{il} + n_{tl}
    % The MMSE estimate is:
    %   hhat_{kl} = sqrt(p_k * tau_p) * R_{kl} * Psi_{tl}^{-1} * y^{pilot}_{tl}

    for t = 1:tau_p
        % Construct the received pilot signal for pilot t at AP l
        usersOnPilot = find(pilotIndex == t);

        % Received pilot signal: y^{pilot}_{tl}
        % Dimensions: N x nbrOfRealizations
        yPilot = (randn(N, nbrOfRealizations) + 1j * randn(N, nbrOfRealizations)) / sqrt(2);

        for idx = 1:length(usersOnPilot)
            i = usersOnPilot(idx);
            yPilot = yPilot + sqrt(tau_p * p(i)) * H(:, :, l, i);
        end

        % Compute MMSE estimates for all UEs on pilot t
        for idx = 1:length(usersOnPilot)
            k = usersOnPilot(idx);
            Rkl = R(:, :, l, k);

            % hhat_{kl} = sqrt(p_k * tau_p) * R_{kl} * Psi_{tl}^{-1} * y^{pilot}_{tl}
            Hhat(:, :, l, k) = sqrt(p(k) * tau_p) * Rkl * PsiInv(:, :, pilotIndex(k)) * yPilot;
        end
    end

end

end
