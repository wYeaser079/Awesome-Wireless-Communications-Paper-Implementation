function [signalCF, interferenceCF, noiseCF, signalSC, interferenceSC, noiseSC] = ...
    functionSINRterms_uplink_ngo(gainOverNoisedB, tau_p, pilotIndex, K, L)
%functionSINRterms_uplink_ngo Extract SINR terms (without power factors)
%for cell-free and small-cell configurations with MR combining (N=1).
%
%   Used for max-min power optimization in Figure 4.
%
%   Inputs:
%       gainOverNoisedB - L x K gain-over-noise in dB
%       tau_p           - Pilot length
%       pilotIndex      - K x 1 pilot indices
%       K, L            - Number of UEs and APs
%
%   Outputs (all without power scaling):
%       signalCF       - K x 1 signal term for cell-free
%       interferenceCF - K x K interference matrix for cell-free
%       noiseCF        - K x 1 noise term for cell-free
%       signalSC       - K x L signal term for small cells
%       interferenceSC - K x K x L interference matrix for small cells
%       noiseSC        - K x L noise term for small cells

sigma2 = 1;
beta = db2pow(gainOverNoisedB);

% Compute gamma_kl for all k, l
gamma = zeros(L, K);
for l = 1:L
    for k = 1:K
        Pk = find(pilotIndex == pilotIndex(k));
        Psi_tl = sigma2;
        for i = Pk'
            Psi_tl = Psi_tl + tau_p * beta(l, i);  % p=1 reference
        end
        gamma(l, k) = tau_p * beta(l, k)^2 / Psi_tl;
    end
end

%% Cell-free terms (Level 2, MR, uniform weights)
signalCF = zeros(K, 1);
interferenceCF = zeros(K, K);
noiseCF = zeros(K, 1);

for k = 1:K
    % Signal: (sum_l gamma_kl)^2
    signalCF(k) = sum(gamma(:, k))^2;

    % Noise
    noiseCF(k) = sigma2 * sum(gamma(:, k));

    % Interference
    Pk = find(pilotIndex == pilotIndex(k));
    for i = 1:K
        if ismember(i, Pk) && i ~= k
            % Coherent contamination
            cross = 0;
            for l = 1:L
                cross = cross + beta(l, i) / (beta(l, k) + 1e-30) * gamma(l, k);
            end
            interferenceCF(k, i) = cross^2;
        end
        % Non-coherent interference
        for l = 1:L
            interferenceCF(k, i) = interferenceCF(k, i) + gamma(l, k) * beta(l, i);
        end
    end
end

%% Small-cell terms (Level 1, per AP)
signalSC = zeros(K, L);
interferenceSC = zeros(K, K, L);
noiseSC = zeros(K, L);

for k = 1:K
    Pk = find(pilotIndex == pilotIndex(k));
    for l = 1:L
        signalSC(k, l) = gamma(l, k);

        noiseSC(k, l) = sigma2;

        for i = 1:K
            if ismember(i, Pk) && i ~= k
                interferenceSC(k, i, l) = gamma(l, i);
            end
            interferenceSC(k, i, l) = interferenceSC(k, i, l) + beta(l, i);
        end
    end
end

end
