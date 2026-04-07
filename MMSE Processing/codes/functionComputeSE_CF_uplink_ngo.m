function SE = functionComputeSE_CF_uplink_ngo(gainOverNoisedB, tau_c, tau_p, pilotIndex, p, K, L)
%functionComputeSE_CF_uplink_ngo Compute closed-form uplink SE for
%cell-free mMIMO with MR combining at Level 2 (N=1).
%
%   SE = functionComputeSE_CF_uplink_ngo(gainOverNoisedB, tau_c, tau_p, pilotIndex, p, K, L)
%
%   Based on the approach in Ngo et al. [16] for Level 2 with MR combining
%   and single-antenna APs.
%
%   Inputs:
%       gainOverNoisedB - L x K matrix of gain-over-noise in dB
%       tau_c           - Coherence block length
%       tau_p           - Pilot length
%       pilotIndex      - K x 1 pilot assignment vector
%       p               - K x 1 transmit powers
%       K               - Number of UEs
%       L               - Number of APs
%
%   Output:
%       SE - K x 1 vector of per-user SE (bit/s/Hz)
%
% Reference:
% H. Q. Ngo et al., "Cell-Free Massive MIMO Versus Small Cells," IEEE TWC, 2017.

sigma2 = 1;
prelogFactor = (tau_c - tau_p) / tau_c;

% Convert to linear scale
beta = db2pow(gainOverNoisedB);  % L x K

SE = zeros(K, 1);

for k = 1:K
    % Pilot sharing set
    Pk = find(pilotIndex == pilotIndex(k));

    % Compute gamma_kl for all APs (MMSE estimation quality)
    gamma_k = zeros(L, 1);
    for l = 1:L
        Psi_tl = sigma2;
        for i = Pk'
            Psi_tl = Psi_tl + p(i) * tau_p * beta(l, i);
        end
        gamma_k(l) = p(k) * tau_p * beta(l, k)^2 / Psi_tl;
    end

    % Level 2 SE with MR combining (closed-form from [16])
    % Numerator: p_k * (sum_l gamma_kl)^2
    numerator = p(k) * (sum(gamma_k))^2;

    % Denominator terms
    denom = 0;

    % Term 1: Interference from co-pilot UEs (coherent)
    for i = Pk'
        if i ~= k
            gamma_i = zeros(L, 1);
            for l = 1:L
                Psi_tl = sigma2;
                for j2 = Pk'
                    Psi_tl = Psi_tl + p(j2) * tau_p * beta(l, j2);
                end
                gamma_i(l) = p(i) * tau_p * beta(l, i)^2 / Psi_tl;
            end
            % Coherent contamination term
            % sum_l sqrt(p_i * gamma_il / gamma_kl * gamma_kl) type terms
            cross_term = 0;
            for l = 1:L
                if gamma_k(l) > 0
                    cross_term = cross_term + sqrt(p(i) / p(k)) * beta(l, i) / beta(l, k) * gamma_k(l);
                end
            end
            denom = denom + p(i) * cross_term^2;
        end
    end

    % Term 2: Non-coherent interference from all UEs
    for i = 1:K
        for l = 1:L
            denom = denom + p(i) * gamma_k(l) * beta(l, i);
        end
    end

    % Term 3: Noise amplification
    for l = 1:L
        denom = denom + sigma2 * gamma_k(l);
    end

    SE(k) = prelogFactor * log2(1 + numerator / max(denom, 1e-20));
end

end
