function [SE_maxBeta, SE_maxSE] = functionComputeSE_SC_uplink_ngo(...
    gainOverNoisedB, tau_c, tau_p, pilotIndex, p, K, L)
%functionComputeSE_SC_uplink_ngo Compute closed-form uplink SE for small
%cells (Level 1) with N=1 using Proposition 3 (Eq. 29).
%
%   Two AP selection strategies:
%   1. Max-beta: UE connects to AP with strongest large-scale fading (Ngo's approach)
%   2. Max-SE: UE connects to AP giving highest SE (this paper's improvement)
%
%   Inputs:
%       gainOverNoisedB - L x K gain-over-noise in dB
%       tau_c           - Coherence block length
%       tau_p           - Pilot length
%       pilotIndex      - K x 1 pilot indices
%       p               - K x 1 transmit powers
%       K               - Number of UEs
%       L               - Number of APs
%
%   Outputs:
%       SE_maxBeta - K x 1 SE using max-beta AP selection
%       SE_maxSE   - K x 1 SE using max-SE AP selection
%
% Reference:
% Emil Bjornson, Luca Sanguinetti, "Making Cell-Free Massive MIMO
% Competitive With MMSE Processing and Centralized Implementation,"
% IEEE Trans. Wireless Commun., 2020. Proposition 3, Eq. (29)-(30).

sigma2 = 1;
prelogFactor = (tau_c - tau_p) / tau_c;

% Convert to linear scale
beta = db2pow(gainOverNoisedB);  % L x K

SE_maxBeta = zeros(K, 1);
SE_maxSE = zeros(K, 1);

for k = 1:K
    Pk = find(pilotIndex == pilotIndex(k));  % Pilot sharing set

    SE_perAP = zeros(L, 1);

    for l = 1:L
        % Compute gamma_kl (Eq. 30)
        Psi_tl = sigma2;
        for i = Pk'
            Psi_tl = Psi_tl + p(i) * tau_p * beta(l, i);
        end
        gamma_kl = p(k)^2 * tau_p * beta(l, k)^2 / Psi_tl;

        % Pilot contamination ratio A_kl
        A_kl = 0;
        for i = Pk'
            if i ~= k
                gamma_il = p(i)^2 * tau_p * beta(l, i)^2 / Psi_tl;
                A_kl = A_kl + p(i) * gamma_il / (p(k) * gamma_kl + 1e-30);
            end
        end

        % Effective SINR denominator terms (non-copilot interference + noise)
        denom_extra = 0;
        for i = 1:K
            if ~ismember(i, Pk)
                C_il = beta(l, i) - p(i) * tau_p * beta(l, i)^2 / Psi_tl;
                denom_extra = denom_extra + p(i) * C_il;
            else
                % Co-pilot UEs: estimation error
                C_il = beta(l, i) - p(i) * tau_p * beta(l, i)^2 / Psi_tl;
                if i ~= k
                    denom_extra = denom_extra + p(i) * C_il;
                end
            end
        end
        denom_extra = denom_extra + sigma2;

        % Compute gamma_kl scaled for SE formula
        gamma_eff = gamma_kl / denom_extra;

        if gamma_eff < 1e-15
            SE_perAP(l) = 0;
            continue;
        end

        % SE using Proposition 3 (Eq. 29) with E1 function
        % SE = (1/ln2) * [e^(1/(gamma*(1+A))) * E1(1/(gamma*(1+A)))
        %                 - e^(1/(gamma*A)) * E1(1/(gamma*A))]
        % where gamma = gamma_eff, A = A_kl

        if A_kl < 1e-10
            % No pilot contamination: simplified formula
            % SE = (1/ln2) * e^(1/gamma) * E1(1/gamma)
            x = 1 / gamma_eff;
            SE_perAP(l) = prelogFactor * computeE1term(x);
        else
            x1 = 1 / (gamma_eff * (1 + A_kl));
            x2 = 1 / (gamma_eff * A_kl);
            SE_perAP(l) = prelogFactor * (computeE1term(x1) - computeE1term(x2));
        end
    end

    % Max-beta AP selection (Ngo's approach)
    [~, bestAP_beta] = max(beta(:, k));
    SE_maxBeta(k) = SE_perAP(bestAP_beta);

    % Max-SE AP selection (this paper's improvement)
    SE_maxSE(k) = max(SE_perAP);
end

end


function val = computeE1term(x)
%computeE1term Compute (1/ln(2)) * exp(x) * E1(x) with numerical stability
%
%   For small x: exp(x)*E1(x) ≈ 1/x (use bound)
%   For large x: use direct computation

if x > 500
    % Use asymptotic approximation: e^x * E1(x) ≈ 1/x for large x
    val = (1/x) / log(2);
elseif x < 1e-10
    % Use bound: e^x * E1(x) ≈ -log(x) - euler_gamma for small x
    % More precisely: use the tight bound x/(1+x) <= e^x*E1(x) <= x for 1/x
    % Actually for small x (i.e., large gamma), e^x*E1(x) ≈ -log(x) - 0.5772
    val = (-log(x) - 0.5772156649) / log(2);
else
    % Direct computation using MATLAB's expint function
    % E1(x) = expint(x) in MATLAB
    val = exp(x) * expint(x) / log(2);
end

end
