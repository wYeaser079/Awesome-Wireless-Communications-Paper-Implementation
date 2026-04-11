function [SE_k, SINR_k] = compute_SE_closedform(eta_mk, beta_mk, gamma_mk, ...
                                                pilotseq, N, rho_d, tau_c, tau_p)
%COMPUTE_SE_CLOSEDFORM  Closed-form downlink SE of every user in a
%  cell-free massive-MIMO system with conjugate beamforming and
%  multi-antenna APs.  Implements Proposition 1 / Eq. (15) of
%  Ngo-Tran-Duong-Matthaiou-Larsson IEEE TGCN 2018:
%
%      SE_k = (1 - tau_p/tau_c) * log2(1 + SINR_k)
%
%  with
%      SINR_k =  rho_d N^2 | gamma_{k,k}^T sqrt(eta_k) |^2
%               / ( rho_d N^2 sum_{k'~=k} | gamma_{k,k'}^T sqrt(eta_{k'}) |^2
%                 + rho_d N  sum_{k'}   || D_{k',k} sqrt(eta_{k'}) ||^2
%                 + 1 )
%
%  where
%      [gamma_{k,k'}]_m   = |phi_{k'}^H phi_k| * gamma_{m,k'} * beta_{m,k}/beta_{m,k'}
%      [D_{k',k}]_{m,m}   = sqrt( beta_{m,k} * gamma_{m,k'} )
%
%  (These are the correct definitions matched to the official reference
%  code; note the ratio beta_{m,k}/beta_{m,k'} is NOT square-rooted.)
%
%  INPUTS:
%     eta_mk     - M-by-K power control coefficients (NOT sqrt).
%                  Must satisfy sum_k eta_{m,k} gamma_{m,k} <= 1/N per m.
%     beta_mk    - M-by-K large-scale fading
%     gamma_mk   - M-by-K MMSE per-component variance
%     pilotseq   - tau_p-by-K pilot matrix
%     N, rho_d   - antennas per AP, normalised DL SNR
%     tau_c,tau_p- coherence and pilot lengths
%
%  OUTPUTS:
%     SE_k   - K-by-1 achievable SE of each user, bit/s/Hz
%     SINR_k - K-by-1 effective SINR

    [M, K] = size(beta_mk);
    sqrt_eta = sqrt(max(eta_mk, 0));                      % M x K
    pilot_abs = abs(pilotseq' * pilotseq);                % K x K, real

    SINR_k = zeros(K, 1);
    for k = 1:K
        % ----- Numerator: |gamma_{k,k}^T sqrt(eta_k)|^2 --------------
        num_vec = gamma_mk(:, k);                         % M x 1
        num     = rho_d * N^2 * (num_vec' * sqrt_eta(:, k))^2;

        % ----- Denominator: pilot contamination sum -----------------
        pilot_sum = 0;
        for kp = 1:K
            if kp == k, continue; end
            % [gamma_{k,kp}]_m = pilot * gamma_{m,kp} * beta_{m,k}/beta_{m,kp}
            gkkp = pilot_abs(kp, k) ...
                 * (gamma_mk(:, kp) .* beta_mk(:, k) ./ beta_mk(:, kp));
            pilot_sum = pilot_sum + (gkkp' * sqrt_eta(:, kp))^2;
        end

        % ----- Denominator: non-coherent interference --------------
        %   sum_{k'} || D_{k',k} sqrt(eta_{k'}) ||^2
        %   = sum_{k'} sum_m beta_{m,k} gamma_{m,k'} eta_{m,k'}
        noncoh_sum = sum( sum( bsxfun(@times, gamma_mk .* eta_mk, ...
                                              beta_mk(:, k)), 1) );

        den = rho_d * N^2 * pilot_sum + rho_d * N * noncoh_sum + 1;

        SINR_k(k) = num / den;
    end

    SE_k = (1 - tau_p/tau_c) * log2(1 + SINR_k);
end
