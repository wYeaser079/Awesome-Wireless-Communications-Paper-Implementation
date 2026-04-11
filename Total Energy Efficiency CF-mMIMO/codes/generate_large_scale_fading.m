function beta_mk = generate_large_scale_fading(M, K, params)
%GENERATE_LARGE_SCALE_FADING  Place M APs and K users uniformly in a
%  D-by-D km square with 3x3 wrap-around and compute the large-scale
%  fading coefficients beta_mk using the three-slope path-loss model
%  of Tang-Sun-Gong 2001 + log-normal shadowing.
%
%  INPUTS:
%     M       - number of APs
%     K       - number of users
%     params  - struct from GET_PARAMS
%
%  OUTPUT:
%     beta_mk - M-by-K matrix of large-scale fading coefficients (linear)
%
%  This is a direct port of GETSLOWFADING.M from
%  https://github.com/tranlenam/cellfreeMIMOenergyefficiency
%  (Ngo-Tran-Duong-Matthaiou-Larsson, IEEE TGCN 2018) with comments and
%  clearer vectorisation.

    D  = params.D;                 % km
    L  = params.L_dB;              % dB constant
    d0 = params.d0_km;             % km - inner cap
    d1 = params.d1_km;             % km - near/far break
    sigma_shd = params.sigma_shd_dB;

    % --- Random positions in the reference square ----------------------
    AP  = unifrnd(-D/2, D/2, M, 2);
    UE  = unifrnd(-D/2, D/2, K, 2);

    % --- 8 periodic copies of the AP positions (3x3 torus) --------------
    shifts = D * [ 1  0;  0  1; -1  0;  0 -1;  1 -1; -1  1;  1  1; -1 -1];

    % Distance matrix (M-by-K) using the shortest wrapped distance
    dist_km = zeros(M, K);
    for m = 1:M
        % M-by-(8+1) copies of this AP
        AP_copies = [AP(m,:); bsxfun(@plus, AP(m,:), shifts)];
        for k = 1:K
            diffs = bsxfun(@minus, AP_copies, UE(k,:));
            dist_km(m,k) = min( sqrt( sum(diffs.^2, 2) ) );
        end
    end

    % --- Three-slope path loss (dB) ------------------------------------
    beta_dB = zeros(M, K);
    for m = 1:M
        for k = 1:K
            d = dist_km(m,k);
            if d < d0
                beta_dB(m,k) = -L - 15*log10(d1) - 20*log10(d0);
            elseif d <= d1
                beta_dB(m,k) = -L - 15*log10(d1) - 20*log10(d);
            else
                beta_dB(m,k) = -L - 35*log10(d) ...
                               + sigma_shd * randn;    % shadow only in far regime
            end
        end
    end

    beta_mk = 10.^(beta_dB / 10);
end
