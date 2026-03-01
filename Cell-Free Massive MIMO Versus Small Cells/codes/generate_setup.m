function [beta, ap_pos, ue_pos] = generate_setup(p, correlated_shadowing)
%GENERATE_SETUP Generate network topology and large-scale fading coefficients.
%   [beta, ap_pos, ue_pos] = generate_setup(p, correlated_shadowing)
%
%   Places M APs and K users randomly in a 1km x 1km area. Computes beta_mk
%   using the three-slope path loss model (Eq. 52-53) with wrap-around
%   topology and shadow fading (correlated or uncorrelated).
%
%   Shadow fading is applied ONLY when d_mk > d1 = 50m (Section VI-A1).
%
%   Inputs:
%       p                      - Parameter struct from params()
%       correlated_shadowing   - true: use two-component model (Eq. 54-55)
%                                false: use i.i.d. shadow fading
%
%   Outputs:
%       beta    - M x K matrix of large-scale fading coefficients (linear)
%       ap_pos  - M x 2 matrix of AP positions [x, y] in meters
%       ue_pos  - K x 2 matrix of user positions [x, y] in meters

    if nargin < 2
        correlated_shadowing = false;
    end

    D = p.area_side;

    % --- Place M APs uniformly at random ---
    ap_pos = D * rand(p.M, 2);

    % --- Place K users uniformly at random ---
    ue_pos = D * rand(p.K, 2);

    % --- Compute distances with wrap-around ---
    dist = zeros(p.M, p.K);
    for m = 1:p.M
        for k = 1:p.K
            dist(m,k) = wrap_around_distance(ap_pos(m,:), ue_pos(k,:), D);
        end
    end

    % --- Compute path loss (dB) ---
    PL_dB = zeros(p.M, p.K);
    for m = 1:p.M
        for k = 1:p.K
            PL_dB(m,k) = three_slope_pathloss(dist(m,k), p);
        end
    end

    % --- Generate shadow fading ---
    if correlated_shadowing
        shadow_dB = generate_correlated_shadow(p, ap_pos, ue_pos, dist);
    else
        shadow_dB = generate_uncorrelated_shadow(p, dist);
    end

    % --- Combine path loss and shadow fading ---
    %   beta_mk = 10^((PL_dB + shadow_dB - noise_dBm) / 10)
    %   We normalize by noise power so that beta includes the path to SNR.
    %   Actually, in the paper, beta is pure channel gain. The normalized SNR
    %   rho_d, rho_u already account for noise. So:
    %   beta_mk = 10^(PL_dB + shadow_dB) / 10   [in linear, relative to noise]
    beta = 10.^((PL_dB + shadow_dB - p.noise_dBm) / 10);

end


function shadow_dB = generate_correlated_shadow(p, ap_pos, ue_pos, dist)
%GENERATE_CORRELATED_SHADOW Two-component shadowing model (Eq. 54-55).
%   z_mk = sqrt(delta)*a_m + sqrt(1-delta)*b_k
%   where a_m and b_k have spatial correlation (Eq. 55).

    D = p.area_side;

    % --- AP-AP correlation matrix (M x M) ---
    dist_AP = zeros(p.M, p.M);
    for m1 = 1:p.M
        for m2 = m1+1:p.M
            dist_AP(m1,m2) = wrap_around_distance(ap_pos(m1,:), ap_pos(m2,:), D);
            dist_AP(m2,m1) = dist_AP(m1,m2);
        end
    end
    C_AP = 2.^(-dist_AP / p.d_decorr);  % Eq. 55

    % --- UE-UE correlation matrix (K x K) ---
    dist_UE = zeros(p.K, p.K);
    for k1 = 1:p.K
        for k2 = k1+1:p.K
            dist_UE(k1,k2) = wrap_around_distance(ue_pos(k1,:), ue_pos(k2,:), D);
            dist_UE(k2,k1) = dist_UE(k1,k2);
        end
    end
    C_UE = 2.^(-dist_UE / p.d_decorr);  % Eq. 55

    % --- Generate correlated Gaussian samples ---
    %   Use Cholesky: a = L_AP * randn(M,1), b = L_UE * randn(K,1)
    L_AP = chol(C_AP, 'lower');
    L_UE = chol(C_UE, 'lower');

    a = L_AP * randn(p.M, 1);   % AP-dependent component ~ N(0, C_AP)
    b = L_UE * randn(p.K, 1);   % UE-dependent component ~ N(0, C_UE)

    % --- Two-component model (Eq. 54): z_mk = sqrt(delta)*a_m + sqrt(1-delta)*b_k ---
    z = sqrt(p.delta) * repmat(a, 1, p.K) + sqrt(1 - p.delta) * repmat(b', p.M, 1);

    % --- Apply shadow fading only where d_mk > d1 ---
    shadow_dB = zeros(p.M, p.K);
    far_mask = dist > p.d1;
    shadow_dB(far_mask) = p.sigma_sf * z(far_mask);

end


function shadow_dB = generate_uncorrelated_shadow(p, dist)
%GENERATE_UNCORRELATED_SHADOW I.i.d. shadow fading (N(0, sigma_sf^2) per link).

    z = randn(p.M, p.K);
    shadow_dB = zeros(p.M, p.K);
    far_mask = dist > p.d1;
    shadow_dB(far_mask) = p.sigma_sf * z(far_mask);

end


function d = wrap_around_distance(pos1, pos2, D)
%WRAP_AROUND_DISTANCE Minimum distance using wrap-around (9-tile method).
%   Treats the square area as a torus to avoid edge effects.

    offsets = [-D, 0, D];
    d_min_sq = inf;

    for dx = offsets
        for dy = offsets
            diff = pos1 - (pos2 + [dx, dy]);
            d_sq = diff(1)^2 + diff(2)^2;
            if d_sq < d_min_sq
                d_min_sq = d_sq;
            end
        end
    end

    d = sqrt(d_min_sq);
end


function PL_dB = three_slope_pathloss(d, p)
%THREE_SLOPE_PATHLOSS Path loss [dB] using the three-slope model (Eq. 52-53).
%   Hata-COST231 propagation model with breakpoints at d0 and d1.

    f = p.freq_GHz * 1000;  % Frequency in MHz

    % Constant L (Eq. 53)
    L = 46.3 + 33.9*log10(f) - 13.82*log10(p.hAP) ...
        - (1.1*log10(f) - 0.7)*p.hUE ...
        + (1.56*log10(f) - 0.8);

    % Enforce minimum distance
    d = max(d, 1);

    % Convert to km (Hata model uses km)
    d_km  = d / 1000;
    d0_km = p.d0 / 1000;
    d1_km = p.d1 / 1000;

    % Three-slope path loss (Eq. 52)
    if d <= p.d0
        PL_dB = -L - 15*log10(d1_km) - 20*log10(d0_km);
    elseif d <= p.d1
        PL_dB = -L - 15*log10(d1_km) - 20*log10(d_km);
    else
        PL_dB = -L - 35*log10(d_km);
    end
end
