function [beta, ap_pos, ue_pos] = generate_setup(p)
%GENERATE_SETUP Generate network topology and large-scale fading coefficients.
%   [beta, ap_pos, ue_pos] = generate_setup(p)
%
%   Places M APs on a grid, K users randomly, computes beta_mk using the
%   three-slope path loss model from Ngo et al. (IEEE TWC, 2017) with
%   wrap-around and log-normal shadow fading.
%
%   Inputs:
%       p       - Parameter struct from params()
%
%   Outputs:
%       beta    - M x K matrix of large-scale fading coefficients (linear)
%       ap_pos  - M x 2 matrix of AP positions [x, y] in meters
%       ue_pos  - K x 2 matrix of user positions [x, y] in meters

    D = p.area_side;  % Area side length [m]

    % --- Place M APs on a uniform grid ---
    M_per_side = ceil(sqrt(p.M));
    spacing = D / M_per_side;
    [grid_x, grid_y] = meshgrid(spacing/2 : spacing : D - spacing/2, ...
                                 spacing/2 : spacing : D - spacing/2);
    grid_x = grid_x(:);
    grid_y = grid_y(:);

    % Take exactly M APs (grid may have more than M if M is not a perfect square)
    ap_pos = [grid_x(1:p.M), grid_y(1:p.M)];

    % --- Place K users uniformly at random ---
    ue_pos = D * rand(p.K, 2);

    % --- Compute large-scale fading beta_mk ---
    beta = zeros(p.M, p.K);

    for m = 1:p.M
        for k = 1:p.K
            % Wrap-around distance (minimum over 9 tiled copies)
            d = wrap_around_distance(ap_pos(m,:), ue_pos(k,:), D);

            % Three-slope path loss [dB]
            PL_dB = three_slope_pathloss(d, p);

            % Add shadow fading [dB]: N(0, sigma_sf^2)
            shadow_dB = p.sigma_sf * randn();

            % Large-scale fading coefficient (linear scale)
            beta(m,k) = 10^((PL_dB + shadow_dB) / 10);
        end
    end
end


function d = wrap_around_distance(ap, ue, D)
%WRAP_AROUND_DISTANCE Minimum distance using wrap-around (9-tile method).
%   Treats the square area as a torus to avoid edge effects.

    % Offsets for the 9 copies: {-D, 0, +D} in both x and y
    offsets = [-D, 0, D];
    d_min_sq = inf;

    for dx = offsets
        for dy = offsets
            diff = ap - (ue + [dx, dy]);
            d_sq = diff(1)^2 + diff(2)^2;
            if d_sq < d_min_sq
                d_min_sq = d_sq;
            end
        end
    end

    d = sqrt(d_min_sq);
end


function PL_dB = three_slope_pathloss(d, p)
%THREE_SLOPE_PATHLOSS Compute path loss [dB] using Ngo et al. 2017 model.
%   Three-slope model (Eq. 46-47 in Ngo et al., IEEE TWC, 2017).
%
%   The model accounts for:
%     - Free-space-like propagation close to the AP
%     - Increasing path loss exponent at larger distances
%     - Antenna heights and carrier frequency

    f = p.freq_GHz * 1000;  % Frequency in MHz

    % Constant term L (from Hata-COST231 model)
    L = 46.3 + 33.9*log10(f) - 13.82*log10(p.hAP) ...
        - (1.1*log10(f) - 0.7)*p.hUE ...
        + (1.56*log10(f) - 0.8);

    % Enforce minimum distance to avoid singularity
    d = max(d, 1);  % At least 1 m

    % Convert distances to km (Hata-COST231 model uses km in log terms)
    d_km  = d / 1000;
    d0_km = p.d0 / 1000;
    d1_km = p.d1 / 1000;

    % Three-slope path loss [dB] (negative value = attenuation)
    if d <= p.d0
        % Region 1: Very close (flat + free-space-like)
        PL_dB = -L - 15*log10(d1_km) - 20*log10(d0_km);
    elseif d <= p.d1
        % Region 2: Near field (slope = -20 dB/decade)
        PL_dB = -L - 15*log10(d1_km) - 20*log10(d_km);
    else
        % Region 3: Far field (slope = -35 dB/decade)
        PL_dB = -L - 35*log10(d_km);
    end
end
