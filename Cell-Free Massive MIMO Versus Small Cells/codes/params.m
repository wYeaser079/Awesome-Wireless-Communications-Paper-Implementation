function p = params()
%PARAMS Default simulation parameters for Ngo et al. (IEEE TWC, 2017).
%   "Cell-Free Massive MIMO versus Small Cells"
%
%   Usage: p = params();
%
%   All parameters match Table I and Section VI of the paper.

    % --- Network geometry ---
    p.area_side = 1000;         % Area side length [m] (1 km x 1 km square)

    % --- Network dimensions ---
    p.M = 100;                  % Number of Access Points (APs)
    p.K = 40;                   % Number of single-antenna users

    % --- Coherence interval ---
    p.tau_c  = 200;             % Coherence interval length [samples]
    p.tau_cf = 20;              % Cell-Free pilot length [samples]

    % --- Transmit powers (Table I) ---
    p.rho_d_mW = 200;           % DL transmit power per AP [mW]
    p.rho_u_mW = 100;           % UL transmit power per user [mW]
    p.rho_p_mW = 100;           % Pilot transmit power [mW]

    % --- Noise power: N0 = kT * B * NF ---
    p.noise_figure_dB = 9;      % Noise figure [dB]
    p.bandwidth_MHz   = 20;     % System bandwidth [MHz]
    kT_dBm_per_Hz = -174;       % Thermal noise floor [dBm/Hz] at T=290K
    p.noise_dBm = kT_dBm_per_Hz + 10*log10(p.bandwidth_MHz * 1e6) ...
                  + p.noise_figure_dB;   % approx -92 dBm

    % --- Normalized SNRs (linear) ---
    noise_mW   = 10^(p.noise_dBm / 10);
    p.rho_d_cf = p.rho_d_mW / noise_mW;  % CF downlink normalized SNR
    p.rho_u_cf = p.rho_u_mW / noise_mW;  % CF uplink normalized SNR
    p.rho_p_cf = p.rho_p_mW / noise_mW;  % CF pilot normalized SNR

    % --- Small-cell normalized SNRs (fair comparison, Section VI-B) ---
    %   rho_d^sc = (M/K) * rho_d^cf  (equal total radiated DL power)
    %   rho_u^sc = rho_u^cf           (same UL power per user)
    %   rho_p^sc = rho_p^cf           (same pilot power)
    p.rho_d_sc = (p.M / p.K) * p.rho_d_cf;
    p.rho_u_sc = p.rho_u_cf;
    p.rho_p_sc = p.rho_p_cf;

    % --- Small-cell pilot lengths ---
    %   tau_d^sc = tau_u^sc = tau_cf (paper uses same pilot length)
    p.tau_d_sc = p.tau_cf;       % SC downlink pilot length
    p.tau_u_sc = p.tau_cf;       % SC uplink pilot length

    % --- Path loss model (Eq. 52-53) ---
    p.d0 = 10;                  % Inner breakpoint [m]
    p.d1 = 50;                  % Outer breakpoint [m]
    p.hAP = 15;                 % AP antenna height [m]
    p.hUE = 1.65;               % User antenna height [m]
    p.freq_GHz = 1.9;           % Carrier frequency [GHz]

    % --- Shadow fading ---
    p.sigma_sf = 8;             % Shadow fading std dev [dB]
    p.d_decorr = 100;           % Decorrelation distance [m] (correlated case)
    p.delta    = 0.5;           % Shadowing balance factor (Eq. 54)

    % --- Pre-log factors for net throughput (Eqs. 56-57) ---
    %   Cell-Free: (1 - tau_cf/tau_c) / 2
    %   Small-Cell: (1 - (tau_d_sc + tau_u_sc)/tau_c) / 2
    p.prelog_cf = (1 - p.tau_cf / p.tau_c) / 2;
    p.prelog_sc = (1 - (p.tau_d_sc + p.tau_u_sc) / p.tau_c) / 2;

    % --- Bandwidth for throughput [Hz] ---
    p.B = p.bandwidth_MHz * 1e6;

    % --- Monte Carlo ---
    p.num_setups = 200;         % Random AP/user realizations (paper uses 200)

    % --- Greedy pilot assignment ---
    p.greedy_iterations = 10;   % Number of greedy iterations (Algorithm 1)

    % --- Power control ---
    p.bisection_tol = 0.01;     % Bisection tolerance
    p.bisection_maxiter = 30;   % Max bisection iterations

end
