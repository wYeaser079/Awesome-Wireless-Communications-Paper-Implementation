function params = get_params()
%GET_PARAMS  Default system parameters for the cell-free EE study.
%
%  Follows Table I and Section VI of:
%    H. Q. Ngo, L.-N. Tran, T. Q. Duong, M. Matthaiou, E. G. Larsson,
%    "On the Total Energy Efficiency of Cell-Free Massive MIMO,"
%    IEEE Trans. Green Commun. and Networking, vol. 2, no. 1,
%    pp. 25-39, Mar. 2018.  arXiv:1702.07601.
%
%  Returns a struct PARAMS packaging every constant that the rest of
%  the simulation needs.  Change the fields here to re-run experiments
%  with different settings.

% ----- Geometry -----
params.D = 1;             % square side length in kilometres (1 or 2 km)

% ----- Bandwidth and noise -----
params.B_MHz = 20;                            % system bandwidth [MHz]
params.NF_dB = 9;                             % receiver noise figure [dB]
params.noise_power = 10^((-203.975 + 10*log10(params.B_MHz*1e6) ...
                          + params.NF_dB) / 10);   % in Watts

% ----- Path loss (three-slope, Tang-Sun-Gong 2001) -----
params.Hb = 15;           % AP height in metres
params.Hm = 1.65;         % user height in metres
params.fc_MHz = 1900;     % carrier frequency in MHz

% Hata-COST231 constants that reproduce the 140.7 dB the paper quotes
aL = (1.1*log10(params.fc_MHz) - 0.7) * params.Hm ...
     - (1.56*log10(params.fc_MHz) - 0.8);
params.L_dB = 46.3 + 33.9*log10(params.fc_MHz) ...
              - 13.82*log10(params.Hb) - aL;

params.d0_km = 0.01;      % 10 m  - inner frozen region
params.d1_km = 0.05;      % 50 m  - near/far breakpoint
params.sigma_shd_dB = 8;  % log-normal shadowing std-dev

% ----- Transmit power / pilots -----
params.P_per_antenna_W = 1;                % physical transmit power per AP antenna [W]
params.pilot_power_W   = 0.2;              % physical UL pilot power per symbol [W]
params.tau_c = 200;                        % coherence block length [samples]
params.tau_p = 20;                         % default pilot length [samples]

% ----- Hardware power consumption (Table I) -----
params.eta_PA   = 0.4;                     % PA efficiency
params.P_tc_W   = 0.2;                     % per-antenna circuit power [W]
params.P_0_W    = 0.825;                   % fixed per-AP backhaul power [W]
params.P_bt_WperGbps = 0.25;               % traffic-dep. backhaul [W/(Gbit/s)]

% ----- QoS and selection -----
params.SE_target_bps_per_Hz = 1;           % per-user QoS S_o,k  [bit/s/Hz]
params.rho_selection_pct    = 95;          % AP-selection percentile

% ----- SCA tolerances (Algorithm 1) -----
params.sca_max_iter = 15;                  % max SCA iterations
params.sca_tol      = 1e-2;                % relative tolerance

end
