function [EE_MbitPerJ, P_total_W, sum_rate_Mbps, breakdown] = ...
        compute_EE(eta_mk, gamma_mk, SE_k, params, N, M)
%COMPUTE_EE  Evaluate the total energy efficiency of the cell-free
%  network given a power-allocation matrix and the achievable SE.
%  Implements Eqs. (17)-(21) of Ngo et al. IEEE TGCN 2018.
%
%  INPUTS:
%     eta_mk   - M-by-K power control matrix (NOT sqrt)
%     gamma_mk - M-by-K MMSE variance (from compute_gamma.m)
%     SE_k     - K-by-1 achievable SE (bit/s/Hz)
%     params   - struct from get_params
%     N        - antennas per AP
%     M        - number of APs
%
%  OUTPUTS:
%     EE_MbitPerJ   - total energy efficiency (Mbit/Joule)
%     P_total_W     - total power consumption (Watt)
%     sum_rate_Mbps - sum-rate delivered (Mbit/s)
%     breakdown     - struct with 'PA', 'circuit', 'fixedBH', 'trafficBH'

    B_MHz   = params.B_MHz;
    eta_PA  = params.eta_PA;
    P_tc    = params.P_tc_W;
    P_0     = params.P_0_W;
    P_bt_per_Mbps = params.P_bt_WperGbps * 1e-3;   % convert W/(Gbit/s) -> W/(Mbit/s)
    P_per_antenna = params.P_per_antenna_W;
    N0            = params.noise_power;
    rho_d         = P_per_antenna * N / N0;        % normalised max TX SNR (per AP: N*P_per_antenna)

    % ---- (1) PA power --------------------------------------------
    % P_PA,m = (1/eta_PA) * rho_d * N0 * N * sum_k eta_mk gamma_mk
    PA_per_AP = (1/eta_PA) * rho_d * N0 * N * sum(eta_mk .* gamma_mk, 2);  % M x 1
    PA_total  = sum(PA_per_AP);

    % ---- (2) Per-antenna circuit power ---------------------------
    circuit_total = M * N * P_tc;

    % ---- (3) Fixed backhaul --------------------------------------
    fixed_BH = M * P_0;

    % ---- (4) Traffic-dependent backhaul -------------------------
    sum_SE        = sum(SE_k);                       % bit/s/Hz
    sum_rate_Mbps = B_MHz * sum_SE;                  % Mbit/s
    traffic_BH    = M * P_bt_per_Mbps * sum_rate_Mbps;  % W

    P_total_W = PA_total + circuit_total + fixed_BH + traffic_BH;

    EE_MbitPerJ = sum_rate_Mbps / P_total_W;         % (Mbit/s)/W = Mbit/J

    breakdown.PA        = PA_total;
    breakdown.circuit   = circuit_total;
    breakdown.fixedBH   = fixed_BH;
    breakdown.trafficBH = traffic_BH;
end
