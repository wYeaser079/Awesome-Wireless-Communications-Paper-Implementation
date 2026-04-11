%% main_Fig3_AP_selection_Pbt.m
% Reproduces Fig. 3 of Ngo et al. IEEE TGCN 2018:
%   Average EE versus the traffic-dependent backhaul coefficient P_bt,
%   for three schemes:
%     (a) no AP selection (Algorithm 1 baseline)
%     (b) largest-LSF-based  selection (Algorithm 3)
%     (c) received-power-based selection (Algorithm 2)
%
% Requires CVX.

clc; clear; close all;

params = get_params();
params.D     = 1;
params.tau_p = 40;            % enough pilots for K = 40 (no contamination)
N            = 1;
M            = 40;            % smaller than paper's 100 for tractability
K            = 20;
NMC          = 3;

P_bt_list_WperGbps = [0.25, 0.5, 1, 2];
rho_pct = params.rho_selection_pct;

EE_noSel   = zeros(size(P_bt_list_WperGbps));
EE_lsfSel  = zeros(size(P_bt_list_WperGbps));
EE_rpSel   = zeros(size(P_bt_list_WperGbps));

rng(11, 'twister');
rho_p    = params.pilot_power_W / params.noise_power;
rho_d    = params.P_per_antenna_W * N / params.noise_power;
epsT     = params.tau_c / (params.tau_c - params.tau_p);
rate_QoS = epsT * params.SE_target_bps_per_Hz * ones(K, 1);

for iP = 1:length(P_bt_list_WperGbps)
    params.P_bt_WperGbps = P_bt_list_WperGbps(iP);
    fprintf('\n=== Fig.3: P_bt = %.2f W/(Gbit/s) ===\n', params.P_bt_WperGbps);

    accNo = 0; accLSF = 0; accRP = 0; nOK = 0;
    for mc = 1:NMC
        beta_mk  = generate_large_scale_fading(M, K, params);
        pilotseq = assign_pilots(K, params.tau_p);
        gamma_mk = compute_gamma(beta_mk, pilotseq, params.tau_p, rho_p);

        % --- No selection: Algorithm 1 -----------------------------
        try
            [eta_noSel, ~, infoN] = algorithm1_EE_SCA(beta_mk, gamma_mk, ...
                                                      pilotseq, N, M, K, ...
                                                      params, rate_QoS);
        catch, infoN.feasible = false;
        end
        if ~infoN.feasible, continue; end

        SE_no = compute_SE_closedform(eta_noSel, beta_mk, gamma_mk, ...
                                       pilotseq, N, rho_d, ...
                                       params.tau_c, params.tau_p);
        EEno  = compute_EE(eta_noSel, gamma_mk, SE_no, params, N, M);

        % --- Algorithm 3: LSF-based selection ----------------------
        gamma_hatLSF = algorithm3_LSF_selection(beta_mk, gamma_mk, rho_pct);
        try
            [eta_lsf, ~, infoL] = algorithm1_EE_SCA(beta_mk, gamma_hatLSF, ...
                                                    pilotseq, N, M, K, ...
                                                    params, rate_QoS);
        catch, infoL.feasible = false;
        end
        if ~infoL.feasible
            eta_lsf = eta_noSel;
        end
        SE_lsf = compute_SE_closedform(eta_lsf, beta_mk, gamma_hatLSF, ...
                                        pilotseq, N, rho_d, ...
                                        params.tau_c, params.tau_p);
        EElsf  = compute_EE(eta_lsf, gamma_hatLSF, SE_lsf, params, N, M);

        % --- Algorithm 2: Received-power-based selection -----------
        gamma_hatRP = algorithm2_RP_selection(eta_noSel, beta_mk, ...
                                              gamma_mk, rho_pct);
        try
            [eta_rp, ~, infoR] = algorithm1_EE_SCA(beta_mk, gamma_hatRP, ...
                                                   pilotseq, N, M, K, ...
                                                   params, rate_QoS);
        catch, infoR.feasible = false;
        end
        if ~infoR.feasible
            eta_rp = eta_noSel;
        end
        SE_rp = compute_SE_closedform(eta_rp, beta_mk, gamma_hatRP, ...
                                       pilotseq, N, rho_d, ...
                                       params.tau_c, params.tau_p);
        EErp  = compute_EE(eta_rp, gamma_hatRP, SE_rp, params, N, M);

        accNo = accNo + EEno;
        accLSF = accLSF + EElsf;
        accRP  = accRP  + EErp;
        nOK = nOK + 1;
    end
    if nOK > 0
        EE_noSel(iP)  = accNo  / nOK;
        EE_lsfSel(iP) = accLSF / nOK;
        EE_rpSel(iP)  = accRP  / nOK;
    end
end

figure('Color','w','Name','Fig.3 AP Selection vs P_bt'); hold on;
plot(P_bt_list_WperGbps, EE_noSel,  '-o', 'LineWidth', 2, 'DisplayName', 'no AP selection');
plot(P_bt_list_WperGbps, EE_lsfSel, '-s', 'LineWidth', 2, 'DisplayName', 'largest-LSF based (Alg. 3)');
plot(P_bt_list_WperGbps, EE_rpSel,  '-d', 'LineWidth', 2, 'DisplayName', 'received-power based (Alg. 2)');
xlabel('P_{bt} [W/(Gbit/s)]'); ylabel('Average EE [Mbit/J]');
title(sprintf('Fig. 3: Impact of AP selection (M=%d, K=%d, N=%d, \\rho=%d%%)', ...
              M, K, N, rho_pct));
grid on; legend('Location','best');
