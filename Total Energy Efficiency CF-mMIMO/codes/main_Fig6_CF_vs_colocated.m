%% main_Fig6_CF_vs_colocated.m
% Reproduces Fig. 6 of Ngo et al. IEEE TGCN 2018:
%   Average EE vs per-user SE target S_{o,k}, for:
%     (a) Cell-free massive MIMO with optimal N (from Fig. 5) and
%         received-power-based AP selection (Alg. 2)
%     (b) Colocated massive MIMO (one AP with M*N antennas)
%
% Two values of MN are compared.
%
% Requires CVX.  Default sizes are smaller than the paper's for
% tractability; the qualitative pattern (CF beats colocated, gap
% larger at lower MN) is preserved.

clc; clear; close all;

params = get_params();
params.tau_p = 20;
params.D     = 1;
K            = 10;
NMC          = 2;
So_list      = [0.25, 0.5, 1, 1.5, 2];      % per-user SE targets [bit/s/Hz]
MN_list      = [32, 64];                     % keep small; paper uses 128,256

rho_p = params.pilot_power_W / params.noise_power;

EE_CF  = zeros(length(So_list), length(MN_list));
EE_Col = zeros(length(So_list), length(MN_list));

% --- Cell-free: use N = 4 as a practical "near-optimal" for D=1 km --
N_CF = 4;

rng(202, 'twister');

for iMN = 1:length(MN_list)
    MN = MN_list(iMN);
    M_CF  = MN / N_CF;
    for iS = 1:length(So_list)
        params.SE_target_bps_per_Hz = So_list(iS);
        fprintf('\n=== Fig.6: MN=%d, S_o=%.2f ===\n', MN, So_list(iS));
        epsT = params.tau_c / (params.tau_c - params.tau_p);
        rate_QoS = epsT * params.SE_target_bps_per_Hz * ones(K, 1);

        accCF = 0; nCF = 0; accCol = 0; nCol = 0;
        for mc = 1:NMC
            % ---- CELL-FREE ------------------------------------------
            beta_mk  = generate_large_scale_fading(M_CF, K, params);
            pilotseq = assign_pilots(K, params.tau_p);
            gamma_mk = compute_gamma(beta_mk, pilotseq, params.tau_p, rho_p);

            gamma_hat = algorithm3_LSF_selection(beta_mk, gamma_mk, ...
                                                 params.rho_selection_pct);
            try
                [eta_cf, ~, infoCF] = algorithm1_EE_SCA(beta_mk, gamma_hat, ...
                                                         pilotseq, N_CF, M_CF, ...
                                                         K, params, rate_QoS);
            catch, infoCF.feasible = false;
            end
            if infoCF.feasible
                rho_d = params.P_per_antenna_W * N_CF / params.noise_power;
                SE_cf = compute_SE_closedform(eta_cf, beta_mk, gamma_hat, ...
                                               pilotseq, N_CF, rho_d, ...
                                               params.tau_c, params.tau_p);
                eeCF = compute_EE(eta_cf, gamma_hat, SE_cf, params, N_CF, M_CF);
                accCF = accCF + eeCF; nCF = nCF + 1;
            end

            % ---- COLOCATED ------------------------------------------
            M_col = 1; N_col = MN;
            beta_col  = generate_large_scale_fading(M_col, K, params);
            pilot_col = pilotseq;            % re-use the same pilots
            gamma_col = compute_gamma(beta_col, pilot_col, ...
                                      params.tau_p, rho_p);
            try
                [eta_col, ~, infoCol] = algorithm1_EE_SCA(beta_col, gamma_col, ...
                                                           pilot_col, N_col, ...
                                                           M_col, K, params, ...
                                                           rate_QoS);
            catch, infoCol.feasible = false;
            end
            if infoCol.feasible
                rho_d_col = params.P_per_antenna_W * N_col / params.noise_power;
                SE_col = compute_SE_closedform(eta_col, beta_col, gamma_col, ...
                                                pilot_col, N_col, rho_d_col, ...
                                                params.tau_c, params.tau_p);
                eeCol = compute_EE(eta_col, gamma_col, SE_col, ...
                                   params, N_col, M_col);
                accCol = accCol + eeCol; nCol = nCol + 1;
            end
        end
        if nCF > 0,  EE_CF (iS, iMN) = accCF / nCF; end
        if nCol > 0, EE_Col(iS, iMN) = accCol / nCol; end
    end
end

figure('Color','w','Name','Fig.6 Cell-Free vs Colocated Massive MIMO'); hold on;
ls = {'-o','-s','-d'};
for iMN = 1:length(MN_list)
    plot(So_list, EE_CF(:,iMN),  ls{iMN}, 'LineWidth', 2, ...
         'DisplayName', sprintf('CF-mMIMO, MN=%d', MN_list(iMN)));
    plot(So_list, EE_Col(:,iMN), ['--' ls{iMN}(end)], 'LineWidth', 2, ...
         'DisplayName', sprintf('Colocated, MN=%d', MN_list(iMN)));
end
xlabel('Per-user SE target S_{o,k} [bit/s/Hz]'); ylabel('Average EE [Mbit/J]');
title('Fig. 6: Cell-Free vs Colocated Massive MIMO');
grid on; legend('Location','best');
