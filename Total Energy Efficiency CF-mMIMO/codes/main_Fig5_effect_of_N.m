%% main_Fig5_effect_of_N.m
% Reproduces Fig. 5 of Ngo et al. IEEE TGCN 2018:
%   Average EE as a function of N (antennas per AP) with the total
%   service antenna count MN fixed.  Three scenarios are shown:
%
%     (i)   D=1 km, P_bt=0.25, S_o=1  (dense + cheap backhaul)
%     (ii)  D=2 km, P_bt=0.50, S_o=1  (sparse + expensive backhaul)
%     (iii) D=2 km, P_bt=0.25, S_o=2  (harder QoS)
%
% Expect an inverted-U shape with a peak at some N*.  AP selection
% (Algorithm 3) is used throughout.
%
% Requires CVX.

clc; clear; close all;

MN   = 128;                      % total service antennas (was 256 in paper)
N_LIST = [1, 2, 4, 8, 16];       % antennas per AP
NMC  = 2;                        % increase for smoother curves
K    = 10;
rng(123, 'twister');

scenarios(1).D = 1; scenarios(1).P_bt = 0.25; scenarios(1).S_o = 1;
scenarios(2).D = 2; scenarios(2).P_bt = 0.50; scenarios(2).S_o = 1;
scenarios(3).D = 2; scenarios(3).P_bt = 0.25; scenarios(3).S_o = 2;

EE = zeros(length(N_LIST), length(scenarios));

for iS = 1:length(scenarios)
    base = get_params();
    base.tau_p = 40;
    base.D     = scenarios(iS).D;
    base.P_bt_WperGbps     = scenarios(iS).P_bt;
    base.SE_target_bps_per_Hz = scenarios(iS).S_o;
    rho_p = base.pilot_power_W / base.noise_power;

    for iN = 1:length(N_LIST)
        N = N_LIST(iN);
        M = MN / N;
        fprintf('\n=== Fig.5: scen %d, N=%d, M=%d ===\n', iS, N, M);

        acc = 0; nOK = 0;
        for mc = 1:NMC
            beta_mk  = generate_large_scale_fading(M, K, base);
            pilotseq = assign_pilots(K, base.tau_p);
            gamma_mk = compute_gamma(beta_mk, pilotseq, base.tau_p, rho_p);

            gamma_hat = algorithm3_LSF_selection(beta_mk, gamma_mk, ...
                                                 base.rho_selection_pct);

            epsT = base.tau_c / (base.tau_c - base.tau_p);
            rate_QoS = epsT * base.SE_target_bps_per_Hz * ones(K, 1);
            try
                [eta_mk, ~, info] = algorithm1_EE_SCA(beta_mk, gamma_hat, ...
                                                       pilotseq, N, M, K, ...
                                                       base, rate_QoS);
            catch, info.feasible = false;
            end
            if ~info.feasible, continue; end

            rho_d = base.P_per_antenna_W * N / base.noise_power;
            SE_k  = compute_SE_closedform(eta_mk, beta_mk, gamma_hat, ...
                                          pilotseq, N, rho_d, ...
                                          base.tau_c, base.tau_p);
            EEiter = compute_EE(eta_mk, gamma_hat, SE_k, base, N, M);
            acc = acc + EEiter; nOK = nOK + 1;
        end
        if nOK > 0, EE(iN, iS) = acc / nOK; end
    end
end

figure('Color','w','Name','Fig.5 Effect of N (antennas per AP)'); hold on;
plot(N_LIST, EE(:,1), '-o', 'LineWidth', 2, ...
     'DisplayName', 'D=1, P_{bt}=0.25, S_o=1');
plot(N_LIST, EE(:,2), '-s', 'LineWidth', 2, ...
     'DisplayName', 'D=2, P_{bt}=0.50, S_o=1');
plot(N_LIST, EE(:,3), '-d', 'LineWidth', 2, ...
     'DisplayName', 'D=2, P_{bt}=0.25, S_o=2');
xlabel('N (antennas per AP)'); ylabel('Average EE [Mbit/J]');
set(gca, 'XScale', 'log');
title(sprintf('Fig. 5: EE vs N (MN = %d, K = %d)', MN, K));
grid on; legend('Location','best');
