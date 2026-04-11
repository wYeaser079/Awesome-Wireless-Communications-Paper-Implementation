%% main_Fig2_power_control.m
% Reproduces Fig. 2 of Ngo et al. IEEE TGCN 2018:
%   "Average EE vs number of APs M, with and without EE power control".
%
% Two curves per K value: Algorithm 1 (optimised) vs. equal-power
% baseline (no power control).  Averaged over NMC Monte Carlo
% large-scale fading realisations.
%
% Requires CVX.

clc; clear; close all;

%% Parameters
params = get_params();
params.D     = 1;
params.tau_p = 20;
N            = 1;
M_LIST       = [20, 40, 60, 80, 100];          % Fig. 2 uses 20..100
K_LIST       = [20, 40];
NMC          = 5;                              % increase for smoother curves
rng(7, 'twister');

EE_opt = zeros(length(M_LIST), length(K_LIST));
EE_eq  = zeros(length(M_LIST), length(K_LIST));

for iK = 1:length(K_LIST)
    K = K_LIST(iK);
    for iM = 1:length(M_LIST)
        M = M_LIST(iM);
        fprintf('\n=== Fig.2: M=%d, K=%d ===\n', M, K);
        ee_opt_acc = 0; ee_eq_acc = 0; nfeas = 0;

        for mc = 1:NMC
            beta_mk  = generate_large_scale_fading(M, K, params);
            pilotseq = assign_pilots(K, params.tau_p);
            rho_p    = params.pilot_power_W / params.noise_power;
            gamma_mk = compute_gamma(beta_mk, pilotseq, params.tau_p, rho_p);
            rho_d    = params.P_per_antenna_W * N / params.noise_power;

            % --- No power control baseline (equal allocation) -------
            eta_eq = equal_power_allocation(gamma_mk, N);
            SE_eq  = compute_SE_closedform(eta_eq, beta_mk, gamma_mk, ...
                                           pilotseq, N, rho_d, ...
                                           params.tau_c, params.tau_p);
            EE_eqiter = compute_EE(eta_eq, gamma_mk, SE_eq, params, N, M);

            % --- Optimised EE power control -----------------------
            % Use the NO-PC SE as the QoS floor to keep comparison fair
            epsT = params.tau_c / (params.tau_c - params.tau_p);
            rate_QoS = epsT * min(SE_eq);       % equal floor across users
            rate_QoS = 0.5 * rate_QoS * ones(K, 1);   % relax slightly
            try
                [eta_opt, ~, info] = algorithm1_EE_SCA(beta_mk, gamma_mk, ...
                                                       pilotseq, N, M, K, ...
                                                       params, rate_QoS);
            catch
                info.feasible = false;
            end
            if info.feasible
                SE_opt = compute_SE_closedform(eta_opt, beta_mk, gamma_mk, ...
                                               pilotseq, N, rho_d, ...
                                               params.tau_c, params.tau_p);
                EE_optiter = compute_EE(eta_opt, gamma_mk, SE_opt, params, N, M);
                ee_opt_acc = ee_opt_acc + EE_optiter;
                ee_eq_acc  = ee_eq_acc  + EE_eqiter;
                nfeas = nfeas + 1;
            end
        end
        if nfeas > 0
            EE_opt(iM, iK) = ee_opt_acc / nfeas;
            EE_eq (iM, iK) = ee_eq_acc  / nfeas;
        end
    end
end

figure('Color','w', 'Name', 'Fig.2 Impact of Power Control'); hold on;
plot(M_LIST, EE_opt(:,1), '-o', 'LineWidth', 2, 'DisplayName', 'Alg.1, K=20');
plot(M_LIST, EE_eq(:,1),  '--o','LineWidth', 2, 'DisplayName', 'no-PC, K=20');
plot(M_LIST, EE_opt(:,2), '-s', 'LineWidth', 2, 'DisplayName', 'Alg.1, K=40');
plot(M_LIST, EE_eq(:,2),  '--s','LineWidth', 2, 'DisplayName', 'no-PC, K=40');
xlabel('Number of APs M'); ylabel('Average EE [Mbit/J]');
title('Fig. 2: Impact of power control (N = 1, D = 1 km)');
grid on; legend('Location','best');
