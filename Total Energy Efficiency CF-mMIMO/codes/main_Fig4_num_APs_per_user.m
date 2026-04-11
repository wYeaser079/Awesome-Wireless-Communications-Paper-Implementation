%% main_Fig4_num_APs_per_user.m
% Reproduces Fig. 4 of Ngo et al. IEEE TGCN 2018:
%   "Average number of APs selected per user vs M" for both selection
%   schemes, at two deployment sizes D = 1 km and D = 2 km.
%
% Note: Algorithm 2 (received-power-based) requires running Algorithm 1
% first.  This is expensive; consider reducing M_LIST or NMC for a
% quick run.
%
% Requires CVX.

clc; clear; close all;

params = get_params();
params.tau_p = 40;
N            = 1;
K            = 20;
NMC          = 3;
M_LIST       = [40, 60, 80, 100];
D_LIST       = [1, 2];                 % km

rho_pct = params.rho_selection_pct;
rho_p   = params.pilot_power_W / params.noise_power;

avg_num_rp  = zeros(length(M_LIST), length(D_LIST));
avg_num_lsf = zeros(length(M_LIST), length(D_LIST));

for iD = 1:length(D_LIST)
    params.D = D_LIST(iD);

    for iM = 1:length(M_LIST)
        M = M_LIST(iM);
        fprintf('\n=== Fig.4: D=%d km, M=%d ===\n', params.D, M);

        acc_rp = 0; acc_lsf = 0; nOK = 0;
        rng(100 + M + 17*iD, 'twister');

        for mc = 1:NMC
            beta_mk  = generate_large_scale_fading(M, K, params);
            pilotseq = assign_pilots(K, params.tau_p);
            gamma_mk = compute_gamma(beta_mk, pilotseq, params.tau_p, rho_p);

            % LSF-based: no need for Algorithm 1
            [~, ~, n_lsf] = algorithm3_LSF_selection(beta_mk, gamma_mk, rho_pct);

            % RP-based: need Algorithm 1 first
            epsT = params.tau_c / (params.tau_c - params.tau_p);
            rate_QoS = epsT * params.SE_target_bps_per_Hz * ones(K, 1);
            try
                [eta_star, ~, info] = algorithm1_EE_SCA(beta_mk, gamma_mk, ...
                                                         pilotseq, N, M, K, ...
                                                         params, rate_QoS);
            catch, info.feasible = false;
            end
            if ~info.feasible, continue; end
            [~, ~, n_rp] = algorithm2_RP_selection(eta_star, beta_mk, ...
                                                    gamma_mk, rho_pct);

            acc_lsf = acc_lsf + mean(n_lsf);
            acc_rp  = acc_rp  + mean(n_rp);
            nOK = nOK + 1;
        end
        if nOK > 0
            avg_num_lsf(iM, iD) = acc_lsf / nOK;
            avg_num_rp (iM, iD) = acc_rp  / nOK;
        end
    end
end

figure('Color','w','Name','Fig.4 # APs selected per user'); hold on;
plot(M_LIST, avg_num_rp(:,1),  '-o', 'LineWidth', 2, ...
     'DisplayName', 'Received-power (D=1 km)');
plot(M_LIST, avg_num_lsf(:,1), '--o','LineWidth', 2, ...
     'DisplayName', 'LSF-based      (D=1 km)');
plot(M_LIST, avg_num_rp(:,2),  '-s', 'LineWidth', 2, ...
     'DisplayName', 'Received-power (D=2 km)');
plot(M_LIST, avg_num_lsf(:,2), '--s','LineWidth', 2, ...
     'DisplayName', 'LSF-based      (D=2 km)');
xlabel('Number of APs M'); ylabel('Average number of APs per user');
title(sprintf('Fig. 4: APs per user (\\rho=%d%%, K=%d, N=%d)', ...
               rho_pct, K, N));
grid on; legend('Location','best');
