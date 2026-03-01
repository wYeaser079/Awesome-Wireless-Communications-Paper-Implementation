%MAIN_FIG11 Reproduce Figure 11 from Ngo et al. (IEEE TWC, 2017).
%   Average downlink net throughput versus the number of APs M,
%   for tau_cf = 5 and 20. Cell-Free vs Small-Cell.
%
%   Parameters: K=20, uncorrelated shadowing, M sweeps from 30 to 100.

clear; clc; close all;

%% Setup
p = params();
p.K = 20;  % Fixed K=20 for this figure

fprintf('=== Figure 11: Average DL Throughput vs M ===\n');

M_values = 30:10:100;
tau_values = [5, 20];

cf_avg = zeros(length(tau_values), length(M_values));
sc_avg = zeros(length(tau_values), length(M_values));

num_setups = 50;

%% Main simulation loop
for t_idx = 1:length(tau_values)
    p_local = p;
    p_local.tau_cf = tau_values(t_idx);
    p_local.tau_d_sc = tau_values(t_idx);
    p_local.tau_u_sc = tau_values(t_idx);
    p_local.prelog_cf = (1 - p_local.tau_cf / p_local.tau_c) / 2;
    p_local.prelog_sc = (1 - (p_local.tau_d_sc + p_local.tau_u_sc) / p_local.tau_c) / 2;

    for m_idx = 1:length(M_values)
        p_local.M = M_values(m_idx);

        % Update SC power for fair comparison
        p_local.rho_d_sc = (p_local.M / p_local.K) * p_local.rho_d_cf;

        fprintf('tau=%d, M=%d\n', tau_values(t_idx), M_values(m_idx));

        cf_sum = 0;
        sc_sum = 0;

        for s = 1:num_setups
            [beta, ~, ~] = generate_setup(p_local, false);

            % --- Cell-Free ---
            pilot_idx = assign_pilots(beta, [], [], p_local, 'greedy');
            gamma = compute_gamma(beta, pilot_idx, p_local);
            eta_dl = power_control_CF_DL(beta, gamma, pilot_idx, p_local);
            [~, tp_cf] = compute_SE_CF_DL(beta, gamma, pilot_idx, p_local, eta_dl);
            cf_sum = cf_sum + min(tp_cf);

            % --- Small-Cell ---
            ap_sel = select_AP(beta);
            pilot_idx_sc = randi(p_local.tau_cf, p_local.K, 1);
            [alpha_d, ~] = power_control_SC(beta, pilot_idx_sc, ap_sel, p_local);
            [~, tp_sc] = compute_SE_SC_DL(beta, pilot_idx_sc, ap_sel, p_local, alpha_d);
            sc_sum = sc_sum + min(tp_sc);
        end

        cf_avg(t_idx, m_idx) = cf_sum / num_setups;
        sc_avg(t_idx, m_idx) = sc_sum / num_setups;
    end
end

%% Plot Figure 11
to_mbps = 1e-6;

figure('Position', [100, 100, 600, 450]);
hold on; grid on; box on;

for t_idx = 1:length(tau_values)
    plot(M_values, cf_avg(t_idx,:) * to_mbps, 'b-o', ...
         'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
    plot(M_values, sc_avg(t_idx,:) * to_mbps, 'r--s', ...
         'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'r');
end

xlabel('Number of APs (M)');
ylabel('Average Downlink Net Throughput (Mbits/s)');
legend(sprintf('Cell-Free, \\tau^{cf}=%d', tau_values(1)), ...
       sprintf('Small-Cell, \\tau^{cf}=%d', tau_values(1)), ...
       sprintf('Cell-Free, \\tau^{cf}=%d', tau_values(2)), ...
       sprintf('Small-Cell, \\tau^{cf}=%d', tau_values(2)), ...
       'Location', 'northwest');
title('Figure 11: Average DL Throughput vs M');
