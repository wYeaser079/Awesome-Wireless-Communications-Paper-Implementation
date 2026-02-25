% main_fig1.m — Reproduce Figure 1: Validation of Theorem 1
%
% Paper: Zhang et al. (JCN, 2019), Section IV, Figure 1
%
% Compares the closed-form ZF SE (Theorem 1, Eq. 16) with the
% Monte Carlo ZF SE (exact SINR, Eq. 15) to validate the
% tight approximation of Theorem 1.
%
% Plots per-user SE as a grouped bar chart, averaged over setups.
% Fixed: M = 60, K = 10, N = 6.

clear; close all; clc;

% --- Configuration ---
p = params();

% Use fewer setups since Monte Carlo is expensive (300 channel realizations each)
num_setups = 20;

R_ZF_cf_avg = zeros(p.K, 1);   % Closed-form (Theorem 1, Eq. 16)
R_ZF_mc_avg = zeros(p.K, 1);   % Monte Carlo (exact SINR, Eq. 15)

% --- Main simulation loop ---
for s = 1:num_setups
    [beta, ~, ~] = generate_setup(p);
    pilot_index  = assign_pilots(p.K, p.tau);
    [alpha, c_coeff] = estimate_channel(beta, pilot_index, p);

    % Closed-form SE (Theorem 1)
    R_ZF_cf = compute_SE_ZF(beta, alpha, pilot_index, p);

    % Monte Carlo SE (Eq. 15 — exact SINR averaged over channel realizations)
    R_ZF_mc = compute_SE_ZF_montecarlo(beta, alpha, c_coeff, pilot_index, p);

    R_ZF_cf_avg = R_ZF_cf_avg + R_ZF_cf;
    R_ZF_mc_avg = R_ZF_mc_avg + R_ZF_mc;

    fprintf('Setup %d/%d done\n', s, num_setups);
end

R_ZF_cf_avg = R_ZF_cf_avg / num_setups;
R_ZF_mc_avg = R_ZF_mc_avg / num_setups;

% Sort by closed-form rate (descending) for cleaner visualization
[R_ZF_cf_sorted, idx] = sort(R_ZF_cf_avg, 'descend');
R_ZF_mc_sorted = R_ZF_mc_avg(idx);

% --- Plot ---
figure; hold on; grid on; box on;

bar_data = [R_ZF_cf_sorted, R_ZF_mc_sorted];
b = bar(bar_data, 'grouped');
b(1).FaceColor = [0 0.45 0.74];
b(2).FaceColor = [0.85 0.33 0.10];

xlabel('User index (sorted by rate)');
ylabel('Spectral efficiency [bits/s/Hz]');
title('Figure 1: Validation of Theorem 1 (ZF)');
legend('Closed-form (Eq. 16)', 'Monte Carlo (Eq. 15)', 'Location', 'northeast');

% --- Summary ---
fprintf('\n--- Validation of Theorem 1 ---\n');
fprintf('Sum SE — Closed-form: %.2f,  Monte Carlo: %.2f\n', ...
    sum(R_ZF_cf_avg), sum(R_ZF_mc_avg));
fprintf('Max per-user error: %.4f bits/s/Hz\n', ...
    max(abs(R_ZF_cf_avg - R_ZF_mc_avg)));
fprintf('Relative error: %.2f%%\n', ...
    100 * norm(R_ZF_cf_avg - R_ZF_mc_avg) / norm(R_ZF_cf_avg));
