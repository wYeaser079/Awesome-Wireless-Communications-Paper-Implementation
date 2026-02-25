% main_fig4.m — Reproduce Figure 4: Sum spectral efficiency vs K
%
% Paper: Zhang et al. (JCN, 2019), Section IV, Figure 4
%
% Plots the sum SE for ZF and CB receivers with equal power control
% as the number of users K varies from 5 to 30, for N = 2, 4, 6.
% Pilot length tau = K (orthogonal pilots, no pilot contamination).
%
% Fixed: M = 60.

clear; close all; clc;

% --- Configuration ---
p = params();

K_vec      = 5:5:30;            % Number of users to sweep
N_vec      = [2, 4, 6];         % Antennas per AP
num_setups = p.num_setups;

% Pre-allocate: |K_vec| x |N_vec|
sum_rate_ZF = zeros(length(K_vec), length(N_vec));
sum_rate_CB = zeros(length(K_vec), length(N_vec));

% --- Main simulation loop ---
for ik = 1:length(K_vec)
    p.K      = K_vec(ik);
    p.tau    = p.K;                  % Orthogonal pilots (tau = K)
    p.prelog = 1 - p.tau / p.T;     % Update pre-log factor

    for s = 1:num_setups
        [beta, ~, ~] = generate_setup(p);
        pilot_index  = assign_pilots(p.K, p.tau);
        [alpha, ~]   = estimate_channel(beta, pilot_index, p);

        for in = 1:length(N_vec)
            p.N = N_vec(in);

            R_ZF = compute_SE_ZF(beta, alpha, pilot_index, p);
            R_CB = compute_SE_CB(beta, alpha, pilot_index, p);

            sum_rate_ZF(ik, in) = sum_rate_ZF(ik, in) + sum(R_ZF);
            sum_rate_CB(ik, in) = sum_rate_CB(ik, in) + sum(R_CB);
        end
    end

    fprintf('K = %d done (%d/%d)\n', K_vec(ik), ik, length(K_vec));
end

% Average over setups
sum_rate_ZF = sum_rate_ZF / num_setups;
sum_rate_CB = sum_rate_CB / num_setups;

% --- Plot ---
figure; hold on; grid on; box on;

markers_ZF = {'-o', '-s', '-d'};
markers_CB = {'--o', '--s', '--d'};
colors     = {[0 0.45 0.74], [0.85 0.33 0.10], [0.47 0.67 0.19]};

for in = 1:length(N_vec)
    plot(K_vec, sum_rate_ZF(:, in), markers_ZF{in}, ...
         'Color', colors{in}, 'LineWidth', 1.5, 'MarkerSize', 7, ...
         'MarkerFaceColor', colors{in});
    plot(K_vec, sum_rate_CB(:, in), markers_CB{in}, ...
         'Color', colors{in}, 'LineWidth', 1.5, 'MarkerSize', 7);
end

xlabel('Number of users (K)');
ylabel('Sum spectral efficiency [bits/s/Hz]');
title('Figure 4: Sum SE vs K (Equal Power Control, \tau = K)');

legend_entries = cell(1, 2*length(N_vec));
for in = 1:length(N_vec)
    legend_entries{2*in-1} = sprintf('ZF, N = %d', N_vec(in));
    legend_entries{2*in}   = sprintf('CB, N = %d', N_vec(in));
end
legend(legend_entries, 'Location', 'northwest');

% --- Summary ---
fprintf('\n--- Results at K=10, N=6 ---\n');
idx_k10 = find(K_vec == 10);
idx_n6  = find(N_vec == 6);
if ~isempty(idx_k10) && ~isempty(idx_n6)
    fprintf('ZF sum rate: %.1f,  CB sum rate: %.1f\n', ...
        sum_rate_ZF(idx_k10, idx_n6), sum_rate_CB(idx_k10, idx_n6));
end
