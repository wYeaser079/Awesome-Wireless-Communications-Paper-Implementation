% main_fig2.m — Reproduce Figure 2: Sum spectral efficiency vs M
%
% Paper: Zhang et al. (JCN, 2019), Section IV, Figure 2
%
% Plots the sum SE for ZF and CB receivers with equal power control
% (eta_k = 1) as M varies from 20 to 80, for N = 2, 4, 6.
%
% The closed-form expressions are used (Theorem 1 for ZF, Lemma 1 for CB).
% Results are averaged over num_setups random network realizations.

clear; close all; clc;

% --- Configuration ---
p = params();

M_vec      = 20:10:80;       % Number of APs to sweep
N_vec      = [2, 4, 6];      % Antennas per AP
num_setups = p.num_setups;    % Network realizations (default 100)

% Pre-allocate: |M_vec| x |N_vec|
sum_rate_ZF = zeros(length(M_vec), length(N_vec));
sum_rate_CB = zeros(length(M_vec), length(N_vec));

% --- Main simulation loop ---
for im = 1:length(M_vec)
    p.M = M_vec(im);

    for s = 1:num_setups
        % Generate one random network setup (beta depends on M, not N)
        [beta, ~, ~] = generate_setup(p);

        % Assign pilots (deterministic when tau >= K)
        pilot_index = assign_pilots(p.K, p.tau);

        % MMSE estimation parameters (depend on beta, not N)
        [alpha, ~] = estimate_channel(beta, pilot_index, p);

        % Compute SE for each N value (reuse same beta and alpha)
        for in = 1:length(N_vec)
            p.N = N_vec(in);

            R_ZF = compute_SE_ZF(beta, alpha, pilot_index, p);
            R_CB = compute_SE_CB(beta, alpha, pilot_index, p);

            sum_rate_ZF(im, in) = sum_rate_ZF(im, in) + sum(R_ZF);
            sum_rate_CB(im, in) = sum_rate_CB(im, in) + sum(R_CB);
        end
    end

    fprintf('M = %d done (%d/%d)\n', M_vec(im), im, length(M_vec));
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
    plot(M_vec, sum_rate_ZF(:, in), markers_ZF{in}, ...
         'Color', colors{in}, 'LineWidth', 1.5, 'MarkerSize', 7, ...
         'MarkerFaceColor', colors{in});
    plot(M_vec, sum_rate_CB(:, in), markers_CB{in}, ...
         'Color', colors{in}, 'LineWidth', 1.5, 'MarkerSize', 7);
end

xlabel('Number of APs (M)');
ylabel('Sum spectral efficiency [bits/s/Hz]');
title('Figure 2: Sum SE vs M (Equal Power Control)');

legend_entries = cell(1, 2*length(N_vec));
for in = 1:length(N_vec)
    legend_entries{2*in-1} = sprintf('ZF, N = %d', N_vec(in));
    legend_entries{2*in}   = sprintf('CB, N = %d', N_vec(in));
end
legend(legend_entries, 'Location', 'northwest');

fprintf('\n--- Results summary ---\n');
for in = 1:length(N_vec)
    fprintf('N = %d:  ZF sum rate at M=60 = %.1f,  CB = %.1f,  ratio = %.2f\n', ...
        N_vec(in), ...
        sum_rate_ZF(M_vec == 60, in), ...
        sum_rate_CB(M_vec == 60, in), ...
        sum_rate_ZF(M_vec == 60, in) / sum_rate_CB(M_vec == 60, in));
end
