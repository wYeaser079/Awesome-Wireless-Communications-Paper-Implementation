% main_fig3.m — Reproduce Figure 3: CDF of per-user spectral efficiency
%
% Paper: Zhang et al. (JCN, 2019), Section IV, Figure 3
%
% Plots the empirical CDF of per-user SE for ZF and CB receivers
% with equal power control (eta_k = 1), comparing:
%   - tau = 10 (orthogonal pilots, no pilot contamination since K = 10)
%   - tau = 8  (pilot contamination, since tau < K)
%
% Fixed: M = 60, K = 10, N = 6.
%
% Key observation (Remark 2): ZF is highly sensitive to channel
% estimation quality — the 5%-outage rate drops sharply when tau < K.

clear; close all; clc;

% --- Configuration ---
p = params();

tau_vec    = [10, 8];           % Pilot lengths to compare
num_setups = p.num_setups;      % Network realizations (default 100)

% Collect per-user rates for CDF
rates_ZF = cell(1, length(tau_vec));
rates_CB = cell(1, length(tau_vec));

% --- Main simulation loop ---
for it = 1:length(tau_vec)
    p.tau    = tau_vec(it);
    p.prelog = 1 - p.tau / p.T;         % Update pre-log factor

    all_ZF = zeros(p.K, num_setups);
    all_CB = zeros(p.K, num_setups);

    for s = 1:num_setups
        % Generate one random network (large-scale fading)
        [beta, ~, ~] = generate_setup(p);

        % Pilot assignment (orthogonal if tau >= K, random otherwise)
        pilot_index = assign_pilots(p.K, p.tau);

        % MMSE estimation parameters (Eq. 5, 8)
        [alpha, ~] = estimate_channel(beta, pilot_index, p);

        % Closed-form SE: ZF (Theorem 1, Eq. 16) and CB (Lemma 1, Eq. 33)
        all_ZF(:, s) = compute_SE_ZF(beta, alpha, pilot_index, p);
        all_CB(:, s) = compute_SE_CB(beta, alpha, pilot_index, p);
    end

    rates_ZF{it} = all_ZF(:);           % K * num_setups x 1
    rates_CB{it} = all_CB(:);

    fprintf('tau = %d done\n', tau_vec(it));
end

% --- Plot CDF ---
figure; hold on; grid on; box on;

colors = {[0 0.45 0.74], [0.85 0.33 0.10]};   % blue for tau=10, red for tau=8

for it = 1:length(tau_vec)
    [f_zf, x_zf] = ecdf(rates_ZF{it});
    [f_cb, x_cb] = ecdf(rates_CB{it});

    plot(x_zf, f_zf, '-',  'Color', colors{it}, 'LineWidth', 1.5);
    plot(x_cb, f_cb, '--', 'Color', colors{it}, 'LineWidth', 1.5);
end

xlabel('Per-user spectral efficiency [bits/s/Hz]');
ylabel('CDF');
title('Figure 3: CDF of Per-User SE (Equal Power Control)');

legend_entries = cell(1, 2 * length(tau_vec));
for it = 1:length(tau_vec)
    legend_entries{2*it-1} = sprintf('ZF, \\tau = %d', tau_vec(it));
    legend_entries{2*it}   = sprintf('CB, \\tau = %d', tau_vec(it));
end
legend(legend_entries, 'Location', 'southeast');

% --- Summary statistics ---
fprintf('\n--- CDF Statistics ---\n');
for it = 1:length(tau_vec)
    zf_5  = prctile(rates_ZF{it}, 5);
    cb_5  = prctile(rates_CB{it}, 5);
    zf_50 = median(rates_ZF{it});
    cb_50 = median(rates_CB{it});

    fprintf('tau = %d:  ZF [5%%=%.2f, 50%%=%.2f]  CB [5%%=%.2f, 50%%=%.2f]\n', ...
        tau_vec(it), zf_5, zf_50, cb_5, cb_50);
end
