% main_fig6.m — Reproduce Figure 6: CDF with power control
%
% Paper: Zhang et al. (JCN, 2019), Section IV, Figure 6
%
% Compares the CDF of per-user SE under:
%   - Equal Power Control (EPC): eta_k = 1 for all k
%   - SCA power control (Algorithms 1 & 3): optimized {eta_k}
% for both ZF and CB receivers.
%
% Requires: CVX toolbox (http://cvxr.com/cvx/)
% Fixed: M = 60, K = 10, N = 6.

clear; close all; clc;

% --- Check CVX ---
if ~exist('cvx_begin', 'file')
    error('This script requires CVX. Install from http://cvxr.com/cvx/');
end

% --- Configuration ---
p = params();
num_setups = p.num_setups;

% Pre-allocate per-user rates
rates_ZF_epc = zeros(p.K, num_setups);
rates_CB_epc = zeros(p.K, num_setups);
rates_ZF_sca = zeros(p.K, num_setups);
rates_CB_sca = zeros(p.K, num_setups);

% --- Main simulation loop ---
for s = 1:num_setups
    [beta, ~, ~] = generate_setup(p);
    pilot_index  = assign_pilots(p.K, p.tau);
    [alpha, ~]   = estimate_channel(beta, pilot_index, p);

    % Equal power control (eta_k = 1)
    rates_ZF_epc(:, s) = compute_SE_ZF(beta, alpha, pilot_index, p);
    rates_CB_epc(:, s) = compute_SE_CB(beta, alpha, pilot_index, p);

    % SCA-optimized power control
    eta_zf = power_control_SCA_ZF(beta, alpha, pilot_index, p);
    eta_cb = power_control_SCA_CB(beta, alpha, pilot_index, p);

    rates_ZF_sca(:, s) = compute_SE_ZF(beta, alpha, pilot_index, p, eta_zf);
    rates_CB_sca(:, s) = compute_SE_CB(beta, alpha, pilot_index, p, eta_cb);

    fprintf('Setup %d/%d done\n', s, num_setups);
end

% --- Plot CDF ---
figure; hold on; grid on; box on;

[f, x] = ecdf(rates_ZF_epc(:));
plot(x, f, '-',  'Color', [0 0.45 0.74], 'LineWidth', 1.5);
[f, x] = ecdf(rates_ZF_sca(:));
plot(x, f, '--', 'Color', [0 0.45 0.74], 'LineWidth', 1.5);
[f, x] = ecdf(rates_CB_epc(:));
plot(x, f, '-',  'Color', [0.85 0.33 0.10], 'LineWidth', 1.5);
[f, x] = ecdf(rates_CB_sca(:));
plot(x, f, '--', 'Color', [0.85 0.33 0.10], 'LineWidth', 1.5);

xlabel('Per-user spectral efficiency [bits/s/Hz]');
ylabel('CDF');
title('Figure 6: CDF of Per-User SE (EPC vs SCA)');
legend('ZF, EPC', 'ZF, SCA', 'CB, EPC', 'CB, SCA', 'Location', 'southeast');

% --- Summary ---
fprintf('\n--- CDF Statistics ---\n');
fprintf('ZF EPC: [5%%=%.2f, 50%%=%.2f]  ZF SCA: [5%%=%.2f, 50%%=%.2f]\n', ...
    prctile(rates_ZF_epc(:), 5), median(rates_ZF_epc(:)), ...
    prctile(rates_ZF_sca(:), 5), median(rates_ZF_sca(:)));
fprintf('CB EPC: [5%%=%.2f, 50%%=%.2f]  CB SCA: [5%%=%.2f, 50%%=%.2f]\n', ...
    prctile(rates_CB_epc(:), 5), median(rates_CB_epc(:)), ...
    prctile(rates_CB_sca(:), 5), median(rates_CB_sca(:)));
fprintf('\nSum SE — ZF EPC: %.1f, ZF SCA: %.1f, CB EPC: %.1f, CB SCA: %.1f\n', ...
    mean(sum(rates_ZF_epc, 1)), mean(sum(rates_ZF_sca, 1)), ...
    mean(sum(rates_CB_epc, 1)), mean(sum(rates_CB_sca, 1)));
