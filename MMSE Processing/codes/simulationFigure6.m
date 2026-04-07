%% Simulation Script: Figure 6
% Reproduces Figure 6 from:
% Emil Bjornson, Luca Sanguinetti, "Making Cell-Free Massive MIMO
% Competitive With MMSE Processing and Centralized Implementation,"
% IEEE Trans. Wireless Commun., vol. 19, no. 2, pp. 1257-1271, Feb. 2020.
%
% Figure 6: Fronthaul signaling load vs coherence block length tau_c
%           L = 100 APs, N = 4, K = 40, tau_p = 10
%
% This is a purely analytical plot (no Monte Carlo simulation needed).
% Shows complex scalars per AP per channel use for:
%   - Level 4: N (constant, independent of tau_c)
%   - Levels 2/3: (1 - tau_p/tau_c) * K (increases with tau_c)

close all; clear; clc;

%% Parameters
L = 100;
N = 4;
K = 40;
tau_p = 10;

% Range of coherence block lengths
tau_c_range = tau_p + 1 : 1 : 500;

fprintf('=== Simulation for Figure 6 ===\n');
fprintf('L = %d, N = %d, K = %d, tau_p = %d\n', L, N, K, tau_p);

%% Compute fronthaul load per AP per channel use

% Level 4: sends raw signal per channel use
% Total per AP per coherence block: tau_c * N
% Per channel use: N (constant)
fronthaul_L4 = N * ones(size(tau_c_range));

% Levels 2/3: send K local estimates per data symbol
% Total per AP per coherence block: (tau_c - tau_p) * K
% Per channel use: (tau_c - tau_p) / tau_c * K
fronthaul_L23 = (tau_c_range - tau_p) ./ tau_c_range * K;

% Crossover point: N = (tau_c - tau_p) / tau_c * K
% => tau_c * N = (tau_c - tau_p) * K
% => tau_c * (K - N) = tau_p * K
% => tau_c_cross = tau_p * K / (K - N) = 10 * 40 / 36 = 11.11
tau_c_crossover = tau_p * K / (K - N);

%% Plot Figure 6
figure('Name', 'Figure 6: Fronthaul Signaling Load', ...
    'Position', [100 100 700 500]);

hold on;
plot(tau_c_range, fronthaul_L4, 'b-', 'LineWidth', 2.5);
plot(tau_c_range, fronthaul_L23, 'r--', 'LineWidth', 2.5);

% Mark the asymptote
yline(K, 'k:', 'LineWidth', 1);
text(400, K + 1.5, sprintf('K = %d', K), 'FontSize', 11);

% Mark the crossover point
if tau_c_crossover > tau_p + 1 && tau_c_crossover < 500
    plot(tau_c_crossover, N, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
    text(tau_c_crossover + 10, N + 2, ...
        sprintf('Crossover: \\tau_c = %.1f', tau_c_crossover), 'FontSize', 10);
end

% Mark the operating point
tau_c_op = 200;
idx_op = find(tau_c_range == tau_c_op);
if ~isempty(idx_op)
    plot(tau_c_op, fronthaul_L4(idx_op), 'bs', 'MarkerSize', 12, 'MarkerFaceColor', 'b');
    plot(tau_c_op, fronthaul_L23(idx_op), 'r^', 'MarkerSize', 12, 'MarkerFaceColor', 'r');

    ratio = fronthaul_L23(idx_op) / fronthaul_L4(idx_op);
    text(tau_c_op + 10, (fronthaul_L4(idx_op) + fronthaul_L23(idx_op))/2, ...
        sprintf('Ratio = %.1fx', ratio), 'FontSize', 11, 'FontWeight', 'bold');
end

hold off;

xlabel('Coherence Block Length \tau_c', 'FontSize', 12);
ylabel('Complex Scalars per AP per Channel Use', 'FontSize', 12);
title(sprintf('Figure 6: Fronthaul Load (L=%d, N=%d, K=%d)', L, N, K), 'FontSize', 13);
legend('Level 4 (Centralized)', 'Levels 2/3 (Distributed)', ...
    'Location', 'east', 'FontSize', 11);
grid on;
xlim([tau_p + 1, 500]);
ylim([0, K + 5]);
set(gca, 'FontSize', 11);

%% Print analysis table
fprintf('\n=== Fronthaul Analysis ===\n');
fprintf('%-10s %-12s %-12s %-12s\n', 'tau_c', 'Level 4', 'Level 2/3', 'Ratio');
fprintf('%s\n', repmat('-', 1, 46));
for tc = [20 50 100 200 500]
    fh_L4 = N;
    fh_L23 = (tc - tau_p) / tc * K;
    fprintf('%-10d %-12.1f %-12.1f %-12.1fx\n', tc, fh_L4, fh_L23, fh_L23/fh_L4);
end

fprintf('\nCrossover point: tau_c = %.1f\n', tau_c_crossover);
fprintf('For tau_c > %.1f, Level 4 requires LESS fronthaul.\n', tau_c_crossover);
fprintf('Asymptotic ratio (tau_c -> inf): K/N = %d/%d = %.1f\n', K, N, K/N);

saveas(gcf, 'Figure6_Fronthaul.fig');
saveas(gcf, 'Figure6_Fronthaul.png');

fprintf('\nFigure 6 completed and saved.\n');
