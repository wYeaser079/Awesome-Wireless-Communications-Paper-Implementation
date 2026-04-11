%% Simulation Script: Figure 5
% Reproduces Figure 5 from:
% Emil Bjornson, Luca Sanguinetti, "Making Cell-Free Massive MIMO
% Competitive With MMSE Processing and Centralized Implementation,"
% IEEE Trans. Wireless Commun., vol. 19, no. 2, pp. 1257-1271, Feb. 2020.
%
% Figure 5: CDF of SUM SE comparing MMSE-SIC vs linear combining
%           L = 400 APs, N = 1, K = 40 UEs
%
% Compares:
%   - Level 4 MMSE-SIC (Proposition 4)
%   - Level 4 MMSE (linear)
%   - Level 3 L-MMSE + LSFD
%   - Level 2 L-MMSE
%   - Level 1 L-MMSE (small cells)

close all; clear; clc;

%% Simulation Parameters
L = 400;
N = 1;
K = 40;
tau_c = 200;
tau_p = 10;
nbrOfBSs = 4;
M = L * N / nbrOfBSs;

nbrOfSetups = 50;           % Random setups
nbrOfRealizations = 200;    % Channel realizations

p_UE = 0.1;
p = p_UE * ones(K, 1);

fprintf('=== Simulation for Figure 5 ===\n');
fprintf('L = %d, N = %d, K = %d\n', L, N, K);
fprintf('Comparing MMSE-SIC vs linear combining (Sum SE)\n');

%% Storage
sumSE_SIC_all = zeros(nbrOfSetups, 1);
sumSE_L4_MMSE_all = zeros(nbrOfSetups, 1);
sumSE_L3_all = zeros(nbrOfSetups, 1);
sumSE_L2_all = zeros(nbrOfSetups, 1);
sumSE_L1_all = zeros(nbrOfSetups, 1);

%% Main loop
for s = 1:nbrOfSetups
    fprintf('[Fig5] Setup %d/%d -- realizations done: %d/%d\n', ...
        s, nbrOfSetups, (s-1)*nbrOfRealizations, nbrOfSetups*nbrOfRealizations);

    [R_AP, ~, ~, ~, pilotIndex, ~, ~, ~] = ...
        generateSetup(L, K, N, tau_p, nbrOfBSs, M);

    [Hhat, H, B, C] = functionChannelEstimates(...
        R_AP, nbrOfRealizations, L, K, N, tau_p, pilotIndex, p);

    [SE_MR, SE_MMSE, sumSE_SIC] = functionComputeSE_AP_uplink(...
        Hhat, H, B, C, R_AP, tau_c, tau_p, nbrOfRealizations, L, K, N, p);

    % Sum SE across all UEs
    sumSE_SIC_all(s) = sumSE_SIC;
    sumSE_L4_MMSE_all(s) = sum(SE_MMSE(:, 1));
    sumSE_L3_all(s) = sum(SE_MMSE(:, 2));
    sumSE_L2_all(s) = sum(SE_MMSE(:, 3));
    sumSE_L1_all(s) = sum(SE_MMSE(:, 4));
end

%% Plot Figure 5
figure('Name', 'Figure 5: Sum SE Comparison', 'Position', [100 100 700 500]);

hold on;
cdfplot_custom(sumSE_SIC_all, 'b-', 2.5);
cdfplot_custom(sumSE_L4_MMSE_all, 'r--', 2);
cdfplot_custom(sumSE_L3_all, 'k-.', 2);
cdfplot_custom(sumSE_L2_all, 'm:', 2);
cdfplot_custom(sumSE_L1_all, 'g-', 1.5);
hold off;

xlabel('Sum Spectral Efficiency [bit/s/Hz]', 'FontSize', 12);
ylabel('CDF', 'FontSize', 12);
title('Figure 5: Sum SE, L=400, N=1', 'FontSize', 13);
legend('Level 4 (MMSE-SIC)', 'Level 4 (MMSE)', ...
    'Level 3 (L-MMSE + LSFD)', 'Level 2 (L-MMSE)', ...
    'Level 1 (L-MMSE, Small Cells)', ...
    'Location', 'northwest', 'FontSize', 10);
grid on;
set(gca, 'FontSize', 11);

% Print summary statistics
fprintf('\n=== Sum SE Summary (bit/s/Hz) ===\n');
fprintf('Level 4 MMSE-SIC:  Mean = %.1f\n', mean(sumSE_SIC_all));
fprintf('Level 4 MMSE:      Mean = %.1f\n', mean(sumSE_L4_MMSE_all));
fprintf('Level 3 L-MMSE:    Mean = %.1f\n', mean(sumSE_L3_all));
fprintf('Level 2 L-MMSE:    Mean = %.1f\n', mean(sumSE_L2_all));
fprintf('Level 1 L-MMSE:    Mean = %.1f\n', mean(sumSE_L1_all));
fprintf('SIC gain over MMSE: %.1f%%\n', ...
    100*(mean(sumSE_SIC_all) - mean(sumSE_L4_MMSE_all)) / mean(sumSE_L4_MMSE_all));

saveas(gcf, 'Figure5_SumSE_SIC.fig');
saveas(gcf, 'Figure5_SumSE_SIC.png');

fprintf('\nFigure 5 completed and saved.\n');


function cdfplot_custom(data, lineStyle, lineWidth)
    sortedData = sort(data);
    n = length(sortedData);
    cdf = (1:n) / n;
    plot(sortedData, cdf, lineStyle, 'LineWidth', lineWidth);
end
