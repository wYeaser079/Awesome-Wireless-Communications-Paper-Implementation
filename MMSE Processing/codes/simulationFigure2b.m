%% Simulation Script: Figure 2(b)
% Reproduces Figure 2(b) from:
% Emil Bjornson, Luca Sanguinetti, "Making Cell-Free Massive MIMO
% Competitive With MMSE Processing and Centralized Implementation,"
% IEEE Trans. Wireless Commun., vol. 19, no. 2, pp. 1257-1271, Feb. 2020.
%
% Figure 2(b): CDF of per-user SE with MMSE/L-MMSE combining
%              L = 100 APs, N = 4 antennas each, K = 40 UEs
%              Compared against cellular mMIMO with 4 BSs, 100 antennas each

close all; clear; clc;

%% Simulation Parameters
L = 100;            % Number of APs
N = 4;              % Antennas per AP
K = 40;             % Number of UEs
tau_c = 200;        % Coherence block length
tau_p = 10;         % Pilot length
nbrOfBSs = 4;      % Number of BSs
M = L * N / nbrOfBSs;  % Antennas per BS = 100

nbrOfSetups = 30;           % Random setups (increase for smoother curves)
nbrOfRealizations = 100;    % Channel realizations per setup

p_UE = 0.1;  % 100 mW
p = p_UE * ones(K, 1);

fprintf('=== Simulation for Figure 2(b) ===\n');
fprintf('L = %d APs, N = %d, K = %d UEs\n', L, N, K);
fprintf('Note: N=4 is computationally heavier than N=1.\n');

%% Storage
SE_MMSE_all = zeros(K, 4, nbrOfSetups);
SE_cellular_MMMSE_all = zeros(K, nbrOfSetups);

%% Main loop
for s = 1:nbrOfSetups
    fprintf('Setup %d/%d...\n', s, nbrOfSetups);

    [R_AP, R_BS, gainOverNoisedB_AP, gainOverNoisedB_BS, ...
        pilotIndex, APpositions, BSpositions, UEpositions] = ...
        generateSetup(L, K, N, tau_p, nbrOfBSs, M);

    [Hhat, H, B, C] = functionChannelEstimates(R_AP, nbrOfRealizations, L, K, N, tau_p, pilotIndex, p);

    [~, SE_MMSE, ~] = functionComputeSE_AP_uplink(...
        Hhat, H, B, C, R_AP, tau_c, tau_p, nbrOfRealizations, L, K, N, p);

    SE_MMSE_all(:, :, s) = SE_MMSE;

    [~, SE_cellular_MMMSE] = functionComputeSE_BS_uplink(...
        R_BS, gainOverNoisedB_BS, tau_c, tau_p, nbrOfRealizations, ...
        nbrOfBSs, M, K, K/nbrOfBSs, p_UE);

    SE_cellular_MMMSE_all(:, s) = SE_cellular_MMMSE;
end

%% Plot Figure 2(b)
figure('Name', 'Figure 2(b): MMSE/L-MMSE (L=100, N=4)', ...
    'Position', [100 100 700 500]);

SE_L4 = reshape(SE_MMSE_all(:, 1, :), [], 1);
SE_L3 = reshape(SE_MMSE_all(:, 2, :), [], 1);
SE_L2 = reshape(SE_MMSE_all(:, 3, :), [], 1);
SE_L1 = reshape(SE_MMSE_all(:, 4, :), [], 1);
SE_cell = reshape(SE_cellular_MMMSE_all, [], 1);

hold on;
cdfplot_custom(SE_L4, 'b-', 2);
cdfplot_custom(SE_L3, 'r--', 2);
cdfplot_custom(SE_L2, 'k-.', 2);
cdfplot_custom(SE_L1, 'm:', 2);
cdfplot_custom(SE_cell, 'g-', 2);
hold off;

xlabel('Spectral Efficiency [bit/s/Hz]', 'FontSize', 12);
ylabel('CDF', 'FontSize', 12);
title('Figure 2(b): L=100, N=4, MMSE/L-MMSE Combining', 'FontSize', 13);
legend('Level 4 (MMSE)', 'Level 3 (L-MMSE + LSFD)', ...
    'Level 2 (L-MMSE)', 'Level 1 (L-MMSE, Small Cells)', ...
    'Cellular (M-MMSE)', 'Location', 'southeast', 'FontSize', 10);
grid on;
xlim([0 12]);
set(gca, 'FontSize', 11);

saveas(gcf, 'Figure2b_MMSE_L100_N4.fig');
saveas(gcf, 'Figure2b_MMSE_L100_N4.png');

fprintf('\nFigure 2(b) completed and saved.\n');


function cdfplot_custom(data, lineStyle, lineWidth)
    sortedData = sort(data);
    n = length(sortedData);
    cdf = (1:n) / n;
    plot(sortedData, cdf, lineStyle, 'LineWidth', lineWidth);
end
