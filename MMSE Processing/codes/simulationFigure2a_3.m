%% Simulation Script: Figures 2(a) and 3
% Reproduces Figures 2(a) and 3 from:
% Emil Bjornson, Luca Sanguinetti, "Making Cell-Free Massive MIMO
% Competitive With MMSE Processing and Centralized Implementation,"
% IEEE Trans. Wireless Commun., vol. 19, no. 2, pp. 1257-1271, Feb. 2020.
%
% Figure 2(a): CDF of per-user SE with MMSE/L-MMSE combining
%              L = 400 APs, N = 1 antenna each, K = 40 UEs
% Figure 3:   CDF of per-user SE with MR combining (same setup)
%
% Both figures compare cell-free Levels 1-4 against cellular M-MMSE.

close all; clear; clc;

%% Simulation Parameters
L = 400;            % Number of APs (cell-free)
N = 1;              % Antennas per AP
K = 40;             % Number of UEs
tau_c = 200;        % Coherence block length (samples)
tau_p = 10;         % Pilot length
nbrOfBSs = 4;      % Number of BSs (cellular)
M = L * N / nbrOfBSs;  % Antennas per BS = 100

nbrOfSetups = 50;          % Number of random setups (increase to 200 for smoother curves)
nbrOfRealizations = 200;   % Channel realizations per setup (increase to 1000)

% Power settings
p_UE = 0.1;         % 100 mW transmit power per UE (in Watts)
p = p_UE * ones(K, 1);

fprintf('=== Simulation for Figures 2(a) and 3 ===\n');
fprintf('L = %d APs, N = %d, K = %d UEs, tau_c = %d, tau_p = %d\n', L, N, K, tau_c, tau_p);
fprintf('Cellular: %d BSs, M = %d antennas each\n', nbrOfBSs, M);
fprintf('Number of setups: %d, realizations per setup: %d\n', nbrOfSetups, nbrOfRealizations);

%% Storage for results
SE_MR_all = zeros(K, 4, nbrOfSetups);      % [L4, L3, L2, L1] x setups
SE_MMSE_all = zeros(K, 4, nbrOfSetups);
SE_cellular_MR_all = zeros(K, nbrOfSetups);
SE_cellular_MMMSE_all = zeros(K, nbrOfSetups);

%% Main simulation loop
for s = 1:nbrOfSetups
    fprintf('Setup %d/%d...\n', s, nbrOfSetups);

    % Generate network topology and channel statistics
    [R_AP, R_BS, gainOverNoisedB_AP, gainOverNoisedB_BS, ...
        pilotIndex, APpositions, BSpositions, UEpositions] = ...
        generateSetup(L, K, N, tau_p, nbrOfBSs, M);

    % Generate channel estimates for cell-free system
    [Hhat, H, B, C] = functionChannelEstimates(R_AP, nbrOfRealizations, L, K, N, tau_p, pilotIndex, p);

    % Compute SE for cell-free system (all 4 levels, MR and MMSE)
    [SE_MR, SE_MMSE, ~] = functionComputeSE_AP_uplink(...
        Hhat, H, B, C, R_AP, tau_c, tau_p, nbrOfRealizations, L, K, N, p);

    SE_MR_all(:, :, s) = SE_MR;
    SE_MMSE_all(:, :, s) = SE_MMSE;

    % Compute SE for cellular Massive MIMO
    [SE_cellular_MR, SE_cellular_MMMSE] = functionComputeSE_BS_uplink(...
        R_BS, gainOverNoisedB_BS, tau_c, tau_p, nbrOfRealizations, ...
        nbrOfBSs, M, K, K/nbrOfBSs, p_UE);

    SE_cellular_MR_all(:, s) = SE_cellular_MR;
    SE_cellular_MMMSE_all(:, s) = SE_cellular_MMMSE;
end

%% Plot Figure 2(a): MMSE/L-MMSE combining
figure('Name', 'Figure 2(a): MMSE/L-MMSE Combining (L=400, N=1)', ...
    'Position', [100 100 700 500]);

% Collect all per-user SEs across setups
SE_L4_MMSE = reshape(SE_MMSE_all(:, 1, :), [], 1);
SE_L3_MMSE = reshape(SE_MMSE_all(:, 2, :), [], 1);
SE_L2_MMSE = reshape(SE_MMSE_all(:, 3, :), [], 1);
SE_L1_MMSE = reshape(SE_MMSE_all(:, 4, :), [], 1);
SE_cell_MMMSE = reshape(SE_cellular_MMMSE_all, [], 1);

hold on;
cdfplot_custom(SE_L4_MMSE, 'b-', 2);
cdfplot_custom(SE_L3_MMSE, 'r--', 2);
cdfplot_custom(SE_L2_MMSE, 'k-.', 2);
cdfplot_custom(SE_L1_MMSE, 'm:', 2);
cdfplot_custom(SE_cell_MMMSE, 'g-', 2);
hold off;

xlabel('Spectral Efficiency [bit/s/Hz]', 'FontSize', 12);
ylabel('CDF', 'FontSize', 12);
title('Figure 2(a): L=400, N=1, MMSE/L-MMSE Combining', 'FontSize', 13);
legend('Level 4 (MMSE)', 'Level 3 (L-MMSE + LSFD)', ...
    'Level 2 (L-MMSE)', 'Level 1 (L-MMSE, Small Cells)', ...
    'Cellular (M-MMSE)', 'Location', 'southeast', 'FontSize', 10);
grid on;
xlim([0 12]);
set(gca, 'FontSize', 11);

saveas(gcf, 'Figure2a_MMSE_L400_N1.fig');
saveas(gcf, 'Figure2a_MMSE_L400_N1.png');

%% Plot Figure 3: MR combining
figure('Name', 'Figure 3: MR Combining (L=400, N=1)', ...
    'Position', [100 100 700 500]);

SE_L4_MR = reshape(SE_MR_all(:, 1, :), [], 1);
SE_L3_MR = reshape(SE_MR_all(:, 2, :), [], 1);
SE_L2_MR = reshape(SE_MR_all(:, 3, :), [], 1);
SE_L1_MR = reshape(SE_MR_all(:, 4, :), [], 1);
SE_cell_MR = reshape(SE_cellular_MR_all, [], 1);

hold on;
cdfplot_custom(SE_L4_MR, 'b-', 2);
cdfplot_custom(SE_L3_MR, 'r--', 2);
cdfplot_custom(SE_L2_MR, 'k-.', 2);
cdfplot_custom(SE_L1_MR, 'm:', 2);
cdfplot_custom(SE_cell_MMMSE, 'g-', 2);  % Cellular always uses M-MMSE
hold off;

xlabel('Spectral Efficiency [bit/s/Hz]', 'FontSize', 12);
ylabel('CDF', 'FontSize', 12);
title('Figure 3: L=400, N=1, MR Combining', 'FontSize', 13);
legend('Level 4 (MR)', 'Level 3 (MR + LSFD)', ...
    'Level 2 (MR)', 'Level 1 (MR, Small Cells)', ...
    'Cellular (M-MMSE)', 'Location', 'southeast', 'FontSize', 10);
grid on;
xlim([0 12]);
set(gca, 'FontSize', 11);

saveas(gcf, 'Figure3_MR_L400_N1.fig');
saveas(gcf, 'Figure3_MR_L400_N1.png');

fprintf('\nFigures 2(a) and 3 completed and saved.\n');


%% Helper function for CDF plotting
function cdfplot_custom(data, lineStyle, lineWidth)
    sortedData = sort(data);
    n = length(sortedData);
    cdf = (1:n) / n;
    plot(sortedData, cdf, lineStyle, 'LineWidth', lineWidth);
end
