%% Simulation Script: Figure 4
% Reproduces Figure 4 from:
% Emil Bjornson, Luca Sanguinetti, "Making Cell-Free Massive MIMO
% Competitive With MMSE Processing and Centralized Implementation,"
% IEEE Trans. Wireless Commun., vol. 19, no. 2, pp. 1257-1271, Feb. 2020.
%
% Figure 4: Revisiting "Cell-free vs Small Cells" using Ngo et al.'s setup
%           Three-slope path loss model, L=100, N=1, K=40, tau_p=20
%
% Compares:
%   - Level 2 (MR) from Ngo et al. [16]
%   - Small cells with max-beta selection [16]
%   - Small cells with improved max-SE selection
%   - Level 1 with corrected SE (Proposition 3)
%   - Level 2 with L-MMSE combining
%
% Sub-figures:
%   (a) Full power transmission
%   (b) Max-min fairness power control

close all; clear; clc;

%% Simulation Parameters
L = 100;            % Number of APs (single-antenna)
K = 40;             % Number of UEs
tau_c = 200;        % Coherence block length
tau_p = 20;         % Pilot length (matching Ngo et al.)

nbrOfSetups = 100;  % Random setups
nbrOfRealizations = 200;  % For Monte Carlo in L-MMSE computation

p_max = 0.1;  % 100 mW max power
p = p_max * ones(K, 1);

fprintf('=== Simulation for Figure 4 ===\n');
fprintf('L = %d, K = %d, tau_p = %d (Ngo et al. setup)\n', L, K, tau_p);
fprintf('Three-slope path loss model\n');

%% Storage
SE_L2_MR_all = [];
SE_SC_maxBeta_all = [];
SE_SC_maxSE_all = [];
SE_L2_LMMSE_all = [];

%% Main simulation loop
for s = 1:nbrOfSetups
    if mod(s, 10) == 0
        fprintf('Setup %d/%d...\n', s, nbrOfSetups);
    end

    % Generate setup with three-slope model
    [gainOverNoisedB, pilotIndex, APpositions, UEpositions] = ...
        generateSetup_threeslope(L, K, tau_p);

    % --- Level 2 with MR combining (Ngo et al.'s approach) ---
    SE_L2_MR = functionComputeSE_CF_uplink_ngo(...
        gainOverNoisedB, tau_c, tau_p, pilotIndex, p, K, L);

    % --- Small cells: max-beta and max-SE selection ---
    [SE_SC_maxBeta, SE_SC_maxSE] = functionComputeSE_SC_uplink_ngo(...
        gainOverNoisedB, tau_c, tau_p, pilotIndex, p, K, L);

    % --- Level 2 with L-MMSE combining ---
    SE_L2_LMMSE = functionComputeSE_L2_LMMSE_ngo(...
        gainOverNoisedB, tau_c, tau_p, pilotIndex, p, K, L, nbrOfRealizations);

    % Store results
    SE_L2_MR_all = [SE_L2_MR_all; SE_L2_MR]; %#ok<AGROW>
    SE_SC_maxBeta_all = [SE_SC_maxBeta_all; SE_SC_maxBeta]; %#ok<AGROW>
    SE_SC_maxSE_all = [SE_SC_maxSE_all; SE_SC_maxSE]; %#ok<AGROW>
    SE_L2_LMMSE_all = [SE_L2_LMMSE_all; SE_L2_LMMSE]; %#ok<AGROW>
end

%% Figure 4(a): Full power transmission
figure('Name', 'Figure 4(a): Full Power', 'Position', [100 100 700 500]);

hold on;
cdfplot_custom(SE_L2_MR_all, 'b-', 2);
cdfplot_custom(SE_SC_maxBeta_all, 'r--', 2);
cdfplot_custom(SE_SC_maxSE_all, 'k-.', 2);
cdfplot_custom(SE_L2_LMMSE_all, 'm-', 2);
hold off;

xlabel('Spectral Efficiency [bit/s/Hz]', 'FontSize', 12);
ylabel('CDF', 'FontSize', 12);
title('Figure 4(a): Three-Slope Model, Full Power', 'FontSize', 13);
legend('L2 (MR) [Ngo et al.]', ...
    'Small Cells, max-\beta [Ngo et al.]', ...
    'Small Cells, max-SE (improved)', ...
    'L2 (L-MMSE)', ...
    'Location', 'southeast', 'FontSize', 10);
grid on;
xlim([0 10]);
set(gca, 'FontSize', 11);

saveas(gcf, 'Figure4a_ThreeSlope_FullPower.fig');
saveas(gcf, 'Figure4a_ThreeSlope_FullPower.png');


%% Figure 4(b): Max-min power control
fprintf('\nComputing max-min power control...\n');

SE_L2_MR_maxmin_all = [];

for s = 1:min(nbrOfSetups, 50)  % Fewer setups for power optimization
    if mod(s, 10) == 0
        fprintf('Power opt setup %d...\n', s);
    end

    [gainOverNoisedB, pilotIndex, ~, ~] = generateSetup_threeslope(L, K, tau_p);

    % Extract SINR terms for power optimization
    [signalCF, interferenceCF, noiseCF, ~, ~, ~] = ...
        functionSINRterms_uplink_ngo(gainOverNoisedB, tau_p, pilotIndex, K, L);

    % Max-min power optimization for Level 2 MR
    p_opt = functionPowerOptimization_maxmin(...
        signalCF, interferenceCF, noiseCF, K, p_max, tau_c, tau_p);

    % Compute SE with optimized powers
    SE_L2_MR_maxmin = functionComputeSE_CF_uplink_ngo(...
        gainOverNoisedB, tau_c, tau_p, pilotIndex, p_opt, K, L);

    SE_L2_MR_maxmin_all = [SE_L2_MR_maxmin_all; SE_L2_MR_maxmin]; %#ok<AGROW>
end

figure('Name', 'Figure 4(b): Max-Min Power Control', 'Position', [100 100 700 500]);

hold on;
cdfplot_custom(SE_L2_MR_all, 'b-', 2);
cdfplot_custom(SE_L2_MR_maxmin_all, 'b--', 2);
cdfplot_custom(SE_SC_maxSE_all, 'k-.', 2);
cdfplot_custom(SE_L2_LMMSE_all, 'm-', 2);
hold off;

xlabel('Spectral Efficiency [bit/s/Hz]', 'FontSize', 12);
ylabel('CDF', 'FontSize', 12);
title('Figure 4(b): Max-Min Power Control Comparison', 'FontSize', 13);
legend('L2 (MR), full power', 'L2 (MR), max-min power', ...
    'Small Cells (max-SE)', 'L2 (L-MMSE), full power', ...
    'Location', 'southeast', 'FontSize', 10);
grid on;
xlim([0 10]);
set(gca, 'FontSize', 11);

saveas(gcf, 'Figure4b_MaxMinPower.fig');
saveas(gcf, 'Figure4b_MaxMinPower.png');

fprintf('\nFigure 4 completed and saved.\n');


function cdfplot_custom(data, lineStyle, lineWidth)
    sortedData = sort(data);
    n = length(sortedData);
    cdf = (1:n) / n;
    plot(sortedData, cdf, lineStyle, 'LineWidth', lineWidth);
end
