%main_figure5 - Reproduce Figure 5 from [1]: UL SE CDF comparison
%
%   This script generates Figure 5(a) and Figure 5(b) of:
%   [1] E. Bjornson, L. Sanguinetti, "Scalable Cell-Free Massive MIMO
%       Systems," IEEE Trans. Commun., vol. 68, no. 7, pp. 4247-4261, 2020.
%
%   Figure 5 shows the CDF of uplink spectral efficiency per UE comparing:
%     - Scalable P-MMSE (centralized, proposed)
%     - Scalable LP-MMSE (distributed, proposed)
%     - Non-scalable MMSE (All) from [7]
%     - Non-scalable L-MMSE (All) from [7]
%     - Non-scalable MR (All) from [5]
%
%   Two setups: (a) L=400, N=1  and  (b) L=100, N=4
%
%   IMPORTANT: This simulation is computationally intensive.
%   Reduce nbrOfSetups or nbrOfRealizations for faster (approximate) results.

close all; clear; clc;

fprintf('=== Reproducing Figure 5: UL SE CDF ===\n');
fprintf('Paper: Scalable Cell-Free Massive MIMO Systems (Bjornson & Sanguinetti, 2020)\n\n');

%% Simulation Parameters (matching the paper)
K = 100;            % Number of UEs
tau_c = 200;        % Coherence block length (samples)
tau_p = 10;         % Pilot sequence length

% UL transmit power: p_k = 100 mW = 20 dBm
% Note: R matrices from generateSetup are already noise-normalized
% (gainOverNoisedB subtracts noise floor), so powers are in raw mW.
p = 100; % 100 mW (raw value, noise normalization is in R)

% Simulation averaging
nbrOfSetups = 25;         % Number of random network layouts
nbrOfRealizations = 1000; % Channel realizations per setup

% NOTE: For faster testing, use nbrOfSetups = 5 and nbrOfRealizations = 100

%% Run simulations for both setups
for setupIdx = 1:2

    if setupIdx == 1
        L = 400; N = 1;
        fprintf('--- Setup (a): L = %d APs, N = %d antenna ---\n', L, N);
    else
        L = 100; N = 4;
        fprintf('--- Setup (b): L = %d APs, N = %d antennas ---\n', L, N);
    end

    % Accumulate SE results across all setups
    SE_MR_all = [];
    SE_LP_MMSE_all = [];
    SE_P_MMSE_all = [];
    SE_MMSE_all_scalable = [];
    SE_MMSE_all_nonscalable = [];

    % Non-scalable benchmarks (all APs serve all UEs)
    SE_MR_nonscalable = [];
    SE_LMMSE_nonscalable = [];
    SE_MMSE_nonscalable_full = [];

    for setup = 1:nbrOfSetups
        fprintf('  Setup %d/%d (L=%d, N=%d)... ', setup, nbrOfSetups, L, N);
        tic;

        %% Generate network layout
        [gainOverNoisedB, R, pilotIndex, D, D_all] = ...
            generateSetup(L, K, N, tau_p, 1, []);

        % Extract single setup
        gainOverNoisedB_s = gainOverNoisedB(:, :, 1);
        R_s = R(:, :, :, :, 1);
        pilotIndex_s = pilotIndex(:, 1);
        D_s = D(:, :, 1);

        % Power vector (same for all UEs)
        p_vec = p * ones(K, 1);

        %% Generate channel estimates
        [Hhat, H, Bmat, Cmat] = functionChannelEstimates(...
            R_s, nbrOfRealizations, L, K, N, tau_p, pilotIndex_s, p_vec);

        %% Compute UL SE with scalable DCC (proposed)
        [SE_MR_s, SE_LP_s, SE_P_s, ~] = functionComputeSE_uplink(...
            Hhat, H, D_s, Bmat, Cmat, tau_c, tau_p, nbrOfRealizations, ...
            N, K, L, p_vec, R_s, pilotIndex_s);

        %% Compute UL SE with all-serve-all (non-scalable benchmarks)
        [SE_MR_ns, SE_LP_ns, ~, SE_MMSE_ns] = functionComputeSE_uplink(...
            Hhat, H, D_all, Bmat, Cmat, tau_c, tau_p, nbrOfRealizations, ...
            N, K, L, p_vec, R_s, pilotIndex_s);

        % Accumulate results
        SE_MR_all = [SE_MR_all; SE_MR_s];
        SE_LP_MMSE_all = [SE_LP_MMSE_all; SE_LP_s];
        SE_P_MMSE_all = [SE_P_MMSE_all; SE_P_s];

        SE_MR_nonscalable = [SE_MR_nonscalable; SE_MR_ns];
        SE_LMMSE_nonscalable = [SE_LMMSE_nonscalable; SE_LP_ns];
        SE_MMSE_nonscalable_full = [SE_MMSE_nonscalable_full; SE_MMSE_ns];

        elapsed = toc;
        fprintf('done (%.1f s)\n', elapsed);
    end

    %% Plot CDF
    figure(setupIdx);
    hold on; box on; grid on;

    % Sort SE values for CDF plotting
    plot(sort(SE_MMSE_nonscalable_full), linspace(0, 1, length(SE_MMSE_nonscalable_full)), ...
        'k-', 'LineWidth', 2, 'DisplayName', 'MMSE (All)');
    plot(sort(SE_P_MMSE_all), linspace(0, 1, length(SE_P_MMSE_all)), ...
        'b--', 'LineWidth', 2, 'DisplayName', 'P-MMSE (Scalable)');
    plot(sort(SE_LMMSE_nonscalable), linspace(0, 1, length(SE_LMMSE_nonscalable)), ...
        'r-', 'LineWidth', 2, 'DisplayName', 'L-MMSE (All)');
    plot(sort(SE_LP_MMSE_all), linspace(0, 1, length(SE_LP_MMSE_all)), ...
        'm--', 'LineWidth', 2, 'DisplayName', 'LP-MMSE (Scalable)');
    plot(sort(SE_MR_nonscalable), linspace(0, 1, length(SE_MR_nonscalable)), ...
        'g-', 'LineWidth', 1.5, 'DisplayName', 'MR (All)');

    xlabel('SE per UE [bit/s/Hz]', 'FontSize', 12);
    ylabel('CDF', 'FontSize', 12);

    if setupIdx == 1
        title('(a) L = 400 APs with N = 1 antenna', 'FontSize', 13);
    else
        title('(b) L = 100 APs with N = 4 antennas', 'FontSize', 13);
    end

    legend('Location', 'southeast', 'FontSize', 10);
    xlim([0 10]);
    ylim([0 1]);
    set(gca, 'FontSize', 11);

    % Print summary statistics
    fprintf('\n  Results Summary (L=%d, N=%d):\n', L, N);
    fprintf('    MMSE (All):        Avg SE = %.2f bit/s/Hz\n', mean(SE_MMSE_nonscalable_full));
    fprintf('    P-MMSE (Scalable): Avg SE = %.2f bit/s/Hz (%.0f%% of MMSE)\n', ...
        mean(SE_P_MMSE_all), 100*mean(SE_P_MMSE_all)/mean(SE_MMSE_nonscalable_full));
    fprintf('    L-MMSE (All):      Avg SE = %.2f bit/s/Hz\n', mean(SE_LMMSE_nonscalable));
    fprintf('    LP-MMSE (Scalable):Avg SE = %.2f bit/s/Hz\n', mean(SE_LP_MMSE_all));
    fprintf('    MR (All):          Avg SE = %.2f bit/s/Hz\n', mean(SE_MR_nonscalable));
    fprintf('    LP-MMSE/MR ratio:  %.1fx\n\n', mean(SE_LP_MMSE_all)/mean(SE_MR_nonscalable));

end

fprintf('=== Figure 5 complete ===\n');
