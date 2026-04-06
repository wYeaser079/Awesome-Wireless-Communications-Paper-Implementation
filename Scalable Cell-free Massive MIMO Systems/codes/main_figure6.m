%main_figure6 - Reproduce Figure 6 from [1]: DL SE CDF comparison
%
%   This script generates Figure 6(a) and Figure 6(b) of:
%   [1] E. Bjornson, L. Sanguinetti, "Scalable Cell-Free Massive MIMO
%       Systems," IEEE Trans. Commun., vol. 68, no. 7, pp. 4247-4261, 2020.
%
%   Figure 6 shows the CDF of downlink spectral efficiency per UE comparing:
%     - Centralized P-MMSE (with and without perfect CSI at UE)
%     - Distributed LP-MMSE (with and without perfect CSI at UE)
%     - Distributed MR (with and without perfect CSI at UE)
%
%   Two setups: (a) L=400, N=1  and  (b) L=100, N=4
%
%   IMPORTANT: This simulation is computationally intensive.
%   Reduce nbrOfSetups or nbrOfRealizations for faster (approximate) results.

close all; clear; clc;

fprintf('=== Reproducing Figure 6: DL SE CDF ===\n');
fprintf('Paper: Scalable Cell-Free Massive MIMO Systems (Bjornson & Sanguinetti, 2020)\n\n');

%% Simulation Parameters
K = 100;            % Number of UEs
tau_c = 200;        % Coherence block length
tau_p = 10;         % Pilot sequence length

% Power parameters (raw mW -- noise normalization is embedded in R)
p = 100;     % UL power: 100 mW per UE
rho = 1000;  % DL power: 1 W = 1000 mW per AP

% Simulation averaging
nbrOfSetups = 25;
nbrOfRealizations = 1000;

%% Run simulations for both setups
for setupIdx = 1:2

    if setupIdx == 1
        L = 400; N = 1;
        fprintf('--- Setup (a): L = %d APs, N = %d antenna ---\n', L, N);
    else
        L = 100; N = 4;
        fprintf('--- Setup (b): L = %d APs, N = %d antennas ---\n', L, N);
    end

    % Accumulate SE results
    SE_MR_all = []; SE_MR_perf_all = [];
    SE_LP_all = []; SE_LP_perf_all = [];
    SE_PM_all = []; SE_PM_perf_all = [];

    for setup = 1:nbrOfSetups
        fprintf('  Setup %d/%d (L=%d, N=%d)... ', setup, nbrOfSetups, L, N);
        tic;

        %% Generate network layout
        [gainOverNoisedB, R, pilotIndex, D, D_all] = ...
            generateSetup(L, K, N, tau_p, 1, []);

        gainOverNoisedB_s = gainOverNoisedB(:, :, 1);
        R_s = R(:, :, :, :, 1);
        pilotIndex_s = pilotIndex(:, 1);
        D_s = D(:, :, 1);

        p_vec = p * ones(K, 1);

        %% Generate channel estimates
        [Hhat, H, Bmat, Cmat] = functionChannelEstimates(...
            R_s, nbrOfRealizations, L, K, N, tau_p, pilotIndex_s, p_vec);

        %% Compute DL power allocation

        % --- Distributed power allocation (Eq. 43 from [1]) ---
        % rho_{kl} = rho * sqrt(beta_{kl}) / sum_{i in D_l} sqrt(beta_{il})
        rho_dist = zeros(L, K);
        for l = 1:L
            servedUEs = find(D_s(l, :) == 1);
            if isempty(servedUEs)
                continue;
            end
            % Compute sqrt(beta) for each served UE
            sqrtBeta = zeros(length(servedUEs), 1);
            for idx = 1:length(servedUEs)
                k = servedUEs(idx);
                sqrtBeta(idx) = sqrt(db2pow(gainOverNoisedB_s(l, k)));
            end
            sumSqrtBeta = sum(sqrtBeta);
            for idx = 1:length(servedUEs)
                k = servedUEs(idx);
                rho_dist(l, k) = rho * sqrtBeta(idx) / sumSqrtBeta;
            end
        end

        % --- Centralized power allocation ---
        % Equal power: rho_k = rho / tau_p for each UE
        rho_central = (rho / tau_p) * ones(K, 1);

        %% Compute DL SE
        [SE_MR_s, SE_LP_s, SE_MR_perf_s, SE_LP_perf_s, SE_PM_s, SE_PM_perf_s] = ...
            functionComputeSE_downlink(Hhat, H, D_s, Bmat, Cmat, tau_c, tau_p, ...
            nbrOfRealizations, N, K, L, p_vec, rho_dist, R_s, pilotIndex_s, rho_central);

        % Accumulate
        SE_MR_all = [SE_MR_all; SE_MR_s];
        SE_MR_perf_all = [SE_MR_perf_all; SE_MR_perf_s];
        SE_LP_all = [SE_LP_all; SE_LP_s];
        SE_LP_perf_all = [SE_LP_perf_all; SE_LP_perf_s];
        SE_PM_all = [SE_PM_all; SE_PM_s];
        SE_PM_perf_all = [SE_PM_perf_all; SE_PM_perf_s];

        elapsed = toc;
        fprintf('done (%.1f s)\n', elapsed);
    end

    %% Plot CDF
    figure(setupIdx);
    hold on; box on; grid on;

    % P-MMSE
    plot(sort(SE_PM_perf_all), linspace(0, 1, length(SE_PM_perf_all)), ...
        'b:', 'LineWidth', 1.5, 'DisplayName', 'P-MMSE (perfect CSI)');
    plot(sort(SE_PM_all), linspace(0, 1, length(SE_PM_all)), ...
        'b-', 'LineWidth', 2, 'DisplayName', 'P-MMSE');

    % LP-MMSE
    plot(sort(SE_LP_perf_all), linspace(0, 1, length(SE_LP_perf_all)), ...
        'r:', 'LineWidth', 1.5, 'DisplayName', 'LP-MMSE (perfect CSI)');
    plot(sort(SE_LP_all), linspace(0, 1, length(SE_LP_all)), ...
        'r-', 'LineWidth', 2, 'DisplayName', 'LP-MMSE');

    % MR
    plot(sort(SE_MR_perf_all), linspace(0, 1, length(SE_MR_perf_all)), ...
        'g:', 'LineWidth', 1.5, 'DisplayName', 'MR (perfect CSI)');
    plot(sort(SE_MR_all), linspace(0, 1, length(SE_MR_all)), ...
        'g-', 'LineWidth', 2, 'DisplayName', 'MR');

    xlabel('SE per UE [bit/s/Hz]', 'FontSize', 12);
    ylabel('CDF', 'FontSize', 12);

    if setupIdx == 1
        title('(a) L = 400 APs with N = 1 antenna', 'FontSize', 13);
    else
        title('(b) L = 100 APs with N = 4 antennas', 'FontSize', 13);
    end

    legend('Location', 'southeast', 'FontSize', 9);
    xlim([0 10]);
    ylim([0 1]);
    set(gca, 'FontSize', 11);

    % Print summary
    fprintf('\n  Results Summary (L=%d, N=%d):\n', L, N);
    fprintf('    P-MMSE:          Avg SE = %.2f bit/s/Hz\n', mean(SE_PM_all));
    fprintf('    P-MMSE (perf):   Avg SE = %.2f bit/s/Hz (bound tightness: %.0f%%)\n', ...
        mean(SE_PM_perf_all), 100*mean(SE_PM_all)/mean(SE_PM_perf_all));
    fprintf('    LP-MMSE:         Avg SE = %.2f bit/s/Hz\n', mean(SE_LP_all));
    fprintf('    LP-MMSE (perf):  Avg SE = %.2f bit/s/Hz (bound tightness: %.0f%%)\n', ...
        mean(SE_LP_perf_all), 100*mean(SE_LP_all)/mean(SE_LP_perf_all));
    fprintf('    MR:              Avg SE = %.2f bit/s/Hz\n', mean(SE_MR_all));
    fprintf('    MR (perf):       Avg SE = %.2f bit/s/Hz (bound tightness: %.0f%%)\n\n', ...
        mean(SE_MR_perf_all), 100*mean(SE_MR_all)/mean(SE_MR_perf_all));
end

fprintf('=== Figure 6 complete ===\n');
