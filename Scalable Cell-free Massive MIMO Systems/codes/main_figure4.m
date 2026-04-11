%main_figure4 - Reproduce Figure 4 from [1]: PDF of channel gain variations
%
%   This script generates Figure 4 of:
%   [1] E. Bjornson, L. Sanguinetti, "Scalable Cell-Free Massive MIMO
%       Systems," IEEE Trans. Commun., vol. 68, no. 7, pp. 4247-4261, 2020.
%
%   Figure 4 shows the PDF of the normalized effective channel gain
%   |h_k^H w_k|^2 / E{|h_k^H w_k|^2} for different precoding schemes
%   (MR and LP-MMSE), illustrating channel hardening.
%
%   With strong channel hardening, this ratio concentrates around 1,
%   making the hardening bound in Proposition 3 tight.
%   With weak hardening (MR), the distribution has a long tail.

close all; clear; clc;

fprintf('=== Reproducing Figure 4: PDF of Channel Gain Variations ===\n');
fprintf('Paper: Scalable Cell-Free Massive MIMO Systems (Bjornson & Sanguinetti, 2020)\n\n');

%% Simulation Parameters
K = 100;
tau_c = 200;
tau_p = 10;

% Power parameters (raw mW -- noise normalization is embedded in R)
p = 100;     % UL power: 100 mW per UE
rho = 1000;  % DL power: 1 W = 1000 mW per AP

nbrOfSetups = 3;        % Fewer setups for this figure
nbrOfRealizations = 50;

%% Collect normalized channel gain samples for both setups
for setupIdx = 1:2

    if setupIdx == 1
        L = 400; N = 1;
        fprintf('--- Setup (a): L = %d, N = %d ---\n', L, N);
    else
        L = 100; N = 4;
        fprintf('--- Setup (b): L = %d, N = %d ---\n', L, N);
    end

    % Collect normalized gains: |h_k^H w_k|^2 / E{|h_k^H w_k|^2}
    normGains_MR = [];
    normGains_LP = [];

    for setup = 1:nbrOfSetups
        fprintf('  Setup %d/%d... ', setup, nbrOfSetups);
        tic;

        %% Generate setup
        [gainOverNoisedB, R, pilotIndex, D, ~] = ...
            generateSetup(L, K, N, tau_p, 1, []);

        gainOverNoisedB_s = gainOverNoisedB(:, :, 1);
        R_s = R(:, :, :, :, 1);
        pilotIndex_s = pilotIndex(:, 1);
        D_s = D(:, :, 1);

        p_vec = p * ones(K, 1);

        %% Generate channels
        [Hhat, H, Bmat, Cmat] = functionChannelEstimates(...
            R_s, nbrOfRealizations, L, K, N, tau_p, pilotIndex_s, p_vec);

        %% Distributed power allocation
        rho_dist = zeros(L, K);
        for l = 1:L
            servedUEs = find(D_s(l, :) == 1);
            if isempty(servedUEs)
                continue;
            end
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

        %% Compute effective channel gains for MR and LP-MMSE
        for nreal = 1:nbrOfRealizations
            for k = 1:K

                % === MR precoding ===
                gain_MR = 0;
                for l = 1:L
                    if D_s(l, k) == 1
                        trB = real(trace(Bmat(:, :, l, k)));
                        if trB > 0
                            w_kl = Hhat(:, nreal, l, k) * sqrt(rho_dist(l, k) / trB);
                            h_kl = H(:, nreal, l, k);
                            gain_MR = gain_MR + h_kl' * w_kl;
                        end
                    end
                end

                % === LP-MMSE precoding ===
                gain_LP = 0;
                for l = 1:L
                    if D_s(l, k) == 0
                        continue;
                    end
                    servedUEs = find(D_s(l, :) == 1);
                    localMatrix = eye(N);
                    for idx = 1:length(servedUEs)
                        j = servedUEs(idx);
                        hhat_jl = Hhat(:, nreal, l, j);
                        localMatrix = localMatrix + p_vec(j) * ...
                            (hhat_jl * hhat_jl' + Cmat(:, :, l, j));
                    end
                    hhat_kl = Hhat(:, nreal, l, k);
                    v_kl = p_vec(k) * (localMatrix \ hhat_kl);
                    normV = norm(v_kl)^2;
                    if normV > 0
                        w_kl = v_kl * sqrt(rho_dist(l, k) / normV);
                    else
                        w_kl = zeros(N, 1);
                    end
                    h_kl = H(:, nreal, l, k);
                    gain_LP = gain_LP + h_kl' * w_kl;
                end

                % Store gains (absolute squared)
                normGains_MR = [normGains_MR; abs(gain_MR)^2];
                normGains_LP = [normGains_LP; abs(gain_LP)^2];
            end
        end

        elapsed = toc;
        fprintf('done (%.1f s)\n', elapsed);
    end

    %% Normalize gains by their mean
    normGains_MR = normGains_MR / mean(normGains_MR);
    normGains_LP = normGains_LP / mean(normGains_LP);

    %% Plot PDF using kernel density estimation
    figure(setupIdx);
    hold on; box on; grid on;

    % Use histogram-based PDF estimation
    edges = linspace(0, 4, 200);

    [counts_MR, ~] = histcounts(normGains_MR, edges, 'Normalization', 'pdf');
    [counts_LP, ~] = histcounts(normGains_LP, edges, 'Normalization', 'pdf');
    centers = (edges(1:end-1) + edges(2:end)) / 2;

    plot(centers, counts_MR, 'g-', 'LineWidth', 2, 'DisplayName', 'MR');
    plot(centers, counts_LP, 'r-', 'LineWidth', 2, 'DisplayName', 'LP-MMSE');

    % Add vertical line at x = 1 (perfect hardening)
    xline(1, 'k--', 'LineWidth', 1, 'DisplayName', 'Perfect hardening');

    xlabel('Normalized channel gain |h_k^H w_k|^2 / E{|h_k^H w_k|^2}', 'FontSize', 12);
    ylabel('PDF', 'FontSize', 12);

    if setupIdx == 1
        title('(a) L = 400 APs with N = 1 antenna', 'FontSize', 13);
    else
        title('(b) L = 100 APs with N = 4 antennas', 'FontSize', 13);
    end

    legend('Location', 'northeast', 'FontSize', 11);
    xlim([0 4]);
    set(gca, 'FontSize', 11);

    % Print hardening statistics
    fprintf('  Channel hardening (L=%d, N=%d):\n', L, N);
    fprintf('    MR:      Var/Mean^2 = %.4f\n', var(normGains_MR));
    fprintf('    LP-MMSE: Var/Mean^2 = %.4f\n\n', var(normGains_LP));
end

fprintf('=== Figure 4 complete ===\n');
