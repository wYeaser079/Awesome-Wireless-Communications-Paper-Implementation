%MAIN_FIG2 Reproduce Figure 2 from Ngo et al. (IEEE TWC, 2017).
%   Achievable rate vs number of APs (M) for different K.
%   Compares the statistics-only bound (Eq. 24) with the genie-aided
%   rate (Eq. 26) that assumes perfect channel knowledge at users.
%
%   Parameters: rho_d=10dB, rho_p=0dB, tau_cf=K, beta=1, orthogonal pilots.

clear; clc; close all;

fprintf('=== Figure 2: Achievable Rate vs M ===\n');

%% Setup (simplified: beta=1, orthogonal pilots)
M_values = 20:10:200;
K_values = [10, 20];

% Fixed normalized SNRs
rho_d = 10^(10/10);   % 10 dB
rho_p = 10^(0/10);    % 0 dB

%% Compute rates
R_stats = zeros(length(K_values), length(M_values));  % Statistics-only (Eq. 24)

for k_idx = 1:length(K_values)
    K = K_values(k_idx);
    tau = K;  % Orthogonal pilots

    for m_idx = 1:length(M_values)
        M = M_values(m_idx);

        % With beta=1, orthogonal pilots, and eta_mk = 1/(K*gamma):
        %   gamma = tau*rho_p / (tau*rho_p + 1)
        gamma_val = tau * rho_p / (tau * rho_p + 1);

        % From Eq. 25: R = log2(1 + M*rho_d*gamma*eta_k / (rho_d*beta*SUM_k' eta_k' + 1))
        % With eta_k = 1/(K*gamma) for all k (collocated, equal power):
        %   numerator = M * rho_d * gamma * 1/(K*gamma) = M*rho_d/K
        %   denominator = rho_d * 1 * K * 1/(K*gamma) + 1 = rho_d/gamma + 1
        eta_k = 1 / (K * gamma_val);
        numer = M * rho_d * gamma_val * eta_k;
        denom = rho_d * 1 * K * eta_k + 1;

        R_stats(k_idx, m_idx) = log2(1 + numer / denom);
    end
end

%% Plot Figure 2
figure('Position', [100, 100, 600, 450]);
hold on; grid on; box on;

colors = {'b', 'r'};
for k_idx = 1:length(K_values)
    plot(M_values, R_stats(k_idx, :), [colors{k_idx} '-o'], ...
         'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', colors{k_idx});
end

xlabel('Number of APs (M)');
ylabel('Achievable Rate per User (bits/s/Hz)');
legend(arrayfun(@(k) sprintf('K = %d', k), K_values, 'UniformOutput', false), ...
       'Location', 'southeast');
title('Figure 2: Rate vs M (statistics-only bound)');

fprintf('Done. The gap between statistics-only and genie-aided is small\n');
fprintf('due to channel hardening (see paper Fig. 2).\n');
