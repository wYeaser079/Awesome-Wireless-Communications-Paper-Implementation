%MAIN_FIG9 Reproduce Figure 9 from Ngo et al. (IEEE TWC, 2017).
%   CDF of the effective number of APs serving each user.
%   The effective number is the minimum number of APs contributing >= 95%
%   of the total power allocated to a given user (Eq. 58).
%
%   Parameters: M=100, K=40, tau_cf = 5 and 20, uncorrelated shadowing.

clear; clc; close all;

%% Setup
p = params();

fprintf('=== Figure 9: Effective Number of APs per User ===\n');
fprintf('M = %d, K = %d\n', p.M, p.K);
fprintf('Number of setups: %d\n\n', p.num_setups);

tau_values = [5, 20];
eff_APs = cell(length(tau_values), 1);

%% Main simulation loop
for t_idx = 1:length(tau_values)
    p.tau_cf = tau_values(t_idx);

    % Update pre-log factor
    p.prelog_cf = (1 - p.tau_cf / p.tau_c) / 2;

    eff_APs_tau = [];

    for s = 1:p.num_setups
        if mod(s, 50) == 0
            fprintf('tau=%d: Setup %d/%d\n', p.tau_cf, s, p.num_setups);
        end

        % Generate setup (uncorrelated shadowing)
        [beta, ~, ~] = generate_setup(p, false);

        % Greedy pilot assignment + max-min DL power control
        pilot_idx = assign_pilots(beta, [], [], p, 'greedy');
        gamma = compute_gamma(beta, pilot_idx, p);
        eta = power_control_CF_DL(beta, gamma, pilot_idx, p);

        % Compute effective number of APs for each user
        for k = 1:p.K
            % Power fraction p(m,k) = eta_mk * gamma_mk / SUM_m' eta_m'k * gamma_m'k (Eq. 58)
            power_frac = eta(:, k) .* gamma(:, k);
            total_power = sum(power_frac);

            if total_power == 0
                eff_APs_tau = [eff_APs_tau; p.M]; %#ok
                continue;
            end

            power_frac = power_frac / total_power;

            % Sort in descending order
            power_sorted = sort(power_frac, 'descend');

            % Find minimum number of APs contributing >= 95%
            cum_power = cumsum(power_sorted);
            n_eff = find(cum_power >= 0.95, 1, 'first');

            if isempty(n_eff)
                n_eff = p.M;
            end

            eff_APs_tau = [eff_APs_tau; n_eff]; %#ok
        end
    end

    eff_APs{t_idx} = eff_APs_tau;
end

%% Plot Figure 9
figure('Position', [100, 100, 600, 450]);
hold on; grid on; box on;

colors = {'b', 'r'};
styles = {'-', '--'};

for t_idx = 1:length(tau_values)
    sorted = sort(eff_APs{t_idx});
    cdf_y = (1:length(sorted))' / length(sorted);
    plot(sorted, cdf_y, [colors{t_idx}, styles{t_idx}], 'LineWidth', 2);
end

xlabel('Effective Number of APs Serving Each User');
ylabel('Cumulative Distribution');
legend(arrayfun(@(t) sprintf('\\tau^{cf} = %d', t), tau_values, 'UniformOutput', false), ...
       'Location', 'southeast');
title('Figure 9: Effective Number of APs');
xlim([0 40]);

% Print statistics
for t_idx = 1:length(tau_values)
    fprintf('tau=%d: median effective APs = %.1f, mean = %.1f\n', ...
        tau_values(t_idx), median(eff_APs{t_idx}), mean(eff_APs{t_idx}));
end
