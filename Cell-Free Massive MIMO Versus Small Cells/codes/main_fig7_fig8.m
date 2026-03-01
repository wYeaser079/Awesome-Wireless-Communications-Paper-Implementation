%MAIN_FIG7_FIG8 Reproduce Figures 7 & 8 from Ngo et al. (IEEE TWC, 2017).
%   CDF of per-user DL/UL net throughput with RANDOM pilot assignment
%   and max-min power control. Cell-Free vs Small-Cell comparison.
%
%   Parameters: M=100, K=40, tau_cf = tau_sc = 20.

clear; clc; close all;

%% Setup
p = params();

fprintf('=== Figures 7 & 8: Random Pilot + Max-Min Power Control ===\n');
fprintf('M = %d, K = %d, tau_cf = %d\n', p.M, p.K, p.tau_cf);
fprintf('Number of setups: %d\n\n', p.num_setups);

% Storage (with max-min PC, all users get same rate per realization)
cf_dl_uncorr = zeros(p.num_setups, 1);  cf_dl_corr = zeros(p.num_setups, 1);
cf_ul_uncorr = zeros(p.num_setups, 1);  cf_ul_corr = zeros(p.num_setups, 1);
sc_dl_uncorr = zeros(p.num_setups, 1);  sc_dl_corr = zeros(p.num_setups, 1);
sc_ul_uncorr = zeros(p.num_setups, 1);  sc_ul_corr = zeros(p.num_setups, 1);

%% Main simulation loop
for s = 1:p.num_setups
    if mod(s, 10) == 0
        fprintf('Setup %d/%d\n', s, p.num_setups);
    end

    for corr = [false, true]
        [beta, ~, ~] = generate_setup(p, corr);

        % --- Cell-Free with RANDOM pilot assignment ---
        pilot_idx_cf = assign_pilots(beta, [], [], p, 'random');
        gamma = compute_gamma(beta, pilot_idx_cf, p);

        % Max-min power control
        eta_dl = power_control_CF_DL(beta, gamma, pilot_idx_cf, p);
        eta_ul = power_control_CF_UL(beta, gamma, pilot_idx_cf, p);

        [~, tp_dl] = compute_SE_CF_DL(beta, gamma, pilot_idx_cf, p, eta_dl);
        [~, tp_ul] = compute_SE_CF_UL(beta, gamma, pilot_idx_cf, p, eta_ul);

        cf_dl = min(tp_dl);
        cf_ul = min(tp_ul);

        % --- Small-Cell with RANDOM pilot assignment ---
        ap_sel = select_AP(beta);
        pilot_idx_sc = randi(p.tau_cf, p.K, 1);  % Random pilots

        [alpha_d, alpha_u] = power_control_SC(beta, pilot_idx_sc, ap_sel, p);

        [~, tp_sc_dl] = compute_SE_SC_DL(beta, pilot_idx_sc, ap_sel, p, alpha_d);
        [~, tp_sc_ul] = compute_SE_SC_UL(beta, pilot_idx_sc, ap_sel, p, alpha_u);

        sc_dl = min(tp_sc_dl);
        sc_ul = min(tp_sc_ul);

        if ~corr
            cf_dl_uncorr(s) = cf_dl;  cf_ul_uncorr(s) = cf_ul;
            sc_dl_uncorr(s) = sc_dl;  sc_ul_uncorr(s) = sc_ul;
        else
            cf_dl_corr(s) = cf_dl;  cf_ul_corr(s) = cf_ul;
            sc_dl_corr(s) = sc_dl;  sc_ul_corr(s) = sc_ul;
        end
    end
end

%% Convert to Mbit/s and plot
to_mbps = 1e-6;

%% Figure 7: Downlink CDF
figure('Position', [100, 100, 600, 450]);
hold on; grid on; box on;

plot_cdf(cf_dl_uncorr * to_mbps, 'b-', 2);
plot_cdf(cf_dl_corr * to_mbps, 'b--', 2);
plot_cdf(sc_dl_uncorr * to_mbps, 'r-', 2);
plot_cdf(sc_dl_corr * to_mbps, 'r--', 2);

xlabel('Per-User Downlink Net Throughput (Mbits/s)');
ylabel('Cumulative Distribution');
legend('Cell-Free, uncorrelated', 'Cell-Free, correlated', ...
       'Small-Cell, uncorrelated', 'Small-Cell, correlated', ...
       'Location', 'southeast');
title('Figure 7: DL, Random Pilot + Max-Min PC');
xlim([0 20]);

fprintf('\n--- Figure 7 (DL) 95%%-likely throughput [Mbit/s] ---\n');
fprintf('CF uncorr: %.2f, CF corr: %.2f\n', prctile(cf_dl_uncorr*to_mbps, 5), prctile(cf_dl_corr*to_mbps, 5));
fprintf('SC uncorr: %.2f, SC corr: %.2f\n', prctile(sc_dl_uncorr*to_mbps, 5), prctile(sc_dl_corr*to_mbps, 5));

%% Figure 8: Uplink CDF
figure('Position', [750, 100, 600, 450]);
hold on; grid on; box on;

plot_cdf(cf_ul_uncorr * to_mbps, 'b-', 2);
plot_cdf(cf_ul_corr * to_mbps, 'b--', 2);
plot_cdf(sc_ul_uncorr * to_mbps, 'r-', 2);
plot_cdf(sc_ul_corr * to_mbps, 'r--', 2);

xlabel('Per-User Uplink Net Throughput (Mbits/s)');
ylabel('Cumulative Distribution');
legend('Cell-Free, uncorrelated', 'Cell-Free, correlated', ...
       'Small-Cell, uncorrelated', 'Small-Cell, correlated', ...
       'Location', 'southeast');
title('Figure 8: UL, Random Pilot + Max-Min PC');
xlim([0 20]);

fprintf('\n--- Figure 8 (UL) 95%%-likely throughput [Mbit/s] ---\n');
fprintf('CF uncorr: %.2f, CF corr: %.2f\n', prctile(cf_ul_uncorr*to_mbps, 5), prctile(cf_ul_corr*to_mbps, 5));
fprintf('SC uncorr: %.2f, SC corr: %.2f\n', prctile(sc_ul_uncorr*to_mbps, 5), prctile(sc_ul_corr*to_mbps, 5));


function plot_cdf(data, linestyle, linewidth)
    sorted = sort(data);
    cdf_y = (1:length(sorted))' / length(sorted);
    plot(sorted, cdf_y, linestyle, 'LineWidth', linewidth);
end
