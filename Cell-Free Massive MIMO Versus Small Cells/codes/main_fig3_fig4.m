%MAIN_FIG3_FIG4 Reproduce Figures 3 & 4 from Ngo et al. (IEEE TWC, 2017).
%   CDF of per-user DL/UL net throughput with greedy pilot assignment
%   and max-min power control, for both correlated and uncorrelated
%   shadow fading. Cell-Free vs Small-Cell comparison.
%
%   Parameters: M=100, K=40, tau_cf = tau_sc = 20.
%
%   Requires: Optimization Toolbox (linprog), CVX (optional, for DL PC)

clear; clc; close all;

%% Setup
p = params();

fprintf('=== Figures 3 & 4: Greedy Pilot + Max-Min Power Control ===\n');
fprintf('M = %d, K = %d, tau_cf = %d\n', p.M, p.K, p.tau_cf);
fprintf('Number of setups: %d\n\n', p.num_setups);

% Storage for throughputs
cf_dl_uncorr = zeros(p.num_setups, 1);
cf_ul_uncorr = zeros(p.num_setups, 1);
sc_dl_uncorr = zeros(p.num_setups, 1);
sc_ul_uncorr = zeros(p.num_setups, 1);

cf_dl_corr = zeros(p.num_setups, 1);
cf_ul_corr = zeros(p.num_setups, 1);
sc_dl_corr = zeros(p.num_setups, 1);
sc_ul_corr = zeros(p.num_setups, 1);

%% Main simulation loop
for s = 1:p.num_setups
    if mod(s, 10) == 0
        fprintf('Setup %d/%d\n', s, p.num_setups);
    end

    for corr = [false, true]
        % Generate network
        [beta, ap_pos, ue_pos] = generate_setup(p, corr);

        % --- Cell-Free Massive MIMO ---
        % Greedy pilot assignment
        pilot_idx_cf = assign_pilots(beta, [], [], p, 'greedy');
        gamma = compute_gamma(beta, pilot_idx_cf, p);

        % Max-min power control
        eta_dl = power_control_CF_DL(beta, gamma, pilot_idx_cf, p);
        eta_ul = power_control_CF_UL(beta, gamma, pilot_idx_cf, p);

        % Compute throughput
        [~, tp_dl] = compute_SE_CF_DL(beta, gamma, pilot_idx_cf, p, eta_dl);
        [~, tp_ul] = compute_SE_CF_UL(beta, gamma, pilot_idx_cf, p, eta_ul);

        % With max-min PC, all users get the same rate
        cf_dl = min(tp_dl);
        cf_ul = min(tp_ul);

        % --- Small-Cell ---
        ap_sel = select_AP(beta);

        % Greedy pilot assignment for small cells
        % (same structure but contamination measured at selected AP only)
        pilot_idx_sc = assign_pilots_sc(beta, ap_sel, p);

        % Max-min power control
        [alpha_d, alpha_u] = power_control_SC(beta, pilot_idx_sc, ap_sel, p);

        % Compute throughput
        [~, tp_sc_dl] = compute_SE_SC_DL(beta, pilot_idx_sc, ap_sel, p, alpha_d);
        [~, tp_sc_ul] = compute_SE_SC_UL(beta, pilot_idx_sc, ap_sel, p, alpha_u);

        sc_dl = min(tp_sc_dl);
        sc_ul = min(tp_sc_ul);

        % Store results
        if ~corr
            cf_dl_uncorr(s) = cf_dl;
            cf_ul_uncorr(s) = cf_ul;
            sc_dl_uncorr(s) = sc_dl;
            sc_ul_uncorr(s) = sc_ul;
        else
            cf_dl_corr(s) = cf_dl;
            cf_ul_corr(s) = cf_ul;
            sc_dl_corr(s) = sc_dl;
            sc_ul_corr(s) = sc_ul;
        end
    end
end

%% Convert to Mbit/s
to_mbps = 1e-6;

%% Figure 3: Downlink CDF
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
title('Figure 3: DL, Greedy Pilot + Max-Min PC');
xlim([0 20]);

% Print 95%-likely values
fprintf('\n--- Figure 3 (DL) 95%%-likely throughput [Mbit/s] ---\n');
fprintf('CF uncorr: %.2f, CF corr: %.2f\n', prctile(cf_dl_uncorr*to_mbps, 5), prctile(cf_dl_corr*to_mbps, 5));
fprintf('SC uncorr: %.2f, SC corr: %.2f\n', prctile(sc_dl_uncorr*to_mbps, 5), prctile(sc_dl_corr*to_mbps, 5));

%% Figure 4: Uplink CDF
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
title('Figure 4: UL, Greedy Pilot + Max-Min PC');
xlim([0 20]);

fprintf('\n--- Figure 4 (UL) 95%%-likely throughput [Mbit/s] ---\n');
fprintf('CF uncorr: %.2f, CF corr: %.2f\n', prctile(cf_ul_uncorr*to_mbps, 5), prctile(cf_ul_corr*to_mbps, 5));
fprintf('SC uncorr: %.2f, SC corr: %.2f\n', prctile(sc_ul_uncorr*to_mbps, 5), prctile(sc_ul_corr*to_mbps, 5));


%% --- Helper: CDF plot ---
function plot_cdf(data, linestyle, linewidth)
    sorted = sort(data);
    cdf_y = (1:length(sorted))' / length(sorted);
    plot(sorted, cdf_y, linestyle, 'LineWidth', linewidth);
end


%% --- Helper: Greedy pilot assignment for small cells ---
function pilot_idx = assign_pilots_sc(beta, ap_sel, p)
%ASSIGN_PILOTS_SC Greedy pilot assignment for small-cell systems.
%   Same as Algorithm 1, but contamination is measured at the selected AP
%   (not summed over all APs).

    K = size(beta, 2);
    pilot_idx = randi(p.tau_cf, K, 1);

    for n = 1:p.greedy_iterations
        % Compute SC DL rates with current assignment and equal power
        [R_dl, ~] = compute_SE_SC_DL(beta, pilot_idx, ap_sel, p);

        [~, k_star] = min(R_dl);
        mk = ap_sel(k_star);

        best_contam = inf;
        best_pilot = pilot_idx(k_star);

        for t = 1:p.tau_cf
            contam = 0;
            for kp = 1:K
                if kp ~= k_star && pilot_idx(kp) == t
                    contam = contam + beta(mk, kp);
                end
            end
            if contam < best_contam
                best_contam = contam;
                best_pilot = t;
            end
        end

        pilot_idx(k_star) = best_pilot;
    end
end
