%MAIN_FIG5_FIG6 Reproduce Figures 5 & 6 from Ngo et al. (IEEE TWC, 2017).
%   CDF of per-user DL/UL net throughput with greedy pilot assignment
%   and WITHOUT power control. Cell-Free vs Small-Cell comparison.
%
%   Parameters: M=100, K=40, tau_cf = tau_sc = 20.

clear; clc; close all;

%% Setup
p = params();

fprintf('=== Figures 5 & 6: Greedy Pilot + NO Power Control ===\n');
fprintf('M = %d, K = %d, tau_cf = %d\n', p.M, p.K, p.tau_cf);
fprintf('Number of setups: %d\n\n', p.num_setups);

% Storage: without PC, each user has different rate -> store all K rates
cf_dl_uncorr = [];  cf_dl_corr = [];
cf_ul_uncorr = [];  cf_ul_corr = [];
sc_dl_uncorr = [];  sc_dl_corr = [];
sc_ul_uncorr = [];  sc_ul_corr = [];

%% Main simulation loop
for s = 1:p.num_setups
    if mod(s, 10) == 0
        fprintf('Setup %d/%d\n', s, p.num_setups);
    end

    for corr = [false, true]
        [beta, ~, ~] = generate_setup(p, corr);

        % --- Cell-Free (no power control) ---
        pilot_idx_cf = assign_pilots(beta, [], [], p, 'greedy');
        gamma = compute_gamma(beta, pilot_idx_cf, p);

        % Equal power: eta_mk = 1/SUM_k' gamma_mk' (DL), eta_k=1 (UL)
        [~, tp_dl] = compute_SE_CF_DL(beta, gamma, pilot_idx_cf, p);
        [~, tp_ul] = compute_SE_CF_UL(beta, gamma, pilot_idx_cf, p);

        % --- Small-Cell (no power control) ---
        ap_sel = select_AP(beta);
        pilot_idx_sc = randi(p.tau_cf, p.K, 1);  % Random pilots for SC

        % Greedy pilot for SC
        for n = 1:p.greedy_iterations
            [R_tmp, ~] = compute_SE_SC_DL(beta, pilot_idx_sc, ap_sel, p);
            [~, k_star] = min(R_tmp);
            mk = ap_sel(k_star);
            best_contam = inf; best_pilot = pilot_idx_sc(k_star);
            for t = 1:p.tau_cf
                contam = 0;
                for kp = 1:p.K
                    if kp ~= k_star && pilot_idx_sc(kp) == t
                        contam = contam + beta(mk, kp);
                    end
                end
                if contam < best_contam
                    best_contam = contam; best_pilot = t;
                end
            end
            pilot_idx_sc(k_star) = best_pilot;
        end

        % Full power (alpha=1)
        [~, tp_sc_dl] = compute_SE_SC_DL(beta, pilot_idx_sc, ap_sel, p);
        [~, tp_sc_ul] = compute_SE_SC_UL(beta, pilot_idx_sc, ap_sel, p);

        % Store all user throughputs
        if ~corr
            cf_dl_uncorr = [cf_dl_uncorr; tp_dl]; %#ok
            cf_ul_uncorr = [cf_ul_uncorr; tp_ul]; %#ok
            sc_dl_uncorr = [sc_dl_uncorr; tp_sc_dl]; %#ok
            sc_ul_uncorr = [sc_ul_uncorr; tp_sc_ul]; %#ok
        else
            cf_dl_corr = [cf_dl_corr; tp_dl]; %#ok
            cf_ul_corr = [cf_ul_corr; tp_ul]; %#ok
            sc_dl_corr = [sc_dl_corr; tp_sc_dl]; %#ok
            sc_ul_corr = [sc_ul_corr; tp_sc_ul]; %#ok
        end
    end
end

%% Convert to Mbit/s
to_mbps = 1e-6;

%% Figure 5: Downlink CDF (no PC)
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
title('Figure 5: DL, Greedy Pilot + No Power Control');
xlim([0 20]);

fprintf('\n--- Figure 5 (DL) 95%%-likely throughput [Mbit/s] ---\n');
fprintf('CF uncorr: %.2f, CF corr: %.2f\n', prctile(cf_dl_uncorr*to_mbps, 5), prctile(cf_dl_corr*to_mbps, 5));
fprintf('SC uncorr: %.2f, SC corr: %.2f\n', prctile(sc_dl_uncorr*to_mbps, 5), prctile(sc_dl_corr*to_mbps, 5));

%% Figure 6: Uplink CDF (no PC)
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
title('Figure 6: UL, Greedy Pilot + No Power Control');
xlim([0 20]);

fprintf('\n--- Figure 6 (UL) 95%%-likely throughput [Mbit/s] ---\n');
fprintf('CF uncorr: %.2f, CF corr: %.2f\n', prctile(cf_ul_uncorr*to_mbps, 5), prctile(cf_ul_corr*to_mbps, 5));
fprintf('SC uncorr: %.2f, SC corr: %.2f\n', prctile(sc_ul_uncorr*to_mbps, 5), prctile(sc_ul_corr*to_mbps, 5));


function plot_cdf(data, linestyle, linewidth)
    sorted = sort(data);
    cdf_y = (1:length(sorted))' / length(sorted);
    plot(sorted, cdf_y, linestyle, 'LineWidth', linewidth);
end
