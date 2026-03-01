%MAIN_TABLE2_TABLE3 Reproduce Tables II & III from Ngo et al. (IEEE TWC, 2017).
%   95%-likely per-user net throughput for Cell-Free and Small-Cell systems
%   under six configurations:
%     1) Greedy pilot + max-min PC (uncorrelated / correlated)
%     2) Greedy pilot + no PC (uncorrelated / correlated)
%     3) Random pilot + max-min PC (uncorrelated / correlated)
%
%   Parameters: M=100, K=40, tau_cf = tau_sc = 20.

clear; clc; close all;

%% Setup
p = params();

fprintf('=== Tables II & III: 95%%-Likely Per-User Net Throughput ===\n');
fprintf('M = %d, K = %d, tau_cf = %d\n', p.M, p.K, p.tau_cf);
fprintf('Number of setups: %d\n\n', p.num_setups);

% Configuration labels
configs = {'Greedy+PC', 'Greedy+NoPC', 'Random+PC'};
shadow_types = {'uncorr', 'corr'};

% Results storage: configs x shadow x [CF_DL, SC_DL, CF_UL, SC_UL]
% For PC cases: store min(rate) per setup -> CDF over setups
% For no-PC cases: store all K rates per setup -> CDF over all

results = struct();

for c = 1:3
    for sh = 1:2
        fname = sprintf('c%d_sh%d', c, sh);
        results.(fname).cf_dl = [];
        results.(fname).cf_ul = [];
        results.(fname).sc_dl = [];
        results.(fname).sc_ul = [];
    end
end

%% Main simulation
for s = 1:p.num_setups
    if mod(s, 10) == 0
        fprintf('Setup %d/%d\n', s, p.num_setups);
    end

    for sh = 1:2
        corr = (sh == 2);
        [beta, ~, ~] = generate_setup(p, corr);

        % --- Small-Cell setup (shared) ---
        ap_sel = select_AP(beta);

        for c = 1:3
            fname = sprintf('c%d_sh%d', c, sh);

            % Determine pilot method
            if c <= 2  % Greedy
                pilot_cf = assign_pilots(beta, [], [], p, 'greedy');
                % SC greedy pilot
                pilot_sc = randi(p.tau_cf, p.K, 1);
                for n = 1:p.greedy_iterations
                    [R_tmp, ~] = compute_SE_SC_DL(beta, pilot_sc, ap_sel, p);
                    [~, k_star] = min(R_tmp);
                    mk = ap_sel(k_star);
                    best_c = inf; best_p = pilot_sc(k_star);
                    for t = 1:p.tau_cf
                        ct = 0;
                        for kp = 1:p.K
                            if kp ~= k_star && pilot_sc(kp) == t
                                ct = ct + beta(mk, kp);
                            end
                        end
                        if ct < best_c; best_c = ct; best_p = t; end
                    end
                    pilot_sc(k_star) = best_p;
                end
            else       % Random
                pilot_cf = assign_pilots(beta, [], [], p, 'random');
                pilot_sc = randi(p.tau_cf, p.K, 1);
            end

            gamma = compute_gamma(beta, pilot_cf, p);

            if c == 1 || c == 3  % With power control
                % CF power control
                eta_dl = power_control_CF_DL(beta, gamma, pilot_cf, p);
                eta_ul = power_control_CF_UL(beta, gamma, pilot_cf, p);
                [~, tp_cf_dl] = compute_SE_CF_DL(beta, gamma, pilot_cf, p, eta_dl);
                [~, tp_cf_ul] = compute_SE_CF_UL(beta, gamma, pilot_cf, p, eta_ul);

                % SC power control
                [alpha_d, alpha_u] = power_control_SC(beta, pilot_sc, ap_sel, p);
                [~, tp_sc_dl] = compute_SE_SC_DL(beta, pilot_sc, ap_sel, p, alpha_d);
                [~, tp_sc_ul] = compute_SE_SC_UL(beta, pilot_sc, ap_sel, p, alpha_u);

                % With PC: all users get same rate -> store min
                results.(fname).cf_dl = [results.(fname).cf_dl; min(tp_cf_dl)];
                results.(fname).cf_ul = [results.(fname).cf_ul; min(tp_cf_ul)];
                results.(fname).sc_dl = [results.(fname).sc_dl; min(tp_sc_dl)];
                results.(fname).sc_ul = [results.(fname).sc_ul; min(tp_sc_ul)];

            else  % No power control (c == 2)
                [~, tp_cf_dl] = compute_SE_CF_DL(beta, gamma, pilot_cf, p);
                [~, tp_cf_ul] = compute_SE_CF_UL(beta, gamma, pilot_cf, p);
                [~, tp_sc_dl] = compute_SE_SC_DL(beta, pilot_sc, ap_sel, p);
                [~, tp_sc_ul] = compute_SE_SC_UL(beta, pilot_sc, ap_sel, p);

                % No PC: each user has different rate -> store all
                results.(fname).cf_dl = [results.(fname).cf_dl; tp_cf_dl];
                results.(fname).cf_ul = [results.(fname).cf_ul; tp_cf_ul];
                results.(fname).sc_dl = [results.(fname).sc_dl; tp_sc_dl];
                results.(fname).sc_ul = [results.(fname).sc_ul; tp_sc_ul];
            end
        end
    end
end

%% Print Table II (Downlink)
to_mbps = 1e-6;

fprintf('\n========================================\n');
fprintf('TABLE II: 95%%-Likely Per-User DL Net Throughput (Mbit/s)\n');
fprintf('M=%d, K=%d, tau=%d\n', p.M, p.K, p.tau_cf);
fprintf('========================================\n');
fprintf('%-20s | %8s %8s | %8s %8s | %8s %8s\n', '', ...
    'GP+PC', '', 'GP+NoPC', '', 'RP+PC', '');
fprintf('%-20s | %8s %8s | %8s %8s | %8s %8s\n', '', ...
    'uncorr', 'corr', 'uncorr', 'corr', 'uncorr', 'corr');
fprintf('%s\n', repmat('-', 1, 76));

for sys = {'cf', 'sc'}
    if strcmp(sys{1}, 'cf'); label = 'Cell-Free'; else; label = 'Small-Cell'; end
    vals = zeros(1, 6);
    idx = 1;
    for c = [1, 2, 3]
        for sh = [1, 2]
            fname = sprintf('c%d_sh%d', c, sh);
            data = results.(fname).([sys{1} '_dl']) * to_mbps;
            vals(idx) = prctile(data, 5);
            idx = idx + 1;
        end
    end
    fprintf('%-20s | %8.2f %8.2f | %8.2f %8.2f | %8.2f %8.2f\n', label, vals);
end

%% Print Table III (Uplink)
fprintf('\n========================================\n');
fprintf('TABLE III: 95%%-Likely Per-User UL Net Throughput (Mbit/s)\n');
fprintf('M=%d, K=%d, tau=%d\n', p.M, p.K, p.tau_cf);
fprintf('========================================\n');
fprintf('%-20s | %8s %8s | %8s %8s | %8s %8s\n', '', ...
    'GP+PC', '', 'GP+NoPC', '', 'RP+PC', '');
fprintf('%-20s | %8s %8s | %8s %8s | %8s %8s\n', '', ...
    'uncorr', 'corr', 'uncorr', 'corr', 'uncorr', 'corr');
fprintf('%s\n', repmat('-', 1, 76));

for sys = {'cf', 'sc'}
    if strcmp(sys{1}, 'cf'); label = 'Cell-Free'; else; label = 'Small-Cell'; end
    vals = zeros(1, 6);
    idx = 1;
    for c = [1, 2, 3]
        for sh = [1, 2]
            fname = sprintf('c%d_sh%d', c, sh);
            data = results.(fname).([sys{1} '_ul']) * to_mbps;
            vals(idx) = prctile(data, 5);
            idx = idx + 1;
        end
    end
    fprintf('%-20s | %8.2f %8.2f | %8.2f %8.2f | %8.2f %8.2f\n', label, vals);
end
