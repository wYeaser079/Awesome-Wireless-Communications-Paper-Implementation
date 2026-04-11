%% main_Fig1_convergence.m
% Reproduces Fig. 1 of Ngo-Tran-Duong-Matthaiou-Larsson, IEEE TGCN 2018:
%   "Convergence of Algorithm 1 (SCA)".
%
% It runs Algorithm 1 on a single large-scale fading realisation for
% several (M, K) pairs and plots the EE value at each SCA iteration.
% The curves rise monotonically and plateau in ~10 iterations.
%
% Requires CVX (http://cvxr.com).  Running the larger (M,K) pairs is
% slow: reduce MK_LIST to e.g. only [8,2; 20,4] for a quick smoke test.

clc; clear; close all;

%% Parameters
params = get_params();
params.D     = 1;         % 1 km square
params.tau_p = 20;
N            = 1;         % single-antenna APs (as in Fig. 1 of the paper)

% (M,K) configurations - the paper shows (8,2),(100,20),(100,40),(200,40)
MK_LIST = [ 8  2;
           20  4;
           40 10];         % kept small for runtime; expand for full repro

rng(42, 'twister');

figure('Name', 'Fig.1 Convergence of Algorithm 1', 'Color', 'w'); hold on;
styles = {'-o', '-s', '-d', '-^'};

for iCfg = 1:size(MK_LIST, 1)
    M = MK_LIST(iCfg, 1);
    K = MK_LIST(iCfg, 2);

    fprintf('\n=== Fig.1: M=%d, K=%d, N=%d ===\n', M, K, N);

    beta_mk  = generate_large_scale_fading(M, K, params);
    pilotseq = assign_pilots(K, params.tau_p);
    rho_p    = params.pilot_power_W / params.noise_power;
    gamma_mk = compute_gamma(beta_mk, pilotseq, params.tau_p, rho_p);

    % Per-user rate target for SCA: (tau_c/(tau_c-tau_p)) * S_{o,k}
    epsT = params.tau_c / (params.tau_c - params.tau_p);
    rate_QoS = epsT * params.SE_target_bps_per_Hz * ones(K, 1);

    [eta_mk, EE_trace, info] = algorithm1_EE_SCA(beta_mk, gamma_mk, ...
                                                 pilotseq, N, M, K, ...
                                                 params, rate_QoS);

    if info.feasible
        plot(1:length(EE_trace), EE_trace, styles{min(iCfg,end)}, ...
             'LineWidth', 1.5, 'DisplayName', ...
             sprintf('M=%d, K=%d', M, K));
    else
        fprintf('  infeasible for (M=%d,K=%d)\n', M, K);
    end
end

xlabel('SCA iteration index'); ylabel('Energy efficiency [Mbit/J]');
title('Fig. 1: Convergence of Algorithm 1 (N = 1, D = 1 km)');
grid on; legend('Location', 'best');
