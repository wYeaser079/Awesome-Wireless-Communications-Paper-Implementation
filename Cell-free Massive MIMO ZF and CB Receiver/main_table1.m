% main_table1.m — Reproduce Table 1: Sum rate comparison (EPC vs SCA)
%
% Paper: Zhang et al. (JCN, 2019), Section IV, Table 1
%
% Compares the sum SE [bits/s/Hz] for:
%   - ZF with Equal Power Control (EPC) vs SCA (Algorithm 1)
%   - CB with Equal Power Control (EPC) vs SCA (Algorithm 3)
% across N = 2, 4, 6 antennas per AP.
%
% Requires: CVX toolbox (http://cvxr.com/cvx/)
% Fixed: M = 60, K = 10.

clear; clc;

% --- Check CVX ---
if ~exist('cvx_begin', 'file')
    error('This script requires CVX. Install from http://cvxr.com/cvx/');
end

% --- Configuration ---
p = params();

N_vec      = [2, 4, 6];
num_setups = p.num_setups;

% Accumulators
sum_ZF_epc = zeros(1, length(N_vec));
sum_ZF_sca = zeros(1, length(N_vec));
sum_CB_epc = zeros(1, length(N_vec));
sum_CB_sca = zeros(1, length(N_vec));

% --- Main simulation loop ---
for in = 1:length(N_vec)
    p.N = N_vec(in);

    for s = 1:num_setups
        [beta, ~, ~] = generate_setup(p);
        pilot_index  = assign_pilots(p.K, p.tau);
        [alpha, ~]   = estimate_channel(beta, pilot_index, p);

        % Equal power control
        sum_ZF_epc(in) = sum_ZF_epc(in) + sum(compute_SE_ZF(beta, alpha, pilot_index, p));
        sum_CB_epc(in) = sum_CB_epc(in) + sum(compute_SE_CB(beta, alpha, pilot_index, p));

        % SCA power control
        eta_zf = power_control_SCA_ZF(beta, alpha, pilot_index, p);
        eta_cb = power_control_SCA_CB(beta, alpha, pilot_index, p);

        sum_ZF_sca(in) = sum_ZF_sca(in) + sum(compute_SE_ZF(beta, alpha, pilot_index, p, eta_zf));
        sum_CB_sca(in) = sum_CB_sca(in) + sum(compute_SE_CB(beta, alpha, pilot_index, p, eta_cb));
    end

    fprintf('N = %d done (%d/%d)\n', N_vec(in), in, length(N_vec));
end

% Average
sum_ZF_epc = sum_ZF_epc / num_setups;
sum_ZF_sca = sum_ZF_sca / num_setups;
sum_CB_epc = sum_CB_epc / num_setups;
sum_CB_sca = sum_CB_sca / num_setups;

% --- Display Table ---
fprintf('\n============================================\n');
fprintf('   Table 1: Sum SE [bits/s/Hz]  (M=%d, K=%d)\n', p.M, p.K);
fprintf('============================================\n');
fprintf('%-12s', '');
for in = 1:length(N_vec), fprintf('  N = %-6d', N_vec(in)); end
fprintf('\n');
fprintf('--------------------------------------------\n');

fprintf('%-12s', 'ZF, EPC');
for in = 1:length(N_vec), fprintf('  %-8.2f', sum_ZF_epc(in)); end
fprintf('\n');

fprintf('%-12s', 'ZF, SCA');
for in = 1:length(N_vec), fprintf('  %-8.2f', sum_ZF_sca(in)); end
fprintf('\n');

fprintf('%-12s', 'CB, EPC');
for in = 1:length(N_vec), fprintf('  %-8.2f', sum_CB_epc(in)); end
fprintf('\n');

fprintf('%-12s', 'CB, SCA');
for in = 1:length(N_vec), fprintf('  %-8.2f', sum_CB_sca(in)); end
fprintf('\n');

fprintf('--------------------------------------------\n');
fprintf('%-12s', 'ZF gain');
for in = 1:length(N_vec)
    fprintf('  %-7.1f%%', 100*(sum_ZF_sca(in)/sum_ZF_epc(in) - 1));
end
fprintf('\n');

fprintf('%-12s', 'CB gain');
for in = 1:length(N_vec)
    fprintf('  %-7.1f%%', 100*(sum_CB_sca(in)/sum_CB_epc(in) - 1));
end
fprintf('\n');
fprintf('============================================\n');
