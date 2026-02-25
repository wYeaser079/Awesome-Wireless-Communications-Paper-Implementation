% main_table2.m — Reproduce Table 2: Single-user max rate (GP)
%
% Paper: Zhang et al. (JCN, 2019), Section IV, Table 2
%
% Maximizes User 1's rate (Eq. 30/44) while all other users
% meet a QoS target of R_bar = 1 bit/s/Hz, using:
%   - Algorithm 2 (GP for ZF, Eq. 31)
%   - Algorithm 4 (GP for CB, Eq. 45)
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
R_bar      = ones(p.K, 1);         % QoS = 1 bit/s/Hz for all users

% Accumulators
rate_ZF_user1 = zeros(1, length(N_vec));
rate_ZF_other = zeros(1, length(N_vec));
rate_CB_user1 = zeros(1, length(N_vec));
rate_CB_other = zeros(1, length(N_vec));

% --- Main simulation loop ---
for in = 1:length(N_vec)
    p.N = N_vec(in);

    for s = 1:num_setups
        [beta, ~, ~] = generate_setup(p);
        pilot_index  = assign_pilots(p.K, p.tau);
        [alpha, ~]   = estimate_channel(beta, pilot_index, p);

        % GP for ZF (Algorithm 2)
        eta_zf = power_control_GP_ZF(beta, alpha, pilot_index, p, R_bar);
        R_zf   = compute_SE_ZF(beta, alpha, pilot_index, p, eta_zf);
        rate_ZF_user1(in) = rate_ZF_user1(in) + R_zf(1);
        rate_ZF_other(in) = rate_ZF_other(in) + mean(R_zf(2:end));

        % GP for CB (Algorithm 4)
        eta_cb = power_control_GP_CB(beta, alpha, pilot_index, p, R_bar);
        R_cb   = compute_SE_CB(beta, alpha, pilot_index, p, eta_cb);
        rate_CB_user1(in) = rate_CB_user1(in) + R_cb(1);
        rate_CB_other(in) = rate_CB_other(in) + mean(R_cb(2:end));
    end

    fprintf('N = %d done (%d/%d)\n', N_vec(in), in, length(N_vec));
end

% Average
rate_ZF_user1 = rate_ZF_user1 / num_setups;
rate_ZF_other = rate_ZF_other / num_setups;
rate_CB_user1 = rate_CB_user1 / num_setups;
rate_CB_other = rate_CB_other / num_setups;

% --- Display Table ---
fprintf('\n==========================================================\n');
fprintf('   Table 2: Single-User Max Rate [bits/s/Hz]  (M=%d, K=%d)\n', p.M, p.K);
fprintf('   QoS for other users: %.1f bits/s/Hz\n', R_bar(2));
fprintf('==========================================================\n');
fprintf('%-18s', '');
for in = 1:length(N_vec), fprintf('  N = %-6d', N_vec(in)); end
fprintf('\n');
fprintf('----------------------------------------------------------\n');

fprintf('%-18s', 'ZF, User 1');
for in = 1:length(N_vec), fprintf('  %-8.2f', rate_ZF_user1(in)); end
fprintf('\n');

fprintf('%-18s', 'ZF, Others (avg)');
for in = 1:length(N_vec), fprintf('  %-8.2f', rate_ZF_other(in)); end
fprintf('\n');

fprintf('%-18s', 'CB, User 1');
for in = 1:length(N_vec), fprintf('  %-8.2f', rate_CB_user1(in)); end
fprintf('\n');

fprintf('%-18s', 'CB, Others (avg)');
for in = 1:length(N_vec), fprintf('  %-8.2f', rate_CB_other(in)); end
fprintf('\n');

fprintf('----------------------------------------------------------\n');
fprintf('%-18s', 'ZF / CB ratio');
for in = 1:length(N_vec)
    fprintf('  %-8.2fx', rate_ZF_user1(in) / rate_CB_user1(in));
end
fprintf('\n');
fprintf('==========================================================\n');
