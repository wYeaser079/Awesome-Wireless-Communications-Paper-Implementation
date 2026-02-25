% main_fig5.m — Reproduce Figure 5: SCA convergence
%
% Paper: Zhang et al. (JCN, 2019), Section IV, Figure 5
%
% Shows the sum rate increasing monotonically with each SCA iteration
% for both ZF (Algorithm 1, Eq. 21-29) and CB (Algorithm 3, Eq. 36-43).
% Iteration 0 corresponds to equal power control (eta_k = 1).
%
% Requires: CVX toolbox (http://cvxr.com/cvx/)
% Fixed: M = 60, K = 10, N = 6.

clear; close all; clc;

% --- Check CVX ---
if ~exist('cvx_begin', 'file')
    error('This script requires CVX. Install from http://cvxr.com/cvx/');
end

% --- Configuration ---
p = params();

num_setups = 50;
max_iter   = 10;

% Pre-allocate: (max_iter+1) x num_setups
%   Row 1 = iteration 0 (EPC), Row 2 = iteration 1, ...
history_ZF = zeros(max_iter + 1, num_setups);
history_CB = zeros(max_iter + 1, num_setups);

% --- Main simulation loop ---
for s = 1:num_setups
    [beta, ~, ~] = generate_setup(p);
    pilot_index  = assign_pilots(p.K, p.tau);
    [alpha, ~]   = estimate_channel(beta, pilot_index, p);

    % ZF SCA with history (Algorithm 1)
    h_zf = sca_track_ZF(beta, alpha, pilot_index, p, max_iter);
    history_ZF(:, s) = h_zf;

    % CB SCA with history (Algorithm 3)
    h_cb = sca_track_CB(beta, alpha, pilot_index, p, max_iter);
    history_CB(:, s) = h_cb;

    fprintf('Setup %d/%d done\n', s, num_setups);
end

% Average over setups
avg_ZF = mean(history_ZF, 2);
avg_CB = mean(history_CB, 2);

% --- Plot ---
figure; hold on; grid on; box on;

plot(0:max_iter, avg_ZF, '-o', 'Color', [0 0.45 0.74], ...
     'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', [0 0.45 0.74]);
plot(0:max_iter, avg_CB, '-s', 'Color', [0.85 0.33 0.10], ...
     'LineWidth', 1.5, 'MarkerSize', 7, 'MarkerFaceColor', [0.85 0.33 0.10]);

xlabel('SCA iteration');
ylabel('Sum spectral efficiency [bits/s/Hz]');
title('Figure 5: SCA Convergence');
legend('ZF (Algorithm 1)', 'CB (Algorithm 3)', 'Location', 'southeast');

% --- Summary ---
fprintf('\n--- Convergence Summary ---\n');
fprintf('ZF: EPC = %.2f → SCA = %.2f  (gain = %.1f%%)\n', ...
    avg_ZF(1), avg_ZF(end), 100*(avg_ZF(end)/avg_ZF(1) - 1));
fprintf('CB: EPC = %.2f → SCA = %.2f  (gain = %.1f%%)\n', ...
    avg_CB(1), avg_CB(end), 100*(avg_CB(end)/avg_CB(1) - 1));


% =====================================================================
%  Local functions (MATLAB R2016b+)
% =====================================================================

function history = sca_track_ZF(beta, alpha, pilot_index, p, max_iter)
%SCA_TRACK_ZF  Run ZF SCA and record sum rate at each iteration.
    [M, K] = size(beta);

    % Precompute ZF parameters (same as power_control_SCA_ZF)
    c_err = max(beta - alpha, [], 1).';
    [~, dom] = max(beta, [], 1);
    a = zeros(K, 1);
    for k = 1:K
        excl = dom([1:k-1, k+1:K]);
        Mk = setdiff(1:M, excl);
        aMk = alpha(Mk, k);
        S1 = sum(aMk); S2 = sum(aMk.^2);
        a(k) = p.rho_u * (p.N * S1^2 - S2) / S1;
    end

    eta = ones(K, 1);
    history = zeros(max_iter + 1, 1);
    history(1) = sum(compute_SE_ZF(beta, alpha, pilot_index, p, eta));

    for iter = 1:max_iter
        D_n = p.rho_u * (c_err' * eta) + 1;
        d   = K * p.rho_u * c_err / D_n;

        cvx_begin quiet
            variable eta_new(K, 1)
            obj = -d' * eta_new;
            for k = 1:K
                obj = obj + log(a(k)*eta_new(k) ...
                            + p.rho_u*(c_err'*eta_new) + 1);
            end
            maximize(obj)
            subject to
                eta_new >= 0;
                eta_new <= 1;
        cvx_end

        if isempty(strfind(cvx_status, 'Solved'))
            history(iter+1:end) = history(iter);
            break;
        end
        eta = eta_new;
        history(iter+1) = sum(compute_SE_ZF(beta, alpha, pilot_index, p, eta));
    end
end


function history = sca_track_CB(beta, alpha, pilot_index, p, max_iter)
%SCA_TRACK_CB  Run CB SCA and record sum rate at each iteration.
    [M, K] = size(beta);

    % Precompute CB parameters (same as power_control_SCA_CB)
    S = sum(alpha.^2, 1).';
    w = zeros(K, K);
    for k = 1:K
        for i = 1:K
            noncoh = sum(beta(:,i) .* alpha(:,k)) / p.N;
            coh = 0;
            if i ~= k && pilot_index(i) == pilot_index(k)
                coh = sum((alpha(:,k) .* beta(:,i) ./ beta(:,k)).^2);
            end
            w(i,k) = coh + noncoh;
        end
    end
    noise = sum(alpha, 1).' / (p.N * p.rho_u);

    eta = ones(K, 1);
    history = zeros(max_iter + 1, 1);
    history(1) = sum(compute_SE_CB(beta, alpha, pilot_index, p, eta));

    for iter = 1:max_iter
        D_n = w' * eta + noise;

        cvx_begin quiet
            variable eta_new(K, 1)
            obj = 0;
            for k = 1:K
                obj = obj + log(eta_new(k)*S(k) + w(:,k)'*eta_new + noise(k));
                obj = obj - w(:,k)'*eta_new / D_n(k);
            end
            maximize(obj)
            subject to
                eta_new >= 0;
                eta_new <= 1;
        cvx_end

        if isempty(strfind(cvx_status, 'Solved'))
            history(iter+1:end) = history(iter);
            break;
        end
        eta = eta_new;
        history(iter+1) = sum(compute_SE_CB(beta, alpha, pilot_index, p, eta));
    end
end
