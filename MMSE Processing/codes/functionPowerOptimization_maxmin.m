function p_opt = functionPowerOptimization_maxmin(signalTerms, interferenceTerms, ...
    noiseTerms, K, p_max, tau_c, tau_p)
%functionPowerOptimization_maxmin Solve max-min fairness power allocation
%using bisection and linear programming.
%
%   p_opt = functionPowerOptimization_maxmin(...)
%
%   Maximizes the minimum SE across all UEs subject to per-user power
%   constraints: 0 <= p_k <= p_max.
%
%   Uses bisection on the target SINR and checks feasibility at each step.
%
%   Inputs:
%       signalTerms      - K x 1 signal terms (without power)
%       interferenceTerms- K x K interference matrix (without power)
%       noiseTerms       - K x 1 noise terms
%       K                - Number of UEs
%       p_max            - Maximum per-user power
%       tau_c            - Coherence block length
%       tau_p            - Pilot length
%
%   Output:
%       p_opt - K x 1 optimal power allocation

% Bisection parameters
maxIter = 50;
tolerance = 1e-4;

% Initial bounds for the target SINR
SINR_lower = 0;
SINR_upper = 100;  % Upper bound on SINR

for iter = 1:maxIter
    SINR_target = (SINR_lower + SINR_upper) / 2;

    % Check feasibility: can we achieve SINR >= SINR_target for all UEs?
    [feasible, p_feasible] = checkFeasibility(signalTerms, interferenceTerms, ...
        noiseTerms, K, p_max, SINR_target);

    if feasible
        SINR_lower = SINR_target;
        p_opt = p_feasible;
    else
        SINR_upper = SINR_target;
    end

    if (SINR_upper - SINR_lower) < tolerance
        break;
    end
end

% If no feasible solution found, use equal power
if ~exist('p_opt', 'var')
    p_opt = p_max * ones(K, 1);
end

end


function [feasible, p_out] = checkFeasibility(signalTerms, interferenceTerms, ...
    noiseTerms, K, p_max, SINR_target)
%checkFeasibility Check if a target SINR is feasible using fixed-point iteration
%
%   For each UE k: SINR_k = p_k * signal_k / (sum_i p_i * interf_ki + noise_k)
%   We need: p_k * signal_k >= SINR_target * (sum_i p_i * interf_ki + noise_k)

maxIter = 100;
p_out = p_max * ones(K, 1);

for iter = 1:maxIter
    p_new = zeros(K, 1);

    for k = 1:K
        % Required power for UE k to achieve SINR_target
        interference_sum = 0;
        for i = 1:K
            interference_sum = interference_sum + p_out(i) * interferenceTerms(k, i);
        end

        p_required = SINR_target * (interference_sum + noiseTerms(k)) / ...
            max(signalTerms(k), 1e-20);

        p_new(k) = min(p_required, p_max);
    end

    % Check convergence
    if max(abs(p_new - p_out)) < 1e-8
        break;
    end

    p_out = p_new;
end

% Check if all UEs achieve the target SINR
feasible = true;
for k = 1:K
    interference_sum = 0;
    for i = 1:K
        interference_sum = interference_sum + p_out(i) * interferenceTerms(k, i);
    end
    achieved_SINR = p_out(k) * signalTerms(k) / max(interference_sum + noiseTerms(k), 1e-20);

    if achieved_SINR < SINR_target * 0.99  % Allow 1% tolerance
        feasible = false;
        break;
    end
end

end
