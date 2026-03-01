function pilot_index = assign_pilots(beta, gamma, pilot_index_init, p, method)
%ASSIGN_PILOTS Assign pilot sequences to users (random or greedy).
%   pilot_index = assign_pilots(beta, gamma, pilot_index_init, p, method)
%
%   Two methods:
%     'random' - Each user randomly picks from tau_cf orthogonal pilots
%     'greedy' - Algorithm 1: iteratively reassign worst user's pilot
%
%   For the greedy algorithm (Algorithm 1, Section IV-A):
%     1) Start with random pilot assignment
%     2) Compute DL rates, find worst user k*
%     3) Reassign k*'s pilot to minimize pilot contamination (Eq. 30)
%     4) Repeat for N iterations
%
%   Inputs:
%       beta             - M x K large-scale fading coefficients
%       gamma            - M x K channel estimate variances (used by greedy)
%                          Can be [] for 'random' method
%       pilot_index_init - K x 1 initial pilot assignment ([] for fresh random)
%       p                - Parameter struct from params()
%       method           - 'random' or 'greedy'
%
%   Output:
%       pilot_index      - K x 1 pilot assignment vector (values in 1:tau_cf)

    K = size(beta, 2);

    % --- Random pilot assignment ---
    if isempty(pilot_index_init)
        pilot_index = randi(p.tau_cf, K, 1);
    else
        pilot_index = pilot_index_init;
    end

    if strcmpi(method, 'random')
        return;
    end

    % --- Greedy pilot assignment (Algorithm 1) ---
    for n = 1:p.greedy_iterations

        % Recompute gamma with current pilot assignment
        gamma_cur = compute_gamma(beta, pilot_index, p);

        % Compute DL rates with equal power allocation (no PC)
        [R_dl, ~] = compute_SE_CF_DL(beta, gamma_cur, pilot_index, p);

        % Step 2: Find the user with the lowest rate (Eq. 31)
        [~, k_star] = min(R_dl);

        % Step 3: Try all tau_cf pilots for user k*, pick the best
        best_contamination = inf;
        best_pilot = pilot_index(k_star);

        for t = 1:p.tau_cf
            % Compute pilot contamination if k* uses pilot t (Eq. 30):
            %   SUM_m SUM_{k' != k*, pilot_index(k')=t} beta_mk'
            contamination = 0;
            for kp = 1:K
                if kp == k_star
                    continue;
                end
                if pilot_index(kp) == t
                    contamination = contamination + sum(beta(:, kp));
                end
            end

            if contamination < best_contamination
                best_contamination = contamination;
                best_pilot = t;
            end
        end

        % Update pilot assignment for k*
        pilot_index(k_star) = best_pilot;
    end

end
