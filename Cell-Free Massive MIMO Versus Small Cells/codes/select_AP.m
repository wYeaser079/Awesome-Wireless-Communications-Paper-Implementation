function ap_selected = select_AP(beta)
%SELECT_AP Select one AP per user for the small-cell system (Eq. 38).
%   ap_selected = select_AP(beta)
%
%   Each user is served by the available AP with the largest beta_mk.
%   If an AP has already been chosen by another user, it becomes unavailable.
%   AP selection is done user by user in a random order (Section V).
%
%   Input:
%       beta        - M x K large-scale fading coefficients
%
%   Output:
%       ap_selected - K x 1 vector, ap_selected(k) = AP index serving user k

    [M, K] = size(beta);
    ap_selected = zeros(K, 1);
    available = true(M, 1);     % Track available APs

    % Process users in random order
    user_order = randperm(K);

    for idx = 1:K
        k = user_order(idx);

        % Find available AP with largest beta_mk
        beta_k = beta(:, k);
        beta_k(~available) = -inf;  % Exclude unavailable APs

        [~, best_m] = max(beta_k);

        ap_selected(k) = best_m;
        available(best_m) = false;  % Mark AP as taken
    end

end
