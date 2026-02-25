function pilot_index = assign_pilots(K, tau)
%ASSIGN_PILOTS Assign pilot sequences to K users from tau orthogonal pilots.
%   pilot_index = assign_pilots(K, tau)
%
%   When tau >= K: each user gets a unique orthogonal pilot (no contamination).
%   When tau <  K: users are randomly assigned one of tau pilots (some share).
%
%   The inner product |phi_i^H * phi_k|^2 equals:
%       1   if pilot_index(i) == pilot_index(k)   (same pilot)
%       0   otherwise                               (orthogonal)
%
%   Inputs:
%       K   - Number of users
%       tau - Number of orthogonal pilot sequences (= pilot length)
%
%   Output:
%       pilot_index - K x 1 vector, pilot_index(k) in {1, ..., tau}

    if tau >= K
        % Enough pilots: assign unique orthogonal pilot to each user
        pilot_index = (1:K).';
    else
        % Not enough pilots: random assignment causes pilot contamination
        pilot_index = randi(tau, K, 1);
    end

end
