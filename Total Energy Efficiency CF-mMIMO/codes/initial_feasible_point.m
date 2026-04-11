function [c0, u0, t0, feasible] = initial_feasible_point(M, N, K, ...
                                                    gamma_mk, beta_mk, pilotseq, ...
                                                    rho_d, rate_QoS)
%INITIAL_FEASIBLE_POINT  Solve the feasibility SOCP that produces a
%  starting point (c0, u0, t0) for Algorithm 1 of Ngo et al. 2018.
%  Mirrors GENERATEINITIALPOINT.M in the reference repo at
%  github.com/tranlenam/cellfreeMIMOenergyefficiency.
%
%  The returned c0 is the square-root of the power control coefficients
%  (c0 = sqrt(eta)), u0 = 2.^t0, and t0 = rate_QoS.
%
%  Requires CVX (http://cvxr.com).

    cvx_quiet true
    cvx_begin
        variable c(M, K)
        c >= 0;                                                         %#ok<VUNUS>
        % (25d) Per-AP power:  sum_k gamma_mk c_mk^2 <= 1/N
        for m = 1:M
            norm( sqrt(gamma_mk(m, :))' .* c(m, :)' ) <= 1/sqrt(N);
        end
        % (27) Per-user SE: SE_k >= rate_QoS_k in SOC form
        for k = 1:K
            interf = interference_vec(M, N, K, c, ...
                                      sqrt(rho_d) * gamma_mk, ...
                                      sqrt(rho_d) * beta_mk, ...
                                      pilotseq, k);
            lhs = (1 / sqrt(2^rate_QoS(k) - 1)) ...
                  * sqrt(rho_d) * (gamma_mk(:, k)' * c(:, k));
            norm( [interf; 1/N] ) <= lhs;
        end
    cvx_end

    if strcmpi(cvx_status, 'Solved') || strcmpi(cvx_status, 'Inaccurate/Solved')
        c0 = c;
        t0 = rate_QoS;
        u0 = 2.^t0;
        feasible = true;
    else
        c0 = NaN; u0 = NaN; t0 = NaN; feasible = false;
    end
end
