function [eta_mk, EE_trace, info] = algorithm1_EE_SCA(beta_mk, gamma_mk, ...
                                                       pilotseq, N, M, K, params, ...
                                                       rate_QoS)
%ALGORITHM1_EE_SCA  Sequential Convex Approximation (SCA) solver for
%  the energy-efficiency maximisation problem (P) / (P^_1) of
%  Ngo-Tran-Duong-Matthaiou-Larsson IEEE TGCN 2018.
%
%  Each iteration solves a Second-Order Cone Program (SOCP) built from
%  the current tangent point (cn, un).  See Eqs. (24)-(37) of the paper.
%  At convergence, the returned ETA_MK is the optimal EE power
%  allocation.
%
%  INPUTS:
%    beta_mk, gamma_mk - M-by-K large-scale fading and MMSE variance
%    pilotseq          - tau_p-by-K pilot matrix
%    N, M, K           - antennas per AP, number of APs, number of users
%    params            - struct from get_params.m
%    rate_QoS          - K-by-1 per-user SE target.  Passed already as
%                        (tau_c/(tau_c-tau_p)) * S_{o,k} so that the
%                        constraint in the SOCP is log2(1+SINR) >= t.
%                        (At the end EE uses the un-scaled S_e.)
%
%  OUTPUTS:
%    eta_mk   - M-by-K optimal power control matrix (actual etas, not sqrt)
%    EE_trace - vector of EE values over iterations [Mbit/J]
%    info     - struct with 'feasible', 'converged', 'num_iter',
%               'SE_k' (per-user SE at the optimum), 'c0'

    % ---- Unpack ------------------------------------------------------
    P_per_antenna = params.P_per_antenna_W;
    N0            = params.noise_power;
    eta_PA        = params.eta_PA;
    P_bt_per_Mbps = params.P_bt_WperGbps * 1e-3;            % W/(Mbit/s)
    P_tc          = params.P_tc_W;
    P_0           = params.P_0_W;
    B_MHz         = params.B_MHz;
    tau_c         = params.tau_c;
    tau_p         = params.tau_p;
    max_iter      = params.sca_max_iter;
    tol           = params.sca_tol;

    rho_d = P_per_antenna * N / N0;        % normalised DL SNR per AP

    alpha_m = (1/eta_PA) * ones(M, 1);
    P_fix_bar = N * sum(P_tc * ones(M, 1)) + sum(P_0 * ones(M, 1));

    info.feasible  = false;
    info.converged = false;
    info.num_iter  = 0;
    info.SE_k      = zeros(K, 1);

    EE_trace = [];
    eta_mk = [];

    % ---- Initialisation: feasibility SOCP ---------------------------
    [cn, un, ~, feasible] = initial_feasible_point(M, N, K, ...
                                                   gamma_mk, beta_mk, pilotseq, ...
                                                   rho_d, rate_QoS);
    if ~feasible
        info.feasible = false;
        return;
    end
    info.feasible = true;
    info.c0       = cn;

    % ---- SCA loop ----------------------------------------------------
    sum_tdot_prev = 0;
    converged = false;
    for iter = 1:max_iter

        cvx_quiet true
        cvx_begin
            variable cdot(M, K)
            variable tdot(K, 1)
            variable udot(K, 1)
            variable theta

            % Objective: maximise sum(tdot)  (equiv. to B*(1-tp/tc)*Se/Aux denom)
            maximise( sum(tdot) )

            subject to
                cdot(:) >= 0;                                           %#ok<VUNUS>
                tdot >= 0;                                              %#ok<VUNUS>
                theta >= 1e-9;                                          %#ok<VUNUS>

                % (36b) Per-AP power:
                %    sum_k gamma_mk cdot_mk^2 <= theta^2/N
                for m = 1:M
                    norm( sqrt(gamma_mk(m, :))' .* cdot(m, :)' ) <= theta/sqrt(N);
                end

                % (36e) Denominator slack (rotated cone form):
                %    P_fix_bar*theta^2 + rho_d*N0*N*sum(alpha_m gamma_mk cdot_mk^2) <= theta
                % Encoded as:  || [sqrt(P_fix_bar)*theta; sqrt(rho_d*N0*N*alpha.*gamma)(:).*cdot(:); (theta-1)/2] ||
                %                                                   <= (theta+1)/2
                Gammaan = sqrt( rho_d * N0 * N * ...
                                repmat(alpha_m, 1, K) .* gamma_mk );    % M x K
                norm( [ sqrt(P_fix_bar) * theta; ...
                        reshape(Gammaan .* cdot, M*K, 1); ...
                        0.5 * (theta - 1) ] ) <= 0.5 * (theta + 1);

                % (36f) Log constraint (SOC form of the hyperbolic  X*Y >= un*theta^2)
                %      X = udot + theta(log(un)+1) - log(2)*tdot
                %      Y = udot - theta(log(un)+1) + log(2)*tdot
                %      X >= || [Y ; 2*sqrt(un)*theta] ||
                logun = log(un);
                for k = 1:K
                    X = udot(k) + theta*(logun(k) + 1) - log(2) * tdot(k);
                    Y = udot(k) - theta*(logun(k) + 1) + log(2) * tdot(k);
                    Z = 2 * sqrt(un(k)) * theta;
                    norm([Y; Z]) <= X;
                end

                % (36d) Per-user SE constraint in SOC form
                % (27)   sqrt(rho_d)(gamma_kk^T cdot_k) / sqrt(2^Rk - 1)
                %        >= || [interf_vec; theta/N] ||
                for k = 1:K
                    interf = interference_vec(M, N, K, cdot, ...
                                              sqrt(rho_d) * gamma_mk, ...
                                              sqrt(rho_d) * beta_mk, ...
                                              pilotseq, k);
                    lhs = (1 / sqrt(2^rate_QoS(k) - 1)) ...
                           * sqrt(rho_d) * (gamma_mk(:, k)' * cdot(:, k));
                    norm( [ interf; theta/N ] ) <= lhs;
                end

                % (36g) SCA convex approximation of SE_k >= tdot_k:
                %    || [ 2*sqrt(rho_d)*N*interf(cdot); 2*theta; theta - approx ] ||
                %        <= theta + approx
                for k = 1:K
                    interf_full = interference_vec(M, N, K, cdot, ...
                                                    gamma_mk, beta_mk, ...
                                                    pilotseq, k);
                    approx = approx_function(M, N, K, gamma_mk, beta_mk, ...
                                              pilotseq, rho_d, cdot, udot, ...
                                              cn, un, theta, k);
                    norm( [ 2*sqrt(rho_d)*N*interf_full; 2*theta; ...
                            theta - approx ] ) <= theta + approx;
                end
        cvx_end

        % --- Check solver status ------------------------------------
        if ~(strcmpi(cvx_status, 'Solved') || strcmpi(cvx_status, 'Inaccurate/Solved'))
            warning('algorithm1_EE_SCA:SolverFailed', ...
                    'CVX status: %s at iter %d', cvx_status, iter);
            break;
        end

        % --- Update tangent point -----------------------------------
        theta_val = max(double(theta), 1e-12);
        cn = double(cdot) / theta_val;
        un = double(udot) / theta_val;
        tn = double(tdot) / theta_val;    %#ok<NASGU>

        % --- Compute EE from this iterate ---------------------------
        eta_current = max(cn.^2, 0);     % since cn = sqrt(eta)
        SE_k = compute_SE_closedform(eta_current, beta_mk, gamma_mk, ...
                                     pilotseq, N, rho_d, tau_c, tau_p);
        [EE_iter, ~, ~, ~] = compute_EE(eta_current, gamma_mk, SE_k, ...
                                        params, N, M);
        EE_trace(end+1) = EE_iter;   %#ok<AGROW>

        % --- Convergence check --------------------------------------
        sum_tdot = sum(double(tdot));
        if iter > 1 && abs(sum_tdot - sum_tdot_prev)/max(abs(sum_tdot_prev),1e-9) < tol
            converged = true;
            info.num_iter = iter;
            break;
        end
        sum_tdot_prev = sum_tdot;
        info.num_iter = iter;
    end

    % ---- Return final power allocation ------------------------------
    eta_mk = max(cn.^2, 0);
    info.converged = converged;
    info.SE_k = compute_SE_closedform(eta_mk, beta_mk, gamma_mk, ...
                                      pilotseq, N, rho_d, tau_c, tau_p);
end
