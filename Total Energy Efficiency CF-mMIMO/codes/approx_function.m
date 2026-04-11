function F = approx_function(M, N, K, gamma_mk, beta_mk, pilotseq, rho_d, ...
                              cdot, udot, cn, un, theta, iUser)
%APPROX_FUNCTION  Convex (affine) Taylor lower bound for the
%  quadratic-over-linear term
%
%      f(c, u_k) = [ rho_d N^2 |gamma_kk^T c_k|^2
%                    + rho_d N^2 sum_{k'!=k}|gamma_kk'^T c_k'|^2
%                    + rho_d N sum_{k'}|| D_{k',k} c_{k'} ||^2
%                    + 1 ] / u_k,
%
%  expanded around the tangent point (CN, UN).  See Eqs. (30)-(33) of
%  Ngo-Tran-Duong-Matthaiou-Larsson IEEE TGCN 2018.
%
%  Inside the SCA-Charnes-Cooper reformulation, all variables are
%  scaled by THETA (so CDOT = c*theta, UDOT = u*theta).  This function
%  returns the Taylor bound expressed in the scaled variables.

    % Scalar term f(cn, un) ------------------------------------------
    interf0 = build_interf(M, N, K, cn, gamma_mk, beta_mk, pilotseq, iUser);
    f1 = rho_d * (N^2) * norm(interf0)^2 ...
       + 1 ...
       + rho_d * (N^2) * (gamma_mk(:, iUser)' * cn(:, iUser))^2;

    % Charnes-Cooper lifted constant term (in scaled variables) -----
    F = theta * f1 / un(iUser) ...
        - (f1 / (un(iUser))^2) * (udot(iUser) - theta * un(iUser));

    % Linearised contribution from every user k ---------------------
    for k = 1:K
        % gamma_bar for the (k, iUser) pair  ---------------------
        % entries: |phi_k^H phi_iUser| * gamma_{m,k} * beta_{m,iUser}/beta_{m,k}
        gbar = (gamma_mk(:, k) ./ beta_mk(:, k)) .* beta_mk(:, iUser) ...
               * abs(pilotseq(:, k)' * pilotseq(:, iUser));
        Dk2  = diag( gamma_mk(:, k) .* beta_mk(:, iUser) );

        M_k  = (N^2) * (gbar * gbar') + N * Dk2;
        grad = cn(:, k)' * M_k;     % 1 x M row

        F = F + (2 * rho_d / un(iUser)) * grad * (cdot(:, k) - theta * cn(:, k));
    end
end

% ----------------------------------------------------------------------
function v = build_interf(M, N, K, c, gamma_mk, beta_mk, pilotseq, iUser)
%BUILD_INTERF  Like interference_vec but for pure numeric CN arguments.
%   Uses un-scaled gamma, beta (the sqrt(rho_d) factor stays outside).

    v = [];
    for jUser = 1:K
        if jUser == iUser, continue; end
        inner = abs(pilotseq(:, jUser)' * pilotseq(:, iUser));
        bargamj = (gamma_mk(:, jUser) ./ beta_mk(:, jUser)) ...
                  .* beta_mk(:, iUser) * inner;
        v = [v; bargamj' * c(:, jUser)];   %#ok<AGROW>
    end

    nc = (1/sqrt(N)) * reshape( c .* sqrt(gamma_mk .* beta_mk(:, iUser)), ...
                                M*K, 1 );
    v = [v; nc];
end
