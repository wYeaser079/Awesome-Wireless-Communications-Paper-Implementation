function v = interference_vec(M, N, K, cdot, gn, bn, pilotseq, iUser)
%INTERFERENCE_VEC  Helper used inside the SCA-SOCP constraints of
%  Algorithm 1.  Builds the vector whose Euclidean norm bounds the
%  "interference + noise" piece of the SINR denominator for user IUSER
%  as a linear function of CDOT.
%
%  The argument GN should be passed already scaled (e.g.  sqrt(rho_d)*gamma
%  or gamma depending on which of constraints (27)/(36d) is being
%  encoded).  The same applies to BN = .*beta or sqrt(rho_d)*beta.
%
%  CDOT is either a numeric M-by-K matrix (at the tangent point) or a
%  CVX variable matrix (inside the SCA loop).
%
%  This mirrors interferencevector.m in the reference repo at
%  github.com/tranlenam/cellfreeMIMOenergyefficiency.

    % --- Pilot-contamination contributions (k' != iUser) -----------
    pilot_terms = [];
    for jUser = 1:K
        if jUser == iUser, continue; end
        inner = abs(pilotseq(:, jUser)' * pilotseq(:, iUser));
        bargamj = (gn(:, jUser) ./ bn(:, jUser)) .* bn(:, iUser) * inner;
        pilot_terms = [pilot_terms; bargamj' * cdot(:, jUser)];   %#ok<AGROW>
    end

    % --- Non-coherent interference contributions -------------------
    nc_terms = (1/sqrt(N)) * reshape( cdot .* sqrt(gn .* bn(:, iUser)), ...
                                      M*K, 1 );

    v = [pilot_terms; nc_terms];
end
