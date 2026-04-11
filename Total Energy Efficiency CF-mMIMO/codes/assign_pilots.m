function pilotseq = assign_pilots(K, tau_p)
%ASSIGN_PILOTS  Random pilot assignment from a pool of TAU_P orthonormal
%  length-TAU_P sequences.
%
%  OUTPUT:
%     pilotseq - tau_p-by-K matrix whose columns are the pilot sequences.
%                Each column has unit norm.  If tau_p >= K, every user
%                gets a distinct orthogonal sequence (no pilot
%                contamination).  If tau_p < K, some users must reuse
%                a pilot and contamination appears.
%
%  This matches Section VI-A of Ngo et al. 2018 and the reference code
%  at github.com/tranlenam/cellfreeMIMOenergyefficiency.

    [U, ~, ~] = svd(randn(tau_p, tau_p));   % tau_p orthonormal cols

    if tau_p >= K
        pilotseq = U(:, 1:K);
    else
        pilotseq = zeros(tau_p, K);
        pilotseq(:, 1:tau_p) = U;
        for iUser = (tau_p+1):K
            pilotseq(:, iUser) = U(:, randi([1, tau_p]));
        end
    end
end
