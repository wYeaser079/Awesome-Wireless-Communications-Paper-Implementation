function [gainOverNoisedB, pilotIndex, APpositions, UEpositions] = ...
    generateSetup_threeslope(L, K, tau_p)
%generateSetup_threeslope Generate network setup using the three-slope
%path loss model from Ngo et al. [16] for Figure 4 comparison.
%
%   Uses the COST-Hata three-slope model (Eq. 39 in the paper).
%
%   Inputs:
%       L     - Number of APs (single antenna, N=1)
%       K     - Number of UEs
%       tau_p - Number of orthogonal pilots
%
%   Outputs:
%       gainOverNoisedB - L x K matrix of channel gain over noise in dB
%       pilotIndex      - K x 1 vector of pilot indices
%       APpositions     - L x 1 complex vector of AP positions
%       UEpositions     - K x 1 complex vector of UE positions
%
% Reference:
% H. Q. Ngo, A. Ashikhmin, H. Yang, E. G. Larsson, and T. L. Marzetta,
% "Cell-Free Massive MIMO Versus Small Cells," IEEE TWC, 2017.

%% Network parameters
squareLength = 1000;  % 1 km x 1 km area
APheight = 15;        % AP height (matching Ngo et al.)
UEheight = 1.65;      % UE height
B = 20e6;             % Bandwidth
noiseFigure = 9;      % Noise figure (matching Ngo et al.)
shadowFadingStd = 8;  % Shadow fading std (matching Ngo et al.)
decorrelationDist = 100; % Decorrelation distance (matching Ngo et al.)

% Noise power
noiseVariancedBm = -174 + 10*log10(B) + noiseFigure;

%% Place APs on a grid
nbrAPsPerDim = sqrt(L);
if nbrAPsPerDim == floor(nbrAPsPerDim)
    APspacing = squareLength / nbrAPsPerDim;
    tmp = APspacing/2 : APspacing : squareLength - APspacing/2;
    [APx, APy] = meshgrid(tmp, tmp);
    APpositions = APx(:) + 1j * APy(:);
else
    APpositions = squareLength * (rand(L, 1) + 1j * rand(L, 1));
end

%% Place UEs randomly
UEpositions = squareLength * (rand(K, 1) + 1j * rand(K, 1));

%% Generate shadow fading
shadowFading = zeros(L, K);
for l = 1:L
    % Correlated shadow fading
    distMatrix = zeros(K, K);
    for i = 1:K
        for j2 = i+1:K
            distMatrix(i, j2) = abs(UEpositions(i) - UEpositions(j2));
            distMatrix(j2, i) = distMatrix(i, j2);
        end
    end
    corrMatrix = shadowFadingStd^2 * 2.^(-distMatrix / decorrelationDist);
    L_chol = chol(corrMatrix + 1e-10 * eye(K), 'lower');
    shadowFading(l, :) = (L_chol * randn(K, 1))';
end

%% Compute path loss using three-slope model
gainOverNoisedB = zeros(L, K);

for l = 1:L
    for k = 1:K
        % Wraparound distance
        dx = abs(real(UEpositions(k)) - real(APpositions(l)));
        dy = abs(imag(UEpositions(k)) - imag(APpositions(l)));
        dx = min(dx, squareLength - dx);
        dy = min(dy, squareLength - dy);
        dist2D = sqrt(dx^2 + dy^2);

        % 3D distance
        dist3D = sqrt(dist2D^2 + (APheight - UEheight)^2);

        % Three-slope path loss (Eq. 39)
        pathLossdB = pathloss_threeslope(dist3D);

        % Add shadow fading (only for d >= 50m in Ngo model)
        if dist3D >= 50
            gainOverNoisedB(l, k) = pathLossdB + shadowFading(l, k) - noiseVariancedBm;
        else
            gainOverNoisedB(l, k) = pathLossdB - noiseVariancedBm;
        end
    end
end

%% Greedy pilot assignment (from Ngo et al. [16, Sec. IV.A])
pilotIndex = zeros(K, 1);

for k = 1:K
    if k <= tau_p
        pilotIndex(k) = k;
    else
        % Assign pilot that minimizes contamination
        minContamination = inf(tau_p, 1);
        for t = 1:tau_p
            usersWithPilot = find(pilotIndex(1:k-1) == t);
            if isempty(usersWithPilot)
                minContamination(t) = 0;
            else
                contamination = 0;
                for u = usersWithPilot'
                    contamination = contamination + ...
                        sum(db2pow(gainOverNoisedB(:, k)) .* db2pow(gainOverNoisedB(:, u)));
                end
                minContamination(t) = contamination;
            end
        end
        [~, pilotIndex(k)] = min(minContamination);
    end
end

end
