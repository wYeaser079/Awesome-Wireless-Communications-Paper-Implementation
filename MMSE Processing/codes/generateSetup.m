function [R_AP, R_BS, gainOverNoisedB_AP, gainOverNoisedB_BS, ...
    pilotIndex, APpositions, BSpositions, UEpositions] = ...
    generateSetup(L, K, N, tau_p, nbrOfBSs, M)
%generateSetup Generate the network topology and channel statistics for
%cell-free mMIMO and cellular mMIMO comparison.
%
%   Uses the 3GPP Urban Microcell propagation model (TR 36.814).
%
%   Inputs:
%       L       - Number of APs in the cell-free system
%       K       - Number of UEs
%       N       - Number of antennas per AP
%       tau_p   - Number of orthogonal pilots
%       nbrOfBSs- Number of BSs in cellular system (default: 4)
%       M       - Number of antennas per BS (default: LN/nbrOfBSs)
%
%   Outputs:
%       R_AP              - N x N x L x K array of AP spatial correlation matrices
%       R_BS              - M x M x nbrOfBSs x K x nbrOfBSs array of BS correlation
%       gainOverNoisedB_AP- L x K matrix of gain-over-noise in dB (APs)
%       gainOverNoisedB_BS- nbrOfBSs x K x nbrOfBSs matrix (BSs)
%       pilotIndex        - K x 1 vector of pilot indices
%       APpositions       - L x 1 complex vector of AP positions
%       BSpositions       - nbrOfBSs x 1 complex vector of BS positions
%       UEpositions       - K x 1 complex vector of UE positions
%
% Reference:
% Emil Bjornson, Luca Sanguinetti, "Making Cell-Free Massive MIMO
% Competitive With MMSE Processing and Centralized Implementation,"
% IEEE Trans. Wireless Commun., vol. 19, no. 2, pp. 1257-1271, Feb. 2020.

if nargin < 5
    nbrOfBSs = 4;
end
if nargin < 6
    M = round(L * N / nbrOfBSs);
end

%% Network parameters
squareLength = 1000;  % Side length of the square area in meters
APheight = 10;        % AP height in meters
BSheight = 10;        % BS height in meters (same as APs for fair comparison)
UEheight = 1.5;       % UE height in meters

% Propagation model parameters (3GPP Urban Microcell)
ASD_deg = 15;           % Angular standard deviation in degrees
shadowFadingStd = 4;    % Shadow fading standard deviation in dB
decorrelationDist = 9;  % Decorrelation distance in meters
B = 20e6;               % Bandwidth in Hz
noiseFigure = 5;        % Noise figure in dB

% Noise power
noiseVariancedBm = -174 + 10*log10(B) + noiseFigure;

%% Place APs on a grid
nbrAPsPerDim = sqrt(L);
if nbrAPsPerDim == floor(nbrAPsPerDim)
    % APs on a regular grid
    APspacing = squareLength / nbrAPsPerDim;
    tmp = APspacing/2 : APspacing : squareLength - APspacing/2;
    [APx, APy] = meshgrid(tmp, tmp);
    APpositions = APx(:) + 1j * APy(:);
else
    % If not perfect square, place randomly
    APpositions = squareLength * (rand(L, 1) + 1j * rand(L, 1));
end

%% Place BSs on a grid (2x2 for 4 BSs)
nbrBSsPerDim = sqrt(nbrOfBSs);
BSspacing = squareLength / nbrBSsPerDim;
tmp = BSspacing/2 : BSspacing : squareLength - BSspacing/2;
[BSx, BSy] = meshgrid(tmp, tmp);
BSpositions = BSx(:) + 1j * BSy(:);

%% Place UEs randomly
UEpositions = squareLength * (rand(K, 1) + 1j * rand(K, 1));

%% Generate shadow fading realizations
% Shadow fading correlation between UEs for same AP/BS
shadowFading_AP = zeros(L, K);
shadowFading_BS = zeros(nbrOfBSs, K);

for l = 1:L
    % Generate correlated shadow fading for all UEs to AP l
    shadowFading_AP(l, :) = generateCorrelatedShadowFading(...
        K, UEpositions, shadowFadingStd, decorrelationDist);
end

for j = 1:nbrOfBSs
    shadowFading_BS(j, :) = generateCorrelatedShadowFading(...
        K, UEpositions, shadowFadingStd, decorrelationDist);
end

%% Compute path loss and gain-over-noise for APs
gainOverNoisedB_AP = zeros(L, K);
R_AP = zeros(N, N, L, K);

for l = 1:L
    for k = 1:K
        % Distance with wraparound
        dist = computeWraparoundDistance(UEpositions(k), APpositions(l), squareLength);

        % 3D distance (including height difference)
        dist3D = sqrt(dist^2 + (APheight - UEheight)^2);

        % 3GPP Urban Microcell path loss (Eq. 37)
        pathLossdB = -30.5 - 36.7 * log10(dist3D);

        % Add shadow fading
        gainOverNoisedB_AP(l, k) = pathLossdB + shadowFading_AP(l, k) - noiseVariancedBm;

        % Angle of arrival
        theta = angle(UEpositions(k) - APpositions(l));

        % Spatial correlation matrix
        if N > 1
            R_AP(:, :, l, k) = db2pow(gainOverNoisedB_AP(l, k)) * ...
                functionRlocalscattering(N, theta, ASD_deg);
        else
            R_AP(:, :, l, k) = db2pow(gainOverNoisedB_AP(l, k));
        end
    end
end

%% Compute path loss and gain-over-noise for BSs
% Each UE belongs to one cell; we need channels to all BSs
gainOverNoisedB_BS = zeros(nbrOfBSs, K, nbrOfBSs);
R_BS = zeros(M, M, nbrOfBSs, K, nbrOfBSs);

for j = 1:nbrOfBSs     % BS index
    for l = 1:nbrOfBSs  % Cell index (which cell the UE is in)
        for k = 1:K
            dist = computeWraparoundDistance(UEpositions(k), BSpositions(j), squareLength);
            dist3D = sqrt(dist^2 + (BSheight - UEheight)^2);

            pathLossdB = -30.5 - 36.7 * log10(dist3D);
            gainOverNoisedB_BS(j, k, l) = pathLossdB + shadowFading_BS(j, k) - noiseVariancedBm;

            theta = angle(UEpositions(k) - BSpositions(j));
            if M > 1
                R_BS(:, :, j, k, l) = db2pow(gainOverNoisedB_BS(j, k, l)) * ...
                    functionRlocalscattering(M, theta, ASD_deg);
            else
                R_BS(:, :, j, k, l) = db2pow(gainOverNoisedB_BS(j, k, l));
            end
        end
    end
end

%% Assign UEs to cells (for cellular system)
% Each UE is assigned to the BS with the strongest average channel
% (This is implicit in the BS indexing)

%% Pilot assignment for cell-free system
% Assign pilots to minimize pilot contamination
pilotIndex = zeros(K, 1);

% Simple greedy assignment based on spatial separation
for k = 1:K
    if k <= tau_p
        % First tau_p UEs get unique pilots
        pilotIndex(k) = k;
    else
        % Remaining UEs: assign pilot that minimizes contamination
        % by choosing the pilot whose existing users are farthest away
        minContamination = inf(tau_p, 1);

        for t = 1:tau_p
            usersWithPilot = find(pilotIndex(1:k-1) == t);
            if isempty(usersWithPilot)
                minContamination(t) = 0;
            else
                % Sum of gain-over-noise products (contamination metric)
                contamination = 0;
                for u = usersWithPilot'
                    contamination = contamination + ...
                        sum(db2pow(gainOverNoisedB_AP(:, k)) .* db2pow(gainOverNoisedB_AP(:, u)));
                end
                minContamination(t) = contamination;
            end
        end

        [~, pilotIndex(k)] = min(minContamination);
    end
end

end


%% Helper function: Compute wraparound distance
function dist = computeWraparoundDistance(pos1, pos2, squareLength)
% Compute distance with wraparound topology

dx = abs(real(pos1) - real(pos2));
dy = abs(imag(pos1) - imag(pos2));

dx = min(dx, squareLength - dx);
dy = min(dy, squareLength - dy);

dist = sqrt(dx^2 + dy^2);
end


%% Helper function: Generate correlated shadow fading
function SF = generateCorrelatedShadowFading(K, UEpositions, std_dB, decorrelationDist)
% Generate spatially correlated shadow fading realizations

% Compute distance matrix between UEs
distMatrix = zeros(K, K);
for i = 1:K
    for j = i+1:K
        distMatrix(i, j) = abs(UEpositions(i) - UEpositions(j));
        distMatrix(j, i) = distMatrix(i, j);
    end
end

% Shadow fading correlation matrix (Eq. 38)
corrMatrix = std_dB^2 * 2.^(-distMatrix / decorrelationDist);

% Generate correlated Gaussian random variables
L_chol = chol(corrMatrix + 1e-10 * eye(K), 'lower');
SF = (L_chol * randn(K, 1))';

end
