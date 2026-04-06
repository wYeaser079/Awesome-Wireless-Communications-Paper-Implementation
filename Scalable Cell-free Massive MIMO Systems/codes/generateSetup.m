function [gainOverNoisedB, R, pilotIndex, D, D_all, APpositions, UEpositions, masterAPs] = ...
    generateSetup(L, K, N, tau_p, nbrOfSetups, seed)
%generateSetup Generate random network layouts for Scalable Cell-Free
%Massive MIMO simulations.
%
%   [gainOverNoisedB, R, pilotIndex, D, D_all, APpositions, UEpositions, masterAPs] = ...
%       generateSetup(L, K, N, tau_p, nbrOfSetups, seed)
%
%   Implements the three-step algorithm for joint initial access, pilot
%   assignment, and cluster formation from Section V-A of [1].
%
%   INPUTS:
%       L            - Number of APs
%       K            - Number of UEs
%       N            - Number of antennas per AP
%       tau_p        - Number of orthogonal pilot sequences
%       nbrOfSetups  - Number of random network realizations
%       seed         - (Optional) Random seed for reproducibility
%
%   OUTPUTS:
%       gainOverNoisedB - L x K x nbrOfSetups array of channel gains
%                         over noise in dB
%       R               - N x N x L x K x nbrOfSetups array of spatial
%                         correlation matrices
%       pilotIndex      - K x nbrOfSetups matrix of pilot assignments
%       D               - L x K x nbrOfSetups binary service indicator
%                         (DCC matrix, scalable)
%       D_all           - L x K matrix of all-ones (non-scalable benchmark)
%       APpositions     - L x nbrOfSetups complex AP positions
%       UEpositions     - K x nbrOfSetups complex UE positions
%       masterAPs       - K x nbrOfSetups master AP indices
%
%   REFERENCES:
%   [1] E. Bjornson, L. Sanguinetti, "Scalable Cell-Free Massive MIMO
%       Systems," IEEE Trans. Commun., vol. 68, no. 7, pp. 4247-4261, 2020.
%   [2] E. Bjornson, J. Hoydis, L. Sanguinetti, "Massive MIMO Networks,"
%       Foundations and Trends in Signal Processing, vol. 11, no. 3-4, 2017.

if nargin >= 6 && ~isempty(seed)
    rng(seed);
end

%% Simulation Parameters

% Coverage area dimensions (meters)
squareLength = 2000; % 2 x 2 km

% Communication bandwidth (Hz)
B = 20e6;

% Noise figure at the APs (dB)
noiseFigure = 7;

% Noise power (dBm)
noiseVariancedBm = -174 + 10*log10(B) + noiseFigure;

% Pathloss parameters (three-slope model from [2, Sec. 4.1.3])
% Using the convention: PL(dB) = constantTerm - alpha * log10(d_3D)
% where alpha = 10 * pathlossExponent. Reference: original code uses 3.76.
alpha = 37.6; % = 10 * 3.76 pathloss exponent
constantTerm = -35.3; % Pathloss at reference distance (dB)
sigma_sf = 10; % Shadow fading standard deviation (dB)

% Note: In the original code by Bjornson, the model uses a simplified
% version. The three-slope pathloss in dB at distance d (km) is:
%   PL = constantTerm - alpha * log10(d_3D)
% where d_3D includes the height difference.

% AP height above UE level (meters)
APheight = 10;

% Antenna spacing (in wavelengths)
antennaSpacing = 0.5;

% Angular standard deviation around the UEs (degrees -> radians)
ASD_deg = 20;
ASD = ASD_deg * pi / 180;

% Cluster formation threshold (dB)
% An AP serves a UE only if its channel gain is within this threshold
% of the master AP's channel gain
clusterThreshold = -40; % dB

%% Non-scalable benchmark (all APs serve all UEs)
D_all = ones(L, K);

%% Pre-allocate output arrays
gainOverNoisedB = zeros(L, K, nbrOfSetups);
R = zeros(N, N, L, K, nbrOfSetups);
pilotIndex = zeros(K, nbrOfSetups);
D = zeros(L, K, nbrOfSetups);
APpositions = zeros(L, nbrOfSetups);
UEpositions = zeros(K, nbrOfSetups);
masterAPs = zeros(K, nbrOfSetups);

%% Generate random network realizations
for n = 1:nbrOfSetups

    %% Step 1: Random placement of APs and UEs
    % APs uniformly distributed in the square
    APpos = squareLength * (rand(L, 1) + 1j * rand(L, 1));
    APpositions(:, n) = APpos;

    % UEs uniformly distributed in the square
    UEpos = squareLength * (rand(K, 1) + 1j * rand(K, 1));
    UEpositions(:, n) = UEpos;

    %% Step 2: Compute channel gains with wrap-around

    % Wrap-around: create 9 copies of the network (3x3 grid)
    wrapPositions = [0; squareLength; -squareLength; ...
        1j*squareLength; -1j*squareLength; ...
        squareLength + 1j*squareLength; ...
        squareLength - 1j*squareLength; ...
        -squareLength + 1j*squareLength; ...
        -squareLength - 1j*squareLength];

    % Compute channel gains for each AP-UE pair
    for l = 1:L
        for k = 1:K

            % Compute distances to all 9 wrap-around copies
            distVec = abs(APpos(l) - UEpos(k) + wrapPositions);

            % Include height difference for 3D distance
            dist3D = sqrt(distVec.^2 + APheight^2);

            % Find the closest copy (wrap-around distance)
            [minDist3D, bestCopy] = min(dist3D);

            % Compute pathloss (dB) using simplified model
            % Based on [2, Sec. 4.1.3]: PL = constantTerm - alpha * log10(d)
            gainOverNoisedB(l, k, n) = constantTerm ...
                - alpha * log10(minDist3D) ...
                + sigma_sf * randn() ...
                - noiseVariancedBm;

            % Compute angle of arrival at AP l from UE k
            % (using the closest wrap-around copy)
            APtoUE = UEpos(k) + wrapPositions(bestCopy) - APpos(l);
            angleOfArrival = angle(APtoUE);

            % Generate spatial correlation matrix
            if N > 1
                R(:, :, l, k, n) = db2pow(gainOverNoisedB(l, k, n)) * ...
                    functionRlocalscattering(N, angleOfArrival, ASD, antennaSpacing);
            else
                % Single antenna: scalar channel gain
                R(:, :, l, k, n) = db2pow(gainOverNoisedB(l, k, n));
            end
        end
    end

    %% Step 3: Three-step algorithm for pilot assignment and cluster formation

    % === Step 1 of Algorithm: Each UE selects its Master AP ===
    % Master AP = AP with largest average channel gain (beta_l)
    [~, masterAPidx] = max(gainOverNoisedB(:, :, n), [], 1);
    masterAPs(:, n) = masterAPidx(:);

    % === Step 2 of Algorithm: Pilot assignment ===
    % First tau_p UEs get orthogonal pilots
    pilotAssignment = zeros(K, 1);

    % Assign first tau_p UEs to orthogonal pilots (1 through tau_p)
    for k = 1:min(K, tau_p)
        pilotAssignment(k) = k;
    end

    % Assign remaining UEs greedily: pick pilot with least interference
    % at the Master AP
    for k = (tau_p+1):K
        masterAP_k = masterAPidx(k);

        % For each pilot, compute the total received pilot power at Master AP
        pilotInterference = zeros(tau_p, 1);
        for t = 1:tau_p
            % Find UEs already assigned to pilot t
            usersOnPilot = find(pilotAssignment == t);
            for u = usersOnPilot'
                % Add the channel gain from this UE to the Master AP
                pilotInterference(t) = pilotInterference(t) + ...
                    db2pow(gainOverNoisedB(masterAP_k, u, n));
            end
        end

        % Assign pilot with minimum interference (Eq. 15)
        [~, bestPilot] = min(pilotInterference);
        pilotAssignment(k) = bestPilot;
    end

    pilotIndex(:, n) = pilotAssignment;

    %% === Step 3 of Algorithm: Cluster Formation (DCC) ===
    % Each AP serves at most one UE per pilot

    % First, each Master AP must serve its UE
    for k = 1:K
        D(masterAPidx(k), k, n) = 1;
    end

    % Then, neighboring APs decide whether to join the cluster
    for l = 1:L
        for t = 1:tau_p
            % Find UEs on pilot t
            usersOnPilot = find(pilotAssignment == t);

            if isempty(usersOnPilot)
                continue;
            end

            % Check if AP l is already the master AP for any UE on pilot t
            isMasterForPilot = false;
            for u = usersOnPilot'
                if masterAPidx(u) == l
                    isMasterForPilot = true;
                    break;
                end
            end

            if isMasterForPilot
                % Already serving as master; D already set
                continue;
            end

            % Find the UE on this pilot with the best channel to AP l
            bestGain = -inf;
            bestUE = -1;
            for u = usersOnPilot'
                if gainOverNoisedB(l, u, n) > bestGain
                    bestGain = gainOverNoisedB(l, u, n);
                    bestUE = u;
                end
            end

            % Check if the channel is within the cluster threshold
            % relative to the master AP's channel
            if bestUE > 0
                masterGain = gainOverNoisedB(masterAPidx(bestUE), bestUE, n);
                if bestGain >= masterGain + clusterThreshold
                    D(l, bestUE, n) = 1;
                end
            end
        end
    end

end

end
