%% Main Runner Script
% =========================================================================
% Making Cell-Free Massive MIMO Competitive With MMSE Processing and
% Centralized Implementation
%
% Authors: Emil Bjornson and Luca Sanguinetti
% Published: IEEE Trans. Wireless Commun., vol. 19, no. 2, Feb. 2020
%
% Code Implementation by: [Your Name]
% =========================================================================
%
% This script provides an interface to run all figure-generating simulations.
% Each figure can be run independently by executing the corresponding script.
%
% IMPORTANT NOTES:
% - Figures 2(a), 2(b), 3, and 5 require significant computation time
%   due to Monte Carlo simulations. Reduce nbrOfSetups and nbrOfRealizations
%   in each script for faster (but less accurate) results.
% - Figure 4 uses the three-slope path loss model from Ngo et al. [16]
% - Figure 6 is purely analytical and runs instantly
% - The default parameters are set lower than the paper for faster execution.
%   For paper-quality results, increase to: nbrOfSetups=200, nbrOfRealizations=1000
%
% FIGURES REPRODUCED:
%   Figure 2(a) - CDF of per-user SE, MMSE/L-MMSE, L=400, N=1
%   Figure 2(b) - CDF of per-user SE, MMSE/L-MMSE, L=100, N=4
%   Figure 3    - CDF of per-user SE, MR combining, L=400, N=1
%   Figure 4    - Revisiting Ngo et al. with three-slope model
%   Figure 5    - Sum SE: MMSE-SIC vs linear combining
%   Figure 6    - Fronthaul signaling load analysis

clear; close all; clc;

fprintf('================================================================\n');
fprintf('  Making Cell-Free Massive MIMO Competitive With MMSE Processing\n');
fprintf('  and Centralized Implementation\n');
fprintf('  Bjornson & Sanguinetti, IEEE TWC, 2020\n');
fprintf('================================================================\n\n');

fprintf('Select which figure(s) to generate:\n');
fprintf('  1 - Figure 2(a) and 3: L=400, N=1 (MMSE + MR comparison)\n');
fprintf('  2 - Figure 2(b): L=100, N=4 (MMSE only)\n');
fprintf('  3 - Figure 4: Ngo et al. comparison (three-slope model)\n');
fprintf('  4 - Figure 5: Sum SE with MMSE-SIC\n');
fprintf('  5 - Figure 6: Fronthaul analysis (instant)\n');
fprintf('  6 - ALL figures (may take several hours)\n');
fprintf('  0 - Exit\n\n');

choice = input('Enter your choice (0-6): ');

switch choice
    case 1
        fprintf('\nRunning Figures 2(a) and 3...\n');
        simulationFigure2a_3;
    case 2
        fprintf('\nRunning Figure 2(b)...\n');
        simulationFigure2b;
    case 3
        fprintf('\nRunning Figure 4...\n');
        simulationFigure4;
    case 4
        fprintf('\nRunning Figure 5...\n');
        simulationFigure5;
    case 5
        fprintf('\nRunning Figure 6...\n');
        simulationFigure6;
    case 6
        fprintf('\nRunning ALL figures...\n');
        fprintf('\n--- Figures 2(a) and 3 ---\n');
        simulationFigure2a_3;
        fprintf('\n--- Figure 2(b) ---\n');
        simulationFigure2b;
        fprintf('\n--- Figure 4 ---\n');
        simulationFigure4;
        fprintf('\n--- Figure 5 ---\n');
        simulationFigure5;
        fprintf('\n--- Figure 6 ---\n');
        simulationFigure6;
        fprintf('\n=== ALL FIGURES COMPLETE ===\n');
    case 0
        fprintf('Exiting.\n');
    otherwise
        fprintf('Invalid choice.\n');
end
