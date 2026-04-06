function R = functionRlocalscattering(N, theta, ASD, antennaSpacing, distribution)
%functionRlocalscattering Generate spatial correlation matrix using local
%scattering model with a Gaussian angular distribution.
%
%   R = functionRlocalscattering(N, theta, ASD, antennaSpacing, distribution)
%
%   Generates the N x N spatial correlation matrix for a Uniform Linear
%   Array (ULA) with N antennas, based on the local scattering model from
%   [1, Sec. 2.6]. The model assumes scatterers are distributed around
%   the UE with a given Angular Standard Deviation (ASD).
%
%   INPUTS:
%       N               - Number of antennas at the AP
%       theta           - Nominal angle of arrival (radians)
%       ASD             - Angular Standard Deviation (radians)
%       antennaSpacing  - Antenna spacing in wavelengths (default: 0.5)
%       distribution    - Angular distribution: 'Gaussian' (default),
%                         'Uniform', or 'Laplace'
%
%   OUTPUT:
%       R   - N x N spatial correlation matrix (Hermitian positive
%             semi-definite)
%
%   REFERENCES:
%   [1] E. Bjornson, J. Hoydis, L. Sanguinetti, "Massive MIMO Networks:
%       Spectral, Energy, and Hardware Efficiency," Foundations and Trends
%       in Signal Processing, vol. 11, no. 3-4, pp. 154-655, 2017.
%
%   This function follows the approach in the original code at:
%   https://github.com/emilbjornson/scalable-cell-free

% Default parameters
if nargin < 4
    antennaSpacing = 0.5; % Half-wavelength spacing
end
if nargin < 5
    distribution = 'Gaussian';
end

% For a single antenna, correlation matrix is just 1
if N == 1
    R = 1;
    return;
end

% Compute the first row of the correlation matrix (Toeplitz structure)
firstRow = zeros(N, 1);

% Set integration limits based on distribution type
switch lower(distribution)
    case 'gaussian'
        % Truncated Gaussian: integrate over [-20*ASD, 20*ASD]
        integrationRange = 20 * ASD;
    case 'uniform'
        % Uniform: angular spread over [-sqrt(3)*ASD, sqrt(3)*ASD]
        integrationRange = sqrt(3) * ASD;
    case 'laplace'
        % Laplace distribution
        integrationRange = 20 * ASD;
    otherwise
        error('Unknown angular distribution: %s', distribution);
end

% Numerical integration using trapezoidal rule
nPoints = 1001; % Number of integration points
deltaGrid = linspace(-integrationRange, integrationRange, nPoints);
dDelta = deltaGrid(2) - deltaGrid(1);

% Compute the angular power spectrum (pdf of the angular distribution)
switch lower(distribution)
    case 'gaussian'
        f_delta = exp(-deltaGrid.^2 / (2 * ASD^2)) / (sqrt(2*pi) * ASD);
    case 'uniform'
        f_delta = ones(size(deltaGrid)) / (2 * sqrt(3) * ASD);
    case 'laplace'
        b = ASD / sqrt(2);
        f_delta = exp(-abs(deltaGrid) / b) / (2 * b);
end

% Compute each element of the first row via numerical integration
for col = 1:N
    % Array response for antenna index (col-1) relative to reference
    integrand = exp(1j * 2 * pi * antennaSpacing * (col-1) ...
        * sin(theta + deltaGrid)) .* f_delta;

    % Trapezoidal integration
    firstRow(col) = dDelta * (sum(integrand) - 0.5*(integrand(1) + integrand(end)));
end

% Build the Toeplitz Hermitian correlation matrix from the first row
R = toeplitz(firstRow);

end
