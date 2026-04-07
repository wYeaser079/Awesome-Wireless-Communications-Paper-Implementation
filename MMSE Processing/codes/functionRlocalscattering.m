function R = functionRlocalscattering(N, theta, ASD_deg, antennaSpacing)
%functionRlocalscattering Generate spatial correlation matrix using the
%Gaussian local scattering model from [1, Sec. 2.6].
%
%   R = functionRlocalscattering(N, theta, ASD_deg, antennaSpacing)
%
%   Inputs:
%       N              - Number of antennas at the ULA
%       theta          - Nominal angle of arrival (radians)
%       ASD_deg        - Angular standard deviation (degrees)
%       antennaSpacing - Antenna spacing in wavelengths (default: 0.5)
%
%   Output:
%       R - N x N spatial correlation matrix
%
%   Reference:
%   [1] E. Bjornson, J. Hoydis, L. Sanguinetti, "Massive MIMO Networks:
%       Spectral, Energy, and Hardware Efficiency," Foundations and Trends
%       in Signal Processing, 2017.
%
% This is used by the paper:
% Emil Bjornson, Luca Sanguinetti, "Making Cell-Free Massive MIMO
% Competitive With MMSE Processing and Centralized Implementation,"
% IEEE Trans. Wireless Commun., vol. 19, no. 2, pp. 1257-1271, Feb. 2020.

if nargin < 4
    antennaSpacing = 0.5;  % Half-wavelength spacing
end

% Convert ASD to radians
ASD_rad = ASD_deg * pi / 180;

% Initialize the correlation matrix
R = zeros(N, N);

% Distance from the first antenna (in wavelengths)
distanceAntenna = antennaSpacing * (0:N-1);

% Numerical integration over the angular domain
% Use high-resolution integration for accuracy
numAngles = 1000;
angles = linspace(-pi, pi, numAngles);
dAngle = angles(2) - angles(1);

% Gaussian angular power spectrum centered at theta
% f(alpha) = (1 / (sqrt(2*pi) * ASD)) * exp(-(alpha - theta)^2 / (2*ASD^2))
f = (1 / (sqrt(2*pi) * ASD_rad)) * exp(-(angles - theta).^2 / (2 * ASD_rad^2));

% Build the correlation matrix
for row = 1:N
    for col = row:N
        d = distanceAntenna(col) - distanceAntenna(row);

        % Integral: R(row,col) = integral of exp(j*2*pi*d*sin(alpha)) * f(alpha) d(alpha)
        integrand = exp(1j * 2 * pi * d * sin(angles)) .* f;
        R(row, col) = sum(integrand) * dAngle;

        % Hermitian symmetry
        R(col, row) = conj(R(row, col));
    end
end

% Ensure the matrix is Hermitian positive semi-definite
R = (R + R') / 2;

end
