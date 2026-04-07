function PL_dB = pathloss_threeslope(dist3D)
%pathloss_threeslope Compute path loss using the three-slope model from
%Ngo et al. [16], based on the COST-Hata model (Eq. 39 in the paper).
%
%   PL_dB = pathloss_threeslope(dist3D)
%
%   Input:
%       dist3D - 3D distance in meters
%
%   Output:
%       PL_dB  - Path loss in dB (negative value)
%
% Three-slope model:
%   d < 10m:       PL = -L0 (constant)
%   10m <= d < 50m: PL = -L0 - 20*log10(d/10)
%   d >= 50m:      PL = -L0 - 20*log10(50/10) - 35*log10(d/50)
%
% Where L0 = 140.7151 dB (reference loss at 1m)
% The constants are chosen for continuity at the breakpoints.

% Reference constant (adjusted for the Ngo et al. setup)
% -61.2 dB at 10m reference => L0 = 81.2 dB at < 10m
% But the three-slope model from [16] uses:

if dist3D < 10
    PL_dB = -81.2;
elseif dist3D < 50
    PL_dB = -61.2 - 20 * log10(dist3D);
else
    PL_dB = -35.7 - 35 * log10(dist3D);
end

end
