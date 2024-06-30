%% GRAINYNOISE
% Adds a grainy noise to a normal map.
%
% Syntax
%   normalMap = grainynoise(normalMap)
%   normalMap = grainynoise(normalMap, strength)
%
% Description
%   normalMap = grainynoise(normalMap)
%       adds a grainy noise to the normal map.
%   normalMap = grainynoise(normalMap, strength)
%       specifies the strength of the noise.
%
% Input Arguments
%   normalMap - Normal map
%       A 3D matrix representing the normal map.
%   strength - Noise strength
%       The value must be in the range of (0, 1]. The default value is 0.1.
%
% Output Arguments
%   normalMap - Normal map
%       A 3D matrix representing the normal map with grainy noise.
%
% Authored by
%   En-Chi Lee <williameclee@gmail.com>, 2024-06-29

function normalMap = grainynoise(normalMap, varargin)
    %% Input parsing
    p = inputParser;
    addRequired(p, 'NormalMap', ...
        @(x) isnumeric(x) && (size(x, 3) == 3));
    addOptional(p, 'Strength', 0.1, ...
        @(x) isscalar(x) && isnumeric(x) && (x > 0));
    parse(p, normalMap, varargin{:});
    normalMap = p.Results.NormalMap;
    noiseStrength = p.Results.Strength;

    % Variable cleaning
    noiseStrength = min(1, noiseStrength);

    %% Generating the new normal map
    % Generate a noise map
    theta = 2 * pi * rand(size(normalMap, 1:2));
    phi = pi * rand(size(normalMap, 1:2));
    % Convert the noise map to a normal map
    noiseMap = cat(3, ...
        cos(theta) .* sin(phi), ...
        sin(theta) .* sin(phi), ...
        cos(phi));

    % Add the noise map to the normal map
    normalMap = normalMap + noiseStrength * noiseMap;
    normalMap = normalMap ./ vecnorm(normalMap, 2, 3);
end
