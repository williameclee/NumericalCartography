%% NORMAL
% Calculates the normal map of a given height map.
%
% Syntax
%   normalMap = normal(Z)
%   normalMap = normal(Z, h)
%	normalMap = normal(__, VE)
%
% Description
%   normalMap = normal(Z)
%       calculates the normal map of the height map Z.
%   normalMap = normal(Z, h)
%       specifies the mesh size h of the height map Z.
%   normalMap = normal(__, VE)
%       specifies the z-factor of the height map Z.
%
% Input Arguments
%   Z - Height map
%       A 2D matrix representing the height map.
%   h - Mesh size
%       The default value is 1.
%   VE - Vertical exaggeration factor
%       The default value is 1 (no exaggeration).
%
% Output Arguments
%   normalMap - Normal map
%       A 3D matrix representing the normal map.
%
% Authored by
%   En-Chi Lee <williameclee@gmail.com>, 2024-06-29

function normalMap = normal(Z, varargin)
    %% Initialisation
    % Input parser
    p = inputParser;
    addRequired(p, 'Z', ...
        @(x) isnumeric(x) && ismatrix(x));
    addRequired(p, 'h', ...
        @(x) isscalar(x) && isnumeric(x) && (x > 0));
    addRequired(p, 'zFactor', ...
        @(x) isscalar(x) && isnumeric(x) && (x > 0));
    parse(p, Z, varargin{:});
    Z = p.Results.Z;
    h = p.Results.h;
    zFactor = p.Results.zFactor;

    % Preallocation
    normalMap = zeros([size(Z), 3]);

    %% Normal calculation
    Z = Z * zFactor / h;
    [normalMap(:, :, 1), normalMap(:, :, 2), normalMap(:, :, 3)] = ...
        surfnorm(Z);
end
