%% SAMPLEDEM2
% Generates a sample DEM
%
% Syntax
%   Z = sampledem2
%   Z = sampledem2(N)
%   [Z, X, Y, h, light] = sampledem2(__)
%
% Description
%   Z = sampledem2
% 		generates a 2D array of a sample DEM.
%  	Z = sampledem2(N)
% 		generates a 2D array of a sample DEM with the size of N x N.
%   [Z, X, Y, h, light] = sampledem2(__)
% 		returns the sample DEM, the X and Y coordinates, the mesh size, and
% 		the light direction.
%
% Input Arguments
%   N - Map size
%       The default value is 512 (2^9).
%
% Output Arguments
%   Z - Sample DEM
%   X - X coordinates
%   Y - Y coordinates
%   h - Mesh size
%   light - Light direction
%
% See also
%   SAMPLEDEM
%
% Authored by
%   En-Chi Lee <williameclee@gmail.com>, 2024-06-30
% Last modified by
%   En-Chi Lee <williameclee@gmail.com>, 2024-07-01

function varargout = sampledem2(varargin)
    %% Initialistion
    p = inputParser;
    addOptional(p, 'MapSize', 2 ^ 9, ...
        @(x) isscalar(x) && isnumeric(x) && (x > 0));
    parse(p, varargin{:});
    N = p.Results.MapSize;

    %% Generating the DEM
    X = linspace(0, 1, N);
    Y = linspace(0, 1, N);
    h = X(2) - X(1);
    Z = peaks(N);
    Z = Z / max(abs(Z(:))) * 0.3;

    light = [-1, 1, 0.5];
    light = LightSource('parallel', light);

    %% Returning requested outputs
    if nargout == 0
        figure('Name', 'Sample DEM')
        surf(X, Y, Z, ...
            'EdgeColor', 'none')
        axis equal
        axis off
        colormap(gray)
        return
    end

    varargout = {Z, X, Y, h, light};
end
