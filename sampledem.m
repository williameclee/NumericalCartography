%% SAMPLEDEM
% Generates a sample DEM
%
% Syntax
%   Z = sampledem
%   Z = sampledem(N)
%   [Z, X, Y, h, light] = sampledem(__)
%
% Description
%   Z = sampledem
% 		generates a 2D array of a sample DEM.
%  	Z = sampledem(N)
% 		generates a 2D array of a sample DEM with the size of N x N.
%   [Z, X, Y, h, light] = sampledem(__)
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
%   light - Light object
%
% See also
%   SAMPLEDEM2
%
% Authored by
%   En-Chi Lee <williameclee@gmail.com>, 2024-06-30
% Last modified by
%   En-Chi Lee <williameclee@gmail.com>, 2024-07-01

function varargout = sampledem(varargin)
    %% Initialistion
    p = inputParser;
    addOptional(p, 'MapSize', 2 ^ 9, ...
        @(x) isscalar(x) && isnumeric(x) && (x > 0));
    parse(p, varargin{:});
    N = p.Results.MapSize;

    %% Generating the DEM
    Z = imread('n24_e120_1arc_v3.tif');
    Z = Z(1201:2:3601, 1201:2:3601);
    Z = flip(Z);
    Z = double(Z);

    % Setup coordinates and interpolate
    XYrange = 111e3;
    X = linspace(0, XYrange, N);
    Y = linspace(0, XYrange, N);
    h = X(2) - X(1);

    Z = interp2(linspace(0, XYrange, size(Z, 2)), ...
        linspace(0, XYrange, size(Z, 1)), Z, ...
        X, Y', 'linear');

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
