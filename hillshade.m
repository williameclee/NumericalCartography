%% HILLSHADE
% Calculates the hillshade of a DEM.
%
% Syntax
%   hillshadeMap = hillshade(Z)
%   hillshadeMap = hillshade(Z, L)
%   hillshadeMap = hillshade(Z, __, h)
%   hillshadeMap = hillshade(Z, __, VE)
%   hillshadeMap = hillshade(Z, __, Name, Value)
%   hillshadeMap = hillshade('demo')
%   [hillshadeMap, L, r, N] = hillshade(__)
%
% Description
%   hillshadeMap = hillshade(Z)
%       calculates the hillshade of the DEM Z with the default light
%       azimuth and altitude of 45 degrees.
%       If Z is a normal map, the normal map is used for the calculation.
%   hillshadeMap = hillshade(Z, L)
%       calculates the hillshade of the DEM Z with the light source.
%   hillshadeMap = hillshade(Z, __, h)
%       calculates the hillshade of the DEM Z with the mesh size h. The DEM
%       is assumed to have the same mesh size in both dimensions.
%   hillshadeMap = hillshade(Z, __, VE)
%       calculates the hillshade of the DEM Z with the vertical
%       exaggeration factor VE.
%   hillshadeMap = hillshade('demo')
%       displays the hillshade map of the peaks function.
%   hillshade(__)
%       displays the hillshade map hs.
%   [hillshadeMap, L, r, N] = hillshade(__)
%       returns the hillshade map hs, the light direction L, the flat
%       reflectance r, and the normal map N.
%
% Input Aarguments
%   Z - DEM matrix (Nx-by-Ny) or normal map (Nx-by-Ny-by-3)
%   L - Light source
%       - LightSource object.
%       - 1-by-2 vector [azimuth, altitude] in degrees.
%       - 1-by-3 vector [Lx, Ly, Lz] of the (invert) light direction.
%       The default value is [-1, 1, 1].
%   h - Mesh size
%       The default value is 1.
%   VE - Vertical exaggeration factor
%       The default value is 1 (no exaggeration).
%   ShadingMethod - Shading method
%       The default value is 'Clamped' (Lambert). The other options are
%       'Half' and 'Soft' (Lambert).
%   IndexScheme - Indexing scheme
%       The default value is 'meshgrid'. The other option is 'ndgrid'.
%
% Output Arguments
%   hillshadeMap - Hillshade map
%   L - Light object
%   r - Flat reflectance
%       The flat reflectance is the reflectance of a flat surface facing
%       upwards.
%   N - Normal map
%
% See also
%   NORMAL, LIGHTSOURCE
%
% Authored by
%   En-Chi Lee <williameclee@gmail.com>, 2024-06-23
% Last modified by
%   En-Chi Lee <williameclee@gmail.com>, 2024-07-01

function varargout = hillshade(Z, varargin)
    %% Initialisation and preallocation
    % Make sure subfunctions in the path
    addpath(fullfile(fileparts(mfilename('fullpath')), ...
    'unit-sphere-transformations'));
    addpath(fullfile(fileparts(mfilename('fullpath')), ...
    'aux'));

    % Input parser
    p = inputParser;
    addRequired(p, 'HeightMap', ...
        @(x) isnumeric(x) || ...
        (ischar(x) && strcmp(x, 'demo')));
    addOptional(p, 'Light', [-1, 1, 1], ...
        @(x) (isvector(x) && (length(x) == 2 || length(x) == 3)) ...
        || isa(x, 'LightSource'));
    addOptional(p, 'MeshSize', 1, ...
        @(x) isscalar(x) && isnumeric(x) && (x > 0));
    addOptional(p, 'ZFactor', 1, ...
        @(x) isscalar(x) && isnumeric(x) && (x > 0));
    addParameter(p, 'ShadingMethod', 'clamped', ...
        @(x) ischar(x) && ...
        ismember(x, {'clamped', 'half', 'soft', 'none'}));
    addParameter(p, 'IndexScheme', 'meshgrid', ...
        @(x) ischar(x) && ismember(x, {'meshgrid', 'ndgrid'}));
    parse(p, Z, varargin{:});
    Z = p.Results.HeightMap;
    h = p.Results.MeshSize;
    light = p.Results.Light;
    zFactor = p.Results.ZFactor;
    shadingMethod = lower(p.Results.ShadingMethod);
    indexScheme = lower(p.Results.IndexScheme);

    % Check for demo mode
    if strcmp(Z, 'demo')
        [hillshadeMap, light, flarReflectance, normalMap] ...
            = demo;

        if nargout == 0
            return
        end

        varargout = ...
            {hillshadeMap, light, flarReflectance, normalMap};
        return
    end

    %% Variable formatting
    % If the input is a normal map
    isNormalMap = (ndims(Z) == 3 && size(Z, 3) == 3);

    % Format the light direction to vector form
    if ~isa(light, 'LightSource')
        light = LightSource('parallel', light);
    end

    % Format the indexing scheme
    if strcmp(indexScheme, 'ndgrid')
        Z = Z';
    end

    %% Hillshade calculation
    % Compute the flat reflectance
    flarReflectance = light.Direction(3);

    if isNormalMap
        normalMap = Z;
    else
        normalMap = normal(Z, h, zFactor);
    end

    normalMap = ...
        [ ...
         reshape(normalMap(:, :, 1), [], 1), ...
         reshape(normalMap(:, :, 2), [], 1), ...
         reshape(normalMap(:, :, 3), [], 1) ...
     ];

    hillshadeMap = normalMap * light.Direction(:);
    hillshadeMap = reshape(hillshadeMap, size(Z(:, :, 1)));

    % Apply different shading methods
    switch shadingMethod
        case 'clamped'
            hillshadeMap = max(hillshadeMap, 0);
        case 'half'
            hillshadeMap = (hillshadeMap + 1) / 2;
        case 'soft'
            hillshadeMap = (hillshadeMap + 1) / 2;
            hillshadeMap = hillshadeMap .^ 2;
    end

    %% Output
    if nargout == 0
        figure('Name', 'Hillshade Map')
        surf(Z, hillshadeMap, ...
            'EdgeColor', 'none');
        colormap(gray);
        axis equal
        axis off
        return
    end

    if strcmp(indexScheme, 'ndgrid')
        hillshadeMap = hillshadeMap';
        normalMap = normalMap';
    end

    varargout = {hillshadeMap, light, flarReflectance, normalMap};

end

%% Demo
function [hillshadeMap, light, flarReflectance, normalMap] = demo
    % Create a simple DEM
    [Z, X, Y, h, light] = sampledem2;

    % Calculate the hillshade map
    [hillshadeMap, light, flarReflectance, normalMap] = ...
        hillshade(Z, light, h, 'ShadingMethod', 'clamped');

    % Display the hillshade map
    figName = ['Demo of ', upper(mfilename)];
    figure('Name', figName);

    surf(X, Y, Z, hillshadeMap, ...
        'EdgeColor', 'none');

    colormap(gray);
    axis equal
    axis off
end
