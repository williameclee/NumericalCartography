%% HILLSHADE
% Calculates the hillshade of a DEM.
%
% Syntax
%   hillshadeMap = hillshade(Z)
%   hillshadeMap = hillshade(Z, [azimuth, altitude])
%   hillshadeMap = hillshade(Z, [Lx, Ly, Lz])
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
%   hillshadeMap = hillshade(Z, [azimuth, altitude])
%       calculates the hillshade of the DEM Z with the light azimuth and
%       altitude specified in the 1-by-2 vector [azimuth, altitude] in 
%       degrees.
%   hillshadeMap = hillshade(Z, [Lx, Ly, Lz])
%       calculates the hillshade of the DEM Z with the light direction
%       specified in the 1-by-3 vector L.
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
%   Z - DEM matrix
%   L - Light direction vector
%       - 1-by-2 vector [azimuth, altitude] in degrees
%       - 1-by-3 vector [Lx, Ly, Lz] of the (invert) light direction
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
%   L - Light direction vector
%   r - Flat reflectance
%       The flat reflectance is the reflectance of a flat surface facing
%       upwards.
%   N - Normal map
%
% Authored by
%   En-Chi Lee <williameclee@gmail.com>, 2024-06-23
% Last modified by
%   En-Chi Lee <williameclee@gmail.com>, 2024-06-28

function varargout = hillshade(Z, varargin)
    %% Initialisation and preallocation
    % Make sure subfunctions in the path
    addpath(fullfile(fileparts(mfilename('fullpath')), ...
    'unit-sphere-transformations'))

    % Input parser
    p = inputParser;
    addRequired(p, 'HeightMap', ...
        @(x) (isnumeric(x) && ismatrix(x)) || ...
        (ischar(x) && strcmp(x, 'demo')));
    addOptional(p, 'LightVector', [-1, 1, 1], ...
        @(x) isvector(x) && (length(x) == 2 || length(x) == 3));
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
    lightVector = p.Results.LightVector;
    zFactor = p.Results.ZFactor;
    shadingMethod = lower(p.Results.ShadingMethod);
    indexScheme = lower(p.Results.IndexScheme);

    % Check for demo mode
    if strcmp(Z, 'demo')
        [hillshadeMap, lightVector, flarReflectance, normalMap] ...
            = demo;

        if nargout == 0
            return
        end

        varargout = ...
            {hillshadeMap, lightVector, flarReflectance, normalMap};
        return
    end

    % Preallocation
    normalMap = zeros([size(Z), 3]);

    %% Variable formatting
    % Rescale the DEM
    Z = Z * zFactor / h;

    if strcmp(indexScheme, 'ndgrid')
        Z = Z';
    end

    % Format the light direction to vector form
    if length(lightVector) == 2
        lightVector = azald2vec(lightVector(1), lightVector(2));
    else
        lightVector = lightVector / norm(lightVector(:));
    end

    % Make sure the light direction is pointing upwards
    % (i.e. towards the light source)
    lightVector = lightVector * sign(lightVector(3));

    % Compute the flat reflectance
    flarReflectance = lightVector(3);

    %% Hillshade calculation
    [normalMap(:, :, 1), normalMap(:, :, 2), normalMap(:, :, 3)] = ...
        surfnorm(Z);
    normalMap = [ ...
                     reshape(normalMap(:, :, 1), [], 1), ...
                     reshape(normalMap(:, :, 2), [], 1), ...
                     reshape(normalMap(:, :, 3), [], 1) ...
                 ];
    hillshadeMap = normalMap * lightVector(:);
    hillshadeMap = reshape(hillshadeMap, size(Z));

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
        figure
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

    lightVector = lightVector(:)';
    varargout = {hillshadeMap, lightVector, flarReflectance, normalMap};

end

%% Demo
function [hillshadeMap, light, flarReflectance, normalMap] = demo
    z = peaks(2 ^ 9) * 20;
    light = [-1, 1, 0.5];
    [hillshadeMap, light, flarReflectance, normalMap] = ...
        hillshade(z, light, 'ShadingMethod', 'clamped');

    figure('Name', ['Demo of ', upper(mfilename)])
    surf(z, hillshadeMap, ...
        'EdgeColor', 'none');
    colormap(gray);
    axis equal
    axis off
end
