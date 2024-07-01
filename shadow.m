%% SHADOW
% Calculates shadow casted on a DEM.
%
% Syntax
%   shadowMap = shadow(Z)
%   shadowMap = shadow(Z, L)
%   shadowMap = shadow(Z, __, h)
%   shadowMap = shadow(Z, __, VE)
%   shadowMap = shadow(Z, __, Name, Value)
%   shadowMap = shadow('demo')
%   [shadowMap, shadowDistMap] = shadow(__)
%
% Description
%   shadowMap = shadow(Z)
%       calculates what pixels are shadowed in the DEM Z with the default
%       light azimuth and altitude of 45 degrees.
%   shadowMap = shadow(Z, L)
%       calculates what pixels are shadowed in the DEM Z with the light
%       source.
%   shadowMap = shadow(Z, __, h)
%       calculates what pixels are shadowed in the DEM Z with the mesh size
%       h. The DEM is assumed to have the same mesh size in both
%       dimensions.
%   shadowMap = shadow(Z, __, VE)
%       calculates what pixels are shadowed in the DEM Z with the vertical
%       exaggeration factor VE.
%   shadowMap = shadow('demo')
%       displays the shadow map of the peaks function.
%   shadow(__)
%       displays the shadow map shadowMap.
%   [shadowMap, shadowDistMap] = shadow(__)
%       returns the shadow map shadowMap and the distance from the shadow
%       shadowDistMap.
%
% Input Arguments
%   Z - DEM matrix
%   L - Light source
%       - LightSource object.
%       - 1-by-2 vector [azimuth, altitude] in degrees.
%       - 1-by-3 vector [Lx, Ly, Lz] of the (invert) light direction.
%       The default value is [-1, 1, 1].
%   h - Mesh size
%       The default value is 1.
%   VE - Vertical exaggeration factor
%       The default value is 1 (no exaggeration).
%   IndexScheme - Indexing scheme
%       The default value is 'meshgrid'. The other option is 'ndgrid'.
%   HillShadeMap - Hillshade map
%       If provided, the function will use the hillshade map to speed up
% 		the calculation.
%
% Output Arguments
%   shadowMap - Shadow map
%   	A logical matrix of the same size as Z. The value is true if the
%       pixel is shadowed.
%   shadowDistMap - Distance from the shadow casting pixel
%       A matrix of the same size as Z. The value is the distance from the
%       shadow in the direction of the light. The value is NaN if the pixel
%       is not shadowed.
%
% See also
%   HILLSHADE, LIGHTSOURCE
%
% Authored by
%   En-Chi Lee <williameclee@gmail.com>, 2024-06-28
% Last modified by
%   En-Chi Lee <williameclee@gmail.com>, 2024-07-01

function varargout = shadow(Z, varargin)
    %% Initilisation
    % Make sure subfunctions in the path
    addpath(fullfile(fileparts(mfilename('fullpath')), ...
    'aux'));

    % Parse the input
    p = inputParser;
    addRequired(p, 'DEM', ...
        @(x) (isnumeric(x) && ismatrix(x)) || ...
        (ischar(x) && strcmp(x, 'demo')));
    addOptional(p, 'Light', [-1, 1, 1], ...
        @(x) (isvector(x) && (length(x) == 2 || length(x) == 3)) ...
        || isa(x, 'LightSource'));
    addOptional(p, 'MeshSize', 1, ...
        @(x) isscalar(x) && isnumeric(x) && (x > 0));
    addOptional(p, 'ZFactor', 1, ...
        @(x) isscalar(x) && isnumeric(x) && (x > 0));
    addParameter(p, 'IndexScheme', 'meshgrid', ...
        @(x) ischar(x) && ismember(x, {'meshgrid', 'ndgrid'}));
    addParameter(p, 'HillShadeMap', [], ...
        @(x) isnumeric(x) && ismatrix(x));
    parse(p, Z, varargin{:});
    Z = p.Results.DEM;
    light = p.Results.Light;
    meshSize = p.Results.MeshSize;
    zFactor = p.Results.ZFactor;
    indexScheme = p.Results.IndexScheme;
    hillshadeMap = p.Results.HillShadeMap;

    % More input checking
    if ~isempty(hillshadeMap) && ~isequal(size(hillshadeMap), size(Z))
        error('The size of the hillshade map must be the same as the DEM.')
    end

    % Check for demo mode
    if strcmp(Z, 'demo')
        [isShadowed, distFromShadow] = demo;

        if nargout == 0

            return
        end

        varargout = {isShadowed, distFromShadow};
        return

    end

    %% Variable formatting
    % Format the indexing scheme
    if strcmp(indexScheme, 'ndgrid')
        Z = Z';
    end

    % Format the light direction to vector form
    if ~isa(light, 'LightSource')
        light = LightSource('parallel', light);
    end

    Z = Z * zFactor / meshSize;

    %% Preallocation
    mapSize = size(Z);
    nPixel = prod(mapSize);

    isShadowed = false([nPixel, 1]);
    isCaster = false([nPixel, 1]);

    if nargout >= 2
        distFromShadow = nan([nPixel, 1]);
        zOfCaster = nan([nPixel, 1]);
    end

    if isempty(hillshadeMap)
        hillshadeMap = zeros([nPixel, 1]);
        hillshadeMap(:) = hillshade(Z, light.Direction);
    else
        hillshadeMap = hillshadeMap(:);
        hillshadeMap = max(0, hillshadeMap);
    end

    %% Rescale and rotate coordinate to light direction
    [I, J] = meshgrid(1:mapSize(2), 1:mapSize(1));
    [I, J, Z] = rotateandtiltmap(I, J, Z, light);

    %% Find shadow-casting pixels
    isShadowed(:) = hillshadeMap(:) <= 0;

    [isCaster(:), nCaster] = findshadowcasters(light, Z);
    IJZcasters = [I(isCaster), J(isCaster), Z(isCaster)];

    [I, J, Z, isShadowed, isCaster, sortId] = ...
        sortpixels(I, J, Z, isShadowed, isCaster);

    % Make an educated guess on where are the pixels that might be shadowed
    beamWidth = dot([1/2, 1/2] .* sign(light.HorizontalDirection), light.HorizontalDirection);

    % Normalise the coordinates to the range of [0, nBin]
    [J, IJZcasters, beamWidth, nBin] = ...
        normalisetobin(J, IJZcasters, beamWidth);

    %% Iterating through shadow casters
    binBounds = findbinbounds(J, nBin);

    for iCaster = 1:nCaster

        if isShadowed(isCaster(iCaster))
            continue
        end

        IJZcaster = IJZcasters(iCaster, :);

        binBound = binBounds(floor(IJZcaster(2)) + 1, :);
        ijCasted = binBound(1):binBound(2);
        % Find pixels that can be shadowed by the current shadow caster
        ijCasted = ijCasted(abs(J(ijCasted) - IJZcaster(2)) <= beamWidth);
        ijCasted = ijCasted(I(ijCasted) > IJZcaster(1));
        % Find if the pixels are actually shadowed
        ijCasted = ijCasted(Z(ijCasted) < IJZcaster(3));
        isShadowed(ijCasted) = true;

        if nargout < 2
            continue
        end

        % First-order approximation: ignore the y difference
        zOfCaster(ijCasted) = max(IJZcaster(3), zOfCaster(ijCasted));
        distFromShadow(ijCasted(IJZcaster(3) == zOfCaster(ijCasted))) = ...
            IJZcasters(iCaster, 1);
    end

    % Unsort the pixels
    if nargout >= 2
        distFromShadow = (I - distFromShadow) * meshSize;
        [isShadowed, distFromShadow] = ...
            unsortpixels(isShadowed, distFromShadow, sortId, mapSize);
    else
        isShadowed = unsortpixels(isShadowed, [], sortId, mapSize);
    end

    %% Output
    if nargout == 0
        figure('Name', 'Shadow Map')
        surf(Z, isShadowed, ...
            'EdgeColor', 'none');
        colormap(gray);
        axis equal
        axis off
        return
    end

    if strcmp(indexScheme, 'ndgrid')
        isShadowed = isShadowed';

        if nargout >= 2
            distFromShadow = distFromShadow';
        end

    end

    if nargout == 1
        varargout = {isShadowed};
    else
        varargout = {isShadowed, distFromShadow};
    end

end

%% Subfunctions

% Rotate and tilt the map to the light direction
function [Irotated, Jrotated, Ztilted] = ...
        rotateandtiltmap(I, J, Z, light)
    % Make sure the transformation functions are found in the path
    addpath(fullfile(fileparts(mfilename('fullpath')), ...
    'unit-sphere-transformations'));

    % Rotate the coordinate system
    angOfRotation = deg2rad(light.Azimuth) + pi / 2;
    RotationMatrix = ...
        [cos(angOfRotation), -sin(angOfRotation); ...
         sin(angOfRotation), cos(angOfRotation)];
    IJ = [I(:), J(:)];
    IJrotated = (RotationMatrix * IJ')';
    Irotated = reshape(IJrotated(:, 1), size(I));
    Jrotated = reshape(IJrotated(:, 2), size(J));

    % Tilt the map
    Ztilted = Z + Irotated * tand(light.Altitude);
end

% Find which pixels are potential shadow casters
function [isCaster, nCaster] = findshadowcasters(light, Z)
    %% Initialisation and preallocation
    light = light.HorizontalDirection;

    % Preallocate the output
    slopeMap = zeros(size(Z));
    isCaster = false(size(Z));
    gradientMapI = zeros([numel(Z), 1]);
    gradientMapJ = zeros([numel(Z), 1]);

    %% Find the slope of Z in the direction of the light direction
    % Calculate the gradient of Z
    [gradientMapI(:), gradientMapJ(:)] = gradient(Z);
    % Project the gradient to the light direction
    slopeMap(:) = [gradientMapI, gradientMapJ] * -light(:);
    % Find the sign of the slope
    slopeMap = sign(slopeMap);

    clear gradientMapI gradientMapJ

    %% Find the edges of the zero-slope region
    % TODO: Is there an economic way to eliminate redundant points?
    if light(1) > 0 % light from left
        isCaster(:, 2:end) = diff(slopeMap, 1, 2) > 0 ...
            | isCaster(:, 2:end);
    elseif light(1) < 0 % light from right
        isCaster(:, 1:end - 1) = diff(slopeMap, 1, 2) < 0 ...
            | isCaster(:, 1:end - 1);
    end

    if light(2) > 0 % light from top
        isCaster(2:end, :) = diff(slopeMap, 1, 1) > 0 ...
            | isCaster(2:end, :);
    elseif light(2) < 0 % light from bottom
        isCaster(1:end - 1, :) = diff(slopeMap, 1, 1) < 0 ...
            | isCaster(1:end - 1, :);
    end

    if light(1) > 0 && light(2) > 0 % light from quadrant 1
        isCaster(2:end, 2:end) = ...
            slopeMap(2:end, 2:end) > slopeMap(1:end - 1, 1:end - 1) ...
            | isCaster(2:end, 2:end);
    elseif light(1) < 0 && light(2) > 0 % light from quadrant 2
        isCaster(2:end, 1:end - 1) = ...
            slopeMap(2:end, 1:end - 1) > slopeMap(1:end - 1, 2:end) ...
            | isCaster(2:end, 1:end - 1);
    elseif light(1) < 0 && light(2) < 0 % light from quadrant 3
        isCaster(1:end - 1, 1:end - 1) = ...
            slopeMap(1:end - 1, 1:end - 1) > slopeMap(2:end, 2:end) ...
            | isCaster(1:end - 1, 1:end - 1);
    elseif light(1) > 0 && light(2) < 0 % light from quadrant 4
        isCaster(1:end - 1, 2:end) = ...
            slopeMap(1:end - 1, 2:end) > slopeMap(2:end, 1:end - 1) ...
            | isCaster(1:end - 1, 2:end);
    end

    %% Add the boundaries
    if light(1) > 0
        isCaster(:, end) = true;
    elseif light(1) < 0
        isCaster(:, 1) = true;
    end

    if light(2) > 0
        isCaster(end, :) = true;
    elseif light(2) < 0
        isCaster(1, :) = true;
    end

    nCaster = sum(isCaster(:));

end

% Sort the pixels to make the search more efficient
function [I, J, Z, shadow, caster, sortId] = ...
        sortpixels(I, J, Z, shadow, caster)
    [IJZ, sortId] = sortrows( ...
        [I(:), J(:), Z(:), shadow(:)], [2, 1]);

    I = IJZ(:, 1);
    J = IJZ(:, 2);
    Z = IJZ(:, 3);

    shadow = logical(IJZ(:, 4));
    caster = caster(sortId);
end

% Reverse the sorting by sortpixels
function [shadow, dist] = ...
        unsortpixels(shadow, dist, sortId, mapSize)
    shadow(sortId) = shadow;
    shadow = reshape(shadow, mapSize);

    if ~isempty(dist)
        dist(sortId) = dist;
        dist = reshape(dist, mapSize);
    else
        dist = [];
    end

end

% Normalise the J coordinate to the range of [0, nBin]
function [J, IJZ, beamWidth, nBin] = ...
        normalisetobin(J, IJZ, beamWidth)
    % Randomly chosen numbers
    nBinScale = 2 ^ 3;
    nBinMax = 2 ^ 12;

    % Normalise the coordinates to the range of [0, nBin]
    Jlim = [min(J), max(J)];
    nBin = min(floor(diff(Jlim) / beamWidth / nBinScale), nBinMax);

    scalingFactor = nBin / diff(Jlim);

    J = (J - Jlim(1)) * scalingFactor;
    J = min(nBin, max(0, J));

    IJZ(:, 2) = (IJZ(:, 2) - Jlim(1)) * scalingFactor;
    IJZ(:, 2) = min(nBin, max(0, IJZ(:, 2)));

    beamWidth = beamWidth * scalingFactor;
end

% Make educated guess on the pixels that might be shadowed
function binBounds = findbinbounds(J, nBin)
    binBounds = zeros([nBin + 1, 2]);

    J = uint16(floor(J));

    binBounds(end, 2) = length(J);

    % Looping backwards to make search range smaller
    for iBinEdge = nBin - 1:-1:0
        binBounds(iBinEdge + 1, 2) = ...
            find(J(1:binBounds(iBinEdge + 2, 2)) <= iBinEdge + 2, ...
            1, 'last');

        if iBinEdge <= 2
            binBounds(iBinEdge + 1, 1) = find( ...
                J >= iBinEdge - 1, ...
                1, 'first');
        end

    end

    binBounds(4:end, 1) = binBounds(1:end - 3, 2) + 1;
end

%% Demo
function [isShadowed, distFromShadow] = demo
    % Create a simple DEM
    [Z, X, Y, h, light] = sampledem2;

    % Calculate the hillshade and shadow map
    [hillshadeMap, ~, ~] = ...
        hillshade(Z, light, h, 'ShadingMethod', 'clamped');
    [isShadowed, distFromShadow] = shadow(Z, light, h, ...
        'HillShadeMap', hillshadeMap);

    % Blur the shadow
    isShadowed = imgaussfilt(double(isShadowed), 3);

    % Display the shadow map
    figName = ['Demo of ', upper(mfilename)];
    figure('Name', figName);

    % Blend the shadow map with the hillshade map
    surf(X, Y, Z, hillshadeMap .* (1 - isShadowed), ...
        'EdgeColor', 'none');

    colormap(gray);
    axis equal
    axis tight
    axis off
end
