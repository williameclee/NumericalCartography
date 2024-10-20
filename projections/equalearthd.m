%% EQUALEARTHD
% Computes the Equal Earth projection of a set of points.
% The Equal Earth projection is an equal-area pseudocylindrical projection,
% making it suitable for assessing surface mass density-related proceses
% (e.g. the GRACE data we are working with).
% This function only differs from EQUALEARTH in that it accepts inputs in
% degrees instead of radians.
%
% Syntax
% 	equalearthd('demo')
%		Runs a demo for the Equal Earth projection.
% 	[XY, X, Y] = equalearthd(lonlatd)
%		Returns the Equal Earth projection of the input longitude-latitude
%       points.
% 	[pd, X, Y] = equalearthd(pd)
%		Returns the Equal Earth projection of the input polyshape object.
% 	[XY, X, Y] = equalearthd(lond, latd)
%		Returns the Equal Earth projection of the input longitude and
%       latitude pairs.
% 	[XY, X, Y] = equalearthd(__, lonOrigind)
%		Returns the Equal Earth projection of the input points with the
%       specified longitude origin.
% 	equalearthd(__)
%		Plots the Equal Earth projection of the input points.
%
% Input arguments
% 	lonlatd - A N-by-2 matrix of longitude-latitude pairs in degrees
% 	pd - A polyshape object, the vertices are in degrees
% 	lond, latd - A vector of longitudes and latitudes in degrees
% 	lonOrigind - The origin of the longitude in degrees
% 		The default value is 0 (i.e. no shifting)
%
% Output arguments
% 	XY - A N-by-2 matrix of the Equal Earth projection of the input points
%	pd - The polyshape object of the Equal Earth projection
% 	X, Y - The X- and Y-coordinates of the Equal Earth projection
%
% Examples
% 	Plot the Equal Earth projection of the world map
% 	>> 	load coastlines
% 	>> 	equalearthd(coastlon, coastlat)
%
% Data source
%	The tranformation are based on the article 'The Equal Earth map
%   projection':
%		Šavrič et al. (2019).
%		doi: 0.1080/13658816.2018.1504949
% Notes
%   The version with radian inputs is EQUALEARTH.
%   This function was originally part of the slepian_ulmo package.
%
% See also
%   EQUALEARTH
%
% Last modified by
%   2024/10/15, williameclee@arizona.edu (@williameclee)

function varargout = equalearthd(varargin)
    %% Initialisation
    % Demos
    if isempty(varargin) || strcmpi(varargin{1}, 'demo')
        fprintf('%s running demo for the equal Earth projection\n', ...
            upper(mfilename))
        equalearth_demo(mfilename)
        return
    end

    % Parse inputs
    [lond, latd, lonOrigind, addAnchors, outputFormat] = parseinputs(varargin{:});

    if addAnchors
        [lond, latd] = addanchors(lond, latd);
    end

    %% Coordinate transformation
    [XY, X, Y] = equalearth(deg2rad(lond), deg2rad(latd), deg2rad(lonOrigind));

    %% Collecting output
    if nargout == 0
        figName = sprintf('Equal Earth Projection (%s)', upper(mfilename));
        figure(999)
        set(gcf, 'Name', figName, 'NumberTitle', 'off')

        plot(X, Y, 'k', 'LineWidth', 0.5)

        axis equal
        xlim([min(X) - 0.1, max(X) + 0.1])
        ylim([min(Y) - 0.1, max(Y) + 0.1])

        return
    end

    if (isa(varargin{1}, 'polyshape') && (nargout == 1 || nargout == 3)) || strcmp(outputFormat, 'polyshape')

        if ~isa(XY, 'polyshape')
            XY = polyshape(XY);
        end

        varargout = {XY, X, Y};
    elseif nargout == 1
        varargout = {XY};
    elseif nargout == 2
        varargout = {X, Y};
    else
        varargout = {XY, X, Y};
    end

end

%% Subfunctions
function varargout = parseinputs(varargin)
    % Find longitude and latitude
    [lon, lat, varargin] = parselonlatinputs(varargin{:});
    % Find the longitude origin
    lonOriginD = 180 * double(min(lon) >= 0);
    p = inputParser;
    addOptional(p, 'LonOrigin', lonOriginD, ...
        @(x) isnumeric(x) || ...
        ((ischar(x) || isstring(x)) && ismember(x, {'center', 'centre', 'c'})));
    addOptional(p, 'Anchors', false, @(x)islogical(x) || isnumeric(x));
    addParameter(p, 'OutputFormat', '', @(x) ischar(x) || isstring(x));
    parse(p, varargin{:});
    lonOrigin = p.Results.LonOrigin;
    addAnchors = logical(p.Results.Anchors);
    outputFormat = char(p.Results.OutputFormat);

    if ischar(lonOrigin) || isstring(lonOrigin)

        if ismember(lonOrigin, {'center', 'centre', 'c'})
            lonOrigin = (max(lon) + min(lon)) / 2;
        else
            error('Invalid lonOrigin argument %s', lonOrigin)
        end

    end

    varargout = {lon, lat, lonOrigin, addAnchors, outputFormat};
end
