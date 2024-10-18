%% GRIDEARTH
% Computes the gridded Earth projection for a set of longitude and latitude coordinates.
%
% Syntax
%	[x, y] = gridearth(lon, lat)
%	[x, y] = gridearth(lon, lat, origin, g)
%	[x, y] = gridearth(__, "Name", Value)
%	xy = gridearth(lonlat, __)
%	p = gridearth(p, __)
%	xy = gridearth(geodomain, __)
%
% Input arguments
%	longitudes and latitudes in radians
%		- lon, lat: numeric arrays of longitudes and latitudes
%		- lonlat: N-by-2 (or 2-by-N) numeric array of longitude and latitude coordinates
%		- p: polyshape object
%		- geodomain: GeoDomain object (in the slepian_ulmo package)
%	origin - which longitude (in radians) to centre the projection on
%	g - geometry of the projection, [height, width, base, e]
%		- height: height of the projection
%		- width: width of the projection
%		- base: base of the projection
%		- e: 'eccentricity' of the projection, controlling how the projection is shaped.
%		The default geometry is [2, 4, 2, 0.5].
%	Anchors - whether to add anchors to the projection
%		The default option is true.
%	OutputFormat - format of the output when there is only one output argument
%		- 'polyshape': output as a polyshape object
%		- 'xy': output as a numeric array
%
% Output arguments
%	x, y - numeric arrays of x and y coordinates
%	xy - numeric array of x and y coordinates
%	p - polyshape object
%
% Notes
%	The version with degree inputs is GRIDEARTHD.
%
% See also
%	GRIDEARTHD
%
% Authored by
%	2024-10-18, En-Chi Lee (williameclee@gmail.com)

function varargout = gridearth(varargin)
    %% Initialisation
    % Demos
    if isempty(varargin) || strcmpi(varargin{1}, 'demo')
        fprintf('%s running demo for the gridded Earth projection\n', ...
            upper(mfilename))
        gridearth_demo(mfilename)
        return
    end

    % Parse inputs
    [lon, lat, lonOrigin, g, addAnchors, outputFormat, typeflag] = ...
        parseinputs(varargin{:});

    if addAnchors
        [lon, lat] = addanchors(rad2deg(lon), rad2deg(lat));
        lon = deg2rad(lon);
        lat = deg2rad(lat);
    end

    lon = lon - lonOrigin;

    %% Coordinate transformation
    Y = gridearthy(lat, g);
    X = gridearthx(lon, Y, g);

    %% Collecting output
    if nargout > 2
        error('Too many output arguments')
    end

    % If no output is requested, plot the result
    if nargout == 0
        plotprojections('Gridded Earth Projection', mfilename)
        return
    end

    % Otherwise, return the output
    if nargout == 2
        varargout = {X, Y};
        return
    end

    % Combine X and Y
    XY = [X(:), Y(:)];

    if strcmp(outputFormat, 'polyshape') || (isempty(outputFormat) && typeflag == 1)
        XY = polyshape(XY);
    end

    varargout = {XY};

end

%% Subfunctions
function varargout = parseinputs(varargin)
    % Find longitude and latitude
    [lon, lat, varargin, typeflag] = parselonlatinputs(varargin{:});

    % Find the longitude origin
    lonOriginD = 0;
    p = inputParser;
    addOptional(p, 'LonOrigin', lonOriginD, ...
        @(x) isnumeric(x) || ...
        ((ischar(x) || isstring(x)) && ismember(x, {'center', 'centre', 'c'})));
    addOptional(p, 'Geometry', [2, 4, 2, 0.5], @(x) isnumeric(x) && numel(x) == 4);
    addParameter(p, 'Anchors', true, @(x) islogical(x) || isnumeric(x));
    addParameter(p, 'OutputFormat', '', ...
        @(x) (ischar(x) || isstring(x)) && ismember(x, {'polyshape', 'xy', ''}));
    parse(p, varargin{:});
    lonOrigin = p.Results.LonOrigin;
    g = p.Results.Geometry(:)';
    addAnchors = logical(p.Results.Anchors);
    outputFormat = char(p.Results.OutputFormat);

    % Make sure the geometry is within the bounds
    g(4) = max(min(g(4), 1), 0);

    if ischar(lonOrigin) || isstring(lonOrigin)

        if ismember(lonOrigin, {'center', 'centre', 'c'})
            lonOrigin = (max(lon) + min(lon)) / 2;
        else
            error('Invalid lonOrigin argument %s', lonOrigin)
        end

    end

    varargout = {lon, lat, lonOrigin, g, addAnchors, outputFormat, typeflag};
end
