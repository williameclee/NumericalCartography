%% PARSELONLATINPUTS
% Parses inputs for any longitude and latitude inputs
%
% Syntax
%   [lon, lat, __, type] = dataPath(lon, lat, __)
%   [lon, lat, __, type] = dataPath(lonlat, __)
%   [lon, lat, __, type] = dataPath(p, __)
%   [lon, lat, __, type] = dataPath(geodomain, __)
%
% Input arguments
%   lon, lat - numeric array of longitudes and latitudes
%       The two arrays must have the same size.
%   lonlat - N-by-2 (or 2-by-N) numeric array of longitude and latitude
%       coordinates
%   p - polyshape object
%   geodomain - GeoDomain object
%       This class is found in the slepian_ulmo package.
%   Any other arguments are not processed and passed back to the caller.
%
% Output arguments
%   lon, lat - numeric array of longitudes and latitudes
%   type - which type of input was used
%       0 - lon and lat
%       1 - polyshape
%       2 - GeoDomain
%
% Notes
%   This function was originally part of the slepian_ulmo package.
%
% Authored by
%   2024-08-12, En-Chi Lee (williameclee@gmail.com)
% Last modified by
%   2024-10-18, En-Chi Lee (williameclee@gmail.com)

function varargout = parselonlatinputs(varargin)
    % Find longitude and latitude
    if isnumeric(varargin{1})
        typeflag = 0;

        if length(varargin) >= 2 && isnumeric(varargin{2}) && ...
                isequal(size(varargin{1}), size(varargin{2}))
            % First input is lon, second is lat
            lon = varargin{1};
            lat = varargin{2};

            % Remove lonlat from varargin
            varargin(1:2) = [];

        elseif any(size(varargin{1}) == 2)
            % First input is lonlat
            lonlat = varargin{1};

            if size(lonlat, 1) == 2 && size(lonlat, 2) ~= 2
                lonlat = lonlat';
            end

            lon = lonlat(:, 1);
            lat = lonlat(:, 2);

            varargin(1) = [];

        else
            error('Invalid input argument size')
        end

    elseif isa(varargin{1}, 'polyshape')
        typeflag = 1;
        % First input is a polyshape
        [lon, lat] = boundary(varargin{1});
        varargin(1) = [];
    elseif isa(varargin{1}, 'GeoDomain')
        % This class is found in the slepian_ulmo package
        typeflag = 2;
        % First input is a GeoDomain
        domain = varargin{1};
        lonlat = domain.Lonlat;
        lon = lonlat(:, 1);
        lat = lonlat(:, 2);
        varargin(1) = [];
    else
        error('Invalid input argument type for the first argument')
    end

    typeflag = uint8(typeflag);

    varargout = {lon, lat, varargin, typeflag};
end
