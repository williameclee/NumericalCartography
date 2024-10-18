%% GRIDEARTHX
% Computes the x-coordinate of a point on the gridearth projection.
%
% Syntax
%   x = gridearthx(lon, y, g)
%
% Input arguments
%   lon - longitudes in radians
%   y - y-coordinates on the gridearth projection
%       This is obtained from the gridearthy function.
%   g - geometry of the gridearth projection
%       g = [height, width, base, 'eccentricity']
%
% Output arguments
%   x - x-coordinates on the gridearth projection
%
% See also
%   GRIDEARTH, GRIDEARTHY
%
% Authored by
%   2024-10-18, En-Chi Lee (williameclee@gmail.com)

function x = gridearthx(lon, y, g)
    hheight = g(1) / 2;
    hwidth = g(2) / 2;
    hbase = g(3) / 2;
    e = g(4);

    Lhalfr = hwidth - (e * abs(y / hheight) .^ 3 - (e - 1) * (y / hheight) .^ 2) * (hwidth - hbase);
    x = Lhalfr .* (lon / pi);
end
