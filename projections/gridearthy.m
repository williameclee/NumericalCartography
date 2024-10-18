%% GRIDEARTHY
% Computes the y-coordinate of a point on the gridearth projection.
%
% Syntax
%   y = gridearthy(lat, g)
%
% Input arguments
%   lat - latitudes in radians
%   g - geometry of the gridearth projection
%       g = [height, width, base, 'eccentricity']
%
% Output arguments
%   y - y-coordinates on the gridearth projection
%
% See also
%   GRIDEARTH, GRIDEARTHX
%
% Authored by
%   2024-10-18, En-Chi Lee (williameclee@gmail.com)

function y = gridearthy(lat, g)
    hheight = 1; % will break if smaller than 1
    hwidth = g(2) / g(1);
    hbase = g(3) / g(1);
    e = g(4);

    l = -abs(sin(lat));
    a = 3 * e * hwidth - 3 * hbase * e;
    b = 4 * e * hheight * hwidth - 4 * hbase * e * hheight - 4 * hheight * hwidth + 4 * hbase * hheight;
    c = 12 * hheight ^ 3 * hwidth;
    d =- 4 * e * hheight * hwidth * l + 4 * hbase * e * hheight * l - 12 * hheight ^ 3 * hwidth * l + 4 * hheight * hwidth * l + 3 * e * hwidth * l - 4 * hbase * hheight * l - 3 * hbase * e * l;

    f1 = 27 * a * c ^ 2 + 27 * b ^ 2 * d;
    f2 = (12 * a * d - 3 * b * c) .^ 3;
    f3 = (sqrt(f1 .^ 2 - 4 * f2) + f1) .^ (1/3);
    f4 = f3 / (3 * 2 ^ (1/3) * a) + (2 ^ (1/3) * (4 * a * d - b * c)) ./ (a * f3);
    y = -1/2 * sqrt(b ^ 2 / (4 * a ^ 2) + f4) ...
        +1/2 * sqrt(b ^ 2 / (2 * a ^ 2) - (-b ^ 3 / a ^ 3 - (8 * c) / a) ./ (4 * sqrt(b ^ 2 / (4 * a ^ 2) + f4)) - f4) ...
        - b / (4 * a);
    y = min(max(real(y), -1), 0);
    y(lat == -90 | lat == 90) = -1;
    y = y .* (-sign(lat));
    y = y * g(1) / 2;
end
