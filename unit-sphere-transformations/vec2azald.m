function [azimuth, altitude] = vec2azald(vec)
    [azimuth, altitude] = vec2azal(vec);
    azimuth = rad2deg(azimuth);
    altitude = rad2deg(altitude);
end