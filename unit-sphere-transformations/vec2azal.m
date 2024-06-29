function [azimuth, altitude] = vec2azal(vec)
    vec = vec / norm(vec(:));
    azimuth = atan2(vec(1), vec(2));
    altitude = asin(vec(3));
end