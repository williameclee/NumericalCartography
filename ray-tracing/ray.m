function R = ray(origin, light_direction)
    % RAY creates a ray object
    R.origin = origin;
    R.originH = R.origin(1:2);
    R.originGrid = round(R.origin);
    R.direction = light_direction;
    R.directionH = R.direction(1:2);

    R.position = R.origin;
    R.positionH = R.position(1:2);
    R.positionGrid = R.originGrid;

    R.aux.originH_sft = R.originH + 0.5;
    R.aux.positionH_sft = R.originH + 0.5;
    R.aux.s(1) = 1 / R.directionH(1); % total ray length within 1 horizontal step
    R.aux.s(2) = 1 / R.directionH(2); % total ray length within 1 vertical step
    R.aux.l(1) = (sign(R.directionH(1)) - mod(R.aux.originH_sft(1), sign(R.directionH(1)))) ...
        * R.aux.s(1); % distance (from origin) to the closest vertical gridline
    R.aux.l(2) = (sign(R.directionH(2)) - mod(R.aux.originH_sft(2), sign(R.directionH(2)))) ...
        * R.aux.s(2); % distance (from origin) to the closest horizontal gridline
    R.aux.closestGrid = R.aux.originH_sft;
end
