function R = ray(origin, light_direction)
    % RAY creates a ray object
    R.status = 0; % Ray at origin
    R.origin = origin;
    R.originH = R.origin(1:2);
    R.originGrid = round(R.origin);
    if norm(light_direction) ~=0
        R.direction = light_direction / norm(light_direction);
    else
        R.direction = light_direction;
    end
    R.directionH = R.direction(1:2);

    R.position = R.origin;
    R.positionH = R.position(1:2);
    R.positionGrid = R.originGrid;

    R.distance = norm(R.position-R.origin);
    R.distanceH = norm(R.positionH-R.originH);

    R.aux.originH_sft = R.originH + 0.5;
    R.aux.positionH_sft = R.originH + 0.5;
    R.aux.s(1) = 1 / R.directionH(1); % total ray length within 1 horizontal step
    R.aux.s(2) = 1 / R.directionH(2); % total ray length within 1 vertical step
    R.aux.l(1) = (sign(R.directionH(1)) - mod(R.aux.originH_sft(1), sign(R.directionH(1)))) ...
        * R.aux.s(1); % distance (from origin) to the closest vertical gridline
    R.aux.l(2) = (sign(R.directionH(2)) - mod(R.aux.originH_sft(2), sign(R.directionH(2)))) ...
        * R.aux.s(2); % distance (from origin) to the closest horizontal gridline
    R.aux.closestGrid = R.aux.originH_sft;

    if length(R.direction) == 3
        R.aux.slope = R.direction(3) / norm(R.direction(1:2));
    else
        R.aux.slope = 0;
    end
end
