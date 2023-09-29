function R = raymarch(R)
    % RAYMARCH marches the ray to the next grid
    if R.directionH(1) == 0 && R.directionH(2) == 0
        % Vertical ray
    elseif R.directionH(1) == 0 && R.directionH(2) ~= 0
        % Horizontal ray
        if R.directionH(2) > 0
            R.aux.positionH_sft(2) = floor(R.aux.positionH_sft(2) + 1);
        else
            R.aux.positionH_sft(2) = ceil(R.aux.positionH_sft(2) - 1);
        end

    elseif R.directionH(1) ~= 0 && R.directionH(2) == 0
        % Vertical ray
        if R.directionH(1) > 0
            R.aux.positionH_sft(1) = floor(R.aux.positionH_sft(1) + 1);
        else
            R.aux.positionH_sft(1) = ceil(R.aux.positionH_sft(1) - 1);
        end

    else
        % General case
        if R.aux.l(1) > R.aux.l(2)
            xy = 2;
            xy_sup = 1;
        else
            xy = 1;
            xy_sup = 2;
        end

        R.aux.positionH_sft(1) = R.aux.closestGrid(1) + R.aux.l(xy) / R.aux.s(1);
        R.aux.positionH_sft(2) = R.aux.closestGrid(2) + R.aux.l(xy) / R.aux.s(2);
        R.aux.positionH_sft(xy) = round(R.aux.positionH_sft(xy));

        if abs(R.aux.positionH_sft(xy_sup) - round(R.aux.positionH_sft(xy_sup))) < 1e-6
            R.aux.positionH_sft(xy_sup) = round(R.aux.positionH_sft(xy_sup));
            R.aux.closestGrid = R.aux.positionH_sft;
            R.aux.l = [abs(R.aux.s(1)), abs(R.aux.s(2))];
        else
            R.aux.l(xy) = norm(R.aux.positionH_sft - R.aux.closestGrid) + abs(R.aux.s(xy));
        end

    end

    R.positionH = R.aux.positionH_sft - 0.5;
    R.position(1:2) = R.positionH;

    if R.directionH(1) >= 0
        R.positionGrid(1) = ceil(R.aux.positionH_sft(1)) - 1;
    else
        R.positionGrid(1) = floor(R.aux.positionH_sft(1));
    end

    if R.directionH(2) >= 0
        R.positionGrid(2) = ceil(R.aux.positionH_sft(2)) - 1;
    else
        R.positionGrid(2) = floor(R.aux.positionH_sft(2));
    end

end

%% An older version of the function
% function [u,v,lx,ly] = raycast(O,lx,ly,sx,sy)
%     if lx > ly
%         u  = O(1) + ly/sx;
%         v  = round(O(2) + ly/sy);
%         ly = ly + abs(sy);
%     else
%         u  = round(O(1) + lx/sx);
%         v  = O(2) + lx/sy;
%         lx = lx + abs(sx);
%     end
% end
