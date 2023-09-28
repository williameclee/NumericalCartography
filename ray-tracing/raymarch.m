function ray = raymarch(ray)
    % RAYMARCH marches the ray to the next grid
    if ray.aux.lx > ray.aux.ly
        ray.aux.u  = ray.aux.originH_sft(1) + ray.aux.ly/ray.aux.sx;
        ray.aux.v  = round(ray.aux.originH_sft(2) + ray.aux.ly/ray.aux.sy);
        ray.aux.ly = ray.aux.ly + abs(ray.aux.sy);
    else
        ray.aux.u  = round(ray.aux.originH_sft(1) + ray.aux.lx/ray.aux.sx);
        ray.aux.v  = ray.aux.originH_sft(2) + ray.aux.lx/ray.aux.sy;
        ray.aux.lx = ray.aux.lx + abs(ray.aux.sx);
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