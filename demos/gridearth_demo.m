%% GRIDEARTH_DEMO
% This is a demo for GRIDEARTH and GRIDEARTHD.
%
% Last modified by
%   2024-10-18, En-Chi Lee (williameclee@gmail.com)

function gridearth_demo(varargin)
    %% Generating data
    lonOrigin = 0;

    % The shape of the oceans
    load coastlines %#ok<LOAD>

    % Flatten if possible
    try
        [coastlat, coastlon] = flatearthpoly(coastlat, coastlon, lonOrigin); %#ok<NODEF>
    catch
    end

    coastLonlat = polyshape(coastlon, coastlat);

    % Load test circles
    ccenterLon = -180:30:150;
    ccenterLat = -60:20:60;
    cRadius = 5;
    [ccenterLonn, ccenterLatt] = meshgrid(ccenterLon, ccenterLat);
    ccenterLonn = reshape([ccenterLonn(:)'; nan([1, numel(ccenterLonn)])], [], 1);
    ccenterLatt = reshape([ccenterLatt(:)'; nan([1, numel(ccenterLatt)])], [], 1);
    [circleslat, circleslon] = bufferm(ccenterLatt, ccenterLonn, cRadius, 'out');
    circlesLonlat = polyshape(circleslon, circleslat);

    % The bounding box
    bbox = polyshape([[-180; 180; 180; -180; -180] + lonOrigin, [-90; -90; 90; 90; -90]]);

    %% Projection
    coastXY = gridearthd(coastLonlat, lonOrigin);
    circlesXY = gridearthd(circlesLonlat, lonOrigin);
    bboxXY = gridearthd(bbox, lonOrigin, "OutputFormat", 'polyshape');

    %% Plotting
    figName = 'Equal Earth Projection';

    if nargin > 0
        funName = varargin{1};
        figName = sprintf('%s (%s)', figName, upper(funName));
    end

    figure(999)
    set(gcf, 'Name', figName, 'NumberTitle', 'off')
    clf

    % Plot the original data
    subplot(2, 1, 1)
    title('Unprojected coastline')

    hold on
    plot(coastLonlat)
    plot(circlesLonlat)
    plot(bbox, "FaceColor", 'none')
    hold off

    formatlonticks
    formatlatticks
    axis equal
    xlim([min(bbox.Vertices(:, 1)), max(bbox.Vertices(:, 1))])
    ylim([min(bbox.Vertices(:, 2)), max(bbox.Vertices(:, 2))])
    set(gca, "Box", 'on', "Layer", 'top', "Color", 'none')

    % Plot the Equal Earth projection
    subplot(2, 1, 2)
    title('Gridded Earth projection')

    hold on
    plot(coastXY)
    plot(circlesXY)
    plot(bboxXY, "FaceColor", 'none')
    hold off

    axis equal
    xlim([min(bboxXY.Vertices(:, 1)), max(bboxXY.Vertices(:, 1))])
    ylim([min(bboxXY.Vertices(:, 2)), max(bboxXY.Vertices(:, 2))])
    set(gca, "XAxisLocation", 'origin', "YAxisLocation", 'origin', ...
        "TickDir", 'both', "Layer", 'top', "Color", 'none')
end
