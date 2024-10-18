%% PLOTPROJECTIONS
% Creates quick and dirty plot of the projections
%
% Authored by
%	2024-10-18, En-Chi Lee (williameclee@gmail.com)

function plotprojections(prjname, funname)
    figName = sprintf('%s (%s)', prjname, upper(funname));
    figure(999)
    set(gcf, 'Name', figName, 'NumberTitle', 'off')

    plot(X, Y, 'k', 'LineWidth', 0.5)

    axis equal
    xlim([min(X) - 0.1, max(X) + 0.1])
    ylim([min(Y) - 0.1, max(Y) + 0.1])
end
