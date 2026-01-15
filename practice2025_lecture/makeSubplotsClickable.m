%--------------------------------------------------------------------------
% makeSubplotsClickable.m
% Author: Yutaka Fuchino
% Created: 2025-02-26
%
% Description:
% This script enables click-to-enlarge functionality for all subplots 
% in the current figure. It supports various plot types, including 
% plot(), imagesc(), and patch().
%
% Features:
% - Click on a subplot to open an enlarged version in a new figure.
% - Supports plot(), imagesc(), patch(), and more.
% - Retains original axis settings (XLim, YLim, XDir, YDir).
% - Preserves original colormap and colorbar.
% - Only makes the axis transparent for patch subplots.
%
% Usage:
% - Run this script after creating subplots to enable click-to-enlarge.
%
% License: MIT License
% See LICENSE file for details.
%--------------------------------------------------------------------------

function makeSubplotsClickable()
    fig = gcf; % Get the current figure
    axs = findall(fig, 'Type', 'axes'); % Get all subplot axes

    for i = 1:length(axs)
        ax = axs(i);
        
        % Make the axis transparent only for patch subplots
        patches = findall(ax, 'Type', 'Patch');
        if ~isempty(patches)
            set(ax, 'XColor', 'none', 'YColor', 'none', 'Box', 'off', 'Color', 'none');
        end
        
        % Disable HitTest for imagesc plots to allow axis click detection
        images = findall(ax, 'Type', 'Image');
        for j = 1:length(images)
            set(images(j), 'HitTest', 'off');
        end
        
        % Apply click event to the subplot
        set(ax, 'ButtonDownFcn', @(src, event) enlarge_subplot(src));
    end
end

% Function to enlarge the subplot when clicked
function enlarge_subplot(ax)
    % Create new figure for enlarged subplot
    newFig = figure;
    newAx = axes(newFig);
    hold(newAx, 'on'); % Allow multiple plots in the new figure

    % Copy title
    if ~isempty(ax.Title.String)
        title(newAx, ax.Title.String);
    end
    xlabel(newAx, ax.XLabel.String);
    ylabel(newAx, ax.YLabel.String);
    
    % Retain original axis settings
    set(newAx, 'XLim', ax.XLim, 'YLim', ax.YLim, ...
               'XDir', ax.XDir, 'YDir', ax.YDir, ...
               'XColor', ax.XColor, 'YColor', ax.YColor, ...
               'Box', ax.Box, 'Color', ax.Color);

    % Copy all line plots (e.g., plot)
    lines = findall(ax, 'Type', 'Line');
    for i = 1:length(lines)
        plot(newAx, lines(i).XData, lines(i).YData, ...
             'Color', lines(i).Color, 'LineStyle', lines(i).LineStyle, ...
             'LineWidth', lines(i).LineWidth, 'Marker', lines(i).Marker);
    end

    % Copy all patches (for patch, etc.)
    patches = findall(ax, 'Type', 'Patch');
    for i = 1:length(patches)
        patch('XData', patches(i).XData, 'YData', patches(i).YData, ...
              'FaceColor', patches(i).FaceColor, 'EdgeColor', patches(i).EdgeColor, ...
              'FaceAlpha', patches(i).FaceAlpha, 'EdgeAlpha', patches(i).EdgeAlpha, ...
              'Parent', newAx);
    end

    % Copy images (for imagesc, imshow, etc.)
    images = findall(ax, 'Type', 'Image');
    for i = 1:length(images)
        xData = images(i).XData;
        yData = images(i).YData;
        cData = images(i).CData;

        % Adjust XData and YData if they are in default [1, size(CData)]
        if isequal(xData, [1 size(cData,2)]) && isequal(yData, [1 size(cData,1)])
            xData = linspace(ax.XLim(1), ax.XLim(2), size(cData,2));
            yData = linspace(ax.YLim(1), ax.YLim(2), size(cData,1));
        end
        
        imagesc('XData', xData, 'YData', yData, 'CData', cData, ...
                'Parent', newAx);
    end

    % Copy colormap
    colormap(newFig, colormap(ax));

    % Copy colorbar if exists
    colorbarObj = findall(ax.Parent, 'Type', 'Colorbar');
    if ~isempty(colorbarObj)
        colorbar(newAx);
    end
end
