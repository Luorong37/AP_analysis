function [rois, traces] = select_ROI(movie, num_rows, num_cols, t, colors, mask, map)
% Initialize variables
rois = {};
traces = [];
its_2D = reshape(movie, num_cols, num_rows, []);
% Set current figure to full screen and add keypress callback
fig = gcf;
set(fig,'Position',get(0,'Screensize'));
set(gcf, 'KeyPressFcn', @(src, event) set_space_pressed(event, fig));

% Display sensitivity map if provided
if ~isempty(map)
    subplot(1,2,1);
    map_axe = gca;
    imagesc(map);
    colorbar;
    title(sprintf('Sensitivity MAP\n\nPress SPACE to continue\nPress ENTER to end'));
    hold on;
    axis image;
end

% Display fluorescent image
subplot(2,2,2);
image_axe = gca;
im_adj = imadjust(uint16(mean(its_2D, 3)));
imshow(im_adj);
hold on;
title('Fluorescent Image');

% Plot Trace
subplot(2,2,4);
trace_axe = gca;
hold on;
xlabel('Time (s)');
ylabel('Intensity');

% Mask provided, overlay ROI and traces
if ~isempty(mask)
    roi_amount = size(mask, 2);
    rois = mask;
    for i = 1:roi_amount
        trace = mean(movie(mask{i}, :));
        traces = [traces trace'];
        boundary = bwboundaries(mask{i});
        plot(boundary{1}(:, 2), boundary{1}(:, 1), 'Color', colors(mod(i - 1, length(colors)) + 1, :), 'LineWidth', 2, 'Parent', map_axe);
        plot(t, trace, 'Color', colors(mod(i - 1, length(colors)) + 1, :), 'Parent', trace_axe);
    end

else
    % ROI selection process
    while true
        if ~isempty(map)
            roi_mask = drawpolygon('Color', colors(mod(length(rois), length(colors)) + 1, :), 'LineWidth', 1, 'Parent', map_axe);
        else
            roi_mask = drawpolygon('Color', colors(mod(length(rois), length(colors)) + 1, :), 'LineWidth', 1, 'Parent', image_axe);
        end
        mask = poly2mask(roi_mask.Position(:, 1), roi_mask.Position(:, 2), size(im_adj, 1), size(im_adj, 2));
        rois{end + 1} = mask;
        trace = mean(movie(mask, :));
        traces = [traces trace'];

        boundary = bwboundaries(mask);
        plot(boundary{1}(:, 2), boundary{1}(:, 1), 'Color', colors(mod(length(rois) - 1, length(colors)) + 1, :), 'LineWidth', 1, 'Parent', image_axe);
        plot(t, trace, 'Color', colors(mod(length(rois) - 1, length(colors)) + 1, :), 'Parent', trace_axe);

        % Wait for user input
        fig.UserData = [];
        waitfor(fig, 'UserData');
        if strcmp(fig.UserData, 'stop')
            break;
        elseif strcmp(fig.UserData, 'spacePressed')
            continue;
        end
    end

end
fprintf('Finished ROI selection\n');
end

function set_space_pressed(event, fig)
if strcmp(event.Key, 'space')
    fig.UserData = 'spacePressed';
elseif strcmp(event.Key, 'return')
    fig.UserData = 'stop';
end
end
