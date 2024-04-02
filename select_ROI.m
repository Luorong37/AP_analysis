function [rois, traces] = select_ROI(movie, nrows, ncols, t, colors, mask, map)
% 
% ----------Write by Liu-Yang Luorong and ChatGPT----------
% ----------POWERED by Zoulab in Peking University----------
% Date: 23.11.16
% MATLAB Version: R2022b
% SELECT_ROI Allows for interactive selection of Regions of Interest (ROIs) in a movie data set.
%
%   This function is designed for interactive selection and visualization of ROIs in movie data. 
%   It displays the movie and an optional sensitivity map, allowing the user to draw polygons 
%   representing ROIs. The function then calculates and returns the mean intensity trace for each ROI.
%
%   Syntax:
%   [rois, traces] = select_ROI(movie, num_rows, num_cols, t, colors, mask, map)
%
%   Parameters:
%   movie - A 2D array representing the movie, with dimensions [ncols*nrows, nframes].
%   nrows - Number of rows in the movie.
%   ncols - Number of columns in the movie.
%   t - Time vector corresponding to the movie frames.
%   colors - A matrix of RGB values for differentiating multiple ROIs.
%   mask - (Optional) A predefined set of masks for ROIs. If provided, the function skips manual selection.
%   map - (Optional) A sensitivity map to be displayed alongside the movie for guidance in ROI selection.
%
%   Returns:
%   rois - A cell array of masks for each selected ROI.
%   traces - A matrix where each column represents the mean intensity trace of a corresponding ROI.
%
%   Description:
%   - The function displays the movie and sensitivity map (if provided) in separate subplots.
%   - Users can draw ROIs on the movie (if map were not provided) or sensitivity map. The ROIs are displayed as colored polygons.
%   - For each ROI, the function calculates the mean intensity trace over time.
%   - The process continues by user presses 'SPACE' until the user presses the 'Enter' key to end selection.
%
%   Example:
%   [rois, traces] = select_ROI(movie_data, ncols, nrows, time_vector, color_matrix, [], sensitivity_map);
%
%   Notes:
%   - The function requires an interactive MATLAB figure environment to work correctly.
%   - The 'colors' parameter should have as many rows as the maximum number of ROIs expected to be selected.
%   - If a 'mask' is provided, the function uses these masks instead of manual ROI selection.
%   - Press 'Space' to continue drawing ROIs; press 'Enter' to end the selection process.
%
% See also IMAGESC, IMSHOW, DRAWPOLYGON, POLY2MASK, MEAN.


% Initialize variables
rois = {};
traces = [];
movie_2D = reshape(movie, ncols, nrows, []);
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
    title(sprintf(['Sensitivity MAP\n\nPress SPACE to continue\nPress ENTER to end' ...
        '\n\n!!!!!!!!!!----------THE LAST ROI MUST BE BACKGROUND----------!!!!!!!!!!']));
    hold on;
    axis image;
end

% Display fluorescent image
subplot(2,2,2);
image_axe = gca;
im_adj = imadjust(uint16(mean(movie_2D, 3)));
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
    bwmask = mask;
    num_roi = max(bwmask(:));
    for i = 1:num_roi
        roi = (bwmask == i);
        trace = mean(movie(roi, :));
        traces = [traces trace'];
        boundary = bwboundaries(roi);
        plot(boundary{1}(:, 2), boundary{1}(:, 1), 'Color', colors(mod(i - 1, length(colors)) + 1, :), 'LineWidth', 2, 'Parent', image_axe);
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
