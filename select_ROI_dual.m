function [bwmask, traces] = select_ROI_dual(movie, movie_ca, ...
    nrows, ncols, nrows_ca, ncols_ca, t, t_ca, colors, mask, map, mask_ca, map_ca)
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
traces = [];
traces_ca = [];
bwmask = zeros(nrows,ncols);
movie_2D = reshape(movie, ncols, nrows, []);
movie_2D_ca = reshape(movie_ca, ncols_ca, nrows_ca, []);
map_merge = [map; map_ca];

% Set current figure to full screen and add keypress callback
fig = gcf;
set(fig,'Position',get(0,'Screensize'));
set(gcf, 'KeyPressFcn', @(src, event) set_space_pressed(event, fig));

% plot figure
if ~isempty(map) && ~isempty(map_ca)
    subplot(2,3,[1,4]);
    imagesc(map_merge);
    title('Sensitivity Map');
    map_axe = gca;
    hold on;

    % Plot Trace
    subplot(2,3,[2,3]);
    trace_axe = gca;
    hold on;
    xlabel('Time (s)');
    ylabel('Voltage Intensity');
    hold on;

    subplot(2,3,[5,6]);
    trace_axe_ca = gca;
    hold on;
    xlabel('Time (s)');
    ylabel('Calcium Intensity');
    hold on;
end

% % Mask provided, overlay ROI and traces
% if ~isempty(mask)
%     bwmask = mask;
%     num_roi = max(bwmask(:));
%     for i = 1:num_roi
%         roi = (bwmask == i);
%         trace = mean(movie(roi, :),1);
%         traces = [traces trace'];
%         boundary = bwboundaries(roi);
%         plot(boundary{1}(:, 2), boundary{1}(:, 1), 'Color', colors(mod(i - 1, length(colors)) + 1, :), 'LineWidth', 2, 'Parent', image_axe);
%         plot(t, trace, 'Color', colors(mod(i - 1, length(colors)) + 1, :), 'Parent', trace_axe);
%     end

% else
    % ROI selection process
    while true

        num_roi = max(bwmask(:)) + 1;
        roi_mask = drawpolygon('Color', colors(mod(num_roi, length(colors)), :), 'LineWidth', 1, 'Parent', map_axe);
        
        if mean(roi_mask.Position(:,2)) > ncols
            % select in calcium map
            mask = poly2mask(roi_mask.Position(:, 1), roi_mask.Position(:, 2)-ncols, ncols, nrows);
            boundary = bwboundaries(mask);
            % draw in voltage map
            plot(boundary{1}(:, 2), boundary{1}(:, 1), ...
                'Color', colors(mod(num_roi, length(colors)), :), 'LineWidth', 1, 'Parent', map_axe);hold on;
        else
            % select in voltage map
            mask = poly2mask(roi_mask.Position(:, 1), roi_mask.Position(:, 2), ncols, nrows);
            boundary = bwboundaries(mask);
            % draw in calcium map
            plot(boundary{1}(:, 2), boundary{1}(:, 1) + ncols, ...
                'Color', colors(mod(num_roi, length(colors)), :), 'LineWidth', 1, 'Parent', map_axe);hold on;
        end

        % save_mask
        bwmask(mask) = num_roi;
        
        % save trace
        trace = mean(movie(mask, :));
        traces = [traces trace'];
        trace_ca = mean(movie_ca(mask, :));
        traces_ca = [traces_ca trace_ca'];

        % plot trace
        plot(t, trace, 'Color', colors(mod(num_roi, length(colors)), :), 'Parent', trace_axe); hold on;
        plot(t_ca, trace_ca, 'Color', colors(mod(num_roi, length(colors)), :), 'Parent', trace_axe_ca); hold on;
        
        % Wait for user input
        fig.UserData.space = [];
        waitfor(fig, 'UserData');
        if strcmp(fig.UserData.space, 'stop')
            break;
        elseif strcmp(fig.UserData.space, 'spacePressed')
            continue;
        end
    end

% end
fprintf('Finished ROI selection\n');
end

function set_space_pressed(event, fig)
if strcmp(event.Key, 'space')
    fig.UserData.space = 'spacePressed';
elseif strcmp(event.Key, 'return')
    fig.UserData.space = 'stop';
end

end
