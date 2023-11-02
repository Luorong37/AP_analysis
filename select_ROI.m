%% Select ROI


function  [rois, traces] = select_ROI(intensity_time_series, num_rows,num_cols,t, colors, mask, map)
% -----------------------------with Mask-----------------------------

if isempty(mask) && isempty(map)
   rois = {}; % Initialize cell array to store the ROIs
    traces = [];
    fig = gcf;
    tolerance = 0.2; % Set a tolerance for double-click detection (in seconds)

    % Visualize Fluorecent Image
    subplot(1,2,1);
    image_axe = gca;
    % draw average z stack
    im_adj = imadjust(uint16(mean(reshape(intensity_time_series, num_cols,num_rows,[]),3)));
    Image = imshow(im_adj);
    hold on;
    title(sprintf('Fluorecent Image\nslow double click to continue, quick double click to end\nthe last trace should be background'));

    % Visualize Trace Image
    subplot(1,2,2);
    trace_axe = gca;
    hold on;
    xlabel('Time (s)');
    ylabel('Intensity');

    %Select ROI
    while true
        %Select ROI in Sensitivity Map aftermatlab2017b
        roi_mask = drawpolygon('Color', colors(mod(length(rois), length(colors))+1,:), 'LineWidth', 1, 'Parent', image_axe); % Select a color based on the number of existing ROIs
        hold on;

        %     %Select ROI in Sensitivity Map before matlab2017b
        %     h = impoly(map_axe);
        %     setColor(h, colors{mod(length(rois), length(colors))+1});
        %     setPositionConstraintFcn(h, @getPositionConstraint);
        %     mask = createMask(h);

        %Store Mask and Trace
        mask = poly2mask(roi_mask.Position(:,1), roi_mask.Position(:,2), size(Image.CData, 1), size(Image.CData, 2));
        rois{end+1} = mask;
        trace = mean(intensity_time_series(mask,:));% Compute mean intensity-time trace for the selected ROI
        traces = [traces trace'];

        % Plot ROI in Fluorescent Image
        boundary = bwboundaries(mask); % get boundary
        plot(boundary{1}(:,2), boundary{1}(:,1), 'Color', colors(mod(length(rois)-1, length(colors))+1,:), 'LineWidth', 1,'Parent',image_axe);
        hold on;

        % Plot the intensity-time trace
        plot(t, trace, 'Color', colors(mod(length(rois)-1, length(colors))+1,:),'Parent',trace_axe);
        hold on;


        % Allow the user to select an ROI or finish
        % Wait for a button press or a click
        t1 = tic; % Start a timer
        waitforbuttonpress;
        t2 = toc(t1); % Get the elapsed time
        waitforbuttonpress;
        % Check if the button press was a double-click
        t3 = toc(t1); % Get the elapsed time again
        % If it's a double-click, exit the loop
        if t3 -t2 < tolerance
            break;
        end
    end
    %


    fprintf('Finished ROI selection\n')

elseif ~isempty(mask)
    rois = {}; % Initialize cell array to store the ROIs
    traces = [];
    fig = gcf;
    set(fig,'Position',get(0,'Screensize'));
    tolerance = 0.2; % Set a tolerance for double-click detection (in seconds)



    % % Visualize Sensitivity Map
    % subplot('Position',position{1});
    % map_axe = gca; % Get the handle of the current axes
    % Map =imagesc(mask_input);
    % hold on;
    % colorbar;
    % title(sprintf('Sensitivity Map \n while double click after each ROI selection: Quick-end; Slow-continue'));
    % axis image;

    % Visualize Fluorecent Image
    subplot(2,1,1);
    image_axe = gca;
    im_adj = imadjust(uint16(mean(reshape(intensity_time_series,num_cols, num_rows,[]),3)));
    Image = imshow(im_adj);
    hold on;
    title('Fluorecent Image');


    % Visualize Trace
    subplot(2,1,2);
    trace_axe = gca;
    hold on;
    xlabel('Time (s)');
    ylabel('Intensity');
    roi_amount = size(mask_input,2);
    rois = mask_input;
    for i = 1: roi_amount
        its = reshape(intensity_time_series,num_cols, num_rows,[]);
        trace = mean(its(mask_input{i},:));
        traces = [traces trace'];
        boundary = bwboundaries(mask_input{i}); % get boundary
        plot(boundary{1}(:,2), boundary{1}(:,1), 'Color', colors(mod(i-1, length(colors))+1,:), 'LineWidth', 1,'Parent',image_axe);
        hold on;
        %     plot(boundary{1}(:,2), boundary{1}(:,1), 'Color', colors(mod(i-1, length(colors))+1,:), 'LineWidth', 1,'Parent',map_axe);
        %     hold on;
        plot(t, trace, 'Color', colors(mod(i-1, length(colors))+1,:),'Parent',trace_axe);
        hold on;
    end


    fprintf('Finished ROI selection\n')
elseif ~isempty(Map)
    % default setting
    figure()
    position = {[0.03 0.05 0.45 0.9] [0.55 0.55 0.4 0.4] [0.55 0.1 0.4 0.4]};
    rois = {}; % Initialize cell array to store the ROIs
    traces = [];
    fig = gcf;
    set(fig,'Position',get(0,'Screensize'));
    % tolerance = 0.2; % Set a tolerance for double-click detection (in seconds)
    % spacePressed = false;
    %     设置一个回调函数来更新标志变量
    set(gcf, 'KeyPressFcn', @(src, event) set_space_pressed(event, fig));



    % Visualize Sensitivity Map
    subplot('Position',position{1});
    map_axe = gca; % Get the handle of the current axes
    Map =imagesc(map);
    hold on;
    colorbar;
    title(sprintf('Sensitivity MAP \n PRESS for continue; ENTER for end'));
    axis image;

    % Visualize Fluorecent Image
    subplot('Position',position{2});
    image_axe = gca;
    %im_adj = imadjust(uint16(mean(intensity_time_series,3)));
    im_adj = imadjust(uint16(mean(reshape(intensity_time_series,num_cols,num_rows,[]),3)));
    Image = imshow(im_adj);
    hold on;
    title('Fluorecent Image');


    % Visualize Trace
    subplot('Position',position{3});
    trace_axe = gca;
    hold on;
    xlabel('Time (s)');
    ylabel('Intensity');

    %Select ROI
    while true
        %Select ROI in Sensitivity Map (after matlab2017b)
        roi_mask = drawpolygon('Color', colors(mod(length(rois), length(colors))+1,:), 'LineWidth', 1, 'Parent', map_axe); % Select a color based on the number of existing ROIs
        hold on;

        %Store Mask and Trace
        mask = poly2mask(roi_mask.Position(:,1), roi_mask.Position(:,2), size(Map.CData, 1), size(Map.CData, 2));
        rois{end+1} = mask;
        trace = mean(intensity_time_series(rois{end},:));% Compute mean intensity-time trace for the selected ROI
        traces = [traces trace'];

        % Plot ROI in Fluorescent Image
        boundary = bwboundaries(mask); % get boundary
        plot(boundary{1}(:,2), boundary{1}(:,1), 'Color', colors(mod(length(rois)-1, length(colors))+1,:), 'LineWidth', 1,'Parent',image_axe);
        hold on;

        % Plot the intensity-time trace
        plot(t, trace, 'Color', colors(mod(length(rois)-1, length(colors))+1,:),'Parent',trace_axe);
        hold on;

        % 等待用户按下空格键或执行其他操作
        fig.UserData = [];
        while true
            waitfor(fig, 'UserData');
            if strcmp(fig.UserData, 'stop')
                break;
            elseif strcmp(fig.UserData, 'spacePressed')
                break;
            end
        end
        %
        if strcmp(fig.UserData, 'stop')
            break;
        end

    end

    fprintf('Finished ROI selection\n')
end
end
