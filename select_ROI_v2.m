%% Find Cell Mask by sensitivity map
%v2.0
% 2.0 update: can input ROI or input without ROI

function  [rois, traces] = select_ROI_v2(intensity_time_series, t, colors, num_rows, num_cols)

% -----------------------------without Mask and Map-----------------------------
% default setting
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

end
