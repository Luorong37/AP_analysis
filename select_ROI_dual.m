function [bwmask, bwmask_ca, traces, traces_ca, offset] = select_ROI_dual(movie, movie_ca, ...
    nrows, ncols, correct, map, map_ca, mask, mask_ca, offset)

colors = lines(100);
if nargin < 7, correct = false; 
end
% movie_3d = reshape(movie,ncols,nrows,[]);
% movie_3d_ca = reshape(movie_ca,ncols,nrows,[]);

traces = []; traces_ca = [];
bwmask = zeros(nrows,ncols); bwmask_ca = zeros(nrows,ncols);
x_offset = 0; y_offset = 0;
key = ''; 
selected = false;

fig = gcf; set(fig,'Position',get(0,'Screensize'));
set(gcf, 'KeyPressFcn', @(src, event) set_key_pressed(event, fig));
[map_axe, map_ca_axe, img_axe, img_ca_axe, trace_axe, trace_ca_axe] = GUI(map, map_ca, movie, movie_ca, ncols, nrows);

if ~isempty(mask) && ~isempty(mask_ca)
    [traces, traces_ca, bwmask, bwmask_ca] = select_by_mask(mask, colors, mask_ca, map_axe, map_ca_axe, img_axe, img_ca_axe, movie, traces, movie_ca, traces_ca, trace_axe, trace_ca_axe);
    selected = true;
end

if correct && ~selected
    [x_offset, y_offset] = handle_manual_correction(map_axe,map_ca_axe, ncols, nrows, img_axe, img_ca_axe);
    offset = [x_offset, y_offset];
    corrected = true;
elseif ~isempty(offset)
    corrected = true;
end

while ~selected
    num_roi = max(bwmask(:)) + 1;
    color = colors(num_roi,:);
    waitforbuttonpress;
    current_axe = gca;
    if isequal(current_axe,img_axe)||isequal(current_axe, img_ca_axe)
        % plot in img
        [mask_v, mask_ca, boundary_v, boundary_ca] = handle_offset_select(x_offset, y_offset, ncols, nrows, img_axe, img_ca_axe,current_axe, color);
    else
        % plot in map
        [mask_v, mask_ca, boundary_v, boundary_ca] = handle_offset_select(x_offset, y_offset, ncols, nrows, map_axe, map_ca_axe,current_axe, color);
    end

    % Plot corrected ROIs on images
    [map_roi_v, map_roi_ca, img_roi_v, img_roi_ca, text_v, text_ca, text_img_v, text_img_ca] = ......
    plotROI(boundary_v, boundary_ca, map_axe, map_ca_axe, img_axe, img_ca_axe, num_roi, color);

    % save trace and mask
    bwmask(mask_v) = num_roi; bwmask_ca(mask_ca) = num_roi;

    trace = mean(movie(mask_v(:),:),1);traces = [traces trace'];
    trace_ca = mean(movie_ca(mask_ca(:),:), 1); traces_ca = [traces_ca trace_ca'];
    cla(trace_axe); plot(trace', 'Color', color, 'Parent', trace_axe);  hold on;
    cla(trace_ca_axe); plot(trace_ca', 'Color', color, 'Parent', trace_ca_axe);  hold on;

    % 
    key = wait_for_key(fig);
    if strcmp(key, 'return')
        break;
    elseif strcmp(key, 'space')
        key = '';
        continue;
    elseif any(strcmp(key, {'v', 'c', 'r'}))
        % reset data
        bwmask(mask_v) = 0; bwmask_ca(mask_ca) = 0;
        if max(bwmask(:)) == 0
            traces = []; traces_ca = [];
        else
            traces = traces(:,end-1); traces_ca = traces_ca(:,end-1);
        end

        % reset plot
        delete(map_roi_v); delete(map_roi_ca); delete(img_roi_v); delete(img_roi_ca); 
        delete(text_v); delete(text_ca); delete(text_img_v); delete(text_img_ca);  
        
        continue;
    end
end
fprintf('Finished ROI selection\n');
end

function set_key_pressed(event, fig)
if any(strcmp(event.Key, {'space', 'return', 'v', 'c', 'r'}))
    fig.UserData.space = event.Key;
end

end
function key = wait_for_key(fig)
% Waits for a keypress event to continue or stop ROI selection

    fig.UserData.space = [];
    waitfor(fig, 'UserData');
    key = fig.UserData.space;
end

function     [map_axe, map_ca_axe, img_axe, img_ca_axe, trace_axe, trace_ca_axe] = ...
    GUI(map, map_ca, movie, movie_ca, ncols, nrows)

subplot(2,3,1); imagesc(map); title(sprintf('Sensitivity Map\nclick to draw ROI on this map\nSelect this first to correct offset if need'));
map_axe = gca; hold on;
subplot(2,3,4); imagesc(map_ca); title(sprintf('Sensitivity Map Ca\nclick to draw ROI on this map\nSelect this second to correct offset if need'));
map_ca_axe = gca; hold on;

img_axe = subplot(2,3,2); hold on;
img = reshape(mean(movie,2), ncols, nrows);
normalized_img = (img - min(img(:))) / (max(img(:)) - min(img(:)));
imshow(normalized_img);
title('Voltage image'); hold on;

img_ca_axe = subplot(2,3,5); hold on;
img_ca = reshape(mean(movie_ca,2), ncols, nrows);
normalized_img_ca = (img_ca - min(img_ca(:))) / (max(img_ca(:)) - min(img_ca(:)));
imshow(normalized_img_ca);
title('Calcium image'); hold on;

subplot(2,3,3); xlabel('Time '); ylabel('Voltage Intensity');
trace_axe = gca; hold on;
subplot(2,3,6);  xlabel('Time'); ylabel('Calcium Intensity');
trace_ca_axe = gca; hold on;
% set(map_axe, 'ButtonDownFcn', @(~,~) set_gca(gcf, map_axe));
% set(map_ca_axe, 'ButtonDownFcn', @(~,~) set_gca(gcf, map_ca_axe));
% set(img_axe, 'ButtonDownFcn', @(~,~) set_gca(gcf, img_axe));
% set(img_ca_axe, 'ButtonDownFcn', @(~,~) set_gca(gcf, img_ca_axe));

end

function [x_offset, y_offset] = handle_manual_correction(map_axe,map_ca_axe, ncols, nrows, img_axe, img_ca_axe)
% Handles manual correction for ROIs in both maps and calculates the offset
    waitforbuttonpress;
    if map_axe == gca || img_axe == gca
        first_axe = 'v'; second_axe = 'c';
    elseif map_ca_axe == gca || img_ca_axe == gca
        first_axe = 'c'; second_axe = 'v';
    end

    % Draw ROI in voltage map
    first_roi = drawpolygon('LineWidth', 2, 'Parent', gca);
    mask = poly2mask(first_roi.Position(:, 1), first_roi.Position(:, 2), ncols, nrows);
    boundary = cell2mat(bwboundaries(mask));
    
    % Draw corrected ROI in calcium map
    if map_axe == gca
        set(gcf, 'CurrentAxes', map_ca_axe);
    elseif map_ca_axe == gca
        set(gcf, 'CurrentAxes', map_axe);
    elseif img_axe == gca
        set(gcf, 'CurrentAxes', img_ca_axe);
    elseif img_ca_axe == gca
        set(gcf, 'CurrentAxes', img_axe);
    end
    ref_roi = plot(boundary(:, 2), boundary(:, 1), 'Color', 'r', 'LineWidth', 1, 'Parent', gca); hold on;
    
    % Select calcium ROI and calculate offset
    second_roi = drawpolygon('Color', 'g', 'LineWidth', 2, 'Parent', gca); hold on;
    % mask_ca = poly2mask(second_roi.Position(:, 1), second_roi.Position(:, 2), ncols, nrows);
    % boundary_ca = cell2mat(bwboundaries(mask_ca));

    % calculate offset
    if first_axe == 'v'
    [cent_v.x, cent_v.y] = centroid(polyshape(first_roi.Position));
    [cent_ca.x, cent_ca.y] = centroid(polyshape(second_roi.Position));
    elseif first_axe == 'c'
    [cent_ca.x, cent_ca.y] = centroid(polyshape(first_roi.Position));
    [cent_v.x, cent_v.y] = centroid(polyshape(second_roi.Position));
    end
    offset = [cent_v.x, cent_v.y] - [cent_ca.x, cent_ca.y];
    x_offset = offset(1);
    y_offset = offset(2);
    plot([cent_v.x, cent_ca.x], [cent_v.y, cent_ca.y], 'black', 'LineWidth', 3, 'Marker', 'o'); hold on;
    fprintf('offset = %f, %f\n', offset);
    
    delete(first_roi);
    delete(second_roi);
    delete(ref_roi);
end

function [mask_v, mask_ca, boundary_v, boundary_ca] = ...
handle_offset_select(x_offset, y_offset, ncols, nrows, v_axe, ca_axe, axe1 , color)

% % Handles ROI selection 
%     peak = true;
% 
%     while peak
%     waitforbuttonpress;
%     axe1 = gca;
%     if isequal(axe1, v_axe) || isequal(axe1, ca_axe)
%         peak = false;
%     end
    % end
    roi1 = drawpolygon('Color', color, 'LineWidth', 2, 'Parent', axe1);
    mask1 = poly2mask(roi1.Position(:, 1), roi1.Position(:, 2), ncols, nrows);
    boundary1 = cell2mat(bwboundaries(mask1));
    
    % Calculate corresponding ROI
    if axe1 == v_axe
        roi2 = translate(polyshape(roi1.Position),[-x_offset, -y_offset]).Vertices;
    elseif axe1 == ca_axe
        roi2 = translate(polyshape(roi1.Position),[x_offset, y_offset]).Vertices;
    end
    mask2 = poly2mask(roi2(:, 1), roi2(:, 2), ncols, nrows);
    boundary2 = cell2mat(bwboundaries(mask2));
    
    % allocate roi to voltage and calcium
    if axe1 == v_axe
        mask_v = mask1; mask_ca = mask2;
        boundary_v = boundary1; boundary_ca = boundary2;
    elseif axe1 == ca_axe
        mask_v = mask2; mask_ca = mask1;
        boundary_v = boundary2; boundary_ca = boundary1;
    end
    delete(roi1);
end

function [map_roi_v, map_roi_ca, img_roi_v, img_roi_ca, text_v, text_ca, text_img_v, text_img_ca] = ......
    plotROI(boundary_v, boundary_ca, map_axe, map_ca_axe, img_axe, img_ca_axe, num_roi, color)
    % Plot corrected ROIs on images
    map_roi_v = plot(boundary_v(:, 2), boundary_v(:, 1), 'Color', color, 'LineWidth', 1, 'Parent', map_axe); hold on;
    map_roi_ca = plot(boundary_ca(:, 2), boundary_ca(:, 1), 'Color', color, 'LineWidth', 1, 'Parent', map_ca_axe); hold on;
    img_roi_v = plot(boundary_v(:, 2), boundary_v(:, 1), 'Color', color, 'LineWidth', 1, 'Parent', img_axe); hold on;
    img_roi_ca = plot(boundary_ca(:, 2), boundary_ca(:, 1), 'Color', color, 'LineWidth', 1, 'Parent', img_ca_axe); hold on;
    
    % 标注ROI编号
    text_v = text(mean(boundary_v(:, 2)), mean(boundary_v(:, 1)), num2str(num_roi), 'Color', 'k', 'FontSize', 12, 'Parent', map_axe); hold on;
    text_ca = text(mean(boundary_ca(:, 2)), mean(boundary_ca(:, 1)), num2str(num_roi), 'Color', 'k', 'FontSize', 12, 'Parent', map_ca_axe); hold on;
    text_img_v = text(mean(boundary_v(:, 2)), mean(boundary_v(:, 1)), num2str(num_roi), 'Color', 'k', 'FontSize', 12, 'Parent', img_axe); hold on;
    text_img_ca = text(mean(boundary_ca(:, 2)), mean(boundary_ca(:, 1)), num2str(num_roi), 'Color', 'k', 'FontSize', 12, 'Parent', img_ca_axe); hold on;
end

function [traces, traces_ca, bwmask, bwmask_ca] = select_by_mask(mask, colors, mask_ca, map_axe, map_ca_axe, img_axe, img_ca_axe, movie, traces, movie_ca, traces_ca, trace_axe, trace_ca_axe)
for i = 1: max(mask(:))
    fprintf('plot ROI %d\n',i);
    color = colors(i,:);
    mask_vi = (mask == i);
    mask_cai = (mask_ca == i);
    boundary_v = cell2mat(bwboundaries(mask_vi));
    boundary_ca = cell2mat(bwboundaries(mask_cai));
    plotROI(boundary_v, boundary_ca, map_axe, map_ca_axe, img_axe, img_ca_axe, i, color);
    trace = mean(movie(mask_vi(:),:),1);
    [trace_corrected, ~] = fit_exp2(trace'); traces = [traces trace'];
    trace_ca = mean(movie_ca(mask_cai(:),:), 1);
    [trace_corrected_ca, ~] = fit_exp2(trace_ca'); traces_ca = [traces_ca trace_ca'];
    cla(trace_axe); plot(trace_corrected, 'Color', color, 'Parent', trace_axe);  hold on;
    cla(trace_ca_axe); plot(trace_corrected_ca, 'Color', color, 'Parent', trace_ca_axe);  hold on;
    % pause(0.5);
end
bwmask = mask;
bwmask_ca = mask_ca;
end