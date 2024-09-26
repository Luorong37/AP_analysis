function [bwmask, traces] = select_ROI(movie, nrows, ncols, mask, map)
% 
% ----------Write by Liu-Yang Luorong and ChatGPT----------
% ----------POWERED by Zoulab in Peking University----------
% Date: 24.09.10
% MATLAB Version: R2023a
% SELECT_ROI 交互式选择兴趣区域（ROIs）
%
%   该函数用于在给定的电影数据集中，允许用户交互式选择兴趣区域（ROIs）。函数显示电影数据
%   和可选的敏感度地图，用户可以通过绘制多边形定义感兴趣的区域。该函数返回每个ROI的二值掩
%   模和对应的平均强度时间轨迹。
%
%   语法：
%   [bwmask, traces] = select_ROI(movie, nrows, ncols, mask, map)
%
%   输入参数：
%   movie - 一个二维矩阵，表示电影数据，尺寸为 [ncols*nrows, nframes]。
%   nrows - 电影的行数。
%   ncols - 电影的列数。
%   mask - （可选）预定义的ROI掩模集，如果提供，函数将跳过手动选择。
%   map - （可选）敏感度地图，将与电影数据一起显示，帮助用户选择ROI。
%
%   输出参数：
%   bwmask - 每个选定ROI的二值掩模矩阵。
%   traces - 一个矩阵，每列表示对应ROI的平均强度轨迹。
%
%   功能描述：
%   - 函数将电影数据和可选的敏感度地图显示在不同的子图中。
%   - 用户可以在电影图像（或敏感度地图）上绘制多边形定义ROI，每个ROI用不同的颜色标记。
%   - 对于每个选定的ROI，函数会计算其随时间变化的平均强度轨迹。
%   - 用户按空格键继续选择ROI，按回车键结束选择过程。
%
%   示例：
%   [bwmask, traces] = select_ROI(movie_data, nrows, ncols, [], sensitivity_map);
%
%   注意事项：
%   - 该函数需要交互式MATLAB图形环境才能正常工作。
%   - 如果提供了' mask' 参数，则跳过手动ROI选择，直接使用预定义的掩模。
%   - 按空格键可以继续选择ROI，按回车键结束选择过程，按 'R' 可以重新选择当前ROI。
%
% See also IMAGESC, IMSHOW, DRAWPOLYGON, POLY2MASK, MEAN.

% Initialize variables
traces = [];
bwmask = zeros(nrows,ncols);
movie_2D = reshape(movie, ncols, nrows, []);
colors = lines(100);
key = ''; previous_key = '';
selected = false;

% Set current figure to full screen and add keypress callback
fig = gcf;
set(fig,'Position',get(0,'Screensize'));
set(gcf, 'KeyPressFcn', @(src, event) set_key_pressed(event, fig));

% Generate GUI for ROI selection
if ~isempty(map)
    [map_axe, image_axe, ~, trace_axe] = GUIwithmap(map, movie_2D);
else
    [image_axe, ~, trace_axe] = GUIwithoutmap(movie_2D);
    map_axe = [];
end

% Mask provided, overlay ROI and traces
if ~isempty(mask)
    [bwmask, traces] = select_by_mask(mask, colors, movie, traces, image_axe, trace_axe);
    selected = true;
end

% ROI selection process
while ~selected
    num_roi = max(bwmask(:)) + 1;
    color = colors(mod(num_roi - 1, length(colors)) + 1, :);
    % if (~isempty(map) && strcmp(key, {''})) || (any(strcmp(previous_key, {'','r'})) && strcmp(key, {'r'}))
    %     current_axe = map_axe;
    % else
    %     current_axe = image_axe;
    % end
    waitforbuttonpress;
    current_axe = gca;
    [mask,boundary] = axe_select(ncols, nrows, current_axe, color);
    bwmask(mask) = num_roi;
    trace = mean(movie(mask(:),:),1);
    traces = [traces trace'];
    
    [image_ROI, image_text, map_ROI, map_text] = ...
        plot_ROI(boundary, color, num_roi, image_axe, map_axe);
    cla(trace_axe); plot(trace', 'Color', color, 'Parent', trace_axe);
    
    previous_key = key;
    key = wait_for_key(fig);
    
    switch key
        case 'return'
            break;
        case  'space'
            key = '';
            continue;
        case  {'v', 'r'}
            % reset data
            bwmask(mask) = 0;
            if num_roi > 1
                traces = traces(:,end-1);
            else
                traces = [];
            end
            % reset plot
            delete(image_ROI);  delete(image_text);
            if ~isempty(map)
                delete(map_ROI); delete(map_text);
            end
            continue;
    end
end
fprintf('Finished ROI selection\n');
end

function set_key_pressed(event, fig)
if any(strcmp(event.Key, {'space', 'return', 'v', 'r'}))
    fig.UserData.space = event.Key;
end

end

function key = wait_for_key(fig)
% Waits for a keypress event to continue or stop ROI selection
    fig.UserData.space = [];
    waitfor(fig, 'UserData');
    key = fig.UserData.space;
end

function [map_axe, image_axe, im_adj, trace_axe] = GUIwithmap(map, movie_2D)
    subplot(2,2,2);
    map_axe = gca;
    imagesc(map);
    colorbar;
    title(sprintf('Sensitivity MAP\n\nPress SPACE to continue\nPress ENTER to end'));
    hold on;
    axis image;
    
    % Display fluorescent image
    subplot(1,2,1);
    image_axe = gca;
    % im_adj = imadjust(uint16(mean(movie_2D, 3)));
    % imshow(im_adj);
    im_adj = mean(movie_2D, 3);
    normalized_img = (im_adj - min(im_adj(:))) / (max(im_adj(:)) - min(im_adj(:)));
    imshow(normalized_img);
    hold on;
    title('Fluorescent Image');
    
    % Plot Trace
    subplot(2,2,4);
    trace_axe = gca;
    hold on;
    xlabel('Time (s)');
    ylabel('Intensity');
end

function [image_axe, im_adj, trace_axe] = GUIwithoutmap(movie_2D)
    % Display fluorescent image without map
    subplot(1,2,1);
    image_axe = gca;
    % im_adj = imadjust(uint16(mean(movie_2D, 3)));
    im_adj = uint16(mean(movie_2D, 3));
    imshow(im_adj,[min(im_adj,[],'all'),max(im_adj,[],'all')]);
    hold on;
    title('Fluorescent Image');
    
    % Plot Trace
    subplot(1,2,2);
    trace_axe = gca;
    hold on;
    xlabel('Time (s)');
    ylabel('Intensity');
end

function [bwmask, traces] = select_by_mask(bwmask, colors, movie, traces, image_axe, trace_axe)
num_roi = max(bwmask(:));
for i = 1:num_roi
    color = colors(mod(i - 1, length(colors)) + 1, :);
    roi = (bwmask == i);
    trace = mean(movie(roi, :),1);
    traces = [traces trace'];
    boundary = cell2mat(bwboundaries(roi));
    plot(boundary(:, 2), boundary(:, 1), 'Color', color, 'LineWidth', 2, 'Parent', image_axe);
    plot(trace', 'Color', color, 'Parent', trace_axe);
    % 标注ROI编号
    text(mean(boundary(:, 2)) + 12, mean(boundary(:, 1)) - 12, num2str(i), ...
        'Color', color, 'FontSize', 12, 'Parent', image_axe); hold on;
    
end
end

function [mask,boundary] = axe_select(ncols, nrows, axe, color)
    roi_mask = drawpolygon('Color', color, 'LineWidth', 1, 'Parent', axe);
    mask = poly2mask(roi_mask.Position(:, 1), roi_mask.Position(:, 2), ncols, nrows);
    boundary = cell2mat(bwboundaries(mask));
    delete(roi_mask)
end

function [image_ROI, image_text, map_ROI, map_text] = ...
plot_ROI(boundary, color, num_roi, image_axe, map_axe)
if ~isempty(map_axe)
    map_ROI = plot(boundary(:, 2), boundary(:, 1), 'Color', color, 'LineWidth', 1, 'Parent', map_axe);
    map_text = text(mean(boundary(:, 2)) + 12, mean(boundary(:, 1)) - 12, num2str(num_roi), 'Color', color, 'FontSize', 12, 'Parent', map_axe); hold on;
else
    map_ROI = [];
    map_text = [];
end
    image_ROI = plot(boundary(:, 2), boundary(:, 1), 'Color', color, 'LineWidth', 1, 'Parent', image_axe);
    image_text = text(mean(boundary(:, 2)) + 12, mean(boundary(:, 1)) - 12, num2str(num_roi), 'Color', color, 'FontSize', 12, 'Parent', image_axe); hold on;
end