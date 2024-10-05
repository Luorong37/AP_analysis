function [rois, traces] = select_ROI(movie, nrows, ncols, mask, map)
% 
% ----------Written by Liu-Yang Luorong and ChatGPT----------
% ----------POWERED by Zoulab in Peking University----------
% Date: 24.10.01
% MATLAB Version: R2023a
% select_ROI 交互式选择兴趣区域（ROIs）
%
%   此函数允许用户在给定的电影数据集中，交互式选择感兴趣区域（ROIs）。
%   用户可以在显示的图像或敏感度图（如有）上手动绘制多边形，定义感兴趣区域。
%   如果提供了预定义的掩码，则跳过手动选择过程。
%   最终输出每个ROI的二值掩码和对应的平均强度时间轨迹。
%
%   语法：
%   [rois, traces] = select_ROI(movie, nrows, ncols, mask, map)
%
%   输入参数：
%   movie - 一个二维矩阵，表示电影数据，尺寸为 [ncols*nrows, nframes]。
%   nrows - 电影的行数。
%   ncols - 电影的列数。
%   mask - （可选）预定义的ROIs掩模集，若提供则跳过手动选择过程。
%   map  - （可选）敏感度图（用于帮助选择ROIs的参考图）。
%
%   输出参数：
%   rois   - 结构体，包含以下字段：
%            - bwmask: 每个选定ROI的二值掩模矩阵。
%            - boundary: 每个选定ROI的边界坐标。
%            - position: 每个选定ROI的顶点坐标。
%   traces - 一个矩阵，每列表示对应ROI的平均强度时间轨迹。
%
%   功能描述：
%   - 显示电影数据和可选的敏感度地图，用户可以通过绘制多边形定义ROIs。
%   - 对于每个选定的ROI，计算其随时间变化的平均强度轨迹。
%   - 用户按空格键继续选择ROIs，按回车键结束选择过程，按 'R' 重新选择当前ROI。
%
%   示例：
%   [rois, traces] = select_ROI(movie_data, nrows, ncols, [], sensitivity_map);
%
%   注意事项：
%   - 函数需要交互式MATLAB图形环境才能正常工作。
%   - 如果提供了'mask'参数，则跳过手动ROIs选择，直接使用预定义的掩模。
%   - 用户可以按空格键继续选择ROIs，按回车键结束选择，按 'R' 重新选择当前ROI。

% 初始化变量
traces = [];
rois.bwmask = zeros(nrows, ncols);  % 初始化二值掩码
rois.boundary = {};  % 存储每个ROI的边界
rois.position = {};  % 存储每个ROI的顶点坐标
movie_2D = reshape(movie, ncols, nrows, []);  % 重塑电影数据为二维图像
colors = lines(100);  % 使用预定义的颜色线来标记不同的ROI
% key = '';  % 用于跟踪用户按键
selected = false;  % 标记是否已完成选择

% 设置当前图形窗口为全屏，并添加按键回调函数
fig = gcf;
set(gcf, 'KeyPressFcn', @(src, event) set_key_pressed(event, fig));  % 设置按键事件
set(fig, 'Position', get(0, 'Screensize'));

% 如果提供了敏感度图，创建带有地图的GUI界面；否则，创建仅显示电影的界面
if ~isempty(map)
    [map_axe, image_axe, ~, trace_axe] = GUIwithmap(map, movie_2D);
else
    [image_axe, ~, trace_axe] = GUIwithoutmap(movie_2D);
    map_axe = [];
end
title(sprintf(['Click an axe to select ROIs.' ...
    '\nAfter an ROI selection, Press SPACE to continue, Press ENTER to end, Press R to reselect' ...
    '\nBefore an ROI selction, Press Q to quit, Press R to reselect.']));

% 如果提供了掩码，跳过选择，直接使用提供的掩码并计算轨迹
if ~isempty(mask)
    traces = select_by_mask(mask, colors, movie, image_axe, trace_axe);
    selected = true;
    rois.bwmask = mask;
end
axe_labels = {};

% 交互式选择ROI的过程
while ~selected
    num_rois = max(rois.bwmask(:)) + 1;  % 当前ROI编号
    color = colors(mod(num_rois - 1, length(colors)) + 1, :);  % 选择颜色

    % 等待用户绘制ROI
    title(trace_axe,sprintf(['Click an axe to select ROIs.' ...
    '\nPress Q/Enter to quit, Press R to reselect.']));
    set(gcf, 'CurrentCharacter', char(0));
    
    waitforbuttonpress;
    current_axe = gca;
    
    % 捕捉用户的按键
    key_pressed = get(gcf, 'CurrentCharacter');
    switch key_pressed
        case {char(13),'q'} % 按下 'q'和 Enter 键退出循环
        disp('ROI selection cancelled by user.');
        break; 

        case 'r'
        key = 'r';

        otherwise
        % 用户在当前轴选择ROI
        [mask, boundary, position] = axe_select(ncols, nrows, current_axe, color);
        rois.bwmask(mask) = num_rois;  % 将掩码添加到总掩码
        rois.boundary = [rois.boundary, boundary];  % 记录边界
        rois.position = [rois.position, position];  % 记录位置
    
        % 计算ROI区域的平均强度轨迹
        trace = mean(movie(mask(:), :), 1);
        traces = [traces, trace'];
        
        % 绘制ROI和轨迹
        [image_roi, image_text, map_roi, map_text] = plot_rois(boundary, color, num_rois, image_axe, map_axe);
        cla(trace_axe); plot(trace', 'Color', color, 'Parent', trace_axe);  % 绘制信号轨迹
        axe_labels = [axe_labels;{image_roi, image_text, map_roi, map_text}];

        % 等待用户的下一步操作
        title(trace_axe,sprintf(['An ROI selected.' ...
    '\nPress SPACE to continue, Press Q/Enter to end, Press R to reselect']));
        key = wait_for_key(fig);
    end

    switch key
        case {'return','q'}
            break;  % 按回车键结束选择
        case 'space'
            % key = '';  % 按空格键继续选择下一个ROI
            continue;
        case {'v', 'r'}
            deleted_nroi = max(rois.bwmask(:));
            % 清除当前ROI的绘图
            delete(axe_labels{deleted_nroi,1});  delete(axe_labels{deleted_nroi,2});
            if ~isempty(map)
                delete(axe_labels{deleted_nroi,3}); delete(axe_labels{deleted_nroi,4});
            end
            axe_labels = axe_labels(1:end-1,:);
            cla(trace_axe);

            % 重新选择当前ROI
            rois.bwmask(rois.bwmask == max(rois.bwmask(:))) = 0;
            if max(rois.bwmask(:)) > 0
                traces = traces(:, 1:end-1);
                rois.boundary = rois.boundary(1:end-1);
                rois.position = rois.position(1:end-1);
                plot(traces(:, end), 'Color',colors(mod(max(rois.bwmask(:))-1, length(colors)) + 1, :),'Parent', trace_axe);
            else
                traces = [];
                rois.boundary = {};
                rois.position = {};
            end

            continue;
    end
end
    disp('Finished ROIs selection');
end

function set_key_pressed(event, fig)
% 处理按键事件，记录按键值
if any(strcmp(event.Key, {'space', 'return', 'v', 'r','q'}))
    fig.UserData.space = event.Key;
end
end

function key = wait_for_key(fig)
% 等待用户按键事件，返回按键值
    fig.UserData.space = [];
    waitfor(fig, 'UserData');
    key = fig.UserData.space;
end

function [map_axe, image_axe, im_adj, trace_axe] = GUIwithmap(map, movie_2D)
% 创建带有敏感度地图和电影的GUI界面
    subplot(2,2,2);
    map_axe = gca;
    imagesc(map);  % 显示敏感度图
    colorbar;
    title(sprintf('Sensitivity MAP'));
    hold on;
    axis image;

    % 显示电影数据
    subplot(1,2,1);
    image_axe = gca;
    im_adj = mean(movie_2D, 3);  % 计算电影的平均图像
    normalized_img = (im_adj - min(im_adj(:))) / (max(im_adj(:)) - min(im_adj(:)));  % 归一化图像
    imshow(normalized_img);
    hold on;
    title('Fluorescent Image');

    % 显示轨迹的轴
    subplot(2,2,4);
    trace_axe = gca;
    hold on;
    xlabel('Time (s)');
    ylabel('Intensity');
end

function [image_axe, im_adj, trace_axe] = GUIwithoutmap(movie_2D)
% 创建仅显示电影的GUI界面
    subplot(1,2,1);
    image_axe = gca;
    im_adj = uint16(mean(movie_2D, 3));  % 显示电影的平均图像
    imshow(im_adj, [min(im_adj,[],'all'), max(im_adj,[],'all')]);
    hold on;
    title('Fluorescent Image');

    % 显示轨迹的轴
    subplot(1,2,2);
    trace_axe = gca;
    hold on;
    xlabel('Time (s)');
    ylabel('Intensity');
end

function traces = select_by_mask(bwmask, colors, movie, image_axe, trace_axe)
% 使用给定掩码选择ROI并计算轨迹
num_rois = max(bwmask(:));  % ROI数量
traces = [];
for i = 1:num_rois
    color = colors(mod(i - 1, length(colors)) + 1, :);  % 选择颜色
    mask = (bwmask == i);  % 获取当前ROI掩码
    trace = mean(movie(mask, :), 1);  % 计算平均强度
    boundary = cell2mat(bwboundaries(mask));  % 提取边界

    traces = [traces, trace'];  % 保存轨迹
    rois.boundary{i} = boundary;

    % 在图像上绘制ROI边界和轨迹
    plot(boundary(:, 2), boundary(:, 1), 'Color', color, 'LineWidth', 2, 'Parent', image_axe);
    plot(trace', 'Color', color, 'Parent', trace_axe);

    % 标注ROI编号
    text(mean(boundary(:, 2)) + 12, mean(boundary(:, 1)) - 12, num2str(i), ...
        'Color', color, 'FontSize', 12, 'Parent', image_axe); hold on;
end
end

function [mask, boundary, position] = axe_select(ncols, nrows, axe, color)
% 在给定的轴上让用户选择ROI并生成掩码
    rois_polygon = drawpolygon('Color', color, 'LineWidth', 1, 'Parent', axe);  % 绘制多边形
    position = rois_polygon.Position;  % 获取顶点坐标
    mask = poly2mask(position(:, 1), position(:, 2), ncols, nrows);  % 生成二值掩码
    boundary = cell2mat(bwboundaries(mask));  % 提取边界
    delete(rois_polygon);  % 删除绘制的多边形
end

function [image_rois, image_text, map_rois, map_text] = plot_rois(boundary, color, num_rois, image_axe, map_axe)
% 绘制ROI的边界和编号
if ~isempty(map_axe)
    map_rois = plot(boundary(:, 2), boundary(:, 1), 'Color', color, 'LineWidth', 1, 'Parent', map_axe);  % 在敏感度图上绘制
    map_text = text(mean(boundary(:, 2)) + 12, mean(boundary(:, 1)) - 12, num2str(num_rois), ...
                    'Color', color, 'FontSize', 12, 'Parent', map_axe); hold on;
else
    map_rois = [];
    map_text = [];
end
    image_rois = plot(boundary(:, 2), boundary(:, 1), 'Color', color, 'LineWidth', 1, 'Parent', image_axe);  % 在图像上绘制
    image_text = text(mean(boundary(:, 2)) + 12, mean(boundary(:, 1)) - 12, num2str(num_rois), ...
                      'Color', color, 'FontSize', 12, 'Parent', image_axe); hold on;
end
