function background_correction(movie, ncols, nrows, rois)
    % 背景校正函数，显示图像、ROI和背景的校正界面
    % movie: 视频数据（时间序列图像）
    % ncols: 图像的列数
    % nrows: 图像的行数
    % rois: 包含每个ROI的位置、边界和掩码数据的结构体
    
    % 将movie重塑为3D图像并计算每帧的均值作为静态图像
    movie3D = reshape(movie, ncols, nrows, []);
    movie2D = mean(movie3D, 3);
    img = (movie2D - min(movie2D(:))) / (max(movie2D(:)) - min(movie2D(:))); % 标准化图像

    % 初始化背景单元格数组，用于存储每个ROI的背景信息
    background = cell(max(rois.bwmask(:)), 1);

    % 创建主窗口界面
    screenSize = get(0, 'ScreenSize');
    fig = uifigure('Position', [0 0 screenSize(3) screenSize(4)], 'Name', 'Background Correction');

    % 创建用于显示原始图像的轴
    axImage = uiaxes('Parent', fig, 'Position', [screenSize(3)/20, screenSize(3)/20, screenSize(4)*0.8, screenSize(4)*0.8]);
    title(axImage, 'Origin Image','HorizontalAlignment', 'center');
    imshow(img, 'Parent', axImage); % 显示图像
    hold(axImage, 'on');

    % 创建用于显示信号的轴
    axSig = uiaxes('Parent', fig, 'Position', [screenSize(3)/20*11, screenSize(3)/8, screenSize(4)*0.3, screenSize(4)*0.6]);
    title(axSig, 'Signal','HorizontalAlignment', 'center');

    % 创建用于显示背景的轴
    axBg = uiaxes('Parent', fig, 'Position', [screenSize(3)/20*15, screenSize(3)/8, screenSize(4)*0.3, screenSize(4)*0.6]);
    title(axBg, 'Background','HorizontalAlignment', 'center');

    % 创建输入框标签
    bottom_edge = screenSize(4)/7;
    uilabel(fig, 'Position', [screenSize(3)/2 - 100 bottom_edge 100 30], 'Text', 'ROI No.');
    uilabel(fig, 'Position', [screenSize(3)/2 + 200 bottom_edge 120 30], 'Text', 'Inner distance (pixel)');
    uilabel(fig, 'Position', [screenSize(3)/2 + 500 bottom_edge 160 30], 'Text', 'Outter distance (pixel)');

    % 创建数值输入框
    ROInumber = uieditfield(fig, 'numeric', 'Position', [screenSize(3)/2 bottom_edge 50 30]);
    Inner_distance = uieditfield(fig, 'numeric', 'Position', [screenSize(3)/2 + 340 bottom_edge 50 30]);
    Outter_distance = uieditfield(fig, 'numeric', 'Position', [screenSize(3)/2 + 660 bottom_edge 50 30]);

    % 设置默认值
    defaultROInumber = 1;
    defaultInner_distance = 2;
    defaultOutter_distance = 5;
    ROInumber.Value = defaultROInumber;
    Inner_distance.Value = defaultInner_distance;
    Outter_distance.Value = defaultOutter_distance;

    % 绘制初始ROI
    h = drawpolygon('Position', rois.position{defaultROInumber}, 'Parent', axImage);
    setappdata(fig, 'h', h);
    setappdata(fig, 'rois', rois);
    setappdata(fig, 'background', background);

    % 设置数值框的回调函数，更新背景信息
    ROInumber.ValueChangedFcn = @(src, event) updatebackground(src.Value, Inner_distance.Value, Outter_distance.Value, movie, img, axImage, axBg, fig);
    Inner_distance.ValueChangedFcn = @(src, event) updatebackground(ROInumber.Value, src.Value, Outter_distance.Value, movie, img, axImage, axBg, fig);
    Outter_distance.ValueChangedFcn = @(src, event) updatebackground(ROInumber.Value, Inner_distance.Value, src.Value, movie, img, axImage, axBg, fig);

    % 初次更新背景
    updatebackground(defaultROInumber, defaultInner_distance, defaultOutter_distance, movie, img, axImage, axBg, fig);

    % 创建按钮：用于绘制轨迹
    plotButton = uibutton(fig, 'Text', 'Plot trace', ...
        'Position', [screenSize(3)/2 + 150 screenSize(4)/8 150 30], ...
        'ButtonPushedFcn', @(src, event) plotbutton(fig, img, ROInumber.Value, Inner_distance.Value, Outter_distance.Value, movie, axImage, axSig, axBg, ncols, nrows));

    % 创建按钮：用于保存背景和ROI数据
    saveButton = uibutton(fig, 'Text', 'Save background', ...
        'Position', [screenSize(3)/2 + 150 screenSize(4)/8 - 50 150 30], ...
        'ButtonPushedFcn', @(src, event) saveROI(fig));
end

% 更新背景区域，使用内外扩展距离来定义背景区域
function updatebackground(nroi, innerdis, outerdis, movie, img, axImage, axBg, fig)
    h = getappdata(fig, 'h');
    rois = getappdata(fig, 'rois');

    % 显示图像
    imshow(img, 'Parent', axImage);
    hold(axImage, 'on');

    % 重新绘制多边形ROI
    h = drawpolygon('Position', rois.position{nroi}, 'Parent', axImage);
    bwmask = rois.bwmask == nroi;

    % 计算内外扩展区域，并显示外扩区域边界
    se1 = strel('disk', innerdis);
    se2 = strel('disk', outerdis);
    expandmask1 = imdilate(bwmask, se1);
    expandmask2 = imdilate(bwmask, se2);
    expandregion = expandmask2 & ~expandmask1;
    expandedBoundaries = bwboundaries(expandregion);
    
    % 绘制扩展区域的边界
    for k = 1:length(expandedBoundaries)
        expandedBoundary = expandedBoundaries{k};
        plot(expandedBoundary(:,2), expandedBoundary(:,1), 'r', 'LineWidth', 0.5, 'Parent', axImage);
    end

    % 根据掩码提取背景信号
    bg = select_by_mask(expandregion, movie, axBg);

    % 保存背景信息
    background = getappdata(fig, 'background');
    background{nroi} = bg;
    setappdata(fig, 'background', background);
    setappdata(fig, 'rois', rois);
    setappdata(fig, 'h', h);
end

% 保存ROI和背景数据到工作区
function saveROI(fig)
    background = getappdata(fig, 'background');
    rois = getappdata(fig, 'rois');

    % 将背景和ROI数据保存到工作区
    assignin('base', 'background', background);
    assignin('base', 'rois', rois);
    
    % 提示保存成功
    disp('ROI updated. Background selected.');
    h = msgbox('ROI updated. Background selected.', 'Done', 'help');
    uiwait(h);
    close(fig);
end

% 绘制ROI的信号轨迹
function plotbutton(fig, img, nroi, innerdis, outerdis, movie, axImage, axSig, axBg, ncols, nrows)
    h = getappdata(fig, 'h');
    rois = getappdata(fig, 'rois');

    % 获取多边形的位置信息
    position = h.Position;
    imshow(img, 'Parent', axImage);
    hold(axImage, 'on');
    h = drawpolygon('Position', position, 'Parent', axImage);

    % 生成ROI掩码并绘制信号
    mask = poly2mask(position(:, 1), position(:, 2), ncols, nrows);
    select_by_mask(mask, movie, axSig);

    % 更新ROI信息
    rois.position{nroi} = position;
    rois.boundary{nroi} = bwboundaries(mask);
    rois.bwmask(rois.bwmask == nroi) = 0;
    rois.bwmask(mask) = nroi;

    % 保存更新的ROI信息
    setappdata(fig, 'rois', rois);
    setappdata(fig, 'h', h);

    % 更新背景
    updatebackground(nroi, innerdis, outerdis, movie, img, axImage, axBg, fig);
end

% 通过掩码提取信号轨迹
function traces = select_by_mask(bwmask, movie, trace_axe)
    num_roi = max(bwmask(:));
    traces = [];

    % 对每个ROI区域提取信号
    for i = 1:num_roi
        roi = (bwmask == i);
        trace = mean(movie(roi, :), 1);
        traces = [traces, trace'];
        plot(trace', 'Parent', trace_axe);  % 绘制信号轨迹
    end
end
