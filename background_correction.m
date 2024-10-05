function fig = background_correction(movie, ncols, nrows, rois)
% BACKGROUND_CORRECTION 交互式工具用于视频数据中的背景校正
%
%   FIG = BACKGROUND_CORRECTION(MOVIE, NCOLS, NROWS, ROIS)
%   创建一个用户界面，允许用户选择和调整感兴趣区域（ROI）的
%   背景校正。
%
%   输入参数：
%     MOVIE  - 视频数据，以 [ncols * nrows, nframes] 的格式表示。
%     NCOLS  - 图像的列数。
%     NROWS  - 图像的行数。
%     ROIS   - 结构体，包含每个ROI的位置、边界和掩码数据。
%
%   输出参数：
%     FIG    - 用于交互式背景校正的图形用户界面（UI Figure 对象）。
%
%   功能：
%     - 显示视频的平均图像。
%     - 允许用户选择ROI并定义背景区域。
%     - 实时更新背景和信号轨迹。
%     - 用户可以调整ROI的内外边界，以定义背景区域。
%
%   示例：
%     % 假设有视频数据 MOVIE 和 ROI 结构 ROIS
%     fig = background_correction(movie, 512, 512, rois);
%
%   参见：PLOT, UIAXES, UIFIGURE

% 将movie重塑为3D图像并计算每帧的均值作为静态图像
    movie2D = mean(reshape(movie, ncols, nrows, []), 3);
    img = (movie2D - min(movie2D(:))) / (max(movie2D(:)) - min(movie2D(:))); % 标准化图像
    num_rois = max(rois.bwmask(:));

    % 初始化背景单元格数组，用于存储每个ROI的背景信息
    background = zeros(size(movie,2),num_rois); 
    background_fitted = zeros(size(movie,2),num_rois); 
    background_mask = zeros(size(rois.bwmask));
    traces_bgcorr = zeros(size(movie,2),num_rois); 
    traces_bgfitcorr = zeros(size(movie,2),num_rois); 

    % 创建主窗口界面
    screenSize = get(0, 'ScreenSize');
    fig = uifigure('Position', [0 0 screenSize(3) screenSize(4)], 'Name', 'Background Correction');

    % 创建用于显示原始图像的轴
    axImage = uiaxes('Parent', fig, 'Position', [screenSize(3)/20, screenSize(3)/20, screenSize(4)*0.8, screenSize(4)*0.8]);
    title(axImage, 'Origin Image','HorizontalAlignment', 'center');
    imshow(img, 'Parent', axImage); % 显示图像
    hold(axImage, 'on');

    % 创建用于显示信号和背景的轴
    axSig = uiaxes('Parent', fig, 'Position', [screenSize(3)/20*10, screenSize(3)/8*2, screenSize(4)*0.7, screenSize(4)*0.2]);
    title(axSig, 'Signal','HorizontalAlignment', 'center');
    axBg = uiaxes('Parent', fig, 'Position', [screenSize(3)/20*10, screenSize(3)/8*3, screenSize(4)*0.7, screenSize(4)*0.2]);
    title(axBg, 'Background','HorizontalAlignment', 'center');
    axSigcorr = uiaxes('Parent', fig, 'Position', [screenSize(3)/20*10, screenSize(3)/8, screenSize(4)*0.7, screenSize(4)*0.2]);
    title(axSigcorr, 'Corrected Signal', 'HorizontalAlignment', 'center');

    % 创建输入框标签和数值输入框
    bottom_edge = screenSize(4)/7;
    uilabel(fig, 'Position', [screenSize(3)/2 - 100 bottom_edge 100 30], 'Text', sprintf('ROI No.(max %d)',num_rois));
    uilabel(fig, 'Position', [screenSize(3)/2 + 200 bottom_edge 120 30], 'Text', 'Inner distance (pixel)');
    uilabel(fig, 'Position', [screenSize(3)/2 + 500 bottom_edge 160 30], 'Text', 'Outer distance (pixel)');
    uilabel(fig, 'Position', [screenSize(3)/2 + 800 bottom_edge 160 30], 'Text', 'background threshold');
    ROInumber = uieditfield(fig, 'numeric', 'Position', [screenSize(3)/2 bottom_edge 50 30]);
    Inner_distance = uieditfield(fig, 'numeric', 'Position', [screenSize(3)/2 + 340 bottom_edge 50 30]);
    Outter_distance = uieditfield(fig, 'numeric', 'Position', [screenSize(3)/2 + 660 bottom_edge 50 30]);
    bgthreshold = uieditfield(fig, 'numeric', 'Position', [screenSize(3)/2 + 980 bottom_edge 50 30]);

    % 设置默认值
    ROInumber.Value = 1;
    lastroi = ROInumber.Value;
    Inner_distance.Value = 2;
    Outter_distance.Value = 32;
    bgthreshold.Value = 0.05;

    % 绘制初始ROI
    h = drawpolygon('Position', rois.position{1}, 'Parent', axImage);

    % 将变量存储在 fig 中
    setappdata(fig, 'movie', movie);
    setappdata(fig, 'ncols', ncols);
    setappdata(fig, 'nrows', nrows);
    setappdata(fig, 'img', img);
    setappdata(fig, 'axImage', axImage);
    setappdata(fig, 'axSig', axSig);
    setappdata(fig, 'axBg', axBg);
    setappdata(fig, 'axSigcorr', axSigcorr);
    setappdata(fig, 'rois', rois);
    setappdata(fig, 'background', background);
    setappdata(fig, 'background_fitted', background_fitted);
    setappdata(fig, 'traces_bgcorr', traces_bgcorr);
    setappdata(fig, 'traces_bgfitcorr', traces_bgfitcorr);
    setappdata(fig, 'h', h);
    setappdata(fig, 'num_rois', num_rois);
    setappdata(fig,'background_mask',background_mask);
    setappdata(fig, 'lastroi', lastroi)

    % 设置数值框的回调函数，更新背景信息
    % ROInumber.ValueChangedFcn = @(src, event) updatebackground(fig, src.Value, Inner_distance.Value, Outter_distance.Value, bgthreshold.Value);
    ROInumber.ValueChangedFcn = @(src, event) plotbutton(fig, src.Value, Inner_distance.Value, Outter_distance.Value, bgthreshold.Value);
    Inner_distance.ValueChangedFcn = @(src, event) updatebackground(fig, ROInumber.Value, src.Value, Outter_distance.Value, bgthreshold.Value);
    Outter_distance.ValueChangedFcn = @(src, event) updatebackground(fig, ROInumber.Value, Inner_distance.Value, src.Value, bgthreshold.Value);
    bgthreshold.ValueChangedFcn = @(src, event) updatebackground(fig, ROInumber.Value, Inner_distance.Value, Outter_distance.Value, src.Value);
    
    % 初次更新背景
    plotbutton(fig, ROInumber.Value, Inner_distance.Value, Outter_distance.Value, bgthreshold.Value);

    % updatebackground(fig, ROInumber.Value, Inner_distance.Value, Outter_distance.Value);

    % 创建按钮：用于绘制轨迹
    plotButton = uibutton(fig, 'Text', 'Plot trace', ...
        'Position', [screenSize(3)/2 + 150 screenSize(4)/8 - 50 150 30], ...
        'ButtonPushedFcn', @(src, event) plotbutton(fig, ROInumber.Value, Inner_distance.Value, Outter_distance.Value, bgthreshold.Value));

    % 创建按钮：用于保存背景和ROI数据
    saveButton = uibutton(fig, 'Text', 'Save background', ...
        'Position', [screenSize(3)/2 + 500 screenSize(4)/8 - 50 150 30], ...
        'ButtonPushedFcn', @(src, event) saveROI(fig, Inner_distance.Value, Outter_distance.Value, bgthreshold.Value));
end

% 更新背景区域，使用内外扩展距离来定义背景区域
function updatebackground(fig, nroi, innerdis, outerdis, bgthreshold)
    movie = getappdata(fig, 'movie');
    img = getappdata(fig, 'img');
    axImage = getappdata(fig, 'axImage');
    axBg = getappdata(fig, 'axBg');
    axSigcorr = getappdata(fig, 'axSigcorr');
    h = getappdata(fig, 'h');
    background = getappdata(fig, 'background');
    background_fitted = getappdata(fig, 'background_fitted');
    traces_bgcorr = getappdata(fig, 'traces_bgcorr');
    traces_bgfitcorr = getappdata(fig, 'traces_bgfitcorr');
    rois = getappdata(fig, 'rois');
    num_rois = getappdata(fig, 'num_rois');
    % sg = getappdata(fig, 'sg');
    background_mask = getappdata(fig,'background_mask');
    lastroi = getappdata(fig, 'lastroi');
    axSig = getappdata(fig, 'axSig');

    % 确认合法性
    if nroi > num_rois
        errordlg('Exceed max ROI number!', 'Error');
        return;
    elseif outerdis <= innerdis
        errordlg('Outer distance must be larger than Inner distance!', 'Error');
        return;
    elseif ~(0<bgthreshold && bgthreshold<=1)
        errordlg('Background threshold must bewteen 0 and 1!', 'Error');
        return;
    end

    % 显示图像
    position = h.Position;
    imshow(img, 'Parent', axImage);
    hold(axImage, 'on');

    % 重新绘制多边形ROI
    h = drawpolygon('Position', position, 'Parent', axImage);

    % 计算并显示背景区域
    se1 = strel('disk', innerdis);
    se2 = strel('disk', outerdis);
    w = msgbox('Rois updating, please wait...', 'Process', 'help');
    for i = 1:num_rois
        if nroi == lastroi || i == nroi
        bwmask = rois.bwmask == i;
        expandmask1 = imdilate(bwmask, se1);
        expandmask2 = imdilate(bwmask, se2);
        expandregion = expandmask2 & ~expandmask1;
        expandedBoundaries = bwboundaries(expandregion);

        % 更新背景掩码
        background_mask(background_mask == nroi) = 0;
        background_mask(expandregion) = nroi;

        % 根据掩码提取信号
        [bg, bg_mask] = selectbg_by_mask(expandregion, movie,bgthreshold);
        bg_fitted = polyval(polyfit(1:length(bg), bg, 1), 1:length(bg));

        sg = select_by_mask(bwmask, movie);
        
        % 保存背景信息
        background(:, i) = bg;
        background_fitted(:, i) = bg_fitted;
        traces_bgcorr(:, i) = sg-bg; 
        traces_bgfitcorr(:, i) = sg-bg_fitted; 
        max_dff = (min(sg-bg_fitted)-mean(sg-bg_fitted))/mean(sg-bg_fitted);
        
        % 绘制区域
        hold(axImage,'on')

        for k = 1:length(expandedBoundaries)
            expandedBoundary = expandedBoundaries{k};
            plot(expandedBoundary(:, 2), expandedBoundary(:, 1), 'r', 'LineWidth', 0.5, 'Parent', axImage);
            if k == length(expandedBoundaries)
            text(mean(expandedBoundary(:, 2)), mean(expandedBoundary(:, 1)), num2str(i), 'Color', 'w', 'FontSize', 18, 'HorizontalAlignment', 'center', 'Parent', axImage);
            end
        end

        % 绘制背景
        bgboundaries = bwboundaries(bg_mask);
        for k = 1:length(bgboundaries)
            bgboundary = bgboundaries{k};
            plot(bgboundary(:, 2), bgboundary(:, 1), 'y', 'LineWidth', 0.5, 'Parent', axImage);
        end

        hold(axImage,'off')

        % 绘制trace
        if i == nroi
            plot(sg, 'Parent', axSig);
            plot(bg, 'Parent', axBg); hold(axBg,'on');
            plot(bg_fitted,'r','LineWidth', 2, 'Parent', axBg); hold(axBg,'off');
            plot(sg-bg, 'Parent', axSigcorr);hold(axSigcorr,'on');
            plot(sg-bg_fitted,'r', 'Parent', axSigcorr);
            plot(zeros(size(bg))','r','LineWidth', 2, 'Parent', axSigcorr);hold(axSigcorr,'off');
            title(axSigcorr,sprintf('Corrected Signal, max dF/F = %f2', max_dff));
        end
        end
    end
    close(w)
    setappdata(fig, 'background', background);
    setappdata(fig, 'h', h);
    setappdata(fig,'traces_bgfitcorr',traces_bgfitcorr);

end

% 保存ROI和背景数据到工作区
function saveROI(fig, innerdis, outerdis, bgthreshold)
    rois = getappdata(fig, 'rois');
    background = getappdata(fig, 'background');
    background_fitted = getappdata(fig, 'background_fitted');
    traces_bgcorr = getappdata(fig, 'traces_bgcorr');
    traces_bgfitcorr = getappdata(fig, 'traces_bgfitcorr');
    background_mask = getappdata(fig,'background_mask');

    % 将背景和ROI数据保存到工作区
    assignin('base', 'background', background);
    assignin('base', 'background_fitted', background_fitted);
    assignin('base', 'background_mask', background_mask);
    assignin('base', 'traces_bgcorr', traces_bgcorr); 
    assignin('base', 'traces_bgfitcorr', traces_bgfitcorr); 
    assignin('base', 'rois_corrected', rois);
    assignin('base', 'innerdis', innerdis);
    assignin('base', 'outerdis', outerdis);
    assignin('base', 'bgthreshold', bgthreshold);

    % 提示保存成功
    disp('ROI updated. Background selected.');
    msgbox('ROI updated. Background selected.', 'Done', 'help');
end

% 绘制ROI的信号轨迹
function plotbutton(fig, nroi, innerdis, outerdis, bgthreshold)
    h = getappdata(fig, 'h');
    rois = getappdata(fig, 'rois');
    % movie = getappdata(fig, 'movie');
    img = getappdata(fig, 'img');
    axImage = getappdata(fig, 'axImage');
    ncols = getappdata(fig, 'ncols');
    nrows = getappdata(fig, 'nrows');
    lastroi = getappdata(fig, 'lastroi');

    % 获取多边形的位置信息
    if nroi == lastroi
        position = h.Position;
    else
        position = rois.position{nroi};
    end
    imshow(img, 'Parent', axImage);
    hold(axImage, 'on');
    h = drawpolygon('Position', position, 'Parent', axImage);

    % 生成ROI掩码并绘制信号
    mask = poly2mask(position(:, 1), position(:, 2), ncols, nrows);

    % 更新ROI信息
    if nroi == lastroi
    rois.boundary = bwboundaries(mask);
    rois.bwmask(rois.bwmask == nroi) = 0;
    rois.bwmask(mask) = nroi;
    rois.position{nroi} = position;
    end

    % 保存更新的ROI信息
    setappdata(fig, 'rois', rois);
    setappdata(fig, 'h', h);
    % setappdata(fig, 'sg', sg);

    % 更新背景
    updatebackground(fig, nroi, innerdis, outerdis, bgthreshold);
    setappdata(fig, 'lastroi', nroi);
end

% 通过掩码提取信号轨迹
function trace = select_by_mask(bwmask, movie)
    trace = mean(movie(bwmask(:), :), 1);
end

function [trace, bg_mask] = selectbg_by_mask(bwmask, movie, bgthreshold)
    % Step 1: 计算掩膜区域在时间维度上的平均图像
    avg_image = zeros(size(bwmask));  % 初始化一个与bwmask大小相同的图像
    avg_image(bwmask) = mean(movie(bwmask(:), :), 2);  % 计算bwmask内像素的平均值

    % Step 2: 找到bwmask掩膜内平均图像中后5%的像素值
    % 提取bwmask区域内的像素值
    masked_pixels = avg_image(bwmask);  % 提取bwmask内像素的平均值
    sorted_pixels = sort(masked_pixels);  % 对这些像素进行排序
    threshold_index = round(length(sorted_pixels) * bgthreshold);  % 找到后5%的阈值索引
    threshold_value = sorted_pixels(threshold_index);  % 阈值对应的像素值

    % Step 3: 生成背景掩膜（在bwmask区域内的低5%像素作为背景）
    bg_mask = (avg_image <= threshold_value) & bwmask;  % 背景掩膜同时考虑bwmask
    
    % Step 4: 根据背景掩膜计算背景的时间轨迹
    % 只在bg_mask内的像素进行背景时间变化轨迹的计算
    trace = mean(movie(bg_mask(:), :), 1);  % 在时间维度上平均背景区域的像素
end


