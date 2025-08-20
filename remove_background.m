function [background, background_fitted, traces_bgcorr, traces_bgfitcorr, background_mask]...
    = remove_background(movie, ncols, nrows, rois, inner_distance, outer_distance, bg_threshold)
% BACKGROUND_CORRECTION 非交互式工具用于视频数据中的背景校正
%
%   [BACKGROUND, BACKGROUND_FITTED, TRACES_BGCORR, TRACES_BGFITCORR, BACKGROUND_MASK]
%   = BACKGROUND_CORRECTION(MOVIE, NCOLS, NROWS, ROIS, INNER_DISTANCE, OUTER_DISTANCE, BG_THRESHOLD)
%   直接计算每个ROI的背景校正。
%
%   输入参数：
%     MOVIE          - 视频数据，以 [ncols * nrows, nframes] 的格式表示。
%     NCOLS          - 图像的列数。
%     NROWS          - 图像的行数。
%     ROIS           - 结构体，包含每个ROI的位置、边界和掩码数据。
%     INNER_DISTANCE - 定义背景区域的内部距离（像素）。
%     OUTER_DISTANCE - 定义背景区域的外部距离（像素）。
%     BG_THRESHOLD   - 背景阈值，用于选择背景区域的低值像素。
%
%   输出参数：
%     BACKGROUND          - 每帧的背景轨迹 [nframes, num_rois]。
%     BACKGROUND_FITTED   - 拟合的背景轨迹 [nframes, num_rois]。
%     TRACES_BGCORR       - 背景校正后的信号轨迹 [nframes, num_rois]。
%     TRACES_BGFITCORR    - 使用拟合背景校正后的信号轨迹 [nframes, num_rois]。
%     BACKGROUND_MASK     - 每个像素的背景归属掩码。

% 
if nargin <5
    inner_distance = 2;
    outer_distance = 32;
    bg_threshold = 0.05;
end

% 预分配输出变量
num_rois = max(rois.bwmask(:));
nframes = size(movie, 2);

background = zeros(nframes, num_rois);
background_fitted = zeros(nframes, num_rois);
traces_bgcorr = zeros(nframes, num_rois);
traces_bgfitcorr = zeros(nframes, num_rois);
background_mask = zeros(nrows, ncols);

% 验证输入
if outer_distance <= inner_distance
    error('Outer distance must be larger than Inner distance.');
end
if ~(0 < bg_threshold && bg_threshold <= 1)
    error('Background threshold must be between 0 and 1.');
end

% 处理每个ROI
se1 = strel('disk', inner_distance);
se2 = strel('disk', outer_distance);

for i = 1:num_rois
    % 获取ROI掩码
    bwmask = rois.bwmask == i;

    % 扩展掩码以定义背景区域
    expandmask1 = imdilate(bwmask, se1);
    expandmask2 = imdilate(bwmask, se2);
    expandregion = expandmask2 & ~expandmask1;

    % 更新背景掩码
    background_mask(expandregion) = i;

    % 根据掩码提取信号
    [bg, bg_mask] = selectbg_by_mask(expandregion, movie, bg_threshold);

    % 拟合背景轨迹
    padlength = round(0.05 * length(bg));
    padded_traces = [repmat(bg(1), padlength, 1)', bg, repmat(bg(end), padlength, 1)'];
    filtered_traces = lowpass(padded_traces, 1/10, 400); % 低通滤波
    bg_fitted = filtered_traces(padlength + 1:end - padlength);

    % 提取原始信号
    sg = select_by_mask(bwmask, movie);

    % 保存结果
    background(:, i) = bg;
    background_fitted(:, i) = bg_fitted;
    traces_bgcorr(:, i) = sg - bg;
    traces_bgfitcorr(:, i) = sg - bg_fitted;
end
fprintf('Background removed\n')
end

function trace = select_by_mask(bwmask, movie)
% 提取ROI区域的信号轨迹
trace = mean(movie(bwmask(:), :), 1);
end

function [trace, bg_mask] = selectbg_by_mask(bwmask, movie, bg_threshold)
% 提取背景区域的信号轨迹
avg_image = zeros(size(bwmask));
avg_image(bwmask) = mean(movie(bwmask(:), :), 2);

masked_pixels = avg_image(bwmask);
sorted_pixels = sort(masked_pixels);
threshold_index = round(length(sorted_pixels) * bg_threshold);
threshold_value = sorted_pixels(threshold_index);

bg_mask = (avg_image <= threshold_value) & bwmask;
trace = mean(movie(bg_mask(:), :), 1);
end
