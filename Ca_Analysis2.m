% AP_ANALYSIS_Auto - Quick Analysis of Voltage Signals
% This script is used for quickly analyzing voltage signals. It requires
% several functions and toolboxes to run properly.
%
% Requirements:
%   Functions: calculate_firing_rate, calculate_FWHM, create_map, calculate_SNR,
%              fit_exp1, highpassfilter, select_ROI
%   Toolboxes: Image Processing Toolbox, Curve Fitting Toolbox, Signal Processing Toolbox
%
% Usage:
%   1. Set the folder_path and file variables to point to your data.
%   2. Specify the frame rate (freq) of your data.

%   3. Run the script.
%
% Example:
%   folder_path = 'F:\20240321\';
%   file = '20240321.tif'; % Include the format in the file name.
%   freq = 400; % Frequency in Hz
%
% Author: Liu-Yang Luorong

% Version: 5.0
% Date: 2025.04.20
% GitHub: https://github.com/Luorong37/AP_analysis
%
% See also calculate_firing_rate, calculate_FWHM, create_map
% , calculate_SNR, fit_exp1, highpassfilter, select_ROI


clear; clc;
%## 写一个画图的段落
%% Loading raw data
nowtime = string(datetime( 'now'));
% Replace colons with hyphens to get the desired output format
nowtime = strrep(nowtime , ':', '-');
fprintf('Loading...\n')


% ↓↓↓↓↓-----------Prompt user for define path-----------↓↓↓↓↓
% support for folder, .tif, .tiff, .bin.
folder_path = 'E:\1_Data\Luorong\26.04.04 dual-colorC192\Methods8_default\Rec1_2026-04-04_22-38-58\Cycle3\\';
file = '\Cam1_Cyan5%_dgod2+1';  % must add format.do not add '\' at last
bin = 1;
downsample = 40;% downsample ratio
transpose = 1;
% ↓↓↓↓↓-----------Prompt user for frame rate------------↓↓↓↓↓
freq = 400; % Hz
gpu = true; % defined gpu open

if exist('movie','var')
    matim =true;
    [folder_path,file,exten] = fileparts(save_file);
    file = [file,exten];
else
    matim = false; % imaging via matlab
end
% loadtype = "mat"; % default "tif" = loading tif file, or "mat" = loading .mat file
% mat_name = 'LEDcyan_rec1_cycles1.mat';     % if loading .mat file

% -----------------------------------------------------------

% Create an analysis folder

file_path = fullfile(folder_path, file);
if isfolder(file_path)
    file_name = file;
    %file_dir = dir(file_path);
    %[~, ~, file_extension] = fileparts(file_dir(3).name);
else
    % [~, file_name, file_extension] = fileparts(file_path)
    [folder_path, file_name, ~] = fileparts(file_path);
end

% create a folder for analysis
save_path = fullfile(folder_path, strcat(file_name, '_Analysis'), nowtime);
mkdir(save_path);

% Load image file
if gpu
    gcp;
end

if ~matim
    fprintf("Start loading movie, please wait...\n");
    tload = tic;
    % if loadtype == "tif"
    [movie, ncols, nrows, nframes] = load_movie(file_path);
    % end
    % if loadtype == "mat"
    %     mat_path = fullfile(file_path, mat_name);
    %     load(mat_path);
    %     [ncols, nrows, nframes] = size(movie);
    % end

    movie = reshape(movie,ncols,nrows, []);

    % 翻转钙的图像
    if transpose
    movie = pagetranspose(movie);
    end
    fprintf("Finished loading and transpose movie after %d s\n", toc(tload));
else
    [ncols, nrows, nframes] = size(movie);
    movie = reshape(movie,ncols*nrows, []);
end

% Presetting
function [dt, colors, t, map, mask, options] = presetting(freq, nframes, movie, ncols, nrows)

% Define parameters
dt = 1 / freq; % Calculate time axis
% colors = [lines(7);hsv(5);spring(3);winter(3);gray(3)];
colors = lines(100);
t = (1:nframes) * dt;
options.colors = colors;
map = [];
mask = [];

end

try
    movie_vol_2D = reshape(mean(movie,2), ncols, nrows, []);
catch ME
    movie_vol_2D = squeeze(mean(movie,3));
end

avg_image  = (movie_vol_2D - min(movie_vol_2D(:))) / (max(movie_vol_2D(:)) - min(movie_vol_2D(:)));

[dt, colors, t, map, mask, options]= presetting(freq, nframes, movie, ncols, nrows);
x = (1:nframes)' * dt;

% Save code
code_path = fullfile(save_path,'Code');
mkdir(code_path);
currentScript = which("Ca_Analysis2.m");
% 获取当前脚本依赖的所有文件
[requiredFiles, ~] = matlab.codetools.requiredFilesAndProducts(currentScript);
% 复制当前脚本和所有依赖文件到目标文件夹
for k = 1:length(requiredFiles)
    [~, name, ext] = fileparts(requiredFiles{k});
    copyfile(requiredFiles{k}, fullfile(code_path, [name, ext]));
end
% 提示完成
fprintf('All codes have been copied to %s\n', code_path);
%% Motion Correction (Support Save/Apply Shifts)
fprintf('Initializing Motion Correction...\n');
% --- 1. 配置参数 ---
apply_only = 1;
loadMC     = 0;
Norigid    = 0;
hp         = 1;      % 是否开启高通滤波（用于辅助估算位移）
template   = [];
dssave     = 1;      % 是否降采样保存

% 确保 movie 是 3 维
if ismatrix(movie)
    movie = reshape(movie, ncols, nrows, []);
end

% NoRMCorre 基础配置
% options_r = NoRMCorreSetParms('d1',ncols,'d2',nrows,'bin_width',200,'max_shift',10,'us_fac',30, ...
%     'grid_size',[128,128],'overlap_pre',[32,32],'mot_uf',4,'max_dev', [5,5],'iter',1,'correct_bidir',false);
options_r = NoRMCorreSetParms('d1',ncols,'d2',nrows,'bin_width',200,'max_shift',30,'us_fac',30,'iter',1,'correct_bidir',false);

% 定义保存路径
[folder_path, file_name, fext] = fileparts(file_path);
shift_res_path = fullfile(save_path, 'motion_shifts_result.mat');
params_save_path = fullfile(save_path, 'motion_correction_para.mat');

% --- 2. 核心处理逻辑 ---
if loadMC
    [Mr, ncols, nrows, ~] = load_movie(mc_path);
else
    t1 = tic;

    % --- 内存优化：数据预处理 ---
    movie = single(movie);
    % movie = movie - min(movie(:)); % 原始数据保留在内存中

    if apply_only
        % --- 功能：直接应用位移 ---
        [shift_filename, shift_foldername] = uigetfile(save_path, '选择位移文件motion_shifts_result.mat');
        if isequal(shift_filename,0), return; end
        shift_file = fullfile(shift_foldername, shift_filename);
        S = load(shift_file);
        copyfile(shift_file, params_save_path);
        fprintf(' -> Backup: Copied shifts to %s\n', params_save_path);
        fprintf(' -> Mode: Apply existing shifts...\n');
        Mr = apply_shifts(movie, S.shifts_r, options_r);
        if Norigid && isfield(S, 'shifts_nr'), Mr = apply_shifts(Mr, S.shifts_nr, S.options_nr); end
        
    else
        % --- 功能：重新估算并保存 ---
        if hp
            % 【高通滤波模式】
            fprintf(' -> High-pass mode: Filtering for better estimation...\n');
            % 内存优化：生成临时滤波数据用于“算位移”，不改变原始 movie
            % 如果文件极大导致内存不足，建议在此处调用 create_highpass_h5
            Y = create_temp_highpass(movie);
        else
            Y = movie;
        end

        fprintf(' -> Mode: Estimating shifts...\n');
        % 1. 刚体校正：用滤波后的图算位移(shifts_r)，但应用到原始 movie 上得到 Mr
        [~, shifts_r, template1] = normcorre_batch(Y, options_r);
        Mr = apply_shifts(movie, shifts_r, options_r);

                clear Y; % 估算完立即释放临时滤波数据

        % 2. 非刚体校正 (根据需要)
        shifts_nr = []; options_nr = [];
                if hp
            % 【高通滤波模式】
            fprintf(' -> High-pass mode: Filtering for better estimation...\n');
            % 内存优化：生成临时滤波数据用于“算位移”，不改变原始 movie
            % 如果文件极大导致内存不足，建议在此处调用 create_highpass_h5
            Y = create_temp_highpass(Mr);
        else
            Y = movie;
                end

        if Norigid
            options_nr = NoRMCorreSetParms('d1',ncols,'d2',nrows,'bin_width',200,'grid_size',[128,128]);
            [~, shifts_nr, ~] = normcorre_batch(Y, options_nr, template1);
            Mr = apply_shifts(Mr, shifts_nr, options_nr);

        end
        clear Y; % 估算完立即释放临时滤波数据
        % 保存
        fprintf(' -> Saving shifts and params...\n');
        save(shift_res_path, 'shifts_r', 'shifts_nr', 'template1', '-v7.3');
        save(params_save_path, 'options_r', 'options_nr', 'hp', 'Norigid');
    end

    tsub = 40;
    % 保存校正后降采样的 TIFF
    fprintf(' -> Downsampling corrected movie (tsub = 40) for saving...\n');

    % 1. 设置下采样倍数 (400Hz -> 10Hz)
    
    if dssave
        % 2. 执行下采样 (使用 NoRMCorre 自带函数)
        % 如果你做了非刚性校正，存 Mpr；否则存 Mr
        if Norigid && exist('Mpr', 'var')
            Mr_ds = downsample_data(Mpr, 'time', tsub);
        else
            Mr_ds = downsample_data(Mr, 'time', tsub);
        end

        % 3. 重新计算显示范围 (Quantile) 确保 TIFF 亮度正常
        nn_ds = quantile(Mr_ds(:), 0.0005);
        mm_ds = quantile(Mr_ds(:), 0.99995);

        % 4. 映射到 uint16 范围并保存
        % 这样做可以确保保存后的 TIFF 在普通播放器里也能看清背景
        Mr_ds = (Mr_ds - nn_ds) / (mm_ds - nn_ds) * 65535;
        Mr_ds(Mr_ds < 0) = 0;
        Mr_ds(Mr_ds > 65535) = 65535;

        % 5. 保存 TIFF
        save_name_ds = fullfile(save_path, ['motion_corrected_ds' mat2str(tsub) '.tif']);
        array2tif(uint16(Mr_ds), save_name_ds);

        % 6. 计算平均图用于后续 ROI 提取
        movie_vol_2D = mean(Mr_ds, 3);
        fprintf(' -> Downsampled movie saved. \n');
    else
        % 保存校正后的 TIFF
        save_name = fullfile(folder_path, [file_name, '_motion_correction.tif']);
        array2tif(uint16(Mr), save_name);
        movie_vol_2D = mean(Mr, 3);
        fprintf('Finished in %d s\n', round(toc(t1)));
    end
end

% --- 3. 后处理 ---
avg_image = (movie_vol_2D - min(movie_vol_2D(:))) ./ (max(movie_vol_2D(:)) - min(movie_vol_2D(:)));

% --- 辅助子函数 ---
function Y_hp = create_temp_highpass(movie)
% 内存中快速创建高通滤波版本
gSig = 7; gSiz = 17;
psf = fspecial('gaussian', round(2*gSiz), gSig);
ind = (psf >= max(psf(:,1)));
psf = psf - mean(psf(ind)); psf(~ind) = 0;
Y_hp = imfilter(movie, psf, 'symmetric');
end

movie = reshape(uint16(Mr), ncols*nrows, []);
%%

% movie = movie-uint16(mean(movie,1));
%% Create a map (optional)
t1 = tic; % Start a timer
fprintf('Creating a map...\n')
% if the map cannot figure out active cells, please large the bin.
mapbin = 4; % defined bin = 4

[quick_map] = create_map(movie, ncols, nrows, mapbin,'calcium');
map = quick_map;

% Visualize correlation coefficients as heatmap
figure()
imagesc(quick_map);
colormap turbo;
title('Sensitivity Map');
axis image;
colorbar;
fig_filename = fullfile(save_path, '0_Sensitivity_Map.fig');
png_filename = fullfile(save_path, '0_Sensitivity_Map.png');
mat_filename = fullfile(save_path, '0_Sensitivity_Map.mat');

saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');
save(mat_filename, 'map');


t2 = toc(t1); % Get the elapsed time
fprintf('Finished mask creating after %d s\n',round(t2))
%% Load Mask (optional)
methods = questdlg('Load previous data?','load data','previous ROIs','cellpose','Cancel');
switch methods
    case 'previous ROIs'
        [roi_filename,roi_foldername] = uigetfile(save_path);
        rois_data = load(fullfile(roi_foldername,roi_filename));
        try
            rois = rois_data.rois;
            mask = rois.bwmask;

            
        catch ME
            rois.bwmask = rois_data.bwmask_ca;
            mask = rois.bwmask;
            boundaries = bwboundaries(mask, 'noholes');
        
            % 取第一个检测到的连通域边界
            % 注意：bwboundaries 返回的是 [row, col]，通常需要转为 [x, y]
            current_boundary = boundaries{1};
            rois.boundary = [current_boundary(:,2), current_boundary(:,1)]; % [X, Y]

            % 2. 计算位置 (Position)
            % 使用 regionprops 提取边界框 (BoundingBox)
            % BoundingBox 格式为 [x_left, y_top, width, height]
            stats = regionprops(mask, 'BoundingBox');
            rois.Position = stats(1).BoundingBox;
        end

        % --- 新增：尺寸检查与自动缩放 ---
        [m_rows, m_cols] = size(mask);

        if m_rows ~= nrows || m_cols ~= ncols
            fprintf('检测到 Mask 尺寸 [%d, %d] 与当前图像 [%d, %d] 不符，正在缩放...\n', ...
                m_rows, m_cols, nrows, ncols);

            % 使用最近邻插值缩放，确保 mask 依然是 logical 类型
            mask = imresize(mask, [nrows, ncols], 'nearest');

            % 更新 rois 结构体中的 mask
            rois.bwmask = mask;
        end



    case 'cellpose'
        fprintf('cellpose running......\n');
        cp = cellpose(ExecutionEnvironment="gpu");
        avgdia = 25;
        % gamma_image = imadjust(map,[],[],1.2); % recommend raise the gamma factor from 1 to 4.
        gamma_image = imadjust(avg_image,[],[],1.2); % recommend raise the gamma factor from 1 to 4.
        % gamma_image = map;
        mask = segmentCells2D(cp, gamma_image , ImageCellDiameter = avgdia,FlowErrorThreshold = 0.4,CellThreshold = 0);% CellThreshold = 0, ,  FlowErrorThreshold = 0.4
        fprintf('%d cells are found.\n',max(mask(:)));
        figure()
        % 使用labeloverlay函数显示图像
        overlayImage = labeloverlay(gamma_image, mask,'Transparency', 0.6);
        % 显示结果
        imshow(overlayImage);
        title('cellpose mask');
        fig_filename = fullfile(save_path, '0_cellpose_mask.fig');
        png_filename = fullfile(save_path, '0_cellpose_mask.png');
        mat_filename = fullfile(save_path, '0_cellpose_mask.mat');

        saveas(gcf, fig_filename, 'fig');
        saveas(gcf, png_filename, 'png');
        save(mat_filename, 'mask');
end
%%

% %翻转mask
% rois.bwmask = rois.bwmask';
% for i = 1:length(rois.boundary)
%     % 假设 rois.boundary{i} 是一个 [N x 2] 的矩阵 [rows, cols]
%     rois.boundary{i} = [rois.boundary{i}(:, 2), rois.boundary{i}(:, 1)];
% end
% for i = 1:length(rois.position)
%     % 假设 rois.position{i} 是 [N x 2] 的矩阵 [x, y]
%     % 对应图像坐标即为 [cols, rows]
%     % 转置后需变为 [rows, cols]，即交换两列
%     rois.position{i} = [rois.position{i}(:, 2), rois.position{i}(:, 1)];
% end
% mask = rois.bwmask;
%% Select ROI
t1 = tic; % Start a timer
figure()
% with or without Mask and Map
% if exist('methods','var')
%     if any(strcmp(methods,{'No',''} ))
%         [rois, traces] = select_ROI(movie, ncols, nrows, mask, map);
%         nrois = max(rois.bwmask,[],'all');
%     else
%         [~, traces] = select_ROI(movie, ncols, nrows, mask, map);
%         nois = max(mask(:));
%         rois.bwmask = mask;
%     end
% else
%     [rois, traces] = select_ROI(movie, ncols, nrows, mask, map);
%     nrois = max(rois.bwmask,[],'all');
% end
movie = reshape(movie, ncols* nrows, []);
[rois, traces] = select_ROI(movie, ncols, nrows, mask, map);
nrois = max(rois.bwmask,[],'all');

try
    bwmask = rois.bwmask;
catch ME

end
traces_original = traces;

% if do not need a map, run the following code:
% map = [];
% [bwmask, traces] = select_ROI(movie, nrows, ncols, t, mask, map);

fig_filename = fullfile(save_path, '1_raw_trace.fig');
png_filename = fullfile(save_path, '1_raw_trace.png');
roi_filename = fullfile(save_path, '1_raw_ROI.mat');

saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');
save(roi_filename, 'rois','avg_image','traces');

h = msgbox('All rois saved.', 'Done', 'help');
uiwait(h);


tif_filename = fullfile(save_path, '1_Averaged.tif');

tiffile = Tiff(tif_filename, 'w');

tiffile.setTag('ImageLength', nrows);
tiffile.setTag('ImageWidth', ncols);
tiffile.setTag('Photometric', Tiff.Photometric.MinIsBlack);
tiffile.setTag('BitsPerSample', 16);
tiffile.setTag('SamplesPerPixel', 1);
tiffile.setTag('RowsPerStrip', 16);
tiffile.setTag('PlanarConfiguration', Tiff.PlanarConfiguration.Chunky);
tiffile.setTag('Compression', Tiff.Compression.None);

tiffile.setTag('Software', 'MATLAB');
tiffile.write(uint16(movie_vol_2D));
tiffile.writeDirectory();
tiffile.close();
%% Mask 全局平移校准
fprintf('-> 开始 Mask 平移校准...\n');

% 1. 让用户指定参考 ROI 编号
roi_idx = input('请输入作为平移参考的 ROI 编号 (例如 1): ');

% 2. 提取该 ROI 并计算原始质心
target_roi = (rois.bwmask == roi_idx);
if ~any(target_roi(:))
    error('找不到编号为 %d 的 ROI，请检查输入。', roi_idx);
end
s_old = regionprops(target_roi, 'Centroid');
old_centroid = s_old.Centroid; % [x, y]

% 3. 在当前图像上手动绘制新的位置
title(sprintf('请在图中为 ROI %d 画出新的位置 (Double click to finish)', roi_idx));
h_poly = drawpolygon(gca);
new_mask = createMask(h_poly);
s_new = regionprops(new_mask, 'Centroid');
new_centroid = s_new.Centroid; % [x, y]

% 4. 计算位移矢量 (dx, dy)
shift_vec = new_centroid - old_centroid;
dx = shift_vec(1);
dy = shift_vec(2);
fprintf('检测到位移: ΔX = %.2f, ΔY = %.2f\n', dx, dy);

% 5. 平移整个 bwmask
% 使用 imtranslate 实现亚像素或像素级平移
% 注意：imtranslate 默认补充 0
shifted_mask = imtranslate(rois.bwmask, [dx, dy], 'FillValues', 0);

% 如果需要保持 ROI 编号不被插值模糊（针对标签矩阵），使用 'nearest'
% 或者重新对每一个 mask 进行平移处理以保证质量
fprintf('正在更新所有 ROI 掩膜...\n');
new_bwmask = zeros(size(rois.bwmask));
num_total = max(rois.bwmask(:));

for m = 1:num_total
    temp_m = (rois.bwmask == m);
    % 对单个 mask 平移
    temp_shifted = imtranslate(temp_m, [dx, dy], 'nearest', 'FillValues', 0);
    % 合并回总 mask (如果有重叠，后者会覆盖前者)
    new_bwmask(temp_shifted > 0) = m;
end

% 6. 更新结果并重新保存
rois.bwmask = new_bwmask;
bwmask = new_bwmask;
mask = new_bwmask;

% 重新保存更新后的文件
save(roi_filename, 'rois', 'avg_image', 'traces','shift_vec', '-append'); 
fprintf('-> Mask 已平移并重新保存至: %s\n', roi_filename);

% 刷新当前绘图查看效果
hold on;
visboundaries(new_bwmask, 'Color', 'y', 'LineWidth', 1);
title('Mask 平移校准完成 (黄色为新边界)');

%close(gcf);
%% Signal Process %%
% Background Correction
fprintf('Removing Background...\n')
[background, background_fitted, traces_bgcorr, traces_bgfitcorr, background_mask]...
    = remove_background(movie, ncols, nrows, rois, freq, bin,1, 16,0.1);
roi_filename = fullfile(save_path, '1_background_ROI.mat');
save(roi_filename, 'background','background_fitted','background_mask','traces_bgcorr','traces_bgfitcorr');
fprintf('Finished\n')

% --- 补全画图功能 ---
fprintf('Generating summary plots...\n')

% 1. 初始化设置
nrois = max(rois.bwmask(:));
colors = lines(nrois); % 生成颜色矩阵，确保每个ROI有唯一颜色

figure('Name', 'Background Correction Summary', 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.7]);

% 2. 绘制左图：原始图像 + ROI边界 + 背景掩码
subplot(1, 2, 1);
movie2D_mean = mean(reshape(movie, ncols, nrows, []), 3);
imshow(movie2D_mean, [], 'InitialMagnification', 'fit');
hold on;

for i = 1:nrois
    % 获取当前 ROI 的颜色
    current_color = colors(i, :);

    % 绘制原始 ROI 边界 (使用实线)
    roi_bw = (rois.bwmask == i);
    roi_boundaries = bwboundaries(roi_bw);
    for k = 1:length(roi_boundaries)
        boundary = roi_boundaries{k};
        plot(boundary(:, 2), boundary(:, 1), 'Color', current_color, 'LineWidth', 1, 'DisplayName', ['ROI ', num2str(i)]);
    end

    % 绘制背景掩码区域 (使用点状或半透明填充)
    bg_bw = (background_mask == i);
    bg_boundaries = bwboundaries(bg_bw);
    for k = 1:length(bg_boundaries)
        boundary = bg_boundaries{k};
        plot(boundary(:, 2), boundary(:, 1), ':', 'Color', current_color, 'LineWidth', 1);
    end

    % 在 ROI 中心标序号
    stats = regionprops(roi_bw, 'Centroid');
    if ~isempty(stats)
        text(stats(1).Centroid(1), stats(1).Centroid(2), num2str(i), ...
            'Color', current_color, 'FontWeight', 'bold', 'FontSize', 12, 'HorizontalAlignment', 'center');
    end
end
title('Image: ROI (Solid) & Background (Dotted)');

% 3. 绘制右图：信号轨迹对比
subplot(1, 2, 2);
hold on;

% 为了避免图例太乱，我们先定义几个占位符用于显示图例
p1 = plot(nan, nan, 'Color', [0.7 0.7 0.7]);
p2 = plot(nan, nan, 'k--', 'LineWidth', 1);
p3 = plot(nan, nan, 'k', 'LineWidth', 1.5);

for i = 1:nrois
    current_color = colors(i, :);

    % 计算原始信号
    raw_signal = traces_bgfitcorr(:, i) + background_fitted(:, i);

    % 绘制原始信号 (浅色背景线)
    plot(t, raw_signal, 'Color', [current_color, 0.3], 'LineWidth', 0.5);

    % 绘制拟合背景 (虚线)
    plot(t, background_fitted(:, i), '--', 'Color', current_color, 'LineWidth', 1);

    % 绘制校正后信号 (深色主线)
    plot(t, traces_bgfitcorr(:, i), 'Color', current_color, 'LineWidth', 1);
end

xlabel('Times(s)');
ylabel('Intensity');
title('Traces: Raw (Faded), BG (Dash), Corrected (Solid)');
grid on;

% 添加统一图例
legend([p1, p2, p3], {'Raw Signal', 'Fitted Background', 'Corrected Signal'}, 'Location', 'northeast');

% 保存图表
saveas(gcf, fullfile(save_path, 'background_correction_summary.png'));
fprintf('Summary plot saved to: %s\n', save_path);

%% Bleaching Correction
bleachmode = 'exp2';% 'linear' 'highpass' 'exp2' 

fprintf('Correcting Bleaching (Mode: %s)...\n', bleachmode);

switch bleachmode
    case 'highpass'
        fc = 0.5/t(end); 
        [traces_corrected, baseline] = highpass_bleach_remove(traces_bgfitcorr, freq, fc);
        
    case 'linear'
        % 使用 detrend 并计算基线
        traces_corrected = detrend(traces_bgfitcorr, 1);
        baseline = traces_bgfitcorr - traces_corrected;
        
    case 'exp2'
        % 双指数拟合提取基线
        [~, baseline] = fit_exp2(traces_bgfitcorr);
        traces_corrected = traces_bgfitcorr - baseline;
end

% 统一调用绘图函数
plot_corrected(traces_corrected, traces_bgfitcorr, baseline, t, save_path);

fprintf('Finished Bleaching Correction.\n');
function plot_corrected(traces_corrected, traces, baseline, t, save_path)
    % 获取 ROI 数量
    nrois = size(traces_corrected, 2);
    
    % 创建画布
    fig = figure('Name', 'Bleaching Correction Overview', 'Color', 'w');
    set(fig, 'Position', get(0, 'Screensize'));
    
    % --- 1. 左侧：原始数据 + 拟合基线 (Stacked) ---
    ax_fit = subplot(1, 2, 1); hold on;
    % 根据数据波动自动计算垂直间距
    spacing_raw = mean(std(traces, 0, 1)) * 5; 
    
    for r = 1:nrois
        offset = (nrois - r) * spacing_raw;
        % 原始数据 (彩色/薄线)
        plot(t, traces(:, r) + offset, 'LineWidth', 0.5, 'DisplayName', 'Original');
        % 拟合基线 (红色/粗线)
        plot(t, baseline(:, r) + offset, 'r', 'LineWidth', 1.2, 'DisplayName', 'Baseline');
        
        if mod(r, 5) == 0 || r == 1 || r == nrois
            text(t(1), offset, [' ROI ', num2str(r)], 'FontSize', 8, 'FontWeight', 'bold');
        end
    end
    title(['Original Traces & Fitted Baselines (N=', num2str(nrois), ')']);
    xlabel('Time (s)'); ylabel('Stacked Magnitude');
    grid on; axis tight;

    % --- 2. 右侧：校正后的信号 (Stacked) ---
    ax_corr = subplot(1, 2, 2); hold on;
    spacing_corr = mean(std(traces_corrected, 0, 1)) * 8;
    
    for r = 1:nrois
        offset = (nrois - r) * spacing_corr;
        % 校正后的信号 (黑色)
        plot(t, traces_corrected(:, r) + offset, 'k', 'LineWidth', 0.5);
        % 零位基准线 (浅蓝色)
        line([t(1) t(end)], [offset offset], 'Color', [0.3 0.7 1, 0.5], 'LineStyle', '--');
    end
    title('Corrected Traces (Residuals)');
    xlabel('Time (s)'); ylabel('Stacked Magnitude');
    grid on; axis tight;
    
    % 联动 X 轴（缩放左侧时右侧同步）
    linkaxes([ax_fit, ax_corr], 'x');

    % 保存图片
    fig_filename = fullfile(save_path, '2_bleach_correction_stacked.fig');
    png_filename = fullfile(save_path, '2_bleach_correction_stacked.png');
    saveas(gcf, fig_filename, 'fig');
    saveas(gcf, png_filename, 'png');
end
%% Wavelet wavelet降噪
Dnmethods = 'FDR';
Dnlevel = 4;
Wavename = 'bior6.8'; % 替代墨西哥帽的离散基
fprintf('Wavelet Denoising...\n')
traces_denoised = wdenoise(traces_corrected, Dnlevel ,DenoisingMethod=Dnmethods,Wavelet = Wavename);

traces_filename = fullfile(save_path, '2_processed_traces.mat');
save(traces_filename,"traces_denoised", "traces_original", 'traces_corrected', 'baseline', 'Dnmethods','Dnlevel')
fprintf('Finished\n');
% 全 ROI 降噪效果可视化
fprintf('Plotting all ROIs comparison...\n');

nrois = size(traces_denoised, 2);
t = (1:size(traces_denoised, 1))';


% --- 堆叠对比图 (Stacked Traces) ---

figure('Name', 'Stacked ROI Comparison', 'Color', 'w', 'Position', [150, 150, 1000, 900]);
max_display = min(10, nrois);
roi_to_show = round(linspace(1, nrois, max_display)); % 均匀选取 10 个 ROI

spacing = max(traces_corrected(:)) * 0.8; % 设置垂直间距

hold on;
for i = 1:nrois

    offset = (i-1) * spacing;

    % 原始信号（灰色）
    plot(t, traces_corrected(:, i) + offset, 'Color', [0.7 0.7 0.7], 'HandleVisibility', 'off');
    % 降噪信号（彩色）
    plot(t, traces_denoised(:, i) + offset, 'LineWidth', 1);
end

set(gca, 'YTick', (0:max_display-1)*spacing, 'YTickLabel', string(roi_to_show));
title('Comparison of Selected ROIs (Stacked View)');
xlabel('Time Points'); ylabel('ROI Index');
legend('Denoised Signal');
grid on; hold off;

fprintf('Finished plotting all requested ROI comparisons.\n');
%% Calculate sensitivity and SNR

traces_sensitivity = traces_corrected./baseline;
noise = traces_corrected-traces_denoised;
traces_SNR = traces_corrected./std(noise);

traces_filename = fullfile(save_path, '3_calculated_traces.mat');
save(traces_filename,"traces_sensitivity", 'traces_corrected', 'baseline', 'noise','traces_SNR')
figure()
%plot raw
subplot(1,3,1);
title('raw');
hold on;
[~] = offset_plot(traces,t);


% plot sensitivity
subplot(1,3,2);
title('Sensitivity');
hold on;
[~] = offset_plot(traces_sensitivity,t);

% plot SNR
subplot(1,3,3);
title('SNR');
hold on;
[~] = offset_plot(traces_SNR,t);

fig_filename = fullfile(save_path, '4_SNR.fig');
png_filename = fullfile(save_path, '4_SNR.png');
trace_filename = fullfile(save_path, '4_SNR.mat');

saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');
%% move mean

traces_ds_sensitivity = movmean(traces_sensitivity,downsample);
traces_ds_SNR = movmean(traces_SNR,downsample);

traces_filename = fullfile(save_path, '3_downsample_traces.mat');
save(traces_filename,"traces_ds_sensitivity", 'traces_corrected', 'baseline', 'noise','traces_ds_SNR','downsample')
figure()
%plot raw
subplot(1,3,1);
title('raw');
hold on;
[~] = offset_plot(movmean(traces,downsample),t);


% plot sensitivity
subplot(1,3,2);
title('Sensitivity');
hold on;
[~] = offset_plot(traces_ds_sensitivity,t);

% plot SNR
subplot(1,3,3);
title('SNR');
hold on;
[~] = offset_plot(traces_ds_SNR,t);

fig_filename = fullfile(save_path, '4_SNR_downsample.fig');
png_filename = fullfile(save_path, '4_SNR_downsample.png');
trace_filename = fullfile(save_path, '4_SNR_downsample.mat');

saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');
%% plot stacked figure

[V_filename,V_foldername] = uigetfile(save_path);
V_data = load(fullfile(V_foldername,V_filename));

SNR_traces = V_data.traces_SNR;
xV =( 0:1/400:(length(SNR_traces)-1)/400)';
%% 

% 创建一个新图形窗口
linemaxroi = 3; % 每行最多绘制3个ROI
% 计算需要绘制的行数
plotlines = floor(nrois/linemaxroi);
if mod(nrois,linemaxroi) == 0 
    plotlines = plotlines;
else
    plotlines = plotlines + 1;
end

xlimit = [0,max(max(x),max(x))+1];
ylimitv = [min(SNR_traces,[],'all'),max(SNR_traces,[],'all')];
ylimitca = [min(traces_ds_SNR,[],'all'),max(traces_ds_SNR,[],'all')];

% 计算 Scale Bar 的长度 (根据数据量调整)
% 假设 x 轴单位是秒，y 轴单位是 SNR 或 荧光强度
x_bar_time = 2 ; %s
x_bar_len = x_bar_time; % 10 
y_bar_len_ca = 10; % 钙信号高度的 30%
y_bar_len_v = 10;   % 电信号高度的 30%


% 1. 预设图形属性提高渲染效率
fig = figure('Color', 'w'); 

% 预定义通用的坐标轴属性，避免在循环中重复设置字符串
axisOpts = {'XTickLabel', [], 'YTickLabel', [], 'TickLength', [0 0], ...
            'XColor', 'none', 'YColor', 'none', 'Box', 'off', 'Color', 'none'};

for i = 0:nrois-1
    % 计算当前 ROI 在网格中的位置
    row = floor(i / linemaxroi);
    col = mod(i, linemaxroi);
    
    % 计算电信号和钙信号的 subplot 索引
    % 电信号在偶数行组，钙信号在奇数行组
    ca_idx = row * 2 * linemaxroi + col + 1;
    v_idx = ca_idx + linemaxroi;

    % --- 绘制电信号 (Voltage SNR) ---
    ax_v = subplot(plotlines * 2, linemaxroi, v_idx);
    plot(xV, SNR_traces(:, i+1), 'r');
    ylim(ylimitv);
    xlim(xlimit);
    set(ax_v, 'YDir', 'reverse', axisOpts{:});

    % --- 绘制钙信号 (Calcium SNR) ---
    ax_ca = subplot(plotlines * 2, linemaxroi, ca_idx);
    plot(x, traces_ds_SNR(:, i+1), 'g');
    ylim(ylimitca);
    xlim(xlimit);
    set(ax_ca, axisOpts{:});

    % --- 绘制 Scale Bar (仅在特定位置绘制，例如第一个或最后一个) ---
    % 建议：只在第一个 ROI 或者最后一个 ROI 绘制，避免遮挡数据
    if i == nrois - 1
        % 钙信号 Scale Bar
        x_start = xlimit(2) - x_bar_len * 1.2;
        y_start_ca = ylimitca(1) + (ylimitca(2) - ylimitca(1)) * 0.1;
        
        hold(ax_ca, 'on');
        plot(ax_ca, [x_start, x_start + x_bar_len], [y_start_ca, y_start_ca], 'k', 'LineWidth', 1.5);
        plot(ax_ca, [x_start + x_bar_len, x_start + x_bar_len], [y_start_ca, y_start_ca + y_bar_len_ca], 'k', 'LineWidth', 1.5);
        text(ax_ca, x_start + x_bar_len/2, y_start_ca, [num2str(x_bar_time),' s'], 'VerticalAlignment','top','HorizontalAlignment','center', 'FontSize', 8);
        text(ax_ca, x_start + x_bar_len, y_start_ca + y_bar_len_ca/2, [' SNR_{ca} = ' ,num2str(y_bar_len_ca) ], 'HorizontalAlignment','left', 'FontSize', 8);

        % 电信号 Scale Bar
        y_start_v = ylimitv(1) + (ylimitv(2) - ylimitv(1)) * 0.1;
        hold(ax_v, 'on');
        plot(ax_v, [x_start, x_start + x_bar_len], [y_start_v, y_start_v], 'k', 'LineWidth', 1.5);
        plot(ax_v, [x_start + x_bar_len, x_start + x_bar_len], [y_start_v, y_start_v + y_bar_len_v], 'k', 'LineWidth', 1.5);
        text(ax_v, x_start + x_bar_len, y_start_v + y_bar_len_v/2, [' SNR_{v} = ' ,num2str(y_bar_len_v) ], 'HorizontalAlignment','left', 'FontSize', 8);
    end
    
    % 标注 ROI 编号 (可选：放在电信号上方)
    title(ax_v, sprintf('ROI %d', i+1), 'FontSize', 7, 'FontWeight', 'normal');
end



fig_filename = fullfile(save_path, '2_stacked_trace.fig');
png_filename = fullfile(save_path, '2_stacked_trace.png');

saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');


%% %% Accumulate Voltage signal
% 显示示例轨迹


% 累积电压信号
% set parameter
g = 20;
k = 0.8;
gthr = 0.5;
kth = 0.21;

ftest = @(x, v) x + (g * max(v - gthr, 0) - k * max(x - kth, 0)) * 1/200;

voltage_accumulated = cell(1, nrois);
calcium_normalized = zeros(size(traces_ds_SNR));

for i = 1:nrois
    % Normalize voltage and calcium signals
    volsele = normalize(-SNR_traces(:, i),'range');
    calsele = normalize(traces_ds_SNR(:, i),'range');
    
    % Initialize accumulation array
    volaccum = zeros(length(volsele), 1);

    % Accumulate voltage signal
    for j = 2:length(volsele)
        volaccum(j) = ftest(volaccum(j-1), volsele(j));
    end

    % 计算滑动窗口大小
    window_size = round(length(volaccum) / length(calsele));
    if window_size < 1
        window_size = 1;
    end
    
    % 滑动窗口平均
    volaccum_smoothed = movmean(volaccum, window_size);
    
    % 插值使得长度相同
    volaccum_resampled = interp1(linspace(1, length(volaccum), length(volaccum)), volaccum_smoothed, linspace(1, length(volaccum), length(calsele)))';

   % 归一化插值后的信号
    voltage_accumulated{i} = normalize(volaccum_resampled, 'range');
    calcium_normalized(:,i) = calsele;
end

voltage_accumulated = cell2mat(voltage_accumulated);
voltage_accumulated = movmean(voltage_accumulated,downsample);
% 保存累积电压信号,和钙对比结果
% for i = 1:length(calculated_set)
%     fig = figure;
%     hold on;
%     plot(t,voltage_accumulated{i}, 'b');
%     plot(t,voltage_accumulated{i} + 1, 'r');
%     legend({'Integral Voltage', 'Calcium'});
%     hold off;
% end

figure()

ylimitv = [min(voltage_accumulated,[],'all'),max(voltage_accumulated,[],'all')];
ylimitca = [min(calcium_normalized,[],'all'),max(calcium_normalized,[],'all')];
for i = 0: nrois-1
    index = floor(i/linemaxroi)*linemaxroi*2 + mod(i,linemaxroi) + 1;
    subplot(plotlines*2,linemaxroi,index)
    plot(voltage_accumulated(:,i+1),'r');ylim(ylimitv);
    if i ~= nrois-1
        set(gca, 'XTickLabel', []);
        set(gca, 'YTickLabel', []);
        set(gca, 'TickLength', [0 0]);
        set(gca, 'XColor', 'none'); % 隐藏 x 轴线
        set(gca, 'YColor', 'none'); % 隐藏
    end

    set(gca, 'Box', 'off');
    set(gca, 'Color', 'none'); % 设置背景为无色
    subplot(plotlines*2,linemaxroi,index+linemaxroi)
    plot(calcium_normalized(:,i+1),'g');ylim(ylimitca);
    if i ~= nrois-1
        set(gca, 'XTickLabel', []);
        set(gca, 'YTickLabel', []);
        set(gca, 'TickLength', [0 0]);
        set(gca, 'XColor', 'none'); % 隐藏 x 轴线
        set(gca, 'YColor', 'none'); % 隐藏
    else
        ylabel(sprintf('ROI %d', floor(i/2)));
    end
    set(gca, 'Box', 'off');
    set(gca, 'Color', 'none'); % 设置背景为无色
end

% axesHandles = findall(gcf, 'type', 'axes');
% for i = 1:length(axesHandles)
%     if i ~= 1 && i ~= 2
%         set(axesHandles(i), 'XTickLabel', []);
%         set(axesHandles(i), 'YTickLabel', []);
%         set(axesHandles(i), 'TickLength', [0 0]);
%         set(axesHandles(i), 'XColor', 'none'); % 隐藏 x 轴线
%         set(axesHandles(i), 'YColor', 'none'); % 隐藏 y 轴线
%         %
%     else
%         ylabel(sprintf('ROI %d', floor(i/2)));
%     end
%     set(axesHandles(i), 'Box', 'off');
%     set(axesHandles(i), 'Color', 'none'); % 设置背景为无色
% end

fig_filename = fullfile(save_path, '3_Integral_trace.fig');
png_filename = fullfile(save_path, '3_Integral_trace.png');
trace_filename = fullfile(save_path, '3_Integral_trace.mat');

saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');
save(trace_filename,'voltage_accumulated','calcium_normalized' );
%% %% %% Calcium Deconvolution and Voltage Comparison
% 使用 FOOPSI 算法对钙信号进行去卷积，并显示去卷积后的信号与电信号对比
% set parameter
deconv_method = 'foopsi'; % 可选 'foopsi', 'constrained_foopsi', 'thresholded'
deconv_type = 'ar1'; % 自动回归模型阶数
% 假设原先是自动估算，现在尝试手动减小
sn_custom = 0.005; % 根据你的信号幅度调整，试着给原估算值的 0.5 - 0.8 倍

calcium_denoised = zeros(size(traces_ds_SNR)); % 存储去噪后的c
calcium_spikes = zeros(size(traces_ds_SNR));   % 存储去卷积后的s (events)
voltage_norm_sync = zeros(size(traces_ds_SNR)); % 存储归一化后的电压信号

for i = 1:nrois
    % 1. 获取并归一化原始信号
    calsele = traces_ds_SNR(:, i);
    % 提取对应ROI的电压信号（假设SNR_traces是高频原始电压）
    volsele = -SNR_traces(:, i); 
    
    % 2. 执行去卷积
    % 注意：根据deconvolveCa的定义，输入y通常是列向量
    % [c, s, options] = deconvolveCa(calsele, 'sn', sn_custom,'type', deconv_type, 'method', deconv_method);
    [c, s, options] = deconvolveCa(calsele,'method', 'foopsi', ...
    'sn', 0.05, ...          % 降低噪声阈值，让c更贴合y
    'pars', 0.992, ...        % 允许信号衰减得稍微快一点
    'optimize_pars', 1, ...
    'optimize_smin', 1); ...   % 建议开启自动优化pars，让算法自己找最合适的衰减率
                   % 确保没有最小脉冲限制

    % 3. 存储去卷积结果并归一化以便显示
    calcium_denoised(:, i) = normalize(c, 'range');
    calcium_spikes(:, i) = normalize(s, 'range');
    
    % 4. 处理电压信号，使其长度与钙信号对齐 (仿照原代码的重采样逻辑)
    vol_norm = normalize(volsele, 'range');
    window_size = round(length(vol_norm) / length(calsele));
    if window_size < 1, window_size = 1; end
    vol_smoothed = movmean(vol_norm, window_size);
    vol_resampled = interp1(linspace(1, length(vol_norm), length(vol_norm)), ...
                            vol_smoothed, ...
                            linspace(1, length(vol_norm), length(calsele)))';
    voltage_norm_sync(:, i) = normalize(vol_resampled, 'range');
end

% 可视化显示
figure('Name', 'Deconvolution Results');
ylimit_v = [0, 1];
ylimit_c = [0, 1];

for i = 0:nrois-1
    % 计算子图索引 (仿照原代码的双行布局)
    index = floor(i/linemaxroi)*linemaxroi*2 + mod(i,linemaxroi) + 1;
    
    % --- 上排：显示归一化电压信号 (红色) ---
    subplot(plotlines*2, linemaxroi, index)
    plot(voltage_norm_sync(:, i+1), 'r'); 
    ylim(ylimit_v);
    title(['ROI ', num2str(i+1)]);
    
    if i < nrois - linemaxroi % 如果不是最后一行，隐藏坐标轴
        set(gca, 'XTickLabel', [], 'XColor', 'none');
    end
    set(gca, 'YTickLabel', [], 'YColor', 'none', 'Box', 'off', 'Color', 'none');

    % --- 下排：显示去卷积后的钙信号 (蓝色为去噪曲线，橙色为脉冲) ---
    subplot(plotlines*2, linemaxroi, index + linemaxroi)
    hold on;
    plot(calcium_denoised(:, i+1), 'b', 'LineWidth', 1); % 去噪后的曲线
    stem(calcium_spikes(:, i+1), 'Marker', 'none', 'Color', [0.85, 0.33, 0.1]); % 脉冲信号
    hold off;
    ylim(ylimit_c);
    
    if i < nrois - linemaxroi
        set(gca, 'XTickLabel', [], 'XColor', 'none');
    end
    set(gca, 'YTickLabel', [], 'YColor', 'none', 'Box', 'off', 'Color', 'none');
end

% 保存结果
fig_deconv_fn = fullfile(save_path, '4_Deconvolution_trace.fig');
png_deconv_fn = fullfile(save_path, '4_Deconvolution_trace.png');
mat_deconv_fn = fullfile(save_path, '4_Deconvolution_trace.mat');

saveas(gcf, fig_deconv_fn, 'fig');
saveas(gcf, png_deconv_fn, 'png');
save(mat_deconv_fn, 'voltage_norm_sync', 'calcium_denoised', 'calcium_spikes','options');

fprintf('去卷积计算与绘图完成。\n');


%% 统计分析
paired_corrtest = zeros(1,nrois);
paired_spcorrtest = zeros(1,nrois);
for i = 1:nrois
    paired_corrtest(i) = corr(voltage_accumulated(:,i),calcium_normalized(:,i));
    paired_spcorrtest(i) = corr(voltage_accumulated(:,i),calcium_normalized(:,i),'Type', 'Spearman');

end

random_corrtest = zeros(nrois);
random_speartest = zeros(nrois);

for i = 1:nrois
    for j = 1:nrois
        if i ~= j
            random_corrtest(i, j) = corr(voltage_accumulated(:,i), calcium_normalized(:,j));
            random_speartest(i, j) = corr(voltage_accumulated(:,i), calcium_normalized(:,j), 'Type', 'Spearman');
        end
    end
end

random_corrtest(random_corrtest == 0) = [];
random_speartest(random_speartest == 0) = [];

% Perform significance tests
[~, p_corr] = ttest2(paired_corrtest, random_corrtest);
[~, p_spcorr] = ttest2(paired_spcorrtest, random_speartest);

% Create figure and boxplots
figure;
subplot(1, 2, 1);
boxplot([paired_corrtest(:); random_corrtest(:)], [repmat({'Paired'}, size(paired_corrtest(:))); repmat({'Random'}, size(random_corrtest(:)))]);
title('Correlation');
text(1.5, max([paired_corrtest(:); random_corrtest(:)])*0.95, sprintf('p = %.3f', p_corr), 'HorizontalAlignment', 'center');

subplot(1, 2, 2);
boxplot([paired_spcorrtest(:); random_speartest(:)], [repmat({'Paired'}, size(paired_spcorrtest(:))); repmat({'Random'}, size(random_speartest(:)))]);
title('Spearman Correlation');
text(1.5, max([paired_spcorrtest(:); random_speartest(:)])*0.95, sprintf('p = %.3f', p_spcorr), 'HorizontalAlignment', 'center');

fig_filename = fullfile(save_path, '4_Correlation_analysis.fig');
png_filename = fullfile(save_path, '4_Correlation_analysis.png');
mat_filename = fullfile(save_path, '4_Correlation_analysis.mat');

saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');
save(mat_filename,'paired_corrtest','random_corrtest','paired_spcorrtest','random_speartest');
%% %% save trace for each ROI
traces_path = fullfile(save_path,'eachROI');
mkdir(traces_path);

for i = 1:nrois
    
    % plot aligned each trace of each ROI
    fig = figure();
    subplot(3,1,1);
    current_SNR = SNR_traces(:,i);
    plot(xV,current_SNR,'r');hold on;xlim(xlimit);
    ylabel('Original Voltage SNR')
    
    title(sprintf('ROI %d',i),sprintf('Correlation effector = %f',paired_corrtest(i)))

    set(gca,'YDir','reverse')

    subplot(3,1,2);
    plot(x,voltage_accumulated(:,i),'r');xlim(xlimit);
    ylabel('Integral Voltage')


    subplot(3,1,3);
    current = traces_ds_SNR(:,i);
    plot(x,current,'g');hold on;xlim(xlimit);

    ylabel('Calcium')
    xlabel('Time')
    
    % save the figure
    fig_filename = fullfile(traces_path, sprintf('ROI %d Corr. = %f.fig', i,paired_corrtest(i)));
    png_filename = fullfile(traces_path, sprintf('ROI %d Corr. = %f.png', i,paired_corrtest(i)));
    saveas(gcf, fig_filename, 'fig');
    saveas(gcf, png_filename, 'png');
    close;
end

fprintf('Traces of each ROI are saved\n')
%% %% %% save trace for each ROI (Updated with Deconvolution)
traces_path = fullfile(save_path, 'eachROI_with_Deconv');
if ~exist(traces_path, 'dir')
    mkdir(traces_path);
end

for i = 1:nrois
    % 创建画布
    fig = figure('Visible', 'off'); % 建议设置Visible为off，大批量保存时不会频繁弹窗
    
    % --- Subplot 1: 钙信号 (Calcium) ---
    subplot(4,1,1);
    current_cal = traces_ds_SNR(:,i);
    plot(x, current_cal, 'g'); 
    hold on; xlim(xlimit);
    ylabel('Calcium');
    % title(sprintf('ROI %d', i), sprintf('Correlation = %f', paired_corrtest(i)));
    
    % --- Subplot 2: 电信号积分 (Integral Voltage) ---
    subplot(4,1,2);
    plot(x, voltage_accumulated(:,i), 'r'); 
    xlim(xlimit);
    ylabel('Integral Voltage');
    
    % --- Subplot 3: 钙信号去卷积 (Deconvolved Calcium) ---
    subplot(4,1,3);
    % 绘制去卷积后的趋势线 (c) 和 脉冲信号 (s)
    hold on;
    % 绘制去噪后的曲线
    plot(x, calcium_denoised(:, i), 'b', 'LineWidth', 1); 
    % 绘制脉冲信号，使用stem或直接用线条表示
    stem(x, calcium_spikes(:, i), 'Marker', 'none', 'Color', [0.85, 0.33, 0.1]); 
    xlim(xlimit);
    ylabel('Deconv Calcium');
    hold off;
    
    % --- Subplot 4: 原始电信号 (Original Voltage SNR) ---
    subplot(4,1,4);
    current_SNR = SNR_traces(:,i);
    plot(xV, current_SNR, 'r'); 
    hold on; xlim(xlimit);
    ylabel('Original Voltage');
    xlabel('Time');
    set(gca, 'YDir', 'reverse'); % 保持你原来的电压反向显示习惯
    
    % 优化布局，防止标签重叠
    set(findall(gcf,'-property','FontSize'),'FontSize', 9);
    
    % --- 保存图片 ---
    fig_filename = fullfile(traces_path, sprintf('ROI_%d.fig', i));
    png_filename = fullfile(traces_path, sprintf('ROI_%d.png', i));
    
    saveas(gcf, fig_filename, 'fig');
    saveas(gcf, png_filename, 'png');
    close(fig);
end
fprintf('Traces (including Deconvolution) of each ROI are saved in: %s\n', traces_path);
%% Save parameter
% 定义保存路径和文件名
save_filename = fullfile(save_path, '-1_workspace_variables.mat');

% % 保存当前工作区中的所有变量到.mat文件
clear movie;
clear Mr;
save(save_filename);



%%
%% 准备数据
roi_ids = 1:nrois;

% 1. 生成整体大图 (Summary)
info_summary.path = save_path;
info_summary.filename = '5_All_ROIs_Comparison';
info_summary.mode = 'summary';

plot_ROI_quad_signals(x, xV, ...
    calcium_normalized, ...   % 钙信号
    voltage_accumulated, ...  % 电信号积分
    calcium_denoised, ...     % 去卷积曲线
    calcium_spikes, ...       % 去卷积脉冲
    SNR_traces, ...           % 原始电信号
    roi_ids, info_summary);

% 2. 生成每个 ROI 的独立详图 (Individual)
info_indiv.path = fullfile(save_path, 'eachROI_Detail');
info_indiv.mode = 'individual';
info_indiv.corr = paired_corrtest; % 传入相关系数用于标题显示

plot_ROI_quad_signals(x, xV, ...
    calcium_normalized, ...
    voltage_accumulated, ...
    calcium_denoised, ...
    calcium_spikes, ...
    SNR_traces, ...
    roi_ids, info_indiv);
function plot_ROI_quad_signals(x, xV, calcium, integral_v, deconv_c, deconv_s, raw_v, titles, save_info)
    % plot_ROI_quad_signals: 绘制四层信号对比图
    % 输入说明:
    % x, xV: 时间轴
    % calcium: 钙信号 (nxM)
    % integral_v: 积分电压 (nxM)
    % deconv_c, deconv_s: 去卷积曲线和脉冲 (nxM)
    % raw_v: 原始电压 (nxM)
    % titles: ROI 编号或名称
    % save_info: 结构体，包含 'path', 'filename', 'mode' ('summary' 或 'individual')

    nrois = size(calcium, 2);
    linemaxroi = 3; 
    plotlines = ceil(nrois / linemaxroi);
    
    if strcmp(save_info.mode, 'summary')
        % --- 模式 A: 整体概览图 (所有 ROI 堆叠在一个 Figure) ---
        fig = figure('Name', 'Summary_Trace', 'Color', 'w', 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.8]);
        axisOpts = {'XTickLabel', [], 'YTickLabel', [], 'TickLength', [0 0], 'XColor', 'none', 'YColor', 'none', 'Box', 'off'};

        for i = 1:nrois
            col = mod(i-1, linemaxroi);
            row = floor((i-1) / linemaxroi);
            % 计算基础索引（每个ROI占4行中的位置）
            % 这里为了在 summary 图中放下 4 层，subplot 的总行数需设为 plotlines*4
            base_idx = row * 4 * linemaxroi + col + 1;

            % 1. Calcium
            ax1 = subplot(plotlines * 4, linemaxroi, base_idx);
            plot(x, calcium(:,i), 'g'); title(['ROI ', num2str(titles(i))]); set(ax1, axisOpts{:});
            
            % 2. Integral V
            ax2 = subplot(plotlines * 4, linemaxroi, base_idx + linemaxroi);
            plot(x, integral_v(:,i), 'r'); set(ax2, axisOpts{:});
            
            % 3. Deconv
            ax3 = subplot(plotlines * 4, linemaxroi, base_idx + 2*linemaxroi);
            plot(x, deconv_c(:,i), 'b'); hold on;
            stem(x, deconv_s(:,i), 'Marker', 'none', 'Color', [0.85 0.33 0.1]); set(ax3, axisOpts{:});
            
            % 4. Raw Voltage
            ax4 = subplot(plotlines * 4, linemaxroi, base_idx + 3*linemaxroi);
            plot(xV, raw_v(:,i), 'r'); set(ax4, 'YDir', 'reverse', axisOpts{:});
        end
        saveas(fig, fullfile(save_info.path, [save_info.filename, '.png']), 'png');
        saveas(fig, fullfile(save_info.path, [save_info.filename, '.fig']), 'fig');
    else
        % --- 模式 B: 循环保存每个 ROI 的独立详图 ---
        if ~exist(save_info.path, 'dir'), mkdir(save_info.path); end
        for i = 1:nrois
            fig = figure('Visible', 'off', 'Color', 'w');
            
            subplot(4,1,1); plot(x, calcium(:,i), 'g'); ylabel('Calcium');
            title(sprintf('ROI %d | Corr: %.3f', titles(i), save_info.corr(i)));
            
            subplot(4,1,2); plot(x, integral_v(:,i), 'r'); ylabel('Integral V');
            
            subplot(4,1,3); plot(x, deconv_c(:,i), 'b'); hold on;
            stem(x, deconv_s(:,i), 'Marker', 'none', 'Color', [0.85 0.33 0.1]); ylabel('Deconv');
            
            subplot(4,1,4); plot(xV, raw_v(:,i), 'r'); ylabel('Raw V');
            xlabel('Time'); set(gca, 'YDir', 'reverse');
            
            saveas(fig, fullfile(save_info.path, sprintf('ROI_%d.png', titles(i))), 'png');
            close(fig);
        end
    end
end

%% 在前面电钙raw trace基础上加上刺激时间窗口 by YHY
% 读取刺激数据
if ~exist('logs', 'var')
    [logs_filename,logs_foldername] = uigetfile(save_path, '在Rec文件夹里选择日志文件 logs.mat');
    load(fullfile(logs_foldername,logs_filename));
end
if ~exist('stimCfg', 'var')
    [stimcfg_filename,stimcfg_foldername] = uigetfile(save_path, '在Method文件夹里选择刺激设定文件 stimcfg.mat');
    load(fullfile(stimcfg_foldername,stimcfg_filename));
end
%% 

% stim_frame = table2array(logs.sync(:,1));
cam_frame = table2array(logs.sync(:,4));
stim_framerate = round(1/stimCfg.ifi);
stim_cycles = stimCfg.numorien;
stim_duration_frame = stimCfg.duration * stim_framerate;
stim_rest_frame = stimCfg.isi * stim_framerate;

stim_range = zeros(stim_cycles,2);
for i = 1:stim_cycles
    stim_range(i,1) = cam_frame((i-1)*(stim_duration_frame + stim_rest_frame) + 1);
    stim_range(i,2) = cam_frame((i-1)*(stim_duration_frame + stim_rest_frame) + stim_duration_frame);  % 得到每次刺激对应的起始/结束帧（每行一组）
end

x_stim = zeros(4,stim_cycles);
y_stim = zeros(4,stim_cycles);
for i = 1:stim_cycles
    x_stim(1,i) = cam_frame((i-1)*(stim_duration_frame + stim_rest_frame) + 1);
    x_stim(4,i) = x_stim(1,i);
    x_stim(2,i) = cam_frame((i-1)*(stim_duration_frame + stim_rest_frame) + stim_duration_frame);
    x_stim(3,i) = x_stim(2,i);
end
x_stim = x_stim/400;

% 创建一个新图形窗口
linemaxroi = 3; % 每行最多绘制3个ROI
% 计算需要绘制的行数
plotlines = floor(nrois/linemaxroi);
if mod(nrois,linemaxroi) == 0 
    plotlines = plotlines;
else
    plotlines = plotlines + 1;
end

xlimit = [0,max(max(x),max(x))+1];
ylimitv = [min(SNR_traces,[],'all'),max(SNR_traces,[],'all')];
ylimitca = [min(traces_ds_SNR,[],'all'),max(traces_ds_SNR,[],'all')];

% 计算 Scale Bar 的长度 (根据数据量调整)
% 假设 x 轴单位是秒，y 轴单位是 SNR 或 荧光强度
x_bar_time = 2 ; %s
x_bar_len = x_bar_time; % 10 
y_bar_len_ca = 10; % 钙信号高度的 30%
y_bar_len_v = 10;   % 电信号高度的 30%


% 1. 预设图形属性提高渲染效率
fig = figure('Color', 'w'); 

% 预定义通用的坐标轴属性，避免在循环中重复设置字符串
axisOpts = {'XTickLabel', [], 'YTickLabel', [], 'TickLength', [0 0], ...
            'XColor', 'none', 'YColor', 'none', 'Box', 'off', 'Color', 'none'};

for i = 0:nrois-1
    % 计算当前 ROI 在网格中的位置
    row = floor(i / linemaxroi);
    col = mod(i, linemaxroi);
    
    % 计算电信号和钙信号的 subplot 索引
    % 电信号在偶数行组，钙信号在奇数行组
    ca_idx = row * 2 * linemaxroi + col + 1;
    v_idx = ca_idx + linemaxroi;

    % --- 绘制电信号 (Voltage SNR) ---
    ax_v = subplot(plotlines * 2, linemaxroi, v_idx);
    plot(xV, SNR_traces(:, i+1), 'r');
    ylim(ylimitv);
    xlim(xlimit);
    hold on;
    for j = 1:stim_cycles
        y_stim(:,j) = [ylimitv(1) ylimitv(1) ylimitv(2) ylimitv(2)]';
    end
    fill(x_stim, y_stim, 'r', 'FaceColor', "#C0E6EB", 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    set(ax_v, 'YDir', 'reverse', axisOpts{:});

    % --- 绘制钙信号 (Calcium SNR) ---
    ax_ca = subplot(plotlines * 2, linemaxroi, ca_idx);
    plot(x, traces_ds_SNR(:, i+1), 'g');
    ylim(ylimitca);
    xlim(xlimit);
    hold on;
    for j = 1:stim_cycles
        y_stim(:,j) = [ylimitca(1) ylimitca(1) ylimitca(2) ylimitca(2)]';
    end
    fill(x_stim, y_stim, 'r', 'FaceColor', "#C0E6EB", 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    set(ax_ca, axisOpts{:});

    % --- 绘制 Scale Bar (仅在特定位置绘制，例如第一个或最后一个) ---
    % 建议：只在第一个 ROI 或者最后一个 ROI 绘制，避免遮挡数据
    if i == nrois - 1
        % 钙信号 Scale Bar
        x_start = xlimit(2) - x_bar_len * 1.2;
        y_start_ca = ylimitca(1) + (ylimitca(2) - ylimitca(1)) * 0.1;
        
        hold(ax_ca, 'on');
        plot(ax_ca, [x_start, x_start + x_bar_len], [y_start_ca, y_start_ca], 'k', 'LineWidth', 1.5);
        plot(ax_ca, [x_start + x_bar_len, x_start + x_bar_len], [y_start_ca, y_start_ca + y_bar_len_ca], 'k', 'LineWidth', 1.5);
        text(ax_ca, x_start + x_bar_len/2, y_start_ca, [num2str(x_bar_time),' s'], 'VerticalAlignment','top','HorizontalAlignment','center', 'FontSize', 8);
        text(ax_ca, x_start + x_bar_len, y_start_ca + y_bar_len_ca/2, [' SNR_{ca} = ' ,num2str(y_bar_len_ca) ], 'HorizontalAlignment','left', 'FontSize', 8);

        % 电信号 Scale Bar
        y_start_v = ylimitv(1) + (ylimitv(2) - ylimitv(1)) * 0.1;
        hold(ax_v, 'on');
        plot(ax_v, [x_start, x_start + x_bar_len], [y_start_v, y_start_v], 'k', 'LineWidth', 1.5);
        plot(ax_v, [x_start + x_bar_len, x_start + x_bar_len], [y_start_v, y_start_v + y_bar_len_v], 'k', 'LineWidth', 1.5);
        text(ax_v, x_start + x_bar_len, y_start_v + y_bar_len_v/2, [' SNR_{v} = ' ,num2str(y_bar_len_v) ], 'HorizontalAlignment','left', 'FontSize', 8);
    end
    
    % 标注 ROI 编号 (可选：放在电信号上方)
    title(ax_v, sprintf('ROI %d', i+1), 'FontSize', 7, 'FontWeight', 'normal');
end



fig_filename = fullfile(save_path, '2_stacked_trace_with_stim.fig');
png_filename = fullfile(save_path, '2_stacked_trace_with_stim.png');

saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');
%% --- 1. 数据准备 (确保已计算出 resp_ca, resp_v, dsi_ca, dsi_v) ---

%% --- 1. 参数初始化与路径准备 ---
% 假设 stim_cycles = stimCfg.numorien
theta = stimCfg.orientations; 
theta_rad = deg2rad(theta);
p_theta = [theta_rad, theta_rad(1)]; % 用于极坐标闭合

% 确保保存路径存在
save_dir = fullfile(save_path, 'Comprehensive_Analysis');
if ~exist(save_dir, 'dir'), mkdir(save_dir); end

% 初始化用于显著性检验的绝对均值矩阵 [ROIs x Trials]
all_stim_ca = zeros(nrois, stim_cycles);
all_base_ca = zeros(nrois, stim_cycles);
all_stim_v  = zeros(nrois, stim_cycles);
all_base_v  = zeros(nrois, stim_cycles);

resp_ca = zeros(nrois, stim_cycles);
resp_v  = zeros(nrois, stim_cycles);

% --- 2. 基于新逻辑确定帧范围 (Stim -> ISI) ---
stim_range = zeros(stim_cycles, 2);
rest_range = zeros(stim_cycles, 2); % 新增：存放刺激之后的 ISI 帧范围

for i = 1:stim_cycles
    % 1. 刺激帧索引 (直接使用您提供的逻辑)
    idx_stim_start = (i-1)*(stim_duration_frame + stim_rest_frame) + 1;
    idx_stim_end   = (i-1)*(stim_duration_frame + stim_rest_frame) + stim_duration_frame;
    
    stim_range(i,1) = cam_frame(idx_stim_start);
    stim_range(i,2) = cam_frame(idx_stim_end);
    
    % 2. 静息(ISI)帧索引：紧接在刺激结束之后
    idx_rest_start = idx_stim_end + 1;
    idx_rest_end   = idx_rest_start + stim_rest_frame - 1;
    
    % 防止超出 cam_frame 总长度
    idx_rest_end = min(idx_rest_end, length(cam_frame));
    
    rest_range(i,1) = cam_frame(idx_rest_start);
    rest_range(i,2) = cam_frame(idx_rest_end);
end

% --- 3. 核心数据提取 ---
fprintf('正在提取响应与静息数据 (ISI在刺激之后)...\n');
for j = 1:stim_cycles
    % 取出当前刺激和其后 ISI 的相机帧
    start_f = stim_range(j, 1);
    end_f   = stim_range(j, 2);
    
    base_start_f = rest_range(j, 1);
    base_end_f   = rest_range(j, 2);
    
    for i = 1:nrois
        % --- 钙信号提取 ---
        ca_seg = traces_sensitivity(start_f:end_f, i);
        ca_base_seg = traces_sensitivity(base_start_f:base_end_f, i);
        
        all_stim_ca(i, j) = mean(ca_seg);
        all_base_ca(i, j) = mean(ca_base_seg);
        resp_ca(i, j) = all_stim_ca(i, j) - all_base_ca(i, j); % 净响应
        
        % --- 电信号提取 (取反，使兴奋向上) ---
        v_seg = -V_data.traces_sensitivity(start_f:end_f, i); 
        v_base_seg = -V_data.traces_sensitivity(base_start_f:base_end_f, i);
        
        all_stim_v(i, j) = mean(v_seg);
        all_base_v(i, j) = mean(v_base_seg);
        resp_v(i, j) = all_stim_v(i, j) - all_base_v(i, j);    % 净响应
    end
end

% 数据清洗：将负净响应归零（仅用于后续调谐计算）
resp_ca(resp_ca < 0) = 0;
resp_v(resp_v < 0) = 0;

% --- 4. 计算 DSI & OSI (向量和法) ---
fprintf('正在计算选择性指数...\n');
dsi_ca = zeros(nrois, 1); osi_ca = zeros(nrois, 1);
dsi_v  = zeros(nrois, 1); osi_v  = zeros(nrois, 1);

for i = 1:nrois
    % 钙信号
    r_c = resp_ca(i,:);
    if sum(r_c) > 0
        dsi_ca(i) = abs(sum(r_c .* exp(1j * theta_rad)) / sum(r_c));
        osi_ca(i) = abs(sum(r_c .* exp(1j * 2 * theta_rad)) / sum(r_c));
    end
    % 电信号
    r_v = resp_v(i,:);
    if sum(r_v) > 0
        dsi_v(i) = abs(sum(r_v .* exp(1j * theta_rad)) / sum(r_v));
        osi_v(i) = abs(sum(r_v .* exp(1j * 2 * theta_rad)) / sum(r_v));
    end
end

% --- 5. 批量生成并保存综合图表 ---
fprintf('正在生成 ROI 分析报告...\n');
for i = 1:nrois
    % 配对 T 检验 (比较刺激期与随后的 ISI 期)
    [~, p_ca] = ttest(all_base_ca(i,:), all_stim_ca(i,:));
    [~, p_v]  = ttest(all_base_v(i,:), all_stim_v(i,:));
    
    h_fig = figure('Color', 'w', 'Position', [100, 50, 1200, 900], 'Visible', 'off'); 
    
    % --- 子图 A: 原始轨迹图 ---
    subplot(3, 2, [1, 2]); 
    yyaxis left
    plot(xV, -V_data.traces_sensitivity(:, i), 'r', 'LineWidth', 0.5); 
    ylabel('Inverted Volt SNR'); set(gca, 'YColor', 'r');
    hold on;
    yyaxis right
    plot(x, traces_sensitivity(:, i), 'g', 'LineWidth', 1.2); 
    ylabel('Calcium SNR'); set(gca, 'YColor', 'g');
    
    % 绘制刺激阴影区 (根据 stim_range 精确绘制)
    yl = ylim;
    for j = 1:stim_cycles
        % 将帧转换为时间作图 (假设 x 和 xV 对应)
        t_start = x(stim_range(j, 1));
        t_end   = x(stim_range(j, 2));
        fill([t_start t_start t_end t_end], [yl(1) yl(2) yl(2) yl(1)], 'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    end
    title(['ROI ' num2str(i) ' Raw Traces (Blue shaded: Stimulus ON)']);
    xlabel('Time (s)'); grid on;

    % --- 子图 B: 极坐标调谐图 ---
    subplot(3, 2, 3);
    p_data_ca = [resp_ca(i,:), resp_ca(i,1)];
    p_data_v  = [resp_v(i,:), resp_v(i,1)];
    polarplot(p_theta, p_data_ca/max(p_data_ca+eps), 'g-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'g'); hold on;
    polarplot(p_theta, p_data_v/max(p_data_v+eps), 'r-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'r');
    title('Direction Preference (Norm.)');
    legend({'Calcium', 'Voltage (Inv)'}, 'Location', 'northeastoutside', 'FontSize', 8);

    % --- 子图 C: 调谐曲线与指数 ---
    subplot(3, 2, 4);
    plot(theta, resp_ca(i,:)/max(resp_ca(i,:)+eps), 'g-s', 'LineWidth', 1.2); hold on;
    plot(theta, resp_v(i,:)/max(resp_v(i,:)+eps), 'r-s', 'LineWidth', 1.2);
    
    txt = {sprintf('DSI (Ca): %.2f', dsi_ca(i)), ...
           sprintf('OSI (Ca): %.2f', osi_ca(i)), ...
           sprintf('DSI (V):  %.2f', dsi_v(i)), ...
           sprintf('OSI (V):  %.2f', osi_v(i))};
    text(0.05, 0.95, txt, 'Units', 'normalized', 'VerticalAlignment', 'top', 'FontSize', 9, 'EdgeColor', 'k', 'BackgroundColor', 'w');
    xlabel('Direction (deg)'); ylabel('Norm. Net Response');
    title('Orientation Tuning Curve');
    xticks(theta); xtickangle(45); grid on;
    
    % --- 子图 D & E: 显著性检验 (Rest vs Stim) ---
    % 钙信号
    subplot(3, 2, 5);
    ca_comp_data = [all_base_ca(i,:)', all_stim_ca(i,:)'];
    plot([1, 2], ca_comp_data', 'Color', [0.7 0.9 0.7], 'Marker', 'o', 'MarkerFaceColor', 'g'); hold on;
    errorbar([1, 2], mean(ca_comp_data), std(ca_comp_data)/sqrt(stim_cycles), 'k', 'LineWidth', 2);
    xlim([0.5 2.5]); xticks([1, 2]); xticklabels({'ISI (Post-Stim)', 'Stim'});
    ylabel('Mean SNR'); 
    sig_str = sprintf('P = %.3f', p_ca);
    if p_ca < 0.05, sig_str = [sig_str ' *']; end; if p_ca < 0.01, sig_str = [sig_str '*']; end
    title(['Calcium Responsiveness (', sig_str, ')']); grid on;

    % 电信号
    subplot(3, 2, 6);
    v_comp_data = [all_base_v(i,:)', all_stim_v(i,:)'];
    plot([1, 2], v_comp_data', 'Color', [0.9 0.7 0.7], 'Marker', 'o', 'MarkerFaceColor', 'r'); hold on;
    errorbar([1, 2], mean(v_comp_data), std(v_comp_data)/sqrt(stim_cycles), 'k', 'LineWidth', 2);
    xlim([0.5 2.5]); xticks([1, 2]); xticklabels({'ISI (Post-Stim)', 'Stim'});
    ylabel('Mean Inverted SNR'); 
    sig_str_v = sprintf('P = %.3f', p_v);
    if p_v < 0.05, sig_str_v = [sig_str_v ' *']; end; if p_v < 0.01, sig_str_v = [sig_str_v '*']; end
    title(['Voltage Responsiveness (', sig_str_v, ')']); grid on;
    
    % 保存并关闭
    save_name = fullfile(save_dir, ['ROI_' num2str(i) '_Analysis.png']);
    saveas(h_fig, save_name);
    close(h_fig); 
end

fprintf('所有 %d 个 ROI 处理完毕！图片已保存至: %s\n', nrois, save_dir);
%% 
%% 高级调谐曲线分析 (基于 Li et al., 2008 & Niell & Stryker, 2008)

%% --- 1. 参数与窗口定义 ---
num_oris = stimCfg.numorien;
theta = stimCfg.orientations; 
theta_rad = deg2rad(theta);
hanning_kernel = [0.5 1 0.5]; % 文献指定的平滑核

% 定义响应计算的帧范围（第2到第8帧）
frame_range = stim_rest_frame;

% 初始化结果
resp_df_f = zeros(nrois, num_oris); % 存储每个方向的平均 Delta F/F
di_index = zeros(nrois, 1);        % 方向指数

% --- 2. 按照文献算法提取响应 ---
fprintf('正在执行 Delta F/F 向量化提取...\n');
for j = 1:num_oris
    % 获取刺激和参考期的帧索引
    s_start = stim_range(j,1);
    r_start = rest_range(j,1); % 假设 rest 在 stim 之后（按你之前的要求）
    
    for i = 1:nrois
        % 提取第 2 到 8 帧
        f_stim = mean(traces_sensitivity(s_start + frame_range - 1, i));
        f_ref  = mean(traces_sensitivity(r_start + frame_range - 1, i));
        
        % 计算百分比增加量 (Equation 1)
        resp_df_f(i, j) = ((f_stim - f_ref) / f_ref) * 100;
    end
end

% --- 3. 平滑与方向指数计算 ---
fprintf('正在进行 Hanning 平滑与 DI 计算...\n');
smoothed_resp = zeros(size(resp_df_f));
for i = 1:nrois
    % 循环平滑 (考虑到 0-360 是周期的)
    raw_curve = resp_df_f(i, :);
    extended_curve = [raw_curve(end), raw_curve, raw_curve(1)];
    temp_smooth = conv(extended_curve, hanning_kernel, 'valid') / sum(hanning_kernel);
    smoothed_resp(i, :) = temp_smooth;
    
    % 寻找最优方向 (Preferred Direction)
    [r_pref, pref_idx] = max(smoothed_resp(i, :));
    theta_pref = theta(pref_idx);
    
    % 寻找相反方向 (Opposite Direction, +180度)
    opp_angle = mod(theta_pref + 180, 360);
    [~, opp_idx] = min(abs(theta - opp_angle));
    r_opp = smoothed_resp(i, opp_idx);
    
    % 计算方向指数 (Direction Index, Eq. 5)
    di_index(i) = (r_pref - r_opp) / (r_pref + r_opp + eps);
end

% --- 4. 向量空间计算 (DSI/OSI 矢量法) ---
dsi_vect = zeros(nrois, 1);
osi_vect = zeros(nrois, 1);
for i = 1:nrois
    r = resp_df_f(i, :);
    % 方向向量 (Eq. 3)
    dsi_vect(i) = abs(sum(r .* exp(1j * theta_rad)) / sum(r));
    % 取向向量 (Eq. 2, 角度翻倍)
    osi_vect(i) = abs(sum(r .* exp(1j * 2 * theta_rad)) / sum(r));
end

% --- 5. 综合绘图展示 ---
for i = 1:nrois
    h_fig = figure('Color', 'w', 'Position', [100, 100, 1000, 450], 'Visible', 'on');
    
    % 子图 1: 原始与平滑调谐曲线
    subplot(1, 2, 1);
    plot(theta, resp_df_f(i,:), 'ko--', 'MarkerFaceColor', 'w', 'DisplayName', 'Raw Delta F/F'); hold on;
    plot(theta, smoothed_resp(i,:), 'r-', 'LineWidth', 2, 'DisplayName', 'Hanning Smoothed');
    xlabel('Direction (deg)'); ylabel('\DeltaF/F (%)');
    title(['ROI ' num2str(i) ' Tuning Curve']);
    xticks(theta); xtickangle(45); grid on;
    legend;

    % 子图 2: 极坐标图与 DI/DSI 标注
    subplot(1, 2, 2);
    p_theta = [theta_rad, theta_rad(1)];
    p_data = [smoothed_resp(i,:), smoothed_resp(i,1)];
    polarplot(p_theta, p_data, 'r-', 'LineWidth', 2);
    
    % 标注算法结果
    res_str = {sprintf('DSI (Vect): %.2f', dsi_vect(i)), ...
               sprintf('OSI (Vect): %.2f', osi_vect(i)), ...
               sprintf('DI (Pref/Opp): %.2f', di_index(i))};
    text(1.2, 1, res_str, 'Units', 'normalized', 'FontSize', 10, 'FontWeight', 'bold');
    title('Spatial Tuning (Smoothed)');
    
    % 自动保存
    saveas(h_fig, fullfile(save_dir, ['ROI_' num2str(i) '_Advanced_Analysis.png']));
end
%% %% 高级调谐曲线分析 - 发放频率版 (Firing Rate Analysis)

%% 高级调谐曲线分析 - 发放频率与显著性对比版

%% --- 1. 参数与路径准备 ---
if ~exist('fs_v', 'var'), fs_v = 400; end % 请根据实际采样率修改
num_oris = stimCfg.numorien;
theta = stimCfg.orientations; 
theta_rad = deg2rad(theta);
cam_frame = table2array(logs.sync(:,4));
stim_framerate = round(1/stimCfg.ifi);
stim_cycles = stimCfg.numorien;
stim_duration_frame = stimCfg.duration * stim_framerate;
stim_rest_frame = stimCfg.isi * stim_framerate;
hanning_kernel = [0.5 1 0.5]; 
stim_range = zeros(stim_cycles,2);
for i = 1:stim_cycles
    stim_range(i,1) = cam_frame((i-1)*(stim_duration_frame + stim_rest_frame) + 1);
    stim_range(i,2) = cam_frame((i-1)*(stim_duration_frame + stim_rest_frame) + stim_duration_frame);  % 得到每次刺激对应的起始/结束帧（每行一组）
end
% 初始化结果矩阵
all_fr_stim = zeros(nrois, num_oris); % 记录每个方向的刺激频率
all_fr_rest = zeros(nrois, num_oris); % 记录每个方向的静息频率
resp_rate   = zeros(nrois, num_oris); % 净增加频率 (Stim - Rest)

% --- 2. 统计每个 Trial 的发放率 ---
fprintf('正在统计峰值并发放频率...\n');
for i = 1:nrois
    roi_peaks = peaks_index{i}; % 当前 ROI 的峰值索引列表
    
    for j = 1:num_oris
        % 获取帧范围
        s_range = stim_range(j,1):stim_range(j,2);
        r_range = rest_range(j,1):rest_range(j,2);
        
        % 计算时间长度 (s)
        dur_s = length(s_range) / fs_v;
        dur_r = length(r_range) / fs_v;
        
        % 统计落在区间内的峰值数
        c_stim = sum(roi_peaks >= s_range(1) & roi_peaks <= s_range(end));
        c_rest = sum(roi_peaks >= r_range(1) & roi_peaks <= r_range(end));
        
        % 计算频率 (Hz)
        all_fr_stim(i, j) = c_stim / dur_s;
        all_fr_rest(i, j) = c_rest / dur_r;
        
        % 净响应 (用于调谐曲线)
        resp_rate(i, j) = max(0, all_fr_stim(i, j) - all_fr_rest(i, j));
    end
end

% --- 3. 批量生成报告 (包含显著性检验) ---
for i = 1:nrois
    % 统计检验：对比该 ROI 在所有方向上的 Rest vs Stim 频率
    % 使用配对 T 检验
    [~, p_val] = ttest(all_fr_rest(i, :), all_fr_stim(i, :));
    
    % 计算 DSI/DI (基于之前平滑后的逻辑)
    raw_r = resp_rate(i, :);
    ext_r = [raw_r(end), raw_r, raw_r(1)];
    sm_r  = conv(ext_r, hanning_kernel, 'valid') / sum(hanning_kernel);
    
    [r_pref, idx] = max(sm_r);
    opp_idx = mod(idx + round(num_oris/2) - 1, num_oris) + 1;
    di_v = (r_pref - sm_r(opp_idx)) / (r_pref + sm_r(opp_idx) + eps);
    
    % --- 绘图 ---
    h_fig = figure('Color', 'w', 'Position', [100, 100, 1200, 500]);
    
    % 子图 1: 频率显著性对比图 (Slope Chart)
    subplot(1, 3, 1);
    plot([1, 2], [all_fr_rest(i,:); all_fr_stim(i,:)], 'Color', [0.8 0.8 0.8], 'Marker', 'o'); 
    hold on;
    % 画出均值线
    errorbar([1, 2], [mean(all_fr_rest(i,:)), mean(all_fr_stim(i,:))], ...
             [std(all_fr_rest(i,:))/sqrt(num_oris), std(all_fr_stim(i,:))/sqrt(num_oris)], ...
             'k-s', 'LineWidth', 2, 'MarkerSize', 10, 'MarkerFaceColor', 'k');
    xlim([0.5 2.5]); xticks([1, 2]); xticklabels({'ISI (Rest)', 'Stimulus'});
    ylabel('Firing Rate (Hz)');
    title(sprintf('Responsiveness (P = %.4f)', p_val));
    grid on;

    % 子图 2: 调谐曲线
    subplot(1, 3, 2);
    bar(theta, resp_rate(i,:), 'FaceColor', [0.9 0.9 0.9], 'EdgeColor', 'none'); hold on;
    plot(theta, sm_r, 'r-o', 'LineWidth', 2);
    xlabel('Direction (deg)'); ylabel('\Delta Firing Rate (Hz)');
    title('Tuning Curve');
    xticks(theta); xtickangle(45);

    % 子图 3: 极坐标
    subplot(1, 3, 3);
    p_theta = [theta_rad, theta_rad(1)];
    p_data  = [sm_r, sm_r(1)];
    polarplot(p_theta, p_data, 'r-o', 'LineWidth', 2);
    title(sprintf('DI: %.2f', di_v));
    
    % 保存
    saveas(h_fig, fullfile(save_dir, ['ROI_' num2str(i) '_Spike_Significance.png']));
end

%% % 每行一个 ROI
linemaxroi = 1;

% 计算需要绘制的行数
plotlines = nrois;

% x 轴范围：兼容 x 和 xV
xlimit = [0, max([max(x(:)), max(xV(:))]) + 1];

% 电信号反转
v_traces_flip = -SNR_traces;

% 电钙共用 y 轴范围
ylimit_all = [ ...
    min([v_traces_flip(:); traces_ds_SNR(:)]), ...
    max([v_traces_flip(:); traces_ds_SNR(:)]) ...
];

% Scale bar
x_bar_time = 2;   % s
x_bar_len = x_bar_time;
y_bar_len_ca = 10;
y_bar_len_v  = 10;

% 创建图形
fig = figure('Color', 'w', 'Position', [100, 100, 1200, 220*nrois]);

% 通用坐标轴属性
axisOpts = {'XTickLabel', [], 'YTickLabel', [], 'TickLength', [0 0], ...
            'XColor', 'none', 'YColor', 'none', 'Box', 'off', 'Color', 'none'};

for i = 0:nrois-1
    ax = subplot(plotlines, 1, i+1);
    hold(ax, 'on');

    % --- 先画电信号 ---
    plot(ax, xV, v_traces_flip(:, i+1), 'r', 'LineWidth', 1);

    % --- 再画钙信号，让钙覆盖在电上面 ---
    plot(ax, x, traces_ds_SNR(:, i+1), 'g', 'LineWidth', 1);

    % 坐标轴设置
    xlim(ax, xlimit);
    ylim(ax, ylimit_all);
    set(ax, axisOpts{:});

    % 只在最后一个 ROI 画 scale bar
    if i == nrois - 1
        x_start = xlimit(2) - x_bar_len * 1.2;
        y_start = ylimit_all(1) + (ylimit_all(2) - ylimit_all(1)) * 0.1;

        % 时间标尺
        plot(ax, [x_start, x_start + x_bar_len], [y_start, y_start], ...
            'k', 'LineWidth', 1.5);

        % 电信号标尺
        plot(ax, [x_start, x_start], ...
            [y_start, y_start + y_bar_len_v], ...
            'r', 'LineWidth', 1.5);

        % 钙信号标尺
        plot(ax, [x_start + x_bar_len, x_start + x_bar_len], ...
            [y_start, y_start + y_bar_len_ca], ...
            'g', 'LineWidth', 1.5);

        % 文字
        text(ax, x_start + x_bar_len/2, y_start, ...
            [num2str(x_bar_time), ' s'], ...
            'VerticalAlignment', 'top', ...
            'HorizontalAlignment', 'center', ...
            'FontSize', 8);

        text(ax, x_start + 0.2, y_start + y_bar_len_v/2, ...
            ['SNR_{v} = ', num2str(y_bar_len_v)], ...
            'Color', 'r', ...
            'HorizontalAlignment', 'left', ...
            'FontSize', 8);

        text(ax, x_start + x_bar_len + 0.2, y_start + y_bar_len_ca/2, ...
            ['SNR_{ca} = ', num2str(y_bar_len_ca)], ...
            'Color', 'g', ...
            'HorizontalAlignment', 'left', ...
            'FontSize', 8);
    end

    % ROI 标题
    title(ax, sprintf('ROI %d', i+1), 'FontSize', 8, 'FontWeight', 'normal');

    hold(ax, 'off');
end

fig_filename = fullfile(save_path, '2_stacked_trace_overlap_ca_on_top.fig');
png_filename = fullfile(save_path, '2_stacked_trace_overlap_ca_on_top.png');

saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');
%% 刺激overlap
% 每行一个 ROI
linemaxroi = 1;

% 计算需要绘制的行数
plotlines = nrois;

% x 轴范围：兼容 x 和 xV
xlimit = [0, max([max(x(:)), max(xV(:))]) + 1];

% 电信号反转
v_traces_flip = -SNR_traces;

% 电钙共用 y 轴范围
ylimit_all = [ ...
    min([v_traces_flip(:); traces_ds_SNR(:)]), ...
    max([v_traces_flip(:); traces_ds_SNR(:)]) ...
];

% Scale bar
x_bar_time = 2;   % s
x_bar_len = x_bar_time;
y_bar_len_ca = 10;
y_bar_len_v  = 10;

% ===== 刺激参数 =====
rest_time = 2;      % 休息 2 s
stim_time = 1;      % 刺激 1 s
cycle_time = rest_time + stim_time;

% 创建图形
fig = figure('Color', 'w', 'Position', [100, 100, 1200, 220*nrois]);

% 通用坐标轴属性
axisOpts = {'XTickLabel', [], 'YTickLabel', [], 'TickLength', [0 0], ...
            'XColor', 'none', 'YColor', 'none', 'Box', 'off', 'Color', 'none'};

for i = 0:nrois-1
    ax = subplot(plotlines, 1, i+1);
    hold(ax, 'on');

    % ===== 先画刺激背景 =====
    % 每个周期: [0,2)休息, [2,3)刺激; [3,5)休息, [5,6)刺激 ...
    stim_starts = rest_time:cycle_time:xlimit(2);
    for k = 1:length(stim_starts)
        stim_start = stim_starts(k);
        stim_end = min(stim_start + stim_time, xlimit(2));

        if stim_start < xlimit(2)
            patch(ax, ...
                [stim_start stim_end stim_end stim_start], ...
                [ylimit_all(1) ylimit_all(1) ylimit_all(2) ylimit_all(2)], ...
                [0 0 1], ...                      % 蓝色
                'FaceAlpha', 0.12, ...           % 透明度
                'EdgeColor', 'none');
        end
    end

    % --- 先画电信号 ---
    plot(ax, xV, v_traces_flip(:, i+1), 'r', 'LineWidth', 1);

    % --- 再画钙信号，让钙覆盖在电上面 ---
    plot(ax, x, traces_ds_SNR(:, i+1), 'g', 'LineWidth', 1);

    % 坐标轴设置
    xlim(ax, xlimit);
    ylim(ax, ylimit_all);
    set(ax, axisOpts{:});

    % 只在最后一个 ROI 画 scale bar
    if i == nrois - 1
        x_start = xlimit(2) - x_bar_len * 1.2;
        y_start = ylimit_all(1) + (ylimit_all(2) - ylimit_all(1)) * 0.1;

        % 时间标尺
        plot(ax, [x_start, x_start + x_bar_len], [y_start, y_start], ...
            'k', 'LineWidth', 1.5);

        % 电信号标尺
        plot(ax, [x_start, x_start], ...
            [y_start, y_start + y_bar_len_v], ...
            'r', 'LineWidth', 1.5);

        % 钙信号标尺
        plot(ax, [x_start + x_bar_len, x_start + x_bar_len], ...
            [y_start, y_start + y_bar_len_ca], ...
            'g', 'LineWidth', 1.5);

        % 文字
        text(ax, x_start + x_bar_len/2, y_start, ...
            [num2str(x_bar_time), ' s'], ...
            'VerticalAlignment', 'top', ...
            'HorizontalAlignment', 'center', ...
            'FontSize', 8);

        text(ax, x_start + 0.2, y_start + y_bar_len_v/2, ...
            ['SNR_{v} = ', num2str(y_bar_len_v)], ...
            'Color', 'r', ...
            'HorizontalAlignment', 'left', ...
            'FontSize', 8);

        text(ax, x_start + x_bar_len + 0.2, y_start + y_bar_len_ca/2, ...
            ['SNR_{ca} = ', num2str(y_bar_len_ca)], ...
            'Color', 'g', ...
            'HorizontalAlignment', 'left', ...
            'FontSize', 8);
    end

    % ROI 标题
    title(ax, sprintf('ROI %d', i+1), 'FontSize', 8, 'FontWeight', 'normal');

    hold(ax, 'off');
end

fig_filename = fullfile(save_path, '2_stacked_trace_overlap_ca_on_top_withStim.fig');
png_filename = fullfile(save_path, '2_stacked_trace_overlap_ca_on_top_withStim.png');

saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');