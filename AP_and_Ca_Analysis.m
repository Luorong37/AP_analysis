% AP_ANALYSIS - Quick Analysis of Voltage Signals
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
% Version: 3.0
% Date: 2024.03
% GitHub: https://github.com/Luorong37/AP_analysis
%
% See also calculate_firing_rate, calculate_FWHM, create_map, calculate_SNR, fit_exp1, highpassfilter, select_ROI

clear;
clc;

%% Loading raw data
clear;clc;
t1 = tic; % Start a timer
nowtime = string(datetime('now'));
% Replace colons with hyphens to get the desired output format
nowtime = strrep(nowtime , ':', '-');
fprintf('Loading...\n') 

% ↓↓↓↓↓-----------Prompt user for define path-----------↓↓↓↓↓
% support for folder, .tif, .tiff, .bin.
folder_path = 'E:\1_Data\LSZ\8.2 HVI2-ST-Cy3b\Exo LplA\ROI1';
file = '';  % must add format.
file_path = fullfile(folder_path, file);
file_path_ca = [file_path,'_Green']; % defined name 
% ↓↓↓↓↓-----------Prompt user for frame rate------------↓↓↓↓↓
freq = 400; % Hz
freq_ca = 10; % Hz
% ↓↓↓↓↓-------------Prompt user for bin-----------------↓↓↓↓↓
bin = 2;
bin_ca = 1;
% -----------------------------------------------------------

% Define parameters
peakfinding = false;
Calcium_analysis = true;
colors = lines(100);
options.colors = colors;

% read path
file_path = fullfile(folder_path, file);
if isfolder(file_path)
    file_name = file;
    file_dir = dir(file_path);
    [~, ~, file_extension] = fileparts(file_dir(3).name);
else
    [~, file_name, file_extension] = fileparts(file_path);
end

% create a folder for analysis
save_path = fullfile(folder_path, [file_name, '_Analysis'], nowtime);
mkdir(save_path);

% Load image file
[movie, ncols, nrows, nframes] = load_movie(file_path);

if Calcium_analysis
    [movie_ca, ncols_ca, nrows_ca, nframes_ca] = load_movie(file_path_ca);
end

% Bin arrange
if bin > bin_ca
    % reshape to high dimension for each bin
    movie_ca = reshape(movie_ca,bin,ncols_ca/bin,bin,nrows_ca/bin,[]);
    % average across x y
    movie_ca = squeeze(mean(mean(movie_ca,1),3));
    % reshape to binned
    movie_ca = reshape(movie_ca,ncols*nrows,nframes_ca);
    nrows_ca = nrows_ca/bin;
    ncols_ca = ncols_ca/bin;
    fprintf('Voltage signal binned\n');
elseif bin < bin_ca
    movie = reshape(movie,bin_ca,ncols/bin_ca,bin_ca,nrows/bin_ca,[]);
    % average across x y
    movie = squeeze(mean(mean(movie,1),3));
    % reshape to binned
    %movie_binned = reshape(movie_ave,ncols,nrows,nframe);
    nrows = nrows/bin_ca;
    ncols = ncols/bin_ca;
    fprintf('Calcium signal binned\n');
end
    
% Define parameters
dt = 1 / freq; % Calculate time axis
dt_ca = 1 / freq_ca;
t = (1:nframes) * dt;
t_ca = (1:nframes_ca) * dt_ca;

% Save code
code_path = fullfile(save_path,'Code');
mkdir(code_path);
currentScript = which("AP_and_Ca_Analysis.m");

% 获取当前脚本依赖的所有文件
[requiredFiles, ~] = matlab.codetools.requiredFilesAndProducts(currentScript);

% 复制当前脚本和所有依赖文件到目标文件夹
for k = 1:length(requiredFiles)
    [~, name, ext] = fileparts(requiredFiles{k});
    copyfile(requiredFiles{k}, fullfile(code_path, [name, ext]));
end

% 提示完成
fprintf('All codes have been copied to %s\n', code_path);

% ----------------------Optional part------------------------

% Save loaded movie (optional)
save_movie = false; % defined as false
if save_movie
    t1 = tic; % Start a timer
    fprintf('Saving...\n')
    raw_filename = fullfile(save_path, '0_Raw_data.mat');
    save(raw_filename,"movie",'ncols','nrows','nframes','freq');
    t2 = toc(t1); % Get the elapsed time
    fprintf('Finished saving movie after %d s\n',round(t2))
end

% Load saved map or mask (optional)
map_path = ''; % define as ''
if isempty(map_path)
    map = [];
    map_ca = [];
else
    map_filename = fullfile(map_path, '0_Sensitivity_Map.mat');
    map = load(map_filename);
    map = map.map;
end

% mask_path = 'E:\0_Code\Luorong\Tools'; % define as ''
mask_path = '';
if isempty(mask_path)
    mask = []; 
    mask_ca = [];
else
    mask_filename = fullfile(mask_path, '1_raw_ROI.mat');
    mask_load = load(mask_filename);
    mask = mask_load.bwmask;
    mask_ca = mask_load.bwmask_ca;
end
t2 = toc(t1); % Get the elapsed time
% -----------------------------------------------------------

%% Create a map 
t1 = tic; % Start a timer
fprintf('Creating...\n')
% if SNR is low, please large the bin.
mapbin = 4; % defined bin = 4
[quick_map] = create_map(movie, nrows, ncols, mapbin);
map = quick_map;
fprintf('Finished voltage map\n')

mapbin_ca = 4; % defined bin_ca = 4
[quick_map] = create_map(movie_ca, nrows_ca, ncols_ca, mapbin_ca,'calcium');
map_ca = quick_map;
fprintf('Finished calcium map\n')

% Visualize correlation coefficients as heatmap
figure()
axe_v = subplot(1,2,1);
imagesc(map);
title('Sensitivity Map');
axis image;
colorbar;

axe_ca = subplot(1,2,2);
imagesc(map_ca);
title('Calcium Map');
axis image;
colorbar;

fig_filename = fullfile(save_path, '0_Sensitivity_Map.fig');
png_filename = fullfile(save_path, '0_Sensitivity_Map.png');
mat_filename = fullfile(save_path, '0_Sensitivity_Map.mat');

saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');
save(mat_filename, 'map', 'map_ca', 'quick_map', 'mapbin', 'mapbin_ca');


t2 = toc(t1); % Get the elapsed time
fprintf('Finished mask creating after %d s\n',round(t2))


%% Photobleaching correction
% 
% mean_movie = mean(movie,1);
% mean_movie_ca = mean(movie_ca,1);
% 
% % fit
% [traces_corrected, fitted_curves] = fit_exp2(mean_movie');
% [traces_corrected_ca, fitted_curves_ca] = fit_exp2(mean_movie_ca');
% 
% % plot
% fig = figure();
% set(fig,'Position',get(0,'Screensize'));
% v_axe = subplot(2,1,1);
% calcium_axe = subplot(2,1,2);
% 
% plot(t,mean_movie,'LineWidth',2,'Parent',v_axe);
% hold(v_axe, 'on');
% plot(t,fitted_curves,'r','LineWidth',2,'Parent',v_axe);
% hold(v_axe, 'on');
% plot(t_ca,mean_movie_ca,'LineWidth',2,'Parent',calcium_axe);
% hold(calcium_axe, 'on');
% plot(t_ca,fitted_curves_ca,'r','LineWidth',2,'Parent',calcium_axe);
% hold(v_axe, 'on');
% 
% % note
% title(v_axe, 'Voltage trace');
% title(calcium_axe, 'Calcium trace');
% 
% fig_filename = fullfile(save_path, '2_fitted_trace.fig');
% png_filename = fullfile(save_path, '2_fitted_trace.png');
% 
% saveas(gcf, fig_filename, 'fig');
% saveas(gcf, png_filename, 'png');
% 
% fprintf('Finished exp2 fit\n')
%% Select ROI
t1 = tic; % Start a timer

% pAce sparaly corrected each ROI
% movie_corrected = movie x`./ fitted_curves';
% movie_ca_corrected = movie_ca ./ fitted_curves_ca';

% with or wihout Mask and Map
correct = false; %true for correct offset
offset_path = '';

if isempty(offset_path)
    offset = [];
else
    offset_filename = fullfile(mask_path, '1_raw_ROI.mat');
    offset_load = load(offset_filename);
    offset = offset_load.offset;
end
[rois, traces, traces_ca, offset] = select_ROI_dual(movie , movie_ca, ...
    nrows, ncols, correct, map, map_ca, mask, mask_ca,offset);

nrois = size(traces,2);
traces_original = traces;
traces_ca_original = traces_ca;


fig_filename = fullfile(save_path, '1_raw_trace.fig');
png_filename = fullfile(save_path, '1_raw_trace.png');
roi_filename = fullfile(save_path, '1_raw_ROI.mat');

saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');
save(roi_filename, 'rois', 'traces', 'traces_ca', 'offset', 'nrows', 'ncols', 'nrois', 'traces_original', 'traces_ca_original');

t2 = toc(t1); % Get the elapsed time
fprintf('Saved ROI figure after %d s\n',round(t2))

%% Background Correction Ca(optional)
rois_ca.bwmask = rois.bwmask_ca;
rois_ca.position = rois.position_ca;
rois_ca.boundary = rois.boundary_ca;
fig = background_correction(movie_ca, ncols_ca, nrows_ca, rois_ca);

waitfor(fig);
% [~,background] = select_ROI(movie, nrows, ncols,background ,[]);
% traces_bgcorr = traces - background;

rois_corr_ca = rois_corrected;
traces_bgfitcorr_ca = traces_bgfitcorr;

figure()
subplot(2,2,1)
imshow(rois_ca.bwmask)
title('Row rois')

subplot(2,2,2)
for i = 1:nrois
    plot(traces_ca(:,i));hold on
end
title('Row traces')

subplot(2,2,3)
imshow(rois_corr_ca.bwmask)
title('Corrected rois')

subplot(2,2,4)
for i = 1:nrois
    plot(traces_bgfitcorr_ca(:,i));hold on
end
title('Background corrected traces')

fig_filename = fullfile(save_path, '1_background_trace_ca.fig');
png_filename = fullfile(save_path, '1_background_trace_ca.png');
save_filename = fullfile(save_path, '1_background_ROI_ca.mat');

saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');
save(save_filename, 'rois', 'rois_corr_ca','background','background_fitted', 'traces_bgcorr', 'traces_bgfitcorr', 'background_mask', 'innerdis', 'outerdis', 'bgthreshold');

traces_ca = traces_bgfitcorr_ca;

%% Background Correction V (optional)
fig = background_correction(movie, ncols, nrows, rois);

waitfor(fig);
% [~,background] = select_ROI(movie, nrows, ncols,background ,[]);
% traces_bgcorr = traces - background;

figure()
subplot(2,2,1)
imshow(rois.bwmask)
title('Row rois')

subplot(2,2,2)
for i = 1:nrois
    plot(traces(:,i));hold on
end
title('Row traces')

subplot(2,2,3)
imshow(rois_corrected.bwmask)
title('Corrected rois')

subplot(2,2,4)
for i = 1:nrois
    plot(traces_bgfitcorr(:,i));hold on
end
title('Background corrected traces')

fig_filename = fullfile(save_path, '1_background_trace.fig');
png_filename = fullfile(save_path, '1_background_trace.png');
roi_filename = fullfile(save_path, '1_background_ROI.mat');

saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');
save(roi_filename, 'rois_corrected','background','background_fitted','background_mask','traces_bgcorr','traces_bgfitcorr', ...
    'innerdis','outerdis','bgthreshold');
traces = traces_bgfitcorr;

%% Correction

% fit
[traces_corrected, fitted_curves] = fit_exp2(traces);

% plot
fig = figure();
set(fig,'Position',get(0,'Screensize'));
fit_axe = subplot(2,1,1);
fited_axe = subplot(2,1,2);

% plot
for i = 1: size(traces_corrected,2)
    plot(traces(:,i),'Color',colors(i,:),'Parent',fit_axe);
    hold(fit_axe, 'on');
    plot(fitted_curves(:,i),'Color',colors(i,:),'LineWidth',2,'Parent',fit_axe);
    hold(fit_axe, 'on');
    plot(traces_corrected(:,i),'Color',colors(i,:),'Parent',fited_axe);
    hold(fited_axe, 'on');
end
hold off;

% note
title(fit_axe, 'Original and Fitted Curves');
title(fited_axe, 'Corrected Traces');
legend(fit_axe, 'Original Trace', 'Fitted Curve');

fig_filename = fullfile(save_path, '2_fitted_trace.fig');
png_filename = fullfile(save_path, '2_fitted_trace.png');

saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');

fprintf('Finished expotential fit\n')

% fit
[traces_corrected_ca, fitted_curves_ca] = fit_exp2(traces_ca);

% plot
fig = figure();
set(fig,'Position',get(0,'Screensize'));
fit_axe = subplot(2,1,1);
fited_axe = subplot(2,1,2);

% plot
for i = 1: size(traces_corrected,2)
    plot(traces_ca(:,i),'Color',colors(i,:),'Parent',fit_axe);
    hold(fit_axe, 'on');
    plot(fitted_curves_ca(:,i),'Color',colors(i,:),'LineWidth',2,'Parent',fit_axe);
    hold(fit_axe, 'on');
    plot(traces_corrected_ca(:,i),'Color',colors(i,:),'Parent',fited_axe);
    hold(fited_axe, 'on');
end
hold off;

% note
title(fit_axe, 'Original and Fitted Curves');
title(fited_axe, 'Corrected Traces');
legend(fit_axe, 'Original Trace', 'Fitted Curve');

fig_filename = fullfile(save_path, '2_fitted_trace_ca.fig');
png_filename = fullfile(save_path, '2_fitted_trace_ca.png');

saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');

fprintf('Finished expotential fit\n')

%% Peak finding (optional)
peakfinding = true; % defined as true
AP_window_width = 40 ; % number of frames to for AP window (defined = 40)
[peaks_polarity, peaks_threshold, peaks_index, peaks_amplitude, peaks_sensitivity] = peak_finding(traces_corrected);


fig_filename = fullfile(save_path, '3_peak_finding.fig');
png_filename = fullfile(save_path, '3_peak_finding.png');
peak_filename = fullfile(save_path, '3_peak_finding.mat');
 
saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');
save(peak_filename, 'peaks_polarity', 'peaks_threshold', 'peaks_index', 'peaks_amplitude', 'peaks_sensitivity');

%% dFF and SNR Analysis
% set up
fig = figure();
set(fig,'Position',get(0,'Screensize'));

SNR_traces = zeros(nframes,nrois);
baselines = zeros(nframes,nrois);

% plot fluorescent image
f_axe = subplot(1,5,1);
f_im = reshape(movie, ncols, nrows, []);hold on;
im_adj = uint16(mean(f_im, 3));
imshow(im_adj,[min(im_adj,[],'all'),max(im_adj,[],'all')]);
for i = 1:nrois
    roi = (rois.bwmask == i);
    boundary = cell2mat(bwboundaries(roi));
    plot(boundary(:, 2), boundary(:, 1), 'Color', colors(i,:), 'LineWidth', 2, 'Parent', f_axe);hold on;
            text(mean(boundary(:, 2)) + 12, mean(boundary(:, 1)) - 12, num2str(i), ...
            'Color', colors(i,:), 'FontSize', 12, 'Parent', f_axe); hold on;
end
hold on;
title('Fluorescent Image');

% plot calcium
calcium_SNR_axe = subplot(1,5,2);
title('calcium');
hold on;

% 设计低通滤波器
fs = 10; % 采样频率 (Hz)
fc = 1; % 截止频率 (Hz)
[b, a] = butter(4, fc/(fs/2)); % 4阶Butterworth低通滤波器
traces_smoothed_ca = zeros(size(traces_corrected_ca));  % 使用双向滤波器进行零相位滤波
SNR_traces_ca = zeros(size(tr aces_corrected_ca));

for i = 1 : nrois
    % Calculate SNR
    traces_smoothed_ca(:,i) = filtfilt(b, a, traces_corrected_ca(:,i));
    SNR_traces_ca(:,i) = calculate_SNR_Ca(traces_corrected_ca(:,i),traces_smoothed_ca(:,i));
end
[~] = offset_plot(SNR_traces_ca,t_ca);
hold off;

% plot SNR
SNR_axe = subplot(1,5,3);
title('SNR');
hold on;
for i = 1 : nrois
    % Calculate SNR
    if peakfinding
        [SNR_traces(:,i),baselines(i)]  = calculate_SNR(traces_corrected(:,i),peaks_index{i},AP_window_width);
    else
        [SNR_traces(:,i),baselines(i)]  = calculate_SNR(traces_corrected(:,i));
    end
% ;
end
[~] = offset_plot(SNR_traces,t);
hold off;

% plot calcium
calcium_dff_axe = subplot(1,5,4);
title('calcium sensitivity');
hold on;

% 设计低通滤波器
fs = 10; % 采样频率 (Hz)
fc = 1; % 截止频率 (Hz)
[b, a] = butter(4, fc/(fs/2)); % 4阶Butterworth低通滤波器
traces_smoothed_ca = zeros(size(traces_corrected_ca));  % 使用双向滤波器进行零相位滤波

for i = 1 : nrois
    % Calculate SNR
    traces_smoothed_ca(:,i) = filtfilt(b, a, traces_corrected_ca(:,i));
end
[~] = offset_plot(traces_smoothed_ca,t_ca);
hold off;

% plot Sensitivity
sensitivity_axe = subplot(1,5,5);
title('Sensitivity');
hold on;
[~] = offset_plot(traces_corrected,t);
hold off;

fig_filename = fullfile(save_path, '4_Signal analysis.fig');
png_filename = fullfile(save_path, '4_Signal analysis.png');
trace_filename = fullfile(save_path, '4_Signal analysis.mat');

save(trace_filename,"SNR_traces_ca",'SNR_traces','traces_smoothed_ca','traces_corrected');
saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');
%% plot stacked figure

% 创建一个新图形窗口
figure('Name','SNR')
linemaxroi = 3; % 每行最多绘制3个ROI
% 计算需要绘制的行数
plotlines = floor(nrois/linemaxroi);
if mod(nrois,linemaxroi) == 0 
    plotlines = plotlines;
else
    plotlines = plotlines + 1;
end

xlimit = [0,max(max(t_ca),max(t))];
ylimitv = [min(SNR_traces,[],'all'),max(SNR_traces,[],'all')];
ylimitca = [min(SNR_traces_ca,[],'all'),max(SNR_traces_ca,[],'all')];

for i = 0: nrois-1
    index = floor(i/linemaxroi)*linemaxroi*2 + mod(i,linemaxroi) + 1;
    tsubplot(plotlines*2,linemaxroi,index,'compact');

    plot(t_ca,SNR_traces_ca(:,i+1),'g');ylim(ylimitca);xlim(xlimit);
    set(gca, 'XTickLabel', []);
    set(gca, 'YTickLabel', []);
    set(gca, 'TickLength', [0 0]);
    set(gca, 'XColor', 'none'); % 隐藏 x 轴线
    set(gca, 'YColor', 'none'); % 隐藏
    set(gca, 'Box', 'off');
    set(gca, 'Color', 'none'); % 设置背景为无色


    tsubplot(plotlines*2,linemaxroi,index+linemaxroi,'compact');
    plot(t,SNR_traces(:,i+1),'r');ylim(ylimitv);xlim(xlimit);
    set(gca, 'YDir', 'reverse');
    set(gca, 'XTickLabel', []);
    set(gca, 'YTickLabel', []);
    set(gca, 'TickLength', [0 0]);
    set(gca, 'XColor', 'none'); % 隐藏 x 轴线
    set(gca, 'YColor', 'none'); % 隐藏
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

fig_filename = fullfile(save_path, '2_stacked_trace.fig');
png_filename = fullfile(save_path, '2_stacked_trace.png');
trace_filename = fullfile(save_path, '2_stacked_trace.mat');

saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');
save(trace_filename, 't','t_ca','SNR_traces','traces_corrected_ca');


% 创建一个新图形窗口
figure('Name','Sensitivity')
linemaxroi = 3; % 每行最多绘制3个ROI
% 计算需要绘制的行数
plotlines = floor(nrois/linemaxroi);
if mod(nrois,linemaxroi) == 0 
    plotlines = plotlines;
else
    plotlines = plotlines + 1;
end

xlimit = [0,max(max(t_ca),max(t))];
ylimitv = [min(traces_corrected,[],'all'),max(traces_corrected,[],'all')];
ylimitca = [min(traces_smoothed_ca,[],'all'),max(traces_smoothed_ca,[],'all')];

for i = 0: nrois-1
    index = floor(i/linemaxroi)*linemaxroi*2 + mod(i,linemaxroi) + 1;
    tsubplot(plotlines*2,linemaxroi,index,'compact');

    plot(t_ca,traces_smoothed_ca(:,i+1),'g');ylim(ylimitca);xlim(xlimit);
    set(gca, 'XTickLabel', []);
    set(gca, 'YTickLabel', []);
    set(gca, 'TickLength', [0 0]);
    set(gca, 'XColor', 'none'); % 隐藏 x 轴线
    set(gca, 'YColor', 'none'); % 隐藏
    set(gca, 'Box', 'off');
    set(gca, 'Color', 'none'); % 设置背景为无色


    tsubplot(plotlines*2,linemaxroi,index+linemaxroi,'compact');
    plot(t,traces_corrected(:,i+1),'r');ylim(ylimitv);xlim(xlimit);
    set(gca, 'YDir', 'reverse');
    set(gca, 'XTickLabel', []);
    set(gca, 'YTickLabel', []);
    set(gca, 'TickLength', [0 0]);
    set(gca, 'XColor', 'none'); % 隐藏 x 轴线
    set(gca, 'YColor', 'none'); % 隐藏
    set(gca, 'Box', 'off');
    set(gca, 'Color', 'none'); % 设置背景为无色
end

fig_filename = fullfile(save_path, '2_stacked_trace_sensitivity.fig');
png_filename = fullfile(save_path, '2_stacked_trace_sensitivity.png');
trace_filename = fullfile(save_path, '2_stacked_trace_sensitivity.mat');

saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');
save(trace_filename, 't','t_ca','traces_corrected','traces_smoothed_ca');


%% Accumulate Voltage signal
% 显示示例轨迹
% sele = 8;
% figure;
% subplot(1, 2, 1);
% trace = -traces(:, sele);
% trace_nol = normalize(trace,'range');
% plot(trace_nol,'Color','r');
% title('Normalized Voltage');
% subplot(1, 2, 2);
% trace_ca = traces_ca(:, sele);
% trace_ca_nol = normalize(trace_ca,'range');
% plot(trace_ca_nol, 'Color', 'g');
% title('Normalized Calcium Signaling');

% 累积电压信号
% set parameter
g = 20;
k = 0.8;
gthr = 0.5;
kth = 0.21;

ftest = @(x, v) x + (g * max(v - gthr, 0) - k * max(x - kth, 0)) * 1/200;

voltage_accumulated = cell(1, nrois);
calcium_normalized = zeros(size(SNR_traces_ca));

for i = 1:nrois
    % Normalize voltage and calcium signals
    volsele = normalize(-SNR_traces(:, i),'range');
    calsele = normalize(SNR_traces_ca(:, i),'range');
    
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

% 保存累积电压信号,和钙对比结果
% for i = 1:length(calculated_set)
%     fig = figure;
%     hold on;
%     plot(t_ca,voltage_accumulated{i}, 'b');
%     plot(t_ca,voltage_accumulated{i} + 1, 'r');
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


%% Statistic AP
AP_list = cell(1, nrois);
each_AP = struct('Trace', [], 'AP_number', [], 'AP_index',[],'AP_amp',[], ...
    'Amplitude', [],'FWHM',[], 'AP_sensitivity',[],'Sensitivity',[], ...
    'AP_SNR', [], 'SNR', []);

% each trace
for i = 1:nrois % i for trace
    peaks_num = length(peaks_index{i});
    peaks_index_i = peaks_index{i};
    peaks_amp_i = peaks_amplitude{i};
    each_trace_amp = traces(:,i);
    each_trace_sensitivity = traces(:,i)-1;
    each_trace_SNR = SNR_traces(:,i);
    AP_list{i} = cell(1, length(peaks_index{i}));

    % each peak
    for j = 1:peaks_num % j for peak
        
        peak_index_ij = peaks_index_i(j);
        peak_amp_ij = peaks_amp_i(j);

        % keep in board
        AP_start_index = max(1, peak_index_ij - AP_window_width);
        AP_end_index = min(nframes, peak_index_ij + AP_window_width);
        AP_index = AP_start_index : AP_end_index;

        % search
        AP_amp = each_trace_amp(AP_start_index:AP_end_index)';
        AP_sensitivity = each_trace_sensitivity(AP_start_index:AP_end_index)';
        AP_SNR = each_trace_SNR(AP_start_index:AP_end_index)';

        % fill NaN
        if 0 > peak_index_ij - AP_window_width
            AP_amp = [NaN(1,0 - (peak_index_ij - AP_window_width)+1), AP_amp];
            AP_sensitivity = [NaN(1,0 - (peak_index_ij - AP_window_width)+1),AP_sensitivity];
            AP_SNR = [NaN(1,0 - (peak_index_ij - AP_window_width)+1),AP_SNR];
        elseif nframes < peak_index_ij + AP_window_width
            AP_amp = [AP_amp, NaN(1,peak_index_ij + AP_window_width - nframes)];
            AP_sensitivity = [AP_sensitivity, NaN(1,peak_index_ij + AP_window_width - nframes)];
            AP_SNR = [AP_SNR, NaN(1,peak_index_ij + AP_window_width - nframes)];
        end

        % Calculate;
        Amplitude = abs(peak_amp_ij);
        Sensitivity = AP_sensitivity(AP_window_width+1) - 1 ;
        SNR = abs(AP_SNR(AP_window_width+1));
        FWHM = calculate_FWHM(AP_amp, dt, peaks_polarity(i));

        % save AP data
        each_AP = struct('Trace', i, 'AP_number', j, 'AP_index',AP_index, ...
            'AP_amp',AP_amp,'Amplitude', Amplitude,'FWHM',FWHM, ...
            'AP_sensitivity',AP_sensitivity,'Sensitivity',Sensitivity, ...
            'AP_SNR', AP_SNR, 'SNR', SNR);
        AP_list{i}{j} = each_AP;
    end
end

% write into excel
% 初始化平均值向量
avg_amp = zeros(length(AP_list), 1);
avg_FWHM = zeros(length(AP_list), 1);
avg_sensitivity = zeros(length(AP_list), 1);
avg_SNR = zeros(length(AP_list), 1);
AP_number = zeros(length(AP_list), 1);
ROI_number = zeros(length(AP_list), 1);

% figure()
% FWHM_axe = subplot(1,3,2);hold on;xlim([0,nrois+1]);
% amp_axe = subplot(1,3,3);hold on;xlim([0,nrois+1]);
% num_axe = subplot(1,3,1);hold on;xlim([0,nrois+1]);
% 计算所有AP的SNR数据并存储在tables中
table_name = fullfile(save_path,'AP_data.xlsx');
for i = 1:length(AP_list)
    
    if cellfun('isempty',AP_list{i}) == 0
        fprintf('Statistic trace %d in %d\n',i, nrois)
        AP_i = AP_list{i}; % 当前trace的所有APs

        % 初始化每个trace的数据向量
        number_i = zeros(length(AP_i), 1);
        amp_i = zeros(length(AP_i), 1);
        FWHM_i = zeros(length(AP_i), 1);
        sensitivity_i = zeros(length(AP_i), 1);
        SNR_i = zeros(length(AP_i), 1);

        for j = 1:length(AP_i)
            
            each_AP = AP_i{j};
            number_i(j) = each_AP.AP_number;
            amp_i(j)  = each_AP.Amplitude;
            FWHM_i(j)  = each_AP.FWHM;
            sensitivity_i(j)  = each_AP.Sensitivity;
            SNR_i(j)  = each_AP.SNR;
        end

        % % plot AP parameters
        % scatter(i,FWHM_i,'filled','Parent',FWHM_axe);hold on;
        % xlabel('ROI number','Parent',FWHM_axe);
        % ylabel('FWHM(ms)','Parent',FWHM_axe);
        % 
        % 
        % scatter(i,amp_i,'filled','Parent',amp_axe);hold on;
        % xlabel('ROI number','Parent',amp_axe);
        % ylabel('Amplitude','Parent',amp_axe);
        % 
        % 
        % scatter(i,number_i,'filled','Parent',num_axe);hold on;
        % xlabel('ROI number','Parent',num_axe);
        % ylabel('AP number','Parent',num_axe);
        

        % 为当前trace创建一个表格
        T = table(number_i, amp_i, FWHM_i, sensitivity_i, SNR_i, ...
            'VariableNames', {'Number', 'Amplitude', 'FWHM (ms)', 'Sensitivity', 'SNR'});

        % 将表格写入Excel的一个新工作表
        sheet_name = string(['ROI ' num2str(i)]);
        writetable(T,table_name, 'Sheet', sheet_name);

        % save average value
        avg_amp(i) = mean(amp_i);
        avg_FWHM(i) = mean(FWHM_i);
        avg_sensitivity(i) = mean(sensitivity_i);
        avg_SNR(i) = mean(SNR_i);
        AP_number(i) = number_i(end);
        ROI_number(i) = i;
    end
end

sgtitle('AP statistic');
fig_filename = fullfile(save_path, '5_AP statistic.fig');
png_filename = fullfile(save_path, '5_AP statistic.png');

saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');

T_ave = table(ROI_number, AP_number, avg_amp, avg_FWHM, avg_sensitivity, avg_SNR, ...
    'VariableNames', {'ROI Number','AP Number', 'Average Amplitude', 'Average FWHM (ms)', 'Average Sensitivity', 'Average SNR'});
writetable(T_ave, table_name, 'Sheet', 'Average');

fprintf('Finished statistic AP\n')

%% save trace for each ROI
traces_path = fullfile(save_path,'eachROI');
mkdir(traces_path);

for i = 1:nrois
    
    % plot aligned each trace of each ROI
    fig = figure();
    subplot(3,1,1);
    current_SNR = SNR_traces(:,i);
    plot(t,current_SNR,'r');hold on;xlim(xlimit);
    ylabel('Original Voltage SNR')
    
    title(sprintf('ROI %d',i),sprintf('Correlation effector = %f',paired_corrtest(i)))
    if peakfinding
        peaks_y = current_SNR(peaks_index{i});
        peaks_x = t(peaks_index{i});
        plot(peaks_x , peaks_y,'v','Color','r','MarkerFaceColor','r','MarkerSize',3);
    end
    set(gca,'YDir','reverse')

    subplot(3,1,2);
    plot(t_ca,voltage_accumulated(:,i),'r');xlim(xlimit);
    ylabel('Integral Voltage')


    subplot(3,1,3);
    current_ca = SNR_traces_ca(:,i);
    plot(t_ca,current_ca,'g');hold on;xlim(xlimit);
    if peakfinding     
        peaks_x = t(peaks_index{i});
        peaks_y = ones(size(peaks_x))*max(current_ca);
        plot(peaks_x , peaks_y,'|','Color','k','MarkerFaceColor','k','MarkerSize',6);
    end

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



%% Heat Map
heatmap_corrtest = zeros(nrois);
heatmap_speartest = zeros(nrois);

for i = 1:nrois
    for j = 1:nrois
        heatmap_corrtest(i, j) = corr(voltage_accumulated(:,i), calcium_normalized(:,j));
        heatmap_speartest(i, j) = corr(voltage_accumulated(:,i), calcium_normalized(:,j), 'Type', 'Spearman');
    end
end

figure;
subplot(1, 2, 1);
imagesc(heatmap_corrtest);
colorbar;
colormap('jet');  % 设定颜色图
axis square;
ylabel('Voltage Index');
xlabel('Calcium Index');
title('Correlation Coefficient Heatmap');

subplot(1, 2, 2);
imagesc(heatmap_speartest);
colorbar;
colormap('jet');  % 设定颜色图
axis square;
ylabel('Voltage Index');
xlabel('Calcium Index');
title('Spearman Correlation Coefficient Heatmap');

fig_filename = fullfile(save_path, '5_Correlation_heatmap.fig');
png_filename = fullfile(save_path, '5_Correlation_heatmap.png');
mat_filename = fullfile(save_path, '5_Correlation_heatmap.mat');

saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');
save(mat_filename,'heatmap_corrtest',"heatmap_speartest");
%% plot fluorescent image
f = figure();
vol_img = subplot(1,2,1);
movie_vol_2D = reshape(mean(movie,2), ncols, nrows, []);
normalized_img = (movie_vol_2D - min(movie_vol_2D(:))) / (max(movie_vol_2D(:)) - min(movie_vol_2D(:)));
imshow(normalized_img);
hold on;

cal_img = subplot(1,2,2);
movie_cal_2D = reshape(mean(movie_ca,2), ncols_ca, nrows_ca, []);
normalized_img_ca = (movie_cal_2D - min(movie_cal_2D(:))) / (max(movie_cal_2D(:)) - min(movie_cal_2D(:)));
imshow(normalized_img_ca);
hold on;

for i = 1:nrois
    roi = (bwmask == i);
    boundary = cell2mat(bwboundaries(roi));
    plot(boundary(:, 2), boundary(:, 1), 'Color', 'r', 'LineWidth', 6, 'Parent',vol_img);
    % 标注ROI编号
    text(mean(boundary(:, 2)) + 6, mean(boundary(:, 1)) - 6, num2str(i), ...
        'Color', 'r', 'FontSize', 36, 'Parent', vol_img); hold on;

    roi_ca = (bwmask_ca == i);
    boundary = cell2mat(bwboundaries(roi_ca));
    plot(boundary(:, 2), boundary(:, 1), 'Color', 'g', 'LineWidth', 6, 'Parent',cal_img);
    % 标注ROI编号
    text(mean(boundary(:, 2)) + 6, mean(boundary(:, 1)) - 6, num2str(i), ...
        'Color', 'g', 'FontSize', 36, 'Parent', cal_img); hold on;
end
fig_filename = fullfile(save_path, '5_fluorescent image.fig');
png_filename = fullfile(save_path, '5_fluorescent image.png');
mat_filename = fullfile(save_path, '5_fluorescent image.mat');

saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');
save(mat_filename,'normalized_img','normalized_img_ca');



%%

%% Save parameter
% 定义保存路径和文件名
save_filename = fullfile(save_path, '-1_workspace_variables.mat');

% 保存当前工作区中的所有变量到.mat文件
clear movie; clear movie_corrected; clear movie_ca; clear movie_corrected_ca;
save(save_filename);


function ax=tsubplot(rows,cols,ind,type)
% @author : slandarer
% gzh  : slandarer随笔

if nargin<4,type='tight';end
sz=[rows,cols];
ratio1=[0,0,1,1];
switch type
    case 'tight'
        ratio1=[0,0,1,1];
        % ratio2=[0.031 0.054 0.9619 0.9254];
    case 'compact'
        ratio1=[0.034 0.0127 0.9256 0.9704];
        % ratio2=[0.065 0.0667 0.8875 0.8958];
    case 'loose'
        ratio1=[0.099 0.056 0.8131 0.8896];
        % ratio2=[0.13 0.11 0.775 0.815];
end
k=1;
posList=zeros(sz(1)*sz(2),4);
for i=1:sz(1)
    for j=1:sz(2)
        tpos=[(j-1)/sz(2),(sz(1)-i)/sz(1),1/sz(2),1/sz(1)];
        posList(k,:)=[tpos(1)+tpos(3).*ratio1(1),tpos(2)+tpos(4).*ratio1(2),...
                      tpos(3).*ratio1(3),tpos(4).*ratio1(4)];
        k=k+1;
    end
end
posSet=posList(ind(:),:);
xmin=min(posSet(:,1));
ymin=min(posSet(:,2));
xmax=max(posSet(:,1)+posSet(:,3));
ymax=max(posSet(:,2)+posSet(:,4));
ax=axes('Parent',gcf,'LooseInset',[0,0,0,0],...
    'OuterPosition',[xmin,ymin,xmax-xmin,ymax-ymin]);
% @author : slandarer
% gzh  : slandarer随笔
end
