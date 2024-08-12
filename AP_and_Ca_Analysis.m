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
t1 = tic; % Start a timer
nowtime = string(datetime('now'));
% Replace colons with hyphens to get the desired output format
nowtime = strrep(nowtime , ':', '-');
fprintf('Loading...\n') 

% ↓↓↓↓↓-----------Prompt user for define path-----------↓↓↓↓↓
% support for folder, .tif, .tiff, .bin.
folder_path = 'G:\2024.05.30_dual_P2A';
file = '20240530-182118POA';  % must add format.
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
colors = [lines(7);hsv(5);spring(3);winter(3);gray(3)];
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
[movie, ncols, nrows, nframes] = load_movie(file_path,file_extension,100000);

if Calcium_analysis
    [movie_ca, ncols_ca, nrows_ca, nframes_ca] = load_movie(file_path_ca,file_extension,100);
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

mask_path = ''; % define as ''
if isempty(mask_path)
    mask = []; 
    mask_ca = [];
else
    mask_path = 'E:\1_Data\Luorong\430_5min\20240430-172447_5min_2_Analysis\2024-05-07 14-56-01';
    mask_filename = fullfile(mask_path, '1_raw_ROI.mat');
    mask = load(mask_filename);
    mask = mask.rois;
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
save(mat_filename, 'map','map_ca');

t2 = toc(t1); % Get the elapsed time
fprintf('Finished mask creating after %d s\n',round(t2))
%% Photobleaching correction

mean_movie = mean(movie,1);
mean_movie_ca = mean(movie_ca,1);

% fit
[traces_corrected, fitted_curves] = fit_exp2(mean_movie');
[traces_corrected_ca, fitted_curves_ca] = fit_exp2(mean_movie_ca');

% plot
fig = figure();
set(fig,'Position',get(0,'Screensize'));
v_axe = subplot(2,1,1);
ca_axe = subplot(2,1,2);

plot(t,mean_movie,'LineWidth',2,'Parent',v_axe);
hold(v_axe, 'on');
plot(t,fitted_curves,'r','LineWidth',2,'Parent',v_axe);
hold(v_axe, 'on');
plot(t_ca,mean_movie_ca,'LineWidth',2,'Parent',ca_axe);
hold(ca_axe, 'on');
plot(t_ca,fitted_curves_ca,'r','LineWidth',2,'Parent',ca_axe);
hold(v_axe, 'on');

% note
title(v_axe, 'Voltage trace');
title(ca_axe, 'Calcium trace');

fig_filename = fullfile(save_path, '2_fitted_trace.fig');
png_filename = fullfile(save_path, '2_fitted_trace.png');

saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');

fprintf('Finished exp2 fit\n')
%% Select ROI
t1 = tic; % Start a timer

% movie = movie ./ fitted_curves';
% movie_ca = movie_ca ./ fitted_curves_ca';

% with or wihout Mask and Map
correct = true; % true for correct offset
[bwmask, bwmask_ca, traces, traces_ca] = select_ROI_dual(movie, movie_ca, ...
    nrows, ncols,nrows_ca, ncols_ca, t, t_ca, colors, mask, map, mask_ca, map_ca, correct);

nrois = size(traces,2);

fig_filename = fullfile(save_path, '1_raw_trace.fig');
png_filename = fullfile(save_path, '1_raw_trace.png');
roi_filename = fullfile(save_path, '1_raw_ROI.mat');

saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');
save(roi_filename, 'bwmask', 'bwmask_ca');

t2 = toc(t1); % Get the elapsed time
fprintf('Saved ROI figure after %d s\n',round(t2))

%% plot stacked figure
% figure()
% subplot(1,2,1)
% stackedplot(t,traces); 
% subplot(1,2,2)
% stackedplot(t_ca,traces_ca);

% 准备列名
% numTraces = size(traces_ca, 2);  % 获取数据列数
% columnNames = cell(1, numTraces);
% for i = 1:numTraces
%     columnNames{i} = sprintf('trace%d', i);
% end
% trace_Time = [t', traces];
% trace_ca_Time = [t_ca', traces_ca];
% 
% trace_table = array2table(trace_Time,'VariableNames', ['Time', columnNames]);
% trace_ca_table = array2table(trace_ca_Time,'VariableNames', ['Time', columnNames]);
% 
% stackedplot(trace_table,trace_ca_table, 'XVariable', 'Time');
% 创建一个新图形窗口
figure()
linemaxroi = 5;
plotlines = floor(nrois/linemaxroi);
if mod(nrois,5) == 0 
    plotlines = plotlines;
else
    plotlines = plotlines + 1;
end
xlimit = [0,max(max(t_ca),max(t))];
ylimitv = [min(traces,[],'all'),max(traces,[],'all')];
ylimitca = [min(traces_ca,[],'all'),max(traces_ca,[],'all')];
for i = 0: nrois-1
    % lines passed + corss index
    index = floor(i/linemaxroi)*linemaxroi*2 + mod(i,linemaxroi) + 1;
    subplot(plotlines*2,linemaxroi,index)
    plot(t,traces(:,i+1),'r');xlim(xlimit);ylim(ylimitv);
    subplot(plotlines*2,linemaxroi,index+linemaxroi)
    plot(t_ca,traces_ca(:,i+1),'g');xlim(xlimit);ylim(ylimitca);
    ylabel(sprintf('ROI %d', i + 1))
end

fig_filename = fullfile(save_path, '2_stacked_trace.fig');
png_filename = fullfile(save_path, '2_stacked_trace.png');
trace_filename = fullfile(save_path, '2_stacked_trace.mat');

saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');
save(trace_filename, 't','t_ca','traces','traces_ca');


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
calcium_normalized = zeros(size(traces_ca));

for i = 1:nrois
    % Normalize voltage and calcium signals
    volsele = normalize(-traces(:, i),'range');
    calsele = normalize(traces_ca(:, i),'range');
    
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
    subplot(plotlines*2,linemaxroi,index+linemaxroi)
    plot(calcium_normalized(:,i+1),'g');ylim(ylimitca);
    ylabel(sprintf('ROI %d', i+1))
end


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

%% save trace for each ROI
traces_path = fullfile(save_path,'eachROI');
mkdir(traces_path);

for i = 1:nrois
    % plot aligned each trace of each ROI
    fig = figure();
    subplot(3,1,1);
    plot(t,traces(:,i),'r');xlim(xlimit);
    ylabel('Original Voltage')
    title(sprintf('ROI %d',i),sprintf('Correlation effector = %f',paired_corrtest(i)))
    subplot(3,1,2);
    plot(t_ca,voltage_accumulated(:,i),'r');xlim(xlimit);
    ylabel('Integral Voltage')
    subplot(3,1,3);
    plot(t_ca,calcium_normalized(:,i),'g');xlim(xlimit);
    ylabel('Calcium')
    xlabel('Time')
    
    % save the figure
    fig_filename = fullfile(traces_path, sprintf('ROI %d.fig', i));
    png_filename = fullfile(traces_path, sprintf('ROI %d.png', i));
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


%% Save parameter
% 定义保存路径和文件名
save_filename = fullfile(save_path, '-1_workspace_variables.mat');

% 保存当前工作区中的所有变量到.mat文件
clear movie; clear movie_corrected; clear movie_ca; clear movie_corrected_ca;
save(save_filename);
