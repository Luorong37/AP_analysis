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
folder_path = 'D:\2024.05.30_dual_P2A';
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
[~, file_name, file_extension] = fileparts(file);
save_path = fullfile(folder_path, [file_name, '_Analysis'], nowtime);
mkdir(save_path);

% when read a folder
if isempty(file_extension)
    file_extension = 'tif';
end

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
currentScript = which("AP_Analysis.m");

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
    index = floor(i/plotlines)*plotlines*2 + mod(i,plotlines) + 1;
    subplot(linemaxroi*2,plotlines,index)
    plot(t,traces(:,i+1),'r');xlim(xlimit);ylim(ylimitv);
    subplot(linemaxroi*2,plotlines,index+plotlines)
    plot(t_ca,traces_ca(:,i+1),'g');xlim(xlimit);ylim(ylimitca);
    ylabel(sprintf('ROI %d', i + 1))

end

fig_filename = fullfile(save_path, '2_stacked_trace.fig');
png_filename = fullfile(save_path, '2_stacked_trace.png');
trace_filename = fullfile(save_path, '2_stacked_trace.mat');

saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');
save(trace_filename, 't','t_ca','traces','traces_ca');

%% Save parameter
% 定义保存路径和文件名
save_filename = fullfile(save_path, '-1_workspace_variables.mat');

% 保存当前工作区中的所有变量到.mat文件
clear movie; clear movie_corrected; clear movie_ca; clear movie_corrected_ca;
save(save_filename);
