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
folder_path = 'E:\1_Data\TEST';
file = '20240322-161555right!';  % must add format.
% ↓↓↓↓↓-----------Prompt user for frame rate------------↓↓↓↓↓
freq = 400; % Hz
% -----------------------------------------------------------

% read path
file_path = fullfile(folder_path, file);
[~, file_name, file_extension] = fileparts(file);
save_path = fullfile(folder_path, [file_name, '_Analysis'], nowtime);
mkdir(save_path);

% when read a folder
if isempty(file_extension)
    file_extension = 'tif';
end

% Load image file
[movie, ncols, nrows, nframes] = load_movie(file_path,file_extension,100000);

% Define parameters
dt = 1 / freq; % Calculate time axis
colors = [lines(7);hsv(5);spring(3);winter(3);gray(3)];
t = (1:nframes) * dt;
peakfinding = false;

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

t2 = toc(t1); % Get the elapsed time
fprintf('Finished loading after %d s\n',round(t2))

% Load saved map or mask (optional)
map_path = ''; % define as ''
if isempty(map_path)
    map = [];
else
    map_filename = fullfile(map_path, '0_Sensitivity_Map.mat');
    map = load(map_filename);
    map = map.map;
end

mask_path = ''; % define as ''
if isempty(mask_path)
    mask = []; 
else
    mask_path = 'E:\Data\20230810\20230810-161506recordSCN_Analysis\2023-11-06 09-12-40';
    mask_filename = fullfile(mask_path, '1_raw_ROI.mat');
    mask = load(mask_filename);
    mask = mask.rois;
end
% -----------------------------------------------------------

%% Create a map (optional)
t1 = tic; % Start a timer
fprintf('Creating...\n')
% if SNR is low, please large the bin.
bin = 4; % defined bin = 4
[quick_map] = create_map(movie, nrows, ncols, bin);
map = quick_map;

% Visualize correlation coefficients as heatmap
figure()
imagesc(quick_map);
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

%% Select ROI
t1 = tic; % Start a timer

% with or wihout Mask and Map
[rois, traces] = select_ROI(movie, nrows, ncols, t, colors, mask, map);

fig_filename = fullfile(save_path, '1_raw_trace.fig');
png_filename = fullfile(save_path, '1_raw_trace.png');
roi_filename = fullfile(save_path, '1_raw_ROI.mat');

saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');
save(roi_filename, 'rois');

traces_input = zeros(size(traces,1),size(traces,2)-1);
background = traces(:,end);

% remove background
for i = 1: size(traces,2)-1
    traces_input(:,i) = traces(:,i) - background;
end

t2 = toc(t1); % Get the elapsed time
fprintf('Saved ROI figure after %d s\n',round(t2))

%% Photobleaching correction

% fit
[traces_corrected, fitted_curves] = fit_exp1(traces_input);

% plot
fig = figure();
set(fig,'Position',get(0,'Screensize'));
fit_axe = subplot(2,1,1);
fited_axe = subplot(2,1,2);

for i = 1: size(traces_corrected,2)
    plot(t,traces_input(:,i),'Color',colors(i,:),'Parent',fit_axe);
    hold(fit_axe, 'on');
    plot(t,fitted_curves(:,i),'Color',colors(i,:),'LineWidth',2,'Parent',fit_axe);
    hold(fit_axe, 'on');
    plot(t,traces_corrected(:,i),'Color',colors(i,:),'Parent',fited_axe);
    hold(fited_axe, 'on');
end

% note
title(fit_axe, 'Original and Fitted Curves');
title(fited_axe, 'Corrected Traces');
legend(fit_axe, 'Original Trace', 'Fitted Curve');

fig_filename = fullfile(save_path, '2_fitted_trace.fig');
png_filename = fullfile(save_path, '2_fitted_trace.png');

saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');

fprintf('Finished exp1 fit\n')

%% Peak finding (optional)
peakfinding = true; % defined as true

if peakfinding

AP_window_width = 40 ; % number of frames to for AP window (defined = 40)
[peak_polarity, peak_threshold, peaks_index, peaks_amplitude, peaks_sensitivity] = peak_finding(traces_corrected, t, colors, rois);

% plot trace
fig = figure();
set(fig,'Position',get(0,'Screensize'));
offset_plot(traces_corrected,t,colors)
% plot label
for i = 1:length(rois)-1
    % plot threshold
    plot(t,ones(1,length(t)).*(1-peak_threshold(i)*peak_polarity(i)) + offset_peak,'Color',colors(i,:),'LineWidth',2); hold on;
    % plot peak
    dt = t(2)-t(1);
    plot(peaks_index{i}.*dt, (1-peaks_sensitivity{i}*peak_polarity(i))+ offset_peak,'v', ...
        'Color',colors(i,:),'MarkerFaceColor',colors(i,:)); hold on;

end
end

fig_filename = fullfile(save_path, '3_peak_finding.fig');
png_filename = fullfile(save_path, '3_peak_finding.png');

saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');

%% Sensitivity and SNR Analysis
% set up
fig = figure();
set(fig,'Position',get(0,'Screensize'));
intensity_trace = zeros(nframes,length(rois));
sentivity_trace = zeros(nframes,length(rois));
SNR_trace = zeros(nframes,length(rois));
fitted_trace = zeros(nframes,length(rois));
options.colors = colors;

 % plot sensitivity
sensitivity_axe = subplot(1,2,1);
title('Sensitivity');
hold on;
for i = 1 : length(rois)-1
    % Calculate Sensitivity
    sentivity_trace(:,i) = traces_corrected(:,i)-1;
    [~] = offset_plot(sentivity_trace,t,options);
end

% plot SNR
SNR_axe = subplot(1,2,2);
title('SNR');
hold on;

for i = 1 : length(rois)-1
    % Calculate SNR
    if peakfinding
        [SNR_trace(:,i),fitted_trace(:,i)]  = calculate_SNR(traces_corrected(:,i), peak_threshold(i), peak_polarity(i));
    else
        [SNR_trace(:,i),fitted_trace(:,i)]  = calculate_SNR(traces_corrected(:,i));
    end
        
    [~] = offset_plot(SNR_trace,t,options);
end

fig_filename = fullfile(save_path, '4_Sensitivity_figure.fig');
png_filename = fullfile(save_path, '4_Sensitivity_figure.png');
trace_filename = fullfile(save_path, '4_Sensitivity_data.mat');

save(trace_filename,"sentivity_trace",'SNR_trace');
saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');
%% Statistic AP
AP_list = cell(1, length(rois)-1);
each_AP = struct('Trace', [], 'AP_number', [], 'AP_index',[],'AP_amp',[], ...
    'Amplitude', [],'FWHM',[], 'AP_sensitivity',[],'Sensitivity',[], ...
    'AP_SNR', [], 'SNR', []);

% each trace
for i = 1:length(rois)-1 % i for trace
    peaks_num = length(peaks_index{i});
    peaks_index_i = peaks_index{i};
    peaks_amp_i = peaks_amplitude{i};
    each_trace_amp = traces_input(:,i)*peak_polarity(i);
    each_trace_sensitivity = sentivity_trace(:,i);
    each_trace_SNR = SNR_trace(:,i);
    each_trace_smooth = fitted_trace(:,i);
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
        AP_smooth = each_trace_smooth(AP_start_index:AP_end_index)';

        % fill NaN
        if 0 > peak_index_ij - AP_window_width
            AP_amp = [NaN(1,0 - (peak_index_ij - AP_window_width)+1), AP_amp];
            AP_sensitivity = [NaN(1,0 - (peak_index_ij - AP_window_width)+1),AP_sensitivity];
            AP_SNR = [NaN(1,0 - (peak_index_ij - AP_window_width)+1),AP_SNR];
            AP_smooth = [NaN(1,0 - (peak_index_ij - AP_window_width)+1),AP_smooth];
        elseif nframes < peak_index_ij + AP_window_width
            AP_amp = [AP_amp, NaN(1,peak_index_ij + AP_window_width - nframes)];
            AP_sensitivity = [AP_sensitivity, NaN(1,peak_index_ij + AP_window_width - nframes)];
            AP_SNR = [AP_SNR, NaN(1,peak_index_ij + AP_window_width - nframes)];
            AP_smooth = [AP_smooth, NaN(1,peak_index_ij + AP_window_width - nframes)];
        end

        % Calculate;
        Amplitude = abs(peak_amp_ij);
        Sensitivity = AP_sensitivity(AP_window_width+1) - 1 ;
        SNR = abs(AP_SNR(AP_window_width+1));
        baseline = mean(AP_smooth,'omitnan');
        FWHM = calculate_FWHM(AP_sensitivity, dt, AP_window_width, baseline, peak_polarity(i));

        % save AP data
        each_AP = struct('Trace', i, 'AP_number', j, 'AP_index',AP_index, ...
            'AP_amp',AP_amp,'Amplitude', Amplitude,'FWHM',FWHM, ...
            'AP_sensitivity',AP_sensitivity,'Sensitivity',Sensitivity, ...
            'AP_SNR', AP_SNR, 'SNR', SNR,'baseline',baseline);
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

% 计算所有AP的SNR数据并存储在tables中
table_name = fullfile(save_path,'AP_data.xlsx');
for i = 1:length(AP_list)
    if cellfun('isempty',AP_list{i}) == 0
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

T_ave = table(ROI_number, AP_number, avg_amp, avg_FWHM, avg_sensitivity, avg_SNR, ...
    'VariableNames', {'ROI Number','AP Number', 'Average Amplitude', 'Average FWHM (ms)', 'Average Sensitivity', 'Average SNR'});
writetable(T_ave, table_name, 'Sheet', 'Average');

fprintf('Finished statistic AP\n')

%% Plot average AP sensitivity
figure();
set(gcf,'Position',get(0,'Screensize'));
% 统计不为空的trace数目
plot_cols = sum(cellfun('isempty',AP_list)==0)+1;
plot_col = 0;

for i = 1:length(rois)-1 % i for trace
    %判断是否为有AP的trace
    peaks_num = length(peaks_index{i});
    if cellfun('isempty',AP_list{i}) == 0
        plot_col = plot_col + 1;
        subplot(2,ceil(plot_cols/2),plot_col);

        % set y axis direction
        if peak_polarity(i) == -1
            set(gca,'YDir','reverse')
            hold on;
        end

        % get each AP
        AP_i = zeros(peaks_num, AP_window_width*2+1);
        for j = 1:peaks_num
            each_AP = AP_list{i}{j};
            AP_i(j,:) = each_AP.AP_sensitivity;
            % plot each AP
            plot((1:AP_window_width*2+1)*dt, each_AP.AP_sensitivity','Color',[0.8 0.8 0.8]);
            hold on;
        end

        % plot average AP for each trace
        subplot(2,ceil(plot_cols/2),plot_col);
        plot((1:AP_window_width*2+1)*dt, mean(AP_i,1,'omitnan'),'Color',colors(i,:),'LineWidth',1);
        hold on;
        title(sprintf('ROI %d\n',i));
        hold on;

        % plot average AP for average trace
        subplot(2,ceil(plot_cols/2),plot_cols);
        plot((1:AP_window_width*2+1)*dt, mean(AP_i,1,'omitnan'),'Color',colors(i,:),'LineWidth',1);
        hold on;
    end
end
title('Averaged of All');
hold on;
sgtitle('Averaged Sensitivity');
hold on;
fig_filename = fullfile(save_path, '5_average_AP_sensitivity.fig');
png_filename = fullfile(save_path, '5_average_AP_sensitivity.png');

saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');
%% Plot average AP SNR
% plot
figure();
title('Statistic AP');
hold on;

% 统计不为空的trace数目
plot_cols = sum(cellfun('isempty',AP_list)==0)+1;
plot_col = 0;
subplot(1,plot_cols,plot_cols);


for i = 1:length(rois)-1
    %判断是否为有AP的trace
    peaks_num = length(peaks_index{i});
    if cellfun('isempty',AP_list{i}) == 0
        plot_col = plot_col + 1;
        subplot(2,ceil(plot_cols/2),plot_col);

        % set y axis direction
        if peak_polarity(i) == -1
            set(gca,'YDir','reverse')
            hold on;
        end

        % extract each AP
        AP_i = zeros(peaks_num, AP_window_width*2+1);
        for j = 1:peaks_num
            each_AP = AP_list{i}{j};
            AP_i(j,:) = each_AP.AP_SNR;
            plot((1:AP_window_width*2+1)*dt, each_AP.AP_SNR','Color',[0.8 0.8 0.8]);
            hold on;
        end

        % plot average AP for trace
        subplot(2,ceil(plot_cols/2),plot_col);
        plot((1:AP_window_width*2+1)*dt, mean(AP_i,1,'omitnan'),'Color',colors(i,:),'LineWidth',1);
        hold on;
        title(sprintf('ROI %d',i));
        hold on;

        % plot average AP for all
        subplot(2,ceil(plot_cols/2),plot_cols);
        plot((1:AP_window_width*2+1)*dt, mean(AP_i,1,'omitnan'),'Color',colors(i,:),'LineWidth',1);
        hold on;
    end
end

sgtitle('Average SNR');
hold on;

fig_filename = fullfile(save_path, '6_average_AP_SNR.fig');
png_filename = fullfile(save_path, '6_average_AP_SNR.png');

saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');
%
%% Firing rate

% set up
window = 1; % second, (def = 1 s)
num_windows = ceil(nframes / window*dt)-1;
window_length = ceil(window / dt);
firing_rate_traces = zeros(num_windows, length(rois));
t_window = (0:window:num_windows-1)+0.5 * window;

% 统计不为空的trace数目
figure();
shift_firingrate = 0;
for i = 1:length(rois)-1
    %判断是否为有AP的trace
    peaks_num = length(peaks_index{i});
    firing_rate_traces(:,i) = calculate_firing_rate(peaks_index{i}, window,num_windows, window_length);
    plot(t_window, firing_rate_traces(:,i) + shift_firingrate,'Color',colors(i,:),'LineWidth',2);
    hold on;
    plot(t_window, shift_firingrate*ones(size(t_window)),'black');
    hold on;
    shift_frstep = max(firing_rate_traces(:,i));
    shift_firingrate = shift_firingrate+ shift_frstep;
end

sgtitle('Firing rate');
hold on;

fig_filename = fullfile(save_path, '7_firing_rate.fig');
png_filename = fullfile(save_path, '7_firing_rate.png');

saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');

%% Save parameter
% 定义保存路径和文件名
save_filename = fullfile(save_path, '-1_workspace_variables.mat');

% 保存当前工作区中的所有变量到.mat文件
clear movie;
save(save_filename);
