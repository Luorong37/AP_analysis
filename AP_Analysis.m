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

clear; clc;

%% Loading raw data
clear; clc;

t1 = tic; % Start a timer
nowtime = string(datetime('now'));
% Replace colons with hyphens to get the desired output format
nowtime = strrep(nowtime , ':', '-');
fprintf('Loading...\n')

% ↓↓↓↓↓-----------Prompt user for define path-----------↓↓↓↓↓
% support for folder, .tif, .tiff, .bin.
folder_path = 'E:\1_Data\Luorong\20240709_optopatch\\';
file = 'sti90%\';  % must add format.
% ↓↓↓↓↓-----------Prompt user for frame rate------------↓↓↓↓↓
freq = 400; % Hz
% -----------------------------------------------------------

% read path
file_path = fullfile(folder_path, file);
if isfolder(file_path)
    file_name = file;
    file_dir = dir(file_path);
    [~, ~, file_extension] = fileparts(file_dir.name(1));
else
    [~, file_name, file_extension] = fileparts(file_path);
end
% % when read a folder
% if isfolder(file_path)
%     file_extension = '.tif';
% end

% create a folder for analysis
save_path = fullfile(folder_path, [file_name, '_Analysis'], nowtime);
mkdir(save_path);

% Load image file

[movie, ncols, nrows, nframes] = load_movie(file_path);


% Define parameters
dt = 1 / freq; % Calculate time axis
% colors = [lines(7);hsv(5);spring(3);winter(3);gray(3)];
colors = lines(100);
t = (1:nframes) * dt;
peakfinding = true;
options.colors = colors;

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
t2 = toc(t1); % Get the elapsed time
fprintf('Finished loading after %d s\n',round(t2))


% ----------------------Optional part------------------------

% Save data
save_stack = false;
if save_stack
    create_tiff_stack(file_path);
end

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
    mask_filename = fullfile(mask_path, '1_raw_ROI.mat');
    mask = load(mask_filename);
    mask = mask.bwmask;
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
[bwmask, traces] = select_ROI(movie, nrows, ncols, t, mask, map);
nrois = max(bwmask,[],'all');

% if do not need a map, run the following code:
% map = [];
% [bwmask, traces] = select_ROI(movie, nrows, ncols, t, mask, map);

fig_filename = fullfile(save_path, '1_raw_trace.fig');
png_filename = fullfile(save_path, '1_raw_trace.png');
roi_filename = fullfile(save_path, '1_raw_ROI.mat');

saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');
save(roi_filename, 'bwmask');

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
    plot(t,traces(:,i),'Color',colors(i,:),'Parent',fit_axe);
    hold(fit_axe, 'on');
    plot(t,fitted_curves(:,i),'Color',colors(i,:),'LineWidth',2,'Parent',fit_axe);
    hold(fit_axe, 'on');
    plot(t,traces_corrected(:,i),'Color',colors(i,:),'Parent',fited_axe);
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

%% Peak finding
peakfinding = true; % defined as true
AP_window_width = 40 ; % number of frames to for AP window (defined = 40)
    
if peakfinding
    [peaks_polarity, peaks_threshold, peaks_index, peaks_amplitude, peaks_sensitivity] = ...
                                            peak_finding(traces_corrected);
    traces_corr_flipped = traces_corrected.*peaks_polarity + 1 - peaks_polarity;

    % plot trace
    fig = figure();
    set(fig,'Position',get(0,'Screensize'));
    offset_array = offset_plot(traces_corr_flipped,t); 
    xlim tight
    sgtitle('Peak finding');
    % plot label
    for i = 1:nrois
        % plot threshold
        plot(t,ones(1,numel(t)).*(peaks_threshold(i)*peaks_polarity(i)+ 1 -peaks_polarity(i)) ...
            + offset_array(i),'Color',colors(i,:),'LineWidth',2); hold on;
        % plot peak
        dt = t(2)-t(1);
        plot(peaks_index{i}.*dt, ( peaks_amplitude{i}.*peaks_polarity(i)+ 1 -peaks_polarity(i)) ...
            + offset_array(i),'v','Color',colors(i,:),'MarkerFaceColor',colors(i,:)); hold on;
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
sentivity_trace = zeros(nframes,nrois);
SNR_traces = zeros(nframes,nrois);
baselines = zeros(nframes,nrois);

% plot fluorescent image
f_axe = subplot(1,3,1);
f_im = reshape(movie, ncols, nrows, []);hold on;
im_adj = uint16(mean(f_im, 3));
imshow(im_adj,[min(im_adj,[],'all'),max(im_adj,[],'all')]);
for i = 1:nrois
    roi = (bwmask == i);
    boundary = cell2mat(bwboundaries(roi));
    plot(boundary(:, 2), boundary(:, 1), 'Color', colors(i,:), 'LineWidth', 2, 'Parent', f_axe);hold on;
            text(mean(boundary(:, 2)) + 12, mean(boundary(:, 1)) - 12, num2str(i), ...
            'Color', colors(i,:), 'FontSize', 12, 'Parent', f_axe); hold on;
end
hold on;
title('Fluorescent Image');

% plot sensitivity
sensitivity_axe = subplot(1,3,2);
title('Sensitivity');
hold on;
for i = 1 : nrois
    % Calculate Sensitivity
    % sentivity_trace(:,i) = traces_corrected(:,i)*peaks_polarity(i)+ 1 -peaks_polarity(i);
    sentivity_trace(:,i) = traces_corrected(:,i) - 1 ;
end
    [~] = offset_plot(sentivity_trace,t);

% plot SNR
SNR_axe = subplot(1,3,3);
title('SNR');
hold on;
for i = 1 : nrois
    % Calculate SNR
    if peakfinding
        [SNR_traces(:,i),baselines(i)]  = calculate_SNR(traces_corrected(:,i),peaks_index{i},AP_window_width);
    else
        [SNR_traces(:,i),baselines(i)]  = calculate_SNR(traces_corrected(:,i));
    end
    [~] = offset_plot(SNR_traces,t);
end

fig_filename = fullfile(save_path, '4_SNR.fig');
png_filename = fullfile(save_path, '4_SNR.png');
trace_filename = fullfile(save_path, '4_SNR.mat');

save(trace_filename,"sentivity_trace",'SNR_traces');
saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');
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
    each_trace_sensitivity = sentivity_trace(:,i);
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

figure()
FWHM_axe = subplot(1,3,2);hold on;xlim([0,nrois+1]);
amp_axe = subplot(1,3,3);hold on;xlim([0,nrois+1]);
num_axe = subplot(1,3,1);hold on;xlim([0,nrois+1]);
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

        % plot AP parameters
        scatter(i,FWHM_i,'filled','Parent',FWHM_axe);hold on;
        xlabel('ROI number','Parent',FWHM_axe);
        ylabel('FWHM(ms)','Parent',FWHM_axe);
       
        
        scatter(i,amp_i,'filled','Parent',amp_axe);hold on;
        xlabel('ROI number','Parent',amp_axe);
        ylabel('Amplitude','Parent',amp_axe);
        

        scatter(i,number_i,'filled','Parent',num_axe);hold on;
        xlabel('ROI number','Parent',num_axe);
        ylabel('AP number','Parent',num_axe);
        

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

%% Plot average AP sensitivity
figure();
% 统计不为空的trace数目
plot_cols = sum(cellfun('isempty',AP_list)==0)+1;
plot_col = 0;

for i = 1:nrois % i for trace
    %判断是否为有AP的trace
    peaks_num = length(peaks_index{i});
    if cellfun('isempty',AP_list{i}) == 0
        plot_col = plot_col + 1;
        subplot(2,ceil(plot_cols/2),plot_col);

        % set y axis direction
        if peaks_polarity(i) == -1
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
fig_filename = fullfile(save_path, '6_average_AP_sensitivity.fig');
png_filename = fullfile(save_path, '6_average_AP_sensitivity.png');

saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');


% Plot average AP SNR
% plot
figure();
title('Statistic AP');
hold on;

% 统计不为空的trace数目
plot_cols = sum(cellfun('isempty',AP_list)==0)+1;
plot_col = 0;
subplot(1,plot_cols,plot_cols);


for i = 1:nrois
    %判断是否为有AP的trace
    peaks_num = length(peaks_index{i});
    if cellfun('isempty',AP_list{i}) == 0
        plot_col = plot_col + 1;
        subplot(2,ceil(plot_cols/2),plot_col);

        % set y axis direction
        if peaks_polarity(i) == -1
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
title('Averaged of All');
hold on;
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
firing_rate_traces = zeros(num_windows, nrois);
t_window = (0:window:num_windows-1)+0.5 * window;

% 统计不为空的trace数目
figure();
shift_firingrate = 0;
for i = 1:nrois
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
xlabel('Time(s)')

fig_filename = fullfile(save_path, '7_firing_rate.fig');
png_filename = fullfile(save_path, '7_firing_rate.png');

saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');

%% Heat Map
heatmap_corrtest = zeros(nrois);
heatmap_speartest = zeros(nrois);

for i = 1:nrois
    for j = 1:nrois
        heatmap_corrtest(i, j) = corr(SNR_traces(:,i), SNR_traces(:,j));
        heatmap_speartest(i, j) = corr(SNR_traces(:,i), SNR_traces(:,j), 'Type', 'Spearman');
    end
end

figure;
subplot(1, 2, 1);
imagesc(heatmap_corrtest);
colorbar;
colormap('jet');  % 设定颜色图
axis square;
ylabel('Voltage Index 1');
xlabel('Voltage Index 2');
title('Correlation Coefficient Heatmap');

subplot(1, 2, 2);
imagesc(heatmap_speartest);
colorbar;
colormap('jet');  % 设定颜色图
axis square;
ylabel('Voltage Index 1');
xlabel('Voltage Index 2');
title('Spearman Correlation Coefficient Heatmap');

fig_filename = fullfile(save_path, '8_Correlation_heatmap.fig');
png_filename = fullfile(save_path, '8_Correlation_heatmap.png');
mat_filename = fullfile(save_path, '8_Correlation_heatmap.mat');

saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');
save(mat_filename,'heatmap_corrtest',"heatmap_speartest");

%% Save parameter
% 定义保存路径和文件名
save_filename = fullfile(save_path, '-1_workspace_variables.mat');

% 保存当前工作区中的所有变量到.mat文件
clear movie;
save(save_filename);
