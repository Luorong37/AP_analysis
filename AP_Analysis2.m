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
% See also calculate_firing_rate, calculate_FWHM, create_map, calculate_SNR, fit_exp1, highpassfilter, select_ROI

%% Loading raw data
nowtime = string(datetime( 'now'));
% Replace colons with hyphens to get the desired output format
nowtime = strrep(nowtime , ':', '-');
fprintf('Loading...\n')

% ↓↓↓↓↓-----------Prompt user for define path-----------↓↓↓↓↓
% support for folder, .tif, .tiff, .bin.
folder_path = 'U:\Luorong\24.07.12_dueplex\1：10\';
file = '50%488-1';  % must add format.do not add '\' at last
% ↓↓↓↓↓-----------Prompt user for frame rate------------↓↓↓↓↓
freq = 400; % Hz
% -----------------------------------------------------------

% Create an analysis folder
function [file_path, save_path] = create_folder(folder_path, file, nowtime)
% read path
file_path = fullfile(folder_path, file);
if isfolder(file_path)
    file_name = file;
    %file_dir = dir(file_path);
    %[~, ~, file_extension] = fileparts(file_dir(3).name);
else
    % [~, file_name, file_extension] = fileparts(file_path)
    [~, file_name, ~] = fileparts(file_path);
end

% create a folder for analysis
save_path = fullfile(folder_path, [file_name, '_Analysis'], nowtime);
mkdir(save_path);
end

[file_path, save_path] = create_folder(folder_path, file, nowtime);

% Load image file
[movie, ncols, nrows, nframes] = load_movie(file_path);

% Presetting
function [dt, colors, t, map, mask, avg_image, options] = presetting(freq, nframes, movie, ncols, nrows, save_path)

% Define parameters
dt = 1 / freq; % Calculate time axis
% colors = [lines(7);hsv(5);spring(3);winter(3);gray(3)];
colors = lines(100);
t = (1:nframes) * dt;
options.colors = colors;
map = [];
mask = [];
movie_vol_2D = reshape(mean(movie,2), ncols, nrows, []);
avg_image  = (movie_vol_2D - min(movie_vol_2D(:))) / (max(movie_vol_2D(:)) - min(movie_vol_2D(:)));


% Save code
code_path = fullfile(save_path,'Code');
mkdir(code_path);
currentScript = which("AP_Analysis_Auto.m");
% 获取当前脚本依赖的所有文件
[requiredFiles, ~] = matlab.codetools.requiredFilesAndProducts(currentScript);
% 复制当前脚本和所有依赖文件到目标文件夹
for k = 1:length(requiredFiles)
    [~, name, ext] = fileparts(requiredFiles{k});
    copyfile(requiredFiles{k}, fullfile(code_path, [name, ext]));
end
end

[dt, colors, t, mask, avg_image, code_path] = presetting(dt, freq, nframes, mask, movie, ncols, nrows, save_path);

% 提示完成
fprintf('All codes have been copied to %s\n', code_path);

%% ----------------------Optional part------------------------
% This part can load previous selected Mask
% Save data
save_stack = false;
if save_stack
    create_tiff_stack(file_path);
end

preload = questdlg('Load previous data?','load data','ROIs','No','Cancel');
switch preload
    case 'ROIs'
        [roi_filename,roi_foldername] = uigetfile(save_path);
        rois_data = load(fullfile(roi_foldername,roi_filename));
        rois = rois_data.rois;
        mask = rois.bwmask;
end
% -----------------------------------------------------------
%% Create a map (optional)
t1 = tic; % Start a timer
fprintf('Creating a map...\n')
% if the map cannot figure out active cells, please large the bin.
bin = 4; % defined bin = 4
[quick_map] = create_map(movie, nrows, ncols, bin);
map = quick_map;

function plot_map(quick_map, save_path, map)
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
end

plot_map(quick_map, save_path, map);

t2 = toc(t1); % Get the elapsed time
fprintf('Finished mask creating after %d s\n',round(t2))

%% Select ROI
t1 = tic; % Start a timer

% with or wihout Mask and Map
if exist('preload','var')
    if any(strcmp(preload,{'No',''} ))
        [rois, traces] = select_ROI(movie, ncols, nrows,  mask, map);
    else
        [~, traces] = select_ROI(movie, ncols, nrows, mask, map);
    end
else
    [rois, traces] = select_ROI(movie, ncols, nrows, mask, map);
end
nrois = max(rois.bwmask,[],'all');
bwmask = rois.bwmask;
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
close(gcf);
%% Signal Process %%

% Background Correction
fprintf('Removing Background...')
[background, background_fitted, traces_bgcorr, traces_bgfitcorr, background_mask]...
    = remove_background(movie, ncols, nrows, rois);
roi_filename = fullfile(save_path, '1_background_ROI.mat');
save(roi_filename, 'background','background_fitted','background_mask','traces_bgcorr','traces_bgfitcorr');
fprintf(' Finished\n')

%% Bleaching Correction

function plot_corrected(traces_corrected, traces, baseline, colors, save_path)
% plot
fig = figure();
set(fig,'Position',get(0,'Screensize'));
fit_axe = subplot(3,1,1);
fited_axe = subplot(3,1,2);
baseline_axe = subplot(3,1,3);

% plot
for i = 1: size(traces_corrected,2)
    plot(traces(:,i),'Color',colors(i,:),'Parent',fit_axe);
    hold(fit_axe, 'on');
    plot(traces_corrected(:,i),'Color',colors(i,:),'Parent',fited_axe);
    hold(fited_axe, 'on');
    plot( baseline(:,i),'Color',colors(i,:),'LineWidth',2,'Parent',baseline_axe);
    hold(baseline_axe, 'on');

end
hold off;

% note
title(fit_axe, 'Original and Fitted Curves');
title(fited_axe, 'Corrected Traces');
title(baseline_axe, 'Baseline');
legend(fit_axe, 'Original Trace', 'Fitted Curve');

fig_filename = fullfile(save_path, '2_fitted_trace.fig');
png_filename = fullfile(save_path, '2_fitted_trace.png');

saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');
end

fprintf('Correcting Bleaching...\n')
[traces_corrected,baseline] = highpass_bleach_remove(traces_bgfitcorr,500);
plot_corrected(traces_corrected, traces, baseline, colors, save_path);
fprintf('Finished\n');

%% Wavelet process

% wavelet降噪
Dnmethods = 'FDR';
Dnlevel = 8;
fprintf('Wavelet Denoising...\')
traces_denoised = wdenoise(traces_corrected, Dnlevel ,DenoisingMethod=Dnmethods);

traces_filename = fullfile(save_path, '2_processed_traces.mat');
save(traces_filename,"traces_denoised", "traces_original", 'traces_corrected', 'baseline', 'Dnmethods','Dnlevel')
fprintf(' Finished\n');

%% AP Processing %%

% Peak finding
[peaks_index, peaks_amplitude, peaks_polarity] = peak_finding_auto(traces_denoised , save_path,'parts',3,'MinPeakProminence_factor',0.4);

%% Calculate sensitivity and SNR
traces_sensitivity = traces_corrected./baseline;

noise = traces_corrected-traces_denoised;
traces_SNR = traces_corrected./std(noise);

peaks_sensitivity = cell(size(peaks_index));
peaks_SNR = cell(size(peaks_index));
for i = 1: nrois
    peaks_sensitivity{i} = traces_sensitivity(peaks_index{i},i);
    peaks_SNR{i} = traces_SNR(peaks_index{i},i);
end

% plot traces

function plot_traces(nrois, traces_sensitivity, traces_SNR, avg_image, bwmask, t, save_path)
% set up
fig = figure();
set(fig,'Position',get(0,'Screensize'));
% plot fluorescent image
imshow(avg_image,[min(avg_image,[],'all'),max(avg_image,[],'all')]);
roiseq = unique(sort(bwmask(:)));
for i = 1:nrois
    roi = (bwmask == roiseq(i+1));
    boundary = cell2mat(bwboundaries(roi));
    plot(boundary(:, 2), boundary(:, 1), 'LineWidth', 2, 'Parent', f_axe);hold on;
    text(mean(boundary(:, 2)) + 12, mean(boundary(:, 1)) - 12, num2str(roiseq(i+1)),'FontSize', 12, 'Parent', f_axe); hold on;
end
hold on;
title('Fluorescent Image');

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

save(trace_filename,"traces_sensitivity",'traces_SNR');
saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');
end

plot_traces(nrois, traces_sensitivity, traces_SNR, avg_image, bwmask, t, save_path);

%% Statistic AP

function [AP_list, AP_data]  = AP_statistic(nrois, peaks_index, peaks_amplitude, traces_corrected, traces_sensitivity, traces_SNR, AP_window_width, nframes, dt, peaks_polarity, save_path)
AP_list = cell(1, nrois);

% each trace
for i = 1:nrois % i for trace
    peaks_num = length(peaks_index{i});
    peaks_index_i = peaks_index{i};
    peaks_amp_i = peaks_amplitude{i};
    each_trace_amp = traces_corrected(:,i);
    each_trace_sensitivity = traces_sensitivity(:,i);
    each_trace_SNR = traces_SNR(:,i);
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
        Sensitivity = AP_sensitivity(AP_window_width+1)*100 ;
        SNR = abs(AP_SNR(AP_window_width+1));
        FWHM = calculate_FWHM(AP_amp, dt, peaks_polarity{i});

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

AP_data.amp = {};
AP_data.FWHM = {};
AP_data.sensitivity = {};
AP_data.SNR = {};


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
        T = table(number_i, amp_i, FWHM_i, sensitivity_i, SNR_i, peaks_index{i}, ...
            'VariableNames', {'Number', 'Amplitude', 'FWHM (ms)', 'Sensitivity', 'SNR', 'Index'});

        % 将表格写入Excel的一个新工作表
        sheet_name = string(['ROI ' num2str(i)]);
        writetable(T,table_name, 'Sheet', sheet_name);

        % save average value
        avg_amp(i) = mean(amp_i,'omitmissing');
        avg_FWHM(i) = mean(FWHM_i,'omitmissing');
        avg_sensitivity(i) = mean(sensitivity_i,'omitmissing');
        avg_SNR(i) = mean(SNR_i,'omitmissing');
        AP_number(i) = number_i(end);
        ROI_number(i) = i;
        AP_data.amp{i} = amp_i;
        AP_data.FWHM{i} = FWHM_i;
        AP_data.sensitivity{i} = sensitivity_i;
        AP_data.SNR{i} = SNR_i;
    end
end


figure()
FWHM_axe = subplot(1,4,2);hold on;xlim([0,nrois+1]);
sensitivity_axe = subplot(1,4,3);hold on;xlim([0,nrois+1]);
SNR_axe = subplot(1,4,4);hold on;xlim([0,nrois+1]);
num_axe = subplot(1,4,1);hold on;xlim([0,nrois+1]);
sgtitle('AP statistic');
bar(avg_FWHM,'Parent',FWHM_axe);hold on;
xlabel('ROI number','Parent',FWHM_axe);
ylabel('FWHM (ms)','Parent',FWHM_axe);

bar(avg_sensitivity,'Parent',sensitivity_axe);hold on;
xlabel('ROI number','Parent',sensitivity_axe);
ylabel('Sensitiviy','Parent',sensitivity_axe);

bar(avg_SNR,'Parent',SNR_axe);hold on;
xlabel('ROI number','Parent',SNR_axe);
ylabel('SNR','Parent',SNR_axe);

bar(AP_number,'Parent',num_axe);hold on;
xlabel('ROI number','Parent',num_axe);
ylabel('AP number','Parent',num_axe);


fig_filename = fullfile(save_path, '5_AP statistic.fig');
png_filename = fullfile(save_path, '5_AP statistic.png');

saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');

T_ave = table(ROI_number, AP_number, avg_amp, avg_FWHM, avg_sensitivity, avg_SNR, ...
    'VariableNames', {'ROI Number','AP Number', 'Average Amplitude', 'Average FWHM (ms)', 'Average Sensitivity', 'Average SNR'});
writetable(T_ave, table_name, 'Sheet', 'Average');
AP_data_filename = fullfile(save_path, 'AP_data.mat');
save(AP_data_filename, "AP_data",'AP_list')
fprintf('Finished statistic AP\n')
end

AP_window_width = 10 ; % number of frames to for AP window (defined = 40)
[AP_list, AP_data] = AP_statistic(nrois, peaks_index, peaks_amplitude, traces_corrected, traces_sensitivity, traces_SNR, AP_window_width, nframes, dt, peaks_polarity, save_path);
fprintf('AP_data.xlsx saved.\n');

%% Plot average AP sensitivity
figure();
% 统计不为空的trace数目
plot_cols = sum(cellfun('isempty',AP_list)==0)+1;
plot_col = 0;
% Initialize arrays to store all traces for final average calculation
trace_AP_mean = [];

for i = 1:nrois % i for trace
    %判断是否为有AP的trace
    peaks_num = length(peaks_index{i});
    if cellfun(['isempt' ...
            'y'],AP_list{i}) == 0
        plot_col = plot_col + 1;
        subplot(2,ceil(plot_cols/2),plot_col);
        set(gca,'color','none');

        % set y axis direction
        if peaks_polarity{i} == -1
            set(gca,'YDir','reverse')
            hold on;
        end

        % get each AP
        AP_i = zeros(peaks_num, AP_window_width*2+1);
        for j = 1:peaks_num
            each_AP = AP_list{i}{j};
            AP_i(j,:) = each_AP.AP_sensitivity;
            % plot each AP
            % plot((1:AP_window_width*2+1)*dt, each_AP.AP_sensitivity','Color',[0.8 0.8 0.8]);
            hold on;
        end

        AP_mean = mean(AP_i, 1, 'omitnan');
        AP_sd = std(AP_i, 0, 1, 'omitnan');

        % plot average AP for each trace
        subplot(2,ceil(plot_cols/2),plot_col);
        plot((1:AP_window_width*2+1)*dt, AP_mean,'Color',colors(i,:),'LineWidth',2);
        hold on;
        title(sprintf('ROI %d\n',i));
        fill([(1:AP_window_width*2+1)*dt, fliplr((1:AP_window_width*2+1)*dt)], ...
            [AP_mean + AP_sd, fliplr(AP_mean - AP_sd)], ...
            colors(i,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        hold on;

        % plot average AP for average trace

        % plot((1:AP_window_width*2+1)*dt, mean(AP_i,1,'omitnan'),'Color',colors(i,:),'LineWidth',1);
        % hold on;
        trace_AP_mean = [trace_AP_mean; AP_mean*peaks_polarity{i}];
    end
end

if ~isempty(trace_AP_mean)
    overall_mean = mean(trace_AP_mean, 1, 'omitnan');
    overall_sem = std(trace_AP_mean, 0, 1, 'omitnan') / sum(~cellfun('isempty', AP_list));
    subplot(2, ceil(plot_cols / 2), plot_cols);set(gca,'color','none');hold on;
    fill([(1:AP_window_width*2+1)*dt, fliplr((1:AP_window_width*2+1)*dt)], ...
        [overall_mean + overall_sem, fliplr(overall_mean - overall_sem)], ...
        [0.8 0.8 0.8], 'EdgeColor', 'none');

    title('Averaged of All');
    plot((1:AP_window_width*2+1) * dt, overall_mean, 'Color', 'k', 'LineWidth', 1);
    hold on;
end

sgtitle('Averaged Sensitivity');
hold on;
fig_filename = fullfile(save_path, '6_average_AP_sensitivity.fig');
png_filename = fullfile(save_path, '6_average_AP_sensitivity.png');
mat_filename = fullfile(save_path, '6_average_AP_sensitivity.mat');
saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');
save(mat_filename, 'AP_list', 'peaks_index', 'nrois', 'AP_window_width', 'dt', ...
    'peaks_polarity', 'colors', 'trace_AP_mean', 'overall_mean', 'overall_sem');

%%
% Plot average AP SNR with SD
figure();
plot_cols = sum(cellfun('isempty', AP_list) == 0) + 1;
plot_col = 0;
% Initialize arrays to store all traces for final average calculation
trace_AP_SNR_mean = [];

for i = 1:nrois
    peaks_num = length(peaks_index{i});
    if ~isempty(AP_list{i})
        plot_col = plot_col + 1;
        subplot(2, ceil(plot_cols / 2), plot_col);
        set(gca,'color','none');

        % Set y axis direction if necessary
        if peaks_polarity{i} == -1
            set(gca, 'YDir', 'reverse');
            hold on;
        end

        % Get each AP and compute mean & SEM
        AP_i = zeros(peaks_num, AP_window_width*2 + 1);
        for j = 1:peaks_num
            each_AP = AP_list{i}{j};
            AP_i(j,:) = each_AP.AP_SNR;
            %plot((1:AP_window_width*2+1) * dt, each_AP.AP_SNR', 'Color', [0.8 0.8 0.8]);
            hold on;
        end

        AP_mean = mean(AP_i, 1, 'omitnan');
        AP_sd = std(AP_i, 0, 1, 'omitnan');

        % Plot mean with Sd
        plot((1:AP_window_width*2+1) * dt, AP_mean, 'Color', colors(i,:), 'LineWidth', 2);
        fill([(1:AP_window_width*2+1)*dt, fliplr((1:AP_window_width*2+1)*dt)], ...
            [AP_mean + AP_sd, fliplr(AP_mean - AP_sd)], ...
            colors(i,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        title(sprintf('ROI %d', i));
        hold on;

        % % Plot average AP for all traces
        % subplot(2, ceil(plot_cols / 2), plot_cols);
        % plot((1:AP_window_width*2+1) * dt, AP_mean*peaks_polarity(i), 'Color', [0.8 0.8 0.8]);
        % hold on;

        % Collect traces for final average calculation
        trace_AP_SNR_mean = [trace_AP_SNR_mean; AP_mean*peaks_polarity{i}];
    end

end

% Plot overall average and SEM in the last subplot
if ~isempty(trace_AP_SNR_mean)
    overall_mean = mean(trace_AP_SNR_mean, 1, 'omitnan');
    overall_sem = std(trace_AP_SNR_mean, 0, 1, 'omitnan') / sum(~cellfun('isempty', AP_list));
    subplot(2, ceil(plot_cols / 2), plot_cols);hold on;
    set(gca,'color','none');
    fill([(1:AP_window_width*2+1)*dt, fliplr((1:AP_window_width*2+1)*dt)], ...
        [overall_mean + overall_sem, fliplr(overall_mean - overall_sem)], ...
        [0.8 0.8 0.8], 'EdgeColor', 'none');

    title('Averaged of All');
    plot((1:AP_window_width*2+1) * dt, overall_mean, 'Color', 'k', 'LineWidth', 1);
    hold on;
end
sgtitle('Average SNR with SD');
fig_filename = fullfile(save_path, '6_average_AP_SNR_with_SD.fig');
png_filename = fullfile(save_path, '6_average_AP_SNR_with_SD.png');

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

%% Save parameter
% 定义保存路径和文件名
save_filename = fullfile(save_path, '-1_workspace_variables.mat');

% % 保存当前工作区中的所有变量到.mat文件
% clear movie;
% save(save_filename);



%%


