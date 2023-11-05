% AP_ANALYSIS POWERED BY LIU-YANG LUORONG.
% version 2

% NEED FUNCTIONS:
% load_movie, create_map, highpassfilter, select_ROI

% NEED TOOLBOX:
% Image Processing Toolbox, Curve Fitting Toolbox, Signal Processing Toolbox

clear;
clc;

%% Loading raw data
t1 = tic; % Start a timer
nowtime = string(datestr(now, 'yyyy-mm-dd_HH-MM-SS'));
fprintf('Loading...\n')

% ↓↓↓↓↓-----------Prompt user for define path-----------↓↓↓↓↓
% support for folder, .tif, .tiff, .bin.
folder_path = 'D:\Temple\20230810-170226recordPVH_2000';
file_name = '';  % must add format.
% ↓↓↓↓↓-----------Prompt user for frame rate------------↓↓↓↓↓
freq = 400; % Hz
% -----------------------------------------------------------

% defined bin 1

% read path
file_path = fullfile(folder_path, file_name);
split_path = split(file_name, '.');
if length(split_path)>1
    % when read a file
    file_extension = string(split_path(end));
    save_path =  fullfile(folder_path, [cell2mat(split_path(1)),'_Analysis'], nowtime);
else
    % when read a folder
    file_extension = 'tif';
    save_path = fullfile([folder_path, '_Analysis'] ,nowtime);
end
mkdir(save_path);

% Load image file
[movie, ncols, nrows, nframes] = load_movie(file_path,file_extension);

% Define parameters
dt = 1 / freq; % Calculate time axis
colors = [lines(7);hsv(5);spring(3);winter(3);gray(3)];
t = (1:nframes) * dt;
map = [];
mask = [];

t2 = toc(t1); % Get the elapsed time
fprintf('Finished loading movie after %d s\n',round(t2))
%% Save loaded movie (optional)
t1 = tic; % Start a timer
fprintf('Saving...\n')

raw_filename = fullfile(save_path, '0_Raw_data.mat');
save(raw_filename,"movie",'ncols','nrows','nframes','freq');

t2 = toc(t1); % Get the elapsed time
fprintf('Finished saving movie after %d s\n',round(t2))
%% Create a map (optional)
t1 = tic; % Start a timer

% if SNR is low, please large the bin.
bin = 2;
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
%% Load saved map (optional)

map_filename = fullfile(folder_path, '0_Sensitivity_Map.mat');
map = load(map_filename);
map = map.map;

fprintf('Finished map loading\n')
%% Load saved ROI (optional)

mask_filename = fullfile(folder_path, '1_raw_ROI.mat');
mask = load(mask_filename);
mask = mask.rois;

%% Select ROI
% with or wihout Mask and Map

[rois, traces] = select_ROI(movie, nrows, ncols, t, colors, mask, map);

fig_filename = fullfile(save_path, '1_raw_trace.fig');
png_filename = fullfile(save_path, '1_raw_trace.png');
roi_filename = fullfile(save_path, '1_raw_ROI.mat');

saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');
save(roi_filename, 'rois');


%% Photobleaching correction (optional)

traces_corrected = highpassfilter(traces, freq);
fprintf('Finished highpass filter\n')

%% Photobleaching correction
traces_input = zeros(size(traces,1),size(traces,2)-1);
background = traces(:,end);

% remove background
for i = 1: size(traces,2)-1
    traces_input(:,i) = traces(:,i) - background;
end

% fit
[traces_corrected, fitted_curves] = fit_exp1(traces_input, freq);

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

%% Peak finding
% AP threshold and polarity
peak_polarity =  zeros(1,length(rois));
peak_threshold = zeros(1,length(rois));

% set peak finding
AP_window_width = 40 ; % number of frames to for AP window (defined = 40)
MinPeakDistance = 20 * dt; % (defined = 20 * dt)
peaks_amplitude = cell(1, length(rois)-1);
peaks_index = cell(1, length(rois)-1);
peaks_sensitivity = cell(1, length(rois)-1);

% Plot mean intensity trace and set threshold to find peak
fig = figure();
set(fig,'Position',get(0,'Screensize'));
for i = 1:length(rois)-1
    clf;
    title(sprintf('ROI %d',i));
    hold on;

    % polarity judge
    if abs(min(traces_corrected(:,i))-mean(traces_corrected(:,i))) < max(abs(traces_corrected(:,i)) - mean(traces_corrected(:,i)))
        peak_polarity(i) = 1;
    else
        peak_polarity(i) = -1;
    end

    % plot trace
    plot_trace = traces_corrected(:,i) * peak_polarity(i);
    plot(t,  plot_trace ,'Color',colors(i,:));

    % set threshold
    [~,peak_threshold(i)] = ginput(1);
    plot(t,ones(1,length(t)).*peak_threshold(i),'Color',colors(i,:),'LineWidth',2);
    hold off;
    pause(0.2);

    % find peak
    MinPeakProminence = (max(plot_trace)-mean(plot_trace))*0.5; % (define factor = 0.7)
    [peak_y,peak_x] =  findpeaks(plot_trace, 'MinPeakProminence', MinPeakProminence ,'MinPeakHeight',peak_threshold(i));
    peaks_index{i} = peak_x;
    current_trace = traces_input(:,i);
    peaks_amplitude{i} = current_trace(peak_x);
    peaks_sensitivity{i} = peak_y;
end
close;

% Plot summary
fig = figure();
set(fig,'Position',get(0,'Screensize'));
for i = 1:length(rois)-1
    if ceil((length(rois)-1)/4) > 1
    subplot(ceil((length(rois)-1)/4),4,i); % each line for 4 ROI
    else
    subplot(1,length(rois)-1,i)
    end
    title(sprintf('ROI %d',i));
    hold on;

    % plot trace
    trace = traces_corrected(:,i) * peak_polarity(i);
    plot(t,  trace ,'Color',colors(i,:));
    hold on;

    % plot threshold
    plot(t,ones(1,length(t)).*peak_threshold(i),'Color',colors(i,:),'LineWidth',2);
    hold on;

    % plot peak
    plot(peaks_index{i}.*dt,peaks_sensitivity{i},'v','Color',colors(i,:),'MarkerFaceColor',colors(i,:));
end
sgtitle('Peak finding');
hold on;

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
smooth_trace = zeros(nframes,length(rois));

% plot
sensitivity_axe = subplot(1,2,1);
title('Sensitivity');
hold on;
SNR_axe = subplot(1,2,2);
title('SNR');
hold on;
shift_sensitivity = zeros(size(rois)); % plot shift
shift_SNR = zeros(size(rois));

for i = 1 : length(rois)-1
    % Calculate Sensitivity
    sentivity_trace(:,i) = traces_corrected(:,i);

    % Calculate SNR
    [SNR_trace(:,i),smooth_trace(:,i)]  = calculate_SNR(traces_corrected(:,i), peak_threshold(i), peak_polarity(i));

    % add plot shift
    if i > 1
        shift_sensitivity(i) = max(sentivity_trace(:,i-1))-min(sentivity_trace(:,i-1));
        shift_SNR(i) = max(SNR_trace(:,i-1))-min(SNR_trace(:,i-1));
    else
        shift_sensitivity(i) = 0;
        shift_SNR(i) = 0;
    end

    % plot sensitivity
    plot(t,sentivity_trace(:,i) + peak_polarity(i)*sum(shift_sensitivity(1:i)),'Color', colors(i,:),'Parent',sensitivity_axe);
    hold(sensitivity_axe,'on');
    % plot smooth curve
    plot(t,smooth_trace(:,i) + peak_polarity(i)*sum(shift_sensitivity(1:i)), ...
        'r','Parent',sensitivity_axe,'LineWidth',2);
    hold(sensitivity_axe,'on');

%     
    % plot SNR
    plot(t,SNR_trace(:,i) + peak_polarity(i)*sum(shift_SNR(1:i)),'Color', colors(i,:),'Parent',SNR_axe);
    hold(SNR_axe,'on');

end

legend(sensitivity_axe, 'Sensitivity', 'Fitted Curve');

fig_filename = fullfile(save_path, '4_Sensitivity_figure.fig');
png_filename = fullfile(save_path, '4_Sensitivity_figure.png');

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
    each_trace_smooth = smooth_trace(:,i);
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
        T = table(number_i, amp_i, FWHM_i.*1000, sensitivity_i, SNR_i, ...
            'VariableNames', {'Number', 'Amplitude', 'FWHM (ms)', 'Sensitivity', 'SNR'});

        % 将表格写入Excel的一个新工作表
        sheet_name = string(['ROI ' num2str(i)]);
        writetable(T,table_name, 'Sheet', sheet_name);

        % save average value
        avg_amp(i) = mean(amp_i);
        avg_FWHM(i) = mean(FWHM_i) .* 1000;
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
%% 
