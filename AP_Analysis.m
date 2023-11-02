clear all;
clc;

%% Loading
% Need Functions:
% load_tif_movie, readBinMov,select_ROI_and_draw_trace, find_cell

% Need Toolbox:
% Image Processing Toolbox, Curve Fitting Toolbox, Signal Processing Toolbox

% Define folder path and file extension
folder_path = 'E:\20231025 zsy\JF552\cell2\CQ1'; 
file_name = '\movie.bin'; % start with '\' if do not save at same folder
% colors = {'r', 'g', 'b', 'm', 'c', 'y', 'k'}; % Define a set of colors
colors = lines;
% Prompt user for frame rate in Hz
Hz = 400;

% read file
file_path = [folder_path file_name];
split_path = split(file_name, '.');
if length(split_path)>1
    file_extension = string(split_path(2));
else
    file_extension = 'tif';
end
% save_path = [folder_path cell2mat(split_path(1)) '_analysis'];
save_path = [folder_path '\Analysis'];
format =string(cell2mat(split_path(end))) ;
mkdir(save_path);

% Calculate time axis
dt = 1 / Hz;

% Load image file
[intensity_time_series,  num_cols,num_rows,num_frames] = load_movie(file_path,file_extension);
t = (1:num_frames) * dt;
Map = [];
% im = cell(num_frames, 1);
% for i = 1:num_frames
%     im{i} = reshape(uint16(intensity_time_series(:,i)), num_rows, num_cols);
% end

fprintf('Finished loading movie\n')

%% Create a map quickly (optional)
t1 = tic; % Start a timer

% if SNR is low, large the bin.
bin = 1;
[quick_map] = quick_find_cell(intensity_time_series, num_rows, num_cols, bin, Hz);
Map = quick_map;

% Visualize correlation coefficients as heatmap
figure()
imagesc(quick_map);
title('Sensitivity Map');
axis image;
colorbar;
fig_filename = fullfile(save_path, '0_Sensitivity_Map.fig');
png_filename = fullfile(save_path, '0_Sensitivity_Map.png');

saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');

t2 = toc(t1); % Get the elapsed time
fprintf('Finished mask creating after %d s\n',round(t2))

%% Select ROI
% with or wihout Mask and Map
% Mask = load('E:\20231025 zsy\JF552\cell2\qucik_map.mat');
% mask_map = {Mask.quick_map};
% if isempty(Cn)
%     Map = [];
% else
%     Map = Cn;
% end
if isempty(Mask) && isempty(Map)
    [rois, traces] = select_ROI_v2(intensity_time_series, num_rows,num_cols,t, colors);
elseif ~isempty(Mask)
    [rois, traces] = selected_ROI_v2(intensity_time_series,num_rows,num_cols,t, colors, mask_map);
elseif ~isempty(Map)
    [rois, traces] = select_ROI_with_map_v2(intensity_time_series,num_rows,num_cols,t, colors, Map);
end

fig_filename = fullfile(save_path, '1_raw_trace.fig');
png_filename = fullfile(save_path, '1_raw_trace.png');

saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');

%% Photobleaching correction
traces_highpassfilted = highpassfilter(traces, Hz);

% polarity judge:
polarity = -1;

% AP threshold
peak_min_height = zeros(1,length(rois));

fig = figure();
set(fig,'Position',get(0,'Screensize'));
% Display ROI selection figure in the same figure
subplot(1,2,1);
imshow(imadjust(uint16(mean(reshape(intensity_time_series, num_cols,num_rows,[]),3))));
hold on;
title(sprintf('Selected ROIs'));
for i = 1:length(rois)
    boundary = bwboundaries(rois{i});
    plot(boundary{1}(:,2), boundary{1}(:,1),'Color', colors(mod(i-1, length(colors))+1,:), 'LineWidth', 2);
end


% Plot mean intensity over time
subplot(1,2,2);
trace_axe = gca;
hold on;
for i = 1:length(rois)
    plot(t, traces_highpassfilted(:,i) * polarity ,'Color',colors(i,:));
    if i < length(rois)
    [~,peak_min_height(i)] = ginput(1);
    plot(t,ones(1,length(t)).*peak_min_height(i),'Color',colors(i,:),'LineWidth',2);
    end
end
hold off;
xlabel('Time (s)');
ylabel('Mean intensity');
title(sprintf('Mean intensity of selected ROI after high-pass filter\nSet threshold of peaks'));

% Create a legend label for each ROI
legend_labels = cell(1, length(rois));
for i = 1:length(rois)
    legend_labels{i} = sprintf('ROI %d', i);
end
legend(legend_labels);

% Save
fig_filename = fullfile(save_path, '2_filtered_trace.fig');
png_filename = fullfile(save_path, '2_filtered_trace.png');

saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');

%% Peak finding
% from stimulated_AP_multicycles_V3.1

figure();
AP_window_size = 100 ; % number of frames to for AP window
MinPeakDistance = 20 * dt;
% MinPeakProminence = 20 ;

title('Peak finding');
hold on;
peaks_amplitude = cell(1, length(rois)-1);
peaks_index = cell(1, length(rois)-1);

% set threshold and find peak
peak_threshold = 0;
peak_min_distance = 5*dt;
% peak_min_height = 10;

for i = 1:length(rois)
    if i < length(rois)
        trace = traces_highpassfilted(:,i) * polarity;

        % set Prominence
        MinPeakProminence = max(trace)*0.3;

%         findpeaks(trace,t, 'MinPeakProminence', MinPeakProminence ,'MinPeakHeight',peak_min_height(i)); 
%         [peaks_amplitude{i}, peaks_index{i}] =  findpeaks(trace, 'MinPeakProminence', MinPeakProminence ,'MinPeakHeight',peak_min_height(i));
        findpeaks(trace,t, 'MinPeakProminence', MinPeakProminence ,'MinPeakHeight',peak_min_height(i)); 
        [peaks_amplitude{i}, peaks_index{i}] =  findpeaks(trace, 'MinPeakProminence', MinPeakProminence ,'MinPeakHeight',peak_min_height(i));

        
        hold on;
    else
        plot(t, traces_highpassfilted(:,i) * polarity ,'Color',colors(i,:));
    end
end

fig_filename = fullfile(save_path, '3_peak_finding.fig');
png_filename = fullfile(save_path, '3_peak_finding.png');

saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');



%% Sensitivity and SNR Analysis
% show fluorescent imageing
figure();
subplot(1,3,1);
title(sprintf('Selected ROIs'));
imshow(imadjust(uint16(mean(reshape(intensity_time_series, num_cols,num_rows,[]),3))));
hold on;
title(sprintf('Selected ROIs'));

for i = 1:length(rois)
    boundary = bwboundaries(rois{i});
    plot(boundary{1}(:,2), boundary{1}(:,1),'Color', colors(mod(i-1, length(colors))+1,:), 'LineWidth', 2);
end

% plot and calculate
intensity_trace = zeros(num_frames,length(rois));
sentivity_trace = zeros(num_frames,length(rois));
SNR_trace = zeros(num_frames,length(rois));

subplot(1,3,2);
title('Sensitivity');
subplot(1,3,3);
title('SNR');
shift_sensitivity = zeros(size(rois)); % plot shift
shift_SNR = zeros(size(rois));

for i = 1 : length(rois)
    intensity_trace(:,i) = traces(:,i) - traces(:,end);
    % calculate sensitivity
    sentivity_trace(:,i) = traces_highpassfilted(:,i)  ./ intensity_trace(:,i);
    [xData, yData] = prepareCurveData( [], traces_highpassfilted(:,i));
    ft = fittype( 'poly1' );
    [fitresult, gof] = fit( xData, yData, ft );
    RMSE_mean=gof.rmse;
    SNR_trace(:,i) = traces_highpassfilted(:,i) ./ RMSE_mean;

    if i > 1
        shift_sensitivity(i) = max(sentivity_trace(:,i-1))-min(sentivity_trace(:,i-1));
        shift_SNR(i) = max(SNR_trace(:,i-1))-min(SNR_trace(:,i-1));
    else
        shift_sensitivity(i) = 0;
        shift_SNR(i) = 0;
    end
    % plot sensitivity
    subplot(1,3,2);
    plot(t,sentivity_trace(:,i) + polarity*sum(shift_sensitivity(1:i)));
    hold on;
    % plot SNR
    subplot(1,3,3);
    plot(t,SNR_trace(:,i) + polarity*sum(shift_SNR(1:i)));
    hold on;

end



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
    each_trace = traces_highpassfilted(:,i) * polarity;
    each_trace_sensitivity = sentivity_trace(:,i);
    each_trace_SNR = SNR_trace(:,i);

    AP_list{i} = cell(1, length(peaks_index{i}));
    
    % each peak
    for j = 1:peaks_num % j for peak
        peak_index_ij = peaks_index_i(j); 
        peak_amp_ij = peaks_amp_i(j);

        % 防止出边缘
        AP_start_index = max(1, peak_index_ij - AP_window_size);
        AP_end_index = min(num_frames, peak_index_ij + AP_window_size);
        AP_index = AP_start_index : AP_end_index;

        % 检索强度，灵敏度，SNR
        AP_amp = each_trace(AP_start_index:AP_end_index)';
        AP_sensitivity = each_trace_sensitivity(AP_start_index:AP_end_index)';
        AP_SNR = each_trace_SNR(AP_start_index:AP_end_index)';
        % calculate SNR
        % noise

        % fill NaN
        if 0 > peak_index_ij - AP_window_size
            AP_amp = [NaN(1,0 - (peak_index_ij - AP_window_size)+1), AP_amp];
            AP_sensitivity = [NaN(1,0 - (peak_index_ij - AP_window_size)+1),AP_sensitivity];
            AP_SNR = [NaN(1,0 - (peak_index_ij - AP_window_size)+1),AP_SNR];
        elseif num_frames < peak_index_ij + AP_window_size
            AP_amp = [AP_amp, NaN(1,peak_index_ij + AP_window_size - num_frames)];
            AP_sensitivity = [AP_sensitivity, NaN(1,peak_index_ij + AP_window_size - num_frames)];
            AP_SNR = [AP_SNR NaN(1,peak_index_ij + AP_window_size - num_frames)];
        end
           
        AP_base = mean(AP_amp,'omitnan');
        HM = 0.5 * (AP_base + peak_amp_ij);
        
        % Calculate FWHM
        % 在峰左侧找到半最大值的位置
        [left_below_idx, left_above_idx] = deal([], []);
        if any(AP_amp(1:AP_window_size + 1) < HM)
            left_below_idx = find(AP_amp(1:AP_window_size+1) < HM, 1,'last');
        end
        if any(AP_amp(1:AP_window_size + 1) > HM)
            left_above_idx = left_below_idx - 1 + find(AP_amp(left_below_idx:AP_window_size+1) > HM, 1, 'first');
        end
        
        % 在峰右侧找到半最大值的位置
        [right_above_idx, right_below_idx] = deal([], []);
        if any(AP_amp(AP_window_size + 1 :end) < HM)
            right_below_idx = find(AP_amp(AP_window_size+1:end) < HM, 1, 'first') + AP_window_size;
        end
        if any(AP_amp(AP_window_size + 1:end) > HM)
            right_above_idx = find(AP_amp(AP_window_size+1:right_below_idx) > HM, 1, 'last') + AP_window_size;
        end

        
        % 使用线性插值估算精确的半最大值位置
        if ~isempty(left_below_idx) && ~isempty(left_above_idx)
            idx_left_exact = left_below_idx + (HM - AP_amp(left_below_idx)) / (AP_amp(left_above_idx) - AP_amp(left_below_idx));
        else
            idx_left_exact = [];
        end
        
        if ~isempty(right_above_idx) && ~isempty(right_below_idx)
            idx_right_exact = right_below_idx - (AP_amp(right_above_idx) - HM) / (AP_amp(right_above_idx) - AP_amp(right_below_idx));
        else
            idx_right_exact = [];
        end
        
        % 计算FWHM
        if isempty(idx_left_exact) || isempty(idx_right_exact)
            FWHM = NaN; % 或者可以设定为某个默认值
        else
            FWHM = (idx_right_exact - idx_left_exact) * dt; 
        end

        % Amplitude = abs(AP_base - peak_amp_ij);
        Amplitude = abs(peak_amp_ij);
        Sensitivity = polarity * max(abs(AP_sensitivity));
        SNR = polarity * max(abs(AP_SNR));

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

% 计算所有AP的SNR数据并存储在tables中
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
        sheet_name = ['ROI ' num2str(i)];
        writetable(T, [save_path '\AP_data.xlsx'], 'Sheet', sheet_name);

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
writetable(T_ave, [save_path '\AP_data.xlsx'], 'Sheet', 'Average');

fprintf('Finished statistic AP\n')

%% Plot average AP sensitivity
% plot
figure();
title('Statistic AP');
hold on;

% 统计不为空的trace数目
plot_cols = sum(cellfun('isempty',AP_list)==0)+1;
plot_col = 0;
subplot(1,plot_cols,plot_cols);
title('All Average');
if polarity == -1
    set(gca,'YDir','reverse') 
end
hold on;

for i = 1:length(rois)-1
    %判断是否为有AP的trace
    peaks_num = length(peaks_index{i});
    if cellfun('isempty',AP_list{i}) == 0
        plot_col = plot_col + 1;
        subplot(1,plot_cols,plot_col);
        title(sprintf('ROI %d\nAverage AP sensitivity',i));
        if polarity == -1
            set(gca,'YDir','reverse')
        end
        hold on;
        AP_i = zeros(peaks_num, AP_window_size*2+1);
        % extract each AP
        for j = 1:peaks_num
            each_AP = AP_list{i}{j};
            % save amp
            AP_i(j,:) = each_AP.AP_sensitivity;
            subplot(1,plot_cols,plot_col);
            plot([1:AP_window_size*2+1]*dt, each_AP.AP_sensitivity','Color',[0.8 0.8 0.8]);
            hold on;
%             subplot(1,plot_cols,plot_cols);
%             plot([1:AP_window_size*2+1]*dt, each_AP.AP_sensitivity','Color',[0.8 0.8 0.8]);
%             hold on;
        end
        %plot average AP
        subplot(1,plot_cols,plot_col);
        plot([1:AP_window_size*2+1]*dt, mean(AP_i,1,'omitnan'),'Color',colors(i,:),'LineWidth',2);
        hold on;
        subplot(1,plot_cols,plot_cols);
        plot([1:AP_window_size*2+1]*dt, mean(AP_i,1,'omitnan'),'Color',colors(i,:),'LineWidth',1);
        hold on;
    end
end


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
title('All Average SNR');
if polarity == -1
    set(gca,'YDir','reverse') 
end
hold on;

for i = 1:length(rois)-1
    %判断是否为有AP的trace
    peaks_num = length(peaks_index{i});
    if cellfun('isempty',AP_list{i}) == 0
        plot_col = plot_col + 1;
        subplot(1,plot_cols,plot_col);
        title(sprintf('ROI %d Average AP SNR',i));
        if polarity == -1
            set(gca,'YDir','reverse')
        end
        hold on;
        AP_i = zeros(peaks_num, AP_window_size*2+1);
        % extract each AP
        for j = 1:peaks_num
            each_AP = AP_list{i}{j};
            % save amp
            AP_i(j,:) = each_AP.AP_SNR;
            subplot(1,plot_cols,plot_col);
            plot([1:AP_window_size*2+1]*dt, each_AP.AP_SNR','Color',[0.8 0.8 0.8]);
            hold on;
%             subplot(1,plot_cols,plot_cols);
%             plot([1:AP_window_size*2+1]*dt, each_AP.AP_SNR','Color',[0.8 0.8 0.8]);
%             hold on;
        end
        %plot average AP
        subplot(1,plot_cols,plot_col);
        plot([1:AP_window_size*2+1]*dt, mean(AP_i,1,'omitnan'),'Color',colors(i,:),'LineWidth',2);
        hold on;
        subplot(1,plot_cols,plot_cols);
        plot([1:AP_window_size*2+1]*dt, mean(AP_i,1,'omitnan'),'Color',colors(i,:),'LineWidth',1);
        hold on;
    end
end


fig_filename = fullfile(save_path, '6_average_AP_SNR.fig');
png_filename = fullfile(save_path, '6_average_AP_SNR.png');

saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');

% %%  SNR calcuation
% 
% % Calculate SNR-time trace for each ROI
% traces_snr = zeros(size(traces_filt));
% shift = zeros(size(rois)); % plot shift
% for i = 1:size(traces_filt, 2)
%     std_dev = std(traces_filt(:,i));
%     % Calculate SNR for each time point
%     traces_snr(:,i) = traces_filt(:,i) ./ std_dev;
%     if i > 1
%         shift(i) = max(traces_snr(:,i))-min(traces_snr(:,i));
%     else
%         shift(i) = 0;
%     end
% end
% % Plot SNR-time trace for each ROI
% figure();
% subplot(1,2,2);
% hold on;
% for i = 1:length(rois)
%     plot(t, traces_snr(:,i) + sum(shift(1:i)), 'Color', colors(mod(i-1, length(colors))+1,:));
% end
% hold off;
% xlabel('Time (s)');
% ylabel('SNR');
% title('SNR of selected ROIs over time');
% 
% % Create a legend label for each ROI
% legend_labels = cell(1, length(rois));
% for i = 1:length(rois)
%     legend_labels{i} = sprintf('ROI %d', i);
% end
% legend(legend_labels);
% 
% % Display ROI selection figure in the same figure
% subplot(1,2,1);
% imshow(imadjust(im{1}));
% title('Selected ROIs');
% hold on;
% for i = 1:length(rois)
%     boundary = bwboundaries(rois{i});
%     plot(boundary{1}(:,2), boundary{1}(:,1),'Color', colors(mod(i-1, length(colors))+1,:), 'LineWidth', 1.5);
% end
% hold off;
% 
% fig_filename = fullfile(save_path, '3_SNR_trace.fig');
% png_filename = fullfile(save_path, '3_SNR_trace.png');
% 
% saveas(gcf, fig_filename, 'fig');
% saveas(gcf, png_filename, 'png');


