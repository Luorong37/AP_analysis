function [peak_polarity, peak_threshold, peaks_index, peaks_amplitude, peaks_sensitivity] = ...
                                            peak_finding(traces_corrected)
    %PEAK_FINDING 寻找给定ROI信号中的峰值。
    %
    %   [peak_polarity, peak_threshold, peaks_index, peaks_amplitude, peaks_sensitivity] = 
    %   peak_finding(traces_corrected, t, colors, rois) 函数通过用户交互的方式设置阈值，来找到
    %   每个ROI信号中的峰值。
    %
    %   输入参数：
    %       traces_corrected - 矫正后的信号矩阵，每一列代表一个ROI的信号。
    %       t - 信号对应的时间轴，与traces_corrected的行数相同。
    %       colors - 用于绘制信号的颜色，每一行代表一个颜色。
    %       rois - ROI的数量。
    %
    %   返回值：
    %       peak_polarity - 每个ROI信号峰值的极性（1或-1）。
    %       peak_threshold - 用于找峰的阈值。
    %       peaks_index - 每个ROI信号峰值的索引。
    %       peaks_amplitude - 每个ROI信号峰值的幅度。
    %       peaks_sensitivity - 每个ROI信号峰值的灵敏度。
    %
    %   请在每个ROI的信号上用鼠标点击一次来设置查找峰值的阈值。
    %   函数将会根据设置的阈值和最小峰值显著性因子来找到信号中的峰值。    


    % 初始化变量
    nrois = size(traces_corrected,2);
    peak_polarity = zeros(1,nrois);
    peak_threshold = zeros(1,nrois);
    peaks_amplitude = cell(1, nrois);
    peaks_index = cell(1, nrois);
    peaks_sensitivity = cell(1, nrois);
    MinPeakProminence_factor = 0.4; % 定义最小峰值显著性因子
    figure; 
    set(gcf,'Position',get(0,'Screensize'));
    % 对每个ROI进行峰值检测
    for i = 1:nrois
        clf;

        % 确定峰值极性
        if abs(min(traces_corrected(:,i)) - mean(traces_corrected(:,i))) < ...
            max(abs(traces_corrected(:,i)) - mean(traces_corrected(:,i)))
            peak_polarity(i) = 1;
        else
            peak_polarity(i) = -1;
        end

        % plot trace
        
        plot_trace = traces_corrected(:,i) * peak_polarity(i);
        plot(plot_trace, 'Color', colors(i,:)); hold on;
        title(sprintf('ROI %d', i));

        % set threshold
        [~, peak_threshold(i)] = ginput(1);
        plot(ones(1,peak_threshold).*peak_threshold(i),'Color',colors(i,:),'LineWidth',2);
        hold off;
        

        % 寻找峰值
        MinPeakProminence = (max(plot_trace)-mean(plot_trace)) * MinPeakProminence_factor;
        [peak_y, peak_x] = findpeaks(plot_trace, 'MinPeakProminence', ...
            MinPeakProminence, 'MinPeakHeight', peak_threshold(i));
        peaks_index{i} = peak_x;
        current_trace = traces_corrected(:,i);
        peaks_amplitude{i} = current_trace(peak_x);
        peaks_sensitivity{i} = peak_y;
        pause(0.5);
    end
    close;
end
