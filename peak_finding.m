function [peak_polarity, peak_threshold, peaks_index, peaks_amplitude, peaks_sensitivity] = peak_finding(traces_corrected, t, colors, rois)
    % 初始化变量
    peak_polarity = zeros(1,length(rois));
    peak_threshold = zeros(1,length(rois));
    peaks_amplitude = cell(1, length(rois)-1);
    peaks_index = cell(1, length(rois)-1);
    peaks_sensitivity = cell(1, length(rois)-1);
    MinPeakProminence_factor = 0.4; % 定义最小峰值显著性因子
    figure; 
    set(gcf,'Position',get(0,'Screensize'));
    % 对每个ROI进行峰值检测
    for i = 1:length(rois)-1
        clf;

        % 确定峰值极性
        if abs(min(traces_corrected(:,i))-mean(traces_corrected(:,i))) < max(abs(traces_corrected(:,i)) - mean(traces_corrected(:,i)))
            peak_polarity(i) = 1;
        else
            peak_polarity(i) = -1;
        end

        % plot trace
        plot_trace = traces_corrected(:,i) * peak_polarity(i);
        plot(t, plot_trace, 'Color', colors(i,:)); hold on;
        title(sprintf('ROI %d', i));

        % set threshold
        [~, peak_threshold(i)] = ginput(1);
        plot(t,ones(1,length(t)).*peak_threshold(i),'Color',colors(i,:),'LineWidth',2);
        hold off;
        pause(0.2);

        % 寻找峰值
        MinPeakProminence = (max(plot_trace)-mean(plot_trace)) * MinPeakProminence_factor;
        [peak_y, peak_x] = findpeaks(plot_trace, 'MinPeakProminence', MinPeakProminence, 'MinPeakHeight', peak_threshold(i));
        peaks_index{i} = peak_x;
        current_trace = traces_corrected(:,i);
        peaks_amplitude{i} = current_trace(peak_x);
        peaks_sensitivity{i} = peak_y;
    end
    close;
end
