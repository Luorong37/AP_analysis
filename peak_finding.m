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
    colors = lines(100);
    peak_polarity = zeros(1,nrois);
    peak_threshold = zeros(1,nrois);
    peaks_amplitude = cell(1, nrois);
    peaks_index = cell(1, nrois);
    peaks_sensitivity = cell(1, nrois);
    MinPeakProminence_factor = 0.4; % 定义最小峰值显著性因子
    fig = figure; 
    set(gcf,'Position',get(0,'Screensize'));

    i = 0;
    set(gcf, 'KeyPressFcn', @(src, event) set_space_pressed(event, fig));
    reselect = false;
    % 对每个ROI进行峰值检测
    while i < nrois
        if reselect
            i = i;
            reselect = false;
        else
            i = i + 1;
        end
        clf;
        current_trace = traces_corrected(:,i);

        % 确定峰值极性
        if abs(min(current_trace) - mean(current_trace)) < max(abs(current_trace) - mean(current_trace))
            peak_polarity(i) = 1;
            polarity = 'Positive';        
        else
            peak_polarity(i) = -1;
            polarity = 'Negative';
        end

        % plot trace
        plot_trace = current_trace * peak_polarity(i)+ 1 -peak_polarity(i) ;
        plot(plot_trace, 'Color', colors(i,:)); hold on;
        title(sprintf(['ROI %d, polarity = %s\n' ...
            'If no peaks, click a point above trace\n' ...
            'Press SPACE for continue, Press ENTER for end, Press R for reselect\n'], i , polarity));

        % set threshold
        [~, peak_threshold(i)] = ginput(1) ;
        plot(ones(1,size(traces_corrected,1)).*peak_threshold(i),'Color',colors(i,:),'LineWidth',2);
        hold on;
        

        % 寻找峰值
        MinPeakProminence = (max(plot_trace)-mean(plot_trace)) * MinPeakProminence_factor;
        [peak_y, peak_x] = findpeaks(plot_trace, 'MinPeakProminence', ...
            MinPeakProminence, 'MinPeakHeight', peak_threshold(i));
        
        peak_threshold(i) = (peak_threshold(i) + peak_polarity(i) -1)/peak_polarity(i);
        peaks_index{i} = peak_x;
        current_trace = current_trace;
        peaks_sensitivity{i} = current_trace(peak_x)-mean(current_trace);

        peaks_amplitude{i}  = peak_y;
        plot(peak_x, (peaks_amplitude{i}),'v','Color',colors(i,:),'MarkerFaceColor',colors(i,:));
        hold off;
        peaks_amplitude{i} = (peaks_amplitude{i} + peak_polarity(i) -1)/peak_polarity(i);
        

        fprintf('%d peaks found in ROI %d, Polarity : %s\n',numel(peak_x),i,polarity)

        % Wait for user input
        fig.UserData = [];
        waitfor(fig, 'UserData');
        if strcmp(fig.UserData, 'stop')
            break;
        elseif strcmp(fig.UserData, 'spacePressed')
            continue;
        elseif strcmp(fig.UserData,'reselect')
            reselect = true;
            continue
        end
    end
    close;
    fprintf('All peaks found.\n')
end

function set_space_pressed(event, fig)
if strcmp(event.Key, 'space')
    fig.UserData = 'spacePressed';
elseif strcmp(event.Key, 'return')
    fig.UserData = 'stop';
elseif strcmp(event.Key, 'r')
    fig.UserData = 'reselect';
end
end