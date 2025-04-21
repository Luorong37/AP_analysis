function [peaks_polarity, peaks_threshold, peaks_index, peaks_amplitude, peaks_sensitivity] = ...
                                            peak_finding(traces,MinPeakProminence_factor)
    % PEAK_FINDING 寻找给定ROI信号中的峰值
    %
    %   此函数通过用户交互方式，帮助检测和识别给定信号中的峰值。用户可以在每个ROI信号上手动
    %   设置阈值，以确定峰值的检测条件。该函数返回每个ROI信号的峰值极性、阈值、峰值索引、峰值
    %   幅度及峰值的灵敏度。
    %
    %   输入参数：
    %       traces - 一个矩阵，其中每列表示一个ROI信号，每行是不同时间点的信号值。
    %       MinPeakProminence_factor - （可选）最小峰值显著性因子，默认为0.5。
    %
    %   输出参数：
    %       peaks_polarity - 每个ROI信号峰值的极性，1代表正峰，-1代表负峰。
    %       peaks_threshold - 检测峰值时使用的阈值。
    %       peaks_index - 每个ROI信号中峰值的索引位置（时间点）。
    %       peaks_amplitude - 每个ROI信号中峰值的幅度。
    %       peaks_sensitivity - 每个ROI信号中峰值的灵敏度。
    %
    %   使用说明：
    %   用户将通过鼠标点击每个ROI信号图来设置阈值，函数根据所设置的阈值和最小峰值显著性因子，
    %   检测信号中的峰值。按空格键继续到下一个ROI信号，按回车键结束，按R键重新选择当前信号的阈值。
    %
    %   示例：
    %   [polarity, threshold, index, amplitude, sensitivity] = peak_finding(traces_corrected);
    %
    %   Notes:
    %   - 此函数采用交互式方式处理数据，请确保图形窗口处于活动状态并能接受用户输入。
    %   - 函数的交互步骤包括用户通过点击来设置阈值和按键控制继续或重选。
    %   - update at 20240910

    % 初始化变量
    if nargin < 2
    MinPeakProminence_factor = 0.5; % 定义最小峰值显著性因子
    end

    nrois = size(traces,2);
    colors = lines(100);
    peaks_polarity = zeros(1,nrois);
    peaks_threshold = zeros(1,nrois);
    peaks_amplitude = cell(1, nrois);
    peaks_index = cell(1, nrois);
    peaks_sensitivity = cell(1, nrois);


    fig = figure; 
    set(gcf,'Position',get(0,'Screensize'));

    i = 0;
    set(gcf, 'KeyPressFcn', @(src, event) set_space_pressed(event, fig));
    reselect = false;
    flipped = false;
    % 对每个ROI进行峰值检测
    while i <= nrois
        if reselect
            i = i;
            reselect = false;
        elseif i ~= nrois
            i = i + 1;
        end
        clf;
        current_trace = traces(:,i);

        % 确定峰值极性
        % abs(min(current_trace) - mean(current_trace)) < max(abs(current_trace) - mean(current_trace))
        if abs(min(current_trace) - mean(current_trace)) < abs(max(current_trace) - mean(current_trace))
            peaks_polarity(i) = 1;
            polarity = 'Positive';        
        else
            peaks_polarity(i) = -1;
            polarity = 'Negative';
        end
        
        % manually flipped polarity
        if flipped
            if strcmp(polarity,'Positive')
                polarity = 'Negative';
                peaks_polarity(i) = -1;
            elseif strcmp(polarity,'Negative')
                polarity = 'Positive';   
                peaks_polarity(i) = 1;
            end
            flipped = false;
        end
        % plot trace
        % plot_trace = current_trace * peaks_polarity(i)+ 1 -peaks_polarity(i) ;
        plot_trace = current_trace * peaks_polarity(i) ;
        plot(plot_trace); hold on;
        title(sprintf(['ROI %d, polarity = %s\n' ...
            'If no peaks, click a point above trace\n' ...
            'Press SPACE for continue, Press ENTER for end, Press R for reselect, Press F for manually flip polarity\n'], i , polarity));

        % set threshold
        [~, peaks_threshold(i)] = ginput(1) ;
        plot(ones(1,size(traces,1)).*peaks_threshold(i),'LineWidth',2);
        hold on;
        

        % 寻找峰值
        MinPeakProminence = (max(plot_trace)-mean(plot_trace)) * MinPeakProminence_factor;
        [peak_y, peak_x] = findpeaks(plot_trace, 'MinPeakProminence', ...
            MinPeakProminence, 'MinPeakHeight', peaks_threshold(i));
        
        % peaks_threshold(i) = (peaks_threshold(i) + peaks_polarity(i) -1)/peaks_polarity(i);
        peaks_threshold(i) = peaks_threshold(i)*peaks_polarity(i);
        peaks_index{i} = peak_x;
        % current_trace = current_trace;
        peaks_sensitivity{i} = current_trace(peak_x)-mean(current_trace);

        peaks_amplitude{i}  = peak_y;
        % plot peak
        plot(peak_x, (peaks_amplitude{i}),'v','MarkerFaceColor',colors(i,:));
        hold off;
        % peaks_amplitude{i} = (peaks_amplitude{i} + peaks_polarity(i) -1)/peaks_polarity(i);
        peaks_amplitude{i} = peaks_amplitude{i} * peaks_polarity(i);

        if numel(peak_x) == 0
            peaks_polarity(i) = 0;
        end

        fprintf('%d peaks found in ROI %d, Polarity : %s\n',numel(peak_x),i,polarity)

        % Wait for user input
        fig.UserData = [];
        waitfor(fig, 'UserData');
        switch  fig.UserData
            case  'stop'
            break;
            case 'spacePressed'
            continue;
            case {'reselect','flipped'}
            reselect = true;
            if strcmp(fig.UserData,'flipped')
                flipped = true;
            end
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
elseif strcmp(event.Key, 'f') 
    fig.UserData = 'flipped';
end
end