function [SNR_trace,avg_baseline] = calculate_SNR(trace,peak_index,AP_window_width)
    % Calculates Signal to Noise Ratio (SNR) using LOWESS smoothing
    % and RMSE calculation.
    %
    % Parameters:
    % trace: the input signal trace for a single ROI
    % threshold: the threshold value to exclude peaks in noise calculation
    if nargin <2
    % extract baseline
        outlierPercentage = 0.05;
        % 计算离群值的数量
        numOutliers = ceil(length(trace) * outlierPercentage);
        
        % 使用标准差方法
        meanData = mean(trace);
        stdData = std(trace);
        
        % 计算每个值的z-score
        zScores = abs((trace - meanData) / stdData);
        
        % 找到z-score最高的numOutliers个值
        [~, sortedIndices] = sort(zScores, 'descend');
        inlierIndices = sortedIndices(numOutliers+1:end);
        
        % 提取离群值
        baseline = trace(inlierIndices);
        avg_baseline = mean(baseline);
        % Calculate std from the baseline
        noise = std(baseline);
    
        % Calculate SNR 
        SNR_trace = (trace-avg_baseline) / noise;
    else
        logicindex = true(size(trace));
        for i = 1:numel(peak_index)
            startindex = max(1, peak_index(i)-AP_window_width);
            endindex = min(length(trace), peak_index(i)+AP_window_width);
            logicindex(startindex:endindex) = false;
        end
        baseline = trace(logicindex);
        noise = std(baseline);
        avg_baseline = mean(baseline);
        % Calculate SNR 
        SNR_trace = (trace-avg_baseline) / noise;
    end
end






