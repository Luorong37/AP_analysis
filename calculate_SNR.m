function SNR_trace = calculate_SNR(trace)
    % Calculates Signal to Noise Ratio (SNR) using LOWESS smoothing
    % and RMSE calculation.
    %
    % Parameters:
    % trace: the input signal trace for a single ROI
    % threshold: the threshold value to exclude peaks in noise calculation

    % extract baseline
    trace_abs = abs(trace);
    baseline = trace_abs(trace_abs < prctile(trace_abs,90)); % defined extract 95%, means max frequence as 40 Hz
    
    % Calculate std from the baseline
    noise = std(baseline);

    % Calculate SNR 
    SNR_trace = (trace-mean(trace)) / noise;
end






