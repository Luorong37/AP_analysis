function [SNR_trace,smooth_trace] = calculate_SNR(trace, threshold, polarity)
    % Calculates Signal to Noise Ratio (SNR) using LOWESS smoothing
    % and RMSE calculation.
    %
    % Parameters:
    % trace: the input signal trace for a single ROI
    % threshold: the threshold value to exclude peaks in noise calculation

    
    % % Apply LOWESS smoothing to the trace
    % smooth_trace = smooth(trace, span, 'lowess');

    % Apply LOWESS smoothing to the trace
    [~, smooth_trace] = fit_exp1(trace);

    % Compute the residuals (difference between the trace and smoothed version)
    residuals = trace - smooth_trace;
    
    % Exclude the residuals that correspond to the peaks above the threshold
    % Only include data below threshold for noise estimation
    if polarity == 1
        residuals = residuals(trace < threshold);
    elseif polarity == -1
        residuals = residuals(trace > threshold);
    end
    
    % Calculate RMSE from the included residuals
    rmse_noise = sqrt(mean(residuals.^2));
    
    % Compute SNR using the peak amplitude divided by RMSE noise    
    SNR_trace = (trace - smooth_trace) / rmse_noise;
end
