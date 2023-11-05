function [SNR_trace,smooth_trace] = calculate_SNR(trace, threshold, polarity, span)
    % Calculates Signal to Noise Ratio (SNR) using LOWESS smoothing
    % and RMSE calculation.
    %
    % Parameters:
    % trace: the input signal trace for a single ROI
    % threshold: the threshold value to exclude peaks in noise calculation
    % span: the fraction of data used in each local regression (0< span <=1)
    
    % Check inputs
    if nargin < 4
        span = 0.1; % Default span if not provided
    end
    
    % Apply LOWESS smoothing to the trace
    smooth_trace = smooth(trace, span, 'lowess');
    
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
    rmse_noise = sqrt(mean(residuals.^2))
    
    % Compute SNR using the peak amplitude divided by RMSE noise    
    SNR_trace = (trace - smooth_trace) / rmse_noise;
end
