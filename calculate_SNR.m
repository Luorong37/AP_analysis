function [SNR_trace,fitted_trace] = calculate_SNR(trace, threshold, polarity)
    % Calculates Signal to Noise Ratio (SNR) using LOWESS smoothing
    % and RMSE calculation.
    %
    % Parameters:
    % trace: the input signal trace for a single ROI
    % threshold: the threshold value to exclude peaks in noise calculation

    
    % % Apply LOWESS smoothing to the trace
    % smooth_trace = smooth(trace, span, 'lowess');

    % Apply LOWESS smoothing to the trace
    [~, fitted_trace] = fit_exp1(trace);

    % Compute the residuals (difference between the trace and smoothed version)
    residuals = trace - fitted_trace;
    
    % Exclude the residuals that correspond to the peaks above the threshold
    % Only include data below threshold for noise estimation
    if nargin > 1
        if polarity == 1
            residuals = residuals(trace < threshold);
        elseif polarity == -1
            residuals = residuals(trace > threshold);
        end
    end
    
    % Calculate RMSE from the included residuals
    rmse_noise = sqrt(mean(residuals.^2));
    
    % Compute SNR using the peak amplitude divided by RMSE noise    
    SNR_trace = (trace - fitted_trace) / rmse_noise;
end
