function traces_corrected = fit_exp1(traces, freq)
% This function performs photobleaching correction on the intensity-time
% traces using a single exponential fit.
% Input arguments:
% traces: a matrix of size num_frames x num_rois containing the intensity-time traces for each ROI
% Hz: the frame rate in Hz

% Calculate time axis
time = (0:size(traces, 1) - 1) / freq;

% Initialize corrected traces matrix
traces_corrected = zeros(size(traces));

% Perform photobleaching correction for each ROI
for i = 1:size(traces, 2)
    % Extract the current trace
    current_trace = traces(:, i);
    
    % Fit the exponential decay model to the trace
    fit_params = fit(time', current_trace, 'exp1'); % 'exp1' is a model type for a single-term exponential
    
    % Calculate the fitted exponential decay curve
    fitted_curve = fit_params.a * exp(fit_params.b * time);
    
    % Correct the original trace by subtracting the fitted curve
    % and adding the offset to maintain non-negative intensity values
    offset = min(fitted_curve);
    traces_corrected(:, i) = current_trace - (fitted_curve - offset);
end

end
