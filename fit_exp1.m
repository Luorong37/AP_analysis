function [traces_corrected, fitted_curves] = fit_exp1(traces, freq)
% This function performs photobleaching correction on the intensity-time
% traces using a single exponential fit and optionally returns the fitted
% curves.
% Input arguments:
% traces: a matrix of size num_frames x num_rois containing the intensity-time traces for each ROI
% freq: the frame rate in Hz
% return_fit: a boolean to decide if the fitted curves should be returned


% Calculate time axis
time = (0:size(traces, 1) - 1) / freq;

% Initialize corrected traces matrix and fitted curves matrix if required
traces_corrected = zeros(size(traces));
fitted_curves = zeros(size(traces));

% Perform photobleaching correction for each ROI
for i = 1:size(traces, 2)
    % Extract the current trace
    current_trace = traces(:, i);
    
    % Fit the exponential decay model to the trace
    [fit_params,~] = fit(time', current_trace, 'exp1'); % 'exp1' is a model type for a single-term exponential
    
    % Generate the fitted curve from the fit parameters
    fitted_curve = fit_params.a * exp(fit_params.b * time);
    fitted_curves(:, i) = fitted_curve;
    
    % Correct the original trace by dividing by the fitted curve to mitigate photobleaching effect
    traces_corrected(:, i) = current_trace ./ fitted_curve';
end

end
