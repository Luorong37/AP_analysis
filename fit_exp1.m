function [traces_corrected, fitted_curves] = fit_exp1(traces, freq)
% ----------Write by Liu-Yang Luorong and ChatGPT----------
% ----------POWERED by Zoulab in Peking University----------
% Date: 23.11.16
% MATLAB Version: R2022b
% 
% FIT_EXP1 Applies exponential fitting for photobleaching correction in fluorescence traces.
%
%   This function corrects for photobleaching in fluorescence traces by fitting a single-term 
%   exponential decay model to each trace. It then divides each original trace by its fitted 
%   exponential curve to mitigate the photobleaching effects.
%
%   Syntax:
%   [traces_corrected, fitted_curves] = fit_exp1(traces, freq)
%
%   Parameters:
%   traces - A matrix containing fluorescence traces, with each column representing a trace.
%   freq - Sampling frequency of the traces.
%
%   Returns:
%   traces_corrected - A matrix of the same size as 'traces', containing the photobleaching-corrected traces.
%   fitted_curves - A matrix of fitted exponential curves corresponding to each trace in 'traces'.
%
%   Description:
%   - For each trace in 'traces', the function fits a single-term exponential decay model.
%   - The fitting is performed using MATLAB's 'fit' function with 'exp1' as the model type.
%   - Each original trace is then divided by its corresponding fitted curve to correct for photobleaching.
%   - The function returns both the corrected traces and the fitted exponential curves.
%
%   Example:
%   [corrected_traces, exp_curves] = fit_exp1(raw_traces, sampling_frequency);
%
%   Notes:
%   - Ensure that 'traces' are correctly preprocessed before applying this function.
%   - The function assumes that photobleaching follows a single-term exponential decay.
%
% See also FIT.


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
