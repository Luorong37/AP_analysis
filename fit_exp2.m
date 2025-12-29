function [traces_corrected, fitted_curves, params, gof] = fit_exp2(traces)
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
time = 1:size(traces, 1);

% Initialize corrected traces matrix and fitted curves matrix if required
traces_corrected = zeros(size(traces));
fitted_curves = zeros(size(traces));
params = cell(size(traces,2),1);
gofs = cell(size(traces,2),1);




% Perform photobleaching correction for each ROI
for i = 1:size(traces, 2)
    % Extract the current trace
    current_trace = traces(:, i);
    
    % --- 动态参数估算 ---
    max_val = max(current_trace);
    min_val = min(current_trace);
    mean_val = mean(current_trace);
    
    ft = fittype('exp2');
    opts = fitoptions(ft);
    
    % --- 自动选取 Upper 和 Lower ---
    % a, c (振幅): 最小为0，最大为信号最大值的2倍（给拟合留空间）
    % b, d (速率): 最小为-1（极速衰减），最大为0（不衰减）。强制为负是防止曲线向上翘。
    opts.Lower = [0, -1, 0, -1];  
    opts.Upper = [max_val*2, 0, max_val*2, 0]; 
    
    % --- 智能设置起始点 (StartPoint) ---
    % 一个快成分，一个慢成分通常效果最好
    opts.StartPoint = [mean_val/2, -1e-3, mean_val/2, -1e-5];
    % Fit the exponential decay model to the trace
    try
    [fit_params,gof] = fit(time', current_trace, ft, opts); 
    catch ME
        ft = fittype('poly1');
        opts = fitoptions(ft);
         [fit_params,gof] = fit(time', current_trace, ft, opts);
    end
    params{i} = fit_params;
    gofs{i} = gof;

    % Generate the fitted curve from the fit parameters
    switch type(ft)
        case 'exp2'
            fitted_curve = fit_params.a * exp(fit_params.b * time) + fit_params.c * exp(fit_params.d * time);
        case 'poly1'
            fitted_curve = fit_params.p1 *time + fit_params.p2;
    end
    fitted_curves(:, i) = fitted_curve;
    
    % Correct the original trace by dividing by the fitted curve to mitigate photobleaching effect
    traces_corrected(:, i) = current_trace ./ fitted_curve';
end

end
