function [smoothed, info] = smooth_trace(trace, spec)
%SMOOTH_TRACE Smooth one explicitly selected trace stage.
%
% Frozen Dual calcium uses method="movmean", window=40, dimension=1 for
% raw, sensitivity, and snr in three separate calls.
% Inputs:
%   trace Numeric [frames x ROI] matrix.
%   spec Scalar struct requiring method, positive window, and dimension.
% Outputs:
%   smoothed Double matrix with the same dimensions as trace.
%   info Scalar struct containing the applied method and parameters.

if ~isstruct(spec) || ~all(isfield(spec, {'method', 'window', 'dimension'}))
    error('AAA:Algorithms:InvalidSmoothingSpec', ...
        'Smoothing spec requires method, window, and dimension.');
end

method = lower(string(spec.method));
switch method
    case "movmean"
        smoothed = movmean(double(trace), spec.window, spec.dimension);
    otherwise
        error('AAA:Algorithms:UnsupportedSmoothingMethod', ...
            'Unsupported smoothing method: %s.', spec.method);
end

info = struct( ...
    'method', char(method), ...
    'parameters', struct( ...
        'window', spec.window, ...
        'dimension', spec.dimension), ...
    'formula', 'smoothed = movmean(input, window, dimension)');
end
