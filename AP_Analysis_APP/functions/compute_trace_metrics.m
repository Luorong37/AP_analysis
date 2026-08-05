function [metrics, info] = compute_trace_metrics(trace, baseline, profile)
%COMPUTE_TRACE_METRICS Compute frozen Dual noise, sensitivity, and SNR.
%
% profile is the complete bound channel profile. Role has one authority at
% profile.role; the effective model is profile.trace.noise_reference.
% The role/method pairs accepted from frozen Dual are:
%   voltage + wdenoise
%   calcium + butterworth_filtfilt
% No smoothing is performed here; use smooth_trace on each
% explicitly selected stage.
% Inputs:
%   trace Numeric [frames x ROI] bleach-removed trace matrix.
%   baseline Numeric [frames x ROI] matrix with exactly the same size.
%   profile Scalar bound profile requiring role and
%           trace.noise_reference. It contains no request.
% Outputs:
%   metrics Struct of [frames x ROI] noise_reference, noise, sensitivity,
%           and snr matrices.
%   info Scalar struct recording the effective role-specific formulas.

require_fields(profile, {'role', 'trace'}, 'channel profile');
require_fields(profile.trace, {'noise_reference'}, ...
    'channel profile.trace');
require_fields(profile.trace.noise_reference, {'method'}, ...
    'channel profile.trace.noise_reference');

role = lower(string(profile.role));
noise_profile = profile.trace.noise_reference;
method = lower(string(noise_profile.method));
trace = double(trace);
baseline = double(baseline);

switch role
    case "voltage"
        if method ~= "wdenoise"
            error('AAA:Algorithms:RoleMethodMismatch', ...
                'Voltage metric profile requires noise_reference.method="wdenoise".');
        end
        require_fields(noise_profile, ...
            {'level', 'denoising_method', 'wavelet_name'}, ...
            'voltage noise-reference profile');
        if exist('wdenoise', 'file') ~= 2
            error('Wavelet Toolbox function wdenoise is required for voltage noise estimation.');
        end
        p = noise_profile;
        noise_reference = wdenoise(trace, p.level, ...
            DenoisingMethod=char(string(p.denoising_method)), ...
            Wavelet=char(string(p.wavelet_name)));
        noise_reference_method = 'wdenoise';
        noise_reference_parameters = struct( ...
            'level', p.level, ...
            'denoising_method', char(string(p.denoising_method)), ...
            'wavelet_name', char(string(p.wavelet_name)));
    case "calcium"
        if method ~= "butterworth_filtfilt"
            error('AAA:Algorithms:RoleMethodMismatch', ...
                'Calcium metric profile requires noise_reference.method="butterworth_filtfilt".');
        end
        require_fields(noise_profile, ...
            {'order', 'normalized_cutoff'}, ...
            'calcium noise-reference profile');
        p = noise_profile;
        [b, a] = butter(p.order, p.normalized_cutoff);
        noise_reference = zeros(size(trace));
        for idx = 1:size(trace, 2)
            noise_reference(:, idx) = filtfilt(b, a, trace(:, idx));
        end
        noise_reference_method = 'butterworth_filtfilt';
        noise_reference_parameters = struct( ...
            'order', p.order, ...
            'normalized_cutoff', p.normalized_cutoff);
    otherwise
        error('AAA:Algorithms:UnsupportedMetricRole', ...
            'Unsupported metric profile role: %s.', profile.role);
end

noise = trace - noise_reference;
baseline_safe = baseline;
baseline_safe(abs(baseline_safe) < eps) = eps;
sensitivity = trace ./ baseline_safe;
noise_std = std(noise, 0, 1);
noise_std(noise_std < eps) = eps;
snr_value = trace ./ noise_std;

metrics = struct( ...
    'noise_reference', noise_reference, ...
    'noise', noise, ...
    'sensitivity', sensitivity, ...
    'snr', snr_value);
info = struct( ...
    'role', role, ...
    'profile', profile, ...
    'noise_reference_method', noise_reference_method, ...
    'noise_reference_parameters', noise_reference_parameters, ...
    'snr_method', 'signal_divided_by_noise_std', ...
    'snr_parameters', struct( ...
        'signal_stage', 'bleach_removed', ...
        'noise_stage', 'noise'), ...
    'formulas', struct( ...
        'noise', 'noise = trace - noise_reference', ...
        'sensitivity', 'sensitivity = trace ./ baseline_safe', ...
        'snr', 'snr = trace ./ std(noise, 0, 1)', ...
        'zero_clamp', 'abs(baseline) < eps -> eps; std(noise) < eps -> eps'));
end

function require_fields(value, names, label)
if ~isstruct(value)
    error('AAA:Algorithms:InvalidMetricProfile', '%s must be a struct.', label);
end
missing = names(~isfield(value, names));
if ~isempty(missing)
    error('AAA:Algorithms:InvalidMetricProfile', ...
        '%s requires: %s.', label, strjoin(missing, ', '));
end
end
