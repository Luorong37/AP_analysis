function [traces_bleach_removed, baseline, parameters] = remove_bleach(traces_input, spec)
%REMOVE_BLEACH Apply one frozen Dual bleaching-removal model.
%
% Required spec fields:
%   mode = "linear" | "highpass" | "exp2"
% For highpass, spec.frame_rate and spec.time_axis are also required. The
% cutoff remains fc = 0.5 / time_axis(end), exactly as in Dual_analysis3.
% Inputs:
%   traces_input Numeric [frames x ROI] trace matrix.
%   spec Scalar algorithm struct. frame_rate is in Hz and time_axis is a
%        [frames x 1] seconds vector when mode is highpass.
% Outputs:
%   traces_bleach_removed [frames x ROI], same size as input.
%   baseline [frames x ROI], same size as input.
%   parameters Scalar struct containing the effective fitted/cutoff values.

if ~isstruct(spec) || ~isfield(spec, 'mode')
    error('AAA:Algorithms:InvalidBleachSpec', ...
        'Bleach spec.mode is required.');
end

bleach_mode = lower(string(spec.mode));
switch bleach_mode
    case "highpass"
        require_fields(spec, {'frame_rate', 'time_axis'}, 'highpass bleach spec');
        freq = spec.frame_rate;
        t_axis = spec.time_axis;
        fc = 0.5 / t_axis(end);
        [traces_bleach_removed, baseline] = ...
            highpass_bleach_remove_frozen(traces_input, freq, fc);
        parameters = struct('fc', fc);
    case "linear"
        traces_bleach_removed = detrend(traces_input, 1);
        baseline = traces_input - traces_bleach_removed;
        parameters = struct();
    case "exp2"
        [~, baseline] = fit_exp2_frozen(traces_input);
        traces_bleach_removed = traces_input - baseline;
        parameters = struct();
    otherwise
        error('Unsupported bleach_mode: %s', spec.mode);
end
end

function require_fields(value, names, label)
missing = names(~isfield(value, names));
if ~isempty(missing)
    error('AAA:Algorithms:InvalidBleachSpec', ...
        '%s requires: %s.', label, strjoin(missing, ', '));
end
end

function [traces_corrected, baseline] = highpass_bleach_remove_frozen(traces, Fs, FcL)
% Exact production body used by the frozen Dual highpass helper, with all
% formerly optional inputs supplied explicitly by remove_bleach.
padlength = round(0.05 * size(traces, 1));
padded_traces = [repmat(traces(1, :), padlength, 1); ...
                 traces; ...
                 repmat(traces(end, :), padlength, 1)];
filtered_traces = zeros(size(padded_traces));
for col = 1:size(traces, 2)
    filtered_traces(:, col) = highpass(padded_traces(:, col), FcL, Fs);
end
traces_corrected = filtered_traces(padlength+1:end-padlength, :);
baseline = smoothdata(traces - traces_corrected);
end

function [traces_corrected, fitted_curves, params, gof] = fit_exp2_frozen(traces)
% Exact numerical fit/fallback sequence used by the frozen fit_exp2 helper.
time = 1:size(traces, 1);
traces_corrected = zeros(size(traces));
fitted_curves = zeros(size(traces));
params = cell(size(traces, 2), 1);
gofs = cell(size(traces, 2), 1); %#ok<NASGU>

for idx = 1:size(traces, 2)
    current_trace = traces(:, idx);
    max_val = max(current_trace);
    mean_val = mean(current_trace);

    ft = fittype('exp2');
    opts = fitoptions(ft);
    opts.Lower = [0, -1, 0, -1];
    opts.Upper = [max_val*2, 0, max_val*2, 0];
    opts.StartPoint = [mean_val/2, -1e-3, mean_val/2, -1e-5];
    try
        [fit_params, gof] = fit(time', current_trace, ft, opts);
    catch
        ft = fittype('poly1');
        opts = fitoptions(ft);
        [fit_params, gof] = fit(time', current_trace, ft, opts);
    end
    params{idx} = fit_params;
    gofs{idx} = gof;

    switch type(ft)
        case 'exp2'
            fitted_curve = fit_params.a * exp(fit_params.b * time) ...
                + fit_params.c * exp(fit_params.d * time);
        case 'poly1'
            fitted_curve = fit_params.p1 * time + fit_params.p2;
    end
    fitted_curves(:, idx) = fitted_curve;
    traces_corrected(:, idx) = current_trace ./ fitted_curve';
end
end
