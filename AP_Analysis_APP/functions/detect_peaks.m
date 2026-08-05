function detected = detect_peaks(trace, profile)
%DETECT_PEAKS Detect peaks using one role-bound effective profile.
% The current shared Analysis3 peak contract supports voltage sensitivity
% traces. AP and Dual bind different public numeric defaults into
% profile.peak before calling this function.
% Inputs:
%   trace Numeric [frames x ROI] voltage sensitivity matrix.
%   profile Scalar bound voltage profile requiring positive frame_rate and
%           the documented profile.peak fields. It contains no request.
% Output:
%   detected Struct containing peak_table plus per-ROI index, amplitude,
%            and polarity cell arrays. Indices are one-based frame indices.

[role, params, frame_rate] = validate_inputs(trace, profile);
if role ~= "voltage"
    error('AAA:Algorithms:UnsupportedPeakRole', ...
        'Peak detection is not defined for role "%s".', role);
end

trace = double(trace);
nrois = size(trace, 2);
peak_table = empty_peak_table();
peaks_polarity = cell(1, nrois);
next_peak_id = 1;
min_distance = max(1, round(params.min_peak_distance_frames));
created_at = datetime("now");

for roi_idx = 1:nrois
    trace_i = trace(:, roi_idx);
    polarity = choose_polarity(trace_i, roi_idx, params);
    peaks_polarity{roi_idx} = polarity;
    plot_trace = trace_i * polarity;
    min_prominence = resolve_min_prominence(plot_trace, params);
    [~, peak_x] = findpeaks(plot_trace, ...
        'MinPeakProminence', min_prominence, ...
        'MinPeakDistance', min_distance, ...
        'MinPeakHeight', params.min_peak_height);
    if isempty(peak_x)
        continue;
    end

    peak_count = numel(peak_x);
    rows = table( ...
        (next_peak_id:next_peak_id + peak_count - 1)', ...
        repmat(roi_idx, peak_count, 1), ...
        peak_x(:), ...
        peak_x(:) ./ frame_rate, ...
        repmat(polarity, peak_count, 1), ...
        trace_i(peak_x(:)), ...
        repmat("accepted", peak_count, 1), ...
        repmat("auto", peak_count, 1), ...
        NaN(peak_count, 1), ...
        repmat("sensitivity_detected", peak_count, 1), ...
        repmat(created_at, peak_count, 1), ...
        'VariableNames', peak_table.Properties.VariableNames);
    peak_table = [peak_table; rows]; %#ok<AGROW>
    next_peak_id = next_peak_id + peak_count;
end

[peaks_index, peaks_amplitude, peaks_polarity] = ...
    peak_table_to_cells(peak_table, nrois, peaks_polarity);
detected = struct( ...
    'peak_table', peak_table, ...
    'index', {peaks_index}, ...
    'amplitude', {peaks_amplitude}, ...
    'polarity', {peaks_polarity});
end

function [role, params, frame_rate] = validate_inputs(trace, profile)
if ~isnumeric(trace) || isempty(trace) || ndims(trace) > 2
    error('AAA:Algorithms:InvalidPeakTrace', ...
        'Peak input must be a nonempty numeric frames-by-ROI matrix.');
end
if ~isstruct(profile) || ~isscalar(profile) ...
        || ~all(isfield(profile, {'role','frame_rate','peak'})) ...
        || ~isstruct(profile.peak) || ~isscalar(profile.peak)
    error('AAA:Algorithms:InvalidPeakProfile', ...
        'Peak profile requires bound role, frame_rate, and peak fields.');
end
role = lower(string(profile.role));
if ~isscalar(role) || strlength(role) == 0
    error('AAA:Algorithms:InvalidPeakProfile', ...
        'profile.role must be one nonempty value.');
end
frame_rate = double(profile.frame_rate);
if ~isscalar(frame_rate) || ~isfinite(frame_rate) || frame_rate <= 0
    error('AAA:Algorithms:InvalidPeakProfile', ...
        'profile.frame_rate must be positive and finite.');
end
params = profile.peak;
required = {'polarity_mode','global_polarity','roi_polarity', ...
    'min_peak_prominence','min_peak_prominence_mode', ...
    'min_peak_distance_frames','min_peak_height'};
missing = required(~isfield(params, required));
if ~isempty(missing)
    error('AAA:Algorithms:InvalidPeakProfile', ...
        'profile.peak is missing: %s.', strjoin(missing, ', '));
end
validateattributes(params.min_peak_prominence, {'numeric'}, ...
    {'scalar','real','finite','nonnegative'});
validateattributes(params.min_peak_distance_frames, {'numeric'}, ...
    {'scalar','real','finite','positive'});
validateattributes(params.min_peak_height, {'numeric'}, ...
    {'scalar','real','finite'});
end

function min_prominence = resolve_min_prominence(plot_trace, params)
mode = lower(string(params.min_peak_prominence_mode));
if mode == "absolute"
    min_prominence = params.min_peak_prominence;
elseif mode == "relative_factor"
    trace_scale = max(plot_trace, [], 'omitnan') ...
        - mean(plot_trace, 'omitnan');
    if ~isfinite(trace_scale) || trace_scale < 0
        trace_scale = 0;
    end
    min_prominence = params.min_peak_prominence * trace_scale;
else
    error('AAA:Algorithms:InvalidPeakProfile', ...
        'Unknown min_peak_prominence_mode: %s.', mode);
end
if ~isfinite(min_prominence) || min_prominence < 0
    min_prominence = 0;
end
end

function polarity = choose_polarity(trace_i, roi_idx, params)
[has_override, polarity] = roi_polarity_override( ...
    params.roi_polarity, roi_idx);
if has_override
    return;
end

mode = lower(string(params.polarity_mode));
switch mode
    case {"global", "voltage", "display"}
        polarity = normalize_polarity(params.global_polarity, ...
            'peak global polarity');
    case "positive"
        polarity = 1;
    case "negative"
        polarity = -1;
    case "auto"
        trace_mean = mean(trace_i, 'omitnan');
        if abs(max(trace_i, [], 'omitnan') - trace_mean) >= ...
                abs(min(trace_i, [], 'omitnan') - trace_mean)
            polarity = 1;
        else
            polarity = -1;
        end
    otherwise
        error('AAA:Algorithms:InvalidPeakProfile', ...
            'Unknown peak polarity_mode: %s.', mode);
end
end

function [has_override, polarity] = roi_polarity_override(values, roi_idx)
has_override = false;
polarity = NaN;
if isempty(values)
    return;
end
if iscell(values)
    if numel(values) < roi_idx || isempty(values{roi_idx})
        return;
    end
    value = values{roi_idx};
else
    if numel(values) < roi_idx || isempty(values(roi_idx))
        return;
    end
    value = values(roi_idx);
end
if isempty(value) || (isnumeric(value) && isscalar(value) ...
        && (isnan(value) || value == 0))
    return;
end
if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value)
    error('AAA:Algorithms:InvalidPeakProfile', ...
        'ROI %d peak polarity must be finite numeric when provided.', roi_idx);
end
polarity = normalize_polarity(value, sprintf('ROI %d peak polarity', roi_idx));
has_override = true;
end

function polarity = normalize_polarity(value, label)
if isempty(value) || ~isnumeric(value) || ~isscalar(value) ...
        || ~isfinite(value) || value == 0
    error('AAA:Algorithms:InvalidPeakProfile', ...
        '%s must be a finite nonzero numeric scalar.', label);
end
polarity = 1;
if value < 0
    polarity = -1;
end
end

function peak_table = empty_peak_table()
peak_table = table( ...
    zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
    strings(0,1), strings(0,1), zeros(0,1), strings(0,1), NaT(0,1), ...
    'VariableNames', {'peak_id', 'roi', 'index', 'time_s', 'polarity', ...
    'amplitude_sensitivity', 'status', 'source', 'parent_peak_id', ...
    'created_stage', 'created_at'});
end
