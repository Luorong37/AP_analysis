function [is_valid, reason] = validate_peak_result(peak_results, expected_nrois, nframes)
%VALIDATE_PEAK_RESULT Validate one accepted_for_events result for reuse.

is_valid = false;
reason = "missing accepted_for_events peak results";
if ~isstruct(peak_results) || ~isscalar(peak_results) ...
        || ~isfield(peak_results, 'accepted_for_events')
    return;
end

accepted = peak_results.accepted_for_events;
required = {'index','amplitude','polarity'};
if ~isstruct(accepted) || ~isfield(accepted, 'data') ...
        || ~isstruct(accepted.data) || ~all(isfield(accepted.data, required))
    reason = "accepted_for_events data is incomplete";
    return;
end

indices = accepted.data.index;
amplitudes = accepted.data.amplitude;
polarities = accepted.data.polarity;
if ~iscell(indices) || ~iscell(amplitudes) || ~iscell(polarities) ...
        || numel(indices) ~= expected_nrois ...
        || numel(amplitudes) ~= expected_nrois ...
        || numel(polarities) ~= expected_nrois
    reason = "saved peak ROI count does not match the current sensitivity traces";
    return;
end

for roi_idx = 1:expected_nrois
    index = double(indices{roi_idx}(:));
    amplitude = double(amplitudes{roi_idx}(:));
    if numel(index) ~= numel(amplitude) || any(~isfinite(index)) ...
            || any(index ~= round(index)) || any(index < 1 | index > nframes)
        reason = sprintf( ...
            'saved peak indices/amplitudes are invalid for ROI %d', roi_idx);
        return;
    end
end

is_valid = true;
reason = "accepted peak ROI count and frame indices match current sensitivity traces";
end
