function edited = edit_peaks(trace, detected, profile)
%EDIT_PEAKS Apply the shared interactive add/delete peak editor.
% The editor changes only the peak table. It does not redirect peaks,
% extract event windows, gate events, or calculate event measurements.
% Inputs:
%   trace Numeric [frames x ROI] voltage sensitivity matrix.
%   detected Struct returned by detect_peaks for the same trace dimensions.
%   profile Scalar bound voltage profile requiring frame_rate and peak
%           fields, including run_manual_edit. It contains no request.
% Output:
%   edited Struct with the accepted peak table, per-ROI cells, edit history,
%          and ROI polarities. When manual editing is disabled it is a
%          noninteractive normalization of detected.

[trace, peak_table, peaks_polarity, params] = ...
    validate_inputs(trace, detected, profile);
nrois = size(trace, 2);
history = empty_history();

if ~params.run_manual_edit
    edited = finish_result(peak_table, history, peaks_polarity, nrois);
    return;
end
if ~usejava('desktop')
    error('AAA:Algorithms:PeakEditorRequiresDesktop', ...
        ['Manual peak editing requires a MATLAB Desktop session. ' ...
         'Disable the public manual-edit control for headless execution.']);
end

next_peak_id = 1;
if height(peak_table) > 0
    next_peak_id = max(peak_table.peak_id) + 1;
end
quit_editor = false;
roi_idx = 1;
while roi_idx <= nrois
    if quit_editor
        break;
    end
    roi_step = 1;
    mode = "delete_box";
    reset_table = peak_table;
    reset_polarity = peaks_polarity{roi_idx};
    fig = figure('Name', ...
        sprintf('ROI %d voltage sensitivity peak editor', roi_idx));
    set(fig, 'Position', get(0, 'Screensize'));
    while ishandle(fig)
        ax = draw_editor(fig, trace(:, roi_idx), peak_table, ...
            peaks_polarity{roi_idx}, roi_idx, mode);
        key = wait_for_key(fig);
        if ~ishandle(fig)
            break;
        end
        switch key
            case "a"
                mode = "add_box";
                [box_position, box_ok] = get_editor_box(ax);
                if box_ok
                    [peak_table, history, next_peak_id] = add_in_box( ...
                        trace(:, roi_idx), peak_table, history, ...
                        next_peak_id, roi_idx, peaks_polarity{roi_idx}, ...
                        profile.frame_rate, params, box_position, mode);
                end
            case "d"
                mode = "delete_box";
                [box_position, box_ok] = get_editor_box(ax);
                if box_ok
                    [peak_table, history] = delete_in_box( ...
                        trace(:, roi_idx), peak_table, history, roi_idx, ...
                        peaks_polarity{roi_idx}, box_position, mode);
                end
            case "c"
                mode = "clear_current_roi";
                [peak_table, history] = clear_roi( ...
                    peak_table, history, roi_idx, mode);
            case "f"
                mode = "flip_polarity";
                old_polarity = normalize_polarity(peaks_polarity{roi_idx});
                new_polarity = -old_polarity;
                peaks_polarity{roi_idx} = new_polarity;
                peak_table.polarity(peak_table.roi == roi_idx) = new_polarity;
                history(end + 1) = history_entry( ...
                    "flip_polarity", roi_idx, NaN, ...
                    old_polarity, new_polarity, mode); %#ok<AGROW>
            case "r"
                peak_table(peak_table.roi == roi_idx, :) = [];
                peak_table = [peak_table; ...
                    reset_table(reset_table.roi == roi_idx, :)];
                peaks_polarity{roi_idx} = reset_polarity;
                history(end + 1) = history_entry( ...
                    "reset", roi_idx, NaN, NaN, NaN, mode); %#ok<AGROW>
            case "q"
                quit_editor = true;
                close(fig);
            case {"n", "space", "return"}
                close(fig);
            case {"p", "leftarrow", "backspace"}
                roi_step = -1;
                close(fig);
        end
    end
    if ~quit_editor
        roi_idx = max(1, roi_idx + roi_step);
    end
end

edited = finish_result(peak_table, history, peaks_polarity, nrois);
end

function [trace, peak_table, polarities, params] = ...
        validate_inputs(trace, detected, profile)
if ~isnumeric(trace) || isempty(trace) || ndims(trace) > 2
    error('AAA:Algorithms:InvalidPeakTrace', ...
        'Peak editor input must be a nonempty numeric matrix.');
end
trace = double(trace);
if ~isstruct(detected) || ~isscalar(detected) ...
        || ~all(isfield(detected, {'peak_table','polarity'})) ...
        || ~istable(detected.peak_table) || ~iscell(detected.polarity)
    error('AAA:Algorithms:InvalidDetectedPeaks', ...
        'Detected peaks must contain peak_table and polarity.');
end
peak_table = detected.peak_table;
polarities = detected.polarity;
if numel(polarities) ~= size(trace, 2)
    error('AAA:Algorithms:PeakRoiCountMismatch', ...
        'Detected peak polarity count does not match trace ROI count.');
end
if ~isstruct(profile) || ~isscalar(profile) ...
        || ~all(isfield(profile, {'role','frame_rate','peak'})) ...
        || lower(string(profile.role)) ~= "voltage" ...
        || ~isstruct(profile.peak) || ~isscalar(profile.peak) ...
        || ~all(isfield(profile.peak, ...
            {'run_manual_edit','min_peak_distance_frames'}))
    error('AAA:Algorithms:InvalidPeakProfile', ...
        'Peak editor requires one bound voltage profile.');
end
params = profile.peak;
end

function edited = finish_result(table_value, history, polarities, nrois)
[index, amplitude, polarity] = peak_table_to_cells( ...
    table_value, nrois, polarities);
edited = struct( ...
    'peak_table', table_value, ...
    'edit_history', history, ...
    'index', {index}, ...
    'amplitude', {amplitude}, ...
    'polarity', {polarity});
end

function ax = draw_editor(fig, trace, table_value, polarity, roi_idx, mode)
polarity = normalize_polarity(polarity);
figure(fig);
clf(fig);
set(fig, 'KeyPressFcn', @(~, event) set_editor_key(event, fig));
ax = axes('Parent', fig, 'Units', 'normalized', ...
    'Position', [0.06 0.08 0.91 0.86]);
plot(ax, trace * polarity, 'k');
hold(ax, 'on');
accepted = table_value.roi == roi_idx & table_value.status == "accepted";
deleted = table_value.roi == roi_idx & table_value.status == "deleted";
accepted_index = table_value.index(accepted);
deleted_index = table_value.index(deleted);
if ~isempty(accepted_index)
    plot(ax, accepted_index, trace(accepted_index) * polarity, ...
        'rv', 'MarkerFaceColor', 'r');
end
if ~isempty(deleted_index)
    plot(ax, deleted_index, trace(deleted_index) * polarity, ...
        'x', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
end
title(ax, sprintf([ ...
    'ROI %d | polarity=%+d | %s | A add, D delete, C clear, ' ...
    'F flip, N/Space next, P/Left previous, R reset, Q quit'], ...
    roi_idx, polarity, upper(mode)));
xlabel(ax, 'Frame');
ylabel(ax, 'Sensitivity x polarity');
grid(ax, 'on');
end

function set_editor_key(event, fig)
allowed = {'a','d','c','f','n','space','return','p', ...
    'leftarrow','backspace','r','q'};
if any(strcmp(event.Key, allowed))
    fig.UserData.key = event.Key;
end
end

function key = wait_for_key(fig)
fig.UserData.key = [];
waitfor(fig, 'UserData');
if ishandle(fig) && isstruct(fig.UserData) && isfield(fig.UserData, 'key')
    key = string(fig.UserData.key);
else
    key = "";
end
end

function [position, ok] = get_editor_box(ax)
position = [NaN NaN NaN NaN];
ok = false;
try
    rectangle_value = drawrectangle(ax);
    wait(rectangle_value);
    if isvalid(rectangle_value)
        position = rectangle_value.Position;
        delete(rectangle_value);
        ok = all(isfinite(position)) && position(3) > 0 && position(4) > 0;
    end
catch ME
    warning('AAA:Algorithms:PeakEditRectangleFailed', ...
        'Rectangle selection failed: %s', ME.message);
end
end

function [table_value, history, next_id] = add_in_box( ...
        trace, table_value, history, next_id, roi_idx, polarity, ...
        frame_rate, params, position, mode)
plot_trace = trace * polarity;
frame_range = max(1, ceil(position(1))): ...
    min(numel(trace), floor(position(1) + position(3)));
if isempty(frame_range)
    return;
end
[peak_y, relative_index] = findpeaks(plot_trace(frame_range), ...
    'MinPeakDistance', max(1, round(params.min_peak_distance_frames)));
candidate_index = frame_range(relative_index);
inside_y = peak_y >= position(2) ...
    & peak_y <= position(2) + position(4);
candidate_index = unique(candidate_index(inside_y));
existing = table_value.index(table_value.roi == roi_idx ...
    & table_value.status == "accepted");

for index = reshape(candidate_index, 1, [])
    if ~isempty(existing) ...
            && any(abs(existing - index) < params.min_peak_distance_frames)
        continue;
    end
    row = table(next_id, roi_idx, index, index / frame_rate, polarity, ...
        trace(index), "accepted", "manual_add", NaN, ...
        "sensitivity_manually_edited", datetime("now"), ...
        'VariableNames', table_value.Properties.VariableNames);
    table_value = [table_value; row]; %#ok<AGROW>
    history(end + 1) = history_entry( ...
        "add_box", roi_idx, next_id, NaN, index, mode); %#ok<AGROW>
    existing(end + 1, 1) = index; %#ok<AGROW>
    next_id = next_id + 1;
end
end

function [table_value, history] = delete_in_box( ...
        trace, table_value, history, roi_idx, polarity, position, mode)
rows = find(table_value.roi == roi_idx ...
    & table_value.status == "accepted");
if isempty(rows)
    return;
end
peak_x = table_value.index(rows);
peak_y = trace(peak_x) * polarity;
inside = peak_x >= position(1) ...
    & peak_x <= position(1) + position(3) ...
    & peak_y >= position(2) ...
    & peak_y <= position(2) + position(4);
for row_idx = reshape(rows(inside), 1, [])
    table_value.status(row_idx) = "deleted";
    history(end + 1) = history_entry( ...
        "delete_box", roi_idx, table_value.peak_id(row_idx), ...
        table_value.index(row_idx), NaN, mode); %#ok<AGROW>
end
end

function [table_value, history] = clear_roi( ...
        table_value, history, roi_idx, mode)
rows = find(table_value.roi == roi_idx ...
    & table_value.status == "accepted");
for row_idx = reshape(rows, 1, [])
    table_value.status(row_idx) = "deleted";
    history(end + 1) = history_entry( ...
        "clear_current_roi", roi_idx, table_value.peak_id(row_idx), ...
        table_value.index(row_idx), NaN, mode); %#ok<AGROW>
end
end

function entry = history_entry(action, roi, peak_id, before, after, mode)
entry = struct( ...
    'action', string(action), ...
    'roi', roi, ...
    'peak_id', peak_id, ...
    'index_before', before, ...
    'index_after', after, ...
    'mode', string(mode), ...
    'created_at', datetime("now"));
end

function history = empty_history()
history = struct( ...
    'action', {}, ...
    'roi', {}, ...
    'peak_id', {}, ...
    'index_before', {}, ...
    'index_after', {}, ...
    'mode', {}, ...
    'created_at', {});
end

function polarity = normalize_polarity(value)
if isempty(value) || ~isnumeric(value) || ~isscalar(value) ...
        || ~isfinite(value) || value == 0
    error('AAA:Algorithms:InvalidPeakPolarity', ...
        'Peak polarity must be a finite nonzero numeric scalar.');
end
polarity = 1;
if value < 0
    polarity = -1;
end
end
