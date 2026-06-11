% check_dual_frame_alignment_from_logs
%
% Check dual-camera frame alignment from multiple acquisition-stage sources:
%   1. NI-DAQ edge-count channels that received camera return/frame signals;
%   2. logs.sync camera frame counters used by the rebuilt visualstim crop;
%   3. saved camera movie files under each Cycle*/Cam* folder.
%
% Rebuilt reference:
%   logs.daq       = read(d.in, "all")
%   Camera counters are logs.daq columns 2..N+1 for N cameras.
%   logs.sync.Camera*_Frame stores those camera counter values sampled at
%   PTB/DAQ stimulus sync events.
%
% Usage:
%   root_path = 'E:\path\to\Methods_or_Rec_folder';
%   run('check_dual_frame_alignment_from_logs.m');
%
% Optional:
%   save_outputs = true;
%   camera_count = 2;
%   count_saved_movies = true;
%   make_plots = true;
%   plot_only_mismatches = true;
%   frame_match_tolerance_s = []; % empty = estimate from DAQ event timing

if ~exist('root_path', 'var') || isempty(root_path)
    root_path = pwd;
end
if ~exist('save_outputs', 'var') || isempty(save_outputs)
    save_outputs = true;
end
if ~exist('camera_count', 'var') || isempty(camera_count)
    camera_count = 2;
end
if ~exist('count_saved_movies', 'var') || isempty(count_saved_movies)
    count_saved_movies = true;
end
if ~exist('make_plots', 'var') || isempty(make_plots)
    make_plots = true;
end
if ~exist('plot_only_mismatches', 'var') || isempty(plot_only_mismatches)
    plot_only_mismatches = true;
end
if ~exist('frame_match_tolerance_s', 'var') || isempty(frame_match_tolerance_s)
    frame_match_tolerance_s = [];
end
if ~exist('plot_output_dir', 'var') || isempty(plot_output_dir)
    plot_output_dir = fullfile(root_path, 'Dual_analysis3_DAQ_frame_alignment_figures');
end

if ~isfolder(root_path)
    error('root_path does not exist or is not a folder: %s', root_path);
end

fprintf('\n============================================================\n');
fprintf('Check Dual Frame Alignment From Logs And Saved Movies\n');
fprintf('============================================================\n');
fprintf('Root: %s\n', root_path);
fprintf('Camera count: %d\n', camera_count);
fprintf('Count saved movie files: %d\n', logical(count_saved_movies));
fprintf('Make plots: %d\n', logical(make_plots));

log_files = find_files_recursive(root_path, 'logs.mat');
fprintf('logs.mat files found: %d\n', numel(log_files));
if isempty(log_files)
    error('No logs.mat files were found under %s.', root_path);
end

rows = repmat(make_empty_summary_row(camera_count), 0, 1);
for file_idx = 1:numel(log_files)
    log_file = char(log_files(file_idx));
    try
        rows(end+1, 1) = summarize_one_logs_file(log_file, camera_count, count_saved_movies);
    catch ME
        failed_row = make_empty_summary_row(camera_count);
        failed_row.log_file = string(log_file);
        failed_row.cycle_path = string(fileparts(log_file));
        [failed_row.rec_name, failed_row.cycle_name] = infer_rec_cycle_from_path(fileparts(log_file));
        failed_row.status = "failed_to_read_log";
        failed_row.message = string(ME.message);
        rows(end+1, 1) = failed_row;
        fprintf(2, 'Failed to read %s: %s\n', log_file, ME.message);
    end
end

daq_frame_alignment_table = struct2table(rows);
daq_frame_alignment_table = sortrows(daq_frame_alignment_table, {'rec_name', 'cycle_number', 'cycle_name'});
mismatch_table = daq_frame_alignment_table( ...
    daq_frame_alignment_table.abs_daq_return_diff > 0 | ...
    daq_frame_alignment_table.abs_sync_predicted_crop_diff > 0 | ...
    daq_frame_alignment_table.abs_saved_movie_frame_diff > 0 | ...
    daq_frame_alignment_table.status ~= "aligned", :);
layer_mismatch_rows = build_layer_mismatch_rows(daq_frame_alignment_table);
if isempty(layer_mismatch_rows)
    layer_mismatch_table = struct2table(make_empty_layer_mismatch_row());
    layer_mismatch_table(1, :) = [];
else
    layer_mismatch_table = struct2table(layer_mismatch_rows);
end

missing_event_rows = repmat(make_empty_missing_event_row(), 0, 1);
if make_plots && ~isfolder(plot_output_dir)
    mkdir(plot_output_dir);
end
for row_idx = 1:height(daq_frame_alignment_table)
    current_row = daq_frame_alignment_table(row_idx, :);
    if current_row.status == "failed_to_read_log" || current_row.status == "missing_logs_daq"
        continue;
    end
    try
        S = load(char(current_row.log_file));
        event_report = analyze_camera_return_events_from_logs(S.logs, camera_count, frame_match_tolerance_s);
        cycle_missing_rows = build_missing_event_rows(current_row, event_report);
        missing_event_rows = [missing_event_rows; cycle_missing_rows(:)]; %#ok<AGROW>

        current_layer_mismatch = has_any_layer_mismatch(current_row);
        should_plot = logical(make_plots) && (~logical(plot_only_mismatches) || ~isempty(cycle_missing_rows) || current_layer_mismatch);
        if should_plot
            plot_file_base = fullfile(plot_output_dir, sanitize_file_component(current_row.rec_name + "_" + current_row.cycle_name + "_DAQ_frame_alignment"));
            plot_daq_frame_alignment_event_report(current_row, event_report, plot_file_base);
        end
    catch ME
        warning('Failed to analyze DAQ return events for %s/%s: %s', ...
            string(current_row.rec_name), string(current_row.cycle_name), ME.message);
    end
end

if isempty(missing_event_rows)
    missing_return_event_table = struct2table(make_empty_missing_event_row());
    missing_return_event_table(1, :) = [];
else
    missing_return_event_table = struct2table(missing_event_rows);
end

fprintf('\nCycles parsed: %d\n', height(daq_frame_alignment_table));
fprintf('Mismatched / incomplete cycles: %d\n', height(mismatch_table));
fprintf('Layer mismatch entries: %d\n', height(layer_mismatch_table));
fprintf('Unmatched DAQ return events: %d\n', height(missing_return_event_table));
if ~isempty(mismatch_table)
    fprintf('\nMismatched or incomplete cycles across DAQ / sync crop / saved movie layers:\n');
    disp(mismatch_table(:, {'rec_name', 'cycle_name', ...
        'daq_camera1_return_count', 'daq_camera2_return_count', ...
        'sync_camera1_predicted_crop_frames', 'sync_camera2_predicted_crop_frames', ...
        'saved_camera1_movie_frames', 'saved_camera2_movie_frames', ...
        'daq_return_diff_cam2_minus_cam1', 'sync_predicted_crop_diff_cam2_minus_cam1', ...
        'saved_movie_frame_diff_cam2_minus_cam1', 'status'}));
else
    fprintf('\nAll parsed cycles have matched DAQ camera-return counts.\n');
end

if ~isempty(missing_return_event_table)
    fprintf('\nUnmatched DAQ return events; these identify where a camera return is missing relative to the other channel:\n');
    disp(missing_return_event_table(:, {'rec_name', 'cycle_name', 'missing_in_camera', ...
        'reference_camera', 'reference_frame_number', 'reference_time_s', ...
        'nearest_other_frame_number', 'nearest_dt_ms'}));
end

if ~isempty(layer_mismatch_table)
    fprintf('\nFrame-count mismatch by layer:\n');
    disp(layer_mismatch_table(:, {'rec_name', 'cycle_name', 'method_layer', ...
        'camera1_count', 'camera2_count', 'diff_cam2_minus_cam1', ...
        'likely_missing_in_camera'}));
end

if save_outputs
    output_mat = fullfile(root_path, 'Dual_analysis3_frame_alignment_all_methods.mat');
    output_csv = fullfile(root_path, 'Dual_analysis3_frame_alignment_all_methods.csv');
    layer_csv = fullfile(root_path, 'Dual_analysis3_frame_alignment_layer_mismatches.csv');
    missing_csv = fullfile(root_path, 'Dual_analysis3_DAQ_missing_return_events_from_logs.csv');
    save(output_mat, 'daq_frame_alignment_table', 'mismatch_table', ...
        'layer_mismatch_table', 'missing_return_event_table', 'log_files', '-v7.3');
    writetable(daq_frame_alignment_table, output_csv);
    writetable(layer_mismatch_table, layer_csv);
    writetable(missing_return_event_table, missing_csv);
    fprintf('\nSaved:\n');
    fprintf('  %s\n', output_mat);
    fprintf('  %s\n', output_csv);
    fprintf('  %s\n', layer_csv);
    fprintf('  %s\n', missing_csv);
    if make_plots
        fprintf('  %s\n', plot_output_dir);
    end
end

function row = summarize_one_logs_file(log_file, camera_count, count_saved_movies)
cycle_path = fileparts(log_file);
[rec_name, cycle_name] = infer_rec_cycle_from_path(cycle_path);

S = load(log_file);
if ~isfield(S, 'logs')
    error('logs.mat does not contain variable logs.');
end
logs = S.logs;

row = make_empty_summary_row(camera_count);
row.log_file = string(log_file);
row.cycle_path = string(cycle_path);
row.rec_name = rec_name;
row.cycle_name = cycle_name;
row.cycle_number = parse_cycle_number(cycle_name);

if isfield(logs, 'daq') && is_table_like(logs.daq)
    daq_table = logs.daq;
    row.daq_row_count = height(daq_table);
    row.daq_duration_s = infer_daq_duration_s(daq_table);
    [return_counts, counter_start, counter_end, counter_names] = extract_camera_counts_from_daq(daq_table, camera_count);
    row.daq_camera_counter_names = strjoin(counter_names, ",");
    row.daq_camera1_return_count = return_counts(1);
    row.daq_camera2_return_count = return_counts(2);
    row.daq_camera1_counter_start = counter_start(1);
    row.daq_camera2_counter_start = counter_start(2);
    row.daq_camera1_counter_end = counter_end(1);
    row.daq_camera2_counter_end = counter_end(2);
else
    row.status = "missing_logs_daq";
    row.message = "logs.daq is missing or is not a table.";
end

if isfield(logs, 'sync') && is_table_like(logs.sync)
    sync_table = logs.sync;
    row.sync_row_count = height(sync_table);
    [sync_start, sync_end, sync_diff, sync_names, sync_crop_start, sync_crop_end, sync_crop_frames] = ...
        extract_camera_sync_summary(sync_table, camera_count);
    row.sync_camera_column_names = strjoin(sync_names, ",");
    row.sync_camera1_start = sync_start(1);
    row.sync_camera2_start = sync_start(2);
    row.sync_camera1_end = sync_end(1);
    row.sync_camera2_end = sync_end(2);
    row.sync_camera1_delta = sync_diff(1);
    row.sync_camera2_delta = sync_diff(2);
    row.sync_delta_diff_cam2_minus_cam1 = sync_diff(2) - sync_diff(1);
    row.sync_camera1_predicted_crop_start = sync_crop_start(1);
    row.sync_camera2_predicted_crop_start = sync_crop_start(2);
    row.sync_camera1_predicted_crop_end = sync_crop_end(1);
    row.sync_camera2_predicted_crop_end = sync_crop_end(2);
    row.sync_camera1_predicted_crop_frames = sync_crop_frames(1);
    row.sync_camera2_predicted_crop_frames = sync_crop_frames(2);
    row.sync_predicted_crop_diff_cam2_minus_cam1 = sync_crop_frames(2) - sync_crop_frames(1);
    row.abs_sync_predicted_crop_diff = abs(row.sync_predicted_crop_diff_cam2_minus_cam1);
else
    row.sync_row_count = 0;
end

if count_saved_movies
    [movie_frames, movie_files, movie_methods, movie_messages] = inspect_saved_camera_movies(cycle_path, camera_count);
    row.saved_camera1_movie_file = movie_files(1);
    row.saved_camera2_movie_file = movie_files(2);
    row.saved_camera1_movie_method = movie_methods(1);
    row.saved_camera2_movie_method = movie_methods(2);
    row.saved_camera1_movie_message = movie_messages(1);
    row.saved_camera2_movie_message = movie_messages(2);
    row.saved_camera1_movie_frames = movie_frames(1);
    row.saved_camera2_movie_frames = movie_frames(2);
    row.saved_movie_frame_diff_cam2_minus_cam1 = movie_frames(2) - movie_frames(1);
    row.abs_saved_movie_frame_diff = abs(row.saved_movie_frame_diff_cam2_minus_cam1);
end

if strlength(row.status) == 0
    row.daq_return_diff_cam2_minus_cam1 = row.daq_camera2_return_count - row.daq_camera1_return_count;
    row.abs_daq_return_diff = abs(row.daq_return_diff_cam2_minus_cam1);
    row.common_daq_return_count = min(row.daq_camera1_return_count, row.daq_camera2_return_count);
    row.drop_cam1_to_common = row.daq_camera1_return_count - row.common_daq_return_count;
    row.drop_cam2_to_common = row.daq_camera2_return_count - row.common_daq_return_count;
    layer_diffs = [row.abs_daq_return_diff, row.abs_sync_predicted_crop_diff, row.abs_saved_movie_frame_diff];
    if all(~isfinite(layer_diffs) | layer_diffs == 0)
        row.status = "aligned";
    else
        row.status = "layer_mismatch";
    end
end
end

function row = make_empty_summary_row(camera_count)
row = struct( ...
    'log_file', "", ...
    'cycle_path', "", ...
    'rec_name', "", ...
    'cycle_name', "", ...
    'cycle_number', NaN, ...
    'status', "", ...
    'message', "", ...
    'daq_row_count', NaN, ...
    'daq_duration_s', NaN, ...
    'daq_camera_counter_names', "", ...
    'daq_camera1_return_count', NaN, ...
    'daq_camera2_return_count', NaN, ...
    'daq_return_diff_cam2_minus_cam1', NaN, ...
    'abs_daq_return_diff', NaN, ...
    'common_daq_return_count', NaN, ...
    'drop_cam1_to_common', NaN, ...
    'drop_cam2_to_common', NaN, ...
    'daq_camera1_counter_start', NaN, ...
    'daq_camera2_counter_start', NaN, ...
    'daq_camera1_counter_end', NaN, ...
    'daq_camera2_counter_end', NaN, ...
    'sync_row_count', NaN, ...
    'sync_camera_column_names', "", ...
    'sync_camera1_start', NaN, ...
    'sync_camera2_start', NaN, ...
    'sync_camera1_end', NaN, ...
    'sync_camera2_end', NaN, ...
    'sync_camera1_delta', NaN, ...
    'sync_camera2_delta', NaN, ...
    'sync_delta_diff_cam2_minus_cam1', NaN, ...
    'sync_camera1_predicted_crop_start', NaN, ...
    'sync_camera2_predicted_crop_start', NaN, ...
    'sync_camera1_predicted_crop_end', NaN, ...
    'sync_camera2_predicted_crop_end', NaN, ...
    'sync_camera1_predicted_crop_frames', NaN, ...
    'sync_camera2_predicted_crop_frames', NaN, ...
    'sync_predicted_crop_diff_cam2_minus_cam1', NaN, ...
    'abs_sync_predicted_crop_diff', NaN, ...
    'saved_camera1_movie_file', "", ...
    'saved_camera2_movie_file', "", ...
    'saved_camera1_movie_method', "", ...
    'saved_camera2_movie_method', "", ...
    'saved_camera1_movie_message', "", ...
    'saved_camera2_movie_message', "", ...
    'saved_camera1_movie_frames', NaN, ...
    'saved_camera2_movie_frames', NaN, ...
    'saved_movie_frame_diff_cam2_minus_cam1', NaN, ...
    'abs_saved_movie_frame_diff', NaN);
if camera_count ~= 2
    error('This first-pass checker currently expects camera_count=2.');
end
end

function row = make_empty_missing_event_row()
row = struct( ...
    'log_file', "", ...
    'cycle_path', "", ...
    'rec_name', "", ...
    'cycle_name', "", ...
    'reference_camera', "", ...
    'missing_in_camera', "", ...
    'reference_frame_number', NaN, ...
    'reference_time_s', NaN, ...
    'nearest_other_frame_number', NaN, ...
    'nearest_other_time_s', NaN, ...
    'nearest_dt_ms', NaN, ...
    'match_tolerance_ms', NaN);
end

function row = make_empty_layer_mismatch_row()
row = struct( ...
    'log_file', "", ...
    'cycle_path', "", ...
    'rec_name', "", ...
    'cycle_name', "", ...
    'method_layer', "", ...
    'camera1_count', NaN, ...
    'camera2_count', NaN, ...
    'diff_cam2_minus_cam1', NaN, ...
    'abs_diff', NaN, ...
    'likely_missing_in_camera', "", ...
    'other_camera_reference', "", ...
    'missing_count', NaN, ...
    'missing_range_start_after_common', NaN, ...
    'missing_range_end_reference', NaN);
end

function layer_rows = build_layer_mismatch_rows(summary_table)
layer_rows = repmat(make_empty_layer_mismatch_row(), 0, 1);
for idx = 1:height(summary_table)
    current_row = summary_table(idx, :);
    layer_rows = [layer_rows; make_layer_mismatch_row_if_needed(current_row, ...
        "daq_return_counter", current_row.daq_camera1_return_count, current_row.daq_camera2_return_count)]; %#ok<AGROW>
    layer_rows = [layer_rows; make_layer_mismatch_row_if_needed(current_row, ...
        "sync_predicted_crop", current_row.sync_camera1_predicted_crop_frames, current_row.sync_camera2_predicted_crop_frames)]; %#ok<AGROW>
    layer_rows = [layer_rows; make_layer_mismatch_row_if_needed(current_row, ...
        "saved_movie_file", current_row.saved_camera1_movie_frames, current_row.saved_camera2_movie_frames)]; %#ok<AGROW>
end
end

function layer_row = make_layer_mismatch_row_if_needed(summary_row, method_layer, camera1_count, camera2_count)
layer_row = repmat(make_empty_layer_mismatch_row(), 0, 1);
if ~isfinite(camera1_count) || ~isfinite(camera2_count) || camera1_count == camera2_count
    return;
end

new_row = make_empty_layer_mismatch_row();
new_row.log_file = string(summary_row.log_file);
new_row.cycle_path = string(summary_row.cycle_path);
new_row.rec_name = string(summary_row.rec_name);
new_row.cycle_name = string(summary_row.cycle_name);
new_row.method_layer = method_layer;
new_row.camera1_count = camera1_count;
new_row.camera2_count = camera2_count;
new_row.diff_cam2_minus_cam1 = camera2_count - camera1_count;
new_row.abs_diff = abs(new_row.diff_cam2_minus_cam1);
common_count = min(camera1_count, camera2_count);
new_row.missing_range_start_after_common = common_count + 1;
if camera1_count < camera2_count
    new_row.likely_missing_in_camera = "Camera1";
    new_row.other_camera_reference = "Camera2";
    new_row.missing_count = camera2_count - camera1_count;
    new_row.missing_range_end_reference = camera2_count;
else
    new_row.likely_missing_in_camera = "Camera2";
    new_row.other_camera_reference = "Camera1";
    new_row.missing_count = camera1_count - camera2_count;
    new_row.missing_range_end_reference = camera1_count;
end
layer_row = new_row;
end

function tf = has_any_layer_mismatch(summary_row)
layer_diffs = [summary_row.abs_daq_return_diff, ...
    summary_row.abs_sync_predicted_crop_diff, ...
    summary_row.abs_saved_movie_frame_diff];
tf = any(isfinite(layer_diffs) & layer_diffs > 0) || summary_row.status ~= "aligned";
end

function event_report = analyze_camera_return_events_from_logs(logs, camera_count, frame_match_tolerance_s)
if ~isfield(logs, 'daq') || ~is_table_like(logs.daq)
    error('logs.daq is missing or is not table-like.');
end
daq_table = logs.daq;
[~, ~, ~, counter_names] = extract_camera_counts_from_daq(daq_table, camera_count);
time_s = get_daq_time_seconds(daq_table);

events = repmat(struct( ...
    'camera_name', "", ...
    'counter_name', "", ...
    'time_s', [], ...
    'frame_number', [], ...
    'sample_index', []), 1, camera_count);

for camera_idx = 1:camera_count
    [event_time_s, frame_number, sample_index] = extract_camera_return_events( ...
        daq_table, counter_names(camera_idx), time_s);
    events(camera_idx).camera_name = "Camera" + camera_idx;
    events(camera_idx).counter_name = counter_names(camera_idx);
    events(camera_idx).time_s = event_time_s;
    events(camera_idx).frame_number = frame_number;
    events(camera_idx).sample_index = sample_index;
end

if isempty(frame_match_tolerance_s)
    tolerance_s = estimate_frame_match_tolerance_s(events, time_s);
else
    tolerance_s = double(frame_match_tolerance_s);
end

matches = match_two_event_trains(events(1), events(2), tolerance_s);
event_report = struct( ...
    'daq_time_s', time_s, ...
    'camera_counter_names', counter_names, ...
    'camera1_counter_values', double(daq_table.(counter_names(1))), ...
    'camera2_counter_values', double(daq_table.(counter_names(2))), ...
    'events', events, ...
    'match_tolerance_s', tolerance_s, ...
    'matches', matches);
end

function [event_time_s, frame_number, sample_index] = extract_camera_return_events(daq_table, counter_name, time_s)
counter_values = double(daq_table.(counter_name));
valid = isfinite(counter_values) & isfinite(time_s);
counter_values = counter_values(valid);
valid_indices = find(valid);
time_s = time_s(valid);

if isempty(counter_values)
    event_time_s = [];
    frame_number = [];
    sample_index = [];
    return;
end

counter_step = diff(counter_values);
step_indices = find(counter_step > 0) + 1;
event_time_s = zeros(0, 1);
frame_number = zeros(0, 1);
sample_index = zeros(0, 1);
for idx = 1:numel(step_indices)
    current_idx = step_indices(idx);
    previous_count = counter_values(current_idx - 1);
    current_count = counter_values(current_idx);
    new_frames = (previous_count + 1):current_count;
    event_time_s = [event_time_s; repmat(time_s(current_idx), numel(new_frames), 1)]; %#ok<AGROW>
    frame_number = [frame_number; new_frames(:)]; %#ok<AGROW>
    sample_index = [sample_index; repmat(valid_indices(current_idx), numel(new_frames), 1)]; %#ok<AGROW>
end
end

function time_s = get_daq_time_seconds(daq_table)
if istimetable(daq_table)
    t = daq_table.Properties.RowTimes;
elseif ismember("Time", string(daq_table.Properties.VariableNames))
    t = daq_table.Time;
else
    t = (0:height(daq_table)-1)';
end

if isduration(t)
    time_s = seconds(t);
elseif isdatetime(t)
    time_s = seconds(t - t(1));
else
    time_s = double(t);
end
time_s = time_s(:);
end

function tolerance_s = estimate_frame_match_tolerance_s(events, time_s)
isi = [];
for idx = 1:numel(events)
    isi = [isi; diff(events(idx).time_s(:))]; %#ok<AGROW>
end
isi = isi(isfinite(isi) & isi > 0);
daq_dt = diff(time_s(:));
daq_dt = daq_dt(isfinite(daq_dt) & daq_dt > 0);
if isempty(isi)
    frame_period_s = 1 / 400;
else
    frame_period_s = median(isi);
end
if isempty(daq_dt)
    daq_sample_s = frame_period_s / 5;
else
    daq_sample_s = median(daq_dt);
end
tolerance_s = min(0.45 * frame_period_s, max(2.5 * daq_sample_s, 0.00075));
end

function matches = match_two_event_trains(cam1_events, cam2_events, tolerance_s)
t1 = cam1_events.time_s(:);
t2 = cam2_events.time_s(:);
i = 1;
j = 1;
paired_cam1_idx = zeros(0, 1);
paired_cam2_idx = zeros(0, 1);
unmatched_cam1_idx = zeros(0, 1);
unmatched_cam2_idx = zeros(0, 1);

while i <= numel(t1) && j <= numel(t2)
    dt = t2(j) - t1(i);
    if abs(dt) <= tolerance_s
        paired_cam1_idx(end+1, 1) = i; %#ok<AGROW>
        paired_cam2_idx(end+1, 1) = j; %#ok<AGROW>
        i = i + 1;
        j = j + 1;
    elseif t1(i) < t2(j)
        unmatched_cam1_idx(end+1, 1) = i; %#ok<AGROW>
        i = i + 1;
    else
        unmatched_cam2_idx(end+1, 1) = j; %#ok<AGROW>
        j = j + 1;
    end
end
if i <= numel(t1)
    unmatched_cam1_idx = [unmatched_cam1_idx; (i:numel(t1))'];
end
if j <= numel(t2)
    unmatched_cam2_idx = [unmatched_cam2_idx; (j:numel(t2))'];
end

matches = struct( ...
    'paired_cam1_idx', paired_cam1_idx, ...
    'paired_cam2_idx', paired_cam2_idx, ...
    'paired_cam1_frame', cam1_events.frame_number(paired_cam1_idx), ...
    'paired_cam2_frame', cam2_events.frame_number(paired_cam2_idx), ...
    'paired_dt_s_cam2_minus_cam1', t2(paired_cam2_idx) - t1(paired_cam1_idx), ...
    'unmatched_cam1_idx', unmatched_cam1_idx, ...
    'unmatched_cam2_idx', unmatched_cam2_idx);
end

function missing_rows = build_missing_event_rows(summary_row, event_report)
missing_rows = repmat(make_empty_missing_event_row(), 0, 1);
cam1_events = event_report.events(1);
cam2_events = event_report.events(2);

for idx = event_report.matches.unmatched_cam1_idx(:)'
    missing_rows(end+1, 1) = make_missing_event_row(summary_row, event_report, ...
        "Camera1", "Camera2", cam1_events, idx, cam2_events); %#ok<AGROW>
end
for idx = event_report.matches.unmatched_cam2_idx(:)'
    missing_rows(end+1, 1) = make_missing_event_row(summary_row, event_report, ...
        "Camera2", "Camera1", cam2_events, idx, cam1_events); %#ok<AGROW>
end
end

function row = make_missing_event_row(summary_row, event_report, reference_camera, missing_in_camera, reference_events, reference_idx, other_events)
row = make_empty_missing_event_row();
row.log_file = string(summary_row.log_file);
row.cycle_path = string(summary_row.cycle_path);
row.rec_name = string(summary_row.rec_name);
row.cycle_name = string(summary_row.cycle_name);
row.reference_camera = reference_camera;
row.missing_in_camera = missing_in_camera;
row.reference_frame_number = reference_events.frame_number(reference_idx);
row.reference_time_s = reference_events.time_s(reference_idx);
row.match_tolerance_ms = event_report.match_tolerance_s * 1000;

if ~isempty(other_events.time_s)
    [nearest_dt, nearest_idx] = min(abs(other_events.time_s(:) - row.reference_time_s));
    row.nearest_other_frame_number = other_events.frame_number(nearest_idx);
    row.nearest_other_time_s = other_events.time_s(nearest_idx);
    row.nearest_dt_ms = nearest_dt * 1000;
end
end

function plot_daq_frame_alignment_event_report(summary_row, event_report, plot_file_base)
fig = figure('Color', 'w', 'Name', char(summary_row.rec_name + " " + summary_row.cycle_name + " DAQ Frame Alignment"));
t = event_report.daq_time_s;
counter_diff = event_report.camera2_counter_values - event_report.camera1_counter_values;
unmatched_cam1 = event_report.matches.unmatched_cam1_idx;
unmatched_cam2 = event_report.matches.unmatched_cam2_idx;
unmatched_times = [event_report.events(1).time_s(unmatched_cam1); event_report.events(2).time_s(unmatched_cam2)];

subplot(3, 1, 1);
paired_frame = event_report.matches.paired_cam1_frame;
paired_dt_ms = event_report.matches.paired_dt_s_cam2_minus_cam1 * 1000;
plot(paired_frame, paired_dt_ms, '.', 'Color', [0.2 0.2 0.2], 'MarkerSize', 4);
hold on;
yline(event_report.match_tolerance_s * 1000, 'r--', 'Tolerance');
yline(-event_report.match_tolerance_s * 1000, 'r--');
grid on;
xlabel('Paired Camera1 frame number');
ylabel('Cam2 - Cam1 time (ms)');
title(sprintf('%s / %s | Paired DAQ return timing', summary_row.rec_name, summary_row.cycle_name), 'Interpreter', 'none');

subplot(3, 1, 2);
plot(t, counter_diff, 'k-', 'LineWidth', 1);
hold on;
if ~isempty(unmatched_times)
    for idx = 1:numel(unmatched_times)
        xline(unmatched_times(idx), 'r-', 'LineWidth', 0.8);
    end
end
grid on;
xlabel('DAQ time (s)');
ylabel('Cam2 count - Cam1 count');
title('Cumulative DAQ return count difference');

subplot(3, 1, 3);
[x_min, x_max] = choose_event_raster_window(event_report, unmatched_times);
cam1_mask = event_report.events(1).time_s >= x_min & event_report.events(1).time_s <= x_max;
cam2_mask = event_report.events(2).time_s >= x_min & event_report.events(2).time_s <= x_max;
plot(event_report.events(1).time_s(cam1_mask), ones(sum(cam1_mask), 1), '.', 'Color', [0 0.45 0.74], 'MarkerSize', 8);
hold on;
plot(event_report.events(2).time_s(cam2_mask), 2 * ones(sum(cam2_mask), 1), '.', 'Color', [0.85 0.33 0.10], 'MarkerSize', 8);
if ~isempty(unmatched_cam1)
    plot(event_report.events(1).time_s(unmatched_cam1), ones(numel(unmatched_cam1), 1), 'o', ...
        'Color', [0 0.45 0.74], 'MarkerSize', 7, 'LineWidth', 1.4);
end
if ~isempty(unmatched_cam2)
    plot(event_report.events(2).time_s(unmatched_cam2), 2 * ones(numel(unmatched_cam2), 1), 'o', ...
        'Color', [0.85 0.33 0.10], 'MarkerSize', 7, 'LineWidth', 1.4);
end
ylim([0.5 2.5]);
yticks([1 2]);
yticklabels({'Camera1', 'Camera2'});
xlim([x_min x_max]);
grid on;
xlabel('DAQ time (s)');
title('Return-event raster around unmatched frames');

annotation_text = build_missing_annotation(summary_row, event_report);
annotation(fig, 'textbox', [0.13 0.01 0.78 0.06], 'String', annotation_text, ...
    'EdgeColor', 'none', 'Interpreter', 'none', 'FontSize', 9);

savefig(fig, [char(plot_file_base) '.fig']);
saveas(fig, [char(plot_file_base) '.png']);
end

function [x_min, x_max] = choose_event_raster_window(event_report, unmatched_times)
if isempty(unmatched_times)
    all_times = [event_report.events(1).time_s(:); event_report.events(2).time_s(:)];
    if isempty(all_times)
        x_min = 0;
        x_max = 1;
    else
        x_min = min(all_times);
        x_max = min(max(all_times), x_min + 1);
    end
else
    x_min = max(0, min(unmatched_times) - 0.05);
    x_max = max(unmatched_times) + 0.05;
    if x_max <= x_min
        x_max = x_min + 0.1;
    end
end
end

function annotation_text = build_missing_annotation(summary_row, event_report)
cam1_missing_in_cam2 = event_report.events(1).frame_number(event_report.matches.unmatched_cam1_idx);
cam2_missing_in_cam1 = event_report.events(2).frame_number(event_report.matches.unmatched_cam2_idx);
layer_text = sprintf('DAQ returns Cam1/Cam2=%g/%g | sync predicted crop=%g/%g | saved movie frames=%g/%g. ', ...
    summary_row.daq_camera1_return_count, summary_row.daq_camera2_return_count, ...
    summary_row.sync_camera1_predicted_crop_frames, summary_row.sync_camera2_predicted_crop_frames, ...
    summary_row.saved_camera1_movie_frames, summary_row.saved_camera2_movie_frames);
missing_text = sprintf('DAQ unmatched: missing in Cam2 relative to Cam1 frames: %s | missing in Cam1 relative to Cam2 frames: %s', ...
    compact_number_list(cam1_missing_in_cam2), compact_number_list(cam2_missing_in_cam1));
annotation_text = [layer_text, missing_text];
end

function text_out = compact_number_list(values)
values = unique(values(:)');
if isempty(values)
    text_out = 'none';
elseif numel(values) <= 12
    text_out = strjoin(string(values), ', ');
else
    text_out = strjoin([string(values(1:10)), "...", string(values(end))], ', ');
end
end

function [return_counts, counter_start, counter_end, counter_names] = extract_camera_counts_from_daq(daq_table, camera_count)
names = string(daq_table.Properties.VariableNames);
names = names(names ~= "Time");
if numel(names) < camera_count + 1
    error('logs.daq has too few columns to contain %d camera counters.', camera_count);
end

% Hamamatsu_Imaging_2_Cameras_rebuilt adds ai0 first, then one EdgeCount
% input per camera, then the PTB/stim counter. Thus camera counters are the
% second and third table variables for two-camera recordings.
counter_names = names(2:(camera_count + 1));
return_counts = NaN(1, camera_count);
counter_start = NaN(1, camera_count);
counter_end = NaN(1, camera_count);
for idx = 1:camera_count
    values = double(daq_table.(counter_names(idx)));
    values = values(isfinite(values));
    if isempty(values)
        continue;
    end
    counter_start(idx) = values(1);
    counter_end(idx) = values(end);
    % EdgeCount channels are reset before acquisition in the rebuilt script,
    % so the final counter value is the DAQ-received camera return count.
    % The start value is kept in the output table for diagnosing unexpected
    % nonzero starts instead of silently subtracting them away.
    return_counts(idx) = counter_end(idx);
end
end

function [sync_start, sync_end, sync_delta, sync_names, sync_crop_start, sync_crop_end, sync_crop_frames] = extract_camera_sync_summary(sync_table, camera_count)
names = string(sync_table.Properties.VariableNames);
sync_names = names(contains(names, "Camera") & contains(names, "Frame"));
if numel(sync_names) < camera_count
    sync_names = strings(1, camera_count);
    sync_start = NaN(1, camera_count);
    sync_end = NaN(1, camera_count);
    sync_delta = NaN(1, camera_count);
    sync_crop_start = NaN(1, camera_count);
    sync_crop_end = NaN(1, camera_count);
    sync_crop_frames = NaN(1, camera_count);
    return;
end
sync_names = sync_names(1:camera_count);
sync_start = NaN(1, camera_count);
sync_end = NaN(1, camera_count);
sync_delta = NaN(1, camera_count);
sync_crop_start = NaN(1, camera_count);
sync_crop_end = NaN(1, camera_count);
sync_crop_frames = NaN(1, camera_count);
for idx = 1:camera_count
    values = double(sync_table.(sync_names(idx)));
    values = values(isfinite(values));
    if isempty(values)
        continue;
    end
    sync_start(idx) = values(1);
    sync_end(idx) = values(end);
    sync_delta(idx) = sync_end(idx) - sync_start(idx);
    if numel(values) >= 2
        sync_crop_start(idx) = values(2);
        sync_crop_end(idx) = values(end) + ceil(mean(diff(values)));
        sync_crop_frames(idx) = sync_crop_end(idx) - sync_crop_start(idx) + 1;
    end
end
end

function [movie_frames, movie_files, movie_methods, movie_messages] = inspect_saved_camera_movies(cycle_path, camera_count)
movie_frames = NaN(1, camera_count);
movie_files = strings(1, camera_count);
movie_methods = strings(1, camera_count);
movie_messages = strings(1, camera_count);
for camera_idx = 1:camera_count
    camera_movie_files = find_saved_camera_movie_files(cycle_path, camera_idx);
    movie_files(camera_idx) = strjoin(camera_movie_files, ";");
    if isempty(camera_movie_files)
        movie_messages(camera_idx) = "No saved movie file found under Cam folder.";
        continue;
    end
    try
        total_frames = 0;
        methods = strings(0, 1);
        for file_idx = 1:numel(camera_movie_files)
            [frame_count, method] = count_movie_file_frames(char(camera_movie_files(file_idx)));
            total_frames = total_frames + frame_count;
            methods(end+1, 1) = method; %#ok<AGROW>
        end
        movie_frames(camera_idx) = total_frames;
        movie_methods(camera_idx) = strjoin(unique(methods, 'stable'), "+");
        movie_messages(camera_idx) = sprintf('%d stack file(s) summed in filename order.', numel(camera_movie_files));
    catch ME
        movie_messages(camera_idx) = string(ME.message);
    end
end
end

function movie_files = find_saved_camera_movie_files(cycle_path, camera_idx)
movie_files = strings(0, 1);
cam_dirs = dir(fullfile(cycle_path, sprintf('Cam%d*', camera_idx)));
cam_dirs = cam_dirs([cam_dirs.isdir]);
if isempty(cam_dirs)
    return;
end

candidate_listing = [];
for dir_idx = 1:numel(cam_dirs)
    cam_path = fullfile(cam_dirs(dir_idx).folder, cam_dirs(dir_idx).name);
    candidate_listing = [candidate_listing; dir(fullfile(cam_path, '*.tif'))]; %#ok<AGROW>
    candidate_listing = [candidate_listing; dir(fullfile(cam_path, '*.tiff'))]; %#ok<AGROW>
    candidate_listing = [candidate_listing; dir(fullfile(cam_path, '*.mat'))]; %#ok<AGROW>
end
candidate_listing = candidate_listing(~[candidate_listing.isdir]);
if isempty(candidate_listing)
    return;
end
candidate_names = string({candidate_listing.name});
[~, file_order] = sort(candidate_names);
candidate_listing = candidate_listing(file_order);
movie_files = string(fullfile({candidate_listing.folder}, {candidate_listing.name}))';
end

function [frame_count, method] = count_movie_file_frames(movie_file)
[~, ~, ext] = fileparts(movie_file);
switch lower(ext)
    case {'.tif', '.tiff'}
        frame_count = count_tiff_directories(movie_file);
        method = "tiff_directory_count";
    case '.mat'
        frame_count = count_mat_movie_frames(movie_file);
        method = "mat_movie_variable_size";
    otherwise
        error('Unsupported movie file extension: %s', ext);
end
end

function frame_count = count_tiff_directories(tif_file)
frame_count = 0;
tiff_obj = Tiff(tif_file, 'r');
cleanup_obj = onCleanup(@() close(tiff_obj));
while true
    frame_count = frame_count + 1;
    if tiff_obj.lastDirectory()
        break;
    end
    tiff_obj.nextDirectory();
end
end

function frame_count = count_mat_movie_frames(mat_file)
vars = whos('-file', mat_file);
if isempty(vars)
    frame_count = NaN;
    return;
end
movie_idx = find(strcmp({vars.name}, 'movie'), 1, 'first');
if isempty(movie_idx)
    numeric_vars = vars(ismember({vars.class}, {'single', 'double', 'uint8', 'uint16', 'uint32', 'int16'}));
    if isempty(numeric_vars)
        frame_count = NaN;
        return;
    end
    [~, movie_idx] = max([numeric_vars.bytes]);
    selected_var = numeric_vars(movie_idx);
else
    selected_var = vars(movie_idx);
end

sz = selected_var.size;
if numel(sz) >= 4
    frame_count = sz(4);
elseif numel(sz) >= 3
    frame_count = sz(3);
elseif numel(sz) >= 2
    frame_count = sz(2);
else
    frame_count = sz(1);
end
end

function duration_s = infer_daq_duration_s(daq_table)
duration_s = NaN;
try
    time_s = get_daq_time_seconds(daq_table);
    if ~isempty(time_s)
        duration_s = time_s(end) - time_s(1);
    end
catch
    duration_s = NaN;
end
end

function tf = is_table_like(x)
tf = istable(x) || istimetable(x);
end

function file_paths = find_files_recursive(root_path, pattern)
file_paths = strings(0, 1);
direct_listing = dir(fullfile(root_path, pattern));
direct_listing = direct_listing(~[direct_listing.isdir]);
for idx = 1:numel(direct_listing)
    file_paths(end+1, 1) = string(fullfile(direct_listing(idx).folder, direct_listing(idx).name)); %#ok<AGROW>
end

child_dirs = dir(root_path);
child_dirs = child_dirs([child_dirs.isdir]);
child_dirs = child_dirs(~ismember({child_dirs.name}, {'.', '..'}));
for idx = 1:numel(child_dirs)
    child_path = fullfile(child_dirs(idx).folder, child_dirs(idx).name);
    child_files = find_files_recursive(child_path, pattern);
    file_paths = [file_paths; child_files]; %#ok<AGROW>
end
file_paths = unique(file_paths, 'stable');
end

function [rec_name, cycle_name] = infer_rec_cycle_from_path(cycle_path)
parts = regexp(char(cycle_path), '[\\/]', 'split');
rec_name = "";
cycle_name = "";
if ~isempty(parts)
    cycle_name = string(parts{end});
end
if numel(parts) >= 2
    rec_name = string(parts{end-1});
end
end

function cycle_number = parse_cycle_number(cycle_name)
tokens = regexp(char(cycle_name), 'Cycle(\d+)', 'tokens', 'once');
if isempty(tokens)
    cycle_number = NaN;
else
    cycle_number = str2double(tokens{1});
end
end
