% AP_analysis3_rec_batch - Run AP_analysis3_rec across all cycles in one record
%
% Workflow:
% 1. Choose one reference cycle.
% 2. Create or reuse one ROI file from that reference cycle.
% 3. Reuse the same ROI file for every other Cycle* under the same record.
% 4. Reopen saved per-cycle traces and build one record-level average summary.
%
% This script keeps AP_analysis3_rec as a normal runnable script.
% It passes parameters in from an outer scope by calling run(...).

clc;

%% Input Setup
if ~exist('rec_path', 'var') || isempty(rec_path)
    rec_path = 'E:\1_Data\CC\20260520_CC_slice_Vega\ROI1\Methods2_default\Rec1_2026-05-21_18-38-41';
end
if ~exist('analysis_backend', 'var') || isempty(analysis_backend)
    analysis_backend = "volpy"; %classic
else
    analysis_backend = string(analysis_backend);
end
if ~exist('source_entry_name', 'var') || isempty(source_entry_name)
    source_entry_name = "";
else
    source_entry_name = string(source_entry_name);
end
if ~exist('reference_cycle_name', 'var') || isempty(reference_cycle_name)
    reference_cycle_name = 'Cycle1';
end
if ~exist('reference_roi_file', 'var') || isempty(reference_roi_file)
    reference_roi_file = "";
else
    reference_roi_file = string(reference_roi_file);
end
if ~exist('reference_mask_method_override', 'var') || isempty(reference_mask_method_override)
    reference_mask_method_override = "";
else
    reference_mask_method_override = string(reference_mask_method_override);
end
if ~exist('reuse_existing_reference_roi', 'var') || isempty(reuse_existing_reference_roi)
    reuse_existing_reference_roi = true;
end
if ~exist('rerun_reference_cycle_for_roi', 'var') || isempty(rerun_reference_cycle_for_roi)
    rerun_reference_cycle_for_roi = false;
end
if ~exist('run_reference_cycle_in_batch', 'var') || isempty(run_reference_cycle_in_batch)
    run_reference_cycle_in_batch = true;
end
if ~exist('cycle_name_filter', 'var') || isempty(cycle_name_filter)
    cycle_name_filter = strings(0, 1);
else
    cycle_name_filter = string(cycle_name_filter(:));
end
if ~exist('skip_cycles_with_existing_results', 'var') || isempty(skip_cycles_with_existing_results)
    skip_cycles_with_existing_results = true;
end
if ~exist('run_record_average_only', 'var') || isempty(run_record_average_only)
    run_record_average_only = false;
end
if ~exist('trace_polarity', 'var') || isempty(trace_polarity)
    trace_polarity = 1;
end
if ~exist('run_manual_peak_refinement', 'var') || isempty(run_manual_peak_refinement)
    run_manual_peak_refinement = false;
end
if ~exist('run_manual_peak_gating', 'var') || isempty(run_manual_peak_gating)
    run_manual_peak_gating = false;
end
if ~exist('show_roi_saved_dialog', 'var') || isempty(show_roi_saved_dialog)
    show_roi_saved_dialog = false;
end
if ~exist('gpu_override', 'var') || isempty(gpu_override)
    gpu_override = false;
end

average_output_dir_name = 'AP_analysis3_rec_average';
batch_summary_file = fullfile(rec_path, 'AP_analysis3_rec_batch_summary.mat');

%% Resolve Cycles
script_dir = fileparts(mfilename('fullpath'));
ap_script_path = fullfile(script_dir, 'AP_analysis3_rec.m');
if ~isfile(ap_script_path)
    error('Cannot find AP_analysis3_rec.m next to this batch script.');
end
if ~isfolder(rec_path)
    error('Record path does not exist: %s', rec_path);
end

cycle_dirs = dir(fullfile(rec_path, 'Cycle*'));
cycle_dirs = cycle_dirs([cycle_dirs.isdir]);
cycle_dirs = sort_cycle_dirs_ap(cycle_dirs);
if isempty(cycle_dirs)
    error('No Cycle* folders found under record path: %s', rec_path);
end

all_cycle_names = string({cycle_dirs.name});
if ~isempty(cycle_name_filter)
    keep_mask = ismember(all_cycle_names, cycle_name_filter);
    cycle_dirs = cycle_dirs(keep_mask);
end
if isempty(cycle_dirs)
    error('No cycles remain after applying cycle_name_filter.');
end

%% Summary-Only Mode
if run_record_average_only
    batch_results = repmat(struct( ...
        'cycle_name', "", ...
        'cycle_path', "", ...
        'status', "", ...
        'save_path', "", ...
        'roi_file', "", ...
        'source_entry_name', "", ...
        'message', ""), 0, 1);

    for i = 1:numel(cycle_dirs)
        current_cycle_name = string(cycle_dirs(i).name);
        current_cycle_path = fullfile(cycle_dirs(i).folder, cycle_dirs(i).name);
        existing_result = find_latest_explicit_result_ap(current_cycle_path);
        if strlength(existing_result) > 0
            batch_results(end+1, 1) = struct( ...
                'cycle_name', current_cycle_name, ...
                'cycle_path', string(current_cycle_path), ...
                'status', "existing_result_detected", ...
                'save_path', fileparts(existing_result), ...
                'roi_file', "", ...
                'source_entry_name', "", ...
                'message', existing_result);
        else
            batch_results(end+1, 1) = struct( ...
                'cycle_name', current_cycle_name, ...
                'cycle_path', string(current_cycle_path), ...
                'status', "missing_result", ...
                'save_path', "", ...
                'roi_file', "", ...
                'source_entry_name', "", ...
                'message', "No -1_explicit_results.mat found under this cycle.");
        end
    end
else
    %% Existing Result Detection
    fprintf('\n============================================================\n');
    fprintf('[Batch] Existing Result Detection\n');
    fprintf('============================================================\n');
    for i = 1:numel(cycle_dirs)
        cycle_path_i = fullfile(cycle_dirs(i).folder, cycle_dirs(i).name);
        existing_run_i = find_latest_explicit_result_ap(cycle_path_i);
        if strlength(existing_run_i) > 0
            fprintf('%s | existing processed run detected:\n  %s\n', cycle_dirs(i).name, existing_run_i);
        else
            fprintf('%s | no processed AP_analysis3_rec run detected.\n', cycle_dirs(i).name);
        end
    end

    reference_cycle_path = fullfile(rec_path, reference_cycle_name);
    if ~isfolder(reference_cycle_path)
        error('Reference cycle folder does not exist: %s', reference_cycle_path);
    end

    %% Resolve Or Create Reference ROI
    fprintf('\n============================================================\n');
    fprintf('[Batch] Resolve Reference ROI\n');
    fprintf('============================================================\n');
    fprintf('Record path: %s\n', rec_path);
    fprintf('Reference cycle: %s\n', reference_cycle_name);

    if strlength(reference_roi_file) == 0 && reuse_existing_reference_roi && ~rerun_reference_cycle_for_roi
        reference_roi_file = string(find_latest_roi_file_ap(reference_cycle_path));
        if strlength(reference_roi_file) > 0
            fprintf('Reusing latest ROI file already present in reference cycle:\n  %s\n', reference_roi_file);
        end
    end

    if strlength(reference_roi_file) == 0 || rerun_reference_cycle_for_roi
        fprintf('Running reference cycle to create/recreate ROI...\n');
        reference_result = run_ap_cycle_with_overrides( ...
            ap_script_path, reference_cycle_path, source_entry_name, "", ...
            reference_mask_method_override, ...
            gpu_override, ...
            analysis_backend, trace_polarity, ...
            run_manual_peak_refinement, run_manual_peak_gating, show_roi_saved_dialog);
        reference_roi_file = string(reference_result.roi_file);
        fprintf('Reference ROI file created:\n  %s\n', reference_roi_file);
    end

    if strlength(reference_roi_file) == 0 || ~isfile(reference_roi_file)
        error('Reference ROI file could not be resolved.');
    end

    %% Run All Cycles
    fprintf('\n============================================================\n');
    fprintf('[Batch] Run Record Cycles\n');
    fprintf('============================================================\n');

    batch_results = repmat(struct( ...
        'cycle_name', "", ...
        'cycle_path', "", ...
        'status', "", ...
        'save_path', "", ...
        'roi_file', "", ...
        'source_entry_name', "", ...
        'message', ""), 0, 1);

    for i = 1:numel(cycle_dirs)
        current_cycle_name = string(cycle_dirs(i).name);
        current_cycle_path = fullfile(cycle_dirs(i).folder, cycle_dirs(i).name);

        if strcmpi(current_cycle_name, string(reference_cycle_name)) && ~run_reference_cycle_in_batch
            batch_results(end+1, 1) = struct( ...
                'cycle_name', current_cycle_name, ...
                'cycle_path', string(current_cycle_path), ...
                'status', "skipped_reference_already_used_for_roi", ...
                'save_path', "", ...
                'roi_file', reference_roi_file, ...
                'source_entry_name', "", ...
                'message', "Reference cycle already used to create/reuse ROI.");
            fprintf('Skipping %s in batch loop because it already served as the ROI reference.\n', current_cycle_name);
            continue;
        end

        if skip_cycles_with_existing_results
            existing_result = find_latest_explicit_result_ap(current_cycle_path);
            if strlength(existing_result) > 0
                batch_results(end+1, 1) = struct( ...
                    'cycle_name', current_cycle_name, ...
                    'cycle_path', string(current_cycle_path), ...
                    'status', "skipped_existing_result", ...
                    'save_path', fileparts(existing_result), ...
                    'roi_file', reference_roi_file, ...
                    'source_entry_name', "", ...
                    'message', existing_result);
                fprintf('Skipping %s because the final output already exists:\n  %s\n', current_cycle_name, existing_result);
                continue;
            end
        end

        fprintf('\n[Batch] Running %s ...\n', current_cycle_name);
        try
            run_result = run_ap_cycle_with_overrides( ...
                ap_script_path, current_cycle_path, source_entry_name, reference_roi_file, ...
                reference_mask_method_override, ...
                gpu_override, ...
                analysis_backend, trace_polarity, ...
                run_manual_peak_refinement, run_manual_peak_gating, show_roi_saved_dialog);

            batch_results(end+1, 1) = struct( ...
                'cycle_name', current_cycle_name, ...
                'cycle_path', string(current_cycle_path), ...
                'status', "completed", ...
                'save_path', string(run_result.save_path), ...
                'roi_file', string(run_result.roi_file), ...
                'source_entry_name', string(run_result.source_entry_name), ...
                'message', "");
        catch ME
            batch_results(end+1, 1) = struct( ...
                'cycle_name', current_cycle_name, ...
                'cycle_path', string(current_cycle_path), ...
                'status', "failed", ...
                'save_path', "", ...
                'roi_file', reference_roi_file, ...
                'source_entry_name', "", ...
                'message', string(ME.message));
            fprintf(2, '[Batch] %s failed: %s\n', current_cycle_name, ME.message);
        end
    end

    save(batch_summary_file, 'batch_results', 'reference_roi_file');
    fprintf('\nBatch summary saved to:\n  %s\n', batch_summary_file);
end

%% Record-Level Average Trace Summary
fprintf('\n============================================================\n');
fprintf('[Batch] Record-Level Average Trace Summary\n');
fprintf('============================================================\n');
record_average_output = build_record_average_summary_ap( ...
    rec_path, batch_results, reference_roi_file, trace_polarity, average_output_dir_name);
fprintf('Record-average summary folder:\n  %s\n', record_average_output.output_dir);
fprintf('Record-average result bundle:\n  %s\n', record_average_output.results_file);

%% Local Functions
function result = run_ap_cycle_with_overrides( ...
    ap_script_path, cycle_path, source_entry_name, reuse_roi_file, reference_mask_method_override, ...
    gpu_override, ...
    analysis_backend, trace_polarity, ...
    run_manual_peak_refinement, run_manual_peak_gating, show_roi_saved_dialog)

[source_entry_path, source_entry_resolved] = resolve_cycle_source_entry_ap(cycle_path, source_entry_name);
[source_parent, source_leaf, source_ext] = fileparts(source_entry_path);
if isfolder(source_entry_path)
    folder_path = cycle_path;
    file = char(source_entry_resolved);
else
    folder_path = source_parent;
    file = [source_leaf, source_ext];
end

[record_path_for_name, cycle_name_for_name] = fileparts(cycle_path);
[~, record_name_for_name] = fileparts(record_path_for_name);
time_tag = string(datetime('now', 'Format', 'yyyy-MM-dd HH-mm-ss-SSS'));
analysis_run_name = sprintf('%s_%s_%s', record_name_for_name, cycle_name_for_name, char(time_tag));

mask_file_override = "";
mask_method_override = "";
if strlength(string(reuse_roi_file)) > 0
    mask_method_override = "previous ROIs";
    mask_file_override = string(reuse_roi_file);
elseif strlength(string(reference_mask_method_override)) > 0
    mask_method_override = string(reference_mask_method_override);
end
gpu = gpu_override;

run(ap_script_path);

result = struct( ...
    'save_path', string(save_path), ...
    'roi_file', string(fullfile(save_path, '1_raw_ROI.mat')), ...
    'source_entry_name', string(source_entry_resolved));
end

function [source_entry_path, source_entry_name] = resolve_cycle_source_entry_ap(cycle_path, source_entry_name)
if strlength(string(source_entry_name)) > 0
    source_entry_path = fullfile(cycle_path, char(source_entry_name));
    if ~isfolder(source_entry_path) && ~isfile(source_entry_path)
        error('source_entry_name does not exist under %s: %s', cycle_path, source_entry_name);
    end
    source_entry_name = string(source_entry_name);
    return;
end

folder_candidates = dir(fullfile(cycle_path, 'Cam*'));
folder_candidates = folder_candidates(~ismember({folder_candidates.name}, {'.', '..'}));
if ~isempty(folder_candidates)
    folder_candidates = folder_candidates([folder_candidates.isdir]);
    folder_names = string({folder_candidates.name});
    folder_candidates = folder_candidates(~endsWith(folder_names, "_Analysis", 'IgnoreCase', true));
end
if numel(folder_candidates) == 1
    source_entry_name = string(folder_candidates(1).name);
    source_entry_path = fullfile(folder_candidates(1).folder, folder_candidates(1).name);
    return;
elseif numel(folder_candidates) > 1
    names = strjoin(string({folder_candidates.name}), ', ');
    error('Multiple Cam* folders found under %s. Set source_entry_name explicitly. Candidates: %s', cycle_path, names);
end

file_candidates = [ ...
    dir(fullfile(cycle_path, '*.tif')); ...
    dir(fullfile(cycle_path, '*.tiff')); ...
    dir(fullfile(cycle_path, '*.btf')); ...
    dir(fullfile(cycle_path, '*.bin'))];
if numel(file_candidates) == 1
    source_entry_name = string(file_candidates(1).name);
    source_entry_path = fullfile(file_candidates(1).folder, file_candidates(1).name);
    return;
elseif numel(file_candidates) > 1
    names = strjoin(string({file_candidates.name}), ', ');
    error('Multiple movie files found under %s. Set source_entry_name explicitly. Candidates: %s', cycle_path, names);
end

error('No usable Cam* folder or movie file was found under cycle path: %s', cycle_path);
end

function cycle_dirs = sort_cycle_dirs_ap(cycle_dirs)
cycle_names = string({cycle_dirs.name});
cycle_numbers = NaN(size(cycle_names));
for idx = 1:numel(cycle_names)
    tokens = regexp(cycle_names(idx), 'Cycle(\d+)', 'tokens', 'once');
    if ~isempty(tokens)
        cycle_numbers(idx) = str2double(tokens{1});
    end
end
[~, order] = sortrows([isnan(cycle_numbers(:)), cycle_numbers(:), (1:numel(cycle_names))']);
cycle_dirs = cycle_dirs(order);
end

function explicit_result = find_latest_explicit_result_ap(cycle_path)
explicit_result = "";
listing = dir(fullfile(cycle_path, '**', '-1_explicit_results.mat'));
if isempty(listing)
    return;
end
[~, newest_idx] = max([listing.datenum]);
explicit_result = string(fullfile(listing(newest_idx).folder, listing(newest_idx).name));
end

function roi_file = find_latest_roi_file_ap(cycle_path)
roi_file = "";
listing = dir(fullfile(cycle_path, '**', '1_raw_ROI.mat'));
if isempty(listing)
    return;
end
[~, newest_idx] = max([listing.datenum]);
roi_file = string(fullfile(listing(newest_idx).folder, listing(newest_idx).name));
end

function output = build_record_average_summary_ap(rec_path, batch_results, reference_roi_file, trace_polarity, average_output_dir_name)
output_dir = fullfile(rec_path, average_output_dir_name);
if ~isfolder(output_dir)
    mkdir(output_dir);
end

[result_dirs, cycle_names] = resolve_successful_result_dirs_ap(batch_results);
if isempty(result_dirs)
    error('No successful cycle result folders are available for record-level averaging.');
end

fprintf('Using %d cycle result folders for record-average analysis.\n', numel(result_dirs));
cycle_entry_list = cell(numel(result_dirs), 1);
for idx = 1:numel(result_dirs)
    fprintf('  Loading %s\n', result_dirs(idx));
    cycle_entry_list{idx} = load_cycle_average_entry_ap(result_dirs(idx), cycle_names(idx));
end
cycle_entries = vertcat(cycle_entry_list{:});

record_average = struct();
record_average.info = struct( ...
    'record_path', string(rec_path), ...
    'reference_roi_file', string(reference_roi_file), ...
    'cycle_names', cycle_names(:), ...
    'result_dirs', result_dirs(:), ...
    'cycle_count', numel(cycle_entries), ...
    'created_at', datetime("now"), ...
    'alignment_mode', "flash_onset_if_available", ...
    'trace_polarity', trace_polarity);

record_average.raw = average_cycle_stage_ap(cycle_entries, 'raw', 'flash_onset_if_available');
record_average.sensitivity = average_cycle_stage_ap(cycle_entries, 'sensitivity', 'flash_onset_if_available');
record_average.snr = average_cycle_stage_ap(cycle_entries, 'snr', 'flash_onset_if_available');
record_average.stim_windows = resolve_record_average_stim_windows_ap(cycle_entries, record_average.raw.time);

raw_summary_png = fullfile(output_dir, '1_record_average_raw_summary.png');
sens_summary_png = fullfile(output_dir, '2_record_average_sensitivity_summary.png');
snr_summary_png = fullfile(output_dir, '3_record_average_snr_summary.png');

plot_record_stage_summary_ap(record_average.raw, record_average.stim_windows, 'Raw', raw_summary_png, [0.85 0.15 0.15], 1, 'Raw');
plot_record_stage_summary_ap(record_average.sensitivity, record_average.stim_windows, 'Sensitivity', sens_summary_png, [0.85 0.15 0.15], trace_polarity, 'Sensitivity');
plot_record_stage_summary_ap(record_average.snr, record_average.stim_windows, 'SNR', snr_summary_png, [0.85 0.15 0.15], trace_polarity, 'SNR');

tf_input = struct( ...
    'trace', trace_polarity * record_average.sensitivity.average, ...
    'time', record_average.sensitivity.time, ...
    'frame_rate', record_average.sensitivity.frame_rate, ...
    'stage_name', "sensitivity_average");
record_average.time_frequency = build_record_average_time_frequency_ap(tf_input, record_average.stim_windows, output_dir);

record_average.visualizations = struct( ...
    'raw_summary_png', string(raw_summary_png), ...
    'sensitivity_summary_png', string(sens_summary_png), ...
    'snr_summary_png', string(snr_summary_png), ...
    'fourier_summary_png', string(record_average.time_frequency.visualizations.fourier_summary_png), ...
    'wavelet_summary_png', string(record_average.time_frequency.visualizations.wavelet_summary_png));

results_file = fullfile(output_dir, 'record_average_results.mat');
save(results_file, 'record_average', '-v7.3');

output = struct( ...
    'output_dir', string(output_dir), ...
    'results_file', string(results_file));
end

function [result_dirs, cycle_names] = resolve_successful_result_dirs_ap(batch_results)
result_dirs = strings(0, 1);
cycle_names = strings(0, 1);
for idx = 1:numel(batch_results)
    status_value = string(batch_results(idx).status);
    result_dir = string(batch_results(idx).save_path);
    if strlength(result_dir) == 0 || ~isfolder(result_dir)
        continue;
    end
    if status_value == "completed" || status_value == "skipped_existing_result" || status_value == "existing_result_detected"
        result_dirs(end+1, 1) = result_dir; %#ok<AGROW>
        cycle_names(end+1, 1) = string(batch_results(idx).cycle_name); %#ok<AGROW>
    end
end
end

function entry = load_cycle_average_entry_ap(result_dir, cycle_name)
movie_info_path = fullfile(result_dir, 'movie_info.mat');
trace_results_path = fullfile(result_dir, 'trace_results.mat');
stim_path = fullfile(result_dir, 'stim_results.mat');

tmp_movie = load(movie_info_path, 'movie_info');
tmp_trace = load(trace_results_path, 'trace_results');
movie_info = tmp_movie.movie_info;
trace_results = tmp_trace.trace_results;

entry = struct();
entry.cycle_name = string(cycle_name);
entry.result_dir = string(result_dir);
entry.raw = extract_cycle_stage_ap(trace_results, movie_info, {'raw'});
entry.sensitivity = extract_cycle_stage_ap(trace_results, movie_info, {'sensitivity'});
entry.snr = extract_cycle_stage_ap(trace_results, movie_info, {'snr'});

entry.stim_windows = struct();
if isfile(stim_path)
    tmp_stim = load(stim_path, 'stim_results');
    if isfield(tmp_stim, 'stim_results') && isstruct(tmp_stim.stim_results) ...
            && isfield(tmp_stim.stim_results, 'windows')
        entry.stim_windows = tmp_stim.stim_results.windows;
    end
end
entry.alignment = struct( ...
    'flash_onset_time', resolve_flash_onset_time_ap(entry.stim_windows));
end

function stage = extract_cycle_stage_ap(trace_results, movie_info, preferred_stages)
[stage_data, stage_name] = resolve_preferred_stage_from_results_ap(trace_results, preferred_stages);
frame_rate = double(movie_info.frame_rate);
nframes = size(stage_data, 1);
stage = struct( ...
    'stage_name', string(stage_name), ...
    'data', double(stage_data), ...
    'frame_rate', frame_rate, ...
    'time', (1:nframes)' / frame_rate);
end

function [stage_data, stage_name] = resolve_preferred_stage_from_results_ap(trace_results, preferred_stages)
for idx = 1:numel(preferred_stages)
    stage_name = preferred_stages{idx};
    if isfield(trace_results, 'trace_results') ...
            && isfield(trace_results.trace_results, stage_name) ...
            && isfield(trace_results.trace_results.(stage_name), 'data')
        stage_data = trace_results.trace_results.(stage_name).data;
        return;
    end
end
error('Required stage is missing. Tried: %s', strjoin(preferred_stages, ', '));
end

function averaged = average_cycle_stage_ap(cycle_entries, stage_name, alignment_mode)
sample = cycle_entries(1).(stage_name).data;
nrois = size(sample, 2);
frame_rate = cycle_entries(1).(stage_name).frame_rate;
stage_labels = strings(numel(cycle_entries), 1);
relative_starts = NaN(numel(cycle_entries), 1);
relative_ends = NaN(numel(cycle_entries), 1);

for idx = 1:numel(cycle_entries)
    current_stage = cycle_entries(idx).(stage_name);
    if size(current_stage.data, 2) ~= nrois
        error('ROI count mismatch while averaging %s across cycles.', stage_name);
    end
    if abs(current_stage.frame_rate - frame_rate) > 1e-9
        error('Frame rate mismatch while averaging %s across cycles.', stage_name);
    end
    stage_labels(idx) = current_stage.stage_name;
    alignment_time = resolve_cycle_alignment_time_ap(cycle_entries(idx), alignment_mode);
    relative_time = current_stage.time(:) - alignment_time;
    relative_starts(idx) = relative_time(1);
    relative_ends(idx) = relative_time(end);
end

common_start = max(relative_starts);
common_end = min(relative_ends);
if ~isfinite(common_start) || ~isfinite(common_end) || common_end <= common_start
    error('No overlapping time range remains after alignment for %s.', stage_name);
end

dt = 1 / frame_rate;
common_time = (common_start:dt:common_end)';
if numel(common_time) < 2
    common_time = [common_start; common_end];
end

per_cycle = NaN(numel(common_time), nrois, numel(cycle_entries));
alignment_time_by_cycle = NaN(numel(cycle_entries), 1);
for idx = 1:numel(cycle_entries)
    current_data = double(cycle_entries(idx).(stage_name).data);
    current_time = cycle_entries(idx).(stage_name).time(:);
    alignment_time = resolve_cycle_alignment_time_ap(cycle_entries(idx), alignment_mode);
    alignment_time_by_cycle(idx) = alignment_time;
    relative_time = current_time - alignment_time;
    for roi_idx = 1:nrois
        per_cycle(:, roi_idx, idx) = interp1(relative_time, current_data(:, roi_idx), common_time, 'linear', NaN);
    end
end

averaged = struct( ...
    'stage_name', string(stage_labels(1)), ...
    'stage_name_by_cycle', stage_labels, ...
    'frame_rate', frame_rate, ...
    'time', common_time, ...
    'per_cycle', per_cycle, ...
    'average', mean(per_cycle, 3, 'omitnan'), ...
    'cycle_names', string({cycle_entries.cycle_name})', ...
    'alignment_mode', string(alignment_mode), ...
    'alignment_time_by_cycle', alignment_time_by_cycle, ...
    'ncycles', numel(cycle_entries), ...
    'nrois', nrois);
end

function alignment_time = resolve_cycle_alignment_time_ap(cycle_entry, alignment_mode)
alignment_time = 0;
if strcmpi(string(alignment_mode), "flash_onset_if_available")
    onset_time = double(cycle_entry.alignment.flash_onset_time);
    if isfinite(onset_time)
        alignment_time = onset_time;
    end
end
end

function onset_time = resolve_flash_onset_time_ap(stim_windows)
onset_time = NaN;
if ~isstruct(stim_windows) || ~isfield(stim_windows, 'stim_type') ...
        || ~strcmpi(string(stim_windows.stim_type), "visualstim_flash") ...
        || ~isfield(stim_windows, 'channel')
    return;
end

channel_windows = stim_windows.channel;
if isfield(channel_windows, 'stim_time_ranges') && ~isempty(channel_windows.stim_time_ranges)
    onset_time = double(channel_windows.stim_time_ranges(1, 1));
elseif isfield(channel_windows, 'stim_frames') && ~isempty(channel_windows.stim_frames)
    onset_time = double(channel_windows.stim_frames(1, 1)) / 400;
end
end

function stim_windows = resolve_record_average_stim_windows_ap(cycle_entries, time_axis)
stim_windows = struct();
for idx = 1:numel(cycle_entries)
    candidate = cycle_entries(idx).stim_windows;
    if isstruct(candidate) && isfield(candidate, 'supported') && candidate.supported
        stim_windows = candidate;
        anchor = cycle_entries(idx).alignment.flash_onset_time;
        break;
    end
end
if isempty(fieldnames(stim_windows))
    return;
end

if isfield(stim_windows, 'channel')
    stim_windows.channel = align_channel_windows_to_flash_ap(stim_windows.channel, anchor, time_axis);
    if isfield(stim_windows, 'flash_windows')
        stim_windows.channel.flash_windows = stim_windows.flash_windows;
    end
end
end

function channel_windows = align_channel_windows_to_flash_ap(channel_windows, anchor_time, time_axis)
if ~isfinite(anchor_time)
    return;
end
frame_rate = 1 / median(diff(time_axis), 'omitnan');
if ~isfinite(frame_rate) || frame_rate <= 0
    frame_rate = 400;
end
fields_to_shift = { ...
    'stim_time_ranges', 'block_time_ranges', 'shading_time_ranges'};
for idx = 1:numel(fields_to_shift)
    field_name = fields_to_shift{idx};
    if isfield(channel_windows, field_name) && ~isempty(channel_windows.(field_name))
        channel_windows.(field_name) = channel_windows.(field_name) - anchor_time;
    end
end
if isfield(channel_windows, 'stim_frames') && ~isempty(channel_windows.stim_frames)
    channel_windows.stim_frames = convert_relative_times_to_frames_ap(channel_windows.stim_time_ranges, time_axis);
end
if isfield(channel_windows, 'baseline_frames') && ~isempty(channel_windows.baseline_frames)
    baseline_time_ranges = channel_windows.baseline_frames / frame_rate;
    baseline_time_ranges = baseline_time_ranges - anchor_time;
    channel_windows.baseline_frames = convert_relative_times_to_frames_ap(baseline_time_ranges, time_axis);
end
if isfield(channel_windows, 'block_frames') && isfield(channel_windows, 'block_time_ranges') && ~isempty(channel_windows.block_time_ranges)
    channel_windows.block_frames = convert_relative_times_to_frames_ap(channel_windows.block_time_ranges, time_axis);
end
if isfield(channel_windows, 'shading_frames') && isfield(channel_windows, 'shading_time_ranges') && ~isempty(channel_windows.shading_time_ranges)
    channel_windows.shading_frames = convert_relative_times_to_frames_ap(channel_windows.shading_time_ranges, time_axis);
end
end

function frame_ranges = convert_relative_times_to_frames_ap(time_ranges, time_axis)
frame_ranges = NaN(size(time_ranges));
for idx = 1:size(time_ranges, 1)
    [~, start_idx] = min(abs(time_axis - time_ranges(idx, 1)));
    [~, end_idx] = min(abs(time_axis - time_ranges(idx, 2)));
    frame_ranges(idx, :) = [start_idx, end_idx];
end
end

function plot_record_stage_summary_ap(stage, stim_windows, metric_name, output_png, trace_color, display_scale, title_text)
fig = figure('Color', 'w', 'Name', sprintf('Record Average %s Summary', metric_name));
ax = axes(fig);
hold(ax, 'on');

average_trace = double(display_scale) * double(stage.average);
per_cycle = double(display_scale) * double(stage.per_cycle);
spacing = median(range(average_trace, 1, 'omitnan'));
if ~isfinite(spacing) || spacing <= 0
    spacing = 1;
end
spacing = spacing * 1.8;

for roi_idx = 1:size(average_trace, 2)
    offset = (size(average_trace, 2) - roi_idx) * spacing;
    for cycle_idx = 1:size(per_cycle, 3)
        light_color = trace_color * 0.28 + [1 1 1] * 0.72;
        plot(ax, stage.time, per_cycle(:, roi_idx, cycle_idx) + offset, ...
            'Color', light_color, 'LineWidth', 0.6);
    end
    plot(ax, stage.time, average_trace(:, roi_idx) + offset, 'Color', trace_color, 'LineWidth', 1.6);
end

if isstruct(stim_windows) && isfield(stim_windows, 'channel')
    add_stim_shading_rec_ap(ax, stim_windows.channel, 0.12);
end

xlabel(ax, 'Time (s)');
ylabel(ax, 'Stacked ROI trace');
title(ax, sprintf('Record-Average %s | %s', metric_name, title_text));
grid(ax, 'on');
save_figure_bundle_rec_ap(fig, replace(output_png, '.png', '.fig'), output_png);
close(fig);
end

function add_stim_shading_rec_ap(ax, channel_windows, alpha_value)
if ~isstruct(channel_windows)
    return;
end
time_ranges = [];
labels = strings(0, 1);
if isfield(channel_windows, 'shading_time_ranges') && ~isempty(channel_windows.shading_time_ranges)
    time_ranges = channel_windows.shading_time_ranges;
    if isfield(channel_windows, 'shading_labels') && numel(channel_windows.shading_labels) == size(time_ranges, 1)
        labels = string(channel_windows.shading_labels(:));
    end
elseif isfield(channel_windows, 'block_time_ranges') && ~isempty(channel_windows.block_time_ranges)
    time_ranges = channel_windows.block_time_ranges;
elseif isfield(channel_windows, 'stim_time_ranges') && ~isempty(channel_windows.stim_time_ranges)
    time_ranges = channel_windows.stim_time_ranges;
end
if isempty(time_ranges)
    return;
end
if isempty(labels)
    labels = repmat("stim", size(time_ranges, 1), 1);
end

yl = ylim(ax);
hold(ax, 'on');
for idx = 1:size(time_ranges, 1)
    [shade_color, shade_alpha] = resolve_stim_shading_style_rec_ap(labels(idx), [0.7 0.7 0.7], alpha_value);
    p = patch(ax, ...
        [time_ranges(idx, 1), time_ranges(idx, 2), time_ranges(idx, 2), time_ranges(idx, 1)], ...
        [yl(1), yl(1), yl(2), yl(2)], ...
        shade_color, 'FaceAlpha', shade_alpha, 'EdgeColor', 'none');
    uistack(p, 'bottom');
end
ylim(ax, yl);
end

function [shade_color, shade_alpha] = resolve_stim_shading_style_rec_ap(label, default_color, default_alpha)
label = lower(char(string(label)));
shade_color = default_color;
shade_alpha = default_alpha;
if contains(label, 'flash')
    shade_color = [1.00 0.82 0.20];
    shade_alpha = max(default_alpha, 0.26);
elseif contains(label, 'black')
    shade_color = [0.05 0.05 0.05];
    shade_alpha = max(default_alpha, 0.18);
elseif contains(label, 'white')
    shade_color = [0.85 0.85 0.85];
    shade_alpha = max(default_alpha, 0.10);
elseif contains(label, 'gray')
    shade_color = [0.55 0.55 0.55];
    shade_alpha = max(default_alpha, 0.14);
elseif contains(label, 'blue')
    shade_color = [0.35 0.55 0.95];
    shade_alpha = max(default_alpha, 0.14);
end
end

function time_frequency = build_record_average_time_frequency_ap(tf_input, stim_windows, output_dir)
params = struct( ...
    'trace_stage', string(tf_input.stage_name), ...
    'min_freq_hz', 0.5, ...
    'max_freq_hz', min(80, tf_input.frame_rate / 2 - eps), ...
    'wavelet_name', "amor");

channel_tf = analyze_population_time_frequency_rec_ap( ...
    tf_input.trace, tf_input.frame_rate, tf_input.time, params, 'Voltage');

[fourier_fig, fourier_png] = plot_population_fourier_summary_rec_ap( ...
    channel_tf, stim_windows, params, output_dir);
[wavelet_fig, wavelet_png] = plot_population_wavelet_summary_rec_ap( ...
    channel_tf, stim_windows, output_dir);

time_frequency = struct( ...
    'parameters', params, ...
    'channel', channel_tf, ...
    'visualizations', struct( ...
        'fourier_summary_fig', string(fourier_fig), ...
        'fourier_summary_png', string(fourier_png), ...
        'wavelet_summary_fig', string(wavelet_fig), ...
        'wavelet_summary_png', string(wavelet_png)));
end

function result = analyze_population_time_frequency_rec_ap(trace_matrix, frame_rate, t_axis, params, channel_name)
trace_matrix = double(trace_matrix);
t_axis = double(t_axis(:));
nrois = size(trace_matrix, 2);
spectrum_frequency = [];
spectrum_amplitude = [];
peak_frequency_hz = NaN(nrois, 1);
peak_amplitude = NaN(nrois, 1);
signal_rms = NaN(nrois, 1);

for roi_idx = 1:nrois
    x = sanitize_trace_for_spectrum_rec_ap(trace_matrix(:, roi_idx));
    signal_rms(roi_idx) = rms(x);
    [frequency_i, amplitude_i] = compute_trace_fft_spectrum_rec_ap(x, frame_rate);
    if isempty(spectrum_frequency)
        spectrum_frequency = frequency_i;
        spectrum_amplitude = NaN(numel(frequency_i), nrois);
    end
    spectrum_amplitude(:, roi_idx) = amplitude_i;
    valid_mask = frequency_i >= params.min_freq_hz & frequency_i <= params.max_freq_hz;
    if any(valid_mask)
        [peak_amplitude(roi_idx), max_idx] = max(amplitude_i(valid_mask));
        valid_frequency = frequency_i(valid_mask);
        peak_frequency_hz(roi_idx) = valid_frequency(max_idx);
    end
end

representative_roi = find(peak_amplitude == max(peak_amplitude, [], 'omitnan'), 1, 'first');
if isempty(representative_roi)
    representative_roi = 1;
end

result = struct( ...
    'channel_name', string(channel_name), ...
    'frame_rate', frame_rate, ...
    'time', t_axis, ...
    'trace_matrix', trace_matrix, ...
    'analysis_parameters', params, ...
    'fft', struct('frequency', spectrum_frequency, 'amplitude', spectrum_amplitude), ...
    'roi', struct( ...
        'peak_frequency_hz', peak_frequency_hz, ...
        'peak_amplitude', peak_amplitude, ...
        'signal_rms', signal_rms), ...
    'representative', struct( ...
        'roi_index', representative_roi, ...
        'trace', sanitize_trace_for_spectrum_rec_ap(trace_matrix(:, representative_roi)), ...
        'peak_frequency_hz', peak_frequency_hz(representative_roi), ...
        'wavelet', compute_wavelet_scalogram_rec_ap(trace_matrix(:, representative_roi), frame_rate, t_axis, params)), ...
    'summary', struct( ...
        'nrois', nrois, ...
        'median_peak_frequency_hz', median(peak_frequency_hz, 'omitnan'), ...
        'median_peak_amplitude', median(peak_amplitude, 'omitnan')));
end

function x = sanitize_trace_for_spectrum_rec_ap(x)
x = double(x(:));
if isempty(x) || all(~isfinite(x))
    x = zeros(size(x));
    return;
end
x = fillmissing(x, 'linear', 'EndValues', 'nearest');
x = x - mean(x, 'omitnan');
end

function [frequency, amplitude] = compute_trace_fft_spectrum_rec_ap(x, frame_rate)
x = sanitize_trace_for_spectrum_rec_ap(x);
nsamples = numel(x);
if nsamples < 8
    frequency = zeros(0, 1);
    amplitude = zeros(0, 1);
    return;
end
nfft = 2 ^ nextpow2(nsamples);
y = fft(x, nfft);
p2 = abs(y / nsamples);
amplitude = p2(1:nfft / 2 + 1);
if numel(amplitude) > 2
    amplitude(2:end-1) = 2 * amplitude(2:end-1);
end
frequency = frame_rate * (0:(nfft / 2))' / nfft;
end

function wavelet_data = compute_wavelet_scalogram_rec_ap(x, frame_rate, t_axis, params)
x = sanitize_trace_for_spectrum_rec_ap(x);
[wt, frequency] = cwt(x, frame_rate);
valid_mask = frequency >= params.min_freq_hz & frequency <= params.max_freq_hz;
wavelet_data = struct( ...
    'time', double(t_axis(:)), ...
    'frequency', frequency(valid_mask), ...
    'power', abs(wt(valid_mask, :)).^2);
end

function [fig_file, png_file] = plot_population_fourier_summary_rec_ap(channel_tf, stim_windows, params, output_dir)
fig = figure('Color', 'w', 'Name', 'Record-Average FFT Summary');
tiledlayout(fig, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
nexttile;
plot_channel_fft_summary_rec_ap(channel_tf, 'r', 'FFT Amplitude', [0, params.max_freq_hz], false);
nexttile;
plot_channel_fft_summary_rec_ap(channel_tf, 'r', 'FFT Amplitude', [0, 10], true);
nexttile;
plot_flash_fft_comparison_rec_ap(channel_tf, resolve_optional_channel_windows_ap(stim_windows), 'r', 'Flash FFT');
nexttile;
plot_roi_peak_summary_rec_ap(channel_tf, [0, 10], 'ROI Peak Frequency (0-10 Hz)');

fig_file = fullfile(output_dir, '4_record_average_fourier_summary.fig');
png_file = fullfile(output_dir, '4_record_average_fourier_summary.png');
save_figure_bundle_rec_ap(fig, fig_file, png_file);
close(fig);
end

function [fig_file, png_file] = plot_population_wavelet_summary_rec_ap(channel_tf, stim_windows, output_dir)
fig = figure('Color', 'w', 'Name', 'Record-Average Wavelet Summary');
nrois = size(channel_tf.trace_matrix, 2);
tiledlayout(fig, nrois, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
for roi_idx = 1:nrois
    nexttile;
    plot_roi_trace_with_stim_rec_ap(gca, channel_tf, roi_idx, resolve_optional_channel_windows_ap(stim_windows), 'r', 'Trace');
    nexttile;
    plot_roi_wavelet_scalogram_rec_ap(gca, channel_tf, roi_idx, resolve_optional_channel_windows_ap(stim_windows), 'Wavelet');
end

fig_file = fullfile(output_dir, '5_record_average_wavelet_summary.fig');
png_file = fullfile(output_dir, '5_record_average_wavelet_summary.png');
save_figure_bundle_rec_ap(fig, fig_file, png_file);
close(fig);
end

function plot_channel_fft_summary_rec_ap(channel_result, trace_color, title_text, display_band_hz, annotate_peak)
frequency = channel_result.fft.frequency;
amplitude = channel_result.fft.amplitude;
display_mask = frequency >= display_band_hz(1) & frequency <= display_band_hz(2);
if isempty(frequency) || isempty(amplitude) || ~any(display_mask)
    text(0.5, 0.5, 'FFT unavailable', 'Units', 'normalized', 'HorizontalAlignment', 'center');
    axis off;
    return;
end
plot(frequency(display_mask), amplitude(display_mask, :), 'Color', [0.82 0.82 0.82], 'LineWidth', 0.8);
hold on;
mean_amplitude = mean(amplitude, 2, 'omitnan');
plot(frequency(display_mask), mean_amplitude(display_mask), trace_color, 'LineWidth', 2);
if annotate_peak
    [peak_hz, peak_amp] = find_band_peak_from_amplitude_rec_ap(frequency, mean_amplitude, display_band_hz);
    if isfinite(peak_hz)
        xline(peak_hz, '--', sprintf('%.2f Hz', peak_hz), 'Color', trace_color, 'LineWidth', 1.2);
        scatter(peak_hz, peak_amp, 42, trace_color, 'filled');
    end
end
xlabel('Frequency (Hz)');
ylabel('Single-Sided Amplitude');
title(title_text);
xlim(display_band_hz);
grid on;
end

function plot_flash_fft_comparison_rec_ap(channel_result, channel_windows, trace_color, title_text)
[baseline_trace, stim_trace, supported] = extract_flash_average_segments_rec_ap(channel_result, channel_windows);
if ~supported
    text(0.5, 0.5, 'Flash pre/post FFT unavailable', 'Units', 'normalized', 'HorizontalAlignment', 'center');
    axis off;
    return;
end
[f_base, a_base] = compute_trace_fft_spectrum_rec_ap(baseline_trace, channel_result.frame_rate);
[f_stim, a_stim] = compute_trace_fft_spectrum_rec_ap(stim_trace, channel_result.frame_rate);
base_mask = f_base >= 0 & f_base <= 10;
stim_mask = f_stim >= 0 & f_stim <= 10;
plot(f_base(base_mask), a_base(base_mask), 'Color', [0.55 0.55 0.55], 'LineWidth', 1.5);
hold on;
plot(f_stim(stim_mask), a_stim(stim_mask), 'Color', trace_color, 'LineWidth', 2);
xlabel('Frequency (Hz)');
ylabel('Single-Sided Amplitude');
title(sprintf('%s | baseline vs response', title_text));
legend({'Baseline', 'Response'}, 'Location', 'best');
xlim([0, 10]);
grid on;
end

function [baseline_trace, stim_trace, supported] = extract_flash_average_segments_rec_ap(channel_result, channel_windows)
baseline_trace = [];
stim_trace = [];
supported = false;
if ~isstruct(channel_windows) || ~isstruct(channel_result) ...
        || ~isfield(channel_result, 'time') || ~isfield(channel_result, 'representative')
    return;
end
baseline_window_sec = 2.0;
response_window_sec = 2.0;
if isfield(channel_windows, 'flash_windows') && isstruct(channel_windows.flash_windows)
    if isfield(channel_windows.flash_windows, 'baseline_window_sec') && ~isempty(channel_windows.flash_windows.baseline_window_sec)
        baseline_window_sec = double(channel_windows.flash_windows.baseline_window_sec);
    end
    if isfield(channel_windows.flash_windows, 'response_window_sec') && ~isempty(channel_windows.flash_windows.response_window_sec)
        response_window_sec = double(channel_windows.flash_windows.response_window_sec);
    end
end

t_axis = double(channel_result.time(:));
trace = double(channel_result.representative.trace(:));
baseline_mask = t_axis >= -baseline_window_sec & t_axis < 0;
stim_mask = t_axis >= 0 & t_axis <= response_window_sec;
if ~any(baseline_mask) || ~any(stim_mask)
    return;
end
baseline_trace = trace(baseline_mask);
stim_trace = trace(stim_mask);
target_length = min(numel(baseline_trace), numel(stim_trace));
if target_length < 8
    baseline_trace = [];
    stim_trace = [];
    return;
end
baseline_trace = baseline_trace(end - target_length + 1:end);
stim_trace = stim_trace(1:target_length);
supported = true;
end

function plot_roi_peak_summary_rec_ap(channel_result, display_band_hz, title_text)
roi_count = size(channel_result.fft.amplitude, 2);
peak_frequency_hz = NaN(roi_count, 1);
peak_amplitude = NaN(roi_count, 1);
for roi_idx = 1:roi_count
    [peak_frequency_hz(roi_idx), peak_amplitude(roi_idx)] = find_band_peak_from_amplitude_rec_ap( ...
        channel_result.fft.frequency, channel_result.fft.amplitude(:, roi_idx), display_band_hz);
end
valid = isfinite(peak_frequency_hz) & isfinite(peak_amplitude);
if ~any(valid)
    text(0.5, 0.5, 'ROI peak summary unavailable', 'Units', 'normalized', 'HorizontalAlignment', 'center');
    axis off;
    return;
end
scatter(find(valid), peak_frequency_hz(valid), 42, peak_amplitude(valid), 'filled', 'MarkerFaceAlpha', 0.8);
hold on;
plot(find(valid), peak_frequency_hz(valid), '-', 'Color', [0.75 0.75 0.75]);
xlabel('ROI');
ylabel('Peak Frequency (Hz)');
ylim(display_band_hz);
title(title_text);
cb = colorbar;
cb.Label.String = 'Peak Amplitude';
grid on;
end

function plot_roi_trace_with_stim_rec_ap(ax, channel_result, roi_idx, channel_windows, trace_color, title_text)
t = channel_result.time(:);
x = double(channel_result.trace_matrix(:, roi_idx));
plot(ax, t, x, 'Color', trace_color, 'LineWidth', 1.0);
apply_trace_axis_limits_rec_ap(ax, x);
if ~isempty(fieldnames(channel_windows))
    add_stim_shading_rec_ap(ax, channel_windows, 0.16);
    hold(ax, 'on');
end
plot(ax, t, x, 'Color', trace_color, 'LineWidth', 1.0);
xlabel(ax, 'Time (s)');
ylabel(ax, resolve_time_frequency_trace_label_rec_ap(channel_result));
title(ax, sprintf('%s | ROI %d | peak %.2f Hz', title_text, roi_idx, channel_result.roi.peak_frequency_hz(roi_idx)));
grid(ax, 'on');
end

function plot_roi_wavelet_scalogram_rec_ap(ax, channel_result, roi_idx, channel_windows, title_text)
trace = double(channel_result.trace_matrix(:, roi_idx));
params = channel_result.analysis_parameters;
spec = compute_wavelet_scalogram_rec_ap(trace, channel_result.frame_rate, channel_result.time, params);
if isempty(spec.frequency) || isempty(spec.power)
    text(ax, 0.5, 0.5, 'Wavelet map unavailable', 'Units', 'normalized', 'HorizontalAlignment', 'center');
    axis(ax, 'off');
    return;
end
[frequency, sort_idx] = sort(spec.frequency(:), 'ascend');
power = spec.power(sort_idx, :);
surface(ax, ...
    repmat(spec.time(:)', numel(frequency), 1), ...
    repmat(frequency, 1, numel(spec.time)), ...
    zeros(size(power)), ...
    power, ...
    'EdgeColor', 'none');
view(ax, 2);
set(ax, 'YDir', 'normal');
set(ax, 'YScale', 'log');
xlabel(ax, 'Time (s)');
ylabel(ax, 'Frequency (Hz)');
title(ax, sprintf('%s | ROI %d | peak %.2f Hz', title_text, roi_idx, channel_result.roi.peak_frequency_hz(roi_idx)));
cb = colorbar(ax);
cb.Label.String = 'Wavelet Power';
hold(ax, 'on');
if isfinite(channel_result.roi.peak_frequency_hz(roi_idx))
    yline(ax, channel_result.roi.peak_frequency_hz(roi_idx), 'w--', 'LineWidth', 1.0);
end
ylim(ax, [max(min(frequency(frequency > 0)), 0.5), max(frequency)]);
yticks(ax, [1 3 10 30 80]);
overlay_stim_boundaries_rec_ap(ax, channel_windows);
end

function y_label = resolve_time_frequency_trace_label_rec_ap(channel_result)
stage_name = "";
if isfield(channel_result, 'analysis_parameters') && isstruct(channel_result.analysis_parameters) ...
        && isfield(channel_result.analysis_parameters, 'trace_stage')
    stage_name = string(channel_result.analysis_parameters.trace_stage);
end
stage_name = lower(stage_name);
switch stage_name
    case {"snr", "snr_average"}
        y_label = 'SNR (signed)';
    case {"sensitivity", "sensitivity_average"}
        y_label = 'Sensitivity (signed)';
    case {"bleach_removed", "bleach_removed_average"}
        y_label = 'Bleach-Removed Signal';
    otherwise
        y_label = 'Trace Value';
end
end

function apply_trace_axis_limits_rec_ap(ax, x)
x = double(x(:));
valid = isfinite(x);
if ~any(valid)
    ylim(ax, [-1, 1]);
    return;
end
ymin = min(x(valid));
ymax = max(x(valid));
if ymax <= ymin
    pad = max(1e-3, abs(ymax) * 0.1 + 1e-3);
else
    pad = 0.08 * (ymax - ymin);
end
ylim(ax, [ymin - pad, ymax + pad]);
end

function overlay_stim_boundaries_rec_ap(ax, channel_windows)
if ~isstruct(channel_windows)
    return;
end
time_ranges = [];
if isfield(channel_windows, 'shading_time_ranges') && ~isempty(channel_windows.shading_time_ranges)
    time_ranges = channel_windows.shading_time_ranges;
elseif isfield(channel_windows, 'block_time_ranges') && ~isempty(channel_windows.block_time_ranges)
    time_ranges = channel_windows.block_time_ranges;
elseif isfield(channel_windows, 'stim_time_ranges') && ~isempty(channel_windows.stim_time_ranges)
    time_ranges = channel_windows.stim_time_ranges;
end
for idx = 1:size(time_ranges, 1)
    xline(ax, time_ranges(idx, 1), 'k:', 'LineWidth', 0.8);
    xline(ax, time_ranges(idx, 2), 'k:', 'LineWidth', 0.8);
end
end

function [peak_frequency_hz, peak_amplitude] = find_band_peak_from_amplitude_rec_ap(frequency, amplitude, display_band_hz)
peak_frequency_hz = NaN;
peak_amplitude = NaN;
if isempty(frequency) || isempty(amplitude)
    return;
end
band_mask = frequency >= display_band_hz(1) & frequency <= display_band_hz(2);
if ~any(band_mask)
    return;
end
band_frequency = frequency(band_mask);
band_amplitude = double(amplitude(band_mask));
if all(~isfinite(band_amplitude))
    return;
end
[peak_amplitude, max_idx] = max(band_amplitude);
peak_frequency_hz = band_frequency(max_idx);
end

function channel_windows = resolve_optional_channel_windows_ap(stim_windows)
channel_windows = struct();
if isstruct(stim_windows) && isfield(stim_windows, 'channel')
    channel_windows = stim_windows.channel;
end
end

function save_figure_bundle_rec_ap(fig_handle, fig_file, png_file)
prepare_figure_for_png_export_rec_ap(fig_handle);
saveas(fig_handle, fig_file, 'fig');
exportgraphics(fig_handle, png_file, 'Resolution', 150);
end

function prepare_figure_for_png_export_rec_ap(fig_handle)
if ~isgraphics(fig_handle, 'figure')
    return;
end
set(fig_handle, 'Units', 'pixels');
try
    set(fig_handle, 'WindowState', 'maximized');
catch
    screen_size = get(groot, 'ScreenSize');
    if isnumeric(screen_size) && numel(screen_size) >= 4
        set(fig_handle, 'Position', screen_size);
    end
end
drawnow;
end
