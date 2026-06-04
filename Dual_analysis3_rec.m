% Dual_analysis3_rec - Run Dual_analysis3 across all cycles in one record
%
% Workflow:
% 1. Choose one reference cycle.
% 2. Create or reuse one ROI file from that reference cycle.
% 3. Reuse the same ROI file for every other Cycle* under the same record.
%
% Analysis-only fork workflow:
%   Set analysis_mode = 'analysis_only' before running this script. The
%   script will find each Cycle*/Dual_analysis3 source result, reuse its
%   saved traces/ROI/info, and write the new analysis outputs under
%   Cycle*/Dual_analysis3_analysis_only/<run_name>. The source folder is
%   read as input and is not used as save_path.
%
% This script keeps Dual_analysis3 as a normal runnable script.
% It passes parameters in from an outer scope by calling run(...).
clear all
%% Input Setup
% Point this to one rebuilt record folder that contains Cycle* subfolders.

% -------------------------------------------------------------------------
% A. Main batch scope
% -------------------------------------------------------------------------
if ~exist('rec_path', 'var') || isempty(rec_path)
    rec_path = 'V:\Luorong\Invivo\26.06.02_dual_color_WT\Methods3_drifting_grating\Rec1_2026-06-02_21-29-09';
end
if ~exist('analysis_backend', 'var') || isempty(analysis_backend)
    analysis_backend = 'default';   % 'default' | 'volpy_voltage_reanalysis'
end
if ~exist('analysis_mode', 'var') || isempty(analysis_mode)
    analysis_mode = 'full';         % 'full' | 'analysis_only'
end
analysis_only_rec_mode = strcmpi(string(analysis_mode), "analysis_only");
source_result_root_dir = resolve_result_root_dir_rec(analysis_backend);
result_root_dir = source_result_root_dir;
average_output_dir_name = resolve_average_output_dir_name_rec(analysis_backend);
batch_summary_file_name = resolve_batch_summary_filename_rec(analysis_backend);
if analysis_only_rec_mode
    result_root_dir = resolve_analysis_only_output_root_dir_rec();
    average_output_dir_name = 'Dual_analysis3_average';
    batch_summary_file_name = 'Dual_analysis3_batch.mat';
end
if ~exist('analysis_only_output_root_dir', 'var') || isempty(analysis_only_output_root_dir)
    analysis_only_output_root_dir = resolve_analysis_only_output_root_dir_rec();
end
if ~exist('analysis_only_run_tag', 'var') || isempty(analysis_only_run_tag)
    analysis_only_run_tag = string(datetime('now', 'Format', 'yyyy-MM-dd HH-mm-ss'));
end

% Choose one cycle as the ROI reference.
if ~exist('reference_cycle_name', 'var') || isempty(reference_cycle_name)
    reference_cycle_name = 'Cycle1';
end

% If you already have a good ROI file, fill it here and the script will
% skip the reference ROI selection step.
if ~exist('reference_roi_file', 'var') || isempty(reference_roi_file)
    reference_roi_file = '';
end

% Reuse the latest ROI result already present in the reference cycle when
% possible. If false, the reference cycle will be rerun to generate ROI.
if ~exist('reuse_existing_reference_roi', 'var') || isempty(reuse_existing_reference_roi)
    reuse_existing_reference_roi = false;
end

% When true, rerun the reference cycle first to create a fresh ROI file.
if ~exist('rerun_reference_cycle_for_roi', 'var') || isempty(rerun_reference_cycle_for_roi)
    rerun_reference_cycle_for_roi = false;
end

% Whether the reference cycle should enter manual inter-camera ROI offset
% correction before ROI drawing.
% Modes:
%   'none'            -> keep reuse_offset / [0 0]
%   'manual_points'   -> click one matching point on the two average images
%   'matlab_register' -> use MATLAB built-in translation registration
% Legacy true/false is still accepted and maps to manual_points/none.
if ~exist('reference_correct_offset_mode', 'var') || isempty(reference_correct_offset_mode)
    if exist('reference_correct_offset', 'var') && ~isempty(reference_correct_offset)
        reference_correct_offset_mode = reference_correct_offset;
    else
        reference_correct_offset_mode = 'none' ;
    end
end
reference_correct_offset_mode = normalize_correct_offset_mode_rec(reference_correct_offset_mode);

% After creating/reusing the reference ROI, also run the full analysis on
% the reference cycle inside the batch loop.
if ~exist('run_reference_cycle_in_batch', 'var') || isempty(run_reference_cycle_in_batch)
    run_reference_cycle_in_batch = true;
end

% Leave empty to run every Cycle*. Otherwise provide a string array like:
% ["Cycle1","Cycle3","Cycle5"] %strings(0, 1)
if ~exist('cycle_name_filter', 'var') || isempty(cycle_name_filter)
    cycle_name_filter = strings(0, 1); %strings(0, 1)
end

% Skip a cycle if it already contains any previous Dual_analysis3 result.
% Leave this on for normal Rec-level reruns so the script can jump straight
% to the final record-average summary without recomputing every cycle.
if ~exist('skip_cycles_with_existing_results', 'var') || isempty(skip_cycles_with_existing_results)
    skip_cycles_with_existing_results = false;
end

% When true, do not revisit the per-cycle workflow. Instead, reopen the
% existing batch summary file and build only the final record-average
% outputs. This is the fast path once all cycles have already been run.
if ~exist('run_record_average_only', 'var') || isempty(run_record_average_only)
    run_record_average_only = false;
end
% Build the cross-cycle record average only when explicitly requested. This
% keeps routine per-cycle batch runs from failing just because no completed
% cycle is available yet, or because the user only wanted per-cycle outputs.
if ~exist('run_record_average', 'var') || isempty(run_record_average)
    run_record_average = false;
end
run_record_average = logical(run_record_average) || logical(run_record_average_only);
if ~exist('stop_on_cycle_error', 'var') || isempty(stop_on_cycle_error)
    stop_on_cycle_error = true;
end

% -------------------------------------------------------------------------
% B. Common per-channel manual settings
% These map onto camera_cfg automatically so routine reruns usually only
% need edits here instead of manual struct rewriting below.
% -------------------------------------------------------------------------
if ~exist('voltage_frame_rate', 'var') || isempty(voltage_frame_rate)
    voltage_frame_rate = 400;
end
if ~exist('calcium_frame_rate', 'var') || isempty(calcium_frame_rate)
    calcium_frame_rate = 400;
end
if ~exist('voltage_transpose_movie', 'var') || isempty(voltage_transpose_movie)
    voltage_transpose_movie = false;
end
if ~exist('calcium_transpose_movie', 'var') || isempty(calcium_transpose_movie)
    calcium_transpose_movie = true;
end
if ~exist('run_motion_correction', 'var') || isempty(run_motion_correction)
    run_motion_correction = true;
end
if ~exist('correct_offset_mode', 'var') || isempty(correct_offset_mode)
    correct_offset_mode = 'none';
end
correct_offset_mode = normalize_correct_offset_mode_rec(correct_offset_mode);

% -------------------------------------------------------------------------
% C. Advanced shared Dual_analysis3 overrides
% -------------------------------------------------------------------------
% Shared Dual_analysis3 overrides.
if ~exist('camera_cfg', 'var') || isempty(camera_cfg)
    camera_cfg(1) = struct( ...
        'camera_index', 1, ...
        'role', "calcium", ...
        'transpose_before_analysis', true, ...
        'frame_rate', 400);
    camera_cfg(2) = struct( ...
        'camera_index', 2, ...
        'role', "voltage", ...
        'transpose_before_analysis', false, ...
        'frame_rate', 400);
end
for camera_idx = 1:numel(camera_cfg)
    current_role = string(camera_cfg(camera_idx).role);
    if strcmpi(current_role, "voltage")
        camera_cfg(camera_idx).transpose_before_analysis = logical(voltage_transpose_movie);
        camera_cfg(camera_idx).frame_rate = double(voltage_frame_rate);
    elseif strcmpi(current_role, "calcium")
        camera_cfg(camera_idx).transpose_before_analysis = logical(calcium_transpose_movie);
        camera_cfg(camera_idx).frame_rate = double(calcium_frame_rate);
    end
end

if ~exist('map_bin', 'var') || isempty(map_bin)
    map_bin = 4;
end
if ~exist('calcium_smoothing_window', 'var') || isempty(calcium_smoothing_window)
    calcium_smoothing_window = 40;
end
if ~exist('bleach_mode_voltage', 'var') || isempty(bleach_mode_voltage)
    bleach_mode_voltage = 'linear';
end
if ~exist('bleach_mode_calcium', 'var') || isempty(bleach_mode_calcium)
    bleach_mode_calcium = 'exp2';
end
if ~exist('run_background_removal', 'var') || isempty(run_background_removal)
    run_background_removal = true;
end
if ~exist('voltage_polarity', 'var') || isempty(voltage_polarity)
    voltage_polarity = -1;
end
if ~exist('calcium_polarity', 'var') || isempty(calcium_polarity)
    calcium_polarity = 1;
end
if ~exist('reuse_offset', 'var') || isempty(reuse_offset)
    reuse_offset = [0,0];
end
if ~exist('gpu', 'var') || isempty(gpu)
    gpu = true;
end
if ~exist('enable_batch_diary', 'var') || isempty(enable_batch_diary)
    enable_batch_diary = true;
end
if ~exist('motion_cfg', 'var') || isempty(motion_cfg)
    motion_cfg = struct( ...
        'enabled', true, ...
        'use_saved_shift', false, ...
        'saved_shift_file', '', ...
        'highpass', true, ...
        'auto_reuse_previous_shift', true);
end
motion_cfg.enabled = logical(run_motion_correction);

%% Resolve Cycles
script_dir = fileparts(mfilename('fullpath'));
dual_script_path = fullfile(script_dir, 'Dual_analysis3.m');
if ~isfile(dual_script_path)
    error('Cannot find Dual_analysis3.m next to this batch script.');
end

if ~isfolder(rec_path)
    error('Record path does not exist: %s', rec_path);
end
batch_summary_file = fullfile(rec_path, batch_summary_file_name);
if enable_batch_diary
    diary_file = fullfile(rec_path, sprintf('Dual_analysis3_rec_%s.log', char(datetime('now', 'Format', 'yyyy-MM-dd HH-mm-ss'))));
    diary(diary_file);
    diary_cleanup = onCleanup(@() diary('off'));
    fprintf('[Batch] Command-window log is being written to:\n  %s\n', diary_file);
end

if run_record_average_only
    cycle_dirs = dir(fullfile(rec_path, 'Cycle*'));
    cycle_dirs = cycle_dirs([cycle_dirs.isdir]);
    cycle_dirs = sort_cycle_dirs(cycle_dirs);

    if isempty(cycle_dirs)
        error('No Cycle* folders found under record path: %s', rec_path);
    end

    all_cycle_names = string({cycle_dirs.name});
    if ~isempty(cycle_name_filter)
        keep_mask = ismember(all_cycle_names, string(cycle_name_filter(:)));
        cycle_dirs = cycle_dirs(keep_mask);
    end

    if isempty(cycle_dirs)
        error('No cycles remain after applying cycle_name_filter.');
    end

    batch_results = repmat(struct( ...
        'cycle_name', "", ...
        'cycle_path', "", ...
        'status', "", ...
        'save_path', "", ...
        'roi_file', "", ...
        'message', ""), 0, 1);

    for i = 1:numel(cycle_dirs)
        current_cycle_name = string(cycle_dirs(i).name);
        current_cycle_path = fullfile(cycle_dirs(i).folder, cycle_dirs(i).name);
        existing_result = find_latest_explicit_result(current_cycle_path, analysis_backend);
        if strlength(existing_result) > 0
            batch_results(end+1, 1) = struct( ...
                'cycle_name', current_cycle_name, ...
                'cycle_path', string(current_cycle_path), ...
                'status', "existing_result_detected", ...
                'save_path', fileparts(existing_result), ...
                'roi_file', "", ...
                'message', existing_result);
        else
            batch_results(end+1, 1) = struct( ...
                'cycle_name', current_cycle_name, ...
                'cycle_path', string(current_cycle_path), ...
                'status', "missing_result", ...
                'save_path', "", ...
                'roi_file', "", ...
                'message', "No -1_explicit_dual_results.mat found under this cycle.");
        end
    end

    if strlength(string(reference_roi_file)) == 0
        reference_cycle_path = fullfile(rec_path, reference_cycle_name);
        reference_roi_file = string(find_latest_dual_roi_file(reference_cycle_path));
    else
        reference_roi_file = string(reference_roi_file);
    end

    fprintf('\n============================================================\n');
    fprintf('[Batch] Summary-Only Mode\n');
    fprintf('============================================================\n');
    fprintf('Scanning existing %s outputs directly from cycle folders.\n', result_root_dir);
    for i = 1:numel(batch_results)
        if strlength(batch_results(i).save_path) > 0
            fprintf('%s | using existing result folder:\n  %s\n', ...
                batch_results(i).cycle_name, batch_results(i).save_path);
        else
            fprintf('%s | missing existing result folder.\n', batch_results(i).cycle_name);
        end
    end
else
    cycle_dirs = dir(fullfile(rec_path, 'Cycle*'));
    cycle_dirs = cycle_dirs([cycle_dirs.isdir]);
    cycle_dirs = sort_cycle_dirs(cycle_dirs);

    if isempty(cycle_dirs)
        error('No Cycle* folders found under record path: %s', rec_path);
    end

    all_cycle_names = string({cycle_dirs.name});
    if ~isempty(cycle_name_filter)
        keep_mask = ismember(all_cycle_names, string(cycle_name_filter(:)));
        cycle_dirs = cycle_dirs(keep_mask);
        all_cycle_names = string({cycle_dirs.name});
    end

    if isempty(cycle_dirs)
        error('No cycles remain after applying cycle_name_filter.');
    end

    fprintf('\n============================================================\n');
    fprintf('[Batch] Existing Result Detection\n');
    fprintf('============================================================\n');
    for i = 1:numel(cycle_dirs)
        cycle_path_i = fullfile(cycle_dirs(i).folder, cycle_dirs(i).name);
        existing_run_i = find_latest_explicit_result(cycle_path_i, analysis_backend);
        if strlength(existing_run_i) > 0
            fprintf('%s | existing result marker detected:\n  %s\n', cycle_dirs(i).name, existing_run_i);
        else
            fprintf('%s | no processed %s run detected.\n', cycle_dirs(i).name, source_result_root_dir);
        end
    end

    if analysis_only_rec_mode
        %% Run Analysis-Only Forks From Existing Cycle Results
        fprintf('\n============================================================\n');
        fprintf('[Batch] Run Analysis-Only Forks\n');
        fprintf('============================================================\n');
        fprintf('Source root: %s\n', resolve_result_root_dir_rec(analysis_backend));
        fprintf('Output root: %s\n', analysis_only_output_root_dir);
        fprintf('Run tag: %s\n', analysis_only_run_tag);

        batch_results = repmat(struct( ...
            'cycle_name', "", ...
            'cycle_path', "", ...
            'status', "", ...
            'save_path', "", ...
            'roi_file', "", ...
            'message', ""), 0, 1);

        for i = 1:numel(cycle_dirs)
            current_cycle_name = string(cycle_dirs(i).name);
            current_cycle_path = fullfile(cycle_dirs(i).folder, cycle_dirs(i).name);
            source_result = find_latest_explicit_result(current_cycle_path, analysis_backend);

            if strlength(source_result) == 0
                batch_results(end+1, 1) = struct( ...
                    'cycle_name', current_cycle_name, ...
                    'cycle_path', string(current_cycle_path), ...
                    'status', "missing_source_result", ...
                    'save_path', "", ...
                    'roi_file', "", ...
                    'message', "No existing Dual_analysis3 source result was found for analysis_only reuse.");
                fprintf('%s | missing source result; skipped.\n', current_cycle_name);
                continue;
            end

            if skip_cycles_with_existing_results
                existing_fork = find_latest_analysis_only_fork_result(current_cycle_path, analysis_only_output_root_dir);
                if strlength(existing_fork) > 0
                    existing_fork_dir = string(fileparts(existing_fork));
                    batch_results(end+1, 1) = struct( ...
                        'cycle_name', current_cycle_name, ...
                        'cycle_path', string(current_cycle_path), ...
                        'status', "skipped_existing_analysis_only_result", ...
                        'save_path', existing_fork_dir, ...
                        'roi_file', string(fullfile(existing_fork_dir, '1_dual_roi_results.mat')), ...
                        'message', existing_fork);
                    fprintf('Skipping %s because an analysis-only fork already exists:\n  %s\n', current_cycle_name, existing_fork);
                    continue;
                end
            end

            fprintf('\n[Batch] Analysis-only fork for %s ...\n', current_cycle_name);
            fprintf('  source: %s\n', fileparts(source_result));
            try
                run_result = run_dual_cycle_analysis_only_fork( ...
                    dual_script_path, current_cycle_path, fileparts(source_result), ...
                    analysis_only_output_root_dir, analysis_only_run_tag, ...
                    camera_cfg, map_bin, calcium_smoothing_window, ...
                    bleach_mode_voltage, bleach_mode_calcium, run_background_removal, ...
                    voltage_polarity, calcium_polarity, ...
                    reuse_offset, gpu, motion_cfg);

                batch_results(end+1, 1) = struct( ...
                    'cycle_name', current_cycle_name, ...
                    'cycle_path', string(current_cycle_path), ...
                    'status', "completed", ...
                    'save_path', string(run_result.save_path), ...
                    'roi_file', string(run_result.roi_file), ...
                    'message', string(run_result.source_save_path));
            catch ME
                error_report = getReport(ME, 'extended', 'hyperlinks', 'off');
                error_log_file = save_cycle_error_report_rec(rec_path, current_cycle_name, error_report);
                batch_results(end+1, 1) = struct( ...
                    'cycle_name', current_cycle_name, ...
                    'cycle_path', string(current_cycle_path), ...
                    'status', "failed", ...
                    'save_path', "", ...
                    'roi_file', "", ...
                    'message', string(error_report));
                fprintf(2, '[Batch] %s analysis-only fork failed:\n%s\n', current_cycle_name, error_report);
                fprintf(2, '[Batch] Error report saved to:\n  %s\n', error_log_file);
                if stop_on_cycle_error
                    rethrow(ME);
                end
            end
        end

        if strlength(string(reference_roi_file)) == 0
            completed_idx = find(strcmpi(string({batch_results.status}), "completed") ...
                & strlength(string({batch_results.roi_file})) > 0, 1, 'first');
            if ~isempty(completed_idx)
                reference_roi_file = string(batch_results(completed_idx).roi_file);
            else
                reference_roi_file = "";
            end
        end
    else
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
        reference_existing_run = find_latest_explicit_result(reference_cycle_path, analysis_backend);
        if strlength(reference_existing_run) > 0
            fprintf('Reference cycle already has a result marker:\n  %s\n', reference_existing_run);
        end

        if strlength(string(reference_roi_file)) == 0
            reference_roi_file = "";
        else
            reference_roi_file = string(reference_roi_file);
        end

        if strlength(reference_roi_file) == 0 && reuse_existing_reference_roi && ~rerun_reference_cycle_for_roi
            reference_roi_file = string(find_latest_dual_roi_file(reference_cycle_path));
            if strlength(reference_roi_file) > 0
                fprintf('Reusing latest ROI file already present in reference cycle:\n  %s\n', reference_roi_file);
            end
        end

        if strlength(reference_roi_file) == 0 || rerun_reference_cycle_for_roi
            fprintf('Running reference cycle to create/recreate ROI...\n');
            reference_result = run_dual_cycle_with_overrides( ...
                dual_script_path, reference_cycle_path, "", reference_correct_offset_mode, ...
                camera_cfg, map_bin, calcium_smoothing_window, ...
                bleach_mode_voltage, bleach_mode_calcium, run_background_removal, ...
                voltage_polarity, calcium_polarity, ...
                reuse_offset, gpu, motion_cfg, analysis_backend);
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
                    'message', "Reference cycle already used to create/reuse ROI.");
                fprintf('Skipping %s in batch loop because it already served as the ROI reference.\n', current_cycle_name);
                continue;
            end

            if skip_cycles_with_existing_results
                existing_result = find_latest_explicit_result(current_cycle_path, analysis_backend);
                if strlength(existing_result) > 0
                    batch_results(end+1, 1) = struct( ...
                        'cycle_name', current_cycle_name, ...
                        'cycle_path', string(current_cycle_path), ...
                        'status', "skipped_existing_result", ...
                        'save_path', fileparts(existing_result), ...
                        'roi_file', reference_roi_file, ...
                        'message', existing_result);
                    fprintf('Skipping %s because the final output already exists:\n  %s\n', current_cycle_name, existing_result);
                    continue;
                end
            end

            fprintf('\n[Batch] Running %s ...\n', current_cycle_name);
            try
                run_result = run_dual_cycle_with_overrides( ...
                    dual_script_path, current_cycle_path, reference_roi_file, correct_offset_mode, ...
                    camera_cfg, map_bin, calcium_smoothing_window, ...
                    bleach_mode_voltage, bleach_mode_calcium, run_background_removal, ...
                    voltage_polarity, calcium_polarity, ...
                    reuse_offset, gpu, motion_cfg, analysis_backend);

                batch_results(end+1, 1) = struct( ...
                    'cycle_name', current_cycle_name, ...
                    'cycle_path', string(current_cycle_path), ...
                    'status', "completed", ...
                    'save_path', string(run_result.save_path), ...
                    'roi_file', string(run_result.roi_file), ...
                    'message', "");
            catch ME
                error_report = getReport(ME, 'extended', 'hyperlinks', 'off');
                error_log_file = save_cycle_error_report_rec(rec_path, current_cycle_name, error_report);
                batch_results(end+1, 1) = struct( ...
                    'cycle_name', current_cycle_name, ...
                    'cycle_path', string(current_cycle_path), ...
                    'status', "failed", ...
                    'save_path', "", ...
                    'roi_file', reference_roi_file, ...
                    'message', string(error_report));
                fprintf(2, '[Batch] %s failed:\n%s\n', current_cycle_name, error_report);
                fprintf(2, '[Batch] Error report saved to:\n  %s\n', error_log_file);
                if stop_on_cycle_error
                    rethrow(ME);
                end
            end
        end
    end

    save(batch_summary_file, 'batch_results', 'reference_roi_file', 'analysis_mode', 'analysis_only_output_root_dir', 'analysis_only_run_tag');
    fprintf('\nBatch summary saved to:\n  %s\n', batch_summary_file);
end

%% Record-Level Average Trace Summary
% After the per-cycle analyses are ready, reopen their saved trace results
% and build one record-level average trace per ROI. This keeps the
% expensive movie-based sections untouched while still letting the user see
% how stable each ROI is across cycles.
if run_record_average
    fprintf('\n============================================================\n');
    fprintf('[Batch] Record-Level Average Trace Summary\n');
    fprintf('============================================================\n');
    record_average_output = build_record_average_summary( ...
        rec_path, batch_results, reference_roi_file, ...
        voltage_polarity, calcium_polarity, calcium_smoothing_window, average_output_dir_name);
    fprintf('Record-average summary folder:\n  %s\n', record_average_output.output_dir);
    fprintf('Record-average result bundle:\n  %s\n', record_average_output.results_file);
else
    fprintf('\n[Batch] Record-level average skipped. Set run_record_average = true when cross-cycle averaging is needed.\n');
end

%% Local Functions
function result = run_dual_cycle_with_overrides( ...
    dual_script_path, cycle_path, reuse_roi_file, correct_offset_mode, ...
    camera_cfg, map_bin, calcium_smoothing_window, ...
    bleach_mode_voltage, bleach_mode_calcium, run_background_removal, ...
    voltage_polarity, calcium_polarity, ...
    reuse_offset, gpu, motion_cfg, analysis_backend)

% All variables defined in this function are visible to Dual_analysis3.m
% when it is executed via run(...). This keeps Dual_analysis3 as a script
% while still letting an outer batch driver provide overrides.
[record_path_for_name, cycle_name_for_name] = fileparts(cycle_path);
[~, record_name_for_name] = fileparts(record_path_for_name);
time_tag = string(datetime('now', 'Format', 'yyyy-MM-dd HH-mm-ss-SSS'));
analysis_run_name = sprintf('%s_%s_%s', record_name_for_name, cycle_name_for_name, char(time_tag));
if strcmpi(string(analysis_backend), "volpy_voltage_reanalysis")
    standard_result = find_latest_explicit_result(cycle_path, 'default');
    if strlength(standard_result) == 0
        error('VolPy voltage re-analysis requires an existing standard Dual_analysis3 result in %s.', cycle_path);
    end
    volpy_source_results_path = string(fileparts(standard_result));
end

run(dual_script_path);

result = struct( ...
    'save_path', save_path, ...
    'roi_file', fullfile(save_path, '1_dual_roi_results.mat'));
end

function result = run_dual_cycle_analysis_only_fork( ...
    dual_script_path, cycle_path, source_results_path, output_root_dir, run_tag, ...
    camera_cfg, map_bin, calcium_smoothing_window, ...
    bleach_mode_voltage, bleach_mode_calcium, run_background_removal, ...
    voltage_polarity, calcium_polarity, ...
    reuse_offset, gpu, motion_cfg)

% Reuse saved ROI/trace/stim context from an existing Dual_analysis3 folder,
% but write the rerun products to a new folder under output_root_dir.
[record_path_for_name, cycle_name_for_name] = fileparts(cycle_path);
[~, record_name_for_name] = fileparts(record_path_for_name);
analysis_mode = 'analysis_only';
analysis_backend = 'default';
analysis_run_name = sprintf('%s_%s_analysis_only_%s', ...
    record_name_for_name, cycle_name_for_name, char(string(run_tag)));
analysis_run_name = sanitize_path_component_rec(analysis_run_name);
reuse_results_path = char(string(source_results_path));
save_path = fullfile(cycle_path, char(string(output_root_dir)), analysis_run_name);
reuse_roi_file = '';
correct_offset_mode = 'none';

run(dual_script_path);

result = struct( ...
    'save_path', save_path, ...
    'roi_file', fullfile(save_path, '1_dual_roi_results.mat'), ...
    'source_save_path', string(source_results_path));
end

function correct_offset_mode = normalize_correct_offset_mode_rec(correct_offset_mode)
if islogical(correct_offset_mode) || (isnumeric(correct_offset_mode) && isscalar(correct_offset_mode))
    if logical(correct_offset_mode)
        correct_offset_mode = "manual_points";
    else
        correct_offset_mode = "none";
    end
    return;
end

correct_offset_mode = lower(strtrim(string(correct_offset_mode)));
if numel(correct_offset_mode) ~= 1
    error('correct_offset_mode must resolve to one scalar mode.');
end

switch correct_offset_mode
    case {"none", "off", "false", "0", "reuse_or_zero"}
        correct_offset_mode = "none";
    case {"manual", "manual_point", "manual_points", "points", "point"}
        correct_offset_mode = "manual_points";
    case {"matlab_register", "register", "imreg", "imregtform"}
        correct_offset_mode = "matlab_register";
    otherwise
        error(['Unsupported correct_offset_mode: %s. Use ''none'', ' ...
            '''manual_points'', or ''matlab_register''.'], char(correct_offset_mode));
end
end

function cycle_dirs = sort_cycle_dirs(cycle_dirs)
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

function result_root_dir = resolve_result_root_dir_rec(analysis_backend)
if strcmpi(string(analysis_backend), "volpy_voltage_reanalysis")
    result_root_dir = 'Dual_analysis3_volpy_voltage';
else
    result_root_dir = 'Dual_analysis3';
end
end

function result_root_dir = resolve_analysis_only_output_root_dir_rec()
result_root_dir = 'Dual_analysis3_analysis_only';
end

function average_output_dir_name = resolve_average_output_dir_name_rec(analysis_backend)
if strcmpi(string(analysis_backend), "volpy_voltage_reanalysis")
    average_output_dir_name = 'Dual_analysis3_rec_average_volpy_voltage';
else
    average_output_dir_name = 'Dual_analysis3_rec_average';
end
end

function batch_summary_file_name = resolve_batch_summary_filename_rec(analysis_backend)
if strcmpi(string(analysis_backend), "volpy_voltage_reanalysis")
    batch_summary_file_name = 'Dual_analysis3_rec_batch_summary_volpy_voltage.mat';
else
    batch_summary_file_name = 'Dual_analysis3_rec_batch_summary.mat';
end
end

function roi_file = find_latest_dual_roi_file(cycle_path)
roi_file = "";
listing = dir(fullfile(cycle_path, 'Dual_analysis3', '**', '1_dual_roi_results.mat'));
if isempty(listing)
    return;
end
[~, newest_idx] = max([listing.datenum]);
roi_file = string(fullfile(listing(newest_idx).folder, listing(newest_idx).name));
end

function explicit_result = find_latest_explicit_result(cycle_path, analysis_backend)
explicit_result = "";
if nargin < 2 || isempty(analysis_backend)
    analysis_backend = 'default';
end
result_root_dir = resolve_result_root_dir_rec(analysis_backend);

% Prefer the final explicit bundle when a cycle completed all Dual_analysis3
% sections. If a run stopped after ROI creation, fall back to the ROI marker
% so rec-level workflows can still discover that partial result folder.
explicit_listing = dir(fullfile(cycle_path, result_root_dir, '**', '-1_explicit_dual_results.mat'));
explicit_result = select_latest_file_from_listing_rec(explicit_listing);
if strlength(explicit_result) > 0
    return;
end

roi_listing = dir(fullfile(cycle_path, result_root_dir, '**', '1_dual_roi_results.mat'));
explicit_result = select_latest_file_from_listing_rec(roi_listing);
end

function file_path = select_latest_file_from_listing_rec(listing)
file_path = "";
if isempty(listing)
    return;
end
[~, newest_idx] = max([listing.datenum]);
file_path = string(fullfile(listing(newest_idx).folder, listing(newest_idx).name));
end

function explicit_result = find_latest_analysis_only_fork_result(cycle_path, output_root_dir)
explicit_result = "";
explicit_listing = dir(fullfile(cycle_path, char(string(output_root_dir)), '**', '-1_explicit_dual_results.mat'));
explicit_result = select_latest_file_from_listing_rec(explicit_listing);
if strlength(explicit_result) > 0
    return;
end
roi_listing = dir(fullfile(cycle_path, char(string(output_root_dir)), '**', '1_dual_roi_results.mat'));
explicit_result = select_latest_file_from_listing_rec(roi_listing);
end

function error_log_file = save_cycle_error_report_rec(rec_path, cycle_name, error_report)
error_dir = fullfile(rec_path, 'Dual_analysis3_rec_error_logs');
if ~isfolder(error_dir)
    mkdir(error_dir);
end
time_tag = char(datetime('now', 'Format', 'yyyy-MM-dd HH-mm-ss'));
safe_cycle_name = sanitize_path_component_rec(cycle_name);
error_log_file = fullfile(error_dir, sprintf('%s_%s_error.txt', safe_cycle_name, time_tag));
fid = fopen(error_log_file, 'w');
if fid < 0
    warning('Dual_analysis3_rec:ErrorLogOpenFailed', 'Could not write cycle error report to %s', error_log_file);
    return;
end
cleanup_obj = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', char(string(error_report)));
end

function path_component = sanitize_path_component_rec(path_component)
path_component = char(string(path_component));
path_component = regexprep(path_component, '[:*?"<>|]', '-');
path_component = strrep(path_component, filesep, '-');
path_component = strrep(path_component, '/', '-');
end

function output = build_record_average_summary(rec_path, batch_results, reference_roi_file, voltage_polarity, calcium_polarity, calcium_smoothing_window, average_output_dir_name)
output_dir = fullfile(rec_path, average_output_dir_name);
if ~isfolder(output_dir)
    mkdir(output_dir);
end

[result_dirs, cycle_names] = resolve_successful_result_dirs(batch_results);
if isempty(result_dirs)
    error('No successful cycle result folders are available for record-level averaging.');
end

fprintf('Using %d cycle result folders for record-average analysis.\n', numel(result_dirs));
cycle_entry_list = cell(numel(result_dirs), 1);
for idx = 1:numel(result_dirs)
    fprintf('  Loading %s\n', result_dirs(idx));
    cycle_entry_list{idx} = load_cycle_average_entry(result_dirs(idx), cycle_names(idx));
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
    'alignment_mode', "flash_onset", ...
    'calcium_smoothing_window', calcium_smoothing_window, ...
    'voltage_polarity', voltage_polarity, ...
    'calcium_polarity', calcium_polarity);

record_average.voltage.raw = average_cycle_stage(cycle_entries, 'voltage', 'raw', 'flash_onset');
record_average.voltage.sensitivity = average_cycle_stage(cycle_entries, 'voltage', 'sensitivity', 'flash_onset');
record_average.voltage.snr = average_cycle_stage(cycle_entries, 'voltage', 'snr', 'flash_onset');
record_average.calcium.raw = average_cycle_stage(cycle_entries, 'calcium', 'raw', 'flash_onset');
record_average.calcium.sensitivity = average_cycle_stage(cycle_entries, 'calcium', 'sensitivity', 'flash_onset');
record_average.calcium.snr = average_cycle_stage(cycle_entries, 'calcium', 'snr', 'flash_onset');

record_average.stim_windows = resolve_record_average_stim_windows( ...
    cycle_entries, ...
    record_average.voltage.raw.time, ...
    record_average.calcium.raw.time);

raw_summary_png = fullfile(output_dir, '1_record_average_raw_summary.png');
sens_summary_png = fullfile(output_dir, '2_record_average_sensitivity_summary.png');
snr_summary_png = fullfile(output_dir, '3_record_average_snr_summary.png');
sens_overlap_png = fullfile(output_dir, '4_record_average_sensitivity_overlap.png');
snr_overlap_png = fullfile(output_dir, '5_record_average_snr_overlap.png');

plot_record_stage_summary( ...
    record_average.voltage.raw, record_average.calcium.raw, ...
    record_average.stim_windows, 'Raw', raw_summary_png, ...
    [0.85 0.15 0.15], [0.10 0.65 0.20], 1, 1, ...
    'Voltage Raw', sprintf('Calcium Raw (%s, window=%d)', record_average.calcium.raw.stage_name, calcium_smoothing_window));
plot_record_stage_summary( ...
    record_average.voltage.sensitivity, record_average.calcium.sensitivity, ...
    record_average.stim_windows, 'Sensitivity', sens_summary_png, ...
    [0.85 0.15 0.15], [0.10 0.65 0.20], voltage_polarity, calcium_polarity, ...
    'Voltage Sensitivity', sprintf('Calcium Sensitivity (%s, window=%d)', record_average.calcium.sensitivity.stage_name, calcium_smoothing_window));
plot_record_stage_summary( ...
    record_average.voltage.snr, record_average.calcium.snr, ...
    record_average.stim_windows, 'SNR', snr_summary_png, ...
    [0.85 0.15 0.15], [0.10 0.65 0.20], voltage_polarity, calcium_polarity, ...
    'Voltage SNR', sprintf('Calcium SNR (%s, window=%d)', record_average.calcium.snr.stage_name, calcium_smoothing_window));

plot_record_metric_overlap( ...
    record_average.voltage.sensitivity, record_average.calcium.sensitivity, ...
    record_average.stim_windows, 'Sensitivity', sens_overlap_png, ...
    [0.85 0.15 0.15], [0.10 0.65 0.20], voltage_polarity, calcium_polarity);
plot_record_metric_overlap( ...
    record_average.voltage.snr, record_average.calcium.snr, ...
    record_average.stim_windows, 'SNR', snr_overlap_png, ...
    [0.85 0.15 0.15], [0.10 0.65 0.20], voltage_polarity, calcium_polarity);

tf_input = struct( ...
    'voltage', struct( ...
        'trace', voltage_polarity * record_average.voltage.sensitivity.average, ...
        'time', record_average.voltage.sensitivity.time, ...
        'frame_rate', record_average.voltage.sensitivity.frame_rate, ...
        'stage_name', "sensitivity_average"), ...
    'calcium', struct( ...
        'trace', calcium_polarity * record_average.calcium.sensitivity.average, ...
        'time', record_average.calcium.sensitivity.time, ...
        'frame_rate', record_average.calcium.sensitivity.frame_rate, ...
        'stage_name', "sensitivity_average"));
record_average.time_frequency = build_record_average_time_frequency( ...
    tf_input, record_average.stim_windows, output_dir);

record_average.visualizations = struct( ...
    'raw_summary_png', string(raw_summary_png), ...
    'sensitivity_summary_png', string(sens_summary_png), ...
    'snr_summary_png', string(snr_summary_png), ...
    'sensitivity_overlap_png', string(sens_overlap_png), ...
    'snr_overlap_png', string(snr_overlap_png), ...
    'fourier_summary_png', string(record_average.time_frequency.visualizations.fourier_summary_png), ...
    'wavelet_summary_png', string(record_average.time_frequency.visualizations.wavelet_summary_png));

results_file = fullfile(output_dir, 'record_average_results.mat');
save(results_file, 'record_average', '-v7.3');

output = struct( ...
    'output_dir', string(output_dir), ...
    'results_file', string(results_file));
end

function [result_dirs, cycle_names] = resolve_successful_result_dirs(batch_results)
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

function entry = load_cycle_average_entry(result_dir, cycle_name)
voltage_path = fullfile(result_dir, 'voltage_results.mat');
calcium_path = fullfile(result_dir, 'calcium_results.mat');
stim_path = fullfile(result_dir, 'stim_results.mat');

tmp_voltage = load(voltage_path, 'voltage_results');
tmp_calcium = load(calcium_path, 'calcium_results');
voltage_results = tmp_voltage.voltage_results;
calcium_results = tmp_calcium.calcium_results;

entry = struct();
entry.cycle_name = string(cycle_name);
entry.result_dir = string(result_dir);
entry.voltage.raw = extract_cycle_stage(voltage_results, {'raw'});
entry.voltage.sensitivity = extract_cycle_stage(voltage_results, {'sensitivity'});
entry.voltage.snr = extract_cycle_stage(voltage_results, {'snr'});
entry.calcium.raw = extract_cycle_stage(calcium_results, {'raw_smoothed', 'raw'});
entry.calcium.sensitivity = extract_cycle_stage(calcium_results, {'sensitivity_smoothed', 'sensitivity'});
entry.calcium.snr = extract_cycle_stage(calcium_results, {'snr_smoothed', 'snr'});

entry.stim_windows = struct();
if isfile(stim_path)
    tmp_stim = load(stim_path, 'stim_results');
    if isfield(tmp_stim, 'stim_results') && isstruct(tmp_stim.stim_results) ...
            && isfield(tmp_stim.stim_results, 'windows')
        entry.stim_windows = tmp_stim.stim_results.windows;
    end
end
entry.alignment = struct( ...
    'voltage_flash_onset_time', resolve_flash_onset_time(entry.stim_windows, 'voltage'), ...
    'calcium_flash_onset_time', resolve_flash_onset_time(entry.stim_windows, 'calcium'));
end

function stage = extract_cycle_stage(results_struct, preferred_stages)
[stage_data, stage_name] = resolve_preferred_stage_from_results(results_struct, preferred_stages);
frame_rate = double(results_struct.movie_info.frame_rate);
nframes = size(stage_data, 1);
stage = struct( ...
    'stage_name', string(stage_name), ...
    'data', double(stage_data), ...
    'frame_rate', frame_rate, ...
    'time', (1:nframes)' / frame_rate);
end

function [stage_data, stage_name] = resolve_preferred_stage_from_results(results_struct, preferred_stages)
for idx = 1:numel(preferred_stages)
    stage_name = preferred_stages{idx};
    if isfield(results_struct, 'trace_results') ...
            && isfield(results_struct.trace_results, stage_name) ...
            && isfield(results_struct.trace_results.(stage_name), 'data')
        stage_data = results_struct.trace_results.(stage_name).data;
        return;
    end
end
error('Required stage is missing. Tried: %s', strjoin(preferred_stages, ', '));
end

function averaged = average_cycle_stage(cycle_entries, channel_name, stage_name, alignment_mode)
sample = cycle_entries(1).(channel_name).(stage_name).data;
nrois = size(sample, 2);
frame_rate = cycle_entries(1).(channel_name).(stage_name).frame_rate;
stage_labels = strings(numel(cycle_entries), 1);
relative_starts = NaN(numel(cycle_entries), 1);
relative_ends = NaN(numel(cycle_entries), 1);
for idx = 1:numel(cycle_entries)
    current_stage = cycle_entries(idx).(channel_name).(stage_name);
    if size(current_stage.data, 2) ~= nrois
        error('ROI count mismatch while averaging %s %s across cycles.', channel_name, stage_name);
    end
    if abs(current_stage.frame_rate - frame_rate) > 1e-9
        error('Frame rate mismatch while averaging %s %s across cycles.', channel_name, stage_name);
    end
    stage_labels(idx) = current_stage.stage_name;
    alignment_time = resolve_cycle_alignment_time(cycle_entries(idx), channel_name, alignment_mode);
    relative_time = current_stage.time(:) - alignment_time;
    relative_starts(idx) = relative_time(1);
    relative_ends(idx) = relative_time(end);
end

common_start = max(relative_starts);
common_end = min(relative_ends);
if ~isfinite(common_start) || ~isfinite(common_end) || common_end <= common_start
    error('No overlapping time range remains after %s alignment for %s %s.', alignment_mode, channel_name, stage_name);
end

dt = 1 / frame_rate;
common_time = (common_start:dt:common_end)';
if numel(common_time) < 2
    common_time = [common_start; common_end];
end

per_cycle = NaN(numel(common_time), nrois, numel(cycle_entries));
alignment_time_by_cycle = NaN(numel(cycle_entries), 1);
for idx = 1:numel(cycle_entries)
    current_data = double(cycle_entries(idx).(channel_name).(stage_name).data);
    current_time = cycle_entries(idx).(channel_name).(stage_name).time(:);
    alignment_time = resolve_cycle_alignment_time(cycle_entries(idx), channel_name, alignment_mode);
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

function stim_windows = resolve_record_average_stim_windows(cycle_entries, voltage_time, calcium_time)
stim_windows = struct();
for idx = 1:numel(cycle_entries)
    candidate = cycle_entries(idx).stim_windows;
    if isstruct(candidate) && isfield(candidate, 'supported') && candidate.supported
        stim_windows = candidate;
        voltage_anchor = cycle_entries(idx).alignment.voltage_flash_onset_time;
        calcium_anchor = cycle_entries(idx).alignment.calcium_flash_onset_time;
        break;
    end
end

if isempty(fieldnames(stim_windows))
    return;
end

if isfield(stim_windows, 'voltage')
    stim_windows.voltage = align_channel_windows_to_flash(stim_windows.voltage, voltage_anchor, voltage_time);
    if isfield(stim_windows, 'flash_windows')
        stim_windows.voltage.flash_windows = stim_windows.flash_windows;
    end
end
if isfield(stim_windows, 'calcium')
    stim_windows.calcium = align_channel_windows_to_flash(stim_windows.calcium, calcium_anchor, calcium_time);
    if isfield(stim_windows, 'flash_windows')
        stim_windows.calcium.flash_windows = stim_windows.flash_windows;
    end
end
end

function channel_windows = align_channel_windows_to_flash(channel_windows, anchor_time, common_time)
if ~isfinite(anchor_time)
    error('Flash onset time is missing, so record-average flash alignment cannot be built.');
end
time_min = min(common_time);
time_max = max(common_time);
field_pairs = { ...
    'baseline_frames', 'baseline_time_ranges'; ...
    'stim_frames', 'stim_time_ranges'; ...
    'block_frames', 'block_time_ranges'; ...
    'shading_frames', 'shading_time_ranges'};

for idx = 1:size(field_pairs, 1)
    frame_field = field_pairs{idx, 1};
    time_field = field_pairs{idx, 2};
    if isfield(channel_windows, time_field) && ~isempty(channel_windows.(time_field))
        time_ranges = double(channel_windows.(time_field));
        time_ranges = time_ranges - anchor_time;
        time_ranges(:, 1) = max(time_min, time_ranges(:, 1));
        time_ranges(:, 2) = min(time_max, time_ranges(:, 2));
        valid = time_ranges(:, 2) >= time_ranges(:, 1);
        channel_windows.(time_field) = time_ranges(valid, :);
        if isfield(channel_windows, frame_field)
            channel_windows.(frame_field) = [];
        end
        if strcmp(time_field, 'stim_time_ranges') && isfield(channel_windows, 'trial_labels') ...
                && numel(channel_windows.trial_labels) == numel(valid)
            channel_windows.trial_labels = channel_windows.trial_labels(valid, :);
        end
        if strcmp(time_field, 'shading_time_ranges') && isfield(channel_windows, 'shading_labels') ...
                && numel(channel_windows.shading_labels) == numel(valid)
            channel_windows.shading_labels = channel_windows.shading_labels(valid, :);
        end
    end
end
end

function onset_time = resolve_flash_onset_time(stim_windows, channel_name)
onset_time = NaN;
if ~isstruct(stim_windows) || ~isfield(stim_windows, 'stim_type') ...
        || ~strcmpi(string(stim_windows.stim_type), "visualstim_flash") ...
        || ~isfield(stim_windows, channel_name)
    return;
end
channel_windows = stim_windows.(channel_name);
if isfield(channel_windows, 'stim_time_ranges') && ~isempty(channel_windows.stim_time_ranges)
    onset_time = double(channel_windows.stim_time_ranges(1, 1));
elseif isfield(channel_windows, 'stim_frames') && ~isempty(channel_windows.stim_frames)
    onset_time = double(channel_windows.stim_frames(1, 1)) / 400;
end
end

function alignment_time = resolve_cycle_alignment_time(cycle_entry, channel_name, alignment_mode)
switch string(alignment_mode)
    case "flash_onset"
        alignment_field = sprintf('%s_flash_onset_time', channel_name);
        alignment_time = cycle_entry.alignment.(alignment_field);
        if ~isfinite(alignment_time)
            error('Missing flash onset time for %s in cycle %s.', channel_name, cycle_entry.cycle_name);
        end
    otherwise
        error('Unsupported record-average alignment mode: %s', alignment_mode);
end
end

function plot_record_stage_summary(voltage_stage, calcium_stage, stim_windows, metric_name, output_png, voltage_color, calcium_color, voltage_display_scale, calcium_display_scale, voltage_title, calcium_title)
fig = figure('Color', 'w', 'Name', sprintf('Record Average %s Summary', metric_name), ...
    'Position', [80, 80, 1500, 900]);
tiledlayout(fig, 2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

ax_voltage = nexttile;
[voltage_ylabel, voltage_scale_suffix] = describe_record_stage_axis_rec(metric_name);
plot_record_stage_stack( ...
    ax_voltage, voltage_stage.time, voltage_display_scale * voltage_stage.per_cycle, ...
    voltage_display_scale * voltage_stage.average, voltage_color, ...
    voltage_title, resolve_optional_channel_windows(stim_windows, 'voltage'), ...
    voltage_ylabel, voltage_scale_suffix);

ax_calcium = nexttile;
[calcium_ylabel, calcium_scale_suffix] = describe_record_stage_axis_rec(metric_name);
plot_record_stage_stack( ...
    ax_calcium, calcium_stage.time, calcium_display_scale * calcium_stage.per_cycle, ...
    calcium_display_scale * calcium_stage.average, calcium_color, ...
    calcium_title, resolve_optional_channel_windows(stim_windows, 'calcium'), ...
    calcium_ylabel, calcium_scale_suffix);

sgtitle(fig, sprintf('Record-Average %s Summary | %d cycles | light = each cycle, bright = cycle average', ...
    metric_name, voltage_stage.ncycles));
save_figure_bundle_rec(fig, strrep(output_png, '.png', '.fig'), output_png);
close(fig);
end

function plot_record_stage_stack(ax, t_axis, per_cycle, average_trace, line_color, title_text, channel_windows, y_axis_label, y_scale_suffix)
[nframes, nrois, ncycles] = size(per_cycle);
hold(ax, 'on');

stack_spacing = compute_stack_spacing_rec(average_trace);
light_color = lighten_color(line_color, 0.72);
x_min = min(t_axis);
x_max = max(t_axis);
x_span = max(eps, x_max - x_min);
left_margin = 0.14 * x_span;
ymin_global = inf;
ymax_global = -inf;

for roi_idx = 1:nrois
    offset = (nrois - roi_idx) * stack_spacing;
    current_cycle = squeeze(per_cycle(:, roi_idx, :));
    ymin_global = min(ymin_global, min(current_cycle(:), [], 'omitnan') + offset);
    ymax_global = max(ymax_global, max(current_cycle(:), [], 'omitnan') + offset);
    ymin_global = min(ymin_global, min(average_trace(:, roi_idx), [], 'omitnan') + offset);
    ymax_global = max(ymax_global, max(average_trace(:, roi_idx), [], 'omitnan') + offset);
end

if ~isfinite(ymin_global) || ~isfinite(ymax_global) || ymax_global <= ymin_global
    ymin_global = -1;
    ymax_global = 1;
else
    ypad = 0.06 * (ymax_global - ymin_global);
    ymin_global = ymin_global - ypad;
    ymax_global = ymax_global + ypad;
end

xlim(ax, [x_min - left_margin, x_max]);
ylim(ax, [ymin_global, ymax_global]);
if ~isempty(fieldnames(channel_windows))
    add_stim_shading_rec(ax, channel_windows, 0.12);
end

for roi_idx = 1:nrois
    offset = (nrois - roi_idx) * stack_spacing;
    for cycle_idx = 1:ncycles
        plot(ax, t_axis, per_cycle(:, roi_idx, cycle_idx) + offset, ...
            'Color', light_color, 'LineWidth', 0.5);
    end
    plot(ax, t_axis, average_trace(:, roi_idx) + offset, ...
        'Color', line_color, 'LineWidth', 1.4);
    y_label = offset + mean(average_trace(:, roi_idx), 'omitnan');
    if ~isfinite(y_label)
        y_label = offset;
    end
    text(ax, x_min - 0.92 * left_margin, y_label, sprintf('ROI %d', roi_idx), ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
        'FontSize', 8, 'FontWeight', 'bold', 'Color', [0.15 0.15 0.15], ...
        'Interpreter', 'none');
end

title(ax, title_text);
xlabel(ax, 'Time (s)');
ylabel(ax, y_axis_label);
set(ax, 'YTick', []);
grid(ax, 'on');
box(ax, 'off');
add_axis_scalebar_rec(ax, t_axis, average_trace(:), [0.1 0.1 0.1], y_scale_suffix);
end

function [y_axis_label, y_scale_suffix] = describe_record_stage_axis_rec(metric_name)
metric_name = lower(string(metric_name));
switch metric_name
    case "sensitivity"
        y_axis_label = 'Sensitivity (signed)';
        y_scale_suffix = 'Sensitivity';
    case "snr"
        y_axis_label = 'SNR (signed)';
        y_scale_suffix = 'SNR';
    otherwise
        y_axis_label = 'Signal (stacked display)';
        y_scale_suffix = 'a.u.';
end
end

function spacing = compute_stack_spacing_rec(trace_matrix)
trace_matrix = double(trace_matrix);
trace_std = std(trace_matrix, 0, 1, 'omitnan');
spacing = mean(trace_std, 'omitnan') * 6;
if ~isfinite(spacing) || spacing <= 0
    trace_range = max(trace_matrix, [], 'all', 'omitnan') - min(trace_matrix, [], 'all', 'omitnan');
    spacing = max(1, trace_range * 0.5);
end
end

function plot_record_metric_overlap(voltage_stage, calcium_stage, stim_windows, metric_name, output_png, voltage_color, calcium_color, voltage_polarity, calcium_polarity)
nrois = voltage_stage.nrois;
fig = figure('Color', 'w', 'Name', sprintf('Record Average %s Overlap', metric_name), ...
    'Position', [80, 80, 1600, max(500, 220 * nrois)]);
tiledlayout(fig, nrois, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

voltage_per_cycle = voltage_polarity * voltage_stage.per_cycle;
voltage_average = voltage_polarity * voltage_stage.average;
calcium_per_cycle = calcium_polarity * calcium_stage.per_cycle;
calcium_average = calcium_polarity * calcium_stage.average;

for roi_idx = 1:nrois
    ax_native = nexttile;
    plot_native_cycle_overlap(ax_native, ...
        voltage_stage.time, squeeze(voltage_per_cycle(:, roi_idx, :)), voltage_average(:, roi_idx), ...
        calcium_stage.time, squeeze(calcium_per_cycle(:, roi_idx, :)), calcium_average(:, roi_idx), ...
        voltage_color, calcium_color, ...
        resolve_optional_channel_windows(stim_windows, 'voltage'), ...
        sprintf('ROI %d | native traces', roi_idx));

    ax_norm = nexttile;
    plot_normalized_average_overlap(ax_norm, ...
        voltage_stage.time, voltage_average(:, roi_idx), ...
        calcium_stage.time, calcium_average(:, roi_idx), ...
        voltage_color, calcium_color, ...
        resolve_optional_channel_windows(stim_windows, 'voltage'), ...
        sprintf('ROI %d | average-only normalized overlap', roi_idx));
end

sgtitle(fig, sprintf(['Record-Average %s Overlap | left: each cycle in light color + bright average | ' ...
    'right: normalize only after averaging across cycles'], metric_name));
save_figure_bundle_rec(fig, strrep(output_png, '.png', '.fig'), output_png);
close(fig);
end

function plot_native_cycle_overlap(ax, t_voltage, voltage_cycles, voltage_average, t_calcium, calcium_cycles, calcium_average, voltage_color, calcium_color, channel_windows, title_text)
hold(ax, 'on');

voltage_light = lighten_color(voltage_color, 0.72);
calcium_light = lighten_color(calcium_color, 0.72);
xlimit = [min([t_voltage(:); t_calcium(:)]), max([t_voltage(:); t_calcium(:)])];
[voltage_ymin, voltage_ymax] = compute_trace_limits_rec([voltage_cycles, voltage_average]);
[calcium_ymin, calcium_ymax] = compute_trace_limits_rec([calcium_cycles, calcium_average]);

yyaxis(ax, 'left');
ylim(ax, [voltage_ymin, voltage_ymax]);
xlim(ax, xlimit);
if ~isempty(fieldnames(channel_windows))
    add_stim_shading_rec(ax, channel_windows, 0.12);
end
for idx = 1:size(voltage_cycles, 2)
    plot(ax, t_voltage, voltage_cycles(:, idx), 'Color', voltage_light, 'LineWidth', 0.5);
end
plot(ax, t_voltage, voltage_average, 'Color', voltage_color, 'LineWidth', 1.5);
ax.YColor = voltage_color;
ylabel(ax, 'Voltage');

yyaxis(ax, 'right');
ylim(ax, [calcium_ymin, calcium_ymax]);
for idx = 1:size(calcium_cycles, 2)
    plot(ax, t_calcium, calcium_cycles(:, idx), 'Color', calcium_light, 'LineWidth', 0.5);
end
plot(ax, t_calcium, calcium_average, 'Color', calcium_color, 'LineWidth', 1.5);
ax.YColor = calcium_color;
ylabel(ax, 'Calcium');

yyaxis(ax, 'left');
xlim(ax, xlimit);
xlabel(ax, 'Time (s)');
title(ax, title_text, 'FontSize', 9, 'FontWeight', 'bold');
grid(ax, 'on');
box(ax, 'off');
end

function plot_normalized_average_overlap(ax, t_voltage, voltage_average, t_calcium, calcium_average, voltage_color, calcium_color, channel_windows, title_text)
hold(ax, 'on');
xlim(ax, [min([t_voltage(:); t_calcium(:)]), max([t_voltage(:); t_calcium(:)])]);
ylim(ax, [0, 1]);
if ~isempty(fieldnames(channel_windows))
    add_stim_shading_rec(ax, channel_windows, 0.12);
end

plot(ax, t_voltage, normalize_trace_to_unit_range(voltage_average), ...
    'Color', voltage_color, 'LineWidth', 1.5);
plot(ax, t_calcium, normalize_trace_to_unit_range(calcium_average), ...
    'Color', calcium_color, 'LineWidth', 1.5);

xlabel(ax, 'Time (s)');
ylabel(ax, 'Normalized average');
title(ax, title_text, 'FontSize', 9, 'FontWeight', 'bold');
legend(ax, {'Voltage avg', 'Calcium avg'}, 'Location', 'best');
grid(ax, 'on');
box(ax, 'off');
end

function x_norm = normalize_trace_to_unit_range(x)
x = double(x(:));
valid = isfinite(x);
if ~any(valid)
    x_norm = zeros(size(x));
    return;
end
x_min = min(x(valid));
x_max = max(x(valid));
if abs(x_max - x_min) < eps
    x_norm = zeros(size(x));
else
    x_norm = (x - x_min) ./ (x_max - x_min);
end
end

function color_out = lighten_color(color_in, amount)
color_out = color_in + (1 - color_in) * amount;
color_out = min(max(color_out, 0), 1);
end

function channel_windows = resolve_optional_channel_windows(stim_windows, channel_name)
channel_windows = struct();
if isstruct(stim_windows) && isfield(stim_windows, channel_name)
    channel_windows = stim_windows.(channel_name);
end
end

function add_stim_shading_rec(ax, channel_windows, alpha_value)
if isempty(channel_windows) || ~isstruct(channel_windows)
    return;
end

if isfield(channel_windows, 'shading_time_ranges') && ~isempty(channel_windows.shading_time_ranges) ...
        && isfield(channel_windows, 'shading_labels') && numel(channel_windows.shading_labels) == size(channel_windows.shading_time_ranges, 1)
    time_ranges = channel_windows.shading_time_ranges;
    labels = string(channel_windows.shading_labels(:));
elseif isfield(channel_windows, 'block_time_ranges') && ~isempty(channel_windows.block_time_ranges) ...
        && isfield(channel_windows, 'block_labels') && numel(channel_windows.block_labels) == size(channel_windows.block_time_ranges, 1)
    time_ranges = channel_windows.block_time_ranges;
    labels = string(channel_windows.block_labels(:));
elseif isfield(channel_windows, 'stim_time_ranges') && ~isempty(channel_windows.stim_time_ranges)
    time_ranges = channel_windows.stim_time_ranges;
    if isfield(channel_windows, 'trial_labels')
        labels = string(channel_windows.trial_labels(:));
    else
        labels = repmat("stim", size(time_ranges, 1), 1);
    end
else
    return;
end

yl = ylim(ax);
for idx = 1:size(time_ranges, 1)
    if any(~isfinite(time_ranges(idx, :)))
        continue;
    end
    [shade_color, shade_alpha] = resolve_stim_shading_style_rec(labels(idx), [0.7 0.7 0.7], alpha_value);
    patch(ax, ...
        [time_ranges(idx, 1), time_ranges(idx, 2), time_ranges(idx, 2), time_ranges(idx, 1)], ...
        [yl(1), yl(1), yl(2), yl(2)], ...
        shade_color, 'FaceAlpha', shade_alpha, 'EdgeColor', 'none');
end
ylim(ax, yl);
end

function [ymin, ymax] = compute_trace_limits_rec(trace_matrix)
trace_values = double(trace_matrix(:));
trace_values = trace_values(isfinite(trace_values));
if isempty(trace_values)
    ymin = -1;
    ymax = 1;
    return;
end
ymin = min(trace_values);
ymax = max(trace_values);
if ymax <= ymin
    pad = max(1e-3, abs(ymax) * 0.1 + 1e-3);
else
    pad = 0.08 * (ymax - ymin);
end
ymin = ymin - pad;
ymax = ymax + pad;
end

function [shade_color, shade_alpha] = resolve_stim_shading_style_rec(label, default_color, default_alpha)
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

function time_frequency = build_record_average_time_frequency(tf_input, stim_windows, output_dir)
params = struct( ...
    'voltage_stage', string(tf_input.voltage.stage_name), ...
    'calcium_stage', string(tf_input.calcium.stage_name), ...
    'min_freq_hz', 0.5, ...
    'max_freq_hz', min([80, tf_input.voltage.frame_rate / 2 - eps, tf_input.calcium.frame_rate / 2 - eps]), ...
    'wavelet_name', "amor");

voltage_tf = analyze_population_time_frequency_rec( ...
    tf_input.voltage.trace, tf_input.voltage.frame_rate, tf_input.voltage.time, params, 'Voltage');
calcium_tf = analyze_population_time_frequency_rec( ...
    tf_input.calcium.trace, tf_input.calcium.frame_rate, tf_input.calcium.time, params, 'Calcium');

[fourier_fig, fourier_png] = plot_population_fourier_summary_rec( ...
    voltage_tf, calcium_tf, stim_windows, params, output_dir);
[wavelet_fig, wavelet_png] = plot_population_wavelet_summary_rec( ...
    voltage_tf, calcium_tf, stim_windows, output_dir);

time_frequency = struct( ...
    'parameters', params, ...
    'voltage', voltage_tf, ...
    'calcium', calcium_tf, ...
    'visualizations', struct( ...
        'fourier_summary_fig', string(fourier_fig), ...
        'fourier_summary_png', string(fourier_png), ...
        'wavelet_summary_fig', string(wavelet_fig), ...
        'wavelet_summary_png', string(wavelet_png)));
end

function result = analyze_population_time_frequency_rec(trace_matrix, frame_rate, t_axis, params, channel_name)
trace_matrix = double(trace_matrix);
t_axis = double(t_axis(:));
nrois = size(trace_matrix, 2);
spectrum_frequency = [];
spectrum_amplitude = [];
peak_frequency_hz = NaN(nrois, 1);
peak_amplitude = NaN(nrois, 1);
signal_rms = NaN(nrois, 1);

for roi_idx = 1:nrois
    x = sanitize_trace_for_spectrum_rec(trace_matrix(:, roi_idx));
    signal_rms(roi_idx) = rms(x);
    [frequency_i, amplitude_i] = compute_trace_fft_spectrum_rec(x, frame_rate);
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
        'trace', sanitize_trace_for_spectrum_rec(trace_matrix(:, representative_roi)), ...
        'peak_frequency_hz', peak_frequency_hz(representative_roi), ...
        'wavelet', compute_wavelet_scalogram_rec(trace_matrix(:, representative_roi), frame_rate, t_axis, params)), ...
    'summary', struct( ...
        'nrois', nrois, ...
        'median_peak_frequency_hz', median(peak_frequency_hz, 'omitnan'), ...
        'median_peak_amplitude', median(peak_amplitude, 'omitnan')));
end

function x = sanitize_trace_for_spectrum_rec(x)
x = double(x(:));
if isempty(x) || all(~isfinite(x))
    x = zeros(size(x));
    return;
end
x = fillmissing(x, 'linear', 'EndValues', 'nearest');
x = x - mean(x, 'omitnan');
end

function [frequency, amplitude] = compute_trace_fft_spectrum_rec(x, frame_rate)
x = sanitize_trace_for_spectrum_rec(x);
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

function wavelet_data = compute_wavelet_scalogram_rec(x, frame_rate, t_axis, params)
x = sanitize_trace_for_spectrum_rec(x);
[wt, frequency] = cwt(x, frame_rate);
valid_mask = frequency >= params.min_freq_hz & frequency <= params.max_freq_hz;
wavelet_data = struct( ...
    'time', double(t_axis(:)), ...
    'frequency', frequency(valid_mask), ...
    'power', abs(wt(valid_mask, :)).^2);
end

function [fig_file, png_file] = plot_population_fourier_summary_rec(voltage_tf, calcium_tf, stim_windows, params, output_dir)
fig = figure('Color', 'w', 'Name', 'Record-Average FFT Summary');
tiledlayout(fig, 4, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
plot_channel_fft_summary_rec(voltage_tf, 'r', 'Voltage FFT Amplitude', [0, params.max_freq_hz], false);
nexttile;
plot_channel_fft_summary_rec(calcium_tf, 'g', 'Calcium FFT Amplitude', [0, params.max_freq_hz], false);
nexttile;
plot_channel_fft_summary_rec(voltage_tf, 'r', 'Voltage FFT Amplitude', [0, 10], true);
nexttile;
plot_channel_fft_summary_rec(calcium_tf, 'g', 'Calcium FFT Amplitude', [0, 10], true);
nexttile;
plot_flash_fft_comparison_rec(voltage_tf, resolve_optional_channel_windows(stim_windows, 'voltage'), 'r', 'Voltage Flash FFT');
nexttile;
plot_flash_fft_comparison_rec(calcium_tf, resolve_optional_channel_windows(stim_windows, 'calcium'), 'g', 'Calcium Flash FFT');
nexttile;
plot_roi_peak_summary_rec(voltage_tf, [0, 10], 'Voltage ROI Peak Frequency (0-10 Hz)');
nexttile;
plot_roi_peak_summary_rec(calcium_tf, [0, 10], 'Calcium ROI Peak Frequency (0-10 Hz)');

fig_file = fullfile(output_dir, '6_record_average_fourier_summary.fig');
png_file = fullfile(output_dir, '6_record_average_fourier_summary.png');
save_figure_bundle_rec(fig, fig_file, png_file);
close(fig);
end

function [fig_file, png_file] = plot_population_wavelet_summary_rec(voltage_tf, calcium_tf, stim_windows, output_dir)
fig = figure('Color', 'w', 'Name', 'Record-Average Wavelet Summary');
nrois = size(voltage_tf.trace_matrix, 2);
tiledlayout(fig, nrois, 4, 'Padding', 'compact', 'TileSpacing', 'compact');

for roi_idx = 1:nrois
    nexttile;
    plot_roi_trace_with_stim_rec(gca, voltage_tf, roi_idx, resolve_optional_channel_windows(stim_windows, 'voltage'), 'r', 'Voltage Trace');
    nexttile;
    plot_roi_wavelet_scalogram_rec(gca, voltage_tf, roi_idx, resolve_optional_channel_windows(stim_windows, 'voltage'), 'Voltage Wavelet');
    nexttile;
    plot_roi_trace_with_stim_rec(gca, calcium_tf, roi_idx, resolve_optional_channel_windows(stim_windows, 'calcium'), 'g', 'Calcium Trace');
    nexttile;
    plot_roi_wavelet_scalogram_rec(gca, calcium_tf, roi_idx, resolve_optional_channel_windows(stim_windows, 'calcium'), 'Calcium Wavelet');
end

fig_file = fullfile(output_dir, '7_record_average_wavelet_summary.fig');
png_file = fullfile(output_dir, '7_record_average_wavelet_summary.png');
save_figure_bundle_rec(fig, fig_file, png_file);
close(fig);
end

function plot_channel_fft_summary_rec(channel_result, trace_color, title_text, display_band_hz, annotate_peak)
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
    [peak_hz, peak_amp] = find_band_peak_from_amplitude_rec(frequency, mean_amplitude, display_band_hz);
    if isfinite(peak_hz)
        xline(peak_hz, '--', sprintf('%.2f Hz', peak_hz), ...
            'Color', trace_color, 'LineWidth', 1.2, 'LabelVerticalAlignment', 'middle');
        scatter(peak_hz, peak_amp, 42, trace_color, 'filled');
    end
end
xlabel('Frequency (Hz)');
ylabel('Single-Sided Amplitude');
title(sprintf('%s | %s', title_text, ternary_rec(annotate_peak, '0-10 Hz zoom', 'full range')));
xlim(display_band_hz);
grid on;
end

function plot_flash_fft_comparison_rec(channel_result, channel_windows, trace_color, title_text)
[baseline_trace, stim_trace, supported] = extract_flash_average_segments_rec(channel_result, channel_windows);
if ~supported
    text(0.5, 0.5, 'Flash pre/post FFT unavailable', 'Units', 'normalized', 'HorizontalAlignment', 'center');
    axis off;
    return;
end

[f_base, a_base] = compute_trace_fft_spectrum_rec(baseline_trace, channel_result.frame_rate);
[f_stim, a_stim] = compute_trace_fft_spectrum_rec(stim_trace, channel_result.frame_rate);
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

function [baseline_trace, stim_trace, supported] = extract_flash_average_segments_rec(channel_result, channel_windows)
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

function trials = extract_aligned_trials_rec(trace, frame_ranges)
trace = double(trace(:));
frame_ranges = round(frame_ranges);
trial_lengths = frame_ranges(:, 2) - frame_ranges(:, 1) + 1;
valid_lengths = trial_lengths(isfinite(trial_lengths) & trial_lengths > 1);
if isempty(valid_lengths)
    trials = [];
    return;
end
target_length = min(valid_lengths);
trials = NaN(target_length, size(frame_ranges, 1));
count = 0;
for idx = 1:size(frame_ranges, 1)
    start_idx = max(1, frame_ranges(idx, 1));
    end_idx = min(numel(trace), frame_ranges(idx, 2));
    if end_idx - start_idx + 1 < target_length
        continue;
    end
    count = count + 1;
    trials(:, count) = trace(start_idx:start_idx + target_length - 1);
end
trials = trials(:, 1:count);
end

function plot_roi_peak_summary_rec(channel_result, display_band_hz, title_text)
roi_count = size(channel_result.fft.amplitude, 2);
peak_frequency_hz = NaN(roi_count, 1);
peak_amplitude = NaN(roi_count, 1);
for roi_idx = 1:roi_count
    [peak_frequency_hz(roi_idx), peak_amplitude(roi_idx)] = find_band_peak_from_amplitude_rec( ...
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

function plot_representative_trace_with_stim_rec(ax, channel_result, channel_windows, trace_color, title_text)
t = channel_result.time(:);
x = channel_result.representative.trace(:);
[y_label, scale_suffix, stage_label] = resolve_time_frequency_trace_label_rec(channel_result);
plot(ax, t, x, 'Color', trace_color, 'LineWidth', 1.2);
apply_trace_axis_limits_rec(ax, x);
if ~isempty(fieldnames(channel_windows))
    add_stim_shading_rec(ax, channel_windows, 0.16);
    hold(ax, 'on');
end
plot(ax, t, x, 'Color', trace_color, 'LineWidth', 1.2);
xlabel(ax, 'Time (s)');
ylabel(ax, y_label);
title(ax, sprintf('%s (%s) | ROI %d | peak %.2f Hz', ...
    title_text, stage_label, channel_result.representative.roi_index, channel_result.representative.peak_frequency_hz));
grid(ax, 'on');
add_axis_scalebar_rec(ax, t, x, trace_color, scale_suffix);
end

function plot_roi_trace_with_stim_rec(ax, channel_result, roi_idx, channel_windows, trace_color, title_text)
t = channel_result.time(:);
x = double(channel_result.trace_matrix(:, roi_idx));
[y_label, scale_suffix, stage_label] = resolve_time_frequency_trace_label_rec(channel_result);
plot(ax, t, x, 'Color', trace_color, 'LineWidth', 1.0);
apply_trace_axis_limits_rec(ax, x);
if ~isempty(fieldnames(channel_windows))
    add_stim_shading_rec(ax, channel_windows, 0.16);
    hold(ax, 'on');
end
plot(ax, t, x, 'Color', trace_color, 'LineWidth', 1.0);
xlabel(ax, 'Time (s)');
ylabel(ax, y_label);
title(ax, sprintf('%s (%s) | ROI %d | peak %.2f Hz', ...
    title_text, stage_label, roi_idx, channel_result.roi.peak_frequency_hz(roi_idx)));
grid(ax, 'on');
add_axis_scalebar_rec(ax, t, x, trace_color, scale_suffix);
end

function [y_label, scale_suffix, stage_label] = resolve_time_frequency_trace_label_rec(channel_result)
stage_name = "";
if isfield(channel_result, 'analysis_parameters') && isstruct(channel_result.analysis_parameters)
    if strcmpi(string(channel_result.channel_name), "Voltage") ...
            && isfield(channel_result.analysis_parameters, 'voltage_stage')
        stage_name = string(channel_result.analysis_parameters.voltage_stage);
    elseif strcmpi(string(channel_result.channel_name), "Calcium") ...
            && isfield(channel_result.analysis_parameters, 'calcium_stage')
        stage_name = string(channel_result.analysis_parameters.calcium_stage);
    end
end
[y_label, scale_suffix, stage_label] = describe_trace_stage_for_display_rec(stage_name);
end

function [y_label, scale_suffix, stage_label] = describe_trace_stage_for_display_rec(stage_name)
stage_name = lower(string(stage_name));
switch stage_name
    case {"snr", "snr_average"}
        y_label = 'SNR (signed)';
        scale_suffix = 'SNR';
        stage_label = 'SNR';
    case {"sensitivity", "sensitivity_average"}
        y_label = 'Sensitivity (signed)';
        scale_suffix = 'Sensitivity';
        stage_label = 'Sensitivity';
    case {"bleach_removed", "bleach_removed_average"}
        y_label = 'Bleach-Removed Signal';
        scale_suffix = 'a.u.';
        stage_label = 'bleach removed';
    case {"bg_removed", "bg_removed_average"}
        y_label = 'BG-Removed Signal';
        scale_suffix = 'a.u.';
        stage_label = 'bg removed';
    otherwise
        y_label = 'Trace Value';
        scale_suffix = 'a.u.';
        if strlength(stage_name) == 0
            stage_label = 'trace';
        else
            stage_label = char(strrep(stage_name, "_", " "));
        end
    end
end

function add_axis_scalebar_rec(ax, x_values, y_values, bar_color, y_suffix)
if nargin < 5
    y_suffix = '';
end

x_values = double(x_values(:));
y_values = double(y_values(:));
y_values = y_values(isfinite(y_values));
if isempty(x_values) || isempty(y_values)
    return;
end

x_limits = xlim(ax);
y_limits = ylim(ax);
x_span = max(eps, diff(x_limits));
y_span = max(eps, diff(y_limits));

x_len = nice_scalebar_value_rec(0.16 * x_span);
y_len = nice_scalebar_value_rec(0.18 * range(y_values));
if ~isfinite(y_len) || y_len <= 0
    y_len = nice_scalebar_value_rec(0.18 * y_span);
end

x_start = x_limits(2) - 0.06 * x_span - x_len;
y_start = y_limits(1) + 0.10 * y_span;

plot(ax, [x_start, x_start + x_len], [y_start, y_start], ...
    'Color', bar_color, 'LineWidth', 1.4, 'Clipping', 'off');
plot(ax, [x_start, x_start], [y_start, y_start + y_len], ...
    'Color', bar_color, 'LineWidth', 1.4, 'Clipping', 'off');

text(ax, x_start + x_len / 2, y_start - 0.04 * y_span, ...
    sprintf('%s s', format_scalebar_value_rec(x_len)), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
    'FontSize', 8, 'Color', bar_color, 'Interpreter', 'none');

if strlength(string(y_suffix)) > 0
    y_label_local = sprintf('%s %s', format_scalebar_value_rec(y_len), y_suffix);
else
    y_label_local = format_scalebar_value_rec(y_len);
end
text(ax, x_start - 0.01 * x_span, y_start + y_len / 2, y_label_local, ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
    'Rotation', 90, 'FontSize', 8, 'Color', bar_color, 'Interpreter', 'none');
end

function value = nice_scalebar_value_rec(target_value)
target_value = abs(double(target_value));
if ~isfinite(target_value) || target_value <= 0
    value = 1;
    return;
end

base = 10 ^ floor(log10(target_value));
candidates = [1, 2, 5, 10] * base;
idx = find(candidates >= target_value, 1, 'first');
if isempty(idx)
    value = candidates(end);
else
    value = candidates(idx);
end
end

function text_value = format_scalebar_value_rec(value)
value = double(value);
if abs(value) >= 10
    text_value = sprintf('%.0f', value);
elseif abs(value) >= 1
    text_value = sprintf('%.1f', value);
elseif abs(value) >= 0.1
    text_value = sprintf('%.2f', value);
else
    text_value = sprintf('%.3f', value);
end
end

function apply_trace_axis_limits_rec(ax, x)
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

function plot_wavelet_scalogram_rec(ax, channel_result, channel_windows, title_text)
spec = channel_result.representative.wavelet;
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
title(ax, sprintf('%s | ROI %d', title_text, channel_result.representative.roi_index));
cb = colorbar(ax);
cb.Label.String = 'Wavelet Power';
hold(ax, 'on');
if isfinite(channel_result.representative.peak_frequency_hz)
    yline(ax, channel_result.representative.peak_frequency_hz, 'w--', 'LineWidth', 1.0);
end
ylim(ax, [max(min(frequency(frequency > 0)), 0.5), max(frequency)]);
yticks(ax, [1 3 10 30 80]);
overlay_stim_boundaries_rec(ax, channel_windows);
end

function plot_roi_wavelet_scalogram_rec(ax, channel_result, roi_idx, channel_windows, title_text)
trace = double(channel_result.trace_matrix(:, roi_idx));
params = channel_result.analysis_parameters;
spec = compute_wavelet_scalogram_rec(trace, channel_result.frame_rate, channel_result.time, params);
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
title(ax, sprintf('%s | ROI %d | peak %.2f Hz', ...
    title_text, roi_idx, channel_result.roi.peak_frequency_hz(roi_idx)));
cb = colorbar(ax);
cb.Label.String = 'Wavelet Power';
hold(ax, 'on');
if isfinite(channel_result.roi.peak_frequency_hz(roi_idx))
    yline(ax, channel_result.roi.peak_frequency_hz(roi_idx), 'w--', 'LineWidth', 1.0);
end
ylim(ax, [max(min(frequency(frequency > 0)), 0.5), max(frequency)]);
yticks(ax, [1 3 10 30 80]);
overlay_stim_boundaries_rec(ax, channel_windows);
end

function overlay_stim_boundaries_rec(ax, channel_windows)
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

function [peak_frequency_hz, peak_amplitude] = find_band_peak_from_amplitude_rec(frequency, amplitude, display_band_hz)
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

function save_figure_bundle_rec(fig_handle, fig_file, png_file)
% Save the editable FIG plus a full-screen PNG snapshot so the exported
% record-level summary is easier to inspect without reopening MATLAB.
prepare_figure_for_png_export_rec(fig_handle);
saveas(fig_handle, fig_file, 'fig');
exportgraphics(fig_handle, png_file, 'Resolution', 150);
end

function prepare_figure_for_png_export_rec(fig_handle)
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

function out = ternary_rec(condition, true_value, false_value)
if condition
    out = true_value;
else
    out = false_value;
end
end
