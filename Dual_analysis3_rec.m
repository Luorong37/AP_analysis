% Dual_analysis3_rec - Run Dual_analysis3 across all cycles in one record
%
% Workflow:
% 1. Choose one reference cycle.
% 2. Create or reuse one ROI file from that reference cycle.
% 3. Reuse the same ROI file for every other Cycle* under the same record.
%
% This script keeps Dual_analysis3 as a normal runnable script.
% It passes parameters in from an outer scope by calling run(...).

clc;

%% Input Setup
% Point this to one rebuilt record folder that contains Cycle* subfolders.
rec_path = 'V:\Luorong\Invivo\26.04.15_dual-color_P195\Methods5_default\Rec1_2026-04-15_20-58-57';

% Choose one cycle as the ROI reference.
reference_cycle_name = 'Cycle1';

% If you already have a good ROI file, fill it here and the script will
% skip the reference ROI selection step.
reference_roi_file = '';

% Reuse the latest ROI result already present in the reference cycle when
% possible. If false, the reference cycle will be rerun to generate ROI.
reuse_existing_reference_roi = true;

% When true, rerun the reference cycle first to create a fresh ROI file.
rerun_reference_cycle_for_roi = false;

% After creating/reusing the reference ROI, also run the full analysis on
% the reference cycle inside the batch loop.
run_reference_cycle_in_batch = false;

% Leave empty to run every Cycle*. Otherwise provide a string array like:
% ["Cycle1","Cycle3","Cycle5"]
cycle_name_filter = strings(0, 1);

% Skip a cycle if it already contains any previous Dual_analysis3 result.
skip_cycles_with_existing_results = false;

% Shared Dual_analysis3 overrides.
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

map_bin = 4;
calcium_smoothing_window = 40;
downsample_window = calcium_smoothing_window;
bleach_mode_voltage = 'linear';
bleach_mode_calcium = 'exp2';
voltage_polarity = -1;
calcium_polarity = 1;
reuse_offset = [];
gpu = true;
motion_cfg = struct( ...
    'enabled', true, ...
    'use_saved_shift', false, ...
    'saved_shift_file', '', ...
    'highpass', true);

%% Resolve Cycles
script_dir = fileparts(mfilename('fullpath'));
dual_script_path = fullfile(script_dir, 'Dual_analysis3.m');
if ~isfile(dual_script_path)
    error('Cannot find Dual_analysis3.m next to this batch script.');
end

if ~isfolder(rec_path)
    error('Record path does not exist: %s', rec_path);
end

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
        dual_script_path, reference_cycle_path, "", true, ...
        camera_cfg, map_bin, calcium_smoothing_window, downsample_window, ...
        bleach_mode_voltage, bleach_mode_calcium, ...
        voltage_polarity, calcium_polarity, ...
        reuse_offset, gpu, motion_cfg);
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
        existing_result = find_latest_explicit_result(current_cycle_path);
        if strlength(existing_result) > 0
            batch_results(end+1, 1) = struct( ...
                'cycle_name', current_cycle_name, ...
                'cycle_path', string(current_cycle_path), ...
                'status', "skipped_existing_result", ...
                'save_path', "", ...
                'roi_file', reference_roi_file, ...
                'message', existing_result);
            fprintf('Skipping %s because an existing Dual_analysis3 result was found.\n', current_cycle_name);
            continue;
        end
    end

    fprintf('\n[Batch] Running %s ...\n', current_cycle_name);
    try
        run_result = run_dual_cycle_with_overrides( ...
            dual_script_path, current_cycle_path, reference_roi_file, false, ...
            camera_cfg, map_bin, calcium_smoothing_window, downsample_window, ...
            bleach_mode_voltage, bleach_mode_calcium, ...
            voltage_polarity, calcium_polarity, ...
            reuse_offset, gpu, motion_cfg);

        batch_results(end+1, 1) = struct( ...
            'cycle_name', current_cycle_name, ...
            'cycle_path', string(current_cycle_path), ...
            'status', "completed", ...
            'save_path', string(run_result.save_path), ...
            'roi_file', string(run_result.roi_file), ...
            'message', "");
    catch ME
        batch_results(end+1, 1) = struct( ...
            'cycle_name', current_cycle_name, ...
            'cycle_path', string(current_cycle_path), ...
            'status', "failed", ...
            'save_path', "", ...
            'roi_file', reference_roi_file, ...
            'message', string(ME.message));
        fprintf(2, '[Batch] %s failed: %s\n', current_cycle_name, ME.message);
    end
end

batch_summary_file = fullfile(rec_path, 'Dual_analysis3_rec_batch_summary.mat');
save(batch_summary_file, 'batch_results', 'reference_roi_file');
fprintf('\nBatch summary saved to:\n  %s\n', batch_summary_file);

%% Local Functions
function result = run_dual_cycle_with_overrides( ...
    dual_script_path, cycle_path, reuse_roi_file, correct_offset, ...
    camera_cfg, map_bin, calcium_smoothing_window, downsample_window, ...
    bleach_mode_voltage, bleach_mode_calcium, ...
    voltage_polarity, calcium_polarity, ...
    reuse_offset, gpu, motion_cfg)

% All variables defined in this function are visible to Dual_analysis3.m
% when it is executed via run(...). This keeps Dual_analysis3 as a script
% while still letting an outer batch driver provide overrides.
[record_path_for_name, cycle_name_for_name] = fileparts(cycle_path);
[~, record_name_for_name] = fileparts(record_path_for_name);
time_tag = string(datetime('now', 'Format', 'yyyy-MM-dd HH-mm-ss-SSS'));
analysis_run_name = sprintf('%s_%s_%s', record_name_for_name, cycle_name_for_name, char(time_tag));

run(dual_script_path);

result = struct( ...
    'save_path', save_path, ...
    'roi_file', fullfile(save_path, '1_dual_roi_results.mat'));
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

function roi_file = find_latest_dual_roi_file(cycle_path)
roi_file = "";
listing = dir(fullfile(cycle_path, 'Dual_analysis3', '**', '1_dual_roi_results.mat'));
if isempty(listing)
    return;
end
[~, newest_idx] = max([listing.datenum]);
roi_file = string(fullfile(listing(newest_idx).folder, listing(newest_idx).name));
end

function explicit_result = find_latest_explicit_result(cycle_path)
explicit_result = "";
listing = dir(fullfile(cycle_path, 'Dual_analysis3', '**', '-1_explicit_dual_results.mat'));
if isempty(listing)
    return;
end
[~, newest_idx] = max([listing.datenum]);
explicit_result = string(fullfile(listing(newest_idx).folder, listing(newest_idx).name));
end
