% Dual_analysis3_volpy - VolPy voltage re-analysis companion for Dual_analysis3
%
% This script converts VolPy output into the existing Dual_analysis3 result
% containers, then delegates all peak/plot/stim/frequency sections to
% Dual_analysis3 in analysis_only mode. The staging folder is removed after
% a successful run; the final MAT variable names and structures stay
% compatible with standard Dual_analysis3 outputs.

nowtime = string(datetime('now', 'Format', 'yyyy-MM-dd HH-mm-ss-SSS'));
fprintf('Initializing Dual_analysis3 VolPy re-analysis...\n');

if ~exist('dual3_volpy_control', 'var') || isempty(dual3_volpy_control)
    dual3_volpy_control = struct();
end
volpy_defaults = struct( ...
    'cycle_path', "", ...
    'source_results_path', "", ...
    'run_name', "", ...
    'output_path', "", ...
    'auto_run', false, ...
    'force_rerun', false, ...
    'use_existing', false, ...
    'flip_signal', false);
dual3_volpy_control = apply_volpy_defaults(dual3_volpy_control, volpy_defaults);

cycle_path = char(string(dual3_volpy_control.cycle_path));
if isempty(cycle_path)
    error('dual3_volpy_control.cycle_path is required.');
end
source_results_path = resolve_volpy_source_results_path( ...
    cycle_path, dual3_volpy_control.source_results_path);

[record_path, cycle_name] = fileparts(cycle_path);
[~, record_name] = fileparts(record_path);
run_name = char(string(dual3_volpy_control.run_name));
if isempty(run_name)
    run_name = sprintf('%s_%s_volpy_%s', record_name, cycle_name, char(nowtime));
end
save_path = char(string(dual3_volpy_control.output_path));
if isempty(save_path)
    save_path = fullfile(cycle_path, 'Dual_analysis3_volpy_voltage', run_name);
end
if isfolder(save_path) && ~isempty(non_dot_dir_entries_volpy(save_path))
    error('VolPy output folder already exists and is not empty: %s', save_path);
end

stage_path = fullfile(cycle_path, 'Dual_analysis3_volpy_voltage', ...
    [run_name '_staging']);
if isfolder(stage_path)
    error('VolPy staging folder already exists: %s', stage_path);
end
mkdir(stage_path);
stage_cleanup = onCleanup(@() cleanup_failed_volpy_stage(stage_path, save_path));

voltage_polarity_local = resolve_workspace_value_volpy('voltage_polarity', -1);
calcium_polarity_local = resolve_workspace_value_volpy('calcium_polarity', 1);
calcium_smoothing_window_local = resolve_workspace_value_volpy( ...
    'calcium_smoothing_window', 40);

[source_dual_info, source_voltage_results, source_calcium_results, ...
    source_dual_results, source_stim_results, source_roi] = ...
    load_standard_volpy_source(source_results_path);
voltage_movie_path = resolve_channel_movie_path_from_dual_info_volpy( ...
    source_dual_info, "voltage");
if ~isfile(voltage_movie_path)
    error('Voltage movie for VolPy re-analysis was not found: %s', voltage_movie_path);
end

backend_dir = fullfile(stage_path, '0_volpy_backend');
mkdir(backend_dir);
roi_mask_file = fullfile(backend_dir, 'volpy_roi_mask.mat');
mask = uint16(source_roi.rois.bwmask);
save(roi_mask_file, 'mask');

result_mat_path = fullfile(backend_dir, 'volpy_results.mat');
voltage_frame_rate = double(source_voltage_results.movie_info.frame_rate);
backend_info = run_volpy_backend_pipeline( ...
    voltage_movie_path, voltage_frame_rate, ...
    logical(dual3_volpy_control.flip_signal), ...
    logical(dual3_volpy_control.auto_run), ...
    logical(dual3_volpy_control.force_rerun), ...
    logical(dual3_volpy_control.use_existing), ...
    backend_dir, result_mat_path, roi_mask_file);
volpy_data = load_volpy_backend_results(char(backend_info.result_mat));

stage_roi_file = fullfile(stage_path, '1_dual_roi_results.mat');
copyfile(fullfile(source_results_path, '1_dual_roi_results.mat'), stage_roi_file);
source_dual_results.registration.info.source_roi_file = ...
    string(fullfile(source_results_path, '1_dual_roi_results.mat'));
source_dual_results.registration.info.roi_file = stage_roi_file;

voltage_results = build_volpy_voltage_results( ...
    source_voltage_results, volpy_data, stage_roi_file, backend_info);
calcium_results = source_calcium_results;
dual_results = source_dual_results;
stim_results = source_stim_results;
dual_info = source_dual_info;
dual_info.analysis_name = 'Dual_analysis3';
dual_info.analysis_backend = 'volpy_voltage_reanalysis';
dual_info.save_path = stage_path;
dual_info.created_at = datetime("now");
dual_info.volpy_source_results_path = source_results_path;
dual_info.volpy_backend = struct( ...
    'output_dir', string(backend_info.output_dir), ...
    'result_mat', string(backend_info.result_mat), ...
    'motion_corrected_file', string(backend_info.motion_corrected_file), ...
    'roi_mask_file', string(roi_mask_file), ...
    'flip_signal', logical(dual3_volpy_control.flip_signal), ...
    'voltage_movie_path', string(voltage_movie_path), ...
    'calcium_smoothing_window', calcium_smoothing_window_local, ...
    'voltage_polarity', voltage_polarity_local, ...
    'calcium_polarity', calcium_polarity_local);
dual_results.volpy_voltage_reanalysis = struct( ...
    'data', struct(), ...
    'info', struct( ...
        'source_results_path', string(source_results_path), ...
        'result_mat', string(backend_info.result_mat), ...
        'motion_corrected_file', string(backend_info.motion_corrected_file), ...
        'roi_mask_file', string(roi_mask_file), ...
        'created_at', datetime("now")));

save(fullfile(stage_path, 'dual_info.mat'), 'dual_info');
save(fullfile(stage_path, 'voltage_results.mat'), 'voltage_results', '-v7.3');
save(fullfile(stage_path, 'calcium_results.mat'), 'calcium_results', '-v7.3');
save(fullfile(stage_path, 'dual_results.mat'), 'dual_results', '-v7.3');
save(fullfile(stage_path, 'stim_results.mat'), 'stim_results', '-v7.3');

dual3_control = struct( ...
    'cycle_path', string(cycle_path), ...
    'preset', "analysis_only", ...
    'run_name', string(run_name), ...
    'output_path', string(save_path), ...
    'source_results_path', string(stage_path), ...
    'offset_mode', "none", ...
    'reuse_offset', [], ...
    'sections', struct('peak', "run"));

dual_script_path = fullfile(fileparts(mfilename('fullpath')), 'Dual_analysis3.m');
run(dual_script_path);

final_backend_dir = fullfile(save_path, '0_volpy_backend');
movefile(backend_dir, final_backend_dir);
update_final_volpy_paths(save_path, stage_path, final_backend_dir, ...
    source_results_path, voltage_movie_path, dual3_volpy_control, ...
    calcium_smoothing_window_local, voltage_polarity_local, calcium_polarity_local);
rmdir(stage_path, 's');
clear stage_cleanup;

fprintf('Dual_analysis3 VolPy re-analysis finished.\n');
fprintf('Results saved to: %s\n', save_path);

function control = apply_volpy_defaults(control, defaults)
if ~isstruct(control) || ~isscalar(control)
    error('dual3_volpy_control must be one scalar struct.');
end
unknown = setdiff(fieldnames(control), fieldnames(defaults));
if ~isempty(unknown)
    error('Unknown dual3_volpy_control field(s): %s', strjoin(unknown, ', '));
end
names = fieldnames(defaults);
for idx = 1:numel(names)
    name = names{idx};
    if ~isfield(control, name) || isempty(control.(name))
        control.(name) = defaults.(name);
    end
end
end

function entries = non_dot_dir_entries_volpy(folder)
entries = dir(folder);
entries = entries(~ismember({entries.name}, {'.', '..'}));
end

function cleanup_failed_volpy_stage(stage_path, save_path)
if isfolder(stage_path) && ~isfolder(save_path)
    fprintf(2, 'VolPy staging folder retained after failure: %s\n', stage_path);
end
end

function value = resolve_workspace_value_volpy(variable_name, default_value)
if evalin('caller', sprintf('exist(''%s'', ''var'')', variable_name))
    value = evalin('caller', variable_name);
else
    value = default_value;
end
end

function source_path = resolve_volpy_source_results_path(cycle_path, source_path)
source_path = string(source_path);
if strlength(source_path) > 0
    if ~isfolder(source_path)
        error('VolPy source_results_path does not exist: %s', source_path);
    end
    return;
end
listing = dir(fullfile(cycle_path, 'Dual_analysis3', '**', ...
    '-1_explicit_dual_results.mat'));
if isempty(listing)
    error('No standard Dual_analysis3 result was found for VolPy re-analysis.');
end
[~, newest_idx] = max([listing.datenum]);
source_path = string(listing(newest_idx).folder);
end

function [dual_info, voltage_results, calcium_results, dual_results, ...
    stim_results, roi_data] = load_standard_volpy_source(source_path)
dual_info = load_required_variable_volpy( ...
    fullfile(source_path, 'dual_info.mat'), 'dual_info');
voltage_results = load_required_variable_volpy( ...
    fullfile(source_path, 'voltage_results.mat'), 'voltage_results');
calcium_results = load_required_variable_volpy( ...
    fullfile(source_path, 'calcium_results.mat'), 'calcium_results');
dual_results = load_required_variable_volpy( ...
    fullfile(source_path, 'dual_results.mat'), 'dual_results');
stim_results = load_required_variable_volpy( ...
    fullfile(source_path, 'stim_results.mat'), 'stim_results');
roi_data = load(fullfile(source_path, '1_dual_roi_results.mat'));
if ~isfield(dual_info, 'stim_context')
    error('VolPy source dual_info is missing stim_context.');
end
if ~isfield(dual_results, 'registration') || ~isfield(roi_data, 'rois') ...
        || ~isfield(roi_data.rois, 'bwmask')
    error('VolPy source is missing the standard dual ROI result.');
end
end

function value = load_required_variable_volpy(path, variable_name)
if ~isfile(path)
    error('Required VolPy source file is missing: %s', path);
end
loaded = load(path, variable_name);
if ~isfield(loaded, variable_name)
    error('Variable %s is missing from %s', variable_name, path);
end
value = loaded.(variable_name);
end

function movie_path = resolve_channel_movie_path_from_dual_info_volpy(dual_info, role)
movie_path = "";
if ~isfield(dual_info, 'camera_source') || ~isfield(dual_info, 'camera_cfg')
    return;
end
for idx = 1:min(numel(dual_info.camera_source), numel(dual_info.camera_cfg))
    if strcmpi(string(dual_info.camera_cfg(idx).role), string(role))
        movie_path = resolve_movie_file_volpy(dual_info.camera_source(idx).path);
        return;
    end
end
end

function movie_file = resolve_movie_file_volpy(source_path)
movie_file = string(source_path);
if isfile(movie_file)
    return;
end
if ~isfolder(movie_file)
    movie_file = "";
    return;
end
listing = [dir(fullfile(movie_file, '*.tif')); ...
    dir(fullfile(movie_file, '*.tiff'))];
if isempty(listing)
    movie_file = "";
    return;
end
[~, newest_idx] = max([listing.datenum]);
movie_file = string(fullfile(listing(newest_idx).folder, ...
    listing(newest_idx).name));
end

function info = run_volpy_backend_pipeline(input_movie_path, frame_rate, ...
    flip_signal, auto_run, force_rerun, use_existing, output_dir, ...
    result_mat_path, roi_mask_file)
python_exe = 'C:\Users\DELL\anaconda3\envs\caiman\python.exe';
script_path = fullfile(fileparts(mfilename('fullpath')), ...
    'python_seg', 'run_volpy_backend.py');
if ~isfile(python_exe)
    error('VolPy backend python interpreter not found: %s', python_exe);
end
if ~isfile(script_path)
    error('VolPy backend script not found: %s', script_path);
end
should_run = force_rerun || ~(isfile(result_mat_path) && use_existing);
command_output = "existing result reused";
if should_run
    if ~auto_run && ~force_rerun
        error('VolPy result is absent and auto_run is false: %s', result_mat_path);
    end
    flip_text = ternary_volpy(flip_signal, 'true', 'false');
    command = sprintf(['"%s" "%s" "%s" --output-dir "%s" ' ...
        '--frame-rate %.12g --flip-signal %s --roi-mask "%s"'], ...
        python_exe, script_path, input_movie_path, output_dir, ...
        frame_rate, flip_text, roi_mask_file);
    [status, raw_output] = system(command);
    command_output = string(raw_output);
    if status ~= 0
        error('VolPy backend command failed:\n%s', raw_output);
    end
end
if ~isfile(result_mat_path)
    error('VolPy backend did not create result file: %s', result_mat_path);
end
manifest_path = fullfile(output_dir, 'run_manifest.json');
motion_corrected_file = "";
memmap_file = "";
if isfile(manifest_path)
    manifest = jsondecode(fileread(manifest_path));
    if isfield(manifest, 'motion_corrected_file')
        motion_corrected_file = string(manifest.motion_corrected_file);
    end
    if isfield(manifest, 'memmap_file')
        memmap_file = string(manifest.memmap_file);
    end
end
info = struct( ...
    'input_movie_path', string(input_movie_path), ...
    'output_dir', string(output_dir), ...
    'result_mat', string(result_mat_path), ...
    'manifest_path', string(manifest_path), ...
    'motion_corrected_file', motion_corrected_file, ...
    'memmap_file', memmap_file, ...
    'frame_rate', frame_rate, ...
    'flip_signal', flip_signal, ...
    'roi_mask_file', string(roi_mask_file), ...
    'script_path', string(script_path), ...
    'command_output', command_output);
end

function data = load_volpy_backend_results(result_path)
s = load(result_path);
data = struct();
data.t = cell_columns_to_matrix_volpy(s.t);
data.ts = cell_columns_to_matrix_volpy(s.ts);
data.t_rec = cell_columns_to_matrix_volpy(s.t_rec);
data.t_sub = cell_columns_to_matrix_volpy(s.t_sub);
if isfield(s, 'F0')
    data.f0 = cell_columns_to_matrix_volpy(s.F0);
else
    data.f0 = ones(size(data.t));
end
if isfield(s, 'dFF')
    data.dff = cell_columns_to_matrix_volpy(s.dFF);
else
    baseline_safe = data.f0;
    baseline_safe(abs(baseline_safe) < eps) = 1;
    data.dff = data.t ./ baseline_safe;
end
data.noise = data.t - data.t_rec;
data.snr_trace = zeros(size(data.t));
for roi_idx = 1:size(data.t, 2)
    noise_std = std(data.noise(:, roi_idx), 0, 1);
    if ~isfinite(noise_std) || noise_std <= eps
        noise_std = 1;
    end
    data.snr_trace(:, roi_idx) = data.t(:, roi_idx) ./ noise_std;
end
end

function matrix = cell_columns_to_matrix_volpy(values)
if ~iscell(values)
    matrix = double(values);
elseif isempty(values)
    matrix = [];
else
    columns = cellfun(@(value) double(value(:)), values, ...
        'UniformOutput', false);
    matrix = cell2mat(columns);
end
end

function results = build_volpy_voltage_results(results, data, roi_file, backend)
stage_names = {'raw', 'bleach_removed', 'baseline', 'noise_reference', ...
    'noise', 'sensitivity', 'snr'};
for idx = 1:numel(stage_names)
    name = stage_names{idx};
    if isfield(results.trace_results, name) ...
            && ~isfield(results.trace_results, ['standard_' name])
        results.trace_results.(['standard_' name]) = results.trace_results.(name);
    end
end
results.movie_info.motion.applied = true;
results.movie_info.motion.method = 'VolPy internal MotionCorrect';
results.movie_info.motion.shift_file = char(backend.motion_corrected_file);
results.movie_info.motion.parameter_file = char(backend.result_mat);
results.movie_info.updated_at = datetime("now");
results = store_trace_stage_volpy(results, 'raw', data.t, {}, roi_file, ...
    'volpy_trace', struct('source', char(backend.result_mat)));
results = store_trace_stage_volpy(results, 'bleach_removed', data.t, ...
    {'raw'}, roi_file, 'volpy_import_compat', struct('source_stage', 'volpy_t'));
results = store_trace_stage_volpy(results, 'baseline', data.f0, ...
    {'volpy_t'}, roi_file, 'volpy_F0', struct());
results = store_trace_stage_volpy(results, 'noise_reference', data.t_rec, ...
    {'volpy_t'}, roi_file, 'volpy_reconstructed_spike_trace', struct());
results = store_trace_stage_volpy(results, 'noise', data.noise, ...
    {'volpy_t', 'noise_reference'}, roi_file, 'volpy_residual', ...
    struct('expression', 'volpy_t - noise_reference'));
results = store_trace_stage_volpy(results, 'sensitivity', data.dff, ...
    {'volpy_t', 'baseline'}, roi_file, 'volpy_dff_import', struct());
results = store_trace_stage_volpy(results, 'snr', data.snr_trace, ...
    {'volpy_t', 'noise'}, roi_file, ...
    'volpy_trace_divided_by_residual_std', struct());
results = store_trace_stage_volpy(results, 'volpy_t', data.t, ...
    {'raw'}, roi_file, 'volpy_trace', struct('source', char(backend.result_mat)));
results = store_trace_stage_volpy(results, 'volpy_ts', data.ts, ...
    {'volpy_t'}, roi_file, 'volpy_matched_filter_trace', struct());
results = store_trace_stage_volpy(results, 'volpy_t_rec', data.t_rec, ...
    {'volpy_t'}, roi_file, 'volpy_reconstructed_spike_trace', struct());
results = store_trace_stage_volpy(results, 'volpy_subthreshold', data.t_sub, ...
    {'volpy_t'}, roi_file, 'volpy_subthreshold', struct());
results = store_trace_stage_volpy(results, 'volpy_dff', data.dff, ...
    {'volpy_t', 'baseline'}, roi_file, 'volpy_dff', struct());
results = store_trace_stage_volpy(results, 'volpy_noise', data.noise, ...
    {'volpy_t', 'volpy_t_rec'}, roi_file, 'volpy_residual', ...
    struct('expression', 'volpy_t - volpy_t_rec'));
end

function results = store_trace_stage_volpy(results, stage_name, data, ...
    parent_results, roi_file, method, parameters)
results.trace_results.(stage_name) = struct( ...
    'data', data, ...
    'time', (0:size(data, 1)-1)' / double(results.movie_info.frame_rate), ...
    'frame_rate', double(results.movie_info.frame_rate), ...
    'info', struct( ...
        'result_name', stage_name, ...
        'parent_results', {parent_results}, ...
        'roi_file', roi_file, ...
        'movie_info', results.movie_info, ...
        'method', method, ...
        'parameters', parameters, ...
        'created_at', datetime("now")));
end

function update_final_volpy_paths(save_path, stage_path, backend_dir, ...
    source_path, movie_path, control, calcium_window, voltage_polarity, ...
    calcium_polarity)
dual_info_file = fullfile(save_path, 'dual_info.mat');
dual_results_file = fullfile(save_path, 'dual_results.mat');
dual_info_data = load(dual_info_file, 'dual_info');
dual_results_data = load(dual_results_file, 'dual_results');
dual_info = dual_info_data.dual_info;
dual_results = dual_results_data.dual_results;
dual_info.analysis_backend = 'volpy_voltage_reanalysis';
dual_info.volpy_source_results_path = string(source_path);
dual_info.volpy_backend.output_dir = string(backend_dir);
dual_info.volpy_backend.result_mat = string(fullfile(backend_dir, 'volpy_results.mat'));
dual_info.volpy_backend.roi_mask_file = string(fullfile(backend_dir, 'volpy_roi_mask.mat'));
dual_info.volpy_backend.voltage_movie_path = string(movie_path);
dual_info.volpy_backend.flip_signal = logical(control.flip_signal);
dual_info.volpy_backend.calcium_smoothing_window = calcium_window;
dual_info.volpy_backend.voltage_polarity = voltage_polarity;
dual_info.volpy_backend.calcium_polarity = calcium_polarity;
if isfield(dual_info.volpy_backend, 'motion_corrected_file')
    dual_info.volpy_backend.motion_corrected_file = replace( ...
        string(dual_info.volpy_backend.motion_corrected_file), ...
        string(stage_path), string(save_path));
end
if isfield(dual_results, 'volpy_voltage_reanalysis')
    info = dual_results.volpy_voltage_reanalysis.info;
    info.result_mat = string(fullfile(backend_dir, 'volpy_results.mat'));
    info.roi_mask_file = string(fullfile(backend_dir, 'volpy_roi_mask.mat'));
    if isfield(info, 'motion_corrected_file')
        info.motion_corrected_file = replace(string(info.motion_corrected_file), ...
            string(stage_path), string(save_path));
    end
    dual_results.volpy_voltage_reanalysis.info = info;
end
save(dual_info_file, 'dual_info');
save(dual_results_file, 'dual_results', '-v7.3');
end

function output = ternary_volpy(condition, true_value, false_value)
if condition
    output = true_value;
else
    output = false_value;
end
end
