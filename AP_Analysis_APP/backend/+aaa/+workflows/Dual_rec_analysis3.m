function receipt = Dual_rec_analysis3(control, runtime)
%DUAL_REC_ANALYSIS3 Analyze and average existing Dual Cycle results.
%
% This script does not process movies or call Dual_analysis3. It can read the
% exact per-Cycle result folders recorded by run_Dual_analysis3_rec in one
% batch summary, or it can scan Cycle*/Dual_analysis3 directly when only
% rec_path is provided.
%
% Control fields:
%   batch_summary_path
%       Full path to Dual_analysis3_rec_batch_summary.mat. This is the
%       optional reproducibility input because it fixes the exact Cycle
%       result versions.
%   rec_path
%       Record folder containing Cycle* folders. When batch_summary_path is
%       empty, the script scans each Cycle*/Dual_analysis3 folder directly
%       and uses the newest complete Cycle result.
%   output_dir_name
%       Output folder under the Rec path. Default: Dual_analysis3_rec_average.
%   voltage_polarity, calcium_polarity, calcium_smoothing_window
%       Optional explicit overrides. Empty first loads values saved in the
%       batch summary, then tries to infer them from Cycle dual_info.mat.
%       calcium_smoothing_window has no safe generic fallback because older
%       datasets may use different windows.
%
% Example:
%   dual_rec_analysis3_control = struct( ...
%       'batch_summary_path', ...
%           'D:\Data\Rec1_2026-01-01_12-00-00\Dual_analysis3_rec_batch_summary.mat');
%   aaa.workflows.Dual_rec_analysis3(control);
%
% Section workflow in the MATLAB Editor:
%   Run %% Control once, then run %% 1 once to create/update rec_avg.mat.
%   Afterwards, run any one of sections %% 2 through %% 7 independently.
%   Plotting sections reload rec_avg.mat and existing Cycle results; they do
%   not recalculate the record average.
if nargin < 1 || isempty(control)
    control = struct();
end
if nargin < 2 || isempty(runtime)
    runtime = struct();
end
dual_rec_analysis3_control = control;
emit_record_average_event(runtime, "running", ...
    "Starting Dual Record-average workflow.", "");
try
%% Control
if ~exist('dual_rec_analysis3_control', 'var') || isempty(dual_rec_analysis3_control)
    dual_rec_analysis3_control = struct();
end
user_supplied_batch_summary_path = isfield(dual_rec_analysis3_control, 'batch_summary_path') ...
    && ~isempty(dual_rec_analysis3_control.batch_summary_path);
user_supplied_rec_path = isfield(dual_rec_analysis3_control, 'rec_path') ...
    && ~isempty(dual_rec_analysis3_control.rec_path);
control_defaults = struct( ...
    'batch_summary_path', "", ...
    'rec_path', "", ...
    'output_dir_name', "Dual_analysis3_rec_average", ...
    'voltage_polarity', [], ...
    'calcium_polarity', [], ...
    'calcium_smoothing_window', []);
if ~isstruct(dual_rec_analysis3_control) || ~isscalar(dual_rec_analysis3_control)
    error('dual_rec_analysis3_control must be one scalar struct.');
end
unknown_fields = setdiff(fieldnames(dual_rec_analysis3_control), fieldnames(control_defaults));
if ~isempty(unknown_fields)
    error('Unknown dual_rec_analysis3_control field(s): %s', strjoin(unknown_fields, ', '));
end
control_fields = fieldnames(control_defaults);
for control_idx = 1:numel(control_fields)
    field_name = control_fields{control_idx};
    if ~isfield(dual_rec_analysis3_control, field_name) || isempty(dual_rec_analysis3_control.(field_name))
        dual_rec_analysis3_control.(field_name) = control_defaults.(field_name);
    end
end
if user_supplied_batch_summary_path && ~user_supplied_rec_path
    dual_rec_analysis3_control.rec_path = "";
end

%% 1. Build or refresh record-average data (no figures; writes rec_avg.mat)
build_input = resolve_rec_average_build_input_rec(dual_rec_analysis3_control);
record_average_output = build_record_average_summary( ...
    build_input.rec_path, build_input.batch_results, build_input.reference_roi_file, ...
    build_input.voltage_polarity, build_input.calcium_polarity, ...
    build_input.calcium_smoothing_window, build_input.average_output_dir_name);
fprintf('Record-average data bundle:\n  %s\n', record_average_output.results_file);

%% 2. Plot ROI-cycle heatmaps and peak rasters (4_roi_cyc_hm_v_sens, 4_roi_cyc_hm_v_snr, 5_roi_cyc_ca_sens_peak)
plot_context = load_record_average_plot_context_rec(dual_rec_analysis3_control);
plot_record_roi_cycle_module_rec(plot_context);

%% 3. Plot flash FFT and wavelet summaries (6_fft, 7_wav; flash records only)
plot_context = load_record_average_plot_context_rec(dual_rec_analysis3_control);
plot_record_flash_time_frequency_module_rec(plot_context);

%% 4. Plot record-average grating tuning summaries (8_grating, 8_grating_roi, 8_grating_v_sens, 8_grating_v_sens_roi)
plot_context = load_record_average_plot_context_rec(dual_rec_analysis3_control);
plot_record_grating_tuning_module_rec(plot_context);

%% 5. Plot random-grating trial stacks (9_rg_stack; includes roi and v_ca stack+tuning figures)
plot_context = load_record_average_plot_context_rec(dual_rec_analysis3_control);
plot_record_random_grating_stack_module_rec(plot_context);

%% 6. Plot random-grating sensitivity heatmaps (10_rg_hm)
plot_context = load_record_average_plot_context_rec(dual_rec_analysis3_control);
plot_record_random_grating_heatmap_module_rec(plot_context);

%% 7. Plot random-grating cycle traces, DSI and per-cycle tuning (11_rg_v_ca)
plot_context = load_record_average_plot_context_rec(dual_rec_analysis3_control);
plot_record_random_grating_cycle_module_rec(plot_context);

output_dir = string(fileparts(record_average_output.results_file));
receipt = struct( ...
    'status', "completed", ...
    'mode', "dual", ...
    'scope', "record_average", ...
    'input_path', string(build_input.rec_path), ...
    'output_path', output_dir, ...
    'manifest_file', "", ...
    'completed_sections', ["record_average";"record_visualization"], ...
    'failed_section', "", ...
    'error', struct(), ...
    'results_file', string(record_average_output.results_file), ...
    'output_files', list_record_average_files(output_dir));
emit_record_average_event(runtime, "completed", ...
    "Dual Record-average workflow completed.", output_dir);
catch ME
    emit_record_average_event(runtime, "failed", string(ME.message), "");
    rethrow(ME);
end
end

%% Local Functions
function build_input = resolve_rec_average_build_input_rec(control)
batch_summary_path = char(string(control.batch_summary_path));
rec_path = char(string(control.rec_path));
saved_summary = struct();
reference_roi_file = "";
if ~isempty(batch_summary_path)
    if ~isfile(batch_summary_path)
        error('Batch summary does not exist: %s', batch_summary_path);
    end
    saved_summary = load(batch_summary_path);
    if ~isfield(saved_summary, 'batch_results')
        error('Batch summary does not contain batch_results: %s', batch_summary_path);
    end
    if isempty(rec_path)
        rec_path = fileparts(batch_summary_path);
    end
    batch_results = saved_summary.batch_results;
    if isfield(saved_summary, 'reference_roi_file')
        reference_roi_file = saved_summary.reference_roi_file;
    end
    input_source = "batch_summary";
else
    if isempty(rec_path) || ~isfolder(rec_path)
        error('Set dual_rec_analysis3_control.rec_path to an existing record folder.');
    end
    [batch_results, reference_roi_file] = scan_record_cycle_results(rec_path);
    input_source = "direct_cycle_scan";
end
build_input = struct( ...
    'rec_path', string(rec_path), ...
    'batch_results', batch_results, ...
    'reference_roi_file', string(reference_roi_file), ...
    'voltage_polarity', resolve_rec_analysis_parameter( ...
        control.voltage_polarity, saved_summary, batch_results, 'voltage_polarity', -1, true), ...
    'calcium_polarity', resolve_rec_analysis_parameter( ...
        control.calcium_polarity, saved_summary, batch_results, 'calcium_polarity', 1, true), ...
    'calcium_smoothing_window', resolve_rec_analysis_parameter( ...
        control.calcium_smoothing_window, saved_summary, batch_results, 'calcium_smoothing_window', 40, false), ...
    'average_output_dir_name', string(control.output_dir_name));
fprintf('Build source: %s\n', input_source);
end

function context = load_record_average_plot_context_rec(control)
rec_path = char(string(control.rec_path));
if isempty(rec_path) && strlength(string(control.batch_summary_path)) > 0
    rec_path = fileparts(char(string(control.batch_summary_path)));
end
if isempty(rec_path)
    error('Set dual_rec_analysis3_control.rec_path before running a plotting section.');
end
output_dir = fullfile(rec_path, char(string(control.output_dir_name)));
results_file = fullfile(output_dir, 'rec_avg.mat');
if ~isfile(results_file)
    error('Missing record-average data bundle: %s\nRun section 1 first.', results_file);
end
saved = load(results_file, 'record_average');
if ~isfield(saved, 'record_average') || ~isfield(saved.record_average, 'info')
    error('Invalid record-average data bundle: %s', results_file);
end
record_average = saved.record_average;
if ~isfield(record_average.info, 'result_dirs') || ~isfield(record_average.info, 'cycle_names')
    error('Record-average data bundle is missing cycle source information: %s', results_file);
end
result_dirs = string(record_average.info.result_dirs(:));
cycle_names = string(record_average.info.cycle_names(:));
if numel(result_dirs) ~= numel(cycle_names) || isempty(result_dirs)
    error('Record-average data bundle has invalid cycle source information: %s', results_file);
end
cycle_entry_list = cell(numel(result_dirs), 1);
for idx = 1:numel(result_dirs)
    cycle_entry_list{idx} = load_cycle_average_entry(result_dirs(idx), cycle_names(idx));
end
context = struct( ...
    'record_average', record_average, ...
    'cycle_entries', vertcat(cycle_entry_list{:}), ...
    'output_dir', string(output_dir), ...
    'voltage_polarity', double(record_average.info.voltage_polarity), ...
    'calcium_polarity', double(record_average.info.calcium_polarity), ...
    'calcium_smoothing_window', double(record_average.info.calcium_smoothing_window), ...
    'stim_type', string(record_average.info.stim_type));
end

function plot_record_roi_cycle_module_rec(context)
record_average = context.record_average;
output_dir = char(context.output_dir);
plot_record_roi_cycle_heatmap_voltage_trace( ...
    record_average.voltage.sensitivity, record_average.calcium.sensitivity, ...
    record_average.stim_windows, fullfile(output_dir, '4_roi_cyc_hm_v_sens'), ...
    fullfile(output_dir, '4_roi_cyc_hm_v_sens.mat'), ...
    context.voltage_polarity, context.calcium_polarity, ...
    context.calcium_smoothing_window, "sensitivity");
plot_record_roi_cycle_heatmap_voltage_trace( ...
    record_average.voltage.snr, record_average.calcium.sensitivity, ...
    record_average.stim_windows, fullfile(output_dir, '4_roi_cyc_hm_v_snr'), ...
    fullfile(output_dir, '4_roi_cyc_hm_v_snr.mat'), ...
    context.voltage_polarity, context.calcium_polarity, ...
    context.calcium_smoothing_window, "snr");
plot_record_roi_cycle_calcium_peak_raster( ...
    record_average.voltage.sensitivity, record_average.calcium.sensitivity, ...
    context.cycle_entries, record_average.stim_windows, ...
    fullfile(output_dir, '5_roi_cyc_ca_sens_peak'), ...
    fullfile(output_dir, '5_roi_cyc_ca_sens_peak.mat'), ...
    context.calcium_polarity, "sensitivity");
end

function plot_record_flash_time_frequency_module_rec(context)
if context.stim_type ~= "flash"
    fprintf('Skipping flash FFT/wavelet section: record stimulus type is %s.\n', context.stim_type);
    return;
end
record_average = context.record_average;
tf_input = struct( ...
    'voltage', struct( ...
        'trace', context.voltage_polarity * record_average.voltage.sensitivity.average, ...
        'time', record_average.voltage.sensitivity.time, ...
        'frame_rate', record_average.voltage.sensitivity.frame_rate, ...
        'stage_name', "sensitivity_average"), ...
    'calcium', struct( ...
        'trace', context.calcium_polarity * record_average.calcium.sensitivity.average, ...
        'time', record_average.calcium.sensitivity.time, ...
        'frame_rate', record_average.calcium.sensitivity.frame_rate, ...
        'stage_name', "sensitivity_average"));
build_record_average_time_frequency(tf_input, record_average.stim_windows, char(context.output_dir));
end

function plot_record_grating_tuning_module_rec(context)
if context.stim_type ~= "grating"
    fprintf('Skipping grating-tuning section: record stimulus type is %s.\n', context.stim_type);
    return;
end
build_record_average_grating_tuning(context.cycle_entries, char(context.output_dir));
end

function [angles, nrois] = resolve_random_grating_plot_inputs_rec(context)
is_random = arrayfun(@(entry) strcmpi( ...
    string(entry.stim_analysis_kind), "random_grating_tuning"), context.cycle_entries);
if ~any(is_random)
    error('Random-grating plotting sections require random_grating_tuning cycles.');
end
if ~all(is_random)
    error('Random and non-random grating cycles cannot share one random-grating plotting section.');
end
nrois = size(context.cycle_entries(1).voltage.sensitivity.data, 2);
angle_cells = arrayfun(@(entry) mod(double(entry.stim_windows.orientations(:)), 360), ...
    context.cycle_entries, 'UniformOutput', false);
angles = unique(vertcat(angle_cells{:}), 'sorted')';
end

function plot_record_random_grating_stack_module_rec(context)
[angles, nrois] = resolve_random_grating_plot_inputs_rec(context);
stack_dir = fullfile(char(context.output_dir), '9_rg_stack');
if ~isfolder(stack_dir)
    mkdir(stack_dir);
end
plot_random_grating_trial_stack_grid_rec( ...
    context.cycle_entries, 'voltage', angles, nrois, context.voltage_polarity, ...
    fullfile(stack_dir, 'v_stack.fig'), fullfile(stack_dir, 'v_stack.png'));
plot_random_grating_trial_stack_grid_rec( ...
    context.cycle_entries, 'calcium', angles, nrois, context.calcium_polarity, ...
    fullfile(stack_dir, 'ca_stack.fig'), fullfile(stack_dir, 'ca_stack.png'));
by_roi_root = fullfile(stack_dir, 'roi');
plot_random_grating_trial_stacks_by_roi_rec( ...
    context.cycle_entries, 'voltage', angles, nrois, context.voltage_polarity, by_roi_root);
plot_random_grating_trial_stacks_by_roi_rec( ...
    context.cycle_entries, 'calcium', angles, nrois, context.calcium_polarity, by_roi_root);
plot_random_grating_roi_voltage_calcium_stack_tuning_rec( ...
    context.cycle_entries, angles, nrois, context.voltage_polarity, ...
    context.calcium_polarity, stack_dir);
end

function plot_record_random_grating_heatmap_module_rec(context)
[angles, nrois] = resolve_random_grating_plot_inputs_rec(context);
heatmap_dir = fullfile(char(context.output_dir), '10_rg_hm');
if ~isfolder(heatmap_dir)
    mkdir(heatmap_dir);
end
for channel_name = ["voltage", "calcium"]
    if channel_name == "voltage"
        polarity = context.voltage_polarity;
        prefix = 'v';
    else
        polarity = context.calcium_polarity;
        prefix = 'ca';
    end
    color_limit = resolve_random_grating_heatmap_clim_rec( ...
        context.cycle_entries, channel_name, angles, nrois, polarity);
    plot_random_grating_heatmap_grid_rec( ...
        context.cycle_entries, channel_name, angles, nrois, polarity, color_limit, ...
        fullfile(heatmap_dir, [prefix, '_hm.fig']), ...
        fullfile(heatmap_dir, [prefix, '_hm.png']));
    plot_random_grating_heatmaps_by_roi_rec( ...
        context.cycle_entries, channel_name, angles, nrois, polarity, color_limit, ...
        fullfile(heatmap_dir, 'roi'));
end
end

function plot_record_random_grating_cycle_module_rec(context)
[angles, nrois] = resolve_random_grating_plot_inputs_rec(context);
subplot_root = fullfile(char(context.output_dir), '11_rg_v_ca');
plot_random_grating_roi_cycle_voltage_calcium_subplots_rec( ...
    context.cycle_entries, nrois, context.voltage_polarity, context.calcium_polarity, subplot_root);
build_random_grating_cycle_voltage_dsi_comparison_rec( ...
    context.cycle_entries, angles, nrois, subplot_root);
plot_random_grating_cycle_tuning_curves_rec( ...
    context.cycle_entries, angles, nrois, subplot_root);
end

function value = resolve_rec_analysis_parameter(control_value, saved_summary, batch_results, field_name, fallback_value, allow_fallback)
if ~isempty(control_value)
    value = control_value;
elseif isfield(saved_summary, field_name) && ~isempty(saved_summary.(field_name))
    value = saved_summary.(field_name);
else
    value = infer_rec_analysis_parameter_from_cycle_results(batch_results, field_name);
    if isempty(value)
        if allow_fallback
            value = fallback_value;
            fprintf(['  %s not found in batch summary or cycle dual_info; ' ...
                'using legacy fallback %g.\n'], field_name, fallback_value);
        else
            error(['Cannot resolve %s from control, batch summary, or cycle dual_info.mat. ' ...
                'Set dual_rec_analysis3_control.%s explicitly.'], field_name, field_name);
        end
    end
end
if ~isnumeric(value) || ~isscalar(value) || ~isfinite(value)
    error('%s must be one finite numeric scalar.', field_name);
end
value = double(value);
end

function [batch_results, reference_roi_file] = scan_record_cycle_results(rec_path)
cycle_dirs = dir(fullfile(rec_path, 'Cycle*'));
cycle_dirs = cycle_dirs([cycle_dirs.isdir]);
cycle_dirs = sort_cycle_dirs_rec(cycle_dirs);
if isempty(cycle_dirs)
    error('No Cycle* folders found under record path: %s', rec_path);
end

batch_results = repmat(struct( ...
    'cycle_name', "", ...
    'cycle_path', "", ...
    'status', "", ...
    'save_path', "", ...
    'roi_file', "", ...
    'message', ""), 0, 1);

for idx = 1:numel(cycle_dirs)
    cycle_name = string(cycle_dirs(idx).name);
    cycle_path = fullfile(cycle_dirs(idx).folder, cycle_dirs(idx).name);
    [result_dir, selection_reason] = find_latest_usable_dual_result_dir(cycle_path);
    if strlength(result_dir) > 0
        batch_results(end+1, 1) = struct( ... %#ok<AGROW>
            'cycle_name', cycle_name, ...
            'cycle_path', string(cycle_path), ...
            'status', "existing_result_detected", ...
            'save_path', result_dir, ...
            'roi_file', string(fullfile(result_dir, '1_dual_roi_results.mat')), ...
            'message', selection_reason);
        fprintf('%s | using latest Rec-analysis-usable result:\n  %s\n  %s\n', ...
            cycle_name, result_dir, selection_reason);
    else
        batch_results(end+1, 1) = struct( ... %#ok<AGROW>
            'cycle_name', cycle_name, ...
            'cycle_path', string(cycle_path), ...
            'status', "missing_result", ...
            'save_path', "", ...
            'roi_file', "", ...
            'message', selection_reason);
        fprintf('%s | no Rec-analysis-usable Dual_analysis3 result folder found. %s\n', ...
            cycle_name, selection_reason);
    end
end

reference_roi_file = "";
for idx = 1:numel(batch_results)
    roi_file = string(batch_results(idx).roi_file);
    if strlength(roi_file) > 0 && isfile(roi_file)
        reference_roi_file = roi_file;
        break;
    end
end
end

function [result_dir, selection_reason] = find_latest_usable_dual_result_dir(cycle_path)
result_dir = "";
selection_reason = "No voltage_results.mat files found under Cycle*/Dual_analysis3.";
listing = dir(fullfile(cycle_path, 'Dual_analysis3', '**', 'voltage_results.mat'));
if isempty(listing)
    return;
end
[~, order] = sort([listing.datenum], 'descend');
listing = listing(order);
required_files = {'dual_info.mat', 'calcium_results.mat', 'stim_results.mat'};
rejected = strings(0, 1);
for idx = 1:numel(listing)
    candidate_dir = string(listing(idx).folder);
    missing_files = required_files(~cellfun(@(name) isfile(fullfile(candidate_dir, name)), required_files));
    if ~isempty(missing_files)
        rejected(end+1, 1) = sprintf('%s: missing %s', char(candidate_dir), strjoin(missing_files, ', ')); %#ok<AGROW>
        continue;
    end
    if ~result_dir_has_rec_average_stages(candidate_dir)
        rejected(end+1, 1) = sprintf('%s: missing required trace stages', char(candidate_dir)); %#ok<AGROW>
        continue;
    end
    [stim_is_usable, stim_reason] = result_dir_has_rec_usable_stim(candidate_dir);
    if ~stim_is_usable
        rejected(end+1, 1) = sprintf('%s: %s', char(candidate_dir), stim_reason); %#ok<AGROW>
        continue;
    end
    result_dir = candidate_dir;
    selection_reason = "Selected because trace stages and Rec-analysis stimulus results are usable.";
    return;
end
if ~isempty(rejected)
    selection_reason = "Rejected candidates: " + strjoin(rejected(1:min(3, numel(rejected))), " | ");
    if numel(rejected) > 3
        selection_reason = selection_reason + sprintf(' | ... %d more', numel(rejected) - 3);
    end
end
end

function is_valid = result_dir_has_rec_average_stages(result_dir)
is_valid = false;
try
    tmp_voltage = load(fullfile(result_dir, 'voltage_results.mat'), 'voltage_results');
    tmp_calcium = load(fullfile(result_dir, 'calcium_results.mat'), 'calcium_results');
catch
    return;
end
if ~isfield(tmp_voltage, 'voltage_results') || ~isfield(tmp_calcium, 'calcium_results')
    return;
end
voltage_results = tmp_voltage.voltage_results;
calcium_results = tmp_calcium.calcium_results;
is_valid = has_stage_for_rec_average(voltage_results, {'raw'}) ...
    && has_stage_for_rec_average(voltage_results, {'sensitivity'}) ...
    && has_stage_for_rec_average(voltage_results, {'snr'}) ...
    && has_stage_for_rec_average(calcium_results, {'raw_smoothed', 'raw'}) ...
    && has_stage_for_rec_average(calcium_results, {'sensitivity_smoothed', 'sensitivity'}) ...
    && has_stage_for_rec_average(calcium_results, {'snr_smoothed', 'snr'});
end

function [is_usable, reason] = result_dir_has_rec_usable_stim(result_dir)
is_usable = false;
reason = "stim_results.mat does not contain usable flash or grating results";
try
    saved = load(fullfile(result_dir, 'stim_results.mat'), 'stim_results');
catch ME
    reason = "stim_results.mat could not be loaded: " + string(ME.message);
    return;
end
if ~isfield(saved, 'stim_results') || ~isstruct(saved.stim_results)
    reason = "stim_results variable is missing or not a struct";
    return;
end
stim_results = saved.stim_results;
if isfield(stim_results, 'analysis_kind') ...
        && is_grating_tuning_kind_rec(stim_results.analysis_kind)
    if isfield(stim_results, 'tuning') && isstruct(stim_results.tuning) ...
            && isfield(stim_results.tuning, 'voltage') ...
            && isfield(stim_results.tuning, 'calcium') ...
            && has_complete_grating_curves_rec(stim_results.tuning.voltage) ...
            && has_complete_grating_curves_rec(stim_results.tuning.calcium)
        is_usable = true;
        reason = "usable " + string(stim_results.analysis_kind) + " found in stim_results.mat";
    else
        reason = string(stim_results.analysis_kind) + " is missing complete voltage/calcium curves";
    end
    return;
end
if isfield(stim_results, 'windows') && isstruct(stim_results.windows) ...
        && isfield(stim_results.windows, 'supported') ...
        && logical(stim_results.windows.supported) ...
        && isfield(stim_results.windows, 'stim_type') ...
        && strcmpi(string(stim_results.windows.stim_type), "visualstim_flash")
    is_usable = true;
    reason = "usable flash windows found in stim_results.mat";
    return;
end
if isfield(stim_results, 'windows') && isstruct(stim_results.windows) ...
        && isfield(stim_results.windows, 'is_grating') ...
        && logical(stim_results.windows.is_grating)
    reason = "grating windows exist but grating/random_grating tuning is missing or incomplete";
end
end

function tf = is_grating_tuning_kind_rec(analysis_kind)
tf = any(strcmpi(string(analysis_kind), ["grating_tuning", "random_grating_tuning"]));
end

function tf = has_stage_for_rec_average(results_struct, preferred_stages)
tf = false;
if ~isstruct(results_struct) || ~isfield(results_struct, 'trace_results') ...
        || ~isstruct(results_struct.trace_results)
    return;
end
for idx = 1:numel(preferred_stages)
    stage_name = preferred_stages{idx};
    if isfield(results_struct.trace_results, stage_name) ...
            && isfield(results_struct.trace_results.(stage_name), 'data') ...
            && ~isempty(results_struct.trace_results.(stage_name).data)
        tf = true;
        return;
    end
end
end

function value = infer_rec_analysis_parameter_from_cycle_results(batch_results, field_name)
value = [];
[result_dirs, ~] = resolve_successful_result_dirs(batch_results);
values = [];
for idx = 1:numel(result_dirs)
    dual_info_file = fullfile(result_dirs(idx), 'dual_info.mat');
    if ~isfile(dual_info_file)
        continue;
    end
    try
        saved = load(dual_info_file, 'dual_info');
    catch
        continue;
    end
    if ~isfield(saved, 'dual_info') || ~isstruct(saved.dual_info)
        continue;
    end
    candidate = extract_parameter_from_dual_info(saved.dual_info, field_name);
    if isnumeric(candidate) && isscalar(candidate) && isfinite(candidate)
        values(end+1, 1) = double(candidate); %#ok<AGROW>
    end
end
if isempty(values)
    return;
end
if any(abs(values - values(1)) > 1e-9)
    error('Inconsistent %s values across Cycle dual_info.mat files: %s', ...
        field_name, strjoin(string(values(:))', ', '));
end
value = values(1);
fprintf('  %s inferred from cycle dual_info.mat: %g\n', field_name, value);
end

function value = extract_parameter_from_dual_info(dual_info, field_name)
value = [];
switch char(field_name)
    case 'voltage_polarity'
        if isfield(dual_info, 'polarity') && isstruct(dual_info.polarity) ...
                && isfield(dual_info.polarity, 'voltage')
            value = dual_info.polarity.voltage;
        end
    case 'calcium_polarity'
        if isfield(dual_info, 'polarity') && isstruct(dual_info.polarity) ...
                && isfield(dual_info.polarity, 'calcium')
            value = dual_info.polarity.calcium;
        end
    case 'calcium_smoothing_window'
        if isfield(dual_info, 'calcium_smoothing_window')
            value = dual_info.calcium_smoothing_window;
        end
end
end

function cycle_dirs = sort_cycle_dirs_rec(cycle_dirs)
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
    'calcium_smoothing_window', calcium_smoothing_window, ...
    'voltage_polarity', voltage_polarity, ...
    'calcium_polarity', calcium_polarity);
stim_groups = strings(numel(cycle_entries), 1);
for idx = 1:numel(cycle_entries)
    stim_groups(idx) = classify_cycle_stim_group_rec(cycle_entries(idx));
end
record_average.info.stim_group_by_cycle = stim_groups;
if any(stim_groups == "unsupported")
    unsupported_cycles = string({cycle_entries(stim_groups == "unsupported").cycle_name})';
    error('Cannot identify flash/grating stimulus type for cycles: %s', ...
        strjoin(unsupported_cycles, ', '));
end
record_stim_types = unique(stim_groups);
if numel(record_stim_types) ~= 1
    error('One Rec must contain one stimulus type. Found: %s', ...
        strjoin(record_stim_types, ', '));
end
record_stim_type = record_stim_types(1);
record_average.info.stim_type = record_stim_type;
fprintf('Building %s record average from %d cycles.\n', record_stim_type, numel(cycle_entries));
stim_average = build_record_average_stim_group(cycle_entries, record_stim_type);
record_average.voltage = stim_average.voltage;
record_average.calcium = stim_average.calcium;
record_average.stim_windows = stim_average.stim_windows;
record_average.visualizations = stim_average.visualizations;
record_average.time_frequency = stim_average.time_frequency;
record_average.grating_tuning = stim_average.grating_tuning;

results_file = fullfile(output_dir, 'rec_avg.mat');
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
    if status_value == "completed" ...
            || status_value == "skipped_existing_result" ...
            || status_value == "skipped_existing_analysis_only_result" ...
            || status_value == "existing_result_detected"
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
entry.voltage.peaks = extract_voltage_accepted_peaks(result_dir, voltage_results);
entry.calcium.raw = extract_cycle_stage(calcium_results, {'raw_smoothed', 'raw'});
entry.calcium.sensitivity = extract_cycle_stage(calcium_results, {'sensitivity_smoothed', 'sensitivity'});
entry.calcium.snr = extract_cycle_stage(calcium_results, {'snr_smoothed', 'snr'});

entry.stim_windows = struct();
entry.grating_tuning = struct();
entry.stim_response = struct();
entry.stim_analysis_kind = "";
if isfile(stim_path)
    tmp_stim = load(stim_path, 'stim_results');
    if isfield(tmp_stim, 'stim_results') && isstruct(tmp_stim.stim_results) ...
            && isfield(tmp_stim.stim_results, 'analysis_kind')
        entry.stim_analysis_kind = string(tmp_stim.stim_results.analysis_kind);
    end
    if isfield(tmp_stim, 'stim_results') && isstruct(tmp_stim.stim_results) ...
            && isfield(tmp_stim.stim_results, 'windows')
        entry.stim_windows = tmp_stim.stim_results.windows;
    end
    if isfield(tmp_stim, 'stim_results') && isstruct(tmp_stim.stim_results) ...
            && isfield(tmp_stim.stim_results, 'response')
        entry.stim_response = tmp_stim.stim_results.response;
    end
    if isfield(tmp_stim, 'stim_results') && isstruct(tmp_stim.stim_results) ...
            && isfield(tmp_stim.stim_results, 'analysis_kind') ...
            && is_grating_tuning_kind_rec(tmp_stim.stim_results.analysis_kind) ...
            && isfield(tmp_stim.stim_results, 'tuning')
        entry.grating_tuning = tmp_stim.stim_results.tuning;
    end
end
entry.grating_tuning = attach_voltage_sensitivity_tuning_from_response_rec( ...
    entry.grating_tuning, entry.stim_response, entry.stim_windows);
entry.alignment = struct( ...
    'voltage_flash_onset_time', resolve_flash_onset_time(entry.stim_windows, 'voltage'), ...
    'calcium_flash_onset_time', resolve_flash_onset_time(entry.stim_windows, 'calcium'));
end

function group_name = classify_cycle_stim_group_rec(entry)
group_name = "unsupported";
if isstruct(entry.grating_tuning) && ~isempty(fieldnames(entry.grating_tuning))
    group_name = "grating";
    return;
end
if isstruct(entry.stim_windows)
    if isfield(entry.stim_windows, 'is_grating') && logical(entry.stim_windows.is_grating)
        group_name = "grating";
        return;
    end
    if isfield(entry.stim_windows, 'stim_type') ...
            && strcmpi(string(entry.stim_windows.stim_type), "visualstim_flash")
        group_name = "flash";
    end
end
end

function tuning = attach_voltage_sensitivity_tuning_from_response_rec(tuning, stim_response, stim_windows)
if ~isstruct(tuning) || isempty(fieldnames(tuning))
    return;
end
if isfield(tuning, 'voltage_sensitivity') && has_complete_grating_curves_rec(tuning.voltage_sensitivity)
    return;
end
if ~isstruct(stim_response) || ~isfield(stim_response, 'sensitivity') ...
        || ~isstruct(stim_response.sensitivity) ...
        || ~isfield(stim_response.sensitivity, 'voltage') ...
        || ~isstruct(stim_response.sensitivity.voltage)
    return;
end
if ~isstruct(stim_windows) || ~isfield(stim_windows, 'orientations')
    return;
end
voltage_metrics = stim_response.sensitivity.voltage;
if ~isfield(voltage_metrics, 'stim_mean') || ~isfield(voltage_metrics, 'baseline_mean')
    return;
end
try
    tuning.voltage_sensitivity = build_grating_tuning_from_trial_response_rec( ...
        double(voltage_metrics.stim_mean), ...
        double(voltage_metrics.baseline_mean), ...
        stim_windows.orientations);
    tuning.trial_response.voltage_mean_sensitivity = struct( ...
        'stim', double(voltage_metrics.stim_mean), ...
        'baseline', double(voltage_metrics.baseline_mean));
    if isfield(tuning, 'curve_definition') && isstruct(tuning.curve_definition)
        tuning.curve_definition.voltage_sensitivity = ...
            "mean sensitivity per paired window; sensitivity baseline is zero";
    end
catch ME
    warning('Dual_rec_analysis3:VoltageSensitivityTuningSkipped', ...
        'Could not build voltage_sensitivity grating tuning from saved response: %s', ME.message);
end
end

function group_average = build_record_average_stim_group(cycle_entries, group_name)
if group_name == "flash"
    alignment_mode = "flash_onset";
else
    alignment_mode = "recording_start";
end

group_average = struct();
group_average.info = struct( ...
    'stim_group', group_name, ...
    'cycle_names', string({cycle_entries.cycle_name})', ...
    'result_dirs', string({cycle_entries.result_dir})', ...
    'cycle_count', numel(cycle_entries), ...
    'alignment_mode', alignment_mode, ...
    'created_at', datetime("now"));
group_average.voltage.raw = average_cycle_stage(cycle_entries, 'voltage', 'raw', alignment_mode);
group_average.voltage.sensitivity = average_cycle_stage(cycle_entries, 'voltage', 'sensitivity', alignment_mode);
group_average.voltage.snr = average_cycle_stage(cycle_entries, 'voltage', 'snr', alignment_mode);
group_average.calcium.raw = average_cycle_stage(cycle_entries, 'calcium', 'raw', alignment_mode);
group_average.calcium.sensitivity = average_cycle_stage(cycle_entries, 'calcium', 'sensitivity', alignment_mode);
group_average.calcium.snr = average_cycle_stage(cycle_entries, 'calcium', 'snr', alignment_mode);
group_average.stim_windows = resolve_record_average_stim_windows( ...
    cycle_entries, group_average.voltage.raw.time, group_average.calcium.raw.time, alignment_mode);

group_average.visualizations = struct('status', "not_generated");
group_average.time_frequency = struct('status', "not_generated");
group_average.grating_tuning = struct('status', "not_generated");

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

function peaks = extract_voltage_accepted_peaks(result_dir, voltage_results)
peaks = struct( ...
    'status', "missing", ...
    'source_file', "", ...
    'index', {{}}, ...
    'amplitude', {{}}, ...
    'polarity', {{}});
peak_results = struct();
if isstruct(voltage_results) && isfield(voltage_results, 'peak_results') ...
        && isstruct(voltage_results.peak_results)
    peak_results = voltage_results.peak_results;
    peaks.source_file = string(fullfile(result_dir, 'voltage_results.mat'));
elseif isfile(fullfile(result_dir, 'voltage_peak_results.mat'))
    try
        saved = load(fullfile(result_dir, 'voltage_peak_results.mat'), 'peak_results');
        if isfield(saved, 'peak_results') && isstruct(saved.peak_results)
            peak_results = saved.peak_results;
            peaks.source_file = string(fullfile(result_dir, 'voltage_peak_results.mat'));
        end
    catch
        peak_results = struct();
    end
end

if ~isstruct(peak_results) || ~isfield(peak_results, 'accepted_for_events') ...
        || ~isstruct(peak_results.accepted_for_events) ...
        || ~isfield(peak_results.accepted_for_events, 'data') ...
        || ~isstruct(peak_results.accepted_for_events.data)
    return;
end
accepted = peak_results.accepted_for_events.data;
if ~isfield(accepted, 'index') || ~iscell(accepted.index)
    return;
end
peaks.index = accepted.index;
if isfield(accepted, 'amplitude') && iscell(accepted.amplitude)
    peaks.amplitude = accepted.amplitude;
else
    peaks.amplitude = cell(size(accepted.index));
end
if isfield(accepted, 'polarity') && iscell(accepted.polarity)
    peaks.polarity = accepted.polarity;
else
    peaks.polarity = cell(size(accepted.index));
end
peaks.status = "loaded";
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

cycle_roi_mean = squeeze(mean(per_cycle, 2, 'omitnan'));
cycle_roi_mean = reshape(cycle_roi_mean, numel(common_time), numel(cycle_entries));
averaged = struct( ...
    'stage_name', string(stage_labels(1)), ...
    'stage_name_by_cycle', stage_labels, ...
    'frame_rate', frame_rate, ...
    'time', common_time, ...
    'per_cycle', per_cycle, ...
    'cycle_roi_mean', cycle_roi_mean, ...
    'average', mean(per_cycle, 3, 'omitnan'), ...
    'cycle_names', string({cycle_entries.cycle_name})', ...
    'alignment_mode', string(alignment_mode), ...
    'alignment_time_by_cycle', alignment_time_by_cycle, ...
    'ncycles', numel(cycle_entries), ...
    'nrois', nrois);
end

function stim_windows = resolve_record_average_stim_windows(cycle_entries, voltage_time, calcium_time, alignment_mode)
stim_windows = struct();
for idx = 1:numel(cycle_entries)
    candidate = cycle_entries(idx).stim_windows;
    if isstruct(candidate) && isfield(candidate, 'supported') && candidate.supported
        stim_windows = candidate;
        if string(alignment_mode) == "recording_start"
            voltage_anchor = cycle_entries(idx).voltage.raw.time(1);
            calcium_anchor = cycle_entries(idx).calcium.raw.time(1);
        else
            voltage_anchor = cycle_entries(idx).alignment.voltage_flash_onset_time;
            calcium_anchor = cycle_entries(idx).alignment.calcium_flash_onset_time;
        end
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
    case "recording_start"
        alignment_time = cycle_entry.(channel_name).raw.time(1);
    otherwise
        error('Unsupported record-average alignment mode: %s', alignment_mode);
end
end

function summary = build_record_average_grating_tuning(cycle_entries, output_dir)
summary = struct( ...
    'status', "not_available", ...
    'reason', "No compatible grating tuning results were found.", ...
    'included_cycles', strings(0, 1), ...
    'excluded_cycles', strings(0, 1), ...
    'orientations', [], ...
    'voltage', struct(), ...
    'voltage_sensitivity', struct(), ...
    'calcium', struct(), ...
    'visualizations', struct('summary_png', "", 'summary_fig', "", 'roi_dir', "", ...
    'voltage_sensitivity_summary_png', "", 'voltage_sensitivity_summary_fig', "", ...
    'voltage_sensitivity_roi_dir', ""), ...
    'result_file', "");

valid = false(numel(cycle_entries), 1);
for idx = 1:numel(cycle_entries)
    tuning = cycle_entries(idx).grating_tuning;
    valid(idx) = isstruct(tuning) ...
        && isfield(tuning, 'voltage') && isfield(tuning, 'calcium') ...
        && has_complete_grating_curves_rec(tuning.voltage) ...
        && has_complete_grating_curves_rec(tuning.calcium);
end
summary.included_cycles = string({cycle_entries(valid).cycle_name})';
summary.excluded_cycles = string({cycle_entries(~valid).cycle_name})';
if ~any(valid)
    return;
end

included_entries = cycle_entries(valid);
orientation_cells = cell(numel(included_entries), 1);
for idx = 1:numel(included_entries)
    orientation_cells{idx} = mod(double(included_entries(idx).grating_tuning.voltage.unique_orientations(:)), 360);
end
orientations = unique(vertcat(orientation_cells{:}));
orientations = sort(orientations(:))';

voltage_roi_count = size(included_entries(1).grating_tuning.voltage.response_by_condition, 1);
calcium_roi_count = size(included_entries(1).grating_tuning.calcium.response_by_condition, 1);
if voltage_roi_count ~= calcium_roi_count
    error('Voltage/calcium ROI count mismatch in grating tuning: voltage=%d calcium=%d.', ...
        voltage_roi_count, calcium_roi_count);
end
for idx = 2:numel(included_entries)
    current_voltage_count = size(included_entries(idx).grating_tuning.voltage.response_by_condition, 1);
    current_calcium_count = size(included_entries(idx).grating_tuning.calcium.response_by_condition, 1);
    if current_voltage_count ~= voltage_roi_count || current_calcium_count ~= calcium_roi_count
        error('ROI count mismatch in grating tuning for cycle %s.', included_entries(idx).cycle_name);
    end
end

summary.status = "completed";
summary.reason = "";
summary.orientations = orientations;
summary.curve_rule = struct( ...
    'orientation_rule', "union after modulo 360; missing directions remain NaN", ...
    'cycle_average_rule', "mean across available cycles with omitnan", ...
    'dsi_mean_rule', "mean of per-cycle DSI", ...
    'dsi_from_mean_curve_rule', "recompute DSI from the cycle-mean grating-present curve", ...
    'nonstim_role', "saved and plotted but excluded from DSI/OSI");
summary.voltage = average_grating_channel_rec(included_entries, 'voltage', orientations, voltage_roi_count);
summary.calcium = average_grating_channel_rec(included_entries, 'calcium', orientations, calcium_roi_count);
has_voltage_sensitivity_tuning = all(arrayfun(@(entry) ...
    isfield(entry.grating_tuning, 'voltage_sensitivity') ...
    && has_complete_grating_curves_rec(entry.grating_tuning.voltage_sensitivity), included_entries));
if has_voltage_sensitivity_tuning
    summary.voltage_sensitivity = average_grating_channel_rec( ...
        included_entries, 'voltage_sensitivity', orientations, voltage_roi_count);
else
    summary.voltage_sensitivity = struct('status', "not_available", ...
        'reason', "At least one included cycle is missing voltage_sensitivity tuning.");
end

summary_png = fullfile(output_dir, '8_grating.png');
summary_fig = fullfile(output_dir, '8_grating.fig');
roi_dir = fullfile(output_dir, '8_grating_roi');
if ~isfolder(roi_dir)
    mkdir(roi_dir);
end
plot_record_grating_summary_rec(summary, summary_fig, summary_png);
plot_record_grating_by_roi_rec(summary, roi_dir);

voltage_sensitivity_summary_png = "";
voltage_sensitivity_summary_fig = "";
voltage_sensitivity_roi_dir = "";
if has_voltage_sensitivity_tuning
    voltage_sensitivity_summary_png = fullfile(output_dir, ...
        '8_grating_v_sens.png');
    voltage_sensitivity_summary_fig = fullfile(output_dir, ...
        '8_grating_v_sens.fig');
    voltage_sensitivity_roi_dir = fullfile(output_dir, ...
        '8_grating_v_sens_roi');
    if ~isfolder(voltage_sensitivity_roi_dir)
        mkdir(voltage_sensitivity_roi_dir);
    end
    plot_record_grating_summary_voltage_sensitivity_rec( ...
        summary, voltage_sensitivity_summary_fig, voltage_sensitivity_summary_png);
    plot_record_grating_by_roi_voltage_sensitivity_rec(summary, voltage_sensitivity_roi_dir);
else
    fprintf('Skipping 8 voltage-sensitivity grating summary: %s\n', summary.voltage_sensitivity.reason);
end

summary.visualizations.summary_png = string(summary_png);
summary.visualizations.summary_fig = string(summary_fig);
summary.visualizations.roi_dir = string(roi_dir);
summary.visualizations.voltage_sensitivity_summary_png = string(voltage_sensitivity_summary_png);
summary.visualizations.voltage_sensitivity_summary_fig = string(voltage_sensitivity_summary_fig);
summary.visualizations.voltage_sensitivity_roi_dir = string(voltage_sensitivity_roi_dir);
summary.result_file = string(fullfile(output_dir, '8_grating.mat'));
summary.metrics_csv = string(fullfile(output_dir, '8_grating_metrics.csv'));
write_record_grating_metrics_table_rec(summary, summary.metrics_csv);
record_grating_tuning = summary;
save(char(summary.result_file), 'record_grating_tuning', '-v7.3');
end

function result = plot_record_random_grating_trial_stacks_rec( ...
    cycle_entries, output_dir, voltage_polarity, calcium_polarity)
is_random = arrayfun(@(entry) strcmpi( ...
    string(entry.stim_analysis_kind), "random_grating_tuning"), cycle_entries);
result = struct( ...
    'status', "not_applicable", ...
    'roi_dir', "", ...
    'voltage_png', "", ...
    'voltage_fig', "", ...
    'calcium_png', "", ...
    'calcium_fig', "", ...
    'voltage_by_roi_dir', "", ...
    'voltage_by_roi_png', strings(0, 1), ...
    'voltage_by_roi_fig', strings(0, 1), ...
    'calcium_by_roi_dir', "", ...
    'calcium_by_roi_png', strings(0, 1), ...
    'calcium_by_roi_fig', strings(0, 1), ...
    'voltage_calcium_stack_tuning_dir', "", ...
    'voltage_calcium_stack_tuning_png', strings(0, 1), ...
    'voltage_calcium_stack_tuning_fig', strings(0, 1), ...
    'voltage_calcium_stack_tuning_rule', "", ...
    'heatmap_dir', "", ...
    'voltage_heatmap_png', "", ...
    'voltage_heatmap_fig', "", ...
    'calcium_heatmap_png', "", ...
    'calcium_heatmap_fig', "", ...
    'voltage_heatmap_by_roi_dir', "", ...
    'voltage_heatmap_by_roi_png', strings(0, 1), ...
    'voltage_heatmap_by_roi_fig', strings(0, 1), ...
    'calcium_heatmap_by_roi_dir', "", ...
    'calcium_heatmap_by_roi_png', strings(0, 1), ...
    'calcium_heatmap_by_roi_fig', strings(0, 1), ...
    'voltage_calcium_subplot_dir', "", ...
    'voltage_calcium_subplot_png', strings(0, 1), ...
    'voltage_calcium_subplot_fig', strings(0, 1), ...
    'voltage_calcium_subplot_emf', strings(0, 1), ...
    'cycle_dsi_comparison_dir', "", ...
    'cycle_dsi_comparison_png', "", ...
    'cycle_dsi_comparison_fig', "", ...
    'cycle_dsi_comparison_csv', "", ...
    'cycle_dsi_comparison_mat', "", ...
    'cycle_dsi_by_roi_dir', "", ...
    'cycle_dsi_by_roi_png', strings(0, 1), ...
    'cycle_dsi_by_roi_fig', strings(0, 1), ...
    'cycle_dsi_by_roi_csv', "", ...
    'cycle_dsi_voltage_sensitivity_dir', "", ...
    'cycle_dsi_voltage_sensitivity_png', "", ...
    'cycle_dsi_voltage_sensitivity_fig', "", ...
    'cycle_dsi_voltage_sensitivity_csv', "", ...
    'cycle_dsi_voltage_sensitivity_mat', "", ...
    'cycle_dsi_voltage_sensitivity_by_roi_dir', "", ...
    'cycle_dsi_voltage_sensitivity_by_roi_png', strings(0, 1), ...
    'cycle_dsi_voltage_sensitivity_by_roi_fig', strings(0, 1), ...
    'cycle_dsi_voltage_sensitivity_by_roi_csv', "", ...
    'cycle_tuning_curve_dir', "", ...
    'cycle_tuning_curve_png', strings(0, 1), ...
    'cycle_tuning_curve_fig', strings(0, 1), ...
    'cycle_tuning_curve_voltage_sensitivity_png', strings(0, 1), ...
    'cycle_tuning_curve_voltage_sensitivity_fig', strings(0, 1), ...
    'cycle_tuning_curve_csv', "", ...
    'cycle_tuning_curve_mat', "", ...
    'strongest_voltage_dsi_cycle', "", ...
    'angles', [], ...
    'trial_count_by_roi_angle', [], ...
    'trace_stage', "sensitivity", ...
    'alignment_rule', "full pre-stimulus ISI baseline_frames plus current stimulus duration stim_frames; stimulus onset at t=0", ...
    'layout_rule', "columns are directions; rows are ROIs; voltage and calcium use separate figures", ...
    'ordering_rule', "within each ROI-direction tile: cycle order, then original trial order", ...
    'display_rule', "native sensitivity values with polarity applied; vertical offsets only");
if ~any(is_random)
    return;
end
if ~all(is_random)
    error('Random and non-random grating cycles cannot share one random-grating trial-stack output.');
end

nrois = size(cycle_entries(1).voltage.sensitivity.data, 2);
for cycle_idx = 1:numel(cycle_entries)
    if size(cycle_entries(cycle_idx).voltage.sensitivity.data, 2) ~= nrois ...
            || size(cycle_entries(cycle_idx).calcium.sensitivity.data, 2) ~= nrois
        error('ROI count mismatch while building random-grating trial stacks for cycle %s.', ...
            cycle_entries(cycle_idx).cycle_name);
    end
end

angle_cells = arrayfun(@(entry) mod(double(entry.stim_windows.orientations(:)), 360), ...
    cycle_entries, 'UniformOutput', false);
angles = unique(vertcat(angle_cells{:}), 'sorted')';
roi_dir = fullfile(output_dir, '9_rg_stack');
if ~isfolder(roi_dir)
    mkdir(roi_dir);
end

voltage_fig = string(fullfile(roi_dir, 'v_stack.fig'));
voltage_png = string(fullfile(roi_dir, 'v_stack.png'));
calcium_fig = string(fullfile(roi_dir, 'ca_stack.fig'));
calcium_png = string(fullfile(roi_dir, 'ca_stack.png'));
voltage_trial_count = plot_random_grating_trial_stack_grid_rec( ...
    cycle_entries, 'voltage', angles, nrois, voltage_polarity, voltage_fig, voltage_png);
calcium_trial_count = plot_random_grating_trial_stack_grid_rec( ...
    cycle_entries, 'calcium', angles, nrois, calcium_polarity, calcium_fig, calcium_png);
if ~isequal(voltage_trial_count, calcium_trial_count)
    error('Voltage/calcium trial-count grids differ in random-grating stack output.');
end
by_roi_root = fullfile(roi_dir, 'roi');
[voltage_by_roi_dir, voltage_by_roi_fig, voltage_by_roi_png] = ...
    plot_random_grating_trial_stacks_by_roi_rec( ...
        cycle_entries, 'voltage', angles, nrois, voltage_polarity, by_roi_root);
[calcium_by_roi_dir, calcium_by_roi_fig, calcium_by_roi_png] = ...
    plot_random_grating_trial_stacks_by_roi_rec( ...
        cycle_entries, 'calcium', angles, nrois, calcium_polarity, by_roi_root);
[voltage_calcium_stack_tuning_dir, voltage_calcium_stack_tuning_fig, ...
    voltage_calcium_stack_tuning_png] = ...
    plot_random_grating_roi_voltage_calcium_stack_tuning_rec( ...
        cycle_entries, angles, nrois, voltage_polarity, calcium_polarity, roi_dir);

heatmap_dir = fullfile(output_dir, '10_rg_hm');
if ~isfolder(heatmap_dir)
    mkdir(heatmap_dir);
end
voltage_heatmap_fig = string(fullfile(heatmap_dir, ...
    'v_hm.fig'));
voltage_heatmap_png = string(fullfile(heatmap_dir, ...
    'v_hm.png'));
calcium_heatmap_fig = string(fullfile(heatmap_dir, ...
    'ca_hm.fig'));
calcium_heatmap_png = string(fullfile(heatmap_dir, ...
    'ca_hm.png'));
if random_grating_heatmap_outputs_complete_rec(heatmap_dir, 'voltage', nrois)
    fprintf('Voltage random-grating heatmap outputs already complete; skipping heatmap color-limit scan.\n');
    voltage_heatmap_clim = [NaN, NaN];
else
    voltage_heatmap_clim = resolve_random_grating_heatmap_clim_rec( ...
        cycle_entries, 'voltage', angles, nrois, voltage_polarity);
end
if random_grating_heatmap_outputs_complete_rec(heatmap_dir, 'calcium', nrois)
    fprintf('Calcium random-grating heatmap outputs already complete; skipping heatmap color-limit scan.\n');
    calcium_heatmap_clim = [NaN, NaN];
else
    calcium_heatmap_clim = resolve_random_grating_heatmap_clim_rec( ...
        cycle_entries, 'calcium', angles, nrois, calcium_polarity);
end
plot_random_grating_heatmap_grid_rec( ...
    cycle_entries, 'voltage', angles, nrois, voltage_polarity, ...
    voltage_heatmap_clim, voltage_heatmap_fig, voltage_heatmap_png);
plot_random_grating_heatmap_grid_rec( ...
    cycle_entries, 'calcium', angles, nrois, calcium_polarity, ...
    calcium_heatmap_clim, calcium_heatmap_fig, calcium_heatmap_png);
heatmap_by_roi_root = fullfile(heatmap_dir, 'roi');
[voltage_heatmap_by_roi_dir, voltage_heatmap_by_roi_fig, voltage_heatmap_by_roi_png] = ...
    plot_random_grating_heatmaps_by_roi_rec( ...
        cycle_entries, 'voltage', angles, nrois, voltage_polarity, ...
        voltage_heatmap_clim, heatmap_by_roi_root);
[calcium_heatmap_by_roi_dir, calcium_heatmap_by_roi_fig, calcium_heatmap_by_roi_png] = ...
    plot_random_grating_heatmaps_by_roi_rec( ...
        cycle_entries, 'calcium', angles, nrois, calcium_polarity, ...
        calcium_heatmap_clim, heatmap_by_roi_root);

subplot_root = fullfile(output_dir, '11_rg_v_ca');
[subplot_dir, subplot_fig, subplot_png, subplot_emf] = ...
    plot_random_grating_roi_cycle_voltage_calcium_subplots_rec( ...
        cycle_entries, nrois, voltage_polarity, calcium_polarity, subplot_root);
dsi_comparison = build_random_grating_cycle_voltage_dsi_comparison_rec( ...
    cycle_entries, angles, nrois, subplot_root);
cycle_tuning = plot_random_grating_cycle_tuning_curves_rec( ...
    cycle_entries, angles, nrois, subplot_root);

result.status = "completed";
result.roi_dir = string(roi_dir);
result.voltage_png = voltage_png;
result.voltage_fig = voltage_fig;
result.calcium_png = calcium_png;
result.calcium_fig = calcium_fig;
result.voltage_by_roi_dir = voltage_by_roi_dir;
result.voltage_by_roi_png = voltage_by_roi_png;
result.voltage_by_roi_fig = voltage_by_roi_fig;
result.calcium_by_roi_dir = calcium_by_roi_dir;
result.calcium_by_roi_png = calcium_by_roi_png;
result.calcium_by_roi_fig = calcium_by_roi_fig;
result.voltage_calcium_stack_tuning_dir = voltage_calcium_stack_tuning_dir;
result.voltage_calcium_stack_tuning_png = voltage_calcium_stack_tuning_png;
result.voltage_calcium_stack_tuning_fig = voltage_calcium_stack_tuning_fig;
result.voltage_calcium_stack_tuning_rule = ...
    "one figure per ROI; all cycles/trials pooled by direction; tuning uses per-trial stim and prestimulus means with SEM";
result.heatmap_dir = string(heatmap_dir);
result.voltage_heatmap_png = voltage_heatmap_png;
result.voltage_heatmap_fig = voltage_heatmap_fig;
result.calcium_heatmap_png = calcium_heatmap_png;
result.calcium_heatmap_fig = calcium_heatmap_fig;
result.voltage_heatmap_by_roi_dir = voltage_heatmap_by_roi_dir;
result.voltage_heatmap_by_roi_png = voltage_heatmap_by_roi_png;
result.voltage_heatmap_by_roi_fig = voltage_heatmap_by_roi_fig;
result.calcium_heatmap_by_roi_dir = calcium_heatmap_by_roi_dir;
result.calcium_heatmap_by_roi_png = calcium_heatmap_by_roi_png;
result.calcium_heatmap_by_roi_fig = calcium_heatmap_by_roi_fig;
result.voltage_calcium_subplot_dir = subplot_dir;
result.voltage_calcium_subplot_png = subplot_png;
result.voltage_calcium_subplot_fig = subplot_fig;
result.voltage_calcium_subplot_emf = subplot_emf;
result.cycle_dsi_comparison_dir = dsi_comparison.output_dir;
result.cycle_dsi_comparison_png = dsi_comparison.summary_png;
result.cycle_dsi_comparison_fig = dsi_comparison.summary_fig;
result.cycle_dsi_comparison_csv = dsi_comparison.summary_csv;
result.cycle_dsi_comparison_mat = dsi_comparison.summary_mat;
result.cycle_dsi_by_roi_dir = dsi_comparison.by_roi_dir;
result.cycle_dsi_by_roi_png = dsi_comparison.by_roi_png;
result.cycle_dsi_by_roi_fig = dsi_comparison.by_roi_fig;
result.cycle_dsi_by_roi_csv = dsi_comparison.roi_summary_csv;
if isfield(dsi_comparison, 'voltage_sensitivity') && isstruct(dsi_comparison.voltage_sensitivity)
    result.cycle_dsi_voltage_sensitivity_dir = dsi_comparison.voltage_sensitivity.output_dir;
    result.cycle_dsi_voltage_sensitivity_png = dsi_comparison.voltage_sensitivity.summary_png;
    result.cycle_dsi_voltage_sensitivity_fig = dsi_comparison.voltage_sensitivity.summary_fig;
    result.cycle_dsi_voltage_sensitivity_csv = dsi_comparison.voltage_sensitivity.summary_csv;
    result.cycle_dsi_voltage_sensitivity_mat = dsi_comparison.voltage_sensitivity.summary_mat;
    result.cycle_dsi_voltage_sensitivity_by_roi_dir = dsi_comparison.voltage_sensitivity.by_roi_dir;
    result.cycle_dsi_voltage_sensitivity_by_roi_png = dsi_comparison.voltage_sensitivity.by_roi_png;
    result.cycle_dsi_voltage_sensitivity_by_roi_fig = dsi_comparison.voltage_sensitivity.by_roi_fig;
    result.cycle_dsi_voltage_sensitivity_by_roi_csv = dsi_comparison.voltage_sensitivity.roi_summary_csv;
end
result.cycle_tuning_curve_dir = cycle_tuning.output_dir;
result.cycle_tuning_curve_png = cycle_tuning.png_files;
result.cycle_tuning_curve_fig = cycle_tuning.fig_files;
result.cycle_tuning_curve_voltage_sensitivity_png = cycle_tuning.voltage_sensitivity_png_files;
result.cycle_tuning_curve_voltage_sensitivity_fig = cycle_tuning.voltage_sensitivity_fig_files;
result.cycle_tuning_curve_csv = cycle_tuning.summary_csv;
result.cycle_tuning_curve_mat = cycle_tuning.summary_mat;
result.strongest_voltage_dsi_cycle = dsi_comparison.strongest_cycle_name;
result.strongest_voltage_dsi_mean = dsi_comparison.strongest_cycle_mean_dsi;
result.heatmap_color_rule = "per-channel global symmetric limits from all extracted trial-window sensitivity values";
result.heatmap_voltage_clim = voltage_heatmap_clim;
result.heatmap_calcium_clim = calcium_heatmap_clim;
result.angles = angles;
result.trial_count_by_roi_angle = voltage_trial_count;
fprintf('Random-grating strongest voltage DSI cycle: %s (mean DSI %.4g).\n', ...
    char(dsi_comparison.strongest_cycle_name), dsi_comparison.strongest_cycle_mean_dsi);
end

function trial_count = plot_random_grating_trial_stack_grid_rec( ...
    cycle_entries, channel_name, angles, nrois, polarity, fig_file, png_file)
nangles = numel(angles);
trial_count = compute_random_grating_trial_count_by_angle_rec(cycle_entries, angles, nrois);
if skip_existing_figure_bundle_rec(fig_file, png_file, ...
        sprintf('%s random-grating trial trace stack grid', channel_name))
    return;
end
fig_width = min(2600, max(1400, 380 * nangles));
fig_height = min(2800, max(850, 300 * nrois));
fig = figure('Color', 'w', ...
    'Name', sprintf('Random Grating %s Trial Trace Stack Grid', channel_name), ...
    'Position', [40, 40, fig_width, fig_height]);
tiledlayout(fig, nrois, nangles, 'Padding', 'compact', 'TileSpacing', 'compact');

for roi_idx = 1:nrois
    all_trials = collect_random_grating_roi_trials_rec( ...
        cycle_entries, channel_name, roi_idx, angles, polarity);
    for angle_idx = 1:nangles
        current_angle = angles(angle_idx);
        tile_trials = all_trials(abs([all_trials.angle] - current_angle) < 1e-9);

        ax = nexttile;
        if roi_idx == 1
            tile_title = sprintf('%g deg', current_angle);
        else
            tile_title = '';
        end
        plot_random_grating_trial_stack_channel_rec(ax, tile_trials, current_angle, tile_title);
        ylabel(ax, sprintf('ROI %03d', roi_idx));
        if roi_idx < nrois
            xlabel(ax, '');
        end
    end
end

sgtitle(fig, sprintf('%s Sensitivity | Random Grating Trials | Columns: Direction, Rows: ROI', ...
    upper_first_rec(channel_name)));
save_figure_bundle_preserve_layout_rec(fig, fig_file, png_file);
close(fig);
end

function trial_count = compute_random_grating_trial_count_by_angle_rec(cycle_entries, angles, nrois)
trial_count_one_roi = zeros(1, numel(angles));
for cycle_idx = 1:numel(cycle_entries)
    cycle_angles = mod(double(cycle_entries(cycle_idx).stim_windows.orientations(:)), 360);
    for angle_idx = 1:numel(angles)
        trial_count_one_roi(angle_idx) = trial_count_one_roi(angle_idx) + ...
            sum(abs(cycle_angles - angles(angle_idx)) < 1e-9);
    end
end
trial_count = repmat(trial_count_one_roi, nrois, 1);
end

function [channel_dir, fig_files, png_files] = plot_random_grating_trial_stacks_by_roi_rec( ...
    cycle_entries, channel_name, angles, nrois, polarity, by_roi_root)
channel_dir = string(fullfile(by_roi_root, char(string(channel_name))));
if ~isfolder(channel_dir)
    mkdir(channel_dir);
end
nangles = numel(angles);
fig_files = strings(nrois, 1);
png_files = strings(nrois, 1);
fig_width = min(2600, max(1400, 380 * nangles));

for roi_idx = 1:nrois
    fig_files(roi_idx) = string(fullfile(channel_dir, sprintf( ...
        'r%03d_%s_stack.fig', roi_idx, channel_name)));
    png_files(roi_idx) = string(fullfile(channel_dir, sprintf( ...
        'r%03d_%s_stack.png', roi_idx, channel_name)));
    if skip_existing_figure_bundle_rec(fig_files(roi_idx), png_files(roi_idx), ...
            sprintf('%s random-grating ROI %03d trial traces', channel_name, roi_idx))
        continue;
    end

    all_trials = collect_random_grating_roi_trials_rec( ...
        cycle_entries, channel_name, roi_idx, angles, polarity);
    fig = figure('Color', 'w', ...
        'Name', sprintf('Random Grating %s ROI %03d Trial Traces', channel_name, roi_idx), ...
        'Position', [40, 80, fig_width, 850]);
    tiledlayout(fig, 1, nangles, 'Padding', 'compact', 'TileSpacing', 'compact');
    for angle_idx = 1:nangles
        current_angle = angles(angle_idx);
        tile_trials = all_trials(abs([all_trials.angle] - current_angle) < 1e-9);
        ax = nexttile;
        plot_random_grating_trial_stack_channel_rec( ...
            ax, tile_trials, current_angle, sprintf('%g deg', current_angle));
        ylabel(ax, 'Trials');
    end
    sgtitle(fig, sprintf('%s Sensitivity | ROI %03d | Columns: Direction', ...
        upper_first_rec(channel_name), roi_idx));

    save_figure_bundle_preserve_layout_rec(fig, fig_files(roi_idx), png_files(roi_idx));
    close(fig);
end
end

function [output_dir, fig_files, png_files] = ...
    plot_random_grating_roi_voltage_calcium_stack_tuning_rec( ...
    cycle_entries, angles, nrois, voltage_polarity, calcium_polarity, stack_dir)
output_dir = string(fullfile(stack_dir, 'v_ca'));
if ~isfolder(output_dir)
    mkdir(output_dir);
end

nangles = numel(angles);
fig_files = strings(nrois, 1);
png_files = strings(nrois, 1);
voltage_color = [1.00 0.20 0.75];
calcium_color = [0.32 0.65 0.12];

for roi_idx = 1:nrois
    fig_files(roi_idx) = string(fullfile(output_dir, sprintf( ...
        'r%03d_v_ca_stack_tuning.fig', roi_idx)));
    png_files(roi_idx) = string(fullfile(output_dir, sprintf( ...
        'r%03d_v_ca_stack_tuning.png', roi_idx)));
    if skip_existing_figure_bundle_rec(fig_files(roi_idx), png_files(roi_idx), ...
            sprintf('ROI %03d random-grating voltage/calcium stack and tuning', roi_idx))
        continue;
    end

    voltage_trials = collect_random_grating_roi_trials_rec( ...
        cycle_entries, 'voltage', roi_idx, angles, voltage_polarity);
    calcium_trials = collect_random_grating_roi_trials_rec( ...
        cycle_entries, 'calcium', roi_idx, angles, calcium_polarity);
    voltage_bar = resolve_trial_stack_scalebar_rec(voltage_trials);
    calcium_bar = resolve_trial_stack_scalebar_rec(calcium_trials);
    max_trial_count = max(numel(voltage_trials), numel(calcium_trials));
    fig_width = min(3000, max(1800, 340 * nangles + 480));
    fig_height = min(2000, max(950, 36 * max_trial_count + 480));
    fig = figure('Color', 'w', ...
        'Name', sprintf('ROI %03d Random Grating Voltage Calcium Stack Tuning', roi_idx), ...
        'Position', [40, 40, fig_width, fig_height]);
    tiledlayout(fig, 2, nangles + 1, 'Padding', 'compact', 'TileSpacing', 'compact');

    for angle_idx = 1:nangles
        target_angle = angles(angle_idx);
        voltage_angle_trials = voltage_trials(abs([voltage_trials.angle] - target_angle) < 1e-9);
        calcium_angle_trials = calcium_trials(abs([calcium_trials.angle] - target_angle) < 1e-9);

        ax_voltage = nexttile(angle_idx);
        plot_random_grating_trial_stack_summary_channel_rec( ...
            ax_voltage, voltage_angle_trials, voltage_color, 1.0, voltage_bar, ...
            sprintf('%g deg', target_angle), angle_idx == 1, false);

        ax_calcium = nexttile(nangles + 1 + angle_idx);
        plot_random_grating_trial_stack_summary_channel_rec( ...
            ax_calcium, calcium_angle_trials, calcium_color, 0.3, calcium_bar, ...
            '', angle_idx == 1, angle_idx == 1);
    end

    ax_voltage_tuning = nexttile(nangles + 1);
    plot_random_grating_pooled_tuning_rec( ...
        ax_voltage_tuning, voltage_trials, angles, voltage_color, 'Voltage tuning');
    ax_calcium_tuning = nexttile(2 * (nangles + 1));
    plot_random_grating_pooled_tuning_rec( ...
        ax_calcium_tuning, calcium_trials, angles, calcium_color, 'Calcium tuning');
    sgtitle(fig, sprintf(['ROI %03d | Random grating | direct pooled trials across cycles | ' ...
        'top: Voltage, bottom: Calcium'], roi_idx));
    save_figure_bundle_preserve_layout_rec(fig, fig_files(roi_idx), png_files(roi_idx));
    close(fig);
end
end

function scale_length = resolve_trial_stack_scalebar_rec(trials)
if isempty(trials)
    scale_length = 1;
    return;
end
values = vertcat(trials.trace);
values = values(isfinite(values));
if isempty(values)
    scale_length = 1;
else
    scale_length = nice_scalebar_value_rec(0.18 * range(values));
end
end

function plot_random_grating_trial_stack_summary_channel_rec( ...
    ax, trials, line_color, spacing_fraction, scale_length, title_text, show_y_scalebar, show_x_scalebar)
if isempty(trials)
    text(ax, 0.5, 0.5, 'No valid trials', ...
        'Units', 'normalized', 'HorizontalAlignment', 'center');
    axis(ax, 'off');
    return;
end

spacing = spacing_fraction * compute_trial_cell_stack_spacing_rec({trials.trace});
ntrials = numel(trials);
offsets = (ntrials - (1:ntrials)) * spacing;
all_time = vertcat(trials.time);
all_trace = vertcat(trials.trace);
x_limits = [min(all_time), max(all_time)];
y_limits = [min(all_trace) - 0.5 * spacing, ...
    max(all_trace) + offsets(1) + 0.5 * spacing];
if y_limits(2) <= y_limits(1)
    y_limits = y_limits + [-0.5, 0.5];
end

hold(ax, 'on');
xlim(ax, x_limits);
ylim(ax, y_limits);
stim_end = max(cellfun(@max, {trials.time}));
if isfinite(stim_end) && stim_end > 0
    patch(ax, [0, stim_end, stim_end, 0], ...
        [y_limits(1), y_limits(1), y_limits(2), y_limits(2)], ...
        [0.55 0.55 0.55], 'FaceAlpha', 0.45, 'EdgeColor', 'none');
end
for trial_idx = 1:ntrials
    plot(ax, trials(trial_idx).time, trials(trial_idx).trace + offsets(trial_idx), ...
        'Color', line_color, 'LineWidth', 0.75);
end
xline(ax, 0, '-k', 'LineWidth', 1.0);
set(ax, 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none', 'TickDir', 'out');
title(ax, title_text, 'FontWeight', 'bold');
box(ax, 'off');

if show_y_scalebar
    add_trial_stack_y_scalebar_rec(ax, scale_length, [0.08 0.08 0.08], 'Sensitivity');
end
if show_x_scalebar
    add_trial_stack_x_scalebar_rec(ax, 2, [0.08 0.08 0.08]);
end
end

function add_trial_stack_y_scalebar_rec(ax, scale_length, bar_color, y_suffix)
x_limits = xlim(ax);
y_limits = ylim(ax);
x_span = max(eps, diff(x_limits));
y_span = max(eps, diff(y_limits));
x_start = x_limits(1) + 0.08 * x_span;
y_start = y_limits(1) + 0.10 * y_span;
plot(ax, [x_start, x_start], [y_start, y_start + scale_length], ...
    'Color', bar_color, 'LineWidth', 1.4, 'Clipping', 'off');
text(ax, x_start - 0.02 * x_span, y_start + scale_length / 2, ...
    sprintf('%s %s', format_scalebar_value_rec(scale_length), y_suffix), ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
    'Rotation', 90, 'FontSize', 8, 'Color', bar_color, 'Interpreter', 'none');
end

function add_trial_stack_x_scalebar_rec(ax, scale_length, bar_color)
x_limits = xlim(ax);
y_limits = ylim(ax);
x_span = max(eps, diff(x_limits));
y_span = max(eps, diff(y_limits));
x_start = x_limits(2) - 0.08 * x_span - scale_length;
y_start = y_limits(1) + 0.10 * y_span;
plot(ax, [x_start, x_start + scale_length], [y_start, y_start], ...
    'Color', bar_color, 'LineWidth', 1.4, 'Clipping', 'off');
text(ax, x_start + scale_length / 2, y_start - 0.04 * y_span, ...
    sprintf('%s s', format_scalebar_value_rec(scale_length)), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
    'FontSize', 8, 'Color', bar_color, 'Interpreter', 'none');
end

function plot_random_grating_pooled_tuning_rec(ax, trials, angles, line_color, title_text)
[stim_mean, stim_sem, nonstim_mean, nonstim_sem] = ...
    summarize_random_grating_trial_tuning_rec(trials, angles);
hold(ax, 'on');
errorbar(ax, angles, nonstim_mean, nonstim_sem, '--o', ...
    'Color', lighten_color(line_color, 0.45), 'MarkerFaceColor', 'w', 'LineWidth', 1.1);
errorbar(ax, angles, stim_mean, stim_sem, '-o', ...
    'Color', line_color, 'MarkerFaceColor', line_color, 'LineWidth', 1.4);
xlabel(ax, 'Direction (deg)');
ylabel(ax, 'Mean sensitivity');
title(ax, title_text);
legend(ax, {'Non-grating', 'Grating'}, 'Location', 'best');
xlim(ax, [-45, 315]);
set(ax, 'XTick', [0, 90, 180, 270], 'TickDir', 'out');
grid(ax, 'off');
box(ax, 'off');
end

function [stim_mean, stim_sem, nonstim_mean, nonstim_sem] = ...
    summarize_random_grating_trial_tuning_rec(trials, angles)
stim_mean = NaN(size(angles));
stim_sem = NaN(size(angles));
nonstim_mean = NaN(size(angles));
nonstim_sem = NaN(size(angles));
for angle_idx = 1:numel(angles)
    angle_trials = trials(abs([trials.angle] - angles(angle_idx)) < 1e-9);
    stim_values = NaN(numel(angle_trials), 1);
    nonstim_values = NaN(numel(angle_trials), 1);
    for trial_idx = 1:numel(angle_trials)
        trace = double(angle_trials(trial_idx).trace(:));
        time = double(angle_trials(trial_idx).time(:));
        stim_values(trial_idx) = mean(trace(time >= 0), 'omitnan');
        nonstim_values(trial_idx) = mean(trace(time < 0), 'omitnan');
    end
    stim_mean(angle_idx) = mean(stim_values, 'omitnan');
    nonstim_mean(angle_idx) = mean(nonstim_values, 'omitnan');
    stim_sem(angle_idx) = standard_error_rec(stim_values);
    nonstim_sem(angle_idx) = standard_error_rec(nonstim_values);
end
end

function value = standard_error_rec(values)
values = double(values(isfinite(values)));
if numel(values) <= 1
    value = 0;
else
    value = std(values, 0) / sqrt(numel(values));
end
end

function color_limit = resolve_random_grating_heatmap_clim_rec( ...
    cycle_entries, channel_name, angles, nrois, polarity)
max_abs_value = 0;
for roi_idx = 1:nrois
    trials = collect_random_grating_roi_trials_rec( ...
        cycle_entries, channel_name, roi_idx, angles, polarity);
    for trial_idx = 1:numel(trials)
        trial_max = max(abs(double(trials(trial_idx).trace)), [], 'omitnan');
        if isfinite(trial_max)
            max_abs_value = max(max_abs_value, trial_max);
        end
    end
end
if ~isfinite(max_abs_value) || max_abs_value <= 0
    max_abs_value = 1;
end
color_limit = [-max_abs_value, max_abs_value];
end

function tf = random_grating_heatmap_outputs_complete_rec(heatmap_dir, channel_name, nrois)
channel_name = char(string(channel_name));
grid_fig = string(fullfile(heatmap_dir, sprintf('%s_hm.fig', channel_name)));
grid_png = string(fullfile(heatmap_dir, sprintf('%s_hm.png', channel_name)));
tf = figure_bundle_complete_rec(grid_fig, grid_png);
if ~tf
    return;
end
channel_dir = fullfile(heatmap_dir, 'roi', channel_name);
for roi_idx = 1:nrois
    roi_fig = string(fullfile(channel_dir, sprintf( ...
        'r%03d_%s_hm.fig', roi_idx, channel_name)));
    roi_png = string(fullfile(channel_dir, sprintf( ...
        'r%03d_%s_hm.png', roi_idx, channel_name)));
    if ~figure_bundle_complete_rec(roi_fig, roi_png)
        tf = false;
        return;
    end
end
end

function plot_random_grating_heatmap_grid_rec( ...
    cycle_entries, channel_name, angles, nrois, polarity, color_limit, fig_file, png_file)
if skip_existing_figure_bundle_rec(fig_file, png_file, ...
        sprintf('%s random-grating sensitivity heatmap grid', channel_name))
    return;
end
nangles = numel(angles);
fig_width = min(2600, max(1400, 380 * nangles));
fig_height = min(2800, max(850, 300 * nrois));
fig = figure('Color', 'w', ...
    'Name', sprintf('Random Grating %s Sensitivity Heatmap Grid', channel_name), ...
    'Position', [40, 40, fig_width, fig_height]);
tiledlayout(fig, nrois, nangles, 'Padding', 'compact', 'TileSpacing', 'compact');
last_ax = [];

for roi_idx = 1:nrois
    all_trials = collect_random_grating_roi_trials_rec( ...
        cycle_entries, channel_name, roi_idx, angles, polarity);
    for angle_idx = 1:nangles
        current_angle = angles(angle_idx);
        tile_trials = all_trials(abs([all_trials.angle] - current_angle) < 1e-9);
        last_ax = nexttile;
        if roi_idx == 1
            tile_title = sprintf('%g deg', current_angle);
        else
            tile_title = '';
        end
        plot_random_grating_trial_heatmap_channel_rec( ...
            last_ax, tile_trials, color_limit, tile_title);
        ylabel(last_ax, sprintf('ROI %03d', roi_idx));
        if roi_idx < nrois
            xlabel(last_ax, '');
        end
    end
end

colormap(fig, redblue_colormap_rec(256));
if ~isempty(last_ax)
    cb = colorbar(last_ax);
    cb.Layout.Tile = 'east';
    cb.Label.String = sprintf('%s sensitivity', upper_first_rec(channel_name));
end
sgtitle(fig, sprintf('%s Sensitivity Heatmap | Columns: Direction, Rows: ROI', ...
    upper_first_rec(channel_name)));
save_figure_bundle_preserve_layout_rec(fig, fig_file, png_file);
close(fig);
end

function [channel_dir, fig_files, png_files] = plot_random_grating_heatmaps_by_roi_rec( ...
    cycle_entries, channel_name, angles, nrois, polarity, color_limit, by_roi_root)
channel_dir = string(fullfile(by_roi_root, char(string(channel_name))));
if ~isfolder(channel_dir)
    mkdir(channel_dir);
end
nangles = numel(angles);
fig_files = strings(nrois, 1);
png_files = strings(nrois, 1);
fig_width = min(2600, max(1400, 380 * nangles));

for roi_idx = 1:nrois
    fig_files(roi_idx) = string(fullfile(channel_dir, sprintf( ...
        'r%03d_%s_hm.fig', roi_idx, channel_name)));
    png_files(roi_idx) = string(fullfile(channel_dir, sprintf( ...
        'r%03d_%s_hm.png', roi_idx, channel_name)));
    if skip_existing_figure_bundle_rec(fig_files(roi_idx), png_files(roi_idx), ...
            sprintf('%s random-grating ROI %03d sensitivity heatmap', channel_name, roi_idx))
        continue;
    end

    all_trials = collect_random_grating_roi_trials_rec( ...
        cycle_entries, channel_name, roi_idx, angles, polarity);
    fig = figure('Color', 'w', ...
        'Name', sprintf('Random Grating %s ROI %03d Sensitivity Heatmap', channel_name, roi_idx), ...
        'Position', [40, 80, fig_width, 850]);
    tiledlayout(fig, 1, nangles, 'Padding', 'compact', 'TileSpacing', 'compact');
    last_ax = [];
    for angle_idx = 1:nangles
        current_angle = angles(angle_idx);
        tile_trials = all_trials(abs([all_trials.angle] - current_angle) < 1e-9);
        last_ax = nexttile;
        plot_random_grating_trial_heatmap_channel_rec( ...
            last_ax, tile_trials, color_limit, sprintf('%g deg', current_angle));
        ylabel(last_ax, 'Trials');
    end
    colormap(fig, redblue_colormap_rec(256));
    if ~isempty(last_ax)
        cb = colorbar(last_ax);
        cb.Layout.Tile = 'east';
        cb.Label.String = sprintf('%s sensitivity', upper_first_rec(channel_name));
    end
    sgtitle(fig, sprintf('%s Sensitivity Heatmap | ROI %03d | Columns: Direction', ...
        upper_first_rec(channel_name), roi_idx));

    save_figure_bundle_preserve_layout_rec(fig, fig_files(roi_idx), png_files(roi_idx));
    close(fig);
end
end

function [subplot_dir, fig_files, png_files, emf_files] = ...
    plot_random_grating_roi_cycle_voltage_calcium_subplots_rec( ...
    cycle_entries, nrois, voltage_polarity, calcium_polarity, subplot_root)
subplot_dir = string(fullfile(subplot_root, 'cyc'));
if ~isfolder(subplot_dir)
    mkdir(subplot_dir);
end

ncycles = numel(cycle_entries);
fig_files = strings(nrois, ncycles);
png_files = strings(nrois, ncycles);
emf_files = strings(nrois, ncycles);

for roi_idx = 1:nrois
    roi_dir = fullfile(subplot_dir, sprintf('r%03d', roi_idx));
    if ~isfolder(roi_dir)
        mkdir(roi_dir);
    end

    for cycle_idx = 1:ncycles
        cycle_name = string(cycle_entries(cycle_idx).cycle_name);
        cycle_tag = sprintf('c%02d', cycle_idx);
        fig_files(roi_idx, cycle_idx) = string(fullfile(roi_dir, sprintf( ...
            'r%03d_%s_v_ca_sens.fig', ...
            roi_idx, cycle_tag)));
        png_files(roi_idx, cycle_idx) = string(fullfile(roi_dir, sprintf( ...
            'r%03d_%s_v_ca_sens.png', ...
            roi_idx, cycle_tag)));
        emf_files(roi_idx, cycle_idx) = replace_file_extension_rec(png_files(roi_idx, cycle_idx), '.emf');
        if skip_existing_figure_bundle_rec(fig_files(roi_idx, cycle_idx), png_files(roi_idx, cycle_idx), ...
                sprintf('ROI %03d %s random-grating continuous trace', roi_idx, char(cycle_name)))
            continue;
        end

        voltage_trace = double(voltage_polarity) * ...
            double(cycle_entries(cycle_idx).voltage.sensitivity.data(:, roi_idx));
        calcium_trace = double(calcium_polarity) * ...
            double(cycle_entries(cycle_idx).calcium.sensitivity.data(:, roi_idx));

        fig = figure('Color', 'w', ...
            'Name', sprintf('ROI %03d %s Random Grating Voltage Calcium Sensitivity', roi_idx, char(cycle_name)), ...
            'Position', [80, 80, 1500, 900]);
        tiledlayout(fig, 2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

        ax_voltage = nexttile;
        plot_random_grating_roi_cycle_continuous_trace_rec( ...
            ax_voltage, cycle_entries(cycle_idx).voltage.sensitivity.time, voltage_trace, ...
            cycle_entries(cycle_idx).stim_windows, 'voltage', [0.85 0.12 0.10], ...
            sprintf('ROI %03d | %s | Voltage sensitivity', roi_idx, char(cycle_name)), ...
            'Voltage sensitivity');

        ax_calcium = nexttile;
        plot_random_grating_roi_cycle_continuous_trace_rec( ...
            ax_calcium, cycle_entries(cycle_idx).calcium.sensitivity.time, calcium_trace, ...
            cycle_entries(cycle_idx).stim_windows, 'calcium', [0.10 0.55 0.18], ...
            sprintf('ROI %03d | %s | Calcium sensitivity', roi_idx, char(cycle_name)), ...
            'Calcium sensitivity');
        xlabel(ax_calcium, 'Time in Cycle (s)');

        sgtitle(fig, sprintf(['Random drifting grating | ROI %03d | %s | ' ...
            'continuous trace with drifting directions'], roi_idx, char(cycle_name)), ...
            'Interpreter', 'none');
        save_figure_bundle_preserve_layout_rec(fig, fig_files(roi_idx, cycle_idx), png_files(roi_idx, cycle_idx));
        close(fig);
    end
end
end

function plot_random_grating_roi_cycle_continuous_trace_rec( ...
    ax, t_axis, trace_data, stim_windows, channel_name, trace_color, title_text, y_axis_label)
t_axis = double(t_axis(:));
trace_data = double(trace_data(:));
hold(ax, 'on');

[ymin, ymax] = compute_trace_limits_rec(trace_data);
ylim(ax, [ymin, ymax]);
if ~isempty(t_axis)
    xlim(ax, [min(t_axis), max(t_axis)]);
end

if isstruct(stim_windows) && isfield(stim_windows, channel_name) ...
        && isfield(stim_windows.(channel_name), 'stim_frames') ...
        && isfield(stim_windows, 'orientations')
    channel_windows = stim_windows.(channel_name);
    stim_frames = round(double(channel_windows.stim_frames));
    orientations = mod(double(stim_windows.orientations(:)), 360);
    frame_rate = resolve_trace_frame_rate_from_time_rec(t_axis);
    if size(stim_frames, 1) == numel(orientations) && isfinite(frame_rate) && frame_rate > 0
        angle_colors = hsv(max(1, numel(unique(orientations))));
        unique_angles = unique(orientations, 'stable');
        for stim_idx = 1:size(stim_frames, 1)
            if any(~isfinite(stim_frames(stim_idx, :)))
                continue;
            end
            start_time = double(stim_frames(stim_idx, 1)) / frame_rate;
            stop_time = double(stim_frames(stim_idx, 2)) / frame_rate;
            if stop_time < start_time
                continue;
            end
            color_idx = find(abs(unique_angles - orientations(stim_idx)) < 1e-9, 1, 'first');
            if isempty(color_idx)
                color_idx = 1;
            end
            shade_color = lighten_color(angle_colors(color_idx, :), 0.70);
            patch(ax, [start_time stop_time stop_time start_time], ...
                [ymin ymin ymax ymax], shade_color, ...
                'FaceAlpha', 0.18, 'EdgeColor', 'none');
            text(ax, mean([start_time, stop_time]), ymax, sprintf('%g deg', orientations(stim_idx)), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
                'FontSize', 7, 'Color', [0.12 0.12 0.12], 'Clipping', 'on');
        end
    end
end

plot(ax, t_axis, trace_data, 'Color', trace_color, 'LineWidth', 0.9);
title(ax, title_text, 'Interpreter', 'none');
ylabel(ax, y_axis_label);
set(ax, 'TickDir', 'out');
box(ax, 'off');
end

function frame_rate = resolve_trace_frame_rate_from_time_rec(t_axis)
frame_rate = NaN;
if numel(t_axis) < 2
    return;
end
dt = median(diff(double(t_axis(:))), 'omitnan');
if isfinite(dt) && dt > 0
    frame_rate = 1 / dt;
end
end

function safe_name = make_safe_filename_rec(name_in)
safe_name = regexprep(char(string(name_in)), '[<>:"/\\|?*]', '_');
safe_name = regexprep(safe_name, '\s+', '_');
if isempty(safe_name)
    safe_name = 'unnamed';
end
end

function comparison = build_random_grating_cycle_voltage_dsi_comparison_rec( ...
    cycle_entries, angles, nrois, subplot_root)
comparison_dir = string(fullfile(subplot_root, 'cycle_direction_preference_comparison'));
if ~isfolder(comparison_dir)
    mkdir(comparison_dir);
end

ncycles = numel(cycle_entries);
cycle_names = string({cycle_entries.cycle_name})';
voltage_dsi_by_roi_cycle = NaN(nrois, ncycles);
voltage_pref_dir_by_roi_cycle = NaN(nrois, ncycles);
source_by_cycle = strings(ncycles, 1);

for cycle_idx = 1:ncycles
    tuning = cycle_entries(cycle_idx).grating_tuning.voltage;
    if isfield(tuning, 'dsi') && numel(tuning.dsi) >= nrois
        voltage_dsi_by_roi_cycle(:, cycle_idx) = double(tuning.dsi(1:nrois));
        if isfield(tuning, 'pref_dir') && numel(tuning.pref_dir) >= nrois
            voltage_pref_dir_by_roi_cycle(:, cycle_idx) = double(tuning.pref_dir(1:nrois));
        end
        source_by_cycle(cycle_idx) = "stim_results.tuning.voltage.dsi";
    else
        response = NaN(nrois, numel(angles));
        if isfield(tuning, 'unique_orientations') && isfield(tuning, 'response_by_condition')
            cycle_angles = mod(double(tuning.unique_orientations(:)'), 360);
            for angle_idx = 1:numel(cycle_angles)
                target_idx = find(abs(angles - cycle_angles(angle_idx)) < 1e-9, 1, 'first');
                if ~isempty(target_idx)
                    response(:, target_idx) = double(tuning.response_by_condition(1:nrois, angle_idx));
                end
            end
            metrics = compute_grating_metrics_from_curve_rec(response, angles);
            voltage_dsi_by_roi_cycle(:, cycle_idx) = metrics.dsi;
            voltage_pref_dir_by_roi_cycle(:, cycle_idx) = metrics.pref_dir;
            source_by_cycle(cycle_idx) = "recomputed_from_voltage_response_by_condition";
        else
            source_by_cycle(cycle_idx) = "missing_voltage_tuning";
        end
    end
end

mean_dsi = mean(voltage_dsi_by_roi_cycle, 1, 'omitnan')';
median_dsi = median(voltage_dsi_by_roi_cycle, 1, 'omitnan')';
sem_dsi = std(voltage_dsi_by_roi_cycle, 0, 1, 'omitnan')' ./ ...
    sqrt(max(1, sum(isfinite(voltage_dsi_by_roi_cycle), 1)'));
roi_count = sum(isfinite(voltage_dsi_by_roi_cycle), 1)';
max_dsi = max(voltage_dsi_by_roi_cycle, [], 1, 'omitnan')';
[strongest_cycle_mean_dsi, strongest_idx] = max(mean_dsi, [], 'omitnan');
if isempty(strongest_idx) || ~isfinite(strongest_cycle_mean_dsi)
    strongest_idx = NaN;
    strongest_cycle_name = "";
else
    strongest_cycle_name = cycle_names(strongest_idx);
end

summary_table = table(cycle_names, mean_dsi, median_dsi, sem_dsi, roi_count, max_dsi, source_by_cycle, ...
    'VariableNames', {'cycle_name', 'mean_voltage_dsi', 'median_voltage_dsi', ...
    'sem_voltage_dsi', 'roi_count', 'max_voltage_dsi', 'dsi_source'});
summary_csv = string(fullfile(comparison_dir, 'cyc_v_dsi.csv'));
writetable(summary_table, summary_csv);
[by_roi_dir, by_roi_fig, by_roi_png, roi_summary_table, ...
    strongest_cycle_by_roi, strongest_dsi_by_roi] = ...
    plot_voltage_dsi_by_roi_cycle_comparison_rec( ...
    comparison_dir, cycle_names, voltage_dsi_by_roi_cycle, voltage_pref_dir_by_roi_cycle);
roi_summary_csv = string(fullfile(comparison_dir, 'roi_cyc_v_dsi.csv'));
writetable(roi_summary_table, roi_summary_csv);

summary_mat = string(fullfile(comparison_dir, 'cyc_v_dsi.mat'));
comparison = struct( ...
    'status', "completed", ...
    'output_dir', comparison_dir, ...
    'summary_csv', summary_csv, ...
    'summary_mat', summary_mat, ...
    'summary_fig', string(fullfile(comparison_dir, 'cyc_v_dsi.fig')), ...
    'summary_png', string(fullfile(comparison_dir, 'cyc_v_dsi.png')), ...
    'roi_summary_csv', roi_summary_csv, ...
    'by_roi_dir', by_roi_dir, ...
    'by_roi_fig', by_roi_fig, ...
    'by_roi_png', by_roi_png, ...
    'strongest_cycle_name', strongest_cycle_name, ...
    'strongest_cycle_mean_dsi', strongest_cycle_mean_dsi, ...
    'strongest_cycle_by_roi', strongest_cycle_by_roi, ...
    'strongest_dsi_by_roi', strongest_dsi_by_roi, ...
    'strongest_rule', "cycle with largest mean voltage DSI across ROIs", ...
    'strongest_roi_rule', "within each ROI, cycle with largest voltage DSI", ...
    'cycle_names', cycle_names, ...
    'angles', angles, ...
    'voltage_dsi_by_roi_cycle', voltage_dsi_by_roi_cycle, ...
    'voltage_pref_dir_by_roi_cycle', voltage_pref_dir_by_roi_cycle, ...
    'summary_table', summary_table);

if ~skip_existing_figure_bundle_rec(comparison.summary_fig, comparison.summary_png, ...
        'random-grating cycle voltage DSI summary')
    fig = figure('Color', 'w', 'Name', 'Random Grating Cycle Voltage DSI Comparison', ...
        'Position', [100, 80, 1500, 850]);
    tiledlayout(fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    ax_bar = nexttile;
    bar(ax_bar, mean_dsi, 'FaceColor', [0.85 0.15 0.15], 'EdgeColor', 'none');
    hold(ax_bar, 'on');
    errorbar(ax_bar, 1:ncycles, mean_dsi, sem_dsi, 'k.', 'LineWidth', 1.2);
    for cycle_idx = 1:ncycles
        roi_values = voltage_dsi_by_roi_cycle(:, cycle_idx);
        roi_values = roi_values(isfinite(roi_values));
        if isempty(roi_values)
            continue;
        end
        scatter(ax_bar, cycle_idx + 0.10 * (rand(size(roi_values)) - 0.5), roi_values, ...
            20, [0.15 0.15 0.15], 'filled', 'MarkerFaceAlpha', 0.35);
    end
    if isfinite(strongest_idx)
        plot(ax_bar, strongest_idx, strongest_cycle_mean_dsi, 'p', ...
            'MarkerSize', 16, 'MarkerFaceColor', [1.0 0.82 0.10], ...
            'MarkerEdgeColor', [0.2 0.2 0.2]);
    end
    set(ax_bar, 'XTick', 1:ncycles, 'XTickLabel', cycle_names, ...
        'TickLabelInterpreter', 'none', 'TickDir', 'out');
    xtickangle(ax_bar, 35);
    ylabel(ax_bar, 'Voltage DSI');
    title(ax_bar, sprintf('Mean voltage DSI by Cycle | strongest: %s', strongest_cycle_name), ...
        'Interpreter', 'none');
    box(ax_bar, 'off');

    ax_heat = nexttile;
    imagesc(ax_heat, 1:ncycles, 1:nrois, voltage_dsi_by_roi_cycle);
    set(ax_heat, 'YDir', 'normal', 'XTick', 1:ncycles, 'XTickLabel', cycle_names, ...
        'TickLabelInterpreter', 'none', 'TickDir', 'out');
    xtickangle(ax_heat, 35);
    xlabel(ax_heat, 'Cycle');
    ylabel(ax_heat, 'ROI');
    title(ax_heat, 'Voltage DSI per ROI and Cycle');
    cb = colorbar(ax_heat);
    ylabel(cb, 'Voltage DSI');
    finite_dsi = voltage_dsi_by_roi_cycle(isfinite(voltage_dsi_by_roi_cycle));
    if ~isempty(finite_dsi)
        clim(ax_heat, [min(0, min(finite_dsi)), max(1, max(finite_dsi))]);
    end
    colormap(ax_heat, turbo(256));
    box(ax_heat, 'off');

    save_figure_bundle_preserve_layout_rec(fig, comparison.summary_fig, comparison.summary_png);
    close(fig);
end
summary_mat_payload = struct( ...
    'comparison', comparison, ...
    'summary_table', summary_table, ...
    'roi_summary_table', roi_summary_table, ...
    'voltage_dsi_by_roi_cycle', voltage_dsi_by_roi_cycle, ...
    'voltage_pref_dir_by_roi_cycle', voltage_pref_dir_by_roi_cycle);
if all(arrayfun(@(entry) isfield(entry.grating_tuning, 'voltage_sensitivity') ...
        && has_complete_grating_curves_rec(entry.grating_tuning.voltage_sensitivity), cycle_entries))
    sensitivity_comparison = build_random_grating_cycle_voltage_dsi_comparison_for_field_rec( ...
        cycle_entries, angles, nrois, comparison_dir, ...
        'voltage_sensitivity', 'Voltage sensitivity-average DSI', ...
        'voltage_sensitivity', 'voltage_sensitivity');
    comparison.voltage_sensitivity = sensitivity_comparison;
    summary_mat_payload.comparison = comparison;
    summary_mat_payload.voltage_sensitivity_dsi_by_roi_cycle = ...
        sensitivity_comparison.voltage_dsi_by_roi_cycle;
    summary_mat_payload.voltage_sensitivity_pref_dir_by_roi_cycle = ...
        sensitivity_comparison.voltage_pref_dir_by_roi_cycle;
else
    fprintf('Skipping 11 voltage-sensitivity DSI comparison: one or more cycles lack voltage_sensitivity tuning.\n');
end
save_mat_file_resilient_rec(summary_mat, summary_mat_payload, '-v7.3');
end

function comparison = build_random_grating_cycle_voltage_dsi_comparison_for_field_rec( ...
    cycle_entries, angles, nrois, comparison_root, tuning_field, plot_label, file_suffix, source_label)
file_suffix_char = char(string(file_suffix));
plot_label_char = char(string(plot_label));
comparison_dir = string(fullfile(comparison_root, file_suffix_char));
if ~isfolder(comparison_dir)
    mkdir(comparison_dir);
end
ncycles = numel(cycle_entries);
cycle_names = string({cycle_entries.cycle_name})';
voltage_dsi_by_roi_cycle = NaN(nrois, ncycles);
voltage_pref_dir_by_roi_cycle = NaN(nrois, ncycles);
source_by_cycle = strings(ncycles, 1);
for cycle_idx = 1:ncycles
    tuning = cycle_entries(cycle_idx).grating_tuning.(tuning_field);
    if isfield(tuning, 'dsi') && numel(tuning.dsi) >= nrois
        voltage_dsi_by_roi_cycle(:, cycle_idx) = double(tuning.dsi(1:nrois));
        if isfield(tuning, 'pref_dir') && numel(tuning.pref_dir) >= nrois
            voltage_pref_dir_by_roi_cycle(:, cycle_idx) = double(tuning.pref_dir(1:nrois));
        end
        source_by_cycle(cycle_idx) = string(source_label);
    end
end
mean_dsi = mean(voltage_dsi_by_roi_cycle, 1, 'omitnan')';
median_dsi = median(voltage_dsi_by_roi_cycle, 1, 'omitnan')';
sem_dsi = std(voltage_dsi_by_roi_cycle, 0, 1, 'omitnan')' ./ ...
    sqrt(max(1, sum(isfinite(voltage_dsi_by_roi_cycle), 1)'));
roi_count = sum(isfinite(voltage_dsi_by_roi_cycle), 1)';
max_dsi = max(voltage_dsi_by_roi_cycle, [], 1, 'omitnan')';
[strongest_cycle_mean_dsi, strongest_idx] = max(mean_dsi, [], 'omitnan');
if isempty(strongest_idx) || ~isfinite(strongest_cycle_mean_dsi)
    strongest_cycle_name = "";
else
    strongest_cycle_name = cycle_names(strongest_idx);
end
summary_table = table(cycle_names, mean_dsi, median_dsi, sem_dsi, roi_count, max_dsi, source_by_cycle, ...
    'VariableNames', {'cycle_name', 'mean_voltage_dsi', 'median_voltage_dsi', ...
    'sem_voltage_dsi', 'roi_count', 'max_voltage_dsi', 'dsi_source'});
summary_csv = string(fullfile(comparison_dir, sprintf('cyc_%s_dsi.csv', file_suffix_char)));
writetable(summary_table, summary_csv);
[by_roi_dir, by_roi_fig, by_roi_png, roi_summary_table, strongest_cycle_by_roi, strongest_dsi_by_roi] = ...
    plot_voltage_dsi_by_roi_cycle_comparison_rec( ...
    comparison_dir, cycle_names, voltage_dsi_by_roi_cycle, voltage_pref_dir_by_roi_cycle, file_suffix, plot_label);
roi_summary_csv = string(fullfile(comparison_dir, sprintf('roi_cyc_%s_dsi.csv', file_suffix_char)));
writetable(roi_summary_table, roi_summary_csv);

summary_fig = string(fullfile(comparison_dir, sprintf('cyc_%s_dsi.fig', file_suffix_char)));
summary_png = string(fullfile(comparison_dir, sprintf('cyc_%s_dsi.png', file_suffix_char)));
if ~skip_existing_figure_bundle_rec(summary_fig, summary_png, ...
        sprintf('random-grating cycle %s summary', plot_label_char))
    fig = figure('Color', 'w', 'Name', sprintf('Random Grating Cycle %s Comparison', plot_label_char), ...
        'Position', [100, 80, 1500, 850]);
    tiledlayout(fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    ax_bar = nexttile;
    bar(ax_bar, mean_dsi, 'FaceColor', [0.85 0.15 0.15], 'EdgeColor', 'none');
    hold(ax_bar, 'on');
    errorbar(ax_bar, 1:ncycles, mean_dsi, sem_dsi, 'k.', 'LineWidth', 1.2);
    set(ax_bar, 'XTick', 1:ncycles, 'XTickLabel', cycle_names, ...
        'TickLabelInterpreter', 'none', 'TickDir', 'out');
    xtickangle(ax_bar, 35);
    ylabel(ax_bar, plot_label_char);
    title(ax_bar, sprintf('Mean %s by Cycle | strongest: %s', plot_label_char, char(strongest_cycle_name)), ...
        'Interpreter', 'none');
    grid(ax_bar, 'off');
    box(ax_bar, 'off');
    ax_heat = nexttile;
    imagesc(ax_heat, 1:ncycles, 1:nrois, voltage_dsi_by_roi_cycle);
    set(ax_heat, 'YDir', 'normal', 'XTick', 1:ncycles, 'XTickLabel', cycle_names, ...
        'TickLabelInterpreter', 'none', 'TickDir', 'out');
    xtickangle(ax_heat, 35);
    xlabel(ax_heat, 'Cycle');
    ylabel(ax_heat, 'ROI');
    title(ax_heat, sprintf('%s per ROI and Cycle', plot_label_char), 'Interpreter', 'none');
    cb = colorbar(ax_heat);
    ylabel(cb, plot_label_char);
    finite_dsi = voltage_dsi_by_roi_cycle(isfinite(voltage_dsi_by_roi_cycle));
    if ~isempty(finite_dsi)
        clim(ax_heat, [min(0, min(finite_dsi)), max(1, max(finite_dsi))]);
    end
    colormap(ax_heat, turbo(256));
    grid(ax_heat, 'off');
    box(ax_heat, 'off');
    save_figure_bundle_preserve_layout_rec(fig, summary_fig, summary_png);
    close(fig);
end

summary_mat = string(fullfile(comparison_dir, sprintf('cyc_%s_dsi.mat', file_suffix_char)));
comparison = struct( ...
    'status', "completed", ...
    'output_dir', comparison_dir, ...
    'summary_csv', summary_csv, ...
    'summary_mat', summary_mat, ...
    'summary_fig', summary_fig, ...
    'summary_png', summary_png, ...
    'roi_summary_csv', roi_summary_csv, ...
    'by_roi_dir', by_roi_dir, ...
    'by_roi_fig', by_roi_fig, ...
    'by_roi_png', by_roi_png, ...
    'strongest_cycle_name', strongest_cycle_name, ...
    'strongest_cycle_mean_dsi', strongest_cycle_mean_dsi, ...
    'strongest_cycle_by_roi', strongest_cycle_by_roi, ...
    'strongest_dsi_by_roi', strongest_dsi_by_roi, ...
    'cycle_names', cycle_names, ...
    'angles', angles, ...
    'voltage_dsi_by_roi_cycle', voltage_dsi_by_roi_cycle, ...
    'voltage_pref_dir_by_roi_cycle', voltage_pref_dir_by_roi_cycle, ...
    'summary_table', summary_table, ...
    'roi_summary_table', roi_summary_table, ...
    'source_label', string(source_label));
payload = struct('comparison', comparison, 'summary_table', summary_table, ...
    'roi_summary_table', roi_summary_table, ...
    'voltage_dsi_by_roi_cycle', voltage_dsi_by_roi_cycle, ...
    'voltage_pref_dir_by_roi_cycle', voltage_pref_dir_by_roi_cycle);
save_mat_file_resilient_rec(summary_mat, payload, '-v7.3');
end

function [roi_dir, fig_files, png_files, roi_summary_table, strongest_cycle_by_roi, strongest_dsi_by_roi] = ...
    plot_voltage_dsi_by_roi_cycle_comparison_rec( ...
    comparison_dir, cycle_names, voltage_dsi_by_roi_cycle, voltage_pref_dir_by_roi_cycle, file_suffix, plot_label)
if nargin < 5 || isempty(file_suffix)
    file_suffix = "voltage";
end
if nargin < 6 || isempty(plot_label)
    plot_label = "Voltage DSI";
end
file_suffix_char = char(string(file_suffix));
plot_label_char = char(string(plot_label));
roi_dir = string(fullfile(comparison_dir, 'roi'));
if ~isfolder(roi_dir)
    mkdir(roi_dir);
end

[nrois, ncycles] = size(voltage_dsi_by_roi_cycle);
fig_files = strings(nrois, 1);
png_files = strings(nrois, 1);
strongest_cycle_by_roi = strings(nrois, 1);
strongest_dsi_by_roi = NaN(nrois, 1);

roi_column = repelem((1:nrois)', ncycles);
cycle_column = repmat(cycle_names(:), nrois, 1);
voltage_dsi_column = reshape(voltage_dsi_by_roi_cycle.', [], 1);
voltage_pref_dir_column = reshape(voltage_pref_dir_by_roi_cycle.', [], 1);
roi_summary_table = table(roi_column, cycle_column, voltage_dsi_column, voltage_pref_dir_column, ...
    'VariableNames', {'roi', 'cycle_name', 'voltage_dsi', 'voltage_pref_dir'});

for roi_idx = 1:nrois
    roi_dsi = voltage_dsi_by_roi_cycle(roi_idx, :);
    roi_pref_dir = voltage_pref_dir_by_roi_cycle(roi_idx, :);
    valid_idx = find(isfinite(roi_dsi));
    if isempty(valid_idx)
        strongest_idx = NaN;
        strongest_cycle_name = "";
        strongest_dsi = NaN;
    else
        [strongest_dsi, local_idx] = max(roi_dsi(valid_idx));
        strongest_idx = valid_idx(local_idx);
        strongest_cycle_name = cycle_names(strongest_idx);
    end
    strongest_cycle_by_roi(roi_idx) = strongest_cycle_name;
    strongest_dsi_by_roi(roi_idx) = strongest_dsi;

    fig_files(roi_idx) = string(fullfile(roi_dir, sprintf( ...
        'r%03d_%s_dsi.fig', roi_idx, file_suffix_char)));
    png_files(roi_idx) = string(fullfile(roi_dir, sprintf( ...
        'r%03d_%s_dsi.png', roi_idx, file_suffix_char)));
    if skip_existing_figure_bundle_rec(fig_files(roi_idx), png_files(roi_idx), ...
            sprintf('ROI %03d %s cycle comparison', roi_idx, plot_label_char))
        continue;
    end

    fig = figure('Color', 'w', ...
        'Name', sprintf('ROI %03d Voltage DSI Cycle Comparison', roi_idx), ...
        'Position', [120, 100, 1100, 650]);
    ax = axes(fig);
    hold(ax, 'on');
    plot(ax, 1:ncycles, roi_dsi, '-o', ...
        'Color', [0.85 0.12 0.10], ...
        'MarkerFaceColor', [0.85 0.12 0.10], ...
        'LineWidth', 1.6);
    if isfinite(strongest_idx)
        plot(ax, strongest_idx, strongest_dsi, 'p', ...
            'MarkerSize', 16, 'MarkerFaceColor', [1.0 0.82 0.10], ...
            'MarkerEdgeColor', [0.20 0.20 0.20]);
    end
    for cycle_idx = 1:ncycles
        if isfinite(roi_dsi(cycle_idx)) && isfinite(roi_pref_dir(cycle_idx))
            text(ax, cycle_idx, roi_dsi(cycle_idx), ...
                sprintf('  %g deg', roi_pref_dir(cycle_idx)), ...
                'VerticalAlignment', 'bottom', ...
                'FontSize', 8, ...
                'Color', [0.15 0.15 0.15], ...
                'Clipping', 'on');
        end
    end
    set(ax, 'XTick', 1:ncycles, 'XTickLabel', cycle_names, ...
        'TickLabelInterpreter', 'none', 'TickDir', 'out');
    xtickangle(ax, 35);
    xlim(ax, [0.75, ncycles + 0.25]);
    finite_dsi = roi_dsi(isfinite(roi_dsi));
    if isempty(finite_dsi)
        ylim(ax, [0, 1]);
    else
        y_min = min(0, min(finite_dsi));
        y_max = max(1, max(finite_dsi));
        if y_max <= y_min
            y_max = y_min + 1;
        end
        ylim(ax, [y_min, y_max + 0.08 * (y_max - y_min)]);
    end
    ylabel(ax, plot_label_char);
    xlabel(ax, 'Cycle');
    title(ax, sprintf('ROI %03d | strongest %s | DSI=%.4g', ...
        roi_idx, char(strongest_cycle_name), strongest_dsi), ...
        'Interpreter', 'none');
    grid(ax, 'off');
    box(ax, 'off');

    save_figure_bundle_preserve_layout_rec(fig, fig_files(roi_idx), png_files(roi_idx));
    close(fig);
end
end

function cycle_tuning = plot_random_grating_cycle_tuning_curves_rec( ...
    cycle_entries, angles, nrois, subplot_root)
output_dir = string(fullfile(subplot_root, 'cyc'));
if ~isfolder(output_dir)
    mkdir(output_dir);
end

ncycles = numel(cycle_entries);
cycle_names = string({cycle_entries.cycle_name})';
fig_files = strings(nrois, ncycles);
png_files = strings(nrois, ncycles);
voltage_sensitivity_fig_files = strings(nrois, ncycles);
voltage_sensitivity_png_files = strings(nrois, ncycles);
rows = cell(0, 1);
skipped_cycles = strings(0, 1);
skipped_reasons = strings(0, 1);

for roi_idx = 1:nrois
    roi_dir = fullfile(output_dir, sprintf('r%03d', roi_idx));
    if ~isfolder(roi_dir)
        mkdir(roi_dir);
    end

    for cycle_idx = 1:ncycles
        cycle_name = string(cycle_entries(cycle_idx).cycle_name);
        cycle_tag = sprintf('c%02d', cycle_idx);
        if ~isfield(cycle_entries(cycle_idx).grating_tuning, 'voltage') ...
                || ~isfield(cycle_entries(cycle_idx).grating_tuning, 'calcium')
            skipped_cycles(end+1, 1) = cycle_name; %#ok<AGROW>
            skipped_reasons(end+1, 1) = "missing voltage or calcium tuning"; %#ok<AGROW>
            warning('Dual_rec_analysis3:CycleTuningMissing', ...
                'Skipping cycle tuning curve for ROI %03d %s: missing voltage or calcium tuning.', ...
                roi_idx, cycle_name);
            continue;
        end

        voltage_tuning = cycle_entries(cycle_idx).grating_tuning.voltage;
        calcium_tuning = cycle_entries(cycle_idx).grating_tuning.calcium;
        if ~has_complete_grating_curves_rec(voltage_tuning) || ~has_complete_grating_curves_rec(calcium_tuning)
            skipped_cycles(end+1, 1) = cycle_name; %#ok<AGROW>
            skipped_reasons(end+1, 1) = "incomplete tuning curves"; %#ok<AGROW>
            warning('Dual_rec_analysis3:CycleTuningIncomplete', ...
                'Skipping cycle tuning curve for ROI %03d %s: incomplete tuning curves.', ...
                roi_idx, cycle_name);
            continue;
        end

        [voltage_response, voltage_nonstim] = align_cycle_tuning_response_rec(voltage_tuning, angles, nrois);
        [calcium_response, calcium_nonstim] = align_cycle_tuning_response_rec(calcium_tuning, angles, nrois);

        fig_files(roi_idx, cycle_idx) = string(fullfile(roi_dir, sprintf( ...
            'r%03d_%s_tuning.fig', roi_idx, cycle_tag)));
        png_files(roi_idx, cycle_idx) = string(fullfile(roi_dir, sprintf( ...
            'r%03d_%s_tuning.png', roi_idx, cycle_tag)));
        if ~skip_existing_figure_bundle_rec(fig_files(roi_idx, cycle_idx), png_files(roi_idx, cycle_idx), ...
                sprintf('ROI %03d %s random-grating tuning curve', roi_idx, char(cycle_name)))
            fig = figure('Color', 'w', ...
                'Name', sprintf('ROI %03d %s Random Grating Tuning Curves', roi_idx, char(cycle_name)), ...
                'Position', [100, 90, 1500, 760]);
            tiledlayout(fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

            plot_roi_cycle_tuning_channel_rec(nexttile, angles, voltage_response(roi_idx, :), voltage_nonstim(roi_idx, :), ...
                [0.85 0.12 0.10], sprintf('ROI %03d | %s | Voltage tuning', roi_idx, char(cycle_name)), ...
                'Voltage response');
            plot_roi_cycle_tuning_channel_rec(nexttile, angles, calcium_response(roi_idx, :), calcium_nonstim(roi_idx, :), ...
                [0.10 0.55 0.18], sprintf('ROI %03d | %s | Calcium tuning', roi_idx, char(cycle_name)), ...
                'Calcium response');
            sgtitle(fig, sprintf('Random drifting grating tuning | ROI %03d | %s', ...
                roi_idx, char(cycle_name)), 'Interpreter', 'none');
            save_figure_bundle_preserve_layout_rec(fig, fig_files(roi_idx, cycle_idx), png_files(roi_idx, cycle_idx));
            close(fig);
        end

        rows{end+1, 1} = build_roi_cycle_tuning_summary_rows_rec( ...
            roi_idx, cycle_name, "voltage", angles, voltage_response(roi_idx, :), voltage_nonstim(roi_idx, :)); %#ok<AGROW>
        rows{end+1, 1} = build_roi_cycle_tuning_summary_rows_rec( ...
            roi_idx, cycle_name, "calcium", angles, calcium_response(roi_idx, :), calcium_nonstim(roi_idx, :)); %#ok<AGROW>

        if isfield(cycle_entries(cycle_idx).grating_tuning, 'voltage_sensitivity') ...
                && has_complete_grating_curves_rec(cycle_entries(cycle_idx).grating_tuning.voltage_sensitivity)
            voltage_sensitivity_tuning = cycle_entries(cycle_idx).grating_tuning.voltage_sensitivity;
            [voltage_sensitivity_response, voltage_sensitivity_nonstim] = ...
                align_cycle_tuning_response_rec(voltage_sensitivity_tuning, angles, nrois);
            voltage_sensitivity_fig_files(roi_idx, cycle_idx) = string(fullfile(roi_dir, sprintf( ...
                'r%03d_%s_v_sens_tuning.fig', roi_idx, cycle_tag)));
            voltage_sensitivity_png_files(roi_idx, cycle_idx) = string(fullfile(roi_dir, sprintf( ...
                'r%03d_%s_v_sens_tuning.png', roi_idx, cycle_tag)));
            if ~skip_existing_figure_bundle_rec( ...
                    voltage_sensitivity_fig_files(roi_idx, cycle_idx), ...
                    voltage_sensitivity_png_files(roi_idx, cycle_idx), ...
                    sprintf('ROI %03d %s random-grating voltage-sensitivity tuning curve', roi_idx, char(cycle_name)))
                fig = figure('Color', 'w', ...
                    'Name', sprintf('ROI %03d %s Random Grating Voltage Sensitivity Tuning Curves', roi_idx, char(cycle_name)), ...
                    'Position', [100, 90, 1500, 760]);
                tiledlayout(fig, 1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
                plot_roi_cycle_tuning_channel_rec(nexttile, angles, ...
                    voltage_sensitivity_response(roi_idx, :), voltage_sensitivity_nonstim(roi_idx, :), ...
                    [0.85 0.12 0.10], sprintf('ROI %03d | %s | Voltage sensitivity tuning', roi_idx, char(cycle_name)), ...
                    'Voltage mean sensitivity');
                plot_roi_cycle_tuning_channel_rec(nexttile, angles, calcium_response(roi_idx, :), calcium_nonstim(roi_idx, :), ...
                    [0.10 0.55 0.18], sprintf('ROI %03d | %s | Calcium tuning', roi_idx, char(cycle_name)), ...
                    'Calcium response');
                sgtitle(fig, sprintf('Random drifting grating tuning | ROI %03d | %s | voltage sensitivity average', ...
                    roi_idx, char(cycle_name)), 'Interpreter', 'none');
                save_figure_bundle_preserve_layout_rec(fig, ...
                    voltage_sensitivity_fig_files(roi_idx, cycle_idx), ...
                    voltage_sensitivity_png_files(roi_idx, cycle_idx));
                close(fig);
            end
            rows{end+1, 1} = build_roi_cycle_tuning_summary_rows_rec( ...
                roi_idx, cycle_name, "voltage_sensitivity", angles, ...
                voltage_sensitivity_response(roi_idx, :), voltage_sensitivity_nonstim(roi_idx, :)); %#ok<AGROW>
        end
    end
end

if isempty(rows)
    summary_table = table();
else
    summary_table = vertcat(rows{:});
end
summary_csv = string(fullfile(output_dir, 'roi_cyc_tuning.csv'));
writetable(summary_table, summary_csv);
summary_mat = string(fullfile(output_dir, 'roi_cyc_tuning.mat'));

cycle_tuning = struct( ...
    'status', "completed", ...
    'output_dir', output_dir, ...
    'fig_files', fig_files, ...
    'png_files', png_files, ...
    'voltage_sensitivity_fig_files', voltage_sensitivity_fig_files, ...
    'voltage_sensitivity_png_files', voltage_sensitivity_png_files, ...
    'summary_csv', summary_csv, ...
    'summary_mat', summary_mat, ...
    'cycle_names', cycle_names, ...
    'angles', angles, ...
    'skipped_cycles', skipped_cycles, ...
    'skipped_reasons', skipped_reasons, ...
    'display_rule', "per ROI per Cycle tuning curve from saved stim_results.tuning; stim and nonstim plotted for the same ROI");

payload = struct( ...
    'cycle_tuning', cycle_tuning, ...
    'summary_table', summary_table);
save_mat_file_resilient_rec(summary_mat, payload, '-v7.3');
end

function [response_aligned, nonstim_aligned] = align_cycle_tuning_response_rec(tuning, angles, nrois)
response_aligned = NaN(nrois, numel(angles));
nonstim_aligned = NaN(nrois, numel(angles));
cycle_angles = mod(double(tuning.unique_orientations(:)'), 360);
for angle_idx = 1:numel(cycle_angles)
    target_idx = find(abs(angles - cycle_angles(angle_idx)) < 1e-9, 1, 'first');
    if isempty(target_idx)
        continue;
    end
    response_aligned(:, target_idx) = double(tuning.response_by_condition(1:nrois, angle_idx));
    nonstim_aligned(:, target_idx) = double(tuning.nonstim_response_by_condition(1:nrois, angle_idx));
end
end

function plot_roi_cycle_tuning_channel_rec(ax, angles, response, nonstim_response, line_color, title_text, y_label)
hold(ax, 'on');
plot(ax, angles, nonstim_response, '--o', ...
    'Color', lighten_color(line_color, 0.45), 'MarkerFaceColor', 'w', ...
    'LineWidth', 1.1);
plot(ax, angles, response, '-o', ...
    'Color', line_color, 'MarkerFaceColor', line_color, 'LineWidth', 1.7);
xlabel(ax, 'Drifting direction (deg)');
ylabel(ax, y_label);
title(ax, title_text, 'Interpreter', 'none');
legend(ax, {'Nonstim response', 'Stim response'}, 'Location', 'best');
set(ax, 'TickDir', 'out');
grid(ax, 'off');
box(ax, 'off');
end

function rows = build_roi_cycle_tuning_summary_rows_rec(roi_idx, cycle_name, channel_name, angles, response, nonstim_response)
rows = table( ...
    repmat(roi_idx, numel(angles), 1), ...
    repmat(string(cycle_name), numel(angles), 1), ...
    repmat(string(channel_name), numel(angles), 1), ...
    double(angles(:)), ...
    double(response(:)), ...
    double(nonstim_response(:)), ...
    'VariableNames', {'roi', 'cycle_name', 'channel', 'angle_deg', ...
    'response', 'nonstim_response'});
end

function plot_random_grating_trial_heatmap_channel_rec(ax, trials, color_limit, title_text)
if isempty(trials)
    text(ax, 0.5, 0.5, 'No valid trials', ...
        'Units', 'normalized', 'HorizontalAlignment', 'center');
    axis(ax, 'off');
    return;
end

[heatmap_data, time_axis, labels] = build_random_grating_trial_heatmap_matrix_rec(trials);
image_handle = imagesc(ax, time_axis, 1:size(heatmap_data, 1), heatmap_data);
set(image_handle, 'AlphaData', isfinite(heatmap_data));
set(ax, 'Color', [0.88 0.88 0.88], 'YDir', 'normal', ...
    'YTick', 1:numel(labels), 'YTickLabel', labels, ...
    'TickLabelInterpreter', 'none', 'TickDir', 'out');
clim(ax, color_limit);
hold(ax, 'on');
xline(ax, 0, '--k', 'Stim onset', 'LineWidth', 1);
xlim(ax, [min(time_axis), max(time_axis)]);
xlabel(ax, 'Time from stimulus onset (s)');
ylabel(ax, 'Trials');
title(ax, title_text);
box(ax, 'off');
end

function [heatmap_data, time_axis, labels] = build_random_grating_trial_heatmap_matrix_rec(trials)
frame_rates = [trials.frame_rate];
reference_rate = frame_rates(1);
if any(abs(frame_rates - reference_rate) > max(1e-9, abs(reference_rate) * 1e-9))
    error('Random-grating heatmap trials have inconsistent frame rates within one channel.');
end

relative_frame_cells = arrayfun(@(trial) ...
    round(double(trial.time(:)) * reference_rate), trials, 'UniformOutput', false);
frame_min = min(cellfun(@min, relative_frame_cells));
frame_max = max(cellfun(@max, relative_frame_cells));
relative_frames = frame_min:frame_max;
heatmap_data = NaN(numel(trials), numel(relative_frames));
labels = strings(numel(trials), 1);
for trial_idx = 1:numel(trials)
    relative_trial_frames = relative_frame_cells{trial_idx};
    column_index = relative_trial_frames - frame_min + 1;
    heatmap_data(trial_idx, column_index) = double(trials(trial_idx).trace(:));
    labels(trial_idx) = sprintf('%s | T%d', ...
        trials(trial_idx).cycle_name, trials(trial_idx).trial_index);
end
time_axis = double(relative_frames) / reference_rate;
end

function text_out = upper_first_rec(text_in)
text_out = char(string(text_in));
if ~isempty(text_out)
    text_out(1) = upper(text_out(1));
end
end

function trials = collect_random_grating_roi_trials_rec( ...
    cycle_entries, channel_name, roi_idx, sorted_angles, polarity)
trial_template = struct( ...
    'trace', [], ...
    'time', [], ...
    'frame_rate', NaN, ...
    'angle', NaN, ...
    'cycle_name', "", ...
    'trial_index', NaN);
trials = repmat(trial_template, 0, 1);

for angle_idx = 1:numel(sorted_angles)
    target_angle = sorted_angles(angle_idx);
    for cycle_idx = 1:numel(cycle_entries)
        entry = cycle_entries(cycle_idx);
        windows = entry.stim_windows;
        if ~isfield(windows, channel_name) || ~isfield(windows.(channel_name), 'baseline_frames') ...
                || ~isfield(windows.(channel_name), 'stim_frames')
            error('Random-grating frame windows are missing for %s in cycle %s.', ...
                channel_name, entry.cycle_name);
        end
        orientations = mod(double(windows.orientations(:)), 360);
        baseline_frames = round(double(windows.(channel_name).baseline_frames));
        stim_frames = round(double(windows.(channel_name).stim_frames));
        if numel(orientations) ~= size(baseline_frames, 1) ...
                || numel(orientations) ~= size(stim_frames, 1)
            error('Random-grating trial labels and frame windows differ in cycle %s.', ...
                entry.cycle_name);
        end

        stage = entry.(channel_name).sensitivity;
        matching_trials = find(abs(orientations - target_angle) < 1e-9);
        for trial_idx = matching_trials(:)'
            frame_start = baseline_frames(trial_idx, 1);
            stim_onset = stim_frames(trial_idx, 1);
            frame_stop = stim_frames(trial_idx, 2);
            if frame_start < 1 || stim_onset < frame_start || frame_stop > size(stage.data, 1)
                error(['Random-grating trial %d in cycle %s has invalid %s frames [%d %d %d] ' ...
                    'for a %d-frame trace.'], trial_idx, entry.cycle_name, channel_name, ...
                    frame_start, stim_onset, frame_stop, size(stage.data, 1));
            end
            frame_index = (frame_start:frame_stop)';
            new_trial = trial_template;
            new_trial.trace = double(polarity) * double(stage.data(frame_index, roi_idx));
            new_trial.time = double(frame_index - stim_onset) / stage.frame_rate;
            new_trial.frame_rate = stage.frame_rate;
            new_trial.angle = target_angle;
            new_trial.cycle_name = entry.cycle_name;
            new_trial.trial_index = trial_idx;
            trials(end+1, 1) = new_trial; %#ok<AGROW>
        end
    end
end
end

function plot_random_grating_trial_stack_channel_rec(ax, trials, angles, title_text)
if isempty(trials)
    text(ax, 0.5, 0.5, 'No valid trials', ...
        'Units', 'normalized', 'HorizontalAlignment', 'center');
    axis(ax, 'off');
    return;
end

hold(ax, 'on');
spacing = compute_trial_cell_stack_spacing_rec({trials.trace});
ntrials = numel(trials);
offsets = (ntrials - (1:ntrials)) * spacing;
angle_colors = hsv(max(1, numel(angles)));
labels = strings(ntrials, 1);
single_direction_tile = numel(unique([trials.angle])) == 1;
for trial_idx = 1:ntrials
    color_idx = find(abs(angles - trials(trial_idx).angle) < 1e-9, 1, 'first');
    plot(ax, trials(trial_idx).time, trials(trial_idx).trace + offsets(trial_idx), ...
        'Color', angle_colors(color_idx, :), 'LineWidth', 0.8);
    if single_direction_tile
        labels(trial_idx) = sprintf('%s | T%d', ...
            trials(trial_idx).cycle_name, trials(trial_idx).trial_index);
    else
        labels(trial_idx) = sprintf('%g deg | %s | T%d', ...
            trials(trial_idx).angle, trials(trial_idx).cycle_name, trials(trial_idx).trial_index);
    end
end
xline(ax, 0, '--k', 'Stim onset', 'LineWidth', 1);

trial_angles = [trials.angle];
group_end = find(diff(trial_angles) ~= 0);
for boundary_idx = group_end(:)'
    separator_y = mean(offsets(boundary_idx:boundary_idx + 1));
    yline(ax, separator_y, ':', 'Color', [0.65 0.65 0.65], 'LineWidth', 0.8);
end

all_time = vertcat(trials.time);
xlim(ax, [min(all_time), max(all_time)]);
ylim(ax, [min(vertcat(trials.trace)) - 0.5 * spacing, ...
    max(vertcat(trials.trace)) + offsets(1) + 0.5 * spacing]);
[tick_values, tick_order] = sort(offsets, 'ascend');
set(ax, 'YTick', tick_values, 'YTickLabel', labels(tick_order), ...
    'TickLabelInterpreter', 'none', 'TickDir', 'out');
xlabel(ax, 'Time from stimulus onset (s)');
ylabel(ax, 'Trials grouped by direction');
title(ax, title_text);
grid(ax, 'off');
box(ax, 'off');
add_axis_scalebar_rec(ax, all_time, vertcat(trials.trace), [0.1 0.1 0.1], 'Sensitivity');
end

function spacing = compute_trial_cell_stack_spacing_rec(trace_cells)
trace_std = cellfun(@(trace) std(double(trace), 0, 'omitnan'), trace_cells);
spacing = median(trace_std, 'omitnan') * 6;
if ~isfinite(spacing) || spacing <= 0
    trace_ranges = cellfun(@(trace) ...
        max(double(trace), [], 'omitnan') - min(double(trace), [], 'omitnan'), trace_cells);
    spacing = max(1, median(trace_ranges, 'omitnan'));
end
end

function tf = has_complete_grating_curves_rec(tuning)
tf = isstruct(tuning) ...
    && isfield(tuning, 'unique_orientations') ...
    && isfield(tuning, 'response_by_condition') ...
    && isfield(tuning, 'nonstim_response_by_condition');
end

function channel = average_grating_channel_rec(cycle_entries, channel_name, orientations, nrois)
ncycles = numel(cycle_entries);
nangles = numel(orientations);
stim_per_cycle = NaN(nrois, nangles, ncycles);
nonstim_per_cycle = NaN(nrois, nangles, ncycles);
metrics_per_cycle = initialize_grating_metric_arrays_rec(nrois, ncycles);

for cycle_idx = 1:ncycles
    tuning = cycle_entries(cycle_idx).grating_tuning.(channel_name);
    cycle_angles = mod(double(tuning.unique_orientations(:)'), 360);
    for angle_idx = 1:numel(cycle_angles)
        target_idx = find(abs(orientations - cycle_angles(angle_idx)) < 1e-9, 1, 'first');
        if isempty(target_idx)
            continue;
        end
        stim_per_cycle(:, target_idx, cycle_idx) = tuning.response_by_condition(:, angle_idx);
        nonstim_per_cycle(:, target_idx, cycle_idx) = tuning.nonstim_response_by_condition(:, angle_idx);
    end
    metric_names = fieldnames(metrics_per_cycle);
    for metric_idx = 1:numel(metric_names)
        metric_name = metric_names{metric_idx};
        if isfield(tuning, metric_name)
            metrics_per_cycle.(metric_name)(:, cycle_idx) = tuning.(metric_name)(:);
        end
    end
end

stim_mean = mean(stim_per_cycle, 3, 'omitnan');
nonstim_mean = mean(nonstim_per_cycle, 3, 'omitnan');
stim_n = sum(isfinite(stim_per_cycle), 3);
nonstim_n = sum(isfinite(nonstim_per_cycle), 3);
stim_sem = std(stim_per_cycle, 0, 3, 'omitnan') ./ sqrt(max(1, stim_n));
nonstim_sem = std(nonstim_per_cycle, 0, 3, 'omitnan') ./ sqrt(max(1, nonstim_n));

metrics_from_mean_curve = compute_grating_metrics_from_curve_rec(stim_mean, orientations);

channel = struct( ...
    'stim_per_cycle', stim_per_cycle, ...
    'nonstim_per_cycle', nonstim_per_cycle, ...
    'stim_mean', stim_mean, ...
    'stim_sem', stim_sem, ...
    'stim_n', stim_n, ...
    'nonstim_mean', nonstim_mean, ...
    'nonstim_sem', nonstim_sem, ...
    'nonstim_n', nonstim_n, ...
    'metrics_per_cycle', metrics_per_cycle, ...
    'metrics_mean_across_cycles', summarize_grating_metrics_across_cycles_rec(metrics_per_cycle), ...
    'metrics_from_mean_curve', metrics_from_mean_curve);
end

function metrics = initialize_grating_metric_arrays_rec(nrois, ncycles)
metric_names = {'pref_dir', 'pref_ori', 'gdsi', 'gosi', 'dsi', 'osi'};
metrics = struct();
for idx = 1:numel(metric_names)
    metrics.(metric_names{idx}) = NaN(nrois, ncycles);
end
end

function summary = summarize_grating_metrics_across_cycles_rec(metrics)
metric_names = fieldnames(metrics);
summary = struct();
for idx = 1:numel(metric_names)
    values = metrics.(metric_names{idx});
    summary.(metric_names{idx}) = struct( ...
        'mean', mean(values, 2, 'omitnan'), ...
        'sem', std(values, 0, 2, 'omitnan') ./ sqrt(max(1, sum(isfinite(values), 2))), ...
        'n', sum(isfinite(values), 2));
end
end

function metrics = compute_grating_metrics_from_curve_rec(response_by_condition, orientations)
nrois = size(response_by_condition, 1);
metrics = struct( ...
    'pref_dir', NaN(nrois, 1), ...
    'pref_ori', NaN(nrois, 1), ...
    'gdsi', NaN(nrois, 1), ...
    'gosi', NaN(nrois, 1), ...
    'dsi', NaN(nrois, 1), ...
    'osi', NaN(nrois, 1));
theta_rad = deg2rad(orientations(:)');
for roi_idx = 1:nrois
    response = response_by_condition(roi_idx, :);
    if all(~isfinite(response)) || sum(response, 'omitnan') <= 0
        continue;
    end
    [~, pref_idx] = max(response, [], 'omitnan');
    metrics.pref_dir(roi_idx) = orientations(pref_idx);
    metrics.pref_ori(roi_idx) = mod(metrics.pref_dir(roi_idx), 180);
    response_sum = sum(response, 'omitnan');
    metrics.gdsi(roi_idx) = abs(sum(response .* exp(1i * theta_rad), 'omitnan') / response_sum);
    metrics.gosi(roi_idx) = abs(sum(response .* exp(1i * 2 * theta_rad), 'omitnan') / response_sum);

    opposite_idx = closest_angle_index_rec(orientations, mod(metrics.pref_dir(roi_idx) + 180, 360));
    orthogonal_idx = closest_angle_index_rec(orientations, mod(metrics.pref_dir(roi_idx) + 90, 360));
    pref_response = response(pref_idx);
    opposite_response = response(opposite_idx);
    orthogonal_response = response(orthogonal_idx);
    metrics.dsi(roi_idx) = (pref_response - opposite_response) / ...
        max(eps, pref_response + opposite_response);
    metrics.osi(roi_idx) = (pref_response - orthogonal_response) / ...
        max(eps, pref_response + orthogonal_response);
end
end

function tuning = build_grating_tuning_from_trial_response_rec(stim_response, nonstim_response, orientations)
stim_tuning = compute_grating_tuning_from_trial_response_rec(stim_response, orientations);
nonstim_tuning = compute_grating_tuning_from_trial_response_rec(nonstim_response, orientations);
tuning = stim_tuning;
tuning.nonstim_trial_response = nonstim_tuning.trial_response;
tuning.nonstim_response_by_condition = nonstim_tuning.response_by_condition;
tuning.nonstim_pref_dir = nonstim_tuning.pref_dir;
tuning.nonstim_pref_ori = nonstim_tuning.pref_ori;
end

function tuning = compute_grating_tuning_from_trial_response_rec(trial_response, orientations)
orientations = mod(double(orientations(:)), 360);
unique_orientations = unique(orientations, 'sorted')';
nrois = size(trial_response, 1);
response_by_condition = NaN(nrois, numel(unique_orientations));
for angle_idx = 1:numel(unique_orientations)
    mask = orientations == unique_orientations(angle_idx);
    response_by_condition(:, angle_idx) = mean(max(double(trial_response(:, mask)), 0), 2, 'omitnan');
end
metrics = compute_grating_metrics_from_curve_rec(response_by_condition, unique_orientations);
tuning = metrics;
tuning.unique_orientations = unique_orientations;
tuning.trial_response = double(trial_response);
tuning.response_by_condition = response_by_condition;
end

function idx = closest_angle_index_rec(angles, target_angle)
wrapped_difference = abs(mod(angles - target_angle + 180, 360) - 180);
[~, idx] = min(wrapped_difference);
end

function plot_record_grating_summary_rec(summary, fig_file, png_file)
if skip_existing_figure_bundle_rec(fig_file, png_file, 'record-average grating tuning summary')
    return;
end
fig = figure('Color', 'w', 'Name', 'Record Average Grating Tuning', ...
    'Position', [100, 80, 1600, 900]);
tiledlayout(fig, 2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
plot_population_grating_curves_rec(nexttile, summary.orientations, summary.voltage, ...
    [0.85 0.15 0.15], 'Voltage Spike Count', 'Spikes / trial');
plot_population_grating_curves_rec(nexttile, summary.orientations, summary.calcium, ...
    [0.10 0.65 0.20], 'Calcium Mean Sensitivity', 'Mean sensitivity');
plot_dsi_direction_radar_rec(nexttile, summary.orientations, summary.voltage, summary.calcium);
plot_osi_orientation_radar_rec(nexttile, summary.orientations, summary.voltage, summary.calcium);
plot_metric_scatter_rec(nexttile, summary.voltage, summary.calcium, 'dsi', 'DSI');
plot_metric_scatter_rec(nexttile, summary.voltage, summary.calcium, 'osi', 'OSI');
sgtitle(fig, sprintf('Record-Average Grating Tuning | %d included cycles', ...
    numel(summary.included_cycles)));
save_figure_bundle_rec(fig, fig_file, png_file);
close(fig);
end

function plot_record_grating_summary_voltage_sensitivity_rec(summary, fig_file, png_file)
if skip_existing_figure_bundle_rec(fig_file, png_file, ...
        'record-average grating tuning voltage-sensitivity summary')
    return;
end
summary_parallel = summary;
summary_parallel.voltage = summary.voltage_sensitivity;
fig = figure('Color', 'w', 'Name', 'Record Average Grating Tuning Voltage Sensitivity', ...
    'Position', [100, 80, 1600, 900]);
tiledlayout(fig, 2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');
plot_population_grating_curves_rec(nexttile, summary_parallel.orientations, summary_parallel.voltage, ...
    [0.85 0.15 0.15], 'Voltage Mean Sensitivity', 'Mean sensitivity');
plot_population_grating_curves_rec(nexttile, summary_parallel.orientations, summary_parallel.calcium, ...
    [0.10 0.65 0.20], 'Calcium Mean Sensitivity', 'Mean sensitivity');
plot_dsi_direction_radar_rec(nexttile, summary_parallel.orientations, summary_parallel.voltage, summary_parallel.calcium);
plot_osi_orientation_radar_rec(nexttile, summary_parallel.orientations, summary_parallel.voltage, summary_parallel.calcium);
plot_metric_scatter_rec(nexttile, summary_parallel.voltage, summary_parallel.calcium, 'dsi', 'DSI');
plot_metric_scatter_rec(nexttile, summary_parallel.voltage, summary_parallel.calcium, 'osi', 'OSI');
sgtitle(fig, sprintf('Record-Average Grating Tuning | Voltage sensitivity average | %d included cycles', ...
    numel(summary_parallel.included_cycles)));
save_figure_bundle_rec(fig, fig_file, png_file);
close(fig);
end

function plot_population_grating_curves_rec(ax, orientations, channel, line_color, title_text, y_label)
stim_population = mean(channel.stim_mean, 1, 'omitnan');
nonstim_population = mean(channel.nonstim_mean, 1, 'omitnan');
stim_sem = std(channel.stim_mean, 0, 1, 'omitnan') ./ ...
    sqrt(max(1, sum(isfinite(channel.stim_mean), 1)));
nonstim_sem = std(channel.nonstim_mean, 0, 1, 'omitnan') ./ ...
    sqrt(max(1, sum(isfinite(channel.nonstim_mean), 1)));
errorbar(ax, orientations, stim_population, stim_sem, '-o', ...
    'Color', line_color, 'MarkerFaceColor', line_color, 'LineWidth', 1.5);
hold(ax, 'on');
errorbar(ax, orientations, nonstim_population, nonstim_sem, '--o', ...
    'Color', lighten_color(line_color, 0.55), 'MarkerFaceColor', 'w', 'LineWidth', 1.2);
xlabel(ax, 'Direction (deg)');
ylabel(ax, y_label);
title(ax, title_text);
legend(ax, {'Grating', 'Non-grating'}, 'Location', 'best');
grid(ax, 'off');
end

function plot_metric_scatter_rec(ax, voltage, calcium, metric_name, metric_label)
voltage_cycle_mean = voltage.metrics_mean_across_cycles.(metric_name).mean;
calcium_cycle_mean = calcium.metrics_mean_across_cycles.(metric_name).mean;
voltage_curve_metric = voltage.metrics_from_mean_curve.(metric_name);
calcium_curve_metric = calcium.metrics_from_mean_curve.(metric_name);
valid_curve = isfinite(voltage_curve_metric) & isfinite(calcium_curve_metric);
valid_cycle = isfinite(voltage_cycle_mean) & isfinite(calcium_cycle_mean);
scatter(ax, voltage_curve_metric(valid_curve), calcium_curve_metric(valid_curve), ...
    42, [0.15 0.25 0.85], 'filled', 'MarkerFaceAlpha', 0.72);
hold(ax, 'on');
scatter(ax, voltage_cycle_mean(valid_cycle), calcium_cycle_mean(valid_cycle), ...
    32, [0.85 0.45 0.10], '^', 'filled', 'MarkerFaceAlpha', 0.62);
plot(ax, [0 1], [0 1], 'k--', 'LineWidth', 1);
xlim(ax, [0 1]);
ylim(ax, [0 1]);
axis(ax, 'square');
xlabel(ax, ['Voltage ', metric_label]);
ylabel(ax, ['Calcium ', metric_label]);
title(ax, sprintf('%s by ROI', metric_label));
legend(ax, {'Metric from mean curve', 'Mean of cycle metrics', 'y=x'}, 'Location', 'best');
grid(ax, 'off');
end

function plot_dsi_direction_radar_rec(ax, orientations, voltage, calcium)
voltage_population = mean(voltage.stim_mean, 1, 'omitnan');
calcium_population = mean(calcium.stim_mean, 1, 'omitnan');
plot_tuning_radar_cartesian_rec(ax, mod(double(orientations(:)'), 360), ...
    voltage_population, calcium_population, 'DSI Direction Radar', 1);
end

function plot_osi_orientation_radar_rec(ax, orientations, voltage, calcium)
voltage_population = mean(voltage.stim_mean, 1, 'omitnan');
calcium_population = mean(calcium.stim_mean, 1, 'omitnan');
[orientation_angles, voltage_orientation] = fold_orientation_response_rec(orientations, voltage_population);
[~, calcium_orientation] = fold_orientation_response_rec(orientations, calcium_population);
plot_tuning_radar_cartesian_rec(ax, orientation_angles, ...
    voltage_orientation, calcium_orientation, 'OSI Orientation Radar', 2);
end

function plot_tuning_radar_cartesian_rec(ax, angles, voltage_response, calcium_response, title_text, angle_multiplier)
angles = double(angles(:)');
voltage_values = normalize_radar_response_rec(voltage_response);
calcium_values = normalize_radar_response_rec(calcium_response);
theta = deg2rad(mod(angle_multiplier * angles, 360));
[theta, order] = sort(theta);
angles = angles(order);
voltage_values = voltage_values(order);
calcium_values = calcium_values(order);
if isempty(theta)
    text(ax, 0.5, 0.5, 'Radar unavailable', 'Units', 'normalized', 'HorizontalAlignment', 'center');
    axis(ax, 'off');
    return;
end
theta_closed = [theta(:); theta(1)];
voltage_closed = [voltage_values(:); voltage_values(1)];
calcium_closed = [calcium_values(:); calcium_values(1)];
cla(ax);
hold(ax, 'on');
axis(ax, 'equal');
axis(ax, 'off');
for radius = 0.25:0.25:1
    plot(ax, radius * cos(linspace(0, 2*pi, 181)), radius * sin(linspace(0, 2*pi, 181)), ...
        '-', 'Color', [0.86 0.86 0.86], 'LineWidth', 0.8);
end
for idx = 1:numel(theta)
    plot(ax, [0 cos(theta(idx))], [0 sin(theta(idx))], '-', 'Color', [0.78 0.78 0.78], 'LineWidth', 0.8);
    text(ax, 1.15 * cos(theta(idx)), 1.15 * sin(theta(idx)), sprintf('%.0f deg', angles(idx)), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 8);
end
h_voltage = plot(ax, voltage_closed .* cos(theta_closed), voltage_closed .* sin(theta_closed), ...
    '-o', 'Color', [0.85 0.15 0.15], 'MarkerFaceColor', [0.85 0.15 0.15], 'LineWidth', 1.7);
h_calcium = plot(ax, calcium_closed .* cos(theta_closed), calcium_closed .* sin(theta_closed), ...
    '-o', 'Color', [0.10 0.65 0.20], 'MarkerFaceColor', [0.10 0.65 0.20], 'LineWidth', 1.7);
text(ax, 0, 0, '0', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'Color', [0.45 0.45 0.45]);
text(ax, 1.02, 0.03, '1', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'Color', [0.45 0.45 0.45]);
legend(ax, [h_voltage, h_calcium], {'Voltage normalized tuning', 'Calcium normalized tuning'}, 'Location', 'southoutside');
title(ax, title_text);
xlim(ax, [-1.35 1.35]);
ylim(ax, [-1.25 1.35]);
end

function [orientation_angles, folded_response] = fold_orientation_response_rec(orientations, response)
orientation_values = mod(double(orientations(:)'), 180);
response = double(response(:)');
orientation_angles = unique(orientation_values, 'sorted');
folded_response = NaN(size(orientation_angles));
for idx = 1:numel(orientation_angles)
    mask = abs(orientation_values - orientation_angles(idx)) < 1e-9;
    folded_response(idx) = mean(response(mask), 'omitnan');
end
end

function values = normalize_radar_response_rec(values)
values = double(values(:)');
values(~isfinite(values)) = NaN;
values(values < 0) = 0;
max_value = max(values, [], 'omitnan');
if isempty(max_value) || ~isfinite(max_value) || max_value <= 0
    values(:) = 0;
else
    values = values ./ max_value;
end
end

function plot_record_grating_by_roi_rec(summary, roi_dir)
nrois = size(summary.voltage.stim_mean, 1);
for roi_idx = 1:nrois
    fig_file = fullfile(roi_dir, sprintf('r%03d_tuning.fig', roi_idx));
    png_file = fullfile(roi_dir, sprintf('r%03d_tuning.png', roi_idx));
    if skip_existing_figure_bundle_rec(fig_file, png_file, ...
            sprintf('ROI %03d record-average grating tuning', roi_idx))
        continue;
    end

    fig = figure('Color', 'w', 'Name', sprintf('Record Average Grating ROI %03d', roi_idx), ...
        'Position', [120, 100, 1100, 760]);
    tiledlayout(fig, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    plot_roi_grating_curves_rec(nexttile, summary.orientations, summary.voltage, roi_idx, ...
        [0.85 0.15 0.15], 'Voltage Spike Count', 'Spikes / trial');
    plot_roi_grating_curves_rec(nexttile, summary.orientations, summary.calcium, roi_idx, ...
        [0.10 0.65 0.20], 'Calcium Mean Sensitivity', 'Mean sensitivity');
    plot_roi_metric_pair_rec(nexttile, summary, roi_idx, 'dsi', 'DSI');
    plot_roi_metric_pair_rec(nexttile, summary, roi_idx, 'osi', 'OSI');
    sgtitle(fig, sprintf('ROI %03d Record-Average Grating Tuning', roi_idx));
    save_figure_bundle_rec(fig, fig_file, png_file);
    close(fig);
end
end

function plot_record_grating_by_roi_voltage_sensitivity_rec(summary, roi_dir)
summary_parallel = summary;
summary_parallel.voltage = summary.voltage_sensitivity;
nrois = size(summary_parallel.voltage.stim_mean, 1);
for roi_idx = 1:nrois
    fig_file = fullfile(roi_dir, sprintf('r%03d_v_sens_tuning.fig', roi_idx));
    png_file = fullfile(roi_dir, sprintf('r%03d_v_sens_tuning.png', roi_idx));
    if skip_existing_figure_bundle_rec(fig_file, png_file, ...
            sprintf('ROI %03d record-average voltage-sensitivity grating tuning', roi_idx))
        continue;
    end

    fig = figure('Color', 'w', 'Name', sprintf('Record Average Grating ROI %03d Voltage Sensitivity', roi_idx), ...
        'Position', [120, 100, 1100, 760]);
    tiledlayout(fig, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');
    plot_roi_grating_curves_rec(nexttile, summary_parallel.orientations, summary_parallel.voltage, roi_idx, ...
        [0.85 0.15 0.15], 'Voltage Mean Sensitivity', 'Mean sensitivity');
    plot_roi_grating_curves_rec(nexttile, summary_parallel.orientations, summary_parallel.calcium, roi_idx, ...
        [0.10 0.65 0.20], 'Calcium Mean Sensitivity', 'Mean sensitivity');
    plot_roi_metric_pair_rec(nexttile, summary_parallel, roi_idx, 'dsi', 'DSI');
    plot_roi_metric_pair_rec(nexttile, summary_parallel, roi_idx, 'osi', 'OSI');
    sgtitle(fig, sprintf('ROI %03d Record-Average Grating Tuning | Voltage sensitivity average', roi_idx));
    save_figure_bundle_rec(fig, fig_file, png_file);
    close(fig);
end
end

function write_record_grating_metrics_table_rec(summary, csv_file)
nrois = size(summary.voltage.stim_mean, 1);
roi = (1:nrois)';
metrics_table = table(roi, ...
    summary.voltage.metrics_mean_across_cycles.dsi.mean, ...
    summary.voltage.metrics_from_mean_curve.dsi, ...
    summary.voltage.metrics_mean_across_cycles.osi.mean, ...
    summary.voltage.metrics_from_mean_curve.osi, ...
    summary.calcium.metrics_mean_across_cycles.dsi.mean, ...
    summary.calcium.metrics_from_mean_curve.dsi, ...
    summary.calcium.metrics_mean_across_cycles.osi.mean, ...
    summary.calcium.metrics_from_mean_curve.osi, ...
    'VariableNames', {'ROI', ...
    'VoltageSpikeCount_DSI_MeanAcrossCycles', 'VoltageSpikeCount_DSI_FromMeanCurve', ...
    'VoltageSpikeCount_OSI_MeanAcrossCycles', 'VoltageSpikeCount_OSI_FromMeanCurve', ...
    'Calcium_DSI_MeanAcrossCycles', 'Calcium_DSI_FromMeanCurve', ...
    'Calcium_OSI_MeanAcrossCycles', 'Calcium_OSI_FromMeanCurve'});
if isfield(summary, 'voltage_sensitivity') && isstruct(summary.voltage_sensitivity) ...
        && isfield(summary.voltage_sensitivity, 'metrics_mean_across_cycles')
    metrics_table.VoltageSensitivity_DSI_MeanAcrossCycles = ...
        summary.voltage_sensitivity.metrics_mean_across_cycles.dsi.mean;
    metrics_table.VoltageSensitivity_DSI_FromMeanCurve = ...
        summary.voltage_sensitivity.metrics_from_mean_curve.dsi;
    metrics_table.VoltageSensitivity_OSI_MeanAcrossCycles = ...
        summary.voltage_sensitivity.metrics_mean_across_cycles.osi.mean;
    metrics_table.VoltageSensitivity_OSI_FromMeanCurve = ...
        summary.voltage_sensitivity.metrics_from_mean_curve.osi;
end
writetable(metrics_table, char(csv_file));
end

function plot_roi_grating_curves_rec(ax, orientations, channel, roi_idx, line_color, title_text, y_label)
errorbar(ax, orientations, channel.stim_mean(roi_idx, :), channel.stim_sem(roi_idx, :), ...
    '-o', 'Color', line_color, 'MarkerFaceColor', line_color, 'LineWidth', 1.5);
hold(ax, 'on');
errorbar(ax, orientations, channel.nonstim_mean(roi_idx, :), channel.nonstim_sem(roi_idx, :), ...
    '--o', 'Color', lighten_color(line_color, 0.55), 'MarkerFaceColor', 'w', 'LineWidth', 1.2);
xlabel(ax, 'Direction (deg)');
ylabel(ax, y_label);
title(ax, title_text);
legend(ax, {'Grating', 'Non-grating'}, 'Location', 'best');
grid(ax, 'off');
end

function plot_roi_metric_pair_rec(ax, summary, roi_idx, metric_name, metric_label)
values = [ ...
    summary.voltage.metrics_mean_across_cycles.(metric_name).mean(roi_idx), ...
    summary.voltage.metrics_from_mean_curve.(metric_name)(roi_idx); ...
    summary.calcium.metrics_mean_across_cycles.(metric_name).mean(roi_idx), ...
    summary.calcium.metrics_from_mean_curve.(metric_name)(roi_idx)];
x_positions = [1, 2];
plot(ax, x_positions, values(:, 1), '-o', ...
    'Color', [0.85 0.45 0.10], 'MarkerFaceColor', [0.85 0.45 0.10], 'LineWidth', 1.3);
hold(ax, 'on');
plot(ax, x_positions, values(:, 2), '-o', ...
    'Color', [0.15 0.25 0.85], 'MarkerFaceColor', [0.15 0.25 0.85], 'LineWidth', 1.5);
xlim(ax, [0.75 2.25]);
ylim(ax, [0 1]);
set(ax, 'XTick', x_positions, 'XTickLabel', {'Voltage', 'Calcium'});
ylabel(ax, metric_label);
title(ax, metric_label);
legend(ax, {'Mean of cycle metrics', 'Metric from mean curve'}, 'Location', 'best');
grid(ax, 'off');
end

function result = plot_record_roi_cycle_heatmap_voltage_trace( ...
    voltage_stage, calcium_stage, stim_windows, roi_dir, mat_file, ...
    voltage_polarity, calcium_polarity, calcium_smoothing_window, voltage_metric_name)
if ~isfolder(roi_dir)
    mkdir(roi_dir);
end
voltage_metric_name = lower(string(voltage_metric_name));
voltage_metric_label = upper(voltage_metric_name);
if voltage_metric_name == "sensitivity"
    voltage_metric_label = "Sensitivity";
end

nrois = min(voltage_stage.nrois, calcium_stage.nrois);
if nrois < 1
    error('No paired ROI columns are available for ROI-cycle heatmap plotting.');
end
if voltage_stage.ncycles ~= calcium_stage.ncycles
    error('Cycle count mismatch between voltage and calcium sensitivity stages.');
end
cycle_names = calcium_stage.cycle_names(:);
ncycles = numel(cycle_names);
voltage_display = double(voltage_polarity) * double(voltage_stage.per_cycle(:, 1:nrois, :));
calcium_display_raw = double(calcium_polarity) * double(calcium_stage.per_cycle(:, 1:nrois, :));
calcium_windows = resolve_optional_channel_windows(stim_windows, 'calcium');

png_files = strings(nrois, 1);
fig_files = strings(nrois, 1);
calcium_zero_percentile = 1;
voltage_trace_row_scale = 1.0;
calcium_cmap = custom_ice_colormap_rec(256);

for roi_idx = 1:nrois
    fig_files(roi_idx) = string(fullfile(roi_dir, ...
        sprintf('r%03d_ca_hm_v_%s.fig', roi_idx, voltage_metric_name)));
    png_files(roi_idx) = string(fullfile(roi_dir, ...
        sprintf('r%03d_ca_hm_v_%s.png', roi_idx, voltage_metric_name)));
    if skip_existing_figure_bundle_rec(fig_files(roi_idx), png_files(roi_idx), ...
            sprintf('ROI %03d cycle calcium heatmap voltage %s trace', roi_idx, voltage_metric_name))
        continue;
    end

    calcium_cycles = squeeze(calcium_display_raw(:, roi_idx, :));
    voltage_cycles = squeeze(voltage_display(:, roi_idx, :));
    if ncycles == 1
        calcium_cycles = reshape(calcium_cycles, [], 1);
        voltage_cycles = reshape(voltage_cycles, [], 1);
    end

    calcium_baseline = prctile(calcium_cycles, calcium_zero_percentile, 1);
    calcium_baseline(~isfinite(calcium_baseline)) = 0;
    calcium_display_zeroed = calcium_cycles - calcium_baseline;
    calcium_display_zeroed(calcium_display_zeroed < 0) = 0;
    finite_calcium = calcium_display_zeroed(isfinite(calcium_display_zeroed));
    if isempty(finite_calcium) || max(finite_calcium) <= 0
        calcium_clim = [0, 1];
    else
        calcium_clim = [0, max(finite_calcium)];
    end

    voltage_cycle_min = min(voltage_cycles, [], 1, 'omitnan');
    voltage_cycle_max = max(voltage_cycles, [], 1, 'omitnan');
    voltage_cycle_range = voltage_cycle_max - voltage_cycle_min;
    voltage_global_range = max(voltage_cycle_range, [], 'omitnan');
    if ~isfinite(voltage_global_range) || voltage_global_range <= 0
        voltage_global_range = 1;
    end
    voltage_cycle_mid = (voltage_cycle_min + voltage_cycle_max) / 2;
    voltage_cycle_mid(~isfinite(voltage_cycle_mid)) = 0;

    fig_height = min(1600, max(520, 42 * ncycles + 180));
    fig = figure('Color', 'w', ...
        'Name', sprintf('Record ROI %03d Cycle Calcium Heatmap With Voltage %s Trace', roi_idx, voltage_metric_label), ...
        'Units', 'pixels', ...
        'Position', [100, 100, 1400, fig_height]);
    ax = axes(fig);
    imagesc(ax, calcium_stage.time, 1:ncycles, calcium_display_zeroed');
    set(ax, 'YDir', 'normal', 'TickDir', 'out', 'Layer', 'top');
    colormap(ax, calcium_cmap);
    clim(ax, calcium_clim);
    hold(ax, 'on');
    for cycle_idx = 1:ncycles
        y_trace = cycle_idx + ...
            (voltage_cycles(:, cycle_idx) - voltage_cycle_mid(cycle_idx)) ./ voltage_global_range .* voltage_trace_row_scale;
        plot(ax, voltage_stage.time, y_trace, 'Color', [1.0, 0.12, 0.02], 'LineWidth', 0.65);
    end
    xlim(ax, [min(calcium_stage.time), max(calcium_stage.time)]);
    ylim(ax, [0.5, ncycles + 0.5]);
    set(ax, 'YTick', 1:ncycles, 'YTickLabel', cycle_names, 'TickLabelInterpreter', 'none');
    if ~isempty(fieldnames(calcium_windows))
        overlay_stim_boundaries_rec(ax, calcium_windows);
    end
    xlabel(ax, 'Time (s)');
    ylabel(ax, 'Cycle');
    title(ax, sprintf(['ROI %03d | Calcium Sensitivity Heatmap (%s, window=%d) ' ...
        '+ Voltage %s Trace'], ...
        roi_idx, strrep(char(calcium_stage.stage_name), '_', '\_'), calcium_smoothing_window, voltage_metric_label));
    cb = colorbar(ax);
    ylabel(cb, 'Calcium sensitivity display value');
    grid(ax, 'off');
    box(ax, 'off');

    save_figure_bundle_preserve_layout_rec(fig, fig_files(roi_idx), png_files(roi_idx));
    close(fig);
end

result = struct( ...
    'status', "completed", ...
    'roi_dir', string(roi_dir), ...
    'mat_file', string(mat_file), ...
    'fig_files', fig_files, ...
    'png_files', png_files, ...
    'roi_count_plotted', nrois, ...
    'cycle_names', cycle_names, ...
    'input_stages', struct('voltage', string(voltage_stage.stage_name), 'calcium', string(calcium_stage.stage_name)), ...
    'voltage_metric_name', voltage_metric_name, ...
    'display_rules', struct( ...
        'calcium_zero_rule', "per-ROI per-cycle 1st percentile subtracted, then values below zero clipped", ...
        'calcium_zero_percentile', calcium_zero_percentile, ...
        'calcium_colormap', "custom_ice", ...
        'voltage_scale_rule', "per-ROI global max-min across cycles; each cycle centered by its own midpoint", ...
        'voltage_trace_metric', voltage_metric_name, ...
        'voltage_trace_row_scale', voltage_trace_row_scale, ...
        'x_axis_source', "calcium time axis; voltage trace uses voltage seconds without resampling"), ...
    'created_at', datetime("now"));
save(mat_file, 'result', 'voltage_display', 'calcium_display_raw', ...
    'voltage_polarity', 'calcium_polarity', 'calcium_smoothing_window', '-v7.3');
end

function result = plot_record_roi_cycle_calcium_peak_raster( ...
    voltage_stage, calcium_stage, cycle_entries, stim_windows, roi_dir, mat_file, calcium_polarity, metric_name)
if ~isfolder(roi_dir)
    mkdir(roi_dir);
end
metric_name = lower(string(metric_name));
metric_label = upper(metric_name);
if metric_name == "sensitivity"
    metric_label = "Sensitivity";
end

nrois = min(voltage_stage.nrois, calcium_stage.nrois);
ncycles = calcium_stage.ncycles;
cycle_names = calcium_stage.cycle_names(:);
calcium_display = double(calcium_polarity) * double(calcium_stage.per_cycle(:, 1:nrois, :));
calcium_windows = resolve_optional_channel_windows(stim_windows, 'calcium');

png_files = strings(nrois, 1);
fig_files = strings(nrois, 1);
peak_time_by_roi_cycle = cell(nrois, ncycles);
missing_peak_cycle = false(ncycles, 1);
trace_row_scale = 0.34;
tick_half_height = 0.34;

for roi_idx = 1:nrois
    fig_files(roi_idx) = string(fullfile(roi_dir, ...
        sprintf('r%03d_ca_%s_peak.fig', roi_idx, metric_name)));
    png_files(roi_idx) = string(fullfile(roi_dir, ...
        sprintf('r%03d_ca_%s_peak.png', roi_idx, metric_name)));
    if skip_existing_figure_bundle_rec(fig_files(roi_idx), png_files(roi_idx), ...
            sprintf('ROI %03d cycle calcium %s peak raster', roi_idx, metric_name))
        continue;
    end

    calcium_cycles = squeeze(calcium_display(:, roi_idx, :));
    if ncycles == 1
        calcium_cycles = reshape(calcium_cycles, [], 1);
    end

    trace_cycle_min = min(calcium_cycles, [], 1, 'omitnan');
    trace_cycle_max = max(calcium_cycles, [], 1, 'omitnan');
    trace_cycle_range = trace_cycle_max - trace_cycle_min;
    trace_global_range = max(trace_cycle_range, [], 'omitnan');
    if ~isfinite(trace_global_range) || trace_global_range <= 0
        trace_global_range = 1;
    end
    trace_cycle_mid = (trace_cycle_min + trace_cycle_max) / 2;
    trace_cycle_mid(~isfinite(trace_cycle_mid)) = 0;

    fig_height = min(1600, max(520, 42 * ncycles + 180));
    fig = figure('Color', 'w', ...
        'Name', sprintf('Record ROI %03d Cycle Calcium %s Peak Raster', roi_idx, metric_label), ...
        'Units', 'pixels', ...
        'Position', [100, 100, 1400, fig_height]);
    ax = axes(fig);
    hold(ax, 'on');
    if ~isempty(fieldnames(calcium_windows))
        add_stim_shading_rec(ax, calcium_windows, 0.08);
    end

    for cycle_idx = 1:ncycles
        current_trace = calcium_cycles(:, cycle_idx);
        y_trace = cycle_idx + ...
            (current_trace - trace_cycle_mid(cycle_idx)) ./ trace_global_range .* trace_row_scale;
        plot(ax, calcium_stage.time, y_trace, 'Color', [0.10 0.55 0.18], 'LineWidth', 0.65);

        peak_times = resolve_cycle_roi_peak_times_rec( ...
            cycle_entries(cycle_idx), roi_idx, voltage_stage.alignment_time_by_cycle(cycle_idx), ...
            [min(calcium_stage.time), max(calcium_stage.time)]);
        peak_time_by_roi_cycle{roi_idx, cycle_idx} = peak_times;
        if cycle_entries(cycle_idx).voltage.peaks.status ~= "loaded"
            missing_peak_cycle(cycle_idx) = true;
        end
        for peak_idx = 1:numel(peak_times)
            line(ax, [peak_times(peak_idx), peak_times(peak_idx)], ...
                [cycle_idx - tick_half_height, cycle_idx + tick_half_height], ...
                'Color', [0.02 0.02 0.02], 'LineWidth', 0.8);
        end
    end

    xlim(ax, [min(calcium_stage.time), max(calcium_stage.time)]);
    ylim(ax, [0.5, ncycles + 0.5]);
    set(ax, 'YTick', 1:ncycles, 'YTickLabel', cycle_names, 'TickLabelInterpreter', 'none');
    if ~isempty(fieldnames(calcium_windows))
        overlay_stim_boundaries_rec(ax, calcium_windows);
    end
    xlabel(ax, 'Time (s)');
    ylabel(ax, 'Cycle');
    title(ax, sprintf('ROI %03d | Calcium %s Stack + Accepted Voltage Peak Raster', roi_idx, metric_label));
    grid(ax, 'off');
    box(ax, 'off');

    save_figure_bundle_preserve_layout_rec(fig, fig_files(roi_idx), png_files(roi_idx));
    close(fig);
end

result = struct( ...
    'status', "completed", ...
    'roi_dir', string(roi_dir), ...
    'mat_file', string(mat_file), ...
    'fig_files', fig_files, ...
    'png_files', png_files, ...
    'roi_count_plotted', nrois, ...
    'cycle_names', cycle_names, ...
    'input_stage', string(calcium_stage.stage_name), ...
    'metric_name', metric_name, ...
    'missing_peak_cycle_names', cycle_names(missing_peak_cycle), ...
    'display_rules', struct( ...
        'peak_source', "voltage_results.peak_results.accepted_for_events, fallback voltage_peak_results.mat", ...
        'peak_mark', "vertical tick line at accepted peak time", ...
        'trace_stage', string(calcium_stage.stage_name), ...
        'trace_metric', metric_name, ...
        'peak_time_rule', "peak frame index / frame_rate - cycle alignment time", ...
        'trace_display', "stacked calcium sensitivity trace behind peak ticks", ...
        'x_axis_source', "record-average aligned calcium time axis"), ...
    'created_at', datetime("now"));
save(mat_file, 'result', 'calcium_display', 'peak_time_by_roi_cycle', ...
    'calcium_polarity', '-v7.3');
end

function peak_times = resolve_cycle_roi_peak_times_rec(cycle_entry, roi_idx, alignment_time, time_limits)
peak_times = [];
peaks = cycle_entry.voltage.peaks;
if ~isstruct(peaks) || peaks.status ~= "loaded" || numel(peaks.index) < roi_idx
    return;
end
peak_idx = double(peaks.index{roi_idx}(:));
peak_idx = peak_idx(isfinite(peak_idx));
if isempty(peak_idx)
    return;
end
frame_rate = cycle_entry.voltage.sensitivity.frame_rate;
peak_times = peak_idx ./ frame_rate - alignment_time;
peak_times = peak_times(peak_times >= time_limits(1) & peak_times <= time_limits(2));
end

function cmap = custom_ice_colormap_rec(ncolors)
if nargin < 1
    ncolors = 256;
end
cmap = zeros(ncolors, 3);
cmap(:, 3) = linspace(0, 1, ncolors);
cyan_start_node = max(1, floor(ncolors * 0.1));
cmap(cyan_start_node:end, 2) = linspace(0, 1, ncolors - cyan_start_node + 1);
white_start_node = max(1, floor(ncolors * 0.9));
cmap(white_start_node:end, 1) = linspace(0, 1, ncolors - white_start_node + 1);
end

function plot_record_stage_summary(voltage_stage, calcium_stage, stim_windows, metric_name, output_png, voltage_color, calcium_color, voltage_display_scale, calcium_display_scale, voltage_title, calcium_title)
if skip_existing_figure_bundle_rec(strrep(output_png, '.png', '.fig'), output_png, ...
        sprintf('record %s stage summary', metric_name))
    return;
end
fig = figure('Color', 'w', 'Name', sprintf('Record Average %s Summary', metric_name), ...
    'Position', [80, 80, 1600, 950]);
tiledlayout(fig, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

voltage_cycle_traces = voltage_display_scale * voltage_stage.cycle_roi_mean;
calcium_cycle_traces = calcium_display_scale * calcium_stage.cycle_roi_mean;

ax_voltage_heatmap = nexttile;
plot_record_cycle_heatmap( ...
    ax_voltage_heatmap, voltage_stage.time, voltage_cycle_traces, ...
    voltage_stage.cycle_names, string(voltage_title) + " Heatmap", metric_name, ...
    resolve_optional_channel_windows(stim_windows, 'voltage'));

ax_voltage_stack = nexttile;
[voltage_ylabel, voltage_scale_suffix] = describe_record_stage_axis_rec(metric_name);
plot_record_cycle_stack( ...
    ax_voltage_stack, voltage_stage.time, voltage_cycle_traces, ...
    voltage_stage.cycle_names, voltage_color, string(voltage_title) + " Stacked Trace", ...
    resolve_optional_channel_windows(stim_windows, 'voltage'), ...
    voltage_ylabel, voltage_scale_suffix);

ax_calcium_heatmap = nexttile;
plot_record_cycle_heatmap( ...
    ax_calcium_heatmap, calcium_stage.time, calcium_cycle_traces, ...
    calcium_stage.cycle_names, string(calcium_title) + " Heatmap", metric_name, ...
    resolve_optional_channel_windows(stim_windows, 'calcium'));

ax_calcium_stack = nexttile;
[calcium_ylabel, calcium_scale_suffix] = describe_record_stage_axis_rec(metric_name);
plot_record_cycle_stack( ...
    ax_calcium_stack, calcium_stage.time, calcium_cycle_traces, ...
    calcium_stage.cycle_names, calcium_color, string(calcium_title) + " Stacked Trace", ...
    resolve_optional_channel_windows(stim_windows, 'calcium'), ...
    calcium_ylabel, calcium_scale_suffix);

sgtitle(fig, sprintf('Record %s Summary | rows/traces = ROI-mean of each of %d cycles', ...
    metric_name, voltage_stage.ncycles));
save_figure_bundle_rec(fig, strrep(output_png, '.png', '.fig'), output_png);
close(fig);
end

function plot_record_cycle_heatmap(ax, t_axis, cycle_traces, cycle_names, title_text, metric_name, channel_windows)
imagesc(ax, t_axis, 1:size(cycle_traces, 2), cycle_traces');
axis(ax, 'xy');
xlabel(ax, 'Time (s)');
ylabel(ax, 'Cycle');
title(ax, title_text);
set(ax, 'YTick', 1:numel(cycle_names), 'YTickLabel', cycle_names, ...
    'TickLabelInterpreter', 'none');
colorbar(ax);
finite_values = cycle_traces(isfinite(cycle_traces));
if ~isempty(finite_values)
    if any(strcmpi(string(metric_name), ["sensitivity", "snr"]))
        color_limit = max(abs(finite_values));
        if isfinite(color_limit) && color_limit > 0
            caxis(ax, [-color_limit, color_limit]);
        end
        colormap(ax, redblue_colormap_rec(256));
    else
        color_limits = [min(finite_values), max(finite_values)];
        if all(isfinite(color_limits)) && color_limits(2) > color_limits(1)
            caxis(ax, color_limits);
        end
        colormap(ax, turbo(256));
    end
end
if ~isempty(fieldnames(channel_windows))
    overlay_stim_boundaries_rec(ax, channel_windows);
end
end

function plot_record_cycle_stack(ax, t_axis, cycle_traces, cycle_names, line_color, title_text, channel_windows, y_axis_label, y_scale_suffix)
[~, ncycles] = size(cycle_traces);
hold(ax, 'on');

stack_spacing = compute_stack_spacing_rec(cycle_traces);
x_min = min(t_axis);
x_max = max(t_axis);
x_span = max(eps, x_max - x_min);
left_margin = 0.14 * x_span;
ymin_global = inf;
ymax_global = -inf;

for cycle_idx = 1:ncycles
    offset = (ncycles - cycle_idx) * stack_spacing;
    current_trace = cycle_traces(:, cycle_idx);
    ymin_global = min(ymin_global, min(current_trace, [], 'omitnan') + offset);
    ymax_global = max(ymax_global, max(current_trace, [], 'omitnan') + offset);
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

for cycle_idx = 1:ncycles
    offset = (ncycles - cycle_idx) * stack_spacing;
    current_trace = cycle_traces(:, cycle_idx);
    plot(ax, t_axis, current_trace + offset, ...
        'Color', line_color, 'LineWidth', 1.4);
    y_label = offset + mean(current_trace, 'omitnan');
    if ~isfinite(y_label)
        y_label = offset;
    end
    text(ax, x_min - 0.92 * left_margin, y_label, char(cycle_names(cycle_idx)), ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
        'FontSize', 8, 'FontWeight', 'bold', 'Color', [0.15 0.15 0.15], ...
        'Interpreter', 'none');
end

title(ax, title_text);
xlabel(ax, 'Time (s)');
ylabel(ax, y_axis_label);
set(ax, 'YTick', []);
grid(ax, 'off');
box(ax, 'off');
add_axis_scalebar_rec(ax, t_axis, cycle_traces(:), [0.1 0.1 0.1], y_scale_suffix);
end

function cmap = redblue_colormap_rec(ncolors)
if nargin < 1
    ncolors = 256;
end
half_count = ceil(ncolors / 2);
blue_to_white = [linspace(0, 1, half_count)', linspace(0, 1, half_count)', ones(half_count, 1)];
red_count = ncolors - half_count;
white_to_red = [ones(red_count, 1), linspace(1, 0, red_count)', linspace(1, 0, red_count)'];
cmap = [blue_to_white; white_to_red];
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
if skip_existing_figure_bundle_rec(strrep(output_png, '.png', '.fig'), output_png, ...
        sprintf('record %s metric overlap', metric_name))
    return;
end
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
grid(ax, 'off');
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
grid(ax, 'off');
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
fig_file = fullfile(output_dir, '6_fft.fig');
png_file = fullfile(output_dir, '6_fft.png');
if skip_existing_figure_bundle_rec(fig_file, png_file, 'record-average FFT summary')
    return;
end
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

save_figure_bundle_rec(fig, fig_file, png_file);
close(fig);
end

function [fig_file, png_file] = plot_population_wavelet_summary_rec(voltage_tf, calcium_tf, stim_windows, output_dir)
fig_file = fullfile(output_dir, '7_wav.fig');
png_file = fullfile(output_dir, '7_wav.png');
if skip_existing_figure_bundle_rec(fig_file, png_file, 'record-average wavelet summary')
    return;
end
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
grid off;
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
grid off;
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
grid off;
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
grid(ax, 'off');
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
grid(ax, 'off');
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
disable_axes_toolbar_for_export_rec(fig_handle);
save_figure_outputs_if_missing_rec(fig_handle, fig_file, png_file);
end

function save_figure_bundle_preserve_layout_rec(fig_handle, fig_file, png_file)
drawnow;
disable_axes_toolbar_for_export_rec(fig_handle);
save_figure_outputs_if_missing_rec(fig_handle, fig_file, png_file);
end

function save_figure_outputs_if_missing_rec(fig_handle, fig_file, png_file)
save_fig_if_missing_rec(fig_handle, fig_file);
export_png_if_missing_rec(fig_handle, png_file);
export_emf_backup_if_missing_rec(fig_handle, png_file);
end

function tf = skip_existing_figure_bundle_rec(fig_file, png_file, description)
tf = figure_bundle_complete_rec(fig_file, png_file);
if tf
    fprintf('Figure bundle already exists, skipping plot generation');
    if nargin >= 3 && strlength(string(description)) > 0
        fprintf(' (%s)', char(string(description)));
    end
    fprintf(':\n  %s\n  %s\n  %s\n', ...
        char(string(fig_file)), char(string(png_file)), ...
        char(replace_file_extension_rec(png_file, '.emf')));
end
end

function tf = figure_bundle_complete_rec(fig_file, png_file)
tf = file_has_content_rec(fig_file) ...
    && file_has_content_rec(png_file) ...
    && file_has_content_rec(replace_file_extension_rec(png_file, '.emf'));
end

function save_fig_if_missing_rec(fig_handle, fig_file)
fig_file = char(string(fig_file));
if file_has_content_rec(fig_file)
    fprintf('FIG already exists, skipping FIG save:\n  %s\n', fig_file);
    return;
end
[fig_folder, fig_name, fig_ext] = fileparts(fig_file);
if ~isfolder(fig_folder)
    mkdir(fig_folder);
end
temp_fig_file = fullfile(fig_folder, sprintf('%s_tmp_%s%s', ...
    fig_name, char(datetime('now', 'Format', 'yyyyMMdd_HHmmss_SSS')), fig_ext));
try
    savefig(fig_handle, temp_fig_file);
    if ~file_has_content_rec(temp_fig_file)
        error('Dual_rec_analysis3:FigSaveEmpty', ...
            'Temporary FIG is missing or zero bytes: %s', temp_fig_file);
    end
    [move_ok, move_message] = movefile(temp_fig_file, fig_file, 'f');
    if ~move_ok
        error('Dual_rec_analysis3:FigMoveFailed', ...
            'Could not replace FIG file %s: %s', fig_file, move_message);
    end
catch ME
    if isfile(temp_fig_file)
        delete(temp_fig_file);
    end
    if isfile(fig_file) && ~file_has_content_rec(fig_file)
        delete(fig_file);
    end
    warning('Dual_rec_analysis3:FigSaveFailed', ...
        'Could not save FIG file, continuing with PNG/EMF export: %s\n%s', ...
        fig_file, ME.message);
end
end

function tf = file_has_content_rec(file_path)
file_path = char(string(file_path));
file_info = dir(file_path);
tf = numel(file_info) == 1 && file_info.bytes > 0;
end

function export_png_if_missing_rec(fig_handle, png_file)
png_file = char(string(png_file));
if file_has_content_rec(png_file)
    fprintf('PNG already exists, skipping PNG export:\n  %s\n', png_file);
    return;
end
try
    exportgraphics(fig_handle, png_file, 'Resolution', 150);
catch ME
    warning('Dual_rec_analysis3:PngExportFailed', ...
        'Could not export PNG file, continuing with EMF export: %s\n%s', ...
        png_file, ME.message);
end
end

function export_emf_backup_if_missing_rec(fig_handle, png_file)
emf_file = char(replace_file_extension_rec(png_file, '.emf'));
if file_has_content_rec(emf_file)
    fprintf('EMF already exists, skipping EMF export:\n  %s\n', emf_file);
    return;
end
[emf_folder, emf_name, emf_ext] = fileparts(emf_file);
if ~isfolder(emf_folder)
    mkdir(emf_folder);
end
temp_emf_file = fullfile(tempdir, sprintf('%s_%s%s', ...
    emf_name, char(datetime('now', 'Format', 'yyyyMMdd_HHmmss_SSS')), emf_ext));
try
    exportgraphics(fig_handle, temp_emf_file, 'ContentType', 'vector');
catch ME
    warning('Dual_rec_analysis3:EmfExportFailed', ...
        'exportgraphics could not export EMF backup; trying print -dmeta fallback: %s\n%s', ...
        emf_file, ME.message);
end
if ~file_has_content_rec(temp_emf_file)
    if isfile(temp_emf_file)
        delete(temp_emf_file);
    end
    try
        drawnow;
        print(fig_handle, temp_emf_file, '-dmeta', '-r150');
    catch ME
        warning('Dual_rec_analysis3:EmfPrintFallbackFailed', ...
            'Could not export EMF backup with print -dmeta, continuing: %s\n%s', ...
            emf_file, ME.message);
    end
end
if file_has_content_rec(temp_emf_file)
    [move_ok, move_message] = movefile(temp_emf_file, emf_file, 'f');
    if ~move_ok
        warning('Dual_rec_analysis3:EmfMoveFailed', ...
            'Could not move temporary EMF backup to target, leaving temp file: %s\nTarget: %s\n%s', ...
            temp_emf_file, emf_file, move_message);
    end
end
if ~file_has_content_rec(emf_file)
    warning('Dual_rec_analysis3:EmfFileMissingAfterExport', ...
        'EMF export finished without creating a file: %s', emf_file);
end
end

function output_file = replace_file_extension_rec(input_file, new_extension)
input_file = char(string(input_file));
new_extension = char(string(new_extension));
if ~startsWith(new_extension, '.')
    new_extension = ['.', new_extension];
end
[input_folder, input_name] = fileparts(input_file);
output_file = string(fullfile(input_folder, [input_name, new_extension]));
end

function save_mat_file_resilient_rec(mat_file, data_struct, varargin)
mat_file = char(string(mat_file));
[mat_folder, mat_name, mat_ext] = fileparts(mat_file);
if isempty(mat_ext)
    mat_ext = '.mat';
end
if ~isfolder(mat_folder)
    mkdir(mat_folder);
end
temp_file = fullfile(mat_folder, sprintf('%s_tmp_%s%s', ...
    mat_name, char(datetime('now', 'Format', 'yyyyMMdd_HHmmss_SSS')), mat_ext));
try
    save(temp_file, '-struct', 'data_struct', varargin{:});
    [move_ok, move_message] = movefile(temp_file, mat_file, 'f');
    if ~move_ok
        warning('Dual_rec_analysis3:MatMoveFailed', ...
            'Could not replace MAT file after temporary save, leaving temp file: %s\n%s', ...
            temp_file, move_message);
    end
catch ME
    if isfile(temp_file)
        delete(temp_file);
    end
    warning('Dual_rec_analysis3:MatSaveFailed', ...
        'Could not save MAT file, continuing: %s\n%s', mat_file, ME.message);
end
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

function disable_axes_toolbar_for_export_rec(fig_handle)
if ~isgraphics(fig_handle, 'figure')
    return;
end
axes_handles = findall(fig_handle, 'Type', 'axes');
for ax_idx = 1:numel(axes_handles)
    ax = axes_handles(ax_idx);
    try
        ax.Toolbar.Visible = 'off';
    catch
    end
    try
        disableDefaultInteractivity(ax);
    catch
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

function files = list_record_average_files(output_dir)
listing = dir(fullfile(char(output_dir), '**', '*'));
listing = listing(~[listing.isdir]);
files = string(arrayfun(@(entry) fullfile(entry.folder, entry.name), ...
    listing, 'UniformOutput', false));
files = files(:);
end

function emit_record_average_event(runtime, state, message, output_path)
event = struct( ...
    'timestamp', datetime('now'), ...
    'type', "section", ...
    'mode', "dual", ...
    'scope', "record", ...
    'input_path', "", ...
    'cycle', "", ...
    'section', "record_average", ...
    'state', string(state), ...
    'message', string(message), ...
    'output_path', string(output_path), ...
    'details', struct());
aaa.io.emit_event(runtime, event);
end
