% Dual_analysis3 - Dual-camera voltage/calcium companion analysis
% Dual_analysis3 notes:
%   - This script is intended for dual-camera acquisitions saved by the
%     rebuilt Hamamatsu workflow under Rec*/Cycle*/Cam*_label paths.
%   - It prefers cycle_manifest.mat / record_manifest.mat when available,
%     then falls back to Cam1_* / Cam2_* folders inside the selected cycle.
%   - It can also pair two raw channel folders directly, for example:
%       ...\20240905-153142POA
%       ...\20240905-153142POA_Green
%     In that case point cycle_path to either folder and keep
%     raw_dual_green_suffix = '_Green'.
%   - Each camera has its own explicit transpose flag. There is no shared
%     transpose variable, so the two camera orientations cannot be mixed
%     up silently across sections.
%   - By default this script follows the old dual-color workflow:
%       Cam1 -> calcium, pagetranspose applied
%       Cam2 -> voltage, no transpose applied
%     If your wiring changed, edit camera_cfg only.
%   - Movies are not saved again. Lightweight movie_info is stored inside
%     voltage_results / calcium_results instead.
%   - The main trace stages use the unified names:
%     raw, bg_removed, bleach_removed, baseline, noise_reference, noise,
%     sensitivity, and snr.
%   - Calcium-only final smoothing stages are stored explicitly as
%     raw_smoothed, sensitivity_smoothed, and snr_smoothed after the main
%     metric stages are computed.
%   - Motion correction is estimated once on the voltage movie and the
%     same shifts are applied to the calcium movie after transpose/resize.
%   - Main output files are dual_info.mat, voltage_results.mat,
%     calcium_results.mat, dual_results.mat, stim_results.mat, and
%     -1_explicit_dual_results.mat.

% Preserve externally provided config variables so this script can still be
% driven by an outer batch script via:
%   cycle_path = ...;
%   reuse_roi_file = ...;
%   correct_offset_mode = 'none';
%   run('Dual_analysis3.m')
%
% Important:
%   this script no longe
% r clears the workspace 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% automatically.
%   - If you run it directly by itself, CLEAR old variables manually first.
%   - If you run it from an outer batch script, do not clear the workspace,
%     because the outer script may be intentionally passing config values in.
% clc;

%% Helper / Readme
% 1. Point cycle_path to one rebuilt acquisition cycle folder.
% 2. Confirm camera_cfg matches your real camera-role mapping.
% 3. Only change transpose flags here; later sections inherit them.
% 4. If motion shifts or ROI files already exist, fill the reuse paths.
% 5. This script focuses on per-channel trace processing plus dual-channel
%    comparison. AP event-level refinement remains in AP_analysis3.
% 6. Every trace-processing section now saves a compact "section record"
%    alongside its usual outputs. The record is designed to answer two
%    practical questions when reopening old results:
%    "What exactly went into this section?" and
%    "What would I need to rerun it?"
%
% Section record interface:
%   record.section_name
%   record.description
%   record.input
%   record.parameters
%   record.output
%   record.excluded_inputs
%   record.rerun
%
% Example:
%   metric_record.input.traces_calcium_bleach
%   metric_record.parameters.calcium_smoothing_parameters.window
%   metric_record.output.calcium_snr_smoothed
%
% Movies are intentionally excluded from these records to avoid duplicating
% large data on disk. When a section depends on movie data, the record says
% so explicitly in excluded_inputs / rerun.note.
% 7. This script does not call clear automatically anymore. For standalone
%    use, clear the workspace manually before running to avoid accidental
%    reuse of old variables. For outer batch scripts, keep the workspace
%    intact so override parameters can be passed in on purpose.
%
% Outer-script override examples:
%   cycle_path = '...\Rec1_...\Cycle3';
%   reuse_roi_file = '...\Cycle1\Dual_analysis3\...\1_dual_roi_results.mat';
%   correct_offset_mode = 'none';
%   run('Dual_analysis3.m');
%
% Quick-start for two-folder raw dual-channel input:
%   cycle_path = 'I:\1_Data\...\20240905-153142POA';
%   raw_dual_green_suffix = '_Green';
%   voltage_transpose_movie = false;
%   calcium_transpose_movie = false;
%   run_motion_correction = false;
%   run('Dual_analysis3.m');
%
% Additional override hooks for legacy datasets:
%   camera_source_override -> explicit movie paths for the two cameras
%   stim_context_override  -> prebuilt stimulus metadata/logs struct

%% Input Setup
% This block collects the main user-tunable analysis settings in one place.
% Most day-to-day reruns should only require edits in this section.
nowtime = string(datetime('now'));
nowtime = strrep(nowtime, ':', '-');
fprintf('Initializing dual analysis...\n');

% -------------------------------------------------------------------------
% A. Main input / output paths
% -------------------------------------------------------------------------
if ~exist('cycle_path', 'var') || isempty(cycle_path)
    cycle_path = 'E:\1_Data\YHY\260520_VIP-POA_NAVI2ST-PG8_sCy5\Methods1_S1R1_1min-test\Rec1_2026-05-20_23-22-53\Cycle1';
end
[record_path_for_name, cycle_name_for_name] = fileparts(cycle_path);
[~, record_name_for_name] = fileparts(record_path_for_name);
if ~exist('analysis_run_name', 'var') || isempty(analysis_run_name)
    analysis_run_name = sprintf('%s_%s_%s', record_name_for_name, cycle_name_for_name, char(nowtime));
end
if ~exist('analysis_mode', 'var') || isempty(analysis_mode)
    analysis_mode = 'full';          % 'full' | 'analysis_only'
end
if ~exist('analysis_backend', 'var') || isempty(analysis_backend)
    analysis_backend = 'default';    % 'default' | 'volpy_voltage_reanalysis'
end
if ~exist('reuse_results_path', 'var') || isempty(reuse_results_path)
    reuse_results_path = '';         % existing Dual_analysis3 output folder for ROI-after analysis reruns
end
if ~exist('volpy_source_results_path', 'var') || isempty(volpy_source_results_path)
    volpy_source_results_path = '';  % existing standard Dual_analysis3 result folder used as calcium/stim/ROI source
end
if ~exist('volpy_auto_run', 'var') || isempty(volpy_auto_run)
    volpy_auto_run = false;
end
if ~exist('volpy_force_rerun', 'var') || isempty(volpy_force_rerun)
    volpy_force_rerun = false;
end
if ~exist('volpy_use_existing', 'var') || isempty(volpy_use_existing)
    volpy_use_existing = false;
end
if ~exist('volpy_flip_signal', 'var') || isempty(volpy_flip_signal)
    volpy_flip_signal = false;
end

% -------------------------------------------------------------------------
% B. Common manual channel settings
% Edit these first when the two movies need different transpose choices,
% frame rates, or when motion correction should be skipped deliberately.
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

% -------------------------------------------------------------------------
% C. Raw two-folder dual-channel compatibility
% When cycle_path is not a Rec*/Cycle* folder, the script can auto-pair one
% "primary" folder with its sibling "..._Green" folder.
% Default meaning:
%   primary folder -> voltage
%   *_Green folder -> calcium
% Change these only if your acquisition naming means something else.
% -------------------------------------------------------------------------
if ~exist('raw_dual_green_suffix', 'var') || isempty(raw_dual_green_suffix)
    raw_dual_green_suffix = '_Green';
end
if ~exist('raw_dual_primary_role', 'var') || isempty(raw_dual_primary_role)
    raw_dual_primary_role = "voltage";
end
if ~exist('raw_dual_green_role', 'var') || isempty(raw_dual_green_role)
    raw_dual_green_role = "calcium";
end
raw_dual_input_cfg = struct( ...
    'green_suffix', string(raw_dual_green_suffix), ...
    'primary_role', string(raw_dual_primary_role), ...
    'green_role', string(raw_dual_green_role));

% -------------------------------------------------------------------------
% D. Advanced per-camera config
% Keep camera_cfg as the single source of truth used downstream, but map
% the simpler role-level settings above onto it so users do not have to
% edit the struct for routine transpose / fps changes.
% -------------------------------------------------------------------------
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

% -------------------------------------------------------------------------
% E. Trace-processing parameters
% -------------------------------------------------------------------------
if ~exist('map_bin', 'var') || isempty(map_bin)
    map_bin = 4;
end
if ~exist('calcium_smoothing_window', 'var') || isempty(calcium_smoothing_window)
    calcium_smoothing_window = 40;
end
if ~exist('bleach_mode_voltage', 'var') || isempty(bleach_mode_voltage)
    bleach_mode_voltage = 'linear';   % 'linear' | 'highpass' | 'exp2'
end
if ~exist('bleach_mode_calcium', 'var') || isempty(bleach_mode_calcium)
    bleach_mode_calcium = 'exp2';     % 'linear' | 'highpass' | 'exp2'
end
if ~exist('run_background_removal', 'var') || isempty(run_background_removal)
    run_background_removal = false;
end
if ~exist('voltage_polarity', 'var') || isempty(voltage_polarity)
    voltage_polarity = -1;            % default display/analysis polarity for voltage traces
end
if ~exist('calcium_polarity', 'var') || isempty(calcium_polarity)
    calcium_polarity = 1;             % default display/analysis polarity for calcium traces
end
% Inter-camera ROI offset mode:
%   'none'            -> do not estimate a new offset; use reuse_offset or [0 0]
%   'manual_points'   -> manually click one matching point on the voltage
%                        and calcium average images
%   'matlab_register' -> use MATLAB built-in translation registration on
%                        the two channel average images
%
% Offset convention used everywhere below:
%   voltage_position = calcium_position + offset
% So both manual_points and matlab_register estimate the shift that moves
% the calcium channel into the voltage channel coordinate system.
%
% Backward compatibility:
%   correct_offset = true  -> 'manual_points'
%   correct_offset = false -> 'none'
if ~exist('correct_offset_mode', 'var') || isempty(correct_offset_mode)
    if exist('correct_offset', 'var') && ~isempty(correct_offset)
        correct_offset_mode = correct_offset;
    else
        correct_offset_mode = 'manual_points';
    end
end
correct_offset_mode = normalize_correct_offset_mode_dual(correct_offset_mode);
if ~exist('reuse_roi_file', 'var') || isempty(reuse_roi_file)
    reuse_roi_file = '';              % e.g. fullfile(save_path, '1_dual_roi_results.mat')
end
reuse_roi_resolution = resolve_dual_reuse_roi_source(reuse_results_path, reuse_roi_file);
reuse_roi_file = char(reuse_roi_resolution.effective_roi_file);
if strlength(reuse_roi_resolution.message) > 0
    fprintf('%s\n', reuse_roi_resolution.message);
end
if ~exist('reuse_offset', 'var') || isempty(reuse_offset)
    reuse_offset = [];
end
if ~exist('gpu', 'var') || isempty(gpu)
    gpu = true;
end

if ~exist('motion_cfg', 'var') || isempty(motion_cfg)
    motion_cfg = struct( ...
        'enabled', true, ...
        'use_saved_shift', false, ...
        'saved_shift_file', '', ...
        'highpass', true, ...
        'auto_reuse_previous_shift', true);
elseif ~isfield(motion_cfg, 'auto_reuse_previous_shift') || isempty(motion_cfg.auto_reuse_previous_shift)
    motion_cfg.auto_reuse_previous_shift = true;
end
motion_cfg.enabled = logical(run_motion_correction);
if ~exist('camera_source_override', 'var') || isempty(camera_source_override)
    camera_source_override = [];
end
if ~exist('stim_context_override', 'var') || isempty(stim_context_override)
    stim_context_override = [];
end

if numel(camera_cfg) ~= 2
    error('Dual_analysis3 currently expects exactly two cameras.');
end

if ~exist('save_path', 'var') || isempty(save_path)
    if strcmpi(string(analysis_mode), "analysis_only") && strlength(string(reuse_results_path)) > 0
        save_path = reuse_results_path;
    else
        if strcmpi(string(analysis_backend), "volpy_voltage_reanalysis")
            save_root_dir = 'Dual_analysis3_volpy_voltage';
        else
            save_root_dir = 'Dual_analysis3';
        end
        save_base_dir = resolve_dual_analysis_save_base_dir(cycle_path, camera_cfg, raw_dual_input_cfg);
        save_path = fullfile(save_base_dir, save_root_dir, analysis_run_name);
    end
end
mkdir(save_path);
analysis_only_mode = strcmpi(string(analysis_mode), "analysis_only");
volpy_reanalysis_mode = strcmpi(string(analysis_backend), "volpy_voltage_reanalysis");
% Parallel pool setup happens once, before any movie loading or map
% generation. Helper functions only reuse an existing pool and fall back to
% serial work when no pool is available.
parallel_pool_info = struct( ...
    'requested', logical(gpu && ~volpy_reanalysis_mode && ~analysis_only_mode), ...
    'ready', false, ...
    'num_workers', 0, ...
    'message', "", ...
    'initialized_before_movie_load', true);
if parallel_pool_info.requested
    print_section('Parallel Pool Setup');
    [~, parallel_pool_ready, parallel_pool_info] = initialize_dual_parallel_pool(parallel_pool_info);
    parallel_pool_info.ready = parallel_pool_ready;
else
    parallel_pool_info.message = "Parallel pool not requested for this mode.";
end

if volpy_reanalysis_mode
    bootstrap_volpy_voltage_reanalysis( ...
        cycle_path, save_path, volpy_source_results_path, ...
        logical(volpy_auto_run), logical(volpy_force_rerun), logical(volpy_use_existing), ...
        logical(volpy_flip_signal), ...
        voltage_polarity, calcium_polarity, calcium_smoothing_window);
end
if analysis_only_mode
    % This rerun mode reuses everything up through ROI selection and then
    % resumes the analysis sections that operate on saved ROI traces.
    [dual_info, voltage_results, calcium_results, dual_results, stim_results, stim_context, stim_windows] = ...
        load_saved_analysis_only_context(reuse_results_path);
    if ~strcmpi(char(string(save_path)), char(string(reuse_results_path)))
        [dual_info, voltage_results, calcium_results, dual_results, stim_results] = ...
            seed_analysis_only_output_folder( ...
            save_path, reuse_results_path, dual_info, voltage_results, calcium_results, dual_results, stim_results);
    end
end

if volpy_reanalysis_mode
    [dual_info, voltage_results, calcium_results, dual_results, stim_results, stim_context, stim_windows] = ...
        load_saved_analysis_only_context(save_path);
    t_voltage = build_time_axis_from_movie_info(voltage_results.movie_info);
    t_calcium = build_time_axis_from_movie_info(calcium_results.movie_info);
    freq_voltage = double(voltage_results.movie_info.frame_rate);
    freq_calcium = double(calcium_results.movie_info.frame_rate);
    nframes_voltage = infer_frame_count_from_channel_results(voltage_results);
    nframes_calcium = infer_frame_count_from_channel_results(calcium_results);
    dual_info_path = fullfile(save_path, 'dual_info.mat');
    voltage_results_path = fullfile(save_path, 'voltage_results.mat');
    calcium_results_path = fullfile(save_path, 'calcium_results.mat');
    dual_results_path = fullfile(save_path, 'dual_results.mat');
    stim_results_path = fullfile(save_path, 'stim_results.mat');
    run_dual_post_trace_sections( ...
        save_path, ...
        dual_info, dual_info_path, ...
        voltage_results, voltage_results_path, ...
        calcium_results, calcium_results_path, ...
        dual_results, dual_results_path, ...
        stim_results, stim_results_path, ...
        stim_context, stim_windows, ...
        t_voltage, t_calcium, ...
        freq_voltage, freq_calcium, ...
        nframes_voltage, nframes_calcium, ...
        voltage_polarity, calcium_polarity, calcium_smoothing_window, true);
    return;
end

if ~analysis_only_mode
%% Resolve Rebuilt Inputs
% Prefer rebuilt manifests so the analysis follows the acquisition layout
% recorded during acquisition instead of guessing from folder names alone.
print_section('Resolve Rebuilt Inputs');
[cycle_manifest, record_manifest, camera_source, input_layout_info] = resolve_dual_camera_sources( ...
    cycle_path, camera_cfg, camera_source_override, raw_dual_input_cfg);

voltage_idx = find(strcmpi(string({camera_cfg.role}), "voltage"), 1, 'first');
calcium_idx = find(strcmpi(string({camera_cfg.role}), "calcium"), 1, 'first');
if isempty(voltage_idx) || isempty(calcium_idx)
    error('camera_cfg must contain one voltage camera and one calcium camera.');
end

stim_context = resolve_stim_context( ...
    cycle_path, cycle_manifest, record_manifest, ...
    camera_cfg(voltage_idx).camera_index, camera_cfg(calcium_idx).camera_index, ...
    stim_context_override, input_layout_info);

dual_info = struct();
dual_info.analysis_name = 'Dual_analysis3';
dual_info.cycle_path = cycle_path;
dual_info.save_path = save_path;
dual_info.created_at = datetime("now");
dual_info.map_bin = map_bin;
dual_info.calcium_smoothing_window = calcium_smoothing_window;
dual_info.bleach_mode = struct( ...
    'voltage', bleach_mode_voltage, ...
    'calcium', bleach_mode_calcium);
dual_info.polarity = struct( ...
    'voltage', voltage_polarity, ...
    'calcium', calcium_polarity);
dual_info.correct_offset = ~strcmpi(string(correct_offset_mode), "none");
dual_info.correct_offset_mode = string(correct_offset_mode);
dual_info.run_background_removal = run_background_removal;
dual_info.camera_cfg = camera_cfg;
dual_info.camera_source = camera_source;
dual_info.input_layout = input_layout_info;
dual_info.input_override_used = ~isempty(camera_source_override) || ~isempty(stim_context_override);
dual_info.role_order = struct( ...
    'voltage_camera_index', camera_cfg(voltage_idx).camera_index, ...
    'calcium_camera_index', camera_cfg(calcium_idx).camera_index);
dual_info.manual_options = struct( ...
    'voltage_frame_rate', voltage_frame_rate, ...
    'calcium_frame_rate', calcium_frame_rate, ...
    'voltage_transpose_movie', logical(voltage_transpose_movie), ...
    'calcium_transpose_movie', logical(calcium_transpose_movie), ...
    'run_motion_correction', logical(run_motion_correction), ...
    'parallel_pool_info', parallel_pool_info, ...
    'correct_offset_mode', string(correct_offset_mode), ...
    'raw_dual_input_cfg', raw_dual_input_cfg, ...
    'reuse_roi_resolution', reuse_roi_resolution);
dual_info.stim_context = rmfield_if_exists(stim_context, {'logs', 'method_manifest'});
dual_info.cycle_manifest_found = ~isempty(cycle_manifest);
dual_info.record_manifest_found = ~isempty(record_manifest);
dual_info.cycle_manifest_path = fullfile(cycle_path, 'cycle_manifest.mat');
dual_info.record_manifest_path = fullfile(fileparts(cycle_path), 'record_manifest.mat');

dual_info_path = fullfile(save_path, 'dual_info.mat');
save(dual_info_path, 'dual_info');

voltage_results = struct('movie_info', struct(), 'trace_results', struct());
calcium_results = struct('movie_info', struct(), 'trace_results', struct());
dual_results = struct();
stim_results = struct();
stim_results.info = rmfield_if_exists(stim_context, {'logs', 'method_manifest'});

voltage_results_path = fullfile(save_path, 'voltage_results.mat');
calcium_results_path = fullfile(save_path, 'calcium_results.mat');
dual_results_path = fullfile(save_path, 'dual_results.mat');
stim_results_path = fullfile(save_path, 'stim_results.mat');
map_results_file = fullfile(save_path, '0_dual_sensitivity_map.mat');
save(voltage_results_path, 'voltage_results', '-v7.3');
save(calcium_results_path, 'calcium_results', '-v7.3');
save(dual_results_path, 'dual_results', '-v7.3');
save(stim_results_path, 'stim_results', '-v7.3');
fprintf('Initialized result containers: voltage, calcium, dual, stim.\n');
fprintf('Record mode detected: %s | Stim supported: %s\n', string(stim_context.recordmode), string(stim_context.supported));

%% Load Movies
% Load both camera movies first, then map them into voltage/calcium roles.
% This is the point where on-disk acquisitions become the in-memory movies
% used by all later processing.
print_section('Load Movies');
fprintf('Loading both camera movies...\n');
camera_data = repmat(struct(), 1, numel(camera_cfg));

for i = 1:numel(camera_cfg)
    [movie_3d, ncols_i, nrows_i, nframes_i, original_ncols_i, original_nrows_i] = load_camera_movie( ...
        camera_source(i).path, camera_cfg(i).transpose_before_analysis);

    camera_data(i).role = camera_cfg(i).role;
    camera_data(i).camera_index = camera_cfg(i).camera_index;
    camera_data(i).label = camera_source(i).label;
    camera_data(i).source_path = camera_source(i).path;
    camera_data(i).movie_3d = movie_3d;
    camera_data(i).ncols = ncols_i;
    camera_data(i).nrows = nrows_i;
    camera_data(i).nframes = nframes_i;
    camera_data(i).frame_rate = camera_cfg(i).frame_rate;
    camera_data(i).transpose_before_analysis = camera_cfg(i).transpose_before_analysis;
    camera_data(i).original_frame_size = [original_ncols_i, original_nrows_i];
end
clear movie_3d;

for i = 1:numel(camera_cfg)
    fprintf('Camera %d | role=%s | label=%s | source=%s | transpose=%d | size=%dx%d | frames=%d | fps=%g\n', ...
        camera_data(i).camera_index, string(camera_data(i).role), string(camera_data(i).label), ...
        camera_data(i).source_path, camera_data(i).transpose_before_analysis, ...
        camera_data(i).ncols, camera_data(i).nrows, camera_data(i).nframes, camera_data(i).frame_rate);
end

[camera_data(voltage_idx), camera_data(calcium_idx), geometry_info] = ...
    match_camera_geometry(camera_data(voltage_idx), camera_data(calcium_idx));
% Force both channels onto one common analysis geometry before any shared
% ROI work. Later sections can then treat voltage and calcium ROIs as
% spatially matched by construction.
fprintf('Geometry match action: %s | common size=%dx%d\n', string(geometry_info.action), geometry_info.common_size(1), geometry_info.common_size(2));

dual_info.geometry = geometry_info;
save(dual_info_path, 'dual_info');

voltage_results.movie_info = make_movie_info(camera_data(voltage_idx));
calcium_results.movie_info = make_movie_info(camera_data(calcium_idx));
save(voltage_results_path, 'voltage_results', '-v7.3');
save(calcium_results_path, 'calcium_results', '-v7.3');

movie_voltage_3d = camera_data(voltage_idx).movie_3d;
ncols = camera_data(voltage_idx).ncols;
nrows = camera_data(voltage_idx).nrows;
nframes_voltage = camera_data(voltage_idx).nframes;
freq_voltage = camera_data(voltage_idx).frame_rate;
t_voltage = (1:nframes_voltage)' / freq_voltage;

movie_calcium_3d = camera_data(calcium_idx).movie_3d;
nframes_calcium = camera_data(calcium_idx).nframes;
freq_calcium = camera_data(calcium_idx).frame_rate;
t_calcium = (1:nframes_calcium)' / freq_calcium;
camera_data = rmfield_if_exists(camera_data, {'movie_3d'});
clear camera_data;

%% Copy Analysis Code
% Save the exact code snapshot used for this run next to the results. This
% makes old result folders easier to interpret after the code evolves.
print_section('Copy Analysis Code');
code_path = fullfile(save_path, 'Code');
mkdir(code_path);
currentScript = '';
if exist('dual_script_path', 'var') && strlength(string(dual_script_path)) > 0 ...
        && isfile(char(string(dual_script_path)))
    currentScript = char(string(dual_script_path));
end
if isempty(currentScript)
    currentScript = mfilename('fullpath');
end
if isempty(currentScript) || ~isfile(currentScript)
    currentScript = which('Dual_analysis3.m');
end
copy_analysis_code(currentScript, code_path);

%% Motion Correction
% Motion is estimated once on voltage and then applied to calcium.
% The processing intent is to keep both channels spatially locked. If the
% two channels were corrected independently, one biological ROI could end
% up drifting to different places in voltage and calcium.
print_section('Motion Correction');
fprintf('Motion correction enabled: %d\n', logical(run_motion_correction));
if run_motion_correction
    fprintf('Applying shared motion correction from voltage to calcium...\n');
    fprintf('Motion model: rigid NoRMCorre estimated on voltage, then apply_shifts to calcium.\n');
    fprintf('Motion QC outputs: save downsampled corrected TIFFs for both channels and save shared motion metrics on the voltage reference movie.\n');
else
    fprintf('Shared motion correction skipped. Downstream ROI and trace analysis will use the geometry-matched but otherwise raw movies.\n');
end
[movie_voltage_3d, movie_calcium_3d, voltage_motion_info, calcium_motion_info] = run_shared_motion_correction( ...
    movie_voltage_3d, movie_calcium_3d, cycle_path, save_path, motion_cfg);

[ncols_v, nrows_v, nframes_voltage] = size(movie_voltage_3d);
[ncols_c, nrows_c, nframes_calcium] = size(movie_calcium_3d);
if ncols_v ~= ncols_c || nrows_v ~= nrows_c
    error('Camera sizes diverged after motion correction. Check geometry matching and transpose settings.');
end

movie_voltage = reshape(movie_voltage_3d, ncols_v * nrows_v, nframes_voltage);
movie_calcium = reshape(movie_calcium_3d, ncols_c * nrows_c, nframes_calcium);
% From here on, movies are stored in [pixels x frames] form because the
% ROI-processing functions work primarily on trace-oriented matrices.
clear movie_voltage_3d movie_calcium_3d;
ncols = ncols_v;
nrows = nrows_v;
t_voltage = (1:nframes_voltage)' / freq_voltage;
t_calcium = (1:nframes_calcium)' / freq_calcium;

voltage_results.movie_info.motion = voltage_motion_info;
voltage_results.movie_info.analysis_frame_size = [ncols, nrows];
voltage_results.movie_info.frame_count = nframes_voltage;
voltage_results.movie_info.updated_at = datetime("now");

calcium_results.movie_info.motion = calcium_motion_info;
calcium_results.movie_info.analysis_frame_size = [ncols, nrows];
calcium_results.movie_info.frame_count = nframes_calcium;
calcium_results.movie_info.updated_at = datetime("now");

save(voltage_results_path, 'voltage_results', '-v7.3');
save(calcium_results_path, 'calcium_results', '-v7.3');
if voltage_motion_info.applied
    fprintf('Motion complete | shared shift file: %s\n', voltage_motion_info.shift_file);
else
    fprintf('Motion correction result: skipped by run_motion_correction=false.\n');
end
if isfield(voltage_motion_info, 'downsampled_tif_file') && strlength(string(voltage_motion_info.downsampled_tif_file)) > 0
    fprintf('Saved voltage downsampled corrected TIFF: %s\n', string(voltage_motion_info.downsampled_tif_file));
end
if isfield(calcium_motion_info, 'downsampled_tif_file') && strlength(string(calcium_motion_info.downsampled_tif_file)) > 0
    fprintf('Saved calcium downsampled corrected TIFF: %s\n', string(calcium_motion_info.downsampled_tif_file));
end
if isfield(voltage_motion_info, 'metrics_png') && strlength(string(voltage_motion_info.metrics_png)) > 0
    fprintf('Saved shared motion metrics plot: %s\n', string(voltage_motion_info.metrics_png));
end

%% Channel Registration
% In matlab_register mode, estimate the calcium-to-voltage translation
% offset, but do not directly move the movie data here. The returned
% offset is used later to map voltage ROI masks into calcium coordinates.
print_section('Channel Registration');
[channel_registration_info, roi_offset_mode_for_selection] = ...
    estimate_dual_channel_registration_offset( ...
        movie_voltage, movie_calcium, nrows, ncols, correct_offset_mode, save_path);
voltage_results.movie_info.channel_registration = channel_registration_info;
calcium_results.movie_info.channel_registration = channel_registration_info;
save(voltage_results_path, 'voltage_results', '-v7.3');
save(calcium_results_path, 'calcium_results', '-v7.3');
fprintf('ROI offset mode used after channel registration: %s\n', string(roi_offset_mode_for_selection));
if isempty(reuse_offset) && isfield(channel_registration_info, 'offset_xy') ...
        && numel(channel_registration_info.offset_xy) == 2
    reuse_offset = double(channel_registration_info.offset_xy);
end

%% Create Sensitivity Maps
% These maps are quick summary images used mainly to guide ROI selection.
% They are visual aids rather than final quantitative outputs.
print_section('Create Sensitivity Maps');
fprintf('Creating voltage/calcium maps...\n');
fprintf('Map method: create_map(movie, nrows, ncols, map_bin), map_bin=%d\n', map_bin);
map_voltage = create_map(movie_voltage, nrows, ncols, map_bin);
map_calcium = create_map(movie_calcium, nrows, ncols, map_bin, 'calcium');

% 暂且用这几行抵消一下相机第一帧第一行65535的bug
if exist('VolMap_Cut', 'var') 
    map_voltage(1,:) = map_voltage(5,:);
    map_voltage(2,:) = map_voltage(5,:);
    map_voltage(3,:) = map_voltage(5,:);
    map_voltage(4,:) = map_voltage(5,:);
end

map_fig = figure('Color', 'w');
subplot(1, 2, 1);
imagesc(map_voltage);
axis image;
title('Voltage Sensitivity Map');
colorbar;
subplot(1, 2, 2);
imagesc(map_calcium);
axis image;
title('Calcium Sensitivity Map');
colorbar;

save(map_results_file, 'map_voltage', 'map_calcium', 'map_bin');
save_figure_bundle(map_fig, fullfile(save_path, '0_dual_sensitivity_map.fig'), fullfile(save_path, '0_dual_sensitivity_map.png'));
close(map_fig);

dual_results.maps = struct( ...
    'data', struct('map_voltage', map_voltage, 'map_calcium', map_calcium), ...
    'info', struct('map_file', map_results_file, 'map_bin', map_bin, 'created_at', datetime("now")));
save(dual_results_path, 'dual_results', '-v7.3');
fprintf('Saved maps to: %s\n', map_results_file);

%% Load Or Create Dual ROI
% Dual ROI selection always starts from the voltage movie and projects to
% the calcium movie through an explicit offset.
print_section('Load Or Create Dual ROI');
fprintf('Selecting dual ROIs...\n');
mask_voltage = [];
mask_calcium = [];
offset = reuse_offset;
roi_map_source = "unknown";

if ~exist('movie_voltage', 'var')
    if exist('movie_voltage_3d', 'var')
        [ncols_v, nrows_v, nframes_voltage] = size(movie_voltage_3d);
        movie_voltage = reshape(movie_voltage_3d, ncols_v * nrows_v, nframes_voltage);
        ncols = ncols_v;
        nrows = nrows_v;
        t_voltage = (1:nframes_voltage)' / freq_voltage;
        fprintf('ROI movie source | voltage: reshaped current 3D movie without requiring Motion Correction section\n');
    else
        error(['ROI selection needs voltage movie data in memory. ' ...
            'Run Load Movies first, then you may skip Motion Correction/Create Sensitivity Maps.']);
    end
end

if ~exist('movie_calcium', 'var')
    if exist('movie_calcium_3d', 'var')
        [ncols_c, nrows_c, nframes_calcium] = size(movie_calcium_3d);
        movie_calcium = reshape(movie_calcium_3d, ncols_c * nrows_c, nframes_calcium);
        t_calcium = (1:nframes_calcium)' / freq_calcium;
        fprintf('ROI movie source | calcium: reshaped current 3D movie without requiring Motion Correction section\n');
    else
        error(['ROI selection needs calcium movie data in memory. ' ...
            'Run Load Movies first, then you may skip Motion Correction/Create Sensitivity Maps.']);
    end
end

if exist('map_voltage', 'var') && exist('map_calcium', 'var')
    roi_map_source = "workspace_maps";
elseif isfile(map_results_file)
    map_cache = load(map_results_file, 'map_voltage', 'map_calcium');
    if isfield(map_cache, 'map_voltage') && isfield(map_cache, 'map_calcium')
        map_voltage = map_cache.map_voltage;
        map_calcium = map_cache.map_calcium;
        roi_map_source = "saved_map_results";
    end
end

if ~(exist('map_voltage', 'var') && exist('map_calcium', 'var'))
    fprintf('Create Sensitivity Maps section was skipped. ROI selection will proceed without maps, following the original select_ROI strategy.\n');
    map_voltage = [];
    map_calcium = [];
    roi_map_source = "no_map";
end

fprintf('ROI prerequisites | motion applied: voltage=%d calcium=%d | map source=%s\n', ...
    voltage_results.movie_info.motion.applied, calcium_results.movie_info.motion.applied, roi_map_source);

if ~isempty(reuse_roi_file) && isfile(reuse_roi_file)
    % ROI reuse is mainly for iteration: it lets later sections be rerun
    % without forcing the user to redraw the same biological regions.
    roi_cache = load(reuse_roi_file);
    offset_for_reuse_mask = offset;
    if isempty(offset_for_reuse_mask) && isfield(roi_cache, 'offset')
        offset_for_reuse_mask = roi_cache.offset;
    end
    [mask_voltage, mask_calcium, reuse_roi_note] = resolve_reuse_dual_masks( ...
        roi_cache, offset_for_reuse_mask, roi_offset_mode_for_selection);
    if strlength(reuse_roi_note) > 0
        fprintf('%s\n', reuse_roi_note);
    end
    if isempty(offset) && isfield(roi_cache, 'offset')
        offset = roi_cache.offset;
    end
end

if isempty(offset)
    offset = [0, 0];
elseif numel(offset) ~= 2
    error('ROI offset must be empty or a 1x2 vector [x_offset, y_offset].');
else
    offset = reshape(double(offset), 1, 2);
end

% Interactive ROI drawing cannot succeed in a headless MATLAB session.
% Offset picking for manual_points is handled earlier in Channel
% Registration, so the ROI section only needs to guard the polygon-drawing
% path here.
requires_interactive_roi_session = isempty(mask_voltage) && isempty(mask_calcium);
if ~usejava('desktop') && requires_interactive_roi_session
    error(['Interactive ROI selection requires a MATLAB Desktop session. ' ...
        'Please rerun Dual_analysis3 or Dual_analysis3_rec inside an existing MATLAB GUI window, ' ...
        'or set reuse_roi_file to an already saved ROI result file.']);
end

roi_figures_before = findall(groot, 'Type', 'figure');
[rois, traces_voltage_raw, traces_calcium_raw, offset] = select_ROI_dual( ...
    movie_voltage, movie_calcium, nrows, ncols, roi_offset_mode_for_selection, ...
    map_voltage, map_calcium, mask_voltage, mask_calcium, offset);
% This is the main transition from movie space to trace space. After this
% point, most analysis steps work on ROI-by-time matrices rather than raw
% image stacks.
roi_figures_after = findall(groot, 'Type', 'figure');
for fig_idx = 1:numel(roi_figures_after)
    if ~any(roi_figures_after(fig_idx) == roi_figures_before) && isgraphics(roi_figures_after(fig_idx), 'figure')
        close(roi_figures_after(fig_idx));
    end
end

if isempty(offset)
    offset = [0, 0];
elseif numel(offset) ~= 2
    error('select_ROI_dual returned an invalid offset. Expected a 1x2 vector.');
else
    offset = reshape(double(offset), 1, 2);
end

nrois = size(traces_voltage_raw, 2);
roi_results_file = fullfile(save_path, '1_dual_roi_results.mat');
% This record is the "interface summary" for ROI selection.
% Natural-language meaning:
%   this section turns the paired movies into paired ROI masks plus the
%   first trace stage ("raw") for both channels.
%
% Interface meaning:
%   input      -> geometry, offset settings, and reuse choices
%   parameters -> which selection function / channel roles were used
%   output     -> ROI masks, offsets, raw traces
%
% Example:
%   roi_record.output.traces_voltage_raw
%   roi_record.output.rois
%   roi_record.excluded_inputs.movie_voltage = 'not saved'
roi_record = build_section_record( ...
    'Load Or Create Dual ROI', ...
    'Select paired voltage/calcium ROIs and extract the raw trace matrices for both channels.', ...
    struct( ...
        'nrows', nrows, ...
        'ncols', ncols, ...
        'correct_offset_mode', string(correct_offset_mode), ...
        'reuse_roi_file', string(reuse_roi_file), ...
        'initial_offset_xy', offset, ...
        'roi_map_source', string(roi_map_source)), ...
    struct( ...
        'selection_function', 'select_ROI_dual', ...
        'reference_role', 'voltage', ...
        'moving_role', 'calcium'), ...
    struct( ...
        'rois', rois, ...
        'offset_xy', offset, ...
        'nrois', nrois, ...
        'traces_voltage_raw', traces_voltage_raw, ...
        'traces_calcium_raw', traces_calcium_raw), ...
    struct( ...
        'movie_voltage', 'not saved', ...
        'movie_calcium', 'not saved', ...
        'map_voltage', 'not saved', ...
        'map_calcium', 'not saved'), ...
    'To rerun ROI selection, reload the movies/maps externally, then call select_ROI_dual with the saved geometry and offset settings from this record.');
save(roi_results_file, 'rois', 'offset', 'nrows', 'ncols', 'nrois', ...
    'traces_voltage_raw', 'traces_calcium_raw', 'roi_record');

dual_results.registration = struct( ...
    'data', struct( ...
        'offset_xy', offset, ...
        'nrois', nrois, ...
        'voltage_mask', rois.bwmask, ...
        'calcium_mask', rois.bwmask_ca), ...
    'info', struct( ...
        'roi_file', roi_results_file, ...
        'reference_role', 'voltage', ...
        'moving_role', 'calcium', ...
        'map_source', roi_map_source, ...
        'offset_mode', string(correct_offset_mode), ...
        'created_at', datetime("now")));
save(dual_results_path, 'dual_results', '-v7.3');

voltage_results = store_trace_stage( ...
    voltage_results, 'raw', traces_voltage_raw, {}, roi_results_file, ...
    voltage_results.movie_info, 'select_ROI_dual', struct('role', 'voltage', 'offset_xy', offset));
calcium_results = store_trace_stage( ...
    calcium_results, 'raw', traces_calcium_raw, {}, roi_results_file, ...
    calcium_results.movie_info, 'select_ROI_dual', struct('role', 'calcium', 'offset_xy', offset));
save(voltage_results_path, 'voltage_results', '-v7.3');
save(calcium_results_path, 'calcium_results', '-v7.3');
raw_trace_fig = fullfile(save_path, '1_dual_raw_trace.fig');
raw_trace_png = fullfile(save_path, '1_dual_raw_trace.png');
map_roi_fig = fullfile(save_path, '0_dual_sensitivity_map_with_roi.fig');
map_roi_png = fullfile(save_path, '0_dual_sensitivity_map_with_roi.png');
merged_average_tif = fullfile(save_path, '0_dual_average_color_merge.tif');
merged_average_png = fullfile(save_path, '0_dual_average_color_merge.png');
merged_average_roi_tif = fullfile(save_path, '0_dual_average_color_merge_with_roi.tif');
merged_average_roi_png = fullfile(save_path, '0_dual_average_color_merge_with_roi.png');
plot_dual_raw_traces( ...
    traces_voltage_raw, traces_calcium_raw, ...
    t_voltage, t_calcium);
raw_trace_handle = gcf;
save_figure_bundle(raw_trace_handle, raw_trace_fig, raw_trace_png);
close(raw_trace_handle);
[map_roi_saved, map_roi_info] = save_dual_map_with_roi( ...
    map_voltage, map_calcium, rois, map_roi_fig, map_roi_png);
[~, merged_average_info] = save_dual_average_color_merge( ...
    movie_voltage, movie_calcium, nrows, ncols, offset, ...
    merged_average_tif, merged_average_png, ...
    merged_average_roi_tif, merged_average_roi_png, rois);
dual_results.visualizations.raw = struct( ...
    'data', struct(), ...
    'info', struct( ...
        'fig_file', raw_trace_fig, ...
        'png_file', raw_trace_png, ...
        'created_at', datetime("now")));
dual_results.visualizations.average_color_merge = struct( ...
    'data', struct(), ...
    'info', merged_average_info);
if map_roi_saved
    dual_results.visualizations.map_with_roi = struct( ...
        'data', struct(), ...
        'info', map_roi_info);
end
save(dual_results_path, 'dual_results', '-v7.3');
fprintf('Raw trace overview saved to: %s\n', raw_trace_png);
if map_roi_saved
    fprintf('Dual sensitivity map with ROI saved to: %s\n', map_roi_png);
end
fprintf('Dual average color merge saved to: %s\n', merged_average_tif);
fprintf('Dual average color merge with ROI saved to: %s\n', merged_average_roi_tif);
fprintf('ROI complete | nrois=%d | offset=[%.3f %.3f] | map source=%s | roi file=%s\n', ...
    nrois, offset(1), offset(2), roi_map_source, roi_results_file);
else
    dual_info_path = fullfile(save_path, 'dual_info.mat');
    voltage_results_path = fullfile(save_path, 'voltage_results.mat');
    calcium_results_path = fullfile(save_path, 'calcium_results.mat');
    dual_results_path = fullfile(save_path, 'dual_results.mat');
    stim_results_path = fullfile(save_path, 'stim_results.mat');
    map_results_file = fullfile(save_path, '0_dual_sensitivity_map.mat');

    cycle_path = char(string(dual_info.cycle_path));
    freq_voltage = resolve_saved_channel_frame_rate(dual_info, "voltage");
    freq_calcium = resolve_saved_channel_frame_rate(dual_info, "calcium");
    nframes_voltage = resolve_saved_channel_frame_count(voltage_results.movie_info);
    nframes_calcium = resolve_saved_channel_frame_count(calcium_results.movie_info);
    [ncols, nrows] = resolve_saved_analysis_frame_size(voltage_results.movie_info, calcium_results.movie_info);
    t_voltage = build_time_axis_from_movie_info(voltage_results.movie_info);
    t_calcium = build_time_axis_from_movie_info(calcium_results.movie_info);

    fprintf('Analysis-only mode: reusing saved ROI/channel results from %s\n', save_path);
    fprintf('Skipped sections: Resolve Inputs, Load Movies, Copy Analysis Code, Motion Correction, Create Sensitivity Maps, Load Or Create Dual ROI\n');
    fprintf(['Analysis-only parameter source: downstream background/bleach/metric/stim/time-frequency settings ' ...
        'still come from the current script variables, while saved ROI/channel/stim files are reused as inputs.\n']);
end

%% Background Removal
% Background removal is optional at the workflow level. If it is skipped,
% later sections automatically fall back to the raw trace stage.
%
% Natural-language meaning:
%   estimate a local/background signal around each ROI and subtract it,
%   so the trace is less affected by slow whole-field intensity changes.
%
% Interface meaning:
%   input      -> raw traces, ROI definitions, frame geometry, movie data
%   parameters -> remove_background settings for each channel
%   output     -> background fits, masks, and bg_removed traces
%
% Note:
%   movie data are required to recompute this section, but are not saved in
%   the background result file. The saved background_record tells you this
%   explicitly so old results remain understandable without duplicating the
%   movie on disk.
print_section('Background Removal');
fprintf('Background removal enabled: %d\n', run_background_removal);
[voltage_results, calcium_results] = load_channel_results( ...
    voltage_results_path, calcium_results_path, voltage_results, calcium_results);
[dual_results, rois, traces_voltage_raw, traces_calcium_raw, roi_results_file, nrois] = ...
    load_dual_roi_context(dual_results_path, dual_results);

rois_voltage = struct('bwmask', rois.bwmask, 'boundary', {rois.boundary}, 'position', {rois.position});
rois_calcium = struct('bwmask', rois.bwmask_ca, 'boundary', {rois.boundary_ca}, 'position', {rois.position_ca});
background_results_file = fullfile(save_path, '1_dual_background_results.mat');

if run_background_removal
    if ~exist('movie_voltage', 'var') || ~exist('movie_calcium', 'var')
        if has_trace_stage(voltage_results, 'bg_removed') && has_trace_stage(calcium_results, 'bg_removed')
            fprintf('Background stage already exists. Reusing saved bg_removed traces without recomputation.\n');
            fprintf('Current bg_removed parent stage: voltage=%s | calcium=%s\n', ...
                strjoin(string(voltage_results.trace_results.bg_removed.info.parent_results), ','), ...
                strjoin(string(calcium_results.trace_results.bg_removed.info.parent_results), ','));
        else
            fprintf('Background removal skipped because movie data are not in memory. Downstream sections will use raw traces.\n');
        end
    else
        [background_voltage, background_fit_voltage, ~, traces_voltage_bg, background_mask_voltage] = ...
            remove_background(movie_voltage, ncols, nrows, rois_voltage, freq_voltage, 1);
        [background_calcium, background_fit_calcium, ~, traces_calcium_bg, background_mask_calcium] = ...
            remove_background(movie_calcium, ncols, nrows, rois_calcium, freq_calcium, 1);

        % Save enough non-movie inputs to explain and reproduce the background
        % calculation later. This is especially useful when reopening only the
        % section result file instead of the whole workspace.
        background_record = build_section_record( ...
            'Background Removal', ...
            'Estimate and subtract ROI-matched background signals for both channels.', ...
            struct( ...
                'voltage_parent_stage', "raw", ...
                'calcium_parent_stage', "raw", ...
                'traces_voltage_raw', traces_voltage_raw, ...
                'traces_calcium_raw', traces_calcium_raw, ...
                'rois_voltage', rois_voltage, ...
                'rois_calcium', rois_calcium, ...
                'nrows', nrows, ...
                'ncols', ncols), ...
            struct( ...
                'voltage', struct('freq', freq_voltage, 'bin', 1), ...
                'calcium', struct('freq', freq_calcium, 'bin', 1), ...
                'function_name', 'remove_background'), ...
            struct( ...
                'background_voltage', background_voltage, ...
                'background_fit_voltage', background_fit_voltage, ...
                'background_mask_voltage', background_mask_voltage, ...
                'traces_voltage_bg', traces_voltage_bg, ...
                'background_calcium', background_calcium, ...
                'background_fit_calcium', background_fit_calcium, ...
                'background_mask_calcium', background_mask_calcium, ...
                'traces_calcium_bg', traces_calcium_bg), ...
            struct( ...
                'movie_voltage', 'not saved', ...
                'movie_calcium', 'not saved'), ...
            'To rerun this section, reload the movies externally and call remove_background with the saved ROI definitions and parameters in this record.');
        save(background_results_file, ...
            'background_voltage', 'background_fit_voltage', 'background_mask_voltage', 'traces_voltage_bg', ...
            'background_calcium', 'background_fit_calcium', 'background_mask_calcium', 'traces_calcium_bg', ...
            'background_record');

        voltage_results = store_trace_stage( ...
            voltage_results, 'bg_removed', traces_voltage_bg, {'raw'}, roi_results_file, ...
            voltage_results.movie_info, 'remove_background', struct('freq', freq_voltage, 'bin', 1));
        calcium_results = store_trace_stage( ...
            calcium_results, 'bg_removed', traces_calcium_bg, {'raw'}, roi_results_file, ...
            calcium_results.movie_info, 'remove_background', struct('freq', freq_calcium, 'bin', 1));
        save(voltage_results_path, 'voltage_results', '-v7.3');
        save(calcium_results_path, 'calcium_results', '-v7.3');
        fprintf('Background removal formula: bg_removed = remove_background(movie, roi)\n');
        fprintf('Background parameters | voltage freq=%g, calcium freq=%g, bin=1\n', freq_voltage, freq_calcium);
    end
else
    fprintf('Background removal skipped. Downstream sections will use raw traces as bleach input.\n');
end

if run_background_removal && isfile(background_results_file)
    % The background summary is a diagnostic check: it helps verify that
    % the chosen background region and fitted background trace look
    % sensible before trusting downstream bg_removed traces.
    bg_cache = load(background_results_file);
    t_voltage_bg = build_time_axis_from_movie_info(voltage_results.movie_info);
    t_calcium_bg = build_time_axis_from_movie_info(calcium_results.movie_info);
    mean_voltage_image = [];
    mean_calcium_image = [];
    if exist('movie_voltage', 'var')
        mean_voltage_image = mean(reshape(movie_voltage, ncols, nrows, []), 3);
    end
    if exist('movie_calcium', 'var')
        mean_calcium_image = mean(reshape(movie_calcium, ncols, nrows, []), 3);
    end
    [bg_fig, bg_png] = plot_dual_background_summary( ...
        rois, mean_voltage_image, mean_calcium_image, ...
        bg_cache.background_mask_voltage, bg_cache.background_mask_calcium, ...
        bg_cache.background_fit_voltage, bg_cache.background_fit_calcium, ...
        bg_cache.traces_voltage_bg, bg_cache.traces_calcium_bg, ...
        t_voltage_bg, t_calcium_bg, save_path);
    dual_results.visualizations.background = struct( ...
        'data', struct(), ...
        'info', struct( ...
            'fig_file', bg_fig, ...
            'png_file', bg_png, ...
            'created_at', datetime("now")));
    save(dual_results_path, 'dual_results', '-v7.3');
    fprintf('Background summary saved to: %s\n', bg_png);
else
    fprintf('Background summary plot skipped.\n');
end

%% Bleaching Removal
% Natural-language meaning:
%   remove the slow baseline drift caused by bleaching and save the fitted
%   baseline itself. The residual trace after this step becomes the common
%   starting point for noise estimation, sensitivity, and SNR.
%
% Interface meaning:
%   input      -> selected parent trace stage for each channel
%   parameters -> bleach mode, fitted parameters, and time axis
%   output     -> bleach_removed trace and fitted baseline
%
% Example:
%   bleach_record.input.voltage_parent_stage   -> "bg_removed"
%   bleach_record.parameters.calcium_mode      -> "exp2"
%   bleach_record.output.baseline_calcium      -> fitted bleach baseline
print_section('Bleaching Removal');
fprintf('Removing bleaching for both channels...\n');
[voltage_results, calcium_results] = load_channel_results( ...
    voltage_results_path, calcium_results_path, voltage_results, calcium_results);
[dual_results, ~, ~, ~, roi_results_file] = load_dual_roi_context(dual_results_path, dual_results);
[traces_voltage_input, voltage_parent_stage] = resolve_preferred_trace_stage(voltage_results, {'bg_removed', 'raw'});
[traces_calcium_input, calcium_parent_stage] = resolve_preferred_trace_stage(calcium_results, {'bg_removed', 'raw'});
freq_voltage_bleach = voltage_results.movie_info.frame_rate;
freq_calcium_bleach = calcium_results.movie_info.frame_rate;
t_voltage_bleach = build_time_axis_from_movie_info(voltage_results.movie_info);
t_calcium_bleach = build_time_axis_from_movie_info(calcium_results.movie_info);
fprintf('Bleach input stage | voltage=%s | calcium=%s\n', voltage_parent_stage, calcium_parent_stage);
fprintf('Bleach mode | voltage=%s | calcium=%s\n', bleach_mode_voltage, bleach_mode_calcium);
fprintf('Bleach formula | voltage=%s\n', describe_bleach_formula(bleach_mode_voltage));
fprintf('Bleach formula | calcium=%s\n', describe_bleach_formula(bleach_mode_calcium));
[traces_voltage_bleach, baseline_voltage, bleach_params_voltage] = ...
    run_bleach_removal(traces_voltage_input, freq_voltage_bleach, t_voltage_bleach, bleach_mode_voltage);
[traces_calcium_bleach, baseline_calcium, bleach_params_calcium] = ...
    run_bleach_removal(traces_calcium_input, freq_calcium_bleach, t_calcium_bleach, bleach_mode_calcium);
[bleach_fig, bleach_png] = plot_dual_bleach_overview( ...
    traces_voltage_bleach, traces_voltage_input, baseline_voltage, t_voltage_bleach, bleach_mode_voltage, ...
    traces_calcium_bleach, traces_calcium_input, baseline_calcium, t_calcium_bleach, bleach_mode_calcium, ...
    save_path);

bleach_results_file = fullfile(save_path, '2_dual_bleach_results.mat');
bleach_record = build_section_record( ...
    'Bleaching Removal', ...
    'Remove slow bleaching trends from the selected input trace stage of each channel and save the fitted baseline.', ...
    struct( ...
        'voltage_parent_stage', string(voltage_parent_stage), ...
        'calcium_parent_stage', string(calcium_parent_stage), ...
        'traces_voltage_input', traces_voltage_input, ...
        'traces_calcium_input', traces_calcium_input, ...
        't_voltage', t_voltage_bleach, ...
        't_calcium', t_calcium_bleach, ...
        'freq_voltage', freq_voltage_bleach, ...
        'freq_calcium', freq_calcium_bleach), ...
    struct( ...
        'voltage_mode', bleach_mode_voltage, ...
        'calcium_mode', bleach_mode_calcium, ...
        'voltage_parameters', bleach_params_voltage, ...
        'calcium_parameters', bleach_params_calcium, ...
        'voltage_formula', describe_bleach_formula(bleach_mode_voltage), ...
        'calcium_formula', describe_bleach_formula(bleach_mode_calcium)), ...
    struct( ...
        'traces_voltage_bleach', traces_voltage_bleach, ...
        'baseline_voltage', baseline_voltage, ...
        'traces_calcium_bleach', traces_calcium_bleach, ...
        'baseline_calcium', baseline_calcium), ...
    struct(), ...
    'To rerun this section, load the saved input traces and parameters from this record, then call run_bleach_removal for each channel.');
save(bleach_results_file, ...
    'traces_voltage_bleach', 'baseline_voltage', 'traces_calcium_bleach', 'baseline_calcium', ...
    'bleach_mode_voltage', 'bleach_mode_calcium', 'bleach_params_voltage', 'bleach_params_calcium', ...
    'bleach_record');

voltage_results = store_trace_stage( ...
    voltage_results, 'bleach_removed', traces_voltage_bleach, {voltage_parent_stage}, roi_results_file, ...
    voltage_results.movie_info, bleach_mode_voltage, bleach_params_voltage);
voltage_results = store_trace_stage( ...
    voltage_results, 'baseline', baseline_voltage, {voltage_parent_stage}, roi_results_file, ...
    voltage_results.movie_info, bleach_mode_voltage, bleach_params_voltage);

calcium_results = store_trace_stage( ...
    calcium_results, 'bleach_removed', traces_calcium_bleach, {calcium_parent_stage}, roi_results_file, ...
    calcium_results.movie_info, bleach_mode_calcium, bleach_params_calcium);
calcium_results = store_trace_stage( ...
    calcium_results, 'baseline', baseline_calcium, {calcium_parent_stage}, roi_results_file, ...
    calcium_results.movie_info, bleach_mode_calcium, bleach_params_calcium);
save(voltage_results_path, 'voltage_results', '-v7.3');
save(calcium_results_path, 'calcium_results', '-v7.3');
dual_results.visualizations.bleach = struct( ...
    'data', struct(), ...
    'info', struct( ...
        'fig_file', bleach_fig, ...
        'png_file', bleach_png, ...
        'voltage_mode', bleach_mode_voltage, ...
        'calcium_mode', bleach_mode_calcium, ...
        'created_at', datetime("now")));
save(dual_results_path, 'dual_results', '-v7.3');
fprintf('Bleach overview saved to: %s\n', bleach_png);
fprintf('Bleach output stages saved: bleach_removed, baseline\n');

%% Noise Estimation, Sensitivity, And SNR
% Natural-language meaning:
%   first build a noise reference for each channel, then convert the
%   bleach-corrected traces into sensitivity and SNR traces that can be
%   compared more meaningfully across ROIs and across channels. For
%   calcium, this section also includes the final smoothing step that later
%   plots and stimulus analyses will preferentially use.
%
% Interface meaning:
%   input      -> bleach_removed traces, fitted baselines, calcium raw
%   parameters -> noise-reference method, SNR rule, calcium smoothing rule
%   output     -> noise_reference, noise, sensitivity, snr, and the saved
%                 calcium final stages raw_smoothed / sensitivity_smoothed /
%                 snr_smoothed
%
% Example:
%   metric_record.input.traces_voltage_bleach
%   metric_record.parameters.calcium_smoothing_parameters.window
%   metric_record.output.calcium_snr_smoothed
print_section('Noise Estimation, Sensitivity, And SNR');
fprintf('Estimating noise references, sensitivity, and SNR...\n');
[voltage_results, calcium_results] = load_channel_results( ...
    voltage_results_path, calcium_results_path, voltage_results, calcium_results);
[~, ~, ~, ~, roi_results_file] = load_dual_roi_context(dual_results_path, dual_results);
traces_voltage_bleach = fetch_trace_stage(voltage_results, 'bleach_removed');
baseline_voltage = fetch_trace_stage(voltage_results, 'baseline');
traces_calcium_bleach = fetch_trace_stage(calcium_results, 'bleach_removed');
baseline_calcium = fetch_trace_stage(calcium_results, 'baseline');
[voltage_noise_reference, voltage_noise, voltage_sensitivity, voltage_snr, voltage_metric_info] = ...
    compute_voltage_metrics(traces_voltage_bleach, baseline_voltage);
[calcium_noise_reference, calcium_noise, calcium_sensitivity, calcium_snr, calcium_metric_info] = ...
    compute_calcium_metrics(traces_calcium_bleach, baseline_calcium, freq_calcium);
fprintf('Voltage noise reference: %s | formula noise = bleach_removed - noise_reference | sensitivity = bleach_removed ./ baseline | snr = bleach_removed ./ std(noise)\n', ...
    voltage_metric_info.noise_reference_method);
fprintf('Calcium noise reference: %s | formula noise = bleach_removed - noise_reference | sensitivity = bleach_removed ./ baseline | snr = bleach_removed ./ std(noise)\n', ...
    calcium_metric_info.noise_reference_method);

metric_results_file = fullfile(save_path, '3_dual_metric_results.mat');

voltage_results = store_trace_stage( ...
    voltage_results, 'noise_reference', voltage_noise_reference, {'bleach_removed'}, roi_results_file, ...
    voltage_results.movie_info, voltage_metric_info.noise_reference_method, voltage_metric_info.noise_reference_parameters);
voltage_results = store_trace_stage( ...
    voltage_results, 'noise', voltage_noise, {'bleach_removed', 'noise_reference'}, roi_results_file, ...
    voltage_results.movie_info, 'residual', struct('expression', 'bleach_removed - noise_reference'));
voltage_results = store_trace_stage( ...
    voltage_results, 'sensitivity', voltage_sensitivity, {'bleach_removed', 'baseline'}, roi_results_file, ...
    voltage_results.movie_info, 'ratio', struct('expression', 'bleach_removed ./ baseline'));
voltage_results = store_trace_stage( ...
    voltage_results, 'snr', voltage_snr, {'bleach_removed', 'noise'}, roi_results_file, ...
    voltage_results.movie_info, voltage_metric_info.snr_method, voltage_metric_info.snr_parameters);

calcium_results = store_trace_stage( ...
    calcium_results, 'noise_reference', calcium_noise_reference, {'bleach_removed'}, roi_results_file, ...
    calcium_results.movie_info, calcium_metric_info.noise_reference_method, calcium_metric_info.noise_reference_parameters);
calcium_results = store_trace_stage( ...
    calcium_results, 'noise', calcium_noise, {'bleach_removed', 'noise_reference'}, roi_results_file, ...
    calcium_results.movie_info, 'residual', struct('expression', 'bleach_removed - noise_reference'));
calcium_results = store_trace_stage( ...
    calcium_results, 'sensitivity', calcium_sensitivity, {'bleach_removed', 'baseline'}, roi_results_file, ...
    calcium_results.movie_info, 'ratio', struct('expression', 'bleach_removed ./ baseline'));
calcium_results = store_trace_stage( ...
    calcium_results, 'snr', calcium_snr, {'bleach_removed', 'noise'}, roi_results_file, ...
    calcium_results.movie_info, calcium_metric_info.snr_method, calcium_metric_info.snr_parameters);

fprintf('Applying final calcium smoothing...\n');
fprintf('  method: movmean\n');
fprintf('  window: %d frames\n', calcium_smoothing_window);
fprintf('  formula: smoothed = movmean(input, window, 1)\n');

% This is treated as the final calcium trace-processing step.
% We save the smoothed traces as separate stages instead of overwriting the
% unsmoothed versions, so the processing history stays explicit and older
% analyses remain reproducible.
calcium_raw_smoothed = movmean(double(fetch_trace_stage(calcium_results, 'raw')), calcium_smoothing_window, 1);
calcium_sensitivity_smoothed = movmean(double(calcium_sensitivity), calcium_smoothing_window, 1);
calcium_snr_smoothed = movmean(double(calcium_snr), calcium_smoothing_window, 1);

calcium_smoothing_info = struct('window', calcium_smoothing_window, 'dimension', 1);
calcium_results = store_trace_stage( ...
    calcium_results, 'raw_smoothed', calcium_raw_smoothed, {'raw'}, roi_results_file, ...
    calcium_results.movie_info, 'movmean', calcium_smoothing_info);
calcium_results = store_trace_stage( ...
    calcium_results, 'sensitivity_smoothed', calcium_sensitivity_smoothed, {'sensitivity'}, roi_results_file, ...
    calcium_results.movie_info, 'movmean', calcium_smoothing_info);
calcium_results = store_trace_stage( ...
    calcium_results, 'snr_smoothed', calcium_snr_smoothed, {'snr'}, roi_results_file, ...
    calcium_results.movie_info, 'movmean', calcium_smoothing_info);

metric_record = build_section_record( ...
    'Noise Estimation, Sensitivity, And SNR', ...
    'Estimate channel noise references, derive sensitivity/SNR traces, and apply the final calcium smoothing used by later sections.', ...
    struct( ...
        'traces_voltage_bleach', traces_voltage_bleach, ...
        'baseline_voltage', baseline_voltage, ...
        'traces_calcium_bleach', traces_calcium_bleach, ...
        'baseline_calcium', baseline_calcium, ...
        'traces_calcium_raw', fetch_trace_stage(calcium_results, 'raw'), ...
        'freq_calcium', freq_calcium), ...
    struct( ...
        'voltage_noise_reference_method', voltage_metric_info.noise_reference_method, ...
        'voltage_noise_reference_parameters', voltage_metric_info.noise_reference_parameters, ...
        'voltage_snr_method', voltage_metric_info.snr_method, ...
        'voltage_snr_parameters', voltage_metric_info.snr_parameters, ...
        'calcium_noise_reference_method', calcium_metric_info.noise_reference_method, ...
        'calcium_noise_reference_parameters', calcium_metric_info.noise_reference_parameters, ...
        'calcium_snr_method', calcium_metric_info.snr_method, ...
        'calcium_snr_parameters', calcium_metric_info.snr_parameters, ...
        'calcium_smoothing_method', 'movmean', ...
        'calcium_smoothing_parameters', struct('window', calcium_smoothing_window, 'dimension', 1)), ...
    struct( ...
        'voltage_noise_reference', voltage_noise_reference, ...
        'voltage_noise', voltage_noise, ...
        'voltage_sensitivity', voltage_sensitivity, ...
        'voltage_snr', voltage_snr, ...
        'calcium_noise_reference', calcium_noise_reference, ...
        'calcium_noise', calcium_noise, ...
        'calcium_sensitivity', calcium_sensitivity, ...
        'calcium_snr', calcium_snr, ...
        'calcium_raw_smoothed', calcium_raw_smoothed, ...
        'calcium_sensitivity_smoothed', calcium_sensitivity_smoothed, ...
        'calcium_snr_smoothed', calcium_snr_smoothed), ...
    struct(), ...
    'To rerun this section, load the saved bleach-stage inputs and parameters from this record, then call compute_voltage_metrics, compute_calcium_metrics, and the final calcium movmean smoothing in the same order.');
save(metric_results_file, ...
    'voltage_noise_reference', 'voltage_noise', 'voltage_sensitivity', 'voltage_snr', 'voltage_metric_info', ...
    'calcium_noise_reference', 'calcium_noise', 'calcium_sensitivity', 'calcium_snr', 'calcium_metric_info', ...
    'calcium_raw_smoothed', 'calcium_sensitivity_smoothed', 'calcium_snr_smoothed', 'calcium_smoothing_info', ...
    'metric_record');

save(voltage_results_path, 'voltage_results', '-v7.3');
save(calcium_results_path, 'calcium_results', '-v7.3');
fprintf('Metric output stages saved: noise_reference, noise, sensitivity, snr\n');
fprintf('Calcium final stages saved: raw_smoothed, sensitivity_smoothed, snr_smoothed\n');
t_voltage_metric = build_time_axis_from_movie_info(voltage_results.movie_info);
t_calcium_metric = build_time_axis_from_movie_info(calcium_results.movie_info);
[noise_ref_fig, noise_ref_png] = plot_dual_noise_reference_comparison( ...
    traces_voltage_bleach, voltage_noise_reference, ...
    traces_calcium_bleach, calcium_noise_reference, ...
    t_voltage_metric, t_calcium_metric, save_path);
dual_results.visualizations.noise_reference = struct( ...
    'data', struct(), ...
    'info', struct( ...
        'fig_file', noise_ref_fig, ...
        'png_file', noise_ref_png, ...
        'created_at', datetime("now")));
save(dual_results_path, 'dual_results', '-v7.3');
fprintf('Noise reference comparison saved to: %s\n', noise_ref_png);

%% Channel Summary Plots
% This figure is a checkpoint only. It does not create new trace stages;
% it simply shows which saved stage version is currently being used for
% each channel, especially on the calcium side.
print_section('Channel Summary Plots');
fprintf('Saving channel summary plots...\n');
[voltage_results, calcium_results] = load_channel_results( ...
    voltage_results_path, calcium_results_path, voltage_results, calcium_results);
traces_voltage_raw = fetch_trace_stage(voltage_results, 'raw');
has_calcium_raw = has_trace_stage(calcium_results, 'raw_smoothed') || has_trace_stage(calcium_results, 'raw');
has_voltage_sensitivity = has_trace_stage(voltage_results, 'sensitivity');
has_voltage_snr = has_trace_stage(voltage_results, 'snr');
has_calcium_sensitivity = has_trace_stage(calcium_results, 'sensitivity_smoothed') || has_trace_stage(calcium_results, 'sensitivity');
has_calcium_snr = has_trace_stage(calcium_results, 'snr_smoothed') || has_trace_stage(calcium_results, 'snr');
if has_voltage_sensitivity
    voltage_sensitivity = fetch_trace_stage(voltage_results, 'sensitivity');
else
    voltage_sensitivity = [];
end
if has_voltage_snr
    voltage_snr = fetch_trace_stage(voltage_results, 'snr');
else
    voltage_snr = [];
end
if has_calcium_raw
    % Prefer the final saved calcium stage when present, but keep a
    % fallback to older unsmoothed results for backward compatibility.
    [traces_calcium_raw, calcium_raw_stage] = resolve_preferred_trace_stage(calcium_results, {'raw_smoothed', 'raw'});
else
    traces_calcium_raw = [];
    calcium_raw_stage = '';
end
if has_calcium_sensitivity
    [calcium_sensitivity, calcium_sensitivity_stage] = resolve_preferred_trace_stage(calcium_results, {'sensitivity_smoothed', 'sensitivity'});
else
    calcium_sensitivity = [];
    calcium_sensitivity_stage = '';
end
if has_calcium_snr
    [calcium_snr, calcium_snr_stage] = resolve_preferred_trace_stage(calcium_results, {'snr_smoothed', 'snr'});
else
    calcium_snr = [];
    calcium_snr_stage = '';
end
summary_fig = figure('Color', 'w');
subplot(2, 3, 1); hold on; title('Voltage Raw'); plot_offset_stage(traces_voltage_raw, t_voltage);
subplot(2, 3, 2); hold on; title('Voltage Sensitivity'); plot_optional_stage(voltage_polarity * voltage_sensitivity, t_voltage, has_voltage_sensitivity, 'Sensitivity stage unavailable');
subplot(2, 3, 3); hold on; title('Voltage SNR'); plot_optional_stage(voltage_polarity * voltage_snr, t_voltage, has_voltage_snr, 'SNR stage unavailable');
subplot(2, 3, 4); hold on; title(sprintf('Calcium Raw (%s, window=%d)', strrep(char(calcium_raw_stage), '_', '\_'), calcium_smoothing_window)); plot_optional_stage(traces_calcium_raw, t_calcium, has_calcium_raw, 'Raw stage unavailable');
subplot(2, 3, 5); hold on; title(sprintf('Calcium Sensitivity (%s, window=%d)', strrep(char(calcium_sensitivity_stage), '_', '\_'), calcium_smoothing_window)); plot_optional_stage(calcium_polarity * calcium_sensitivity, t_calcium, has_calcium_sensitivity, 'Sensitivity stage unavailable');
subplot(2, 3, 6); hold on; title(sprintf('Calcium SNR (%s, window=%d)', strrep(char(calcium_snr_stage), '_', '\_'), calcium_smoothing_window)); plot_optional_stage(calcium_polarity * calcium_snr, t_calcium, has_calcium_snr, 'SNR stage unavailable');
save_figure_bundle(summary_fig, fullfile(save_path, '4_dual_trace_summary.fig'), fullfile(save_path, '4_dual_trace_summary.png'));
close(summary_fig);
fprintf('Summary plot stage availability | voltage sensitivity=%d snr=%d | calcium sensitivity=%d snr=%d\n', ...
    has_voltage_sensitivity, has_voltage_snr, has_calcium_sensitivity, has_calcium_snr);
fprintf('Channel summary polarity | voltage=%d | calcium=%d\n', voltage_polarity, calcium_polarity);
fprintf('Channel summary calcium stages | raw=%s | sensitivity=%s | snr=%s\n', ...
    string(calcium_raw_stage), string(calcium_sensitivity_stage), string(calcium_snr_stage));
fprintf('Channel summary calcium final smoothing window=%d\n', calcium_smoothing_window);

%% ROI Calcium Heatmap With Voltage Trace
% One ROI-aligned population view for inspecting voltage/calcium coupling.
%
% Display contract:
%   - each y-row is one paired ROI;
%   - calcium sensitivity is drawn first as the heatmap background;
%   - voltage sensitivity is then overlaid on top of the corresponding ROI
%     row as a red trace;
%   - the x coordinates keep each channel's saved physical time axis in
%     seconds. Voltage is not resampled onto calcium frames; this preserves
%     the original high-rate voltage timing while the heatmap x-limits come
%     from the calcium axis;
%   - the ROI order is the saved ROI column order from select_ROI_dual.
%
% Voltage vertical scaling rule:
%   - for every ROI, compute voltage max-min after display polarity is
%     applied;
%   - find the ROI with the largest max-min range;
%   - use that single largest range as the global scale for all ROIs;
%   - center every ROI trace by its own midpoint, then divide by the global
%     range. Therefore, the reference ROI's min lies at the bottom edge of
%     its heatmap row and its max lies at the top edge; all other ROIs are
%     shown on the same comparable scale.
%
% Calcium heatmap rule:
%   - prefer the final saved calcium stage sensitivity_smoothed, then fall
%     back to sensitivity for old result folders;
%   - after polarity, subtract each ROI's 1st percentile as that ROI's
%     displayed sensitivity zero;
%   - clip values below that percentile-zero baseline to 0 for display;
%   - use the black-blue-white "ice" colormap from heatmap_sensitivity.mlx;
%   - use [0 max] on the percentile-zeroed display matrix.
print_section('ROI Calcium Heatmap With Voltage Trace');
fprintf('Saving ROI-aligned calcium heatmap with overlaid voltage sensitivity traces...\n');
[voltage_results, calcium_results] = load_channel_results( ...
    voltage_results_path, calcium_results_path, voltage_results, calcium_results);

% Initialize a result record before any validation. This keeps
% dual_results.visualizations.roi_calcium_heatmap_voltage_trace present
% even if the section has to skip because an old result folder lacks one of
% the required sensitivity stages.
roi_heatmap_trace_result = struct( ...
    'status', "skipped", ...
    'reason', "", ...
    'fig_file', "", ...
    'png_file', "", ...
    'mat_file', "", ...
    'created_at', datetime("now"));

% The figure requires one voltage sensitivity matrix and one calcium
% sensitivity matrix. Calcium has a preferred final display stage because
% earlier sections intentionally create sensitivity_smoothed for calcium
% after the main metric calculation.
has_voltage_sensitivity = has_trace_stage(voltage_results, 'sensitivity');
has_calcium_sensitivity = has_trace_stage(calcium_results, 'sensitivity_smoothed') ...
    || has_trace_stage(calcium_results, 'sensitivity');
if ~has_voltage_sensitivity || ~has_calcium_sensitivity
    roi_heatmap_trace_result.reason = sprintf('Missing required sensitivity stage: voltage=%d calcium=%d.', ...
        has_voltage_sensitivity, has_calcium_sensitivity);
    fprintf('Skipping ROI calcium heatmap/voltage trace plot: %s\n', roi_heatmap_trace_result.reason);
else
    % Fetch saved trace stages and apply display polarity immediately.
    % From this point onward, "display" variables are exactly what will be
    % plotted and saved in the diagnostic MAT file.
    voltage_heatmap_stage = 'sensitivity';
    voltage_heatmap_display = double(voltage_polarity) * double(fetch_trace_stage(voltage_results, voltage_heatmap_stage));
    [calcium_heatmap_display, calcium_heatmap_stage] = resolve_preferred_trace_stage(calcium_results, {'sensitivity_smoothed', 'sensitivity'});
    calcium_heatmap_display = double(calcium_polarity) * double(calcium_heatmap_display);

    % Make the time vectors match their own trace matrices. This is only a
    % defensive repair for older or hand-edited result files:
    %   - if the saved time vector is too long, truncate it;
    %   - if it is too short, extend it using the median frame interval;
    %   - do not interpolate trace values.
    t_voltage_heatmap = double(t_voltage(:));
    voltage_time_info = struct('channel', "voltage", 'input_time_points', numel(t_voltage_heatmap), ...
        'trace_frames', size(voltage_heatmap_display, 1), 'rule', "unchanged");
    if numel(t_voltage_heatmap) > size(voltage_heatmap_display, 1)
        t_voltage_heatmap = t_voltage_heatmap(1:size(voltage_heatmap_display, 1));
        voltage_time_info.rule = "truncated time axis to trace frame count";
    elseif numel(t_voltage_heatmap) < size(voltage_heatmap_display, 1)
        if numel(t_voltage_heatmap) >= 2
            dt_heatmap = median(diff(t_voltage_heatmap), 'omitnan');
            if ~isfinite(dt_heatmap) || dt_heatmap <= 0
                dt_heatmap = 1;
            end
            first_t_heatmap = t_voltage_heatmap(1);
        else
            dt_heatmap = 1;
            first_t_heatmap = 0;
        end
        t_voltage_heatmap = first_t_heatmap + (0:size(voltage_heatmap_display, 1)-1)' * dt_heatmap;
        voltage_time_info.rule = "extended time axis using median dt";
    end

    t_calcium_heatmap = double(t_calcium(:));
    calcium_time_info = struct('channel', "calcium", 'input_time_points', numel(t_calcium_heatmap), ...
        'trace_frames', size(calcium_heatmap_display, 1), 'rule', "unchanged");
    if numel(t_calcium_heatmap) > size(calcium_heatmap_display, 1)
        t_calcium_heatmap = t_calcium_heatmap(1:size(calcium_heatmap_display, 1));
        calcium_time_info.rule = "truncated time axis to trace frame count";
    elseif numel(t_calcium_heatmap) < size(calcium_heatmap_display, 1)
        if numel(t_calcium_heatmap) >= 2
            dt_heatmap = median(diff(t_calcium_heatmap), 'omitnan');
            if ~isfinite(dt_heatmap) || dt_heatmap <= 0
                dt_heatmap = 1;
            end
            first_t_heatmap = t_calcium_heatmap(1);
        else
            dt_heatmap = 1;
            first_t_heatmap = 0;
        end
        t_calcium_heatmap = first_t_heatmap + (0:size(calcium_heatmap_display, 1)-1)' * dt_heatmap;
        calcium_time_info.rule = "extended time axis using median dt";
    end

    % Pair ROIs by column index. Dual ROI selection saves voltage and
    % calcium traces in matching ROI order, so column 1 is ROI 1 in both
    % channels. If a legacy file has unequal ROI counts, keep the shared
    % prefix and report the mismatch rather than guessing a remapping.
    nrois_voltage_heatmap = size(voltage_heatmap_display, 2);
    nrois_calcium_heatmap = size(calcium_heatmap_display, 2);
    nrois_heatmap = min(nrois_voltage_heatmap, nrois_calcium_heatmap);
    if nrois_heatmap == 0
        roi_heatmap_trace_result.reason = 'No shared ROI columns are available in the selected sensitivity stages.';
        fprintf('Skipping ROI calcium heatmap/voltage trace plot: %s\n', roi_heatmap_trace_result.reason);
    else
        if nrois_voltage_heatmap ~= nrois_calcium_heatmap
            fprintf('ROI heatmap/trace ROI mismatch: voltage=%d calcium=%d; plotting first %d paired ROIs.\n', ...
                nrois_voltage_heatmap, nrois_calcium_heatmap, nrois_heatmap);
        end
        voltage_heatmap_display = voltage_heatmap_display(:, 1:nrois_heatmap);
        calcium_heatmap_display_raw = calcium_heatmap_display(:, 1:nrois_heatmap);

        % Voltage overlay scaling:
        %   y_trace = roi_index + (trace - roi_midpoint) / global_range
        %
        % Because global_range is the largest ROI max-min value, the ROI
        % with the largest sensitivity swing fills exactly one heatmap row
        % from bottom edge (roi-0.5) to top edge (roi+0.5). Smaller ROIs
        % keep their relative amplitude on the same global scale.
        voltage_roi_min = min(voltage_heatmap_display, [], 1, 'omitnan');
        voltage_roi_max = max(voltage_heatmap_display, [], 1, 'omitnan');
        voltage_roi_range = voltage_roi_max - voltage_roi_min;
        [voltage_global_range, voltage_reference_roi] = max(voltage_roi_range);
        if ~isfinite(voltage_global_range) || voltage_global_range <= 0
            voltage_global_range = 1;
            voltage_reference_roi = 1;
        end
        voltage_roi_mid = (voltage_roi_min + voltage_roi_max) / 2;
        voltage_roi_mid(~isfinite(voltage_roi_mid)) = 0;

        % Calcium heatmap baseline rule:
        %   each ROI uses its own 1st percentile as displayed sensitivity 0.
        % Values below that baseline are clipped to 0 for display, while
        % the raw polarity-adjusted calcium matrix is still saved below.
        calcium_heatmap_zero_percentile = 1;
        calcium_heatmap_baseline = prctile(calcium_heatmap_display_raw, calcium_heatmap_zero_percentile, 1);
        calcium_heatmap_baseline(~isfinite(calcium_heatmap_baseline)) = 0;
        calcium_heatmap_display_zeroed = calcium_heatmap_display_raw - calcium_heatmap_baseline;
        calcium_heatmap_display_zeroed(calcium_heatmap_display_zeroed < 0) = 0;

        % After percentile-zeroing, color limits are intentionally [0 max].
        % A flat or empty heatmap falls back to [0 1] so clim remains valid.
        finite_calcium = calcium_heatmap_display_zeroed(isfinite(calcium_heatmap_display_zeroed));
        if isempty(finite_calcium)
            calcium_heatmap_clim = [0, 1];
            calcium_heatmap_clim_rule = "fallback_empty_to_0_1";
        else
            calcium_heatmap_max = max(finite_calcium);
            if calcium_heatmap_max > 0
                calcium_heatmap_clim = [0, calcium_heatmap_max];
                calcium_heatmap_clim_rule = "roi_1st_percentile_zero_to_max";
            else
                calcium_heatmap_clim = [0, 1];
                calcium_heatmap_clim_rule = "fallback_flat_to_0_1";
            end
        end
        if calcium_heatmap_clim(1) == calcium_heatmap_clim(2)
            calcium_heatmap_clim = calcium_heatmap_clim + [-0.5, 0.5];
            calcium_heatmap_clim_rule = calcium_heatmap_clim_rule + "_expanded_flat_range";
        end

        % Inline version of heatmap_sensitivity.mlx custom_ice_adjust:
        % black background, blue/cyan midrange, white strongest values.
        % gamma controls how much low activity stays dark; pivot controls
        % how quickly blue transitions toward brighter cyan/white.
        gamma_val_calcium = 0.6;
        color_pivot_calcium = 0.4;
        cmap_n = 256;
        base_ice = zeros(cmap_n, 3);
        base_ice(:, 3) = linspace(0, 1, cmap_n);
        cyan_start_node = max(1, floor(cmap_n * 0.1));
        base_ice(cyan_start_node:end, 2) = linspace(0, 1, cmap_n - cyan_start_node + 1);
        white_start_node = max(1, floor(cmap_n * 0.9));
        base_ice(white_start_node:end, 1) = linspace(0, 1, cmap_n - white_start_node + 1);
        x_old = linspace(0, 1, cmap_n);
        x_new = linspace(0, 1, cmap_n) .^ gamma_val_calcium;
        warped = interp1([0, color_pivot_calcium, 1], [0, 0.5, 1], x_new, 'linear', 'extrap');
        warped = min(max(warped, 0), 1);
        calcium_heatmap_cmap = interp1(x_old, base_ice, warped);

        % Diagnostic transfer-curve plot for tuning the calcium heatmap.
        % normalized input x means the calcium value after clim maps into
        % [0, 1]. The dashed curve is x^gamma; the solid curve is the final
        % colormap index after applying pivot. The vertical pivot line marks
        % the raw x value where x^gamma reaches color_pivot_calcium.
        gamma_curve_input = linspace(0, 1, cmap_n);
        gamma_curve_after_gamma = gamma_curve_input .^ gamma_val_calcium;
        gamma_curve_warped = interp1([0, color_pivot_calcium, 1], [0, 0.5, 1], ...
            gamma_curve_after_gamma, 'linear', 'extrap');
        gamma_curve_warped = min(max(gamma_curve_warped, 0), 1);
        if gamma_val_calcium > 0
            gamma_curve_pivot_input = color_pivot_calcium .^ (1 / gamma_val_calcium);
        else
            gamma_curve_pivot_input = NaN;
        end

        colormap_curve_fig = figure( ...
            'Name', 'Calcium Heatmap Colormap Transfer Curve', ...
            'Color', 'w', ...
            'Units', 'pixels', ...
            'Position', [120, 120, 920, 620]);
        curve_layout = tiledlayout(colormap_curve_fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
        ax_curve = nexttile(curve_layout, 1);
        plot(ax_curve, gamma_curve_input, gamma_curve_after_gamma, '--', 'Color', [0.25 0.25 0.25], 'LineWidth', 1.3);
        hold(ax_curve, 'on');
        plot(ax_curve, gamma_curve_input, gamma_curve_warped, 'b-', 'LineWidth', 1.8);
        if isfinite(gamma_curve_pivot_input)
            xline(ax_curve, gamma_curve_pivot_input, ':', sprintf('pivot x=%.3g', gamma_curve_pivot_input), ...
                'Color', [0.2 0.2 0.2], 'LabelVerticalAlignment', 'bottom');
        end
        yline(ax_curve, 0.5, ':', 'colormap midpoint', 'Color', [0.45 0.45 0.45]);
        xlabel(ax_curve, 'Normalized calcium value after clim');
        ylabel(ax_curve, 'Colormap lookup index');
        title(ax_curve, sprintf('Calcium colormap transfer | gamma=%.3g, pivot=%.3g', ...
            gamma_val_calcium, color_pivot_calcium));
        legend(ax_curve, {'x^{gamma}', 'warped final index'}, 'Location', 'southeast');
        grid(ax_curve, 'on');
        ylim(ax_curve, [0 1]);

        ax_strip = nexttile(curve_layout, 2);
        image(ax_strip, gamma_curve_input, 1, reshape(calcium_heatmap_cmap, [1, cmap_n, 3]));
        set(ax_strip, 'YTick', [], 'TickDir', 'out');
        xlabel(ax_strip, 'Normalized calcium value after clim');
        title(ax_strip, 'Resulting black-blue-white color strip');
        xlim(ax_strip, [0 1]);

        colormap_curve_fig_file = fullfile(save_path, '4_dual_roi_calcium_heatmap_colormap_curve.fig');
        colormap_curve_png_file = fullfile(save_path, '4_dual_roi_calcium_heatmap_colormap_curve.png');
        save_figure_bundle_preserve_layout(colormap_curve_fig, colormap_curve_fig_file, colormap_curve_png_file);
        close(colormap_curve_fig);

        % Figure height grows with ROI count but is capped so export remains
        % manageable. save_figure_bundle_preserve_layout is used here so
        % this tall layout is not maximized and distorted before PNG export.
        fig_height = min(2200, max(650, 24 * nrois_heatmap + 180));
        roi_heatmap_fig = figure( ...
            'Name', 'ROI Calcium Heatmap With Voltage Sensitivity Trace', ...
            'Color', 'w', ...
            'Units', 'pixels', ...
            'Position', [80, 80, 1500, fig_height]);
        ax = axes(roi_heatmap_fig);
        imagesc(ax, t_calcium_heatmap, 1:nrois_heatmap, calcium_heatmap_display_zeroed');
        set(ax, 'YDir', 'normal', 'TickDir', 'out', 'Layer', 'top');
        colormap(ax, calcium_heatmap_cmap);
        clim(ax, calcium_heatmap_clim);
        hold(ax, 'on');
        for roi_idx = 1:nrois_heatmap
            % The centered and globally scaled trace is drawn directly in
            % ROI-row coordinates, so y=roi_idx is that row's midline.
            y_trace = roi_idx + (voltage_heatmap_display(:, roi_idx) - voltage_roi_mid(roi_idx)) / voltage_global_range *1.5;
            plot(ax, t_voltage_heatmap, y_trace, 'Color', [1.0, 0.12, 0.02], 'LineWidth', 0.45);
        end
        xlim(ax, [min(t_calcium_heatmap), max(t_calcium_heatmap)]);
        ylim(ax, [0.5, nrois_heatmap + 0.5]);
        if nrois_heatmap <= 20
            yticks(ax, 1:nrois_heatmap);
        else
            roi_tick_step = max(1, ceil(nrois_heatmap / 20));
            yticks(ax, unique([1:roi_tick_step:nrois_heatmap, nrois_heatmap]));
        end
        xlabel(ax, 'Time (s)');
        ylabel(ax, 'ROI');
        title(ax, sprintf('Calcium Sensitivity Heatmap (%s) + Voltage Sensitivity Trace', ...
            strrep(char(calcium_heatmap_stage), '_', '\_')));
        cb = colorbar(ax);
        ylabel(cb, 'Calcium sensitivity display value');
        grid(ax, 'on');
        box(ax, 'off');

        roi_heatmap_fig_file = fullfile(save_path, '4_dual_roi_calcium_heatmap_voltage_trace.fig');
        roi_heatmap_png_file = fullfile(save_path, '4_dual_roi_calcium_heatmap_voltage_trace.png');
        roi_heatmap_mat_file = fullfile(save_path, '4_dual_roi_calcium_heatmap_voltage_trace.mat');
        save_figure_bundle_preserve_layout(roi_heatmap_fig, roi_heatmap_fig_file, roi_heatmap_png_file);
        close(roi_heatmap_fig);

        % Store enough provenance for later interpretation without forcing
        % the user to infer which stage, polarity, scale, ROI truncation,
        % and color-limit rules were active when the figure was generated.
        roi_heatmap_trace_result = struct( ...
            'status', "completed", ...
            'reason', "", ...
            'fig_file', string(roi_heatmap_fig_file), ...
            'png_file', string(roi_heatmap_png_file), ...
            'mat_file', string(roi_heatmap_mat_file), ...
            'colormap_curve_fig_file', string(colormap_curve_fig_file), ...
            'colormap_curve_png_file', string(colormap_curve_png_file), ...
            'input_stages', struct('voltage', string(voltage_heatmap_stage), 'calcium', string(calcium_heatmap_stage)), ...
            'parameters', struct( ...
                'voltage_polarity', voltage_polarity, ...
                'calcium_polarity', calcium_polarity, ...
                'calcium_smoothing_window', calcium_smoothing_window, ...
                'calcium_colormap', "custom_ice_adjust_inline", ...
                'calcium_colormap_gamma', gamma_val_calcium, ...
                'calcium_colormap_pivot', color_pivot_calcium, ...
                'calcium_zero_rule', "per-ROI 1st percentile subtracted, then values below zero clipped", ...
                'calcium_zero_percentile', calcium_heatmap_zero_percentile, ...
                'calcium_clim', calcium_heatmap_clim, ...
                'calcium_clim_rule', string(calcium_heatmap_clim_rule), ...
                'voltage_scale_rule', "global max-min range; each ROI centered by its own midpoint", ...
                'voltage_global_range', voltage_global_range, ...
                'voltage_reference_roi', voltage_reference_roi, ...
                'roi_count_plotted', nrois_heatmap, ...
                'roi_count_voltage', nrois_voltage_heatmap, ...
                'roi_count_calcium', nrois_calcium_heatmap, ...
                'x_axis_source', "calcium time axis", ...
                'voltage_trace_time_rule', "overlay voltage trace using voltage seconds without resampling"), ...
            'time_axis_info', struct('voltage', voltage_time_info, 'calcium', calcium_time_info), ...
            'created_at', datetime("now"));

        % Save the exact displayed matrices rather than the raw stage data:
        % polarity has already been applied and the ROI columns have already
        % been restricted to the paired shared set. This makes the MAT file
        % match the figure pixel-for-pixel.
        save(roi_heatmap_mat_file, ...
            'roi_heatmap_trace_result', ...
            'voltage_heatmap_display', ...
            'calcium_heatmap_display_raw', ...
            'calcium_heatmap_display_zeroed', ...
            'calcium_heatmap_baseline', ...
            't_voltage_heatmap', ...
            't_calcium_heatmap', ...
            'voltage_roi_min', ...
            'voltage_roi_max', ...
            'voltage_roi_range', ...
            'voltage_reference_roi', ...
            'gamma_curve_input', ...
            'gamma_curve_after_gamma', ...
            'gamma_curve_warped', ...
            'gamma_curve_pivot_input', ...
            '-v7.3');
    end
end
if ~isfield(dual_results, 'visualizations') || ~isstruct(dual_results.visualizations)
    dual_results.visualizations = struct();
end
dual_results.visualizations.roi_calcium_heatmap_voltage_trace = roi_heatmap_trace_result;
save(dual_results_path, 'dual_results', '-v7.3');
fprintf('ROI heatmap/trace status: %s\n', string(roi_heatmap_trace_result.status));
if isfield(roi_heatmap_trace_result, 'png_file') && strlength(string(roi_heatmap_trace_result.png_file)) > 0
    fprintf('ROI heatmap/trace saved to: %s\n', roi_heatmap_trace_result.png_file);
end

%% Dual Comparison
% This section asks whether the processed voltage and calcium traces tell a
% consistent ROI-by-ROI story after each channel finishes its own
% preprocessing.
print_section('Dual Comparison');
fprintf('Building dual-channel comparison results...\n');
stim_windows = struct();
if stim_context.supported
    if analysis_only_mode && isstruct(stim_results) && isfield(stim_results, 'windows') && ~isempty(stim_results.windows)
        stim_windows = stim_results.windows;
        fprintf('Analysis-only mode: reusing saved stim windows instead of rebuilding them from logs.\n');
    elseif analysis_only_mode && ~(isfield(stim_context, 'logs') && isstruct(stim_context.logs))
        fprintf(['Analysis-only mode: saved stim windows are unavailable and saved stim_context has no logs. ' ...
            'Rebuilding stim metadata from the current cycle manifests/logs.\n']);
        voltage_idx_local = find(strcmpi(string({camera_cfg.role}), "voltage"), 1, 'first');
        calcium_idx_local = find(strcmpi(string({camera_cfg.role}), "calcium"), 1, 'first');
        [cycle_manifest_reload, record_manifest_reload, ~, input_layout_reload] = resolve_dual_camera_sources( ...
            cycle_path, camera_cfg, camera_source_override, raw_dual_input_cfg);
        stim_context = resolve_stim_context( ...
            cycle_path, cycle_manifest_reload, record_manifest_reload, ...
            camera_cfg(voltage_idx_local).camera_index, camera_cfg(calcium_idx_local).camera_index, ...
            stim_context_override, input_layout_reload);
        stim_windows = build_visual_stim_windows(stim_context, freq_voltage, freq_calcium, nframes_voltage, nframes_calcium);
    else
        stim_windows = build_visual_stim_windows(stim_context, freq_voltage, freq_calcium, nframes_voltage, nframes_calcium);
    end
    print_stim_context_summary(stim_context, stim_windows, 'Dual Comparison');
else
    print_stim_context_summary(stim_context, struct(), 'Dual Comparison');
end
[voltage_results, calcium_results] = load_channel_results( ...
    voltage_results_path, calcium_results_path, voltage_results, calcium_results);
has_sensitivity_pair = has_trace_stage(voltage_results, 'sensitivity') ...
    && (has_trace_stage(calcium_results, 'sensitivity_smoothed') || has_trace_stage(calcium_results, 'sensitivity'));
has_snr_pair = has_trace_stage(voltage_results, 'snr') ...
    && (has_trace_stage(calcium_results, 'snr_smoothed') || has_trace_stage(calcium_results, 'snr'));
if has_sensitivity_pair
    % Voltage keeps its native saved stage, while calcium may switch to its
    % final-smoothed stage if that stage exists.
    voltage_sensitivity = fetch_trace_stage(voltage_results, 'sensitivity');
    [calcium_sensitivity, calcium_sensitivity_stage] = resolve_preferred_trace_stage(calcium_results, {'sensitivity_smoothed', 'sensitivity'});
else
    voltage_sensitivity = [];
    calcium_sensitivity = [];
    calcium_sensitivity_stage = '';
end
if has_snr_pair
    voltage_snr = fetch_trace_stage(voltage_results, 'snr');
    [calcium_snr, calcium_snr_stage] = resolve_preferred_trace_stage(calcium_results, {'snr_smoothed', 'snr'});
else
    voltage_snr = [];
    calcium_snr = [];
    calcium_snr_stage = '';
end
fprintf('Dual comparison timing and polarity:\n');
fprintf('  voltage fps used: %g\n', voltage_results.movie_info.frame_rate);
fprintf('  calcium fps used: %g\n', calcium_results.movie_info.frame_rate);
fprintf('  fps source: movie_info.frame_rate <- camera_cfg.frame_rate\n');
fprintf('  voltage polarity/display: %d\n', voltage_polarity);
fprintf('  calcium polarity/display: %d\n', calcium_polarity);
fprintf('  overlap plot uses t_voltage for voltage and t_calcium for calcium\n');
fprintf('Dual comparison available stage pairs: sensitivity=%d | snr=%d\n', has_sensitivity_pair, has_snr_pair);
fprintf('Dual comparison input stages | voltage sensitivity=sensitivity | voltage snr=snr | calcium sensitivity=%s | calcium snr=%s\n', ...
    string(calcium_sensitivity_stage), string(calcium_snr_stage));
fprintf('Dual comparison formulas | display traces = saved stage values | calcium final smoothing window=%d | integral uses leaky accumulator | deconvolution uses FOOPSI when available\n', ...
    calcium_smoothing_window);
fprintf('Dual comparison statistics | correlation figure = paired ROI vs shuffled ROI distributions (Pearson and Spearman)\n');

comparison_data = struct();
comparison_info = struct();

if has_sensitivity_pair
    % Build one comparison package per metric representation so sensitivity
    % and SNR can be inspected independently.
    comparison_sensitivity = build_dual_metric_comparison( ...
    'sensitivity', voltage_sensitivity, calcium_sensitivity, ...
    t_voltage, t_calcium, calcium_smoothing_window, save_path, stim_windows, voltage_polarity, calcium_polarity);
    comparison_data.sensitivity = comparison_sensitivity.data;
    comparison_info.sensitivity = comparison_sensitivity.info;
else
    fprintf('Skipping dual sensitivity comparison: sensitivity stage missing in saved channel results\n');
end

if has_snr_pair
    comparison_snr = build_dual_metric_comparison( ...
        'snr', voltage_snr, calcium_snr, ...
        t_voltage, t_calcium, calcium_smoothing_window, save_path, stim_windows, voltage_polarity, calcium_polarity);
    comparison_data.snr = comparison_snr.data;
    comparison_info.snr = comparison_snr.info;
else
    fprintf('Skipping dual SNR comparison: snr stage missing in saved channel results\n');
end

if ~has_sensitivity_pair && ~has_snr_pair
    error('Dual comparison cannot run because neither sensitivity nor snr stage is available in both saved channel results.');
end

dual_results.comparison = struct( ...
    'data', comparison_data, ...
    'info', struct( ...
        'method', 'dual_metric_comparison_with_quad_and_overlap_plots', ...
        'parameters', struct('calcium_smoothing_window', calcium_smoothing_window), ...
        'input_stages', struct( ...
            'voltage_sensitivity', 'sensitivity', ...
            'voltage_snr', 'snr', ...
            'calcium_sensitivity', string(calcium_sensitivity_stage), ...
            'calcium_snr', string(calcium_snr_stage)), ...
        'available_stage_pairs', struct('sensitivity', has_sensitivity_pair, 'snr', has_snr_pair), ...
        'comparison_files', comparison_info, ...
        'created_at', datetime("now")));
save(dual_results_path, 'dual_results', '-v7.3');
fprintf('Dual comparison saved to: %s\n', dual_results_path);

%% Stim-Aware Analysis
% Add stimulus shading and trial-based stimulation metrics when the cycle
% includes visual stimulation logs.
% This reconnects the processed traces to stimulus time and turns them into
% trial-based response summaries that can be compared across conditions.
print_section('Stim-Aware Analysis');
fprintf('Running stimulus-specific analysis after dual comparison...\n');
[voltage_results, calcium_results] = load_channel_results( ...
    voltage_results_path, calcium_results_path, voltage_results, calcium_results);
has_voltage_sensitivity = has_trace_stage(voltage_results, 'sensitivity');
has_voltage_snr = has_trace_stage(voltage_results, 'snr');
has_calcium_sensitivity = has_trace_stage(calcium_results, 'sensitivity_smoothed') || has_trace_stage(calcium_results, 'sensitivity');
has_calcium_snr = has_trace_stage(calcium_results, 'snr_smoothed') || has_trace_stage(calcium_results, 'snr');
if has_voltage_sensitivity
    voltage_sensitivity = fetch_trace_stage(voltage_results, 'sensitivity');
else
    voltage_sensitivity = [];
end
if has_voltage_snr
    voltage_snr = fetch_trace_stage(voltage_results, 'snr');
else
    voltage_snr = [];
end
if has_calcium_sensitivity
    [calcium_sensitivity, calcium_sensitivity_stage] = resolve_preferred_trace_stage(calcium_results, {'sensitivity_smoothed', 'sensitivity'});
else
    calcium_sensitivity = [];
    calcium_sensitivity_stage = '';
end
if has_calcium_snr
    [calcium_snr, calcium_snr_stage] = resolve_preferred_trace_stage(calcium_results, {'snr_smoothed', 'snr'});
else
    calcium_snr = [];
    calcium_snr_stage = '';
end

stim_results.info = rmfield_if_exists(stim_context, {'logs', 'method_manifest'});
stim_results.windows = stim_windows;

if stim_context.supported
    % The same saved stages used for plotting are also used for numeric
    % stimulus metrics, so visualization and quantification stay aligned.
    print_stim_context_summary(stim_context, stim_windows, 'Stim-Aware Analysis');
    stim_results.status = ternary(isstruct(stim_windows) && isfield(stim_windows, 'supported') && stim_windows.supported, ...
        'supported', 'unsupported_window_layout');

    if stim_windows.supported
        if ~(has_voltage_snr && has_calcium_snr)
            error('Stim-aware analysis requires snr stage in both voltage_results and calcium_results.');
        end
        fprintf('Stim response metrics:\n');
        fprintf('  trace stages | voltage snr=snr | calcium snr=%s | calcium sensitivity=%s\n', ...
            string(calcium_snr_stage), string(calcium_sensitivity_stage));
        fprintf('  voltage polarity: %d\n', voltage_polarity);
        fprintf('  calcium polarity: %d\n', calcium_polarity);
        fprintf('  formulas: delta_mean, delta_peak, delta_auc = stim - baseline\n');
        fprintf('  calcium final smoothing window: %d\n', calcium_smoothing_window);
        fprintf('  plotted stages: sensitivity and snr only\n');
        if strcmpi(string(stim_context.stim_type), "visualstim_flash")
            fprintf('  flash window rule: baseline = gray 2 s before flash onset | response = 2 s from flash onset\n');
            fprintf('  flash shading rule: actual flash pulse only\n');
        end

        stim_trace_fig = figure('Color', 'w');
        subplot(2, 2, 1); hold on; title('Voltage Sensitivity'); plot_optional_stage(voltage_polarity * voltage_sensitivity, t_voltage, has_voltage_sensitivity, 'Sensitivity stage unavailable'); add_stim_shading(gca, stim_windows.voltage, stim_windows.condition_index, stim_windows.condition_colors, 0.14, stim_windows.trial_labels, stim_windows.block_labels);
        subplot(2, 2, 2); hold on; title('Voltage SNR'); plot_optional_stage(voltage_polarity * voltage_snr, t_voltage, has_voltage_snr, 'SNR stage unavailable'); add_stim_shading(gca, stim_windows.voltage, stim_windows.condition_index, stim_windows.condition_colors, 0.14, stim_windows.trial_labels, stim_windows.block_labels);
        subplot(2, 2, 3); hold on; title(sprintf('Calcium Sensitivity (%s)', strrep(char(calcium_sensitivity_stage), '_', '\_'))); plot_optional_stage(calcium_polarity * calcium_sensitivity, t_calcium, has_calcium_sensitivity, 'Sensitivity stage unavailable'); add_stim_shading(gca, stim_windows.calcium, stim_windows.condition_index, stim_windows.condition_colors, 0.14, stim_windows.trial_labels, stim_windows.block_labels);
        subplot(2, 2, 4); hold on; title(sprintf('Calcium SNR (%s)', strrep(char(calcium_snr_stage), '_', '\_'))); plot_optional_stage(calcium_polarity * calcium_snr, t_calcium, has_calcium_snr, 'SNR stage unavailable'); add_stim_shading(gca, stim_windows.calcium, stim_windows.condition_index, stim_windows.condition_colors, 0.14, stim_windows.trial_labels, stim_windows.block_labels);
        save_figure_bundle(stim_trace_fig, fullfile(save_path, '4_dual_trace_summary_with_stim.fig'), fullfile(save_path, '4_dual_trace_summary_with_stim.png'));
        close(stim_trace_fig);

        voltage_stim_metrics_snr = compute_stim_trial_metrics(voltage_snr, stim_windows.voltage, voltage_polarity, freq_voltage);
        calcium_stim_metrics_snr = compute_stim_trial_metrics(calcium_snr, stim_windows.calcium, calcium_polarity, freq_calcium);

        stim_results.response = struct( ...
            'analysis_trace_stage', struct( ...
                'snr', struct('voltage', "snr", 'calcium', string(calcium_snr_stage)), ...
                'sensitivity', struct('voltage', "sensitivity", 'calcium', string(calcium_sensitivity_stage))), ...
            'plot_trace_stage', struct('voltage_sensitivity', "sensitivity", 'voltage_snr', "snr", ...
                'calcium_sensitivity', string(calcium_sensitivity_stage), 'calcium_snr', string(calcium_snr_stage)), ...
            'calcium_smoothing_window', calcium_smoothing_window, ...
            'voltage_polarity', voltage_polarity, ...
            'calcium_polarity', calcium_polarity, ...
            'snr', struct('voltage', voltage_stim_metrics_snr, 'calcium', calcium_stim_metrics_snr));

        if has_voltage_sensitivity && has_calcium_sensitivity
            voltage_stim_metrics_sensitivity = compute_stim_trial_metrics(voltage_sensitivity, stim_windows.voltage, voltage_polarity, freq_voltage);
            calcium_stim_metrics_sensitivity = compute_stim_trial_metrics(calcium_sensitivity, stim_windows.calcium, calcium_polarity, freq_calcium);
            stim_results.response.sensitivity = struct( ...
                'voltage', voltage_stim_metrics_sensitivity, ...
                'calcium', calcium_stim_metrics_sensitivity);
            fprintf('  condition-response metrics will be generated for both sensitivity and snr\n');
        else
            voltage_stim_metrics_sensitivity = [];
            calcium_stim_metrics_sensitivity = [];
            fprintf('  sensitivity metrics skipped because one or both sensitivity stages are unavailable\n');
        end

        if stim_windows.is_grating
            fprintf('Stim analysis branch: grating tuning\n');
            fprintf('  OSI/DSI: enabled\n');
            fprintf('  outputs: pref_dir, pref_ori, gOSI, gDSI, OSI, DSI\n');
            [voltage_tuning, calcium_tuning, per_roi_tuning] = compute_dual_grating_tuning_by_roi( ...
                voltage_stim_metrics_snr.delta_mean, calcium_stim_metrics_snr.delta_mean, stim_windows.orientations);
            stim_results.analysis_kind = 'grating_tuning';
            stim_results.tuning = struct( ...
                'voltage', voltage_tuning, ...
                'calcium', calcium_tuning, ...
                'per_roi', per_roi_tuning);

            per_roi_tuning = save_grating_tuning_by_roi(per_roi_tuning, save_path, '7_stim_tuning_by_roi');
            stim_results.tuning.per_roi = per_roi_tuning;
            save(fullfile(save_path, '7_stim_tuning_summary.mat'), 'voltage_tuning', 'calcium_tuning', 'per_roi_tuning');

            stim_tuning_fig = figure('Color', 'w');
            subplot(2, 2, 1);
            plot_tuning_population(voltage_tuning.unique_orientations, voltage_tuning.response_by_condition, 'r', 'Voltage');
            subplot(2, 2, 2);
            plot_tuning_population(calcium_tuning.unique_orientations, calcium_tuning.response_by_condition, 'g', 'Calcium');
            subplot(2, 2, 3);
            scatter(voltage_tuning.gosi, calcium_tuning.gosi, 28, 'filled');
            xlabel('Voltage gOSI'); ylabel('Calcium gOSI'); title('Orientation Selectivity'); grid on;
            subplot(2, 2, 4);
            scatter(voltage_tuning.pref_dir, calcium_tuning.pref_dir, 28, 'filled');
            hold on; plot([0 360], [0 360], 'k--');
            xlim([0 360]); ylim([0 360]);
            xlabel('Voltage Pref. Dir (deg)'); ylabel('Calcium Pref. Dir (deg)');
            title('Preferred Direction'); grid on;
            save_figure_bundle(stim_tuning_fig, fullfile(save_path, '7_stim_tuning_summary.fig'), fullfile(save_path, '7_stim_tuning_summary.png'));
            close(stim_tuning_fig);
        else
            fprintf('Stim analysis branch: condition response\n');
            fprintf('  OSI/DSI: skipped\n');
            fprintf('  reason: non-grating stimulus type = %s\n', stim_context.stim_type);
            if strcmpi(string(stim_context.stim_type), "visualstim_flash")
                fprintf('  flash condition summary uses post-flash 2 s response windows\n');
                fprintf('  flash delta summary compares response 2 s against the preceding gray 2 s baseline\n');
            end
            voltage_block_summary_snr = summarize_block_by_condition( ...
                voltage_snr, stim_windows.voltage, stim_windows.block_labels, voltage_polarity, stim_windows.trial_labels);
            calcium_block_summary_snr = summarize_block_by_condition( ...
                calcium_snr, stim_windows.calcium, stim_windows.block_labels, calcium_polarity, stim_windows.trial_labels);
            voltage_delta_summary_snr = summarize_stim_by_condition(voltage_stim_metrics_snr, stim_windows);
            calcium_delta_summary_snr = summarize_stim_by_condition(calcium_stim_metrics_snr, stim_windows);
            stim_results.analysis_kind = 'condition_response';
            stim_results.condition_summary = struct( ...
                'snr', struct( ...
                    'block', struct('voltage', voltage_block_summary_snr, 'calcium', calcium_block_summary_snr), ...
                    'delta', struct('voltage', voltage_delta_summary_snr, 'calcium', calcium_delta_summary_snr)));

            plot_stim_condition_summary( ...
                voltage_block_summary_snr, calcium_block_summary_snr, ...
                stim_context, save_path, '7_stim_block_summary_snr', 'SNR Block Response', 'SNR');
            plot_stim_condition_summary( ...
                voltage_delta_summary_snr, calcium_delta_summary_snr, ...
                stim_context, save_path, '7_stim_delta_summary_snr', 'SNR Delta Response', '\Delta SNR (stim - baseline)');

            if ~isempty(voltage_stim_metrics_sensitivity) && ~isempty(calcium_stim_metrics_sensitivity)
                voltage_block_summary_sensitivity = summarize_block_by_condition( ...
                    voltage_sensitivity, stim_windows.voltage, stim_windows.block_labels, voltage_polarity, stim_windows.trial_labels);
                calcium_block_summary_sensitivity = summarize_block_by_condition( ...
                    calcium_sensitivity, stim_windows.calcium, stim_windows.block_labels, calcium_polarity, stim_windows.trial_labels);
                voltage_delta_summary_sensitivity = summarize_stim_by_condition(voltage_stim_metrics_sensitivity, stim_windows);
                calcium_delta_summary_sensitivity = summarize_stim_by_condition(calcium_stim_metrics_sensitivity, stim_windows);

                stim_results.condition_summary.sensitivity = struct( ...
                    'block', struct('voltage', voltage_block_summary_sensitivity, 'calcium', calcium_block_summary_sensitivity), ...
                    'delta', struct('voltage', voltage_delta_summary_sensitivity, 'calcium', calcium_delta_summary_sensitivity));

                plot_stim_condition_summary( ...
                    voltage_block_summary_sensitivity, calcium_block_summary_sensitivity, ...
                    stim_context, save_path, '7_stim_block_summary_sensitivity', 'Sensitivity Block Response', 'Sensitivity');
                plot_stim_condition_summary( ...
                    voltage_delta_summary_sensitivity, calcium_delta_summary_sensitivity, ...
                    stim_context, save_path, '7_stim_delta_summary_sensitivity', 'Sensitivity Delta Response', '\Delta Sensitivity (stim - baseline)');
            end
        end
    else
        fprintf('Stim analysis skipped after classification:\n');
        fprintf('  reason: stim windows could not be reconstructed\n');
        fprintf('  stim type: %s\n', stim_context.stim_type);
        fprintf('  selected program: %s\n', stim_context.selected_program);
    end
else
    stim_results.status = 'not_applicable_or_missing_inputs';
    fprintf('Stim-aware analysis skipped:\n');
    fprintf('  reason: required stimulation metadata incomplete\n');
    fprintf('  record mode: %s\n', string(stim_context.recordmode));
    fprintf('  has logs: %d\n', stim_context.has_logs);
    fprintf('  has stimSpec: %d\n', stim_context.has_stimSpec);
end
save(stim_results_path, 'stim_results', '-v7.3');

%% Population Time-Frequency Analysis
% Replace the previous "oscillation event" summary with a more direct
% frequency-domain view:
%   1. Direct FFT amplitude spectra for dominant frequency content across ROIs
%   2. Wavelet scalograms for how that content evolves over time
% Important:
%   use unsmoothed channel stages for spectral estimation here. The
%   calcium movmean stage is still useful for time-domain display in other
%   sections, but it can imprint comb-like ripple structure onto the PSD.
print_section('Population Time-Frequency Analysis');
fprintf('Running Fourier and wavelet analysis after dual/stim summaries...\n');
[voltage_results, calcium_results] = load_channel_results( ...
    voltage_results_path, calcium_results_path, voltage_results, calcium_results);

time_frequency_results = run_dual_time_frequency_section( ...
    voltage_results, calcium_results, ...
    t_voltage, t_calcium, ...
    freq_voltage, freq_calcium, ...
    stim_windows, ...
    save_path, ...
    voltage_polarity, calcium_polarity);

time_frequency_record = build_section_record( ...
    'Population Time-Frequency Analysis', ...
    'Summarize processed dual-channel traces with direct FFT spectra and wavelet time-frequency maps, while reusing the saved processed trace stages as section inputs.', ...
    struct( ...
        'voltage_results_file', string(voltage_results_path), ...
        'calcium_results_file', string(calcium_results_path), ...
        'stim_results_file', string(stim_results_path), ...
        'voltage_stage', string(time_frequency_results.parameters.voltage_stage), ...
        'calcium_stage', string(time_frequency_results.parameters.calcium_stage)), ...
    time_frequency_results.parameters, ...
    struct( ...
        'fourier_summary_png', string(time_frequency_results.visualizations.fourier_summary.png_file), ...
        'wavelet_summary_png', string(time_frequency_results.visualizations.wavelet_summary.png_file)), ...
    struct(), ...
    "To rerun the ROI-after analysis pipeline, set analysis_mode = 'analysis_only', point reuse_results_path to an existing Dual_analysis3 output folder, and run Dual_analysis3 again.");

time_frequency_results.record = time_frequency_record;

time_frequency_results_path = fullfile(save_path, '8_time_frequency_results.mat');
save(time_frequency_results_path, 'time_frequency_results', '-v7.3');
fprintf('Fourier summary saved to: %s\n', time_frequency_results.visualizations.fourier_summary.png_file);
fprintf('Wavelet summary saved to: %s\n', time_frequency_results.visualizations.wavelet_summary.png_file);
fprintf('Time-frequency result bundle saved to: %s\n', time_frequency_results_path);

%% Save Explicit Results
% Save one bundled snapshot for manual inspection. The modular section
% files remain the main working outputs, while this file is a convenient
% "open everything at once" archive.
print_section('Save Explicit Results');
save_explicit_dual_results_summary( ...
    save_path, ...
    dual_info, dual_info_path, ...
    voltage_results, voltage_results_path, ...
    calcium_results, calcium_results_path, ...
    dual_results, dual_results_path, ...
    stim_results, stim_results_path, ...
    time_frequency_results, time_frequency_results_path);

fprintf('Dual_analysis3 finished.\n');
fprintf('Results saved to: %s\n', save_path);

function [poolObj, poolReady, poolInfo] = initialize_dual_parallel_pool(poolInfo)
poolObj = [];
poolReady = false;

if exist('gcp', 'file') ~= 2
    poolInfo.message = "Parallel Computing Toolbox is not available.";
    warning('Dual_analysis3:ParallelUnavailable', '%s Continuing in serial mode.', poolInfo.message);
    return;
end

try
    poolObj = gcp('nocreate');
    if isempty(poolObj)
        fprintf('Opening parallel pool before movie loading...\n');
        poolObj = gcp;
    else
        fprintf('Reusing existing parallel pool before movie loading.\n');
    end

    poolReady = ~isempty(poolObj);
    poolInfo.ready = poolReady;
    if poolReady && isprop(poolObj, 'NumWorkers')
        poolInfo.num_workers = poolObj.NumWorkers;
        poolInfo.message = sprintf('Parallel pool ready with %d workers.', poolObj.NumWorkers);
    elseif poolReady
        poolInfo.message = "Parallel pool ready.";
    else
        poolInfo.message = "Parallel pool was requested but is empty.";
    end
    fprintf('%s\n', poolInfo.message);
catch ME
    poolInfo.ready = false;
    poolInfo.num_workers = 0;
    poolInfo.message = string(ME.message);
    warning('Dual_analysis3:ParpoolUnavailable', ...
        ['Parallel pool could not be started (%s). ' ...
        'Continuing in serial mode.'], ME.message);
end
end

function [cycle_manifest, record_manifest, camera_source, input_layout_info] = resolve_dual_camera_sources(cycle_path, camera_cfg, camera_source_override, raw_dual_input_cfg)
% Resolve where each camera movie came from and which label belongs to it.
% The intent is to keep analysis inputs tied to acquisition metadata when
% available, while still supporting older folder-only datasets.
cycle_manifest = [];
record_manifest = [];
camera_source = repmat(struct('path', '', 'label', "", 'original_frame_size', [NaN, NaN]), 1, numel(camera_cfg));
input_layout_info = struct( ...
    'mode', "rebuild_or_legacy_cycle_folder", ...
    'description', "Resolve from manifests first, then Cam*_ folders.", ...
    'primary_path', "", ...
    'green_path', "", ...
    'green_suffix', "");

if nargin < 4 || isempty(raw_dual_input_cfg)
    raw_dual_input_cfg = struct();
end
raw_dual_input_cfg = normalize_raw_dual_input_cfg(raw_dual_input_cfg);

if nargin >= 3 && ~isempty(camera_source_override)
    camera_source = normalize_camera_source_override(camera_source_override, camera_cfg);
    input_layout_info.mode = "explicit_camera_source_override";
    input_layout_info.description = "Caller provided camera_source_override directly.";
    return;
end

[raw_pair_source, raw_pair_found] = try_resolve_raw_dual_pair_sources(cycle_path, camera_cfg, raw_dual_input_cfg);
if raw_pair_found
    camera_source = raw_pair_source;
    input_layout_info.mode = "paired_raw_channel_folders";
    input_layout_info.description = "Resolved from one primary raw folder plus one sibling green raw folder.";
    input_layout_info.primary_path = string(camera_source(find(strcmpi(string({camera_cfg.role}), string(raw_dual_input_cfg.primary_role)), 1, 'first')).path);
    input_layout_info.green_path = string(camera_source(find(strcmpi(string({camera_cfg.role}), string(raw_dual_input_cfg.green_role)), 1, 'first')).path);
    input_layout_info.green_suffix = string(raw_dual_input_cfg.green_suffix);
    return;
end

cycle_manifest_file = fullfile(cycle_path, 'cycle_manifest.mat');
record_manifest_file = fullfile(fileparts(cycle_path), 'record_manifest.mat');

if isfile(cycle_manifest_file)
    tmp = load(cycle_manifest_file);
    if isfield(tmp, 'manifest')
        cycle_manifest = tmp.manifest;
    end
end

if isfile(record_manifest_file)
    tmp = load(record_manifest_file);
    if isfield(tmp, 'manifest')
        record_manifest = tmp.manifest;
    end
end

for i = 1:numel(camera_cfg)
    cam_idx = camera_cfg(i).camera_index;
    cam_path = '';
    cam_label = "";

    if ~isempty(cycle_manifest) && isfield(cycle_manifest, 'actual') ...
            && isfield(cycle_manifest.actual, 'movie_paths')
        movie_paths = cycle_manifest.actual.movie_paths;
        if numel(movie_paths) >= cam_idx
            if iscell(movie_paths)
                candidate_path = string(movie_paths{cam_idx});
            else
                candidate_path = string(movie_paths(cam_idx));
            end
            candidate_path = char(candidate_path);
            if exist(candidate_path, 'file') || exist(candidate_path, 'dir')
                cam_path = candidate_path;
            end
        end
    end

    if cam_label == "" && ~isempty(cycle_manifest) && isfield(cycle_manifest, 'spec') ...
            && isfield(cycle_manifest.spec, 'labels')
        labels = cycle_manifest.spec.labels;
        if numel(labels) >= cam_idx
            if iscell(labels)
                cam_label = string(labels{cam_idx});
            else
                cam_label = string(labels(cam_idx));
            end
        end
    end

    if isempty(cam_path)
        % Fall back to folder scanning only when manifest information is
        % missing, so older recordings remain usable.
        listing = dir(fullfile(cycle_path, sprintf('Cam%d_*', cam_idx)));
        listing = listing(~startsWith({listing.name}, '.'));
        if isempty(listing)
            error('Cannot resolve input for camera %d under %s', cam_idx, cycle_path);
        end
        cam_path = fullfile(listing(1).folder, listing(1).name);
        if cam_label == ""
            cam_label = string(erase(listing(1).name, sprintf('Cam%d_', cam_idx)));
        end
    end

    camera_source(i).path = cam_path;
    camera_source(i).label = cam_label;
    camera_source(i).original_frame_size = [NaN, NaN];
end
end

function raw_dual_input_cfg = normalize_raw_dual_input_cfg(raw_dual_input_cfg)
if ~isstruct(raw_dual_input_cfg)
    error('raw_dual_input_cfg must be a struct when provided.');
end
if ~isfield(raw_dual_input_cfg, 'green_suffix') || strlength(string(raw_dual_input_cfg.green_suffix)) == 0
    raw_dual_input_cfg.green_suffix = "_Green";
end
if ~isfield(raw_dual_input_cfg, 'primary_role') || strlength(string(raw_dual_input_cfg.primary_role)) == 0
    raw_dual_input_cfg.primary_role = "voltage";
end
if ~isfield(raw_dual_input_cfg, 'green_role') || strlength(string(raw_dual_input_cfg.green_role)) == 0
    raw_dual_input_cfg.green_role = "calcium";
end
end

function [camera_source, found_pair] = try_resolve_raw_dual_pair_sources(cycle_path, camera_cfg, raw_dual_input_cfg)
camera_source = repmat(struct('path', '', 'label', "", 'original_frame_size', [NaN, NaN]), 1, numel(camera_cfg));
found_pair = false;

if ~isfolder(cycle_path)
    return;
end

folder_name = string(get_last_path_part(cycle_path));
parent_dir = fileparts(cycle_path);
green_suffix = string(raw_dual_input_cfg.green_suffix);
primary_path = "";
green_path = "";

if strlength(folder_name) == 0 || strlength(green_suffix) == 0
    return;
end

if endsWith(folder_name, green_suffix, 'IgnoreCase', true)
    green_path = string(cycle_path);
    primary_name = extractBefore(folder_name, strlength(folder_name) - strlength(green_suffix) + 1);
    primary_candidate = fullfile(parent_dir, char(primary_name));
    if isfolder(primary_candidate)
        primary_path = string(primary_candidate);
    else
        return;
    end
else
    green_candidate = fullfile(parent_dir, char(folder_name + green_suffix));
    if isfolder(green_candidate)
        primary_path = string(cycle_path);
        green_path = string(green_candidate);
    else
        return;
    end
end

role_names = string({camera_cfg.role});
primary_role_idx = find(strcmpi(role_names, string(raw_dual_input_cfg.primary_role)), 1, 'first');
green_role_idx = find(strcmpi(role_names, string(raw_dual_input_cfg.green_role)), 1, 'first');
if isempty(primary_role_idx) || isempty(green_role_idx)
    error(['Raw dual-folder pairing requires camera_cfg to contain roles matching ' ...
        'raw_dual_input_cfg.primary_role and raw_dual_input_cfg.green_role.']);
end

camera_source(primary_role_idx).path = char(primary_path);
camera_source(primary_role_idx).label = get_last_path_part(char(primary_path));
camera_source(primary_role_idx).original_frame_size = [NaN, NaN];

camera_source(green_role_idx).path = char(green_path);
camera_source(green_role_idx).label = get_last_path_part(char(green_path));
camera_source(green_role_idx).original_frame_size = [NaN, NaN];

for idx = 1:numel(camera_source)
    if strlength(string(camera_source(idx).path)) == 0
        error('Raw dual-folder pairing did not assign a source path for camera_cfg entry %d.', idx);
    end
end

found_pair = true;
end

function part_name = get_last_path_part(input_path)
[~, part_name, ext_name] = fileparts(char(string(input_path)));
part_name = [part_name, ext_name];
end

function camera_source = normalize_camera_source_override(camera_source_override, camera_cfg)
% Normalize explicit camera-path overrides into the same struct layout used
% by rebuilt inputs. This keeps the rest of the analysis blind to where
% the movies originally came from.
camera_source = repmat(struct('path', '', 'label', "", 'original_frame_size', [NaN, NaN]), 1, numel(camera_cfg));

if ~isstruct(camera_source_override) || isempty(camera_source_override)
    error('camera_source_override must be a non-empty struct array.');
end

for i = 1:numel(camera_cfg)
    cam_idx = camera_cfg(i).camera_index;
    match_idx = [];

    if isfield(camera_source_override, 'camera_index')
        match_idx = find([camera_source_override.camera_index] == cam_idx, 1, 'first');
    end
    if isempty(match_idx) && numel(camera_source_override) >= i
        match_idx = i;
    end
    if isempty(match_idx)
        error('camera_source_override is missing camera %d.', cam_idx);
    end

    source_i = camera_source_override(match_idx);
    if ~isfield(source_i, 'path') || strlength(string(source_i.path)) == 0
        error('camera_source_override for camera %d is missing a movie path.', cam_idx);
    end

    camera_source(i).path = char(string(source_i.path));
    if ~exist(camera_source(i).path, 'file')
        error('Explicit movie path does not exist for camera %d: %s', cam_idx, camera_source(i).path);
    end

    if isfield(source_i, 'label') && strlength(string(source_i.label)) > 0
        camera_source(i).label = string(source_i.label);
    else
        [~, fallback_name, fallback_ext] = fileparts(camera_source(i).path);
        camera_source(i).label = string([fallback_name, fallback_ext]);
    end

    if isfield(source_i, 'original_frame_size') && numel(source_i.original_frame_size) == 2
        camera_source(i).original_frame_size = double(source_i.original_frame_size);
    end
end
end

function [movie_3d, ncols, nrows, nframes, original_ncols, original_nrows] = load_camera_movie(movie_path, do_transpose)
% Load one camera movie and apply only the orientation normalization needed
% for consistent later ROI analysis. No biological processing happens yet.
[movie_loaded, ncols_raw, nrows_raw, nframes] = load_movie(movie_path);
original_ncols = ncols_raw;
original_nrows = nrows_raw;

if ismatrix(movie_loaded)
    movie_3d = reshape(movie_loaded, ncols_raw, nrows_raw, []);
else
    movie_3d = movie_loaded;
end

if do_transpose
    % Each camera keeps its own transpose flag because the two optical
    % paths may not share the same orientation.
    movie_3d = pagetranspose(movie_3d);
end

[ncols, nrows, nframes] = size(movie_3d);
end

function save_base_dir = resolve_dual_analysis_save_base_dir(cycle_path, camera_cfg, raw_dual_input_cfg)
save_base_dir = cycle_path;
try
    [~, raw_pair_found] = try_resolve_raw_dual_pair_sources(cycle_path, camera_cfg, raw_dual_input_cfg);
    if raw_pair_found
        parent_dir = fileparts(cycle_path);
        if strlength(string(parent_dir)) > 0
            save_base_dir = parent_dir;
        end
    end
catch
    % Save-path selection should stay conservative. If raw-pair detection
    % fails here, fall back to the historical cycle_path-local save root
    % and let the later input-resolution stage raise the real error.
end
end

function [mask_voltage, mask_calcium, note_text] = resolve_reuse_dual_masks(roi_cache, offset_xy, correct_offset_mode)
mask_voltage = [];
mask_calcium = [];
note_text = "";

if isfield(roi_cache, 'rois') && isstruct(roi_cache.rois)
    rois_cache = roi_cache.rois;
else
    rois_cache = struct();
end

if isfield(rois_cache, 'bwmask') && ~isempty(rois_cache.bwmask)
    mask_voltage = rois_cache.bwmask;
elseif isfield(roi_cache, 'bwmask') && ~isempty(roi_cache.bwmask)
    mask_voltage = roi_cache.bwmask;
end

if isfield(rois_cache, 'bwmask_ca') && ~isempty(rois_cache.bwmask_ca)
    mask_calcium = rois_cache.bwmask_ca;
elseif isfield(roi_cache, 'bwmask_ca') && ~isempty(roi_cache.bwmask_ca)
    mask_calcium = roi_cache.bwmask_ca;
end

if isempty(mask_voltage) && isempty(mask_calcium)
    return;
end

if isempty(mask_voltage) && ~isempty(mask_calcium)
    if has_valid_dual_offset(offset_xy)
        mask_voltage = translate_label_mask_by_offset(mask_calcium, offset_xy(1), offset_xy(2));
        note_text = "Reuse ROI fallback: only calcium mask was found, so voltage mask was reconstructed from offset.";
    else
        mask_voltage = mask_calcium;
        note_text = "Reuse ROI fallback: only calcium mask was found, so it was copied to voltage.";
    end
    return;
end

if ~isempty(mask_voltage) && isempty(mask_calcium)
    if has_valid_dual_offset(offset_xy)
        mask_calcium = translate_label_mask_by_offset(mask_voltage, -offset_xy(1), -offset_xy(2));
        if strcmpi(string(correct_offset_mode), "matlab_register")
            note_text = ['Reuse ROI fallback: old ROI file only contains one bwmask. ' ...
                'bwmask_ca was reconstructed from the matlab_register offset without moving the movie data.'];
        elseif strcmpi(string(correct_offset_mode), "manual_points")
            note_text = ['Reuse ROI fallback: old ROI file only contains one bwmask. ' ...
                'bwmask_ca was reconstructed from the manually selected offset without moving the movie data.'];
        else
            note_text = ['Reuse ROI fallback: old ROI file only contains one bwmask. ' ...
                'bwmask_ca was reconstructed from the available offset.'];
        end
    elseif strcmpi(string(correct_offset_mode), "manual_points")
        note_text = ['Reuse ROI fallback: old ROI file only contains one bwmask. ' ...
            'manual_points mode will estimate a fresh offset first, then rebuild bwmask_ca from the saved voltage ROI mask.'];
    else
        error(['The provided reuse_roi_file only contains one mask (bwmask) and no bwmask_ca. ' ...
            'Provide a valid offset or a dual ROI file that contains both masks.']);
    end
end
end

function tf = has_valid_dual_offset(offset_xy)
tf = ~isempty(offset_xy) && isnumeric(offset_xy) && numel(offset_xy) == 2 && all(isfinite(offset_xy));
end

function translated_mask = translate_label_mask_by_offset(mask_in, dx, dy)
mask_in = double(mask_in);
if exist('imtranslate', 'file') == 2
    translated_mask = imtranslate(mask_in, [dx, dy], 'nearest', 'FillValues', 0);
else
    translated_mask = shift_label_mask_integer(mask_in, dx, dy);
end
translated_mask = round(translated_mask);
translated_mask(translated_mask < 0) = 0;
end

function shifted_mask = shift_label_mask_integer(mask_in, dx, dy)
shift_x = round(dx);
shift_y = round(dy);
shifted_mask = zeros(size(mask_in));
[nrows_mask, ncols_mask] = size(mask_in);

src_rows = max(1, 1 - shift_y):min(nrows_mask, nrows_mask - shift_y);
src_cols = max(1, 1 - shift_x):min(ncols_mask, ncols_mask - shift_x);
dst_rows = src_rows + shift_y;
dst_cols = src_cols + shift_x;
if isempty(src_rows) || isempty(src_cols)
    return;
end
shifted_mask(dst_rows, dst_cols) = mask_in(src_rows, src_cols);
end

function [registration_info, roi_offset_mode_for_selection] = estimate_dual_channel_registration_offset( ...
    movie_voltage, movie_calcium, nrows, ncols, correct_offset_mode, save_path)
registration_info = struct( ...
    'mode', string(correct_offset_mode), ...
    'applied', false, ...
    'offset_xy', [], ...
    'preview_fig', '', ...
    'preview_png', '', ...
    'note', '');
roi_offset_mode_for_selection = string(correct_offset_mode);

if strcmpi(string(correct_offset_mode), "none")
    registration_info.note = 'No automatic registration offset estimated in this mode.';
    return;
end

mean_voltage = reshape(mean(double(movie_voltage), 2), ncols, nrows);
mean_calcium = reshape(mean(double(movie_calcium), 2), ncols, nrows);

switch lower(string(correct_offset_mode))
    case "manual_points"
        [offset_xy, preview_fig, preview_png] = estimate_dual_offset_from_manual_points_main( ...
            mean_voltage, mean_calcium, save_path);
        registration_info.applied = true;
        registration_info.offset_xy = offset_xy;
        registration_info.preview_fig = preview_fig;
        registration_info.preview_png = preview_png;
        registration_info.note = ['Offset was picked manually from the two average images before ROI selection. ' ...
            'Movie data are not moved here.'];
    case "matlab_register"
        [offset_xy, preview_fig, preview_png] = estimate_dual_offset_with_matlab_register_main(mean_voltage, mean_calcium, save_path);
        registration_info.applied = true;
        registration_info.offset_xy = offset_xy;
        registration_info.preview_fig = preview_fig;
        registration_info.preview_png = preview_png;
        registration_info.note = 'Offset maps calcium coordinates into voltage coordinates. Movie data are not moved here.';
    otherwise
        error('Unsupported correct_offset_mode during channel registration: %s', char(string(correct_offset_mode)));
end

% Offset estimation is handled at the script level so ROI selection can
% stay focused on ROI replay/drawing instead of opening a second point-pick
% path inside select_ROI_dual.
roi_offset_mode_for_selection = "none";
end

function [offset_xy, fig_file, png_file] = estimate_dual_offset_from_manual_points_main(mean_voltage, mean_calcium, save_path)
if ~usejava('desktop')
    error(['manual_points requires a MATLAB Desktop session because the offset is picked ' ...
        'directly on the voltage and calcium average images.']);
end

display_voltage = normalize_dual_registration_image_main(mean_voltage);
display_calcium = normalize_dual_registration_image_main(mean_calcium);

fig = figure('Color', 'w', 'Name', 'Dual Manual Offset Selection', ...
    'Position', [100, 100, 1200, 520]);
ax_voltage = subplot(1, 2, 1, 'Parent', fig);
imshow(display_voltage, 'Parent', ax_voltage);
title(ax_voltage, 'Voltage Average: click reference point');
ax_calcium = subplot(1, 2, 2, 'Parent', fig);
imshow(display_calcium, 'Parent', ax_calcium);
title(ax_calcium, 'Calcium Average: click matching point');

disp(['Manual dual offset estimation: first click one point on the Voltage Average image, ' ...
    'then click the matching point on the Calcium Average image.']);

axes(ax_voltage); %#ok<LAXES>
[x_voltage, y_voltage] = ginput(1);
hold(ax_voltage, 'on');
plot(ax_voltage, x_voltage, y_voltage, 'ro', 'MarkerSize', 10, 'LineWidth', 1.5);
text(ax_voltage, x_voltage, y_voltage, ' Voltage point', 'Color', [0.8 0.1 0.1], ...
    'FontWeight', 'bold', 'Interpreter', 'none');

axes(ax_calcium); %#ok<LAXES>
[x_calcium, y_calcium] = ginput(1);
hold(ax_calcium, 'on');
plot(ax_calcium, x_calcium, y_calcium, 'go', 'MarkerSize', 10, 'LineWidth', 1.5);
text(ax_calcium, x_calcium, y_calcium, ' Calcium point', 'Color', [0.1 0.6 0.1], ...
    'FontWeight', 'bold', 'Interpreter', 'none');

offset_xy = [x_voltage - x_calcium, y_voltage - y_calcium];
sgtitle(fig, sprintf('Manual offset [x y] = [%.3f %.3f]', offset_xy(1), offset_xy(2)));

fig_file = fullfile(save_path, '0_dual_manual_points_preview.fig');
png_file = fullfile(save_path, '0_dual_manual_points_preview.png');
save_figure_bundle_preserve_layout(fig, fig_file, png_file);
close(fig);
end

function [offset_xy, fig_file, png_file] = estimate_dual_offset_with_matlab_register_main(mean_voltage, mean_calcium, save_path)
if exist('imregconfig', 'file') ~= 2 || exist('imregtform', 'file') ~= 2
    error(['MATLAB built-in registration requires imregconfig and imregtform ' ...
        '(Image Processing Toolbox).']);
end

fixed_image = normalize_dual_registration_image_main(mean_voltage);
moving_image = normalize_dual_registration_image_main(mean_calcium);
[optimizer, metric] = imregconfig('multimodal');
tform = imregtform(moving_image, fixed_image, 'translation', optimizer, metric);
offset_xy = [double(tform.T(3, 1)), double(tform.T(3, 2))];

registered_moving = apply_xy_shift_to_image(moving_image, offset_xy);
fig = figure('Color', 'w', 'Name', 'Dual MATLAB Register Preview');
subplot(1, 3, 1); imshowpair(fixed_image, moving_image); title('Voltage vs Raw Calcium');
subplot(1, 3, 2); imshowpair(fixed_image, registered_moving); title(sprintf('Voltage vs Registered Calcium [%.2f %.2f]', offset_xy(1), offset_xy(2)));
subplot(1, 3, 3); imagesc(cat(3, fixed_image, registered_moving, zeros(size(fixed_image)))); axis image off; title('Color Merge Preview');
fig_file = fullfile(save_path, '0_dual_matlab_register_preview.fig');
png_file = fullfile(save_path, '0_dual_matlab_register_preview.png');
save_figure_bundle_preserve_layout(fig, fig_file, png_file);
close(fig);
end

function image_out = normalize_dual_registration_image_main(image_in)
image_out = double(image_in);
image_out(~isfinite(image_out)) = 0;
if exist('mat2gray', 'file') == 2
    image_out = mat2gray(image_out);
else
    finite_values = image_out(isfinite(image_out));
    if isempty(finite_values)
        image_out = zeros(size(image_out));
    else
        min_value = min(finite_values);
        max_value = max(finite_values);
        if max_value > min_value
            image_out = (image_out - min_value) ./ (max_value - min_value);
        else
            image_out = zeros(size(image_out));
        end
    end
end
end

function shifted_image = apply_xy_shift_to_image(image_in, offset_xy)
if exist('imtranslate', 'file') == 2
    shifted_image = imtranslate(double(image_in), offset_xy, 'FillValues', 0);
else
    shifted_image = shift_label_mask_integer(double(image_in), offset_xy(1), offset_xy(2));
end
end

function [rgb_preview, info] = save_dual_average_color_merge( ...
    movie_voltage, movie_calcium, nrows, ncols, offset_xy, ...
    tif_path, png_path, tif_roi_path, png_roi_path, rois)
mean_voltage = reshape(mean(double(movie_voltage), 2), ncols, nrows);
mean_calcium = reshape(mean(double(movie_calcium), 2), ncols, nrows);
mean_calcium_aligned = mean_calcium;
if has_valid_dual_offset(offset_xy)
    mean_calcium_aligned = apply_xy_shift_to_image(mean_calcium, offset_xy);
end

voltage_uint16 = scale_image_to_uint16_display(mean_voltage);
calcium_uint16 = scale_image_to_uint16_display(mean_calcium_aligned);
rgb_preview = cat(3, voltage_uint16, calcium_uint16, zeros(size(voltage_uint16), 'uint16'));
imwrite(rgb_preview, tif_path, 'tif', 'Compression', 'none');
imwrite(rgb_preview, png_path, 'png');

rgb_with_roi = rgb_preview;
if nargin >= 10 && isstruct(rois) && isfield(rois, 'boundary') && ~isempty(rois.boundary)
    rgb_with_roi = draw_roi_boundaries_on_rgb_merge(rgb_with_roi, rois.boundary);
    imwrite(rgb_with_roi, tif_roi_path, 'tif', 'Compression', 'none');
    imwrite(rgb_with_roi, png_roi_path, 'png');
end

info = struct( ...
    'tif_file', tif_path, ...
    'png_file', png_path, ...
    'tif_roi_file', tif_roi_path, ...
    'png_roi_file', png_roi_path, ...
    'offset_xy_used_for_alignment', offset_xy, ...
    'alignment_rule', 'calcium average image shifted by offset into voltage coordinates before RGB merge', ...
    'roi_overlay_rule', 'Voltage ROI boundaries are overlaid on the aligned RGB merge image', ...
    'created_at', datetime("now"));
end

function rgb_out = draw_roi_boundaries_on_rgb_merge(rgb_in, boundary_cells)
rgb_out = rgb_in;
if isempty(boundary_cells)
    return;
end

[image_height, image_width, ~] = size(rgb_out);
overlay_color = uint16(65535);
for roi_idx = 1:numel(boundary_cells)
    boundary = boundary_cells{roi_idx};
    if isempty(boundary) || size(boundary, 2) < 2
        continue;
    end
    row_idx = round(boundary(:, 1));
    col_idx = round(boundary(:, 2));
    valid = row_idx >= 1 & row_idx <= image_height & col_idx >= 1 & col_idx <= image_width;
    row_idx = row_idx(valid);
    col_idx = col_idx(valid);
    if isempty(row_idx)
        continue;
    end
    linear_idx = sub2ind([image_height, image_width], row_idx, col_idx);
    for channel_idx = 1:3
        channel_plane = rgb_out(:, :, channel_idx);
        channel_plane(linear_idx) = overlay_color;
        rgb_out(:, :, channel_idx) = channel_plane;
    end
end
end

function [saved_ok, info] = save_dual_map_with_roi(map_voltage, map_calcium, rois, fig_path, png_path)
saved_ok = false;
info = struct();
if isempty(map_voltage) || isempty(map_calcium) || ~isstruct(rois) ...
        || ~isfield(rois, 'boundary') || ~isfield(rois, 'boundary_ca')
    return;
end

fig = figure('Color', 'w');
ax1 = subplot(1, 2, 1, 'Parent', fig);
imagesc(ax1, map_voltage);
axis(ax1, 'image');
title(ax1, 'Voltage Sensitivity Map With ROI');
colorbar(ax1);
hold(ax1, 'on');
overlay_roi_boundaries(ax1, rois.boundary, 'w');

ax2 = subplot(1, 2, 2, 'Parent', fig);
imagesc(ax2, map_calcium);
axis(ax2, 'image');
title(ax2, 'Calcium Sensitivity Map With ROI');
colorbar(ax2);
hold(ax2, 'on');
overlay_roi_boundaries(ax2, rois.boundary_ca, 'w');

save_figure_bundle_preserve_layout(fig, fig_path, png_path);
close(fig);
saved_ok = true;
info = struct( ...
    'fig_file', fig_path, ...
    'png_file', png_path, ...
    'overlay_rule', 'Voltage ROI boundaries on voltage map and calcium ROI boundaries on calcium map', ...
    'created_at', datetime("now"));
end

function overlay_roi_boundaries(ax, boundary_cells, edge_color)
if isempty(boundary_cells)
    return;
end
for roi_idx = 1:numel(boundary_cells)
    boundary = boundary_cells{roi_idx};
    if isempty(boundary) || size(boundary, 2) < 2
        continue;
    end
    plot(ax, boundary(:, 2), boundary(:, 1), 'Color', edge_color, 'LineWidth', 0.8);
end
end

function image_uint16 = scale_image_to_uint16_display(image_in)
image_in = double(image_in);
image_in(~isfinite(image_in)) = 0;
low_q = quantile(image_in(:), 0.0005);
high_q = quantile(image_in(:), 0.9995);
if ~isfinite(low_q) || ~isfinite(high_q) || high_q <= low_q
    image_uint16 = zeros(size(image_in), 'uint16');
    return;
end
image_scaled = (image_in - low_q) ./ (high_q - low_q);
image_scaled(image_scaled < 0) = 0;
image_scaled(image_scaled > 1) = 1;
image_uint16 = uint16(round(image_scaled * 65535));
end

function correct_offset_mode = normalize_correct_offset_mode_dual(correct_offset_mode)
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

function resolution = resolve_dual_reuse_roi_source(reuse_results_path, reuse_roi_file)
resolution = struct( ...
    'effective_roi_file', "", ...
    'source', "none", ...
    'results_candidate', "", ...
    'message', "");

results_path_text = strtrim(string(reuse_results_path));
roi_file_text = strtrim(string(reuse_roi_file));

results_candidate = "";
if strlength(results_path_text) > 0
    if isfolder(char(results_path_text))
        results_candidate = fullfile(char(results_path_text), '1_dual_roi_results.mat');
    elseif isfile(char(results_path_text))
        [~, results_name, results_ext] = fileparts(char(results_path_text));
        if strcmpi([results_name results_ext], '1_dual_roi_results.mat')
            results_candidate = results_path_text;
        end
    end
end
resolution.results_candidate = string(results_candidate);

results_candidate_exists = strlength(string(results_candidate)) > 0 && isfile(char(results_candidate));
roi_file_exists = strlength(roi_file_text) > 0 && isfile(char(roi_file_text));
has_conflict = results_candidate_exists && roi_file_exists ...
    && ~strcmpi(char(string(results_candidate)), char(roi_file_text));

if results_candidate_exists
    resolution.effective_roi_file = string(results_candidate);
    resolution.source = "reuse_results_path";
    if has_conflict
        resolution.message = sprintf(['ROI reuse conflict: reuse_results_path has higher priority, so this run will use\n' ...
            '  %s\ninstead of\n  %s'], char(string(results_candidate)), char(roi_file_text));
    elseif roi_file_exists
        resolution.message = sprintf('ROI reuse: both reuse_results_path and reuse_roi_file resolve to the same file: %s', ...
            char(string(results_candidate)));
    elseif strlength(results_path_text) > 0
        resolution.message = sprintf('ROI reuse: using 1_dual_roi_results.mat from reuse_results_path: %s', ...
            char(string(results_candidate)));
    end
    return;
end

if strlength(results_path_text) > 0 && strlength(string(results_candidate)) > 0 && ~results_candidate_exists
    if roi_file_exists
        resolution.effective_roi_file = roi_file_text;
        resolution.source = "reuse_roi_file_fallback";
        resolution.message = sprintf(['ROI reuse warning: reuse_results_path was given higher priority, but its expected ROI file was not found:\n' ...
            '  %s\nFalling back to reuse_roi_file:\n  %s'], ...
            char(string(results_candidate)), char(roi_file_text));
        return;
    else
        resolution.message = sprintf(['ROI reuse warning: reuse_results_path was provided, but its expected ROI file was not found:\n' ...
            '  %s'], char(string(results_candidate)));
        return;
    end
end

if roi_file_exists
    resolution.effective_roi_file = roi_file_text;
    resolution.source = "reuse_roi_file";
    resolution.message = sprintf('ROI reuse: using explicit reuse_roi_file: %s', char(roi_file_text));
    return;
end

if strlength(roi_file_text) > 0
    resolution.effective_roi_file = roi_file_text;
    resolution.source = "reuse_roi_file_missing";
    resolution.message = sprintf('ROI reuse warning: explicit reuse_roi_file was provided but not found: %s', char(roi_file_text));
end
end

function [voltage_data, calcium_data, geometry_info] = match_camera_geometry(voltage_data, calcium_data)
% Force both channels onto one common spatial grid before any shared ROI
% selection or dual-channel trace comparison.
geometry_info = struct();
geometry_info.voltage_before = [voltage_data.ncols, voltage_data.nrows];
geometry_info.calcium_before = [calcium_data.ncols, calcium_data.nrows];
geometry_info.action = "none";

if voltage_data.ncols == calcium_data.ncols && voltage_data.nrows == calcium_data.nrows
    geometry_info.common_size = [voltage_data.ncols, voltage_data.nrows];
    return;
end

if mod(calcium_data.ncols, voltage_data.ncols) == 0 && mod(calcium_data.nrows, voltage_data.nrows) == 0
    ratio_x = calcium_data.ncols / voltage_data.ncols;
    ratio_y = calcium_data.nrows / voltage_data.nrows;
    if ratio_x ~= ratio_y
        error('Calcium resize ratios are inconsistent across x and y.');
    end
    calcium_data.movie_3d = average_pool_movie(calcium_data.movie_3d, ratio_x);
    [calcium_data.ncols, calcium_data.nrows, calcium_data.nframes] = size(calcium_data.movie_3d);
    geometry_info.action = "downsample_calcium_to_voltage";
elseif mod(voltage_data.ncols, calcium_data.ncols) == 0 && mod(voltage_data.nrows, calcium_data.nrows) == 0
    ratio_x = voltage_data.ncols / calcium_data.ncols;
    ratio_y = voltage_data.nrows / calcium_data.nrows;
    if ratio_x ~= ratio_y
        error('Voltage resize ratios are inconsistent across x and y.');
    end
    voltage_data.movie_3d = average_pool_movie(voltage_data.movie_3d, ratio_x);
    [voltage_data.ncols, voltage_data.nrows, voltage_data.nframes] = size(voltage_data.movie_3d);
    geometry_info.action = "downsample_voltage_to_calcium";
else
    error('Voltage and calcium camera sizes are incompatible: [%d %d] vs [%d %d].', ...
        voltage_data.ncols, voltage_data.nrows, calcium_data.ncols, calcium_data.nrows);
end

geometry_info.common_size = [voltage_data.ncols, voltage_data.nrows];
end

function movie_out = average_pool_movie(movie_in, ratio)
% Use simple average pooling when one channel needs to be downsampled to
% match the spatial sampling of the other channel.
ratio = round(ratio);
[ncols, nrows, nframes] = size(movie_in);
movie_out = reshape(movie_in, ratio, ncols / ratio, ratio, nrows / ratio, nframes);
movie_out = squeeze(mean(mean(movie_out, 1), 3));
movie_out = reshape(movie_out, ncols / ratio, nrows / ratio, nframes);
end

function movie_info = make_movie_info(camera_data)
movie_info = struct( ...
    'role', camera_data.role, ...
    'camera_index', camera_data.camera_index, ...
    'label', camera_data.label, ...
    'source_path', camera_data.source_path, ...
    'frame_rate', camera_data.frame_rate, ...
    'frame_count', camera_data.nframes, ...
    'original_frame_size', camera_data.original_frame_size, ...
    'analysis_frame_size', [camera_data.ncols, camera_data.nrows], ...
    'transpose_before_analysis', camera_data.transpose_before_analysis, ...
    'motion', struct( ...
        'applied', false, ...
        'method', '', ...
        'shift_file', '', ...
        'parameter_file', ''), ...
    'created_at', datetime("now"), ...
    'updated_at', datetime("now"));
end

function t_axis = build_time_axis_from_movie_info(movie_info)
t_axis = (1:movie_info.frame_count)' / movie_info.frame_rate;
end

function results = store_trace_stage(results, stage_name, data, parent_results, roi_file, movie_info, method, parameters)
% Save one named trace stage together with the minimum provenance needed
% to understand how that stage was produced.
%
% Natural-language meaning:
%   this is the common bookkeeping step used throughout the script after a
%   new trace representation has been computed.
%
% Interface example:
%   results.trace_results.snr.data
%   results.trace_results.snr.info.parent_results   -> {'bleach_removed','noise'}
%   results.trace_results.snr.info.method           -> 'signal_divided_by_noise_std'
results.trace_results.(stage_name).data = data;
results.trace_results.(stage_name).info = struct( ...
    'result_name', stage_name, ...
    'parent_results', {parent_results}, ...
    'roi_file', roi_file, ...
    'movie_source', movie_info.source_path, ...
    'transpose_before_analysis', movie_info.transpose_before_analysis, ...
    'motion_applied', movie_info.motion.applied, ...
    'method', method, ...
    'parameters', parameters, ...
    'created_at', datetime("now"));
end

function copy_analysis_code(currentScript, code_path)
currentScript = char(string(currentScript));
if isempty(currentScript) || ~isfile(currentScript)
    warning('Dual_analysis3:CodeSnapshotSkipped', ...
        'Analysis code snapshot skipped because Dual_analysis3.m could not be resolved as a file: %s', currentScript);
    return;
end
try
    [requiredFiles, ~] = matlab.codetools.requiredFilesAndProducts(currentScript);
catch ME
    warning('Dual_analysis3:CodeDependencyScanFailed', ...
        'Dependency scan failed (%s). Copying only the main script.', ME.message);
    requiredFiles = {currentScript};
end
for k = 1:length(requiredFiles)
    [~, name, ext] = fileparts(requiredFiles{k});
    if isfile(requiredFiles{k})
        copyfile(requiredFiles{k}, fullfile(code_path, [name, ext]));
    end
end
fprintf('All required code files were copied to %s\n', code_path);
end

function [voltage_corrected, calcium_corrected, voltage_motion_info, calcium_motion_info] = run_shared_motion_correction(voltage_movie, calcium_movie, cycle_path, save_path, cfg)
% Estimate one shared motion model on voltage, then apply the same shifts
% to calcium so both channels remain spatially coupled after correction.
voltage_corrected = voltage_movie;
calcium_corrected = calcium_movie;

voltage_motion_info = struct( ...
    'applied', false, ...
    'method', '', ...
    'shift_file', '', ...
    'source_shift_file', '', ...
    'parameter_file', '', ...
    'downsample_factor', NaN, ...
    'downsampled_tif_file', '', ...
    'metrics_file', '', ...
    'metrics_fig', '', ...
    'metrics_png', '', ...
    'highpass', cfg.highpass, ...
    'use_saved_shift', cfg.use_saved_shift, ...
    'auto_reused_previous_shift', false, ...
    'shared_with_role', 'calcium');
calcium_motion_info = struct( ...
    'applied', false, ...
    'method', '', ...
    'shift_file', '', ...
    'source_shift_file', '', ...
    'parameter_file', '', ...
    'downsample_factor', NaN, ...
    'downsampled_tif_file', '', ...
    'highpass', cfg.highpass, ...
    'use_saved_shift', cfg.use_saved_shift, ...
    'auto_reused_previous_shift', false, ...
    'source_role', 'voltage');

if ~cfg.enabled
    return;
end

[ncols_v, nrows_v, nframes_v] = size(voltage_movie);
[ncols_c, nrows_c, nframes_c] = size(calcium_movie);
if ncols_v ~= ncols_c || nrows_v ~= nrows_c
    error('Shared motion correction requires matched voltage/calcium frame size before applying shifts.');
end
if nframes_v ~= nframes_c
    error('Shared motion correction requires the same frame count for voltage and calcium movies.');
end

options_r = NoRMCorreSetParms('d1', ncols_v, 'd2', nrows_v, 'bin_width', 200, ...
    'max_shift', 30, 'us_fac', 30, 'iter', 1, 'correct_bidir', false);

shift_res_path = fullfile(save_path, 'shared_motion_shifts_result.mat');
params_save_path = fullfile(save_path, 'shared_motion_correction_para.mat');
voltage_single = single(voltage_movie);
calcium_single = single(calcium_movie);

shift_source_file = string(cfg.saved_shift_file);
auto_reused_previous_shift = false;
if strlength(shift_source_file) == 0 && cfg.auto_reuse_previous_shift
    shift_source_file = find_latest_previous_motion_shift(cycle_path, save_path);
    auto_reused_previous_shift = strlength(shift_source_file) > 0;
end

if (cfg.use_saved_shift && strlength(shift_source_file) > 0) || auto_reused_previous_shift
    % Reusing a previously estimated shift field is useful when iterating
    % on later sections without recomputing motion every time.
    fprintf('Reusing existing motion shift file:\n  %s\n', shift_source_file);
    shift_data = load(char(shift_source_file));
    if isfield(shift_data, 'options_r')
        options_r = shift_data.options_r;
    end
    if ~isfield(shift_data, 'shifts_r')
        error('Saved motion shift file does not contain shifts_r: %s', shift_source_file);
    end
    shifts_r = shift_data.shifts_r;
    if ~strcmpi(char(shift_source_file), shift_res_path)
        copyfile(char(shift_source_file), shift_res_path);
    end
else
    % Only build the high-pass temporary movie when NoRMCorre actually has
    % to estimate shifts. When reusing a saved shift field, this full-size
    % temporary copy is unnecessary and can add tens of GB to the memory
    % peak before apply_shifts.
    if cfg.highpass
        movie_for_estimation = create_temp_highpass(voltage_single);
    else
        movie_for_estimation = voltage_single;
    end
    [~, shifts_r, ~] = normcorre_batch(movie_for_estimation, options_r);
    clear movie_for_estimation;
    save(shift_res_path, 'shifts_r', 'options_r', '-v7.3');
end

voltage_corrected = apply_shifts(voltage_single, shifts_r, options_r);
calcium_corrected = apply_shifts(calcium_single, shifts_r, options_r);
save(params_save_path, 'options_r', 'cfg');

% Save the same tsub=40 motion-corrected QC TIFF that AP_analysis2 /
% AP_analysis3 produce, but do it for both channels in the shared-motion
% dual workflow.
downsample_factor = 40;
voltage_ds_path = fullfile(save_path, sprintf('voltage_motion_corrected_ds%d.tif', downsample_factor));
calcium_ds_path = fullfile(save_path, sprintf('calcium_motion_corrected_ds%d.tif', downsample_factor));
save_downsampled_motion_tif(voltage_corrected, downsample_factor, voltage_ds_path);
save_downsampled_motion_tif(calcium_corrected, downsample_factor, calcium_ds_path);
clear calcium_single;

% Quantify motion-correction quality on a temporally downsampled voltage
% reference movie. Running NoRMCorre motion_metrics on the full dual-camera
% movie allocates large cropped copies and can exceed memory on long runs.
motion_metrics_file = fullfile(save_path, 'shared_motion_metrics.mat');
motion_metrics_fig = fullfile(save_path, 'shared_motion_metrics.fig');
motion_metrics_png = fullfile(save_path, 'shared_motion_metrics.png');
save_shared_motion_metrics( ...
    voltage_single, voltage_corrected, ...
    shifts_r, options_r.max_shift, ...
    motion_metrics_file, motion_metrics_fig, motion_metrics_png, downsample_factor);

voltage_motion_info.applied = true;
voltage_motion_info.method = 'NoRMCorre_rigid_shared_voltage_reference';
voltage_motion_info.shift_file = shift_res_path;
voltage_motion_info.source_shift_file = char(shift_source_file);
voltage_motion_info.parameter_file = params_save_path;
voltage_motion_info.downsample_factor = downsample_factor;
voltage_motion_info.downsampled_tif_file = voltage_ds_path;
voltage_motion_info.metrics_file = motion_metrics_file;
voltage_motion_info.metrics_fig = motion_metrics_fig;
voltage_motion_info.metrics_png = motion_metrics_png;
voltage_motion_info.use_saved_shift = cfg.use_saved_shift || auto_reused_previous_shift;
voltage_motion_info.auto_reused_previous_shift = auto_reused_previous_shift;

calcium_motion_info.applied = true;
calcium_motion_info.method = 'reuse_voltage_motion_shifts';
calcium_motion_info.shift_file = shift_res_path;
calcium_motion_info.source_shift_file = char(shift_source_file);
calcium_motion_info.parameter_file = params_save_path;
calcium_motion_info.downsample_factor = downsample_factor;
calcium_motion_info.downsampled_tif_file = calcium_ds_path;
calcium_motion_info.use_saved_shift = cfg.use_saved_shift || auto_reused_previous_shift;
calcium_motion_info.auto_reused_previous_shift = auto_reused_previous_shift;
end

function save_downsampled_motion_tif(movie_3d, downsample_factor, save_path_tif)
if exist('array2tif', 'file') ~= 2
    error('array2tif helper is required to save downsampled motion-corrected TIFFs.');
end
movie_ds = downsample_movie_time_mean(movie_3d, downsample_factor);
movie_ds = scale_movie_to_uint16_display(movie_ds);
array2tif(uint16(movie_ds), save_path_tif);
end

function movie_ds = downsample_movie_time_mean(movie_3d, downsample_factor)
% Downsample by averaging contiguous temporal blocks without reshaping the
% full movie. NoRMCorre's downsample_data is convenient but can allocate a
% very large temporary array for long dual-camera movies.
movie_size = size(movie_3d);
if numel(movie_size) < 3
    movie_size(3) = 1;
end
ncols = movie_size(1);
nrows = movie_size(2);
nframes = movie_size(3);
downsample_factor = max(1, round(double(downsample_factor)));
ndownsampled = floor(nframes / downsample_factor);
if ndownsampled < 1
    movie_ds = single(movie_3d);
    return;
end

movie_ds = zeros(ncols, nrows, ndownsampled, 'single');
for idx = 1:ndownsampled
    frame_start = (idx - 1) * downsample_factor + 1;
    frame_stop = idx * downsample_factor;
    movie_ds(:, :, idx) = mean(single(movie_3d(:, :, frame_start:frame_stop)), 3, 'omitnan');
end
end

function movie_scaled = scale_movie_to_uint16_display(movie_3d)
nn = quantile(movie_3d(:), 0.0005);
mm = quantile(movie_3d(:), 0.99995);
if ~isfinite(nn) || ~isfinite(mm) || mm <= nn
    movie_scaled = zeros(size(movie_3d), 'single');
    return;
end
movie_scaled = (single(movie_3d) - nn) / (mm - nn) * 65535;
movie_scaled(movie_scaled < 0) = 0;
movie_scaled(movie_scaled > 65535) = 65535;
end

function save_shared_motion_metrics(full_raw, full_corrected, shifts_r, max_shift, mat_file, fig_file, png_file, metric_downsample_factor)
if exist('motion_metrics', 'file') ~= 2
    error('NoRMCorre helper motion_metrics is required for motion-correction QC metrics.');
end
if nargin < 8 || isempty(metric_downsample_factor)
    metric_downsample_factor = 40;
end
metric_downsample_factor = max(1, round(double(metric_downsample_factor)));

full_raw_ds = downsample_movie_time_mean(full_raw, metric_downsample_factor);
[c_full_raw, m_full_raw, v_full_raw] = motion_metrics(full_raw_ds, max_shift);
clear full_raw_ds;

full_corrected_ds = downsample_movie_time_mean(full_corrected, metric_downsample_factor);
[c_full_corrected, m_full_corrected, v_full_corrected] = motion_metrics(full_corrected_ds, max_shift);
clear full_corrected_ds;

shifts_plot = squeeze(cat(3, shifts_r(:).shifts));
metrics_data = struct( ...
    'rigid_shifts', shifts_plot, ...
    'metric_downsample_factor', metric_downsample_factor, ...
    'filtered', struct( ...
        'status', 'skipped', ...
        'reason', 'High-pass filtered metrics are skipped to avoid an extra full-movie apply_shifts and motion_metrics memory peak.'), ...
    'full', struct( ...
        'correlation_raw', c_full_raw, ...
        'correlation_corrected', c_full_corrected, ...
        'mean_raw', m_full_raw, ...
        'mean_corrected', m_full_corrected, ...
        'variance_raw', v_full_raw, ...
        'variance_corrected', v_full_corrected));
save(mat_file, 'metrics_data', '-v7.3');

fig = figure('Color', 'w', 'Name', 'Shared Motion Metrics');
subplot(3, 1, 1);
plot(shifts_plot, 'LineWidth', 1);
title('Shared Rigid Shifts');
legend('y-shifts', 'x-shifts');
grid on;

subplot(3, 1, 2);
text(0.5, 0.5, sprintf('High-pass filtered metrics skipped\nFull-voltage metrics use time downsample x%d', metric_downsample_factor), ...
    'Units', 'normalized', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
axis off;

subplot(3, 1, 3);
plot(c_full_raw, 'Color', [0.45 0.45 0.45], 'LineWidth', 1);
hold on;
plot(c_full_corrected, 'r', 'LineWidth', 1.2);
title('Correlation Coefficients On Full Voltage Movie');
legend('raw', 'corrected');
ylim([0.8, 1]);
grid on;
xlabel(sprintf('Downsampled frame (x%d)', metric_downsample_factor));

save_figure_bundle(fig, fig_file, png_file);
close(fig);
end

function [traces_bleach_removed, baseline, parameters] = run_bleach_removal(traces_bg_removed, freq, t_axis, bleach_mode)
% Treat bleaching as a slow baseline component and separate it from the
% faster activity of interest. Different modes reflect different baseline
% assumptions.
switch lower(bleach_mode)
    case 'highpass'
        fc = 0.5 / t_axis(end);
        [traces_bleach_removed, baseline] = highpass_bleach_remove(traces_bg_removed, freq, fc);
        parameters = struct('fc', fc);
    case 'linear'
        traces_bleach_removed = detrend(traces_bg_removed, 1);
        baseline = traces_bg_removed - traces_bleach_removed;
        parameters = struct();
    case 'exp2'
        [~, baseline] = fit_exp2(traces_bg_removed);
        traces_bleach_removed = traces_bg_removed - baseline;
        parameters = struct();
    otherwise
        error('Unsupported bleach_mode: %s', bleach_mode);
end
end

function shift_file = find_latest_previous_motion_shift(cycle_path, current_save_path)
shift_file = "";
dual_root = fullfile(cycle_path, 'Dual_analysis3');
if isfolder(dual_root)
    listing = dir(fullfile(dual_root, '**', 'shared_motion_shifts_result.mat'));
    shift_file = select_latest_shift_candidate(listing, current_save_path);
    if strlength(shift_file) > 0
        return;
    end
end

% Fall back to legacy single-channel Cam*_Analysis outputs when no modern
% Dual_analysis3 shared shift file exists yet. In rebuilt dual-color data,
% Cam2 is the voltage channel by convention, so it is checked before other
% camera analysis folders.
legacy_voltage_listing = dir(fullfile(cycle_path, 'Cam2*_Analysis', '**', 'motion_shifts_result.mat'));
shift_file = select_latest_shift_candidate(legacy_voltage_listing, current_save_path);
if strlength(shift_file) > 0
    return;
end

legacy_listing = dir(fullfile(cycle_path, 'Cam*_Analysis', '**', 'motion_shifts_result.mat'));
shift_file = select_latest_shift_candidate(legacy_listing, current_save_path);
end

function shift_file = select_latest_shift_candidate(listing, current_save_path)
shift_file = "";
if isempty(listing)
    return;
end
current_save_path = string(current_save_path);
candidate_paths = strings(0, 1);
candidate_times = zeros(0, 1);
for idx = 1:numel(listing)
    candidate = string(fullfile(listing(idx).folder, listing(idx).name));
    if strncmpi(char(candidate), char(current_save_path), strlength(current_save_path))
        continue;
    end
    candidate_paths(end+1, 1) = candidate; %#ok<AGROW>
    candidate_times(end+1, 1) = listing(idx).datenum; %#ok<AGROW>
end

if isempty(candidate_paths)
    return;
end

[~, newest_idx] = max(candidate_times);
shift_file = candidate_paths(newest_idx);
end

function [noise_reference, noise, sensitivity, snr_value, info] = compute_voltage_metrics(traces_bleach_removed, baseline)
% Build voltage-side derived traces from the bleach-corrected signal.
% The denoised reference is used to estimate noise, not to replace the main
% signal trace.
if exist('wdenoise', 'file') ~= 2
    error('Wavelet Toolbox function wdenoise is required for voltage noise estimation.');
end

traces_bleach_removed = double(traces_bleach_removed);
baseline = double(baseline);

noise_reference = wdenoise(traces_bleach_removed, 8, DenoisingMethod='FDR', Wavelet='bior6.8');
noise_reference_method = 'wdenoise';
noise_reference_parameters = struct('level', 8, 'denoising_method', 'FDR', 'wavelet_name', 'bior6.8');

noise = traces_bleach_removed - noise_reference;
baseline_safe = baseline;
baseline_safe(abs(baseline_safe) < eps) = eps;
sensitivity = traces_bleach_removed ./ baseline_safe;
noise_std = std(noise, 0, 1);
noise_std(noise_std < eps) = eps;
snr_value = traces_bleach_removed ./ noise_std;

info = struct( ...
    'noise_reference_method', noise_reference_method, ...
    'noise_reference_parameters', noise_reference_parameters, ...
    'snr_method', 'signal_divided_by_noise_std', ...
    'snr_parameters', struct('signal_stage', 'bleach_removed', 'noise_stage', 'noise'));
end

function [noise_reference, noise, sensitivity, snr_value, info] = compute_calcium_metrics(traces_bleach_removed, baseline, freq)
% Build calcium-side derived traces from the bleach-corrected signal. The
% low-pass reference estimates the slow component, and the residual is used
% to measure noise before final calcium smoothing.
traces_bleach_removed = double(traces_bleach_removed);
baseline = double(baseline);

% Use a 2nd-order low-pass Butterworth filter as the calcium noise
% reference. MATLAB's butter expects a cutoff normalized by Nyquist:
%   Wn = fc / (freq / 2)
% Here Wn = min(0.49, 20 / freq), so for the common case freq >= 40 Hz,
% the effective cutoff is:
%   fc = (20 / freq) * (freq / 2) = 10 Hz
% The 0.49 cap only prevents invalid normalized cutoffs when freq is low.
[b, a] = butter(2, min(0.49, 20 / max(freq, 1)));
noise_reference = zeros(size(traces_bleach_removed));

for i = 1:size(traces_bleach_removed, 2)
    noise_reference(:, i) = filtfilt(b, a, traces_bleach_removed(:, i));
end

noise = traces_bleach_removed - noise_reference;
baseline_safe = baseline;
baseline_safe(abs(baseline_safe) < eps) = eps;
sensitivity = traces_bleach_removed ./ baseline_safe;
noise_std = std(noise, 0, 1);
noise_std(noise_std < eps) = eps;
snr_value = traces_bleach_removed ./ noise_std;

info = struct( ...
    'noise_reference_method', 'butterworth_filtfilt', ...
    'noise_reference_parameters', struct('order', 2, 'normalized_cutoff', min(0.49, 20 / max(freq, 1))), ...
    'snr_method', 'signal_divided_by_noise_std', ...
    'snr_parameters', struct('signal_stage', 'bleach_removed', 'noise_stage', 'noise'));
end

function stim_context = resolve_stim_context(cycle_path, cycle_manifest, record_manifest, voltage_camera_index, calcium_camera_index, stim_context_override, input_layout_info)
% Gather all stimulus-related metadata into one struct so later sections
% can reason about stimulation using one consistent source of truth.
stim_context = struct( ...
    'recordmode', "unknown", ...
    'supported', false, ...
    'stim_type', "none", ...
    'selected_program', "unknown", ...
    'selected_label', "unknown", ...
    'has_logs', false, ...
    'has_stimSpec', false, ...
    'logs_path', '', ...
    'method_manifest_path', '', ...
    'voltage_camera_index', voltage_camera_index, ...
    'calcium_camera_index', calcium_camera_index, ...
    'voltage_sync_var', "", ...
    'calcium_sync_var', "", ...
    'orientation_count', 0, ...
    'logs', struct(), ...
    'stimSpec', struct(), ...
    'stimRuntime', struct(), ...
    'method_manifest', struct());

if nargin < 7 || isempty(input_layout_info)
    input_layout_info = struct();
end

if isstruct(input_layout_info) && isfield(input_layout_info, 'mode') ...
        && strcmpi(string(input_layout_info.mode), "paired_raw_channel_folders")
    stim_context.recordmode = "raw_dual_pair";
    stim_context.supported = false;
    stim_context.stim_type = "none";
    return;
end

if nargin >= 6 && ~isempty(stim_context_override)
    stim_context = apply_stim_context_override(stim_context, stim_context_override);
    if stim_context.has_stimSpec && isfield(stim_context.stimSpec, 'orientations')
        stim_context.orientation_count = numel(stim_context.stimSpec.orientations);
    end
    if stim_context.supported && (strlength(string(stim_context.stim_type)) == 0 || strcmpi(string(stim_context.stim_type), "none"))
        stim_context.stim_type = classify_visual_stim_type(stim_context.stimSpec);
    end
    return;
end

record_path = fileparts(cycle_path);
method_path = fileparts(record_path);

method_manifest_path = '';
if ~isempty(cycle_manifest) && isfield(cycle_manifest, 'refs') && isfield(cycle_manifest.refs, 'method_manifest')
    method_manifest_path = char(string(cycle_manifest.refs.method_manifest));
elseif ~isempty(record_manifest) && isfield(record_manifest, 'refs') && isfield(record_manifest.refs, 'method_manifest')
    method_manifest_path = char(string(record_manifest.refs.method_manifest));
else
    method_manifest_candidate = fullfile(method_path, 'method_manifest.mat');
    if isfile(method_manifest_candidate)
        method_manifest_path = method_manifest_candidate;
    end
end

if isempty(method_manifest_path) || ~isfile(method_manifest_path)
    fallback_candidates = { ...
        fullfile(method_path, 'method_manifest.mat'), ...
        fullfile(record_path, 'method_manifest.mat'), ...
        fullfile(cycle_path, 'method_manifest.mat')};
    for i = 1:numel(fallback_candidates)
        if isfile(fallback_candidates{i})
            method_manifest_path = fallback_candidates{i};
            break;
        end
    end
end

if ~isempty(method_manifest_path) && isfile(method_manifest_path)
    tmp = load(method_manifest_path);
    if isfield(tmp, 'manifest')
        stim_context.method_manifest = tmp.manifest;
        stim_context.method_manifest_path = method_manifest_path;
    end
end

if ~isempty(cycle_manifest) && isfield(cycle_manifest, 'spec') && isfield(cycle_manifest.spec, 'recordmode')
    stim_context.recordmode = string(cycle_manifest.spec.recordmode);
elseif ~isempty(record_manifest) && isfield(record_manifest, 'spec') && isfield(record_manifest.spec, 'recordmode')
    stim_context.recordmode = string(record_manifest.spec.recordmode);
elseif ~isempty(fieldnames(stim_context.method_manifest)) && isfield(stim_context.method_manifest, 'spec') ...
        && isfield(stim_context.method_manifest.spec, 'recordmode')
    stim_context.recordmode = string(stim_context.method_manifest.spec.recordmode);
end

if ~isempty(fieldnames(stim_context.method_manifest)) && isfield(stim_context.method_manifest, 'spec') ...
        && isfield(stim_context.method_manifest.spec, 'stimSpec')
    stim_context.stimSpec = stim_context.method_manifest.spec.stimSpec;
    stim_context.has_stimSpec = ~isempty(fieldnames(stim_context.stimSpec));
    if isfield(stim_context.stimSpec, 'selectedProgram')
        stim_context.selected_program = string(stim_context.stimSpec.selectedProgram);
    end
    if isfield(stim_context.stimSpec, 'selectedLabel')
        stim_context.selected_label = string(stim_context.stimSpec.selectedLabel);
    end
end

if ~isempty(fieldnames(stim_context.method_manifest)) && isfield(stim_context.method_manifest, 'actual') ...
        && isfield(stim_context.method_manifest.actual, 'stimRuntime')
    stim_context.stimRuntime = stim_context.method_manifest.actual.stimRuntime;
end

logs_path = fullfile(cycle_path, 'logs.mat');
if ~isempty(cycle_manifest) && isfield(cycle_manifest, 'artifacts') && isfield(cycle_manifest.artifacts, 'logs_mat')
    candidate_logs = char(string(cycle_manifest.artifacts.logs_mat));
    if isfile(candidate_logs)
        logs_path = candidate_logs;
    end
end

if isfile(logs_path)
    tmp = load(logs_path);
    if isfield(tmp, 'logs')
        stim_context.logs = tmp.logs;
    else
        stim_context.logs = tmp;
    end
    stim_context.logs_path = logs_path;
    stim_context.has_logs = true;
end

if stim_context.has_logs && isfield(stim_context.logs, 'sync') && istable(stim_context.logs.sync)
    stim_context.voltage_sync_var = find_camera_sync_variable(stim_context.logs.sync, voltage_camera_index);
    stim_context.calcium_sync_var = find_camera_sync_variable(stim_context.logs.sync, calcium_camera_index);
end

if stim_context.has_stimSpec && isfield(stim_context.stimSpec, 'orientations')
    stim_context.orientation_count = numel(stim_context.stimSpec.orientations);
end

stim_context.supported = strcmpi(stim_context.recordmode, "visualstim") ...
    && stim_context.has_logs ...
    && stim_context.has_stimSpec ...
    && isfield(stim_context.logs, 'sync') ...
    && istable(stim_context.logs.sync) ...
    && strlength(stim_context.voltage_sync_var) > 0 ...
    && strlength(stim_context.calcium_sync_var) > 0;

if stim_context.supported
    stim_context.stim_type = classify_visual_stim_type(stim_context.stimSpec);
end
end

function stim_context = apply_stim_context_override(stim_context, stim_context_override)
% Merge a caller-provided stimulus description into the standard
% stim_context layout expected by later dual-analysis sections.
override_fields = fieldnames(stim_context_override);
for i = 1:numel(override_fields)
    field_name = override_fields{i};
    stim_context.(field_name) = stim_context_override.(field_name);
end

if isfield(stim_context_override, 'logs') && isfield(stim_context.logs, 'sync') && istable(stim_context.logs.sync)
    stim_context.has_logs = true;
end
if isfield(stim_context_override, 'stimSpec') && isstruct(stim_context.stimSpec)
    stim_context.has_stimSpec = ~isempty(fieldnames(stim_context.stimSpec));
end
if ~isfield(stim_context_override, 'supported')
    stim_context.supported = strcmpi(string(stim_context.recordmode), "visualstim") ...
        && stim_context.has_logs ...
        && stim_context.has_stimSpec ...
        && isfield(stim_context.logs, 'sync') ...
        && istable(stim_context.logs.sync) ...
        && strlength(string(stim_context.voltage_sync_var)) > 0 ...
        && strlength(string(stim_context.calcium_sync_var)) > 0;
end
end

function sync_var = find_camera_sync_variable(sync_table, camera_index)
sync_var = "";
var_names = string(sync_table.Properties.VariableNames);
camera_vars = var_names(contains(lower(var_names), "camera"));

if numel(camera_vars) >= camera_index
    sync_var = camera_vars(camera_index);
    return;
end

fallback_idx = 3 + camera_index;
if width(sync_table) >= fallback_idx
    sync_var = var_names(fallback_idx);
end
end

function stim_windows = build_visual_stim_windows(stim_context, freq_voltage, freq_calcium, nframes_voltage, nframes_calcium)
% Convert stimulus metadata into frame ranges for each channel. This is
% the bridge between acquisition logs and trace-domain stimulus analysis.
stim_windows = struct( ...
    'supported', false, ...
    'is_grating', false, ...
    'stim_type', stim_context.stim_type, ...
    'selected_program', stim_context.selected_program, ...
    'trial_count', 0, ...
    'orientations', [], ...
    'unique_orientations', [], ...
    'condition_index', [], ...
    'condition_labels', strings(0, 1), ...
    'trial_labels', strings(0, 1), ...
    'condition_colors', [], ...
    'voltage', struct(), ...
    'calcium', struct());

if ~stim_context.supported
    return;
end

switch char(stim_context.stim_type)
    case 'visualstim_grating'
        stim_windows = build_grating_stim_windows(stim_context, freq_voltage, freq_calcium, nframes_voltage, nframes_calcium);
    case 'visualstim_flash'
        stim_windows = build_flash_stim_windows(stim_context, freq_voltage, freq_calcium, nframes_voltage, nframes_calcium);
    case {'visualstim_blue', 'visualstim_luminance'}
        stim_windows = build_block_sequence_stim_windows(stim_context, freq_voltage, freq_calcium, nframes_voltage, nframes_calcium);
    case 'visualstim_flicker'
        stim_windows = build_flicker_stim_windows(stim_context, freq_voltage, freq_calcium, nframes_voltage, nframes_calcium);
    otherwise
        if isfield(stim_context.stimSpec, 'blockSequence')
            stim_windows = build_block_sequence_stim_windows(stim_context, freq_voltage, freq_calcium, nframes_voltage, nframes_calcium);
        elseif isfield(stim_context.stimSpec, 'flicker') || isfield(stim_context.stimSpec, 'contrastReverse')
            stim_windows = build_flicker_stim_windows(stim_context, freq_voltage, freq_calcium, nframes_voltage, nframes_calcium);
        else
            stim_windows = build_grating_stim_windows(stim_context, freq_voltage, freq_calcium, nframes_voltage, nframes_calcium);
        end
end
end

function stim_windows = build_flash_stim_windows(stim_context, freq_voltage, freq_calcium, nframes_voltage, nframes_calcium)
% Flash analysis uses equal-length comparison windows rather than the very
% short flash pulse itself:
%   baseline -> 2 s gray immediately before flash onset
%   response -> 2 s starting at flash onset
% The actual flash pulse is still saved separately for trace shading.
stim_windows = initialize_stim_windows(stim_context);
sync_table = stim_context.logs.sync;
ifi = resolve_stim_ifi(stim_context, sync_table);

if ~isfield(stim_context.stimSpec, 'blockSequence')
    return;
end

block_sequence = stim_context.stimSpec.blockSequence;
durations = double(block_sequence.durations(:));
labels = derive_block_labels(block_sequence);
if isempty(durations)
    return;
end
if numel(labels) ~= numel(durations)
    labels = repmat("block", numel(durations), 1);
end

repeat_count = 1;
if isfield(block_sequence, 'repeatCount') && ~isempty(block_sequence.repeatCount)
    repeat_count = max(1, round(double(block_sequence.repeatCount)));
end

[block_rows, block_labels] = build_repeated_block_rows(sync_table, durations, labels, repeat_count, ifi);
if isempty(block_rows)
    return;
end

baseline_window_sec = 2.0;
response_window_sec = 2.0;
baseline_window_frames = max(1, round(baseline_window_sec / ifi));
response_window_frames = max(1, round(response_window_sec / ifi));

flash_mask = ~arrayfun(@is_baseline_like_label, block_labels);
flash_indices = find(flash_mask);
if isempty(flash_indices)
    return;
end

baseline_rows = NaN(numel(flash_indices), 2);
stim_rows = NaN(numel(flash_indices), 2);
event_labels = strings(numel(flash_indices), 1);
shading_rows = NaN(numel(flash_indices), 2);
shading_labels = strings(numel(flash_indices), 1);

for i = 1:numel(flash_indices)
    flash_idx = flash_indices(i);
    flash_row_start = block_rows(flash_idx, 1);
    flash_row_end = block_rows(flash_idx, 2);

    base_row_end = flash_row_start - 1;
    base_row_start = max(2, base_row_end - baseline_window_frames + 1);
    stim_row_start = flash_row_start;
    stim_row_end = min(height(sync_table), stim_row_start + response_window_frames - 1);

    if base_row_end < base_row_start || stim_row_end < stim_row_start
        continue;
    end

    current_label = make_flash_condition_label(block_labels(flash_idx));
    baseline_rows(i, :) = [base_row_start, base_row_end];
    stim_rows(i, :) = [stim_row_start, stim_row_end];
    event_labels(i) = current_label;
    shading_rows(i, :) = [flash_row_start, flash_row_end];
    shading_labels(i) = current_label;
end

valid_trials = all(isfinite(baseline_rows), 2) & all(isfinite(stim_rows), 2);
baseline_rows = baseline_rows(valid_trials, :);
stim_rows = stim_rows(valid_trials, :);
event_labels = event_labels(valid_trials);
shading_rows = shading_rows(valid_trials, :);
shading_labels = shading_labels(valid_trials);
if isempty(event_labels)
    return;
end

[condition_index, condition_labels] = encode_condition_labels(event_labels);
stim_windows = finalize_stim_windows_from_rows( ...
    stim_windows, stim_context, sync_table, baseline_rows, stim_rows, ...
    condition_index, condition_labels, event_labels, ...
    [], [], freq_voltage, freq_calcium, nframes_voltage, nframes_calcium, ...
    lines(max(1, numel(condition_labels))));

[voltage_shading_frames, calcium_shading_frames, valid_shading_mask] = convert_block_rows_to_channel_frames( ...
    sync_table, stim_context.voltage_sync_var, stim_context.calcium_sync_var, shading_rows, nframes_voltage, nframes_calcium);
stim_windows.shading_labels = shading_labels(valid_shading_mask);
stim_windows.voltage.shading_labels = stim_windows.shading_labels;
stim_windows.calcium.shading_labels = stim_windows.shading_labels;
stim_windows.voltage.shading_time_ranges = voltage_shading_frames / freq_voltage;
stim_windows.calcium.shading_time_ranges = calcium_shading_frames / freq_calcium;
stim_windows.voltage.shading_frames = voltage_shading_frames;
stim_windows.calcium.shading_frames = calcium_shading_frames;
stim_windows.block_labels = strings(0, 1);
stim_windows.voltage.block_frames = [];
stim_windows.voltage.block_time_ranges = [];
stim_windows.calcium.block_frames = [];
stim_windows.calcium.block_time_ranges = [];
stim_windows.flash_windows = struct( ...
    'baseline_window_sec', baseline_window_sec, ...
    'response_window_sec', response_window_sec);
end

function stim_windows = build_grating_stim_windows(stim_context, freq_voltage, freq_calcium, nframes_voltage, nframes_calcium)
stim_windows = initialize_stim_windows(stim_context);
sync_table = stim_context.logs.sync;
orientations = double(stim_context.stimSpec.orientations(:)');
anchor_row = 2;

ifi = resolve_stim_ifi(stim_context, sync_table);

num_base_frames = max(1, round(double(stim_context.stimSpec.isi) / ifi));
num_stim_frames = max(1, round(double(stim_context.stimSpec.duration) / ifi));
frames_per_trial = num_base_frames + num_stim_frames;
available_trials = floor((height(sync_table) - 1) / frames_per_trial);
ntrials = min(numel(orientations), available_trials);
if ntrials < 1
    error('No valid stimulation trials were reconstructed from logs.sync.');
end

orientations = orientations(1:ntrials);
voltage_base_rows = NaN(ntrials, 2);
voltage_stim_rows = NaN(ntrials, 2);

for i = 1:ntrials
    base_row_start = anchor_row + (i - 1) * frames_per_trial;
    base_row_end = base_row_start + num_base_frames - 1;
    stim_row_start = base_row_end + 1;
    stim_row_end = stim_row_start + num_stim_frames - 1;

    if stim_row_end > height(sync_table)
        break;
    end

    voltage_base_rows(i, :) = [base_row_start, base_row_end];
    voltage_stim_rows(i, :) = [stim_row_start, stim_row_end];
end

valid_trials = all(isfinite(voltage_stim_rows), 2);
orientations = orientations(valid_trials);
base_rows = voltage_base_rows(valid_trials, :);
stim_rows = voltage_stim_rows(valid_trials, :);
ntrials = numel(orientations);
unique_orientations = unique(orientations, 'stable');
condition_index = zeros(ntrials, 1);
for i = 1:ntrials
    condition_index(i) = find(unique_orientations == orientations(i), 1, 'first');
end
condition_labels = string(arrayfun(@(x) sprintf('%.0f deg', x), unique_orientations, 'UniformOutput', false));
trial_labels = string(arrayfun(@(x) sprintf('%.0f deg', x), orientations(:)', 'UniformOutput', false))';
stim_windows = finalize_stim_windows_from_rows( ...
    stim_windows, stim_context, sync_table, base_rows, stim_rows, ...
    condition_index(:), condition_labels(:), trial_labels(:), ...
    orientations(:), unique_orientations(:)', ...
    freq_voltage, freq_calcium, nframes_voltage, nframes_calcium, hsv(max(1, numel(unique_orientations))));
stim_windows.is_grating = stim_context.orientation_count > 1;
end

function stim_windows = build_block_sequence_stim_windows(stim_context, freq_voltage, freq_calcium, nframes_voltage, nframes_calcium)
stim_windows = initialize_stim_windows(stim_context);
sync_table = stim_context.logs.sync;
ifi = resolve_stim_ifi(stim_context, sync_table);

if ~isfield(stim_context.stimSpec, 'blockSequence')
    return;
end

block_sequence = stim_context.stimSpec.blockSequence;
durations = double(block_sequence.durations(:));
labels = derive_block_labels(block_sequence);
if isempty(durations)
    return;
end

if numel(labels) ~= numel(durations)
    labels = repmat("block", numel(durations), 1);
end

repeat_count = 1;
if isfield(block_sequence, 'repeatCount') && ~isempty(block_sequence.repeatCount)
    repeat_count = max(1, round(double(block_sequence.repeatCount)));
end

[block_rows, block_labels] = build_repeated_block_rows(sync_table, durations, labels, repeat_count, ifi);
[baseline_rows, stim_rows, event_labels] = build_block_events(block_rows, block_labels);
if isempty(event_labels)
    return;
end

[condition_index, condition_labels] = encode_condition_labels(event_labels);
stim_windows = finalize_stim_windows_from_rows( ...
    stim_windows, stim_context, sync_table, baseline_rows, stim_rows, ...
    condition_index, condition_labels, event_labels, ...
    [], [], freq_voltage, freq_calcium, nframes_voltage, nframes_calcium, ...
    lines(max(1, numel(condition_labels))));

[voltage_block_frames, calcium_block_frames, valid_block_mask] = convert_block_rows_to_channel_frames( ...
    sync_table, stim_context.voltage_sync_var, stim_context.calcium_sync_var, block_rows, nframes_voltage, nframes_calcium);
stim_windows.block_labels = block_labels(valid_block_mask);
stim_windows.voltage.block_frames = voltage_block_frames;
stim_windows.voltage.block_time_ranges = voltage_block_frames / freq_voltage;
stim_windows.calcium.block_frames = calcium_block_frames;
stim_windows.calcium.block_time_ranges = calcium_block_frames / freq_calcium;
end

function stim_windows = build_flicker_stim_windows(stim_context, freq_voltage, freq_calcium, nframes_voltage, nframes_calcium)
stim_windows = initialize_stim_windows(stim_context);
sync_table = stim_context.logs.sync;
ifi = resolve_stim_ifi(stim_context, sync_table);

if isfield(stim_context.stimSpec, 'flicker')
    half_cycle_frames = max(1, round((1 / (2 * double(stim_context.stimSpec.flicker.frequencyHz))) / ifi));
    block_count = max(1, round(double(stim_context.stimSpec.flicker.duration) * double(stim_context.stimSpec.flicker.frequencyHz) * 2));
    block_labels = repmat(["high"; "low"], ceil(block_count / 2), 1);
    block_labels = block_labels(1:block_count);
elseif isfield(stim_context.stimSpec, 'contrastReverse')
    half_cycle_frames = max(1, round((1 / (2 * double(stim_context.stimSpec.contrastReverse.frequencyHz))) / ifi));
    block_count = max(1, round(double(stim_context.stimSpec.contrastReverse.duration) * double(stim_context.stimSpec.contrastReverse.frequencyHz) * 2));
    block_labels = repmat(["phase_pos"; "phase_neg"], ceil(block_count / 2), 1);
    block_labels = block_labels(1:block_count);
else
    return;
end

[block_rows, block_labels] = build_fixed_length_blocks(sync_table, half_cycle_frames, block_count, block_labels);
if size(block_rows, 1) < 2
    return;
end

baseline_rows = block_rows(1:end-1, :);
stim_rows = block_rows(2:end, :);
event_labels = block_labels(2:end);
event_labels = string(event_labels(:));

[condition_index, condition_labels] = encode_condition_labels(event_labels);
stim_windows = finalize_stim_windows_from_rows( ...
    stim_windows, stim_context, sync_table, baseline_rows, stim_rows, ...
    condition_index, condition_labels, event_labels, ...
    [], [], freq_voltage, freq_calcium, nframes_voltage, nframes_calcium, ...
    lines(max(1, numel(condition_labels))));
end

function stim_windows = initialize_stim_windows(stim_context)
stim_windows = struct( ...
    'supported', false, ...
    'is_grating', false, ...
    'stim_type', stim_context.stim_type, ...
    'selected_program', stim_context.selected_program, ...
    'trial_count', 0, ...
    'orientations', [], ...
    'unique_orientations', [], ...
    'condition_index', [], ...
    'condition_labels', strings(0, 1), ...
    'trial_labels', strings(0, 1), ...
    'condition_colors', [], ...
    'voltage', struct(), ...
    'calcium', struct(), ...
    'block_labels', strings(0, 1), ...
    'shading_labels', strings(0, 1));
end

function ifi = resolve_stim_ifi(stim_context, sync_table)
ifi = NaN;
if isfield(stim_context.stimRuntime, 'ifi') && ~isempty(stim_context.stimRuntime.ifi)
    ifi = double(stim_context.stimRuntime.ifi);
elseif ismember('PTB_VBL_Time', sync_table.Properties.VariableNames) && height(sync_table) > 2
    ifi = median(diff(sync_table.PTB_VBL_Time), 'omitnan');
end
if ~isfinite(ifi) || ifi <= 0
    error('Cannot determine stimulus frame interval for visual stimulation.');
end
end

function stim_windows = finalize_stim_windows_from_rows( ...
    stim_windows, stim_context, sync_table, baseline_rows, stim_rows, ...
    condition_index, condition_labels, trial_labels, orientations, unique_orientations, ...
    freq_voltage, freq_calcium, nframes_voltage, nframes_calcium, condition_colors)

[voltage_base, voltage_stim] = convert_sync_rows_to_frame_ranges( ...
    sync_table, stim_context.voltage_sync_var, baseline_rows, stim_rows, nframes_voltage);
[calcium_base, calcium_stim] = convert_sync_rows_to_frame_ranges( ...
    sync_table, stim_context.calcium_sync_var, baseline_rows, stim_rows, nframes_calcium);

valid_trials = all(isfinite(voltage_stim), 2) & all(isfinite(calcium_stim), 2);
valid_trials = valid_trials ...
    & voltage_base(:, 1) >= 1 & voltage_stim(:, 2) <= nframes_voltage ...
    & calcium_base(:, 1) >= 1 & calcium_stim(:, 2) <= nframes_calcium;

voltage_base = voltage_base(valid_trials, :);
voltage_stim = voltage_stim(valid_trials, :);
calcium_base = calcium_base(valid_trials, :);
calcium_stim = calcium_stim(valid_trials, :);
trial_labels = trial_labels(valid_trials, :);

if ~isempty(orientations)
    orientations = orientations(valid_trials, :);
end

[condition_index, condition_labels] = encode_condition_labels(trial_labels);
condition_colors = condition_colors(1:max(1, numel(condition_labels)), :);

stim_windows.supported = ~isempty(condition_index);
stim_windows.trial_count = size(voltage_stim, 1);
stim_windows.condition_index = condition_index(:);
stim_windows.condition_labels = string(condition_labels(:));
stim_windows.trial_labels = string(trial_labels(:));
stim_windows.condition_colors = condition_colors;
stim_windows.orientations = orientations;
stim_windows.unique_orientations = unique_orientations;
stim_windows.voltage = struct( ...
    'baseline_frames', voltage_base, ...
    'stim_frames', voltage_stim, ...
    'stim_time_ranges', voltage_stim / freq_voltage, ...
    'block_frames', [], ...
    'block_time_ranges', []);
stim_windows.calcium = struct( ...
    'baseline_frames', calcium_base, ...
    'stim_frames', calcium_stim, ...
    'stim_time_ranges', calcium_stim / freq_calcium, ...
    'block_frames', [], ...
    'block_time_ranges', []);
end

function [baseline_frames, stim_frames] = convert_sync_rows_to_frame_ranges(sync_table, sync_var, baseline_rows, stim_rows, nframes_channel)
% Translate stimulus timing rows into concrete frame ranges so the same
% trial logic can be applied directly to ROI traces.
series = double(sync_table.(char(sync_var)));
anchor_row = 2;
offset = series(anchor_row) - 1;

baseline_frames = NaN(size(baseline_rows));
stim_frames = NaN(size(stim_rows));
for i = 1:size(baseline_rows, 1)
    base_idx = baseline_rows(i, :);
    stim_idx = stim_rows(i, :);
    if any(base_idx < 1) || any(stim_idx < 1) || base_idx(2) > numel(series) || stim_idx(2) > numel(series)
        continue;
    end
    baseline_frames(i, :) = round([series(base_idx(1)), series(base_idx(2))] - offset);
    stim_frames(i, :) = round([series(stim_idx(1)), series(stim_idx(2))] - offset);
end

baseline_frames = max(1, baseline_frames);
stim_frames(:, 2) = min(nframes_channel, stim_frames(:, 2));
end

function [voltage_block_frames, calcium_block_frames, valid_blocks] = convert_block_rows_to_channel_frames(sync_table, voltage_sync_var, calcium_sync_var, block_rows, nframes_voltage, nframes_calcium)
% Convert full block rows into channel-specific frame ranges so block-based
% shading can follow the actual step sequence rather than only the stim
% windows used for statistics.
voltage_series = double(sync_table.(char(voltage_sync_var)));
calcium_series = double(sync_table.(char(calcium_sync_var)));
anchor_row = 2;
voltage_offset = voltage_series(anchor_row) - 1;
calcium_offset = calcium_series(anchor_row) - 1;

voltage_block_frames = NaN(size(block_rows));
calcium_block_frames = NaN(size(block_rows));
for i = 1:size(block_rows, 1)
    block_idx = block_rows(i, :);
    if any(block_idx < 1) || block_idx(2) > numel(voltage_series) || block_idx(2) > numel(calcium_series)
        continue;
    end
    voltage_block_frames(i, :) = round([voltage_series(block_idx(1)), voltage_series(block_idx(2))] - voltage_offset);
    calcium_block_frames(i, :) = round([calcium_series(block_idx(1)), calcium_series(block_idx(2))] - calcium_offset);
end

valid_blocks = all(isfinite(voltage_block_frames), 2) & all(isfinite(calcium_block_frames), 2);
valid_blocks = valid_blocks ...
    & voltage_block_frames(:, 1) >= 1 & voltage_block_frames(:, 2) <= nframes_voltage ...
    & calcium_block_frames(:, 1) >= 1 & calcium_block_frames(:, 2) <= nframes_calcium;

voltage_block_frames = voltage_block_frames(valid_blocks, :);
calcium_block_frames = calcium_block_frames(valid_blocks, :);
end

function [block_rows, block_labels] = build_repeated_block_rows(sync_table, durations, labels, repeat_count, ifi)
total_blocks = numel(durations) * repeat_count;
block_rows = NaN(total_blocks, 2);
block_labels = strings(total_blocks, 1);
row_cursor = 2;
block_counter = 0;

for rep_idx = 1:repeat_count
    for block_idx = 1:numel(durations)
        frames_this_block = max(1, round(durations(block_idx) / ifi));
        row_start = row_cursor;
        row_end = row_start + frames_this_block - 1;
        if row_end > height(sync_table)
            return;
        end
        block_counter = block_counter + 1;
        block_rows(block_counter, :) = [row_start, row_end];
        block_labels(block_counter) = labels(block_idx);
        row_cursor = row_end + 1;
    end
end

block_rows = block_rows(1:block_counter, :);
block_labels = block_labels(1:block_counter);
end

function [block_rows, block_labels] = build_fixed_length_blocks(sync_table, frames_per_block, block_count, labels)
block_rows = NaN(block_count, 2);
block_labels = strings(block_count, 1);
row_cursor = 2;
valid_count = 0;

for block_idx = 1:block_count
    row_start = row_cursor;
    row_end = row_start + frames_per_block - 1;
    if row_end > height(sync_table)
        break;
    end
    valid_count = valid_count + 1;
    block_rows(valid_count, :) = [row_start, row_end];
    block_labels(valid_count) = labels(block_idx);
    row_cursor = row_end + 1;
end

block_rows = block_rows(1:valid_count, :);
block_labels = block_labels(1:valid_count);
end

function [baseline_rows, stim_rows, event_labels] = build_block_events(block_rows, block_labels)
baseline_rows = NaN(0, 2);
stim_rows = NaN(0, 2);
event_labels = strings(0, 1);

for block_idx = 2:numel(block_labels)
    current_label = block_labels(block_idx);
    if is_baseline_like_label(current_label)
        continue;
    end

    baseline_idx = find_previous_baseline_block(block_labels, block_idx);
    if isempty(baseline_idx)
        baseline_idx = block_idx - 1;
    end

    baseline_rows(end+1, :) = block_rows(baseline_idx, :); %#ok<AGROW>
    stim_rows(end+1, :) = block_rows(block_idx, :); %#ok<AGROW>
    event_labels(end+1, 1) = current_label; %#ok<AGROW>
end
end

function baseline_idx = find_previous_baseline_block(block_labels, current_idx)
baseline_idx = [];
for idx = current_idx-1:-1:1
    if is_baseline_like_label(block_labels(idx))
        baseline_idx = idx;
        return;
    end
end
end

function labels = derive_block_labels(block_sequence)
labels = strings(0, 1);
if isfield(block_sequence, 'labels') && ~isempty(block_sequence.labels)
    raw_labels = block_sequence.labels;
    labels = string(raw_labels(:));
end

if ~isempty(labels) && all(strlength(labels) > 0)
    return;
end

if isfield(block_sequence, 'colors') && ~isempty(block_sequence.colors)
    colors = double(block_sequence.colors);
    if max(colors(:)) > 1
        colors = colors / 255;
    end
    labels = strings(size(colors, 1), 1);
    for i = 1:size(colors, 1)
        labels(i) = classify_color_block(colors(i, :));
    end
end
end

function label = make_flash_condition_label(raw_label)
raw_label = lower(char(string(raw_label)));
if contains(raw_label, 'white')
    label = "white_flash";
elseif contains(raw_label, 'black')
    label = "black_flash";
elseif contains(raw_label, 'blue')
    label = "blue_flash";
else
    label = "flash";
end
end

function label = classify_color_block(rgb)
rgb = double(rgb(:)');
if numel(rgb) ~= 3
    label = "block";
    return;
end

if max(rgb) > 1
    rgb = rgb / 255;
end

if max(abs(rgb - mean(rgb))) < 0.08
    if mean(rgb) < 0.2
        label = "black";
    elseif mean(rgb) > 0.8
        label = "white";
    else
        label = "gray";
    end
elseif rgb(3) > rgb(1) + 0.15 && rgb(3) > rgb(2) + 0.15
    label = "blue";
else
    label = "color";
end
end

function tf = is_baseline_like_label(label)
label = lower(char(string(label)));
tf = contains(label, 'gray') || contains(label, 'baseline') || contains(label, 'rest') || contains(label, 'isi');
end

function [condition_index, condition_labels] = encode_condition_labels(event_labels)
event_labels = string(event_labels(:));
condition_labels = unique(event_labels, 'stable');
condition_index = zeros(numel(event_labels), 1);
for i = 1:numel(event_labels)
    condition_index(i) = find(condition_labels == event_labels(i), 1, 'first');
end
end

function stim_type = classify_visual_stim_type(stimSpec)
selected_program = lower(char(get_struct_string(stimSpec, 'selectedProgram', "unknown")));
selected_label = lower(char(get_struct_string(stimSpec, 'selectedLabel', "unknown")));

if strcmp(selected_program, 'flash') || contains(selected_program, 'flash') || contains(selected_label, 'flash')
    stim_type = "visualstim_flash";
elseif strcmp(selected_program, 'gray_blue_gray') || contains(selected_program, 'blue') || contains(selected_label, 'blue')
    stim_type = "visualstim_blue";
elseif strcmp(selected_program, 'gray_white_gray_black') ...
        || contains(selected_program, 'white') || contains(selected_program, 'black') ...
        || contains(selected_label, 'white') || contains(selected_label, 'black')
    stim_type = "visualstim_luminance";
elseif strcmp(selected_program, 'white_black_flicker') || strcmp(selected_program, 'contrast_reverse') ...
        || contains(selected_program, 'flicker') || contains(selected_program, 'contrast')
    stim_type = "visualstim_flicker";
elseif strcmp(selected_program, 'drifting_grating') || contains(selected_program, 'grating')
    stim_type = "visualstim_grating";
elseif (isfield(stimSpec, 'orientations') && numel(stimSpec.orientations) > 1)
    stim_type = "visualstim_grating";
else
    stim_type = "visualstim_generic";
end
end

function value = get_struct_string(s, field_name, default_value)
value = default_value;
if isstruct(s) && isfield(s, field_name) && ~isempty(s.(field_name))
    value = string(s.(field_name));
end
end

function print_stim_context_summary(stim_context, stim_windows, caller_name)
sync_available = false;
if isfield(stim_context, 'has_logs') && stim_context.has_logs ...
        && isfield(stim_context, 'logs') && isstruct(stim_context.logs) ...
        && isfield(stim_context.logs, 'sync') && istable(stim_context.logs.sync)
    sync_available = true;
end
fprintf('Stim context summary [%s]:\n', caller_name);
fprintf('  record mode: %s\n', string(stim_context.recordmode));
fprintf('  selected program: %s\n', string(stim_context.selected_program));
fprintf('  selected label: %s\n', string(stim_context.selected_label));
fprintf('  classified stim type: %s\n', string(stim_context.stim_type));
fprintf('  logs available: %d | stimSpec available: %d | sync available: %d\n', ...
    stim_context.has_logs, stim_context.has_stimSpec, sync_available);
fprintf('  voltage sync var: %s | calcium sync var: %s\n', ...
    string(stim_context.voltage_sync_var), string(stim_context.calcium_sync_var));
fprintf('  OSI/DSI eligible: %s\n', ternary(strcmpi(string(stim_context.stim_type), "visualstim_grating"), 'yes', 'no'));

if nargin < 2 || isempty(stim_windows) || ~isstruct(stim_windows) || ~isfield(stim_windows, 'supported')
    fprintf('  stim windows: not built in this section\n');
    return;
end

fprintf('  stim windows supported: %d\n', stim_windows.supported);
if ~stim_windows.supported
    return;
end

fprintf('  trial count: %d\n', stim_windows.trial_count);
fprintf('  condition count: %d\n', numel(stim_windows.condition_labels));
fprintf('  condition labels: %s\n', preview_string_list(stim_windows.condition_labels, 8));
fprintf('  trial labels preview: %s\n', preview_string_list(stim_windows.trial_labels, 10));
end

function out = preview_string_list(values, max_items)
values = string(values(:));
if isempty(values)
    out = '(none)';
    return;
end

if numel(values) <= max_items
    out = strjoin(cellstr(values), ', ');
else
    head_values = values(1:max_items);
    out = sprintf('%s, ... (%d total)', strjoin(cellstr(head_values), ', '), numel(values));
end
end

function plot_optional_stage(stage_data, t_axis, is_available, missing_text)
if is_available && ~isempty(stage_data)
    plot_offset_stage(stage_data, t_axis);
else
    axis tight;
    text(0.5, 0.5, missing_text, 'Units', 'normalized', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    set(gca, 'XTick', [], 'YTick', []);
    box on;
end
end

function plot_offset_stage(stage_data, t_axis)
offset_array = offset_plot(stage_data, t_axis);
annotate_offset_plot(gca, stage_data, t_axis, offset_array);
end

function annotate_offset_plot(ax, trace_data, t_axis, offset_array)
trace_data = double(trace_data);
t_axis = double(t_axis(:));
offset_array = double(offset_array(:)');

if isempty(trace_data) || isempty(t_axis)
    return;
end

nrois = size(trace_data, 2);
x_min = min(t_axis);
x_max = max(t_axis);
x_span = max(eps, x_max - x_min);
left_margin = 0.14 * x_span;
right_margin = 0.10 * x_span;
xlim(ax, [x_min - left_margin, x_max + right_margin]);

for roi_idx = 1:nrois
    y_anchor = offset_array(roi_idx) + mean(trace_data(:, roi_idx), 'omitnan');
    if ~isfinite(y_anchor)
        y_anchor = offset_array(roi_idx);
    end
    text(ax, x_min - left_margin * 0.92, y_anchor, sprintf('ROI %d', roi_idx), ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
        'FontSize', 8, 'FontWeight', 'bold', 'Color', [0.15 0.15 0.15], ...
        'Interpreter', 'none');
end

add_axis_scalebar(ax, t_axis, trace_data(:), [0.1 0.1 0.1], '');
end

function add_axis_scalebar(ax, x_values, y_values, bar_color, y_suffix)
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

x_len = nice_scalebar_value(0.16 * x_span);
y_len = nice_scalebar_value(0.18 * range(y_values));
if ~isfinite(y_len) || y_len <= 0
    y_len = nice_scalebar_value(0.18 * y_span);
end

x_start = x_limits(2) - 0.06 * x_span - x_len;
y_start = y_limits(1) + 0.10 * y_span;

plot(ax, [x_start, x_start + x_len], [y_start, y_start], ...
    'Color', bar_color, 'LineWidth', 1.4, 'Clipping', 'off');
plot(ax, [x_start, x_start], [y_start, y_start + y_len], ...
    'Color', bar_color, 'LineWidth', 1.4, 'Clipping', 'off');

text(ax, x_start + x_len / 2, y_start - 0.04 * y_span, ...
    sprintf('%s s', format_scalebar_value(x_len)), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
    'FontSize', 8, 'Color', bar_color, 'Interpreter', 'none');

if strlength(string(y_suffix)) > 0
    y_label = sprintf('%s %s', format_scalebar_value(y_len), y_suffix);
else
    y_label = format_scalebar_value(y_len);
end
text(ax, x_start - 0.01 * x_span, y_start + y_len / 2, y_label, ...
    'HorizontalAlignment', 'right', 'VerticalAlignment', 'middle', ...
    'Rotation', 90, 'FontSize', 8, 'Color', bar_color, 'Interpreter', 'none');
end

function text_value = format_p_value(p_value)
p_value = double(p_value);
if ~isfinite(p_value)
    text_value = 'NaN';
elseif p_value < 1e-3
    text_value = sprintf('%.2e', p_value);
else
    text_value = sprintf('%.3f', p_value);
end
end

function value = nice_scalebar_value(target_value)
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

function text_value = format_scalebar_value(value)
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

function limits = compute_shared_ylim(trace_matrix)
trace_matrix = double(trace_matrix);
trace_values = trace_matrix(isfinite(trace_matrix));
if isempty(trace_values)
    limits = [-1, 1];
    return;
end

trace_min = min(trace_values);
trace_max = max(trace_values);
if abs(trace_max - trace_min) < eps
    pad = max(1, abs(trace_max) * 0.2);
else
    pad = 0.08 * (trace_max - trace_min);
end
limits = [trace_min - pad, trace_max + pad];
end

function add_trace_badge(ax, badge_text, badge_color, x_pos, y_pos)
if nargin < 4
    x_pos = 0.02;
end
if nargin < 5
    y_pos = 0.84;
end

text(ax, x_pos, y_pos, badge_text, 'Units', 'normalized', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', ...
    'FontSize', 8, 'FontWeight', 'bold', 'Color', badge_color, ...
    'Interpreter', 'none');
end

function add_stim_shading(ax, channel_windows, condition_index, condition_colors, alpha_value, trial_labels, block_labels)
if isempty(channel_windows) || ~isstruct(channel_windows)
    return;
end
if nargin < 6
    trial_labels = strings(size(condition_index));
else
    trial_labels = string(trial_labels(:));
end
if nargin < 7
    block_labels = strings(0, 1);
else
    block_labels = string(block_labels(:));
end

if isfield(channel_windows, 'shading_time_ranges') && ~isempty(channel_windows.shading_time_ranges) ...
        && isfield(channel_windows, 'shading_labels') && numel(channel_windows.shading_labels) == size(channel_windows.shading_time_ranges, 1)
    time_ranges = channel_windows.shading_time_ranges;
    shading_labels = string(channel_windows.shading_labels(:));
    use_default_colors = false;
elseif isfield(channel_windows, 'block_time_ranges') && ~isempty(channel_windows.block_time_ranges) && numel(block_labels) == size(channel_windows.block_time_ranges, 1)
    time_ranges = channel_windows.block_time_ranges;
    shading_labels = block_labels;
    use_default_colors = false;
else
    if ~isfield(channel_windows, 'stim_time_ranges') || isempty(channel_windows.stim_time_ranges)
        return;
    end
    time_ranges = channel_windows.stim_time_ranges;
    shading_labels = trial_labels;
    use_default_colors = true;
end

yl = ylim(ax);
hold(ax, 'on');
for i = 1:size(time_ranges, 1)
    if any(~isfinite(time_ranges(i, :)))
        continue;
    end
    if use_default_colors
        color = condition_colors(condition_index(i), :);
    else
        color = [0.7 0.7 0.7];
    end
    alpha_i = alpha_value;
    if numel(shading_labels) >= i
        [color, alpha_i] = resolve_stim_shading_style(shading_labels(i), color, alpha_value);
    end
    p = patch(ax, ...
        [time_ranges(i, 1), time_ranges(i, 2), time_ranges(i, 2), time_ranges(i, 1)], ...
        [yl(1), yl(1), yl(2), yl(2)], ...
        color, 'FaceAlpha', alpha_i, 'EdgeColor', 'none');
    uistack(p, 'bottom');
end
ylim(ax, yl);
end

function [shade_color, shade_alpha] = resolve_stim_shading_style(trial_label, default_color, default_alpha)
label = lower(char(string(trial_label)));
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

function metrics = compute_stim_trial_metrics(trace_matrix, channel_windows, polarity, frame_rate)
% Reduce each trial into baseline-vs-stim summary numbers. These metrics
% are the main numeric outputs used for condition comparison across ROIs.
trace_matrix = polarity * trace_matrix;
nrois = size(trace_matrix, 2);
ntrials = size(channel_windows.stim_frames, 1);

metrics = struct();
metrics.frame_rate = frame_rate;
metrics.polarity = polarity;
metrics.baseline_mean = NaN(nrois, ntrials);
metrics.stim_mean = NaN(nrois, ntrials);
metrics.delta_mean = NaN(nrois, ntrials);
metrics.baseline_peak = NaN(nrois, ntrials);
metrics.stim_peak = NaN(nrois, ntrials);
metrics.delta_peak = NaN(nrois, ntrials);
metrics.baseline_auc = NaN(nrois, ntrials);
metrics.stim_auc = NaN(nrois, ntrials);
metrics.delta_auc = NaN(nrois, ntrials);
metrics.summary = repmat(struct(), nrois, 1);

for trial_idx = 1:ntrials
    base_range = channel_windows.baseline_frames(trial_idx, 1):channel_windows.baseline_frames(trial_idx, 2);
    stim_range = channel_windows.stim_frames(trial_idx, 1):channel_windows.stim_frames(trial_idx, 2);

    base_seg = trace_matrix(base_range, :);
    stim_seg = trace_matrix(stim_range, :);

    metrics.baseline_mean(:, trial_idx) = mean(base_seg, 1, 'omitnan')';
    metrics.stim_mean(:, trial_idx) = mean(stim_seg, 1, 'omitnan')';
    metrics.delta_mean(:, trial_idx) = metrics.stim_mean(:, trial_idx) - metrics.baseline_mean(:, trial_idx);

    metrics.baseline_peak(:, trial_idx) = max(base_seg, [], 1)';
    metrics.stim_peak(:, trial_idx) = max(stim_seg, [], 1)';
    metrics.delta_peak(:, trial_idx) = metrics.stim_peak(:, trial_idx) - metrics.baseline_peak(:, trial_idx);

    metrics.baseline_auc(:, trial_idx) = trapz(base_seg, 1)' / frame_rate;
    metrics.stim_auc(:, trial_idx) = trapz(stim_seg, 1)' / frame_rate;
    metrics.delta_auc(:, trial_idx) = metrics.stim_auc(:, trial_idx) - metrics.baseline_auc(:, trial_idx);
end

for roi_idx = 1:nrois
    metrics.summary(roi_idx).delta_mean = mean(metrics.delta_mean(roi_idx, :), 'omitnan');
    metrics.summary(roi_idx).delta_peak = mean(metrics.delta_peak(roi_idx, :), 'omitnan');
    metrics.summary(roi_idx).delta_auc = mean(metrics.delta_auc(roi_idx, :), 'omitnan');
    metrics.summary(roi_idx).p_mean = paired_test_p(metrics.baseline_mean(roi_idx, :), metrics.stim_mean(roi_idx, :));
    metrics.summary(roi_idx).p_peak = paired_test_p(metrics.baseline_peak(roi_idx, :), metrics.stim_peak(roi_idx, :));
    metrics.summary(roi_idx).p_auc = paired_test_p(metrics.baseline_auc(roi_idx, :), metrics.stim_auc(roi_idx, :));
end
end

function p_value = paired_test_p(x, y)
x = x(:);
y = y(:);
valid = isfinite(x) & isfinite(y);
x = x(valid);
y = y(valid);

if numel(x) < 2
    p_value = NaN;
    return;
end

try
    p_value = signrank(x, y);
catch
    [~, p_value] = ttest(x, y);
end
end

function tuning = compute_grating_tuning(delta_mean, orientations)
orientations = orientations(:);
unique_orientations = unique(orientations, 'stable');
nrois = size(delta_mean, 1);
nconditions = numel(unique_orientations);
response_by_condition = NaN(nrois, nconditions);

for i = 1:nconditions
    mask = orientations == unique_orientations(i);
    response_by_condition(:, i) = mean(max(delta_mean(:, mask), 0), 2, 'omitnan');
end

tuning = struct();
tuning.unique_orientations = unique_orientations(:)';
tuning.response_by_condition = response_by_condition;
tuning.pref_dir = NaN(nrois, 1);
tuning.pref_ori = NaN(nrois, 1);
tuning.gdsi = NaN(nrois, 1);
tuning.gosi = NaN(nrois, 1);
tuning.dsi = NaN(nrois, 1);
tuning.osi = NaN(nrois, 1);

theta_rad = deg2rad(unique_orientations(:)');
for roi_idx = 1:nrois
    r = response_by_condition(roi_idx, :);
    if all(~isfinite(r)) || sum(r, 'omitnan') <= 0
        continue;
    end

    [~, pref_idx] = max(r);
    tuning.pref_dir(roi_idx) = unique_orientations(pref_idx);
    tuning.pref_ori(roi_idx) = mod(tuning.pref_dir(roi_idx), 180);
    tuning.gdsi(roi_idx) = abs(nansum(r .* exp(1i * theta_rad)) / nansum(r));
    tuning.gosi(roi_idx) = abs(nansum(r .* exp(1i * 2 * theta_rad)) / nansum(r));

    opp_angle = mod(tuning.pref_dir(roi_idx) + 180, 360);
    orth_angle = mod(tuning.pref_dir(roi_idx) + 90, 360);
    opp_idx = find_closest_angle_index(unique_orientations, opp_angle);
    orth_idx = find_closest_angle_index(unique_orientations, orth_angle);
    pref_resp = r(pref_idx);
    opp_resp = r(opp_idx);
    orth_resp = r(orth_idx);
    tuning.dsi(roi_idx) = (pref_resp - opp_resp) / max(eps, pref_resp + opp_resp);
    tuning.osi(roi_idx) = (pref_resp - orth_resp) / max(eps, pref_resp + orth_resp);
end
end

function [voltage_tuning, calcium_tuning, per_roi_tuning] = compute_dual_grating_tuning_by_roi(voltage_delta_mean, calcium_delta_mean, orientations)
% Compute grating tuning once per ROI, then assemble population arrays.
nrois = max(size(voltage_delta_mean, 1), size(calcium_delta_mean, 1));
voltage_tuning = compute_grating_tuning(voltage_delta_mean, orientations);
calcium_tuning = compute_grating_tuning(calcium_delta_mean, orientations);
per_roi_tuning = repmat(empty_grating_roi_tuning(), nrois, 1);

for roi_idx = 1:nrois
    voltage_roi_delta = NaN(1, numel(orientations));
    calcium_roi_delta = NaN(1, numel(orientations));
    if roi_idx <= size(voltage_delta_mean, 1)
        voltage_roi_delta = voltage_delta_mean(roi_idx, :);
    end
    if roi_idx <= size(calcium_delta_mean, 1)
        calcium_roi_delta = calcium_delta_mean(roi_idx, :);
    end

    per_roi_tuning(roi_idx).roi_index = roi_idx;
    per_roi_tuning(roi_idx).voltage = compute_grating_tuning(voltage_roi_delta, orientations);
    per_roi_tuning(roi_idx).calcium = compute_grating_tuning(calcium_roi_delta, orientations);
end
end

function entry = empty_grating_roi_tuning()
entry = struct( ...
    'roi_index', NaN, ...
    'voltage', struct(), ...
    'calcium', struct(), ...
    'mat_file', "", ...
    'fig_file', "", ...
    'png_file', "");
end

function per_roi_tuning = save_grating_tuning_by_roi(per_roi_tuning, save_path, output_dir_name)
% Save each ROI grating tuning result as its own MAT/FIG/PNG bundle.
roi_tuning_dir = fullfile(save_path, output_dir_name);
if ~isfolder(roi_tuning_dir)
    mkdir(roi_tuning_dir);
end

for roi_idx = 1:numel(per_roi_tuning)
    roi_entry = per_roi_tuning(roi_idx);
    voltage_roi_tuning = roi_entry.voltage;
    calcium_roi_tuning = roi_entry.calcium;
    roi_label = sprintf('roi_%03d', roi_entry.roi_index);
    mat_file = fullfile(roi_tuning_dir, sprintf('%s_stim_tuning.mat', roi_label));
    fig_file = fullfile(roi_tuning_dir, sprintf('%s_stim_tuning.fig', roi_label));
    png_file = fullfile(roi_tuning_dir, sprintf('%s_stim_tuning.png', roi_label));

    roi_entry.mat_file = string(mat_file);
    roi_entry.fig_file = string(fig_file);
    roi_entry.png_file = string(png_file);
    save(mat_file, 'roi_entry', 'voltage_roi_tuning', 'calcium_roi_tuning');

    fig = figure('Color', 'w', 'Name', sprintf('ROI %03d Grating Tuning', roi_entry.roi_index), ...
        'Position', [120, 120, 1000, 750]);
    subplot(2, 2, 1);
    plot_tuning_single(voltage_roi_tuning, 'r', sprintf('ROI %03d Voltage', roi_entry.roi_index));
    subplot(2, 2, 2);
    plot_tuning_single(calcium_roi_tuning, 'g', sprintf('ROI %03d Calcium', roi_entry.roi_index));
    subplot(2, 2, 3);
    plot_roi_tuning_metric_pair(voltage_roi_tuning, calcium_roi_tuning, 'gosi', 'gOSI');
    subplot(2, 2, 4);
    plot_roi_tuning_metric_pair(voltage_roi_tuning, calcium_roi_tuning, 'pref_dir', 'Preferred direction');
    save_figure_bundle(fig, fig_file, png_file);
    close(fig);

    per_roi_tuning(roi_idx) = roi_entry;
end
end

function plot_tuning_single(tuning, line_color, panel_title)
angles = tuning.unique_orientations(:)';
response = tuning.response_by_condition;
if isempty(response)
    response = NaN(1, numel(angles));
elseif size(response, 1) > 1
    response = response(1, :);
end
plot(angles, response, '-o', 'Color', line_color, 'MarkerFaceColor', line_color, 'LineWidth', 1.5);
xlabel('Direction (deg)');
ylabel('Net Response');
title(panel_title);
grid on;
end

function plot_roi_tuning_metric_pair(voltage_tuning, calcium_tuning, field_name, metric_label)
voltage_value = extract_first_finite_field(voltage_tuning, field_name);
calcium_value = extract_first_finite_field(calcium_tuning, field_name);
bar([voltage_value, calcium_value], 'FaceColor', [0.45, 0.55, 0.70]);
set(gca, 'XTickLabel', {'Voltage', 'Calcium'});
ylabel(metric_label);
title(metric_label);
grid on;
end

function value = extract_first_finite_field(s, field_name)
value = NaN;
if isstruct(s) && isfield(s, field_name)
    candidate = s.(field_name);
    candidate = candidate(isfinite(candidate));
    if ~isempty(candidate)
        value = candidate(1);
    end
end
end

function idx = find_closest_angle_index(angles, target_angle)
wrapped_diff = abs(mod(angles - target_angle + 180, 360) - 180);
[~, idx] = min(wrapped_diff);
end

function plot_tuning_population(angles, response_by_condition, line_color, panel_title)
mean_response = mean(response_by_condition, 1, 'omitnan');
sem_response = std(response_by_condition, 0, 1, 'omitnan') ./ sqrt(max(1, sum(isfinite(response_by_condition), 1)));
errorbar(angles, mean_response, sem_response, 'Color', line_color, 'LineWidth', 1.5);
xlabel('Direction (deg)');
ylabel('Mean Net Response');
title([panel_title ' Tuning']);
grid on;
end

function plot_dual_raw_traces(traces_voltage_raw, traces_calcium_raw, t_voltage, t_calcium)
figure('Color', 'w');
subplot(2, 1, 1);
hold on;
title('Voltage Raw');
plot_offset_stage(traces_voltage_raw, t_voltage);
xlabel('Time (s)');

subplot(2, 1, 2);
hold on;
title('Calcium Raw');
plot_offset_stage(traces_calcium_raw, t_calcium);
xlabel('Time (s)');
end

function [fig_filename, png_filename] = plot_dual_background_summary( ...
    rois, mean_voltage_image, mean_calcium_image, ...
    background_mask_voltage, background_mask_calcium, ...
    background_fit_voltage, background_fit_calcium, ...
    traces_voltage_bg, traces_calcium_bg, ...
    t_voltage, t_calcium, save_path)

nrois = max(rois.bwmask(:));
colors = lines(max(1, nrois));
fig = figure('Name', 'Dual Background Correction Summary', 'Units', 'normalized', ...
    'Position', [0.05 0.08 0.9 0.8], 'Color', 'w');

subplot(2, 2, 1);
plot_background_mask_panel(mean_voltage_image, rois.bwmask, background_mask_voltage, colors, 'Voltage: ROI And Background');

subplot(2, 2, 2);
plot_background_trace_panel(traces_voltage_bg, background_fit_voltage, t_voltage, colors, 'Voltage: Raw/Fitted BG/Corrected');

subplot(2, 2, 3);
plot_background_mask_panel(mean_calcium_image, rois.bwmask_ca, background_mask_calcium, colors, 'Calcium: ROI And Background');

subplot(2, 2, 4);
plot_background_trace_panel(traces_calcium_bg, background_fit_calcium, t_calcium, colors, 'Calcium: Raw/Fitted BG/Corrected');

fig_filename = fullfile(save_path, '1_dual_background_correction_summary.fig');
png_filename = fullfile(save_path, '1_dual_background_correction_summary.png');
save_figure_bundle(fig, fig_filename, png_filename);
close(fig);
end

function plot_background_mask_panel(mean_image, roi_mask, background_mask, colors, panel_title)
if isempty(mean_image)
    axis off;
    title(panel_title);
    text(0.5, 0.5, 'Mean image unavailable in this section run', 'HorizontalAlignment', 'center');
    return;
end

imshow(mean_image, [], 'InitialMagnification', 'fit');
hold on;
nrois = max(roi_mask(:));
for i = 1:nrois
    current_color = colors(i, :);
    roi_bw = (roi_mask == i);
    roi_boundaries = bwboundaries(roi_bw);
    for k = 1:length(roi_boundaries)
        boundary = roi_boundaries{k};
        plot(boundary(:, 2), boundary(:, 1), 'Color', current_color, 'LineWidth', 1);
    end

    bg_bw = (background_mask == i);
    bg_boundaries = bwboundaries(bg_bw);
    for k = 1:length(bg_boundaries)
        boundary = bg_boundaries{k};
        plot(boundary(:, 2), boundary(:, 1), ':', 'Color', current_color, 'LineWidth', 1);
    end

    stats = regionprops(roi_bw, 'Centroid');
    if ~isempty(stats)
        text(stats(1).Centroid(1), stats(1).Centroid(2), num2str(i), ...
            'Color', current_color, 'FontWeight', 'bold', 'FontSize', 10, 'HorizontalAlignment', 'center');
    end
end
title(panel_title);
end

function plot_background_trace_panel(corrected_traces, background_fit, t_axis, colors, panel_title)
hold on;
nrois = size(corrected_traces, 2);
p1 = plot(nan, nan, 'Color', [0.7 0.7 0.7]);
p2 = plot(nan, nan, 'k--', 'LineWidth', 1);
p3 = plot(nan, nan, 'k', 'LineWidth', 1);

for i = 1:nrois
    current_color = colors(i, :);
    raw_signal = corrected_traces(:, i) + background_fit(:, i);
    faded_color = current_color * 0.4 + [1 1 1] * 0.6;
    plot(t_axis, raw_signal, 'Color', faded_color, 'LineWidth', 0.5);
    plot(t_axis, background_fit(:, i), '--', 'Color', current_color, 'LineWidth', 1);
    plot(t_axis, corrected_traces(:, i), 'Color', current_color, 'LineWidth', 1);
end

xlabel('Time (s)');
ylabel('Intensity');
title(panel_title);
grid on;
legend([p1, p2, p3], {'Raw Signal', 'Fitted Background', 'Corrected Signal'}, 'Location', 'northeast');
end

function summary = summarize_stim_by_condition(metrics, stim_windows)
nconditions = numel(stim_windows.condition_labels);
nrois = size(metrics.delta_mean, 1);

summary = struct( ...
    'condition_labels', string(stim_windows.condition_labels(:)), ...
    'response_by_condition', NaN(nrois, nconditions), ...
    'p_by_condition', NaN(nrois, nconditions), ...
    'population_mean', NaN(1, nconditions), ...
    'population_sem', NaN(1, nconditions));

for cond_idx = 1:nconditions
    mask = stim_windows.condition_index == cond_idx;
    summary.response_by_condition(:, cond_idx) = mean(metrics.delta_mean(:, mask), 2, 'omitnan');
    for roi_idx = 1:nrois
        summary.p_by_condition(roi_idx, cond_idx) = paired_test_p( ...
            metrics.baseline_mean(roi_idx, mask), metrics.stim_mean(roi_idx, mask));
    end
end

summary.population_mean = mean(summary.response_by_condition, 1, 'omitnan');
summary.population_sem = std(summary.response_by_condition, 0, 1, 'omitnan') ...
    ./ sqrt(max(1, sum(isfinite(summary.response_by_condition), 1)));
end

function summary = summarize_block_by_condition(trace_matrix, channel_windows, block_labels, polarity, trial_labels)
trace_matrix = polarity * double(trace_matrix);

if isfield(channel_windows, 'block_frames') && ~isempty(channel_windows.block_frames) ...
        && numel(block_labels) == size(channel_windows.block_frames, 1)
    labels = string(block_labels(:));
    [condition_index, condition_labels] = encode_condition_labels(labels);
    block_frames = channel_windows.block_frames;
    nrois = size(trace_matrix, 2);
    nblocks = size(block_frames, 1);
    block_mean = NaN(nrois, nblocks);

    for block_idx = 1:nblocks
        frame_range = block_frames(block_idx, 1):block_frames(block_idx, 2);
        block_mean(:, block_idx) = mean(trace_matrix(frame_range, :), 1, 'omitnan')';
    end

    response_by_condition = NaN(nrois, numel(condition_labels));
    for cond_idx = 1:numel(condition_labels)
        mask = condition_index == cond_idx;
        response_by_condition(:, cond_idx) = mean(block_mean(:, mask), 2, 'omitnan');
    end
else
    labels = string(trial_labels(:));
    [condition_index, condition_labels] = encode_condition_labels(labels);
    nrois = size(trace_matrix, 2);
    ntrials = size(channel_windows.stim_frames, 1);
    stim_mean = NaN(nrois, ntrials);

    for trial_idx = 1:ntrials
        frame_range = channel_windows.stim_frames(trial_idx, 1):channel_windows.stim_frames(trial_idx, 2);
        stim_mean(:, trial_idx) = mean(trace_matrix(frame_range, :), 1, 'omitnan')';
    end

    response_by_condition = NaN(nrois, numel(condition_labels));
    for cond_idx = 1:numel(condition_labels)
        mask = condition_index == cond_idx;
        response_by_condition(:, cond_idx) = mean(stim_mean(:, mask), 2, 'omitnan');
    end
end

summary = struct( ...
    'condition_labels', string(condition_labels(:)), ...
    'response_by_condition', response_by_condition, ...
    'population_mean', mean(response_by_condition, 1, 'omitnan'), ...
    'population_sem', std(response_by_condition, 0, 1, 'omitnan') ./ sqrt(max(1, sum(isfinite(response_by_condition), 1))));
end

function plot_stim_condition_summary(voltage_summary, calcium_summary, stim_context, save_path, file_stem, figure_label, y_label_text)
nrois = max(size(voltage_summary.response_by_condition, 1), size(calcium_summary.response_by_condition, 1));
fig = figure('Color', 'w', 'Position', [120, 80, 1300, max(420, 190 * nrois)]);

voltage_labels = cellstr(voltage_summary.condition_labels(:));
calcium_labels = cellstr(calcium_summary.condition_labels(:));
xv = 1:numel(voltage_labels);
xc = 1:numel(calcium_labels);

for roi_idx = 1:nrois
    ax_v = subplot(nrois, 2, 2 * roi_idx - 1);
    if roi_idx <= size(voltage_summary.response_by_condition, 1)
        bar(ax_v, xv, voltage_summary.response_by_condition(roi_idx, :), 'FaceColor', [0.85 0.3 0.3]);
    end
    title(ax_v, sprintf('Voltage ROI %d', roi_idx));
    ylabel(ax_v, y_label_text);
    xticks(ax_v, xv); xticklabels(ax_v, voltage_labels); xtickangle(ax_v, 25);
    grid(ax_v, 'on');
    if roi_idx < nrois
        ax_v.XTickLabel = [];
    else
        xlabel(ax_v, 'Condition');
    end

    ax_c = subplot(nrois, 2, 2 * roi_idx);
    if roi_idx <= size(calcium_summary.response_by_condition, 1)
        bar(ax_c, xc, calcium_summary.response_by_condition(roi_idx, :), 'FaceColor', [0.2 0.7 0.3]);
    end
    title(ax_c, sprintf('Calcium ROI %d', roi_idx));
    ylabel(ax_c, y_label_text);
    xticks(ax_c, xc); xticklabels(ax_c, calcium_labels); xtickangle(ax_c, 25);
    grid(ax_c, 'on');
    if roi_idx < nrois
        ax_c.XTickLabel = [];
    else
        xlabel(ax_c, 'Condition');
    end
end

sgtitle(sprintf('%s (%s)', figure_label, stim_context.selected_label));
save_figure_bundle(fig, fullfile(save_path, [file_stem '.fig']), fullfile(save_path, [file_stem '.png']));
close(fig);
end

function comparison = build_dual_metric_comparison(metric_name, voltage_traces, calcium_traces, t_voltage, t_calcium, calcium_smoothing_window, save_path, stim_windows, voltage_polarity, calcium_polarity)
% Build the main dual-channel comparison package for one metric type
% (currently sensitivity or SNR). The output combines overlap,
% accumulated-voltage comparison, calcium deconvolution, and ROI-matched
% correlation statistics.
voltage_metric = double(voltage_traces);
calcium_metric = double(calcium_traces);
voltage_display = voltage_polarity * voltage_metric;
calcium_display = calcium_polarity * calcium_metric;

[voltage_integral, calcium_normalized, paired_corr, paired_spearman] = ...
    build_integral_comparison(voltage_metric, calcium_metric, voltage_polarity, calcium_polarity);
[calcium_deconv, calcium_spikes, deconv_info] = build_calcium_deconvolution(calcium_metric);
calcium_deconv_display = calcium_polarity * calcium_deconv;
correlation_stats = compute_dual_correlation_statistics(voltage_integral, calcium_normalized);

quad_fig = fullfile(save_path, sprintf('5_dual_%s_quad_summary.fig', metric_name));
quad_png = fullfile(save_path, sprintf('5_dual_%s_quad_summary.png', metric_name));
quad_roi_dir = fullfile(save_path, sprintf('5_dual_%s_quad_rois', metric_name));
overlap_fig = fullfile(save_path, sprintf('5_dual_%s_overlap.fig', metric_name));
overlap_png = fullfile(save_path, sprintf('5_dual_%s_overlap.png', metric_name));
overlap_roi_dir = fullfile(save_path, sprintf('5_dual_%s_overlap_rois', metric_name));
correlation_fig = fullfile(save_path, sprintf('6_dual_%s_correlation.fig', metric_name));
correlation_png = fullfile(save_path, sprintf('6_dual_%s_correlation.png', metric_name));
mat_file = fullfile(save_path, sprintf('6_dual_%s_comparison.mat', metric_name));

plot_dual_quad_summary( ...
    t_calcium, t_voltage, calcium_normalized, calcium_deconv_display, voltage_integral, voltage_display, ...
    size(calcium_metric, 2), metric_name, stim_windows);
quad_handle = gcf;
save_figure_bundle(quad_handle, quad_fig, quad_png);
close(quad_handle);
save_dual_quad_per_roi( ...
    t_calcium, t_voltage, calcium_normalized, calcium_deconv_display, voltage_integral, voltage_display, ...
    size(calcium_metric, 2), metric_name, stim_windows, quad_roi_dir);

plot_dual_overlap_summary( ...
    t_calcium, t_voltage, calcium_display, voltage_display, ...
    size(calcium_metric, 2), metric_name, stim_windows);
overlap_handle = gcf;
save_figure_bundle_preserve_layout(overlap_handle, overlap_fig, overlap_png);
close(overlap_handle);
save_dual_overlap_per_roi( ...
    t_calcium, t_voltage, calcium_display, voltage_display, ...
    size(calcium_metric, 2), metric_name, stim_windows, overlap_roi_dir);

plot_dual_correlation_summary(correlation_stats, metric_name);
correlation_handle = gcf;
save_figure_bundle(correlation_handle, correlation_fig, correlation_png);
close(correlation_handle);

save(mat_file, ...
    'voltage_metric', 'calcium_metric', 'voltage_display', 'calcium_display', ...
    'voltage_integral', 'calcium_normalized', 'calcium_deconv', 'calcium_spikes', ...
    'paired_corr', 'paired_spearman', 'correlation_stats', 'deconv_info');

fprintf('Dual %s comparison saved | quad=%s | overlap=%s | corr=%s\n', ...
    metric_name, quad_png, overlap_png, correlation_png);

comparison = struct( ...
    'data', struct( ...
        'voltage_metric', voltage_metric, ...
        'calcium_metric', calcium_metric, ...
        'voltage_display', voltage_display, ...
        'calcium_display', calcium_display, ...
        'voltage_integral', voltage_integral, ...
        'calcium_normalized', calcium_normalized, ...
        'calcium_deconv', calcium_deconv, ...
        'calcium_spikes', calcium_spikes, ...
        'paired_corr', paired_corr, ...
        'paired_spearman', paired_spearman, ...
        'correlation_stats', correlation_stats), ...
    'info', struct( ...
        'metric_name', metric_name, ...
        'deconvolution', deconv_info, ...
        'display_smoothing', 'none', ...
        'calcium_smoothing_window', calcium_smoothing_window, ...
        'quad_fig', quad_fig, ...
        'quad_png', quad_png, ...
        'quad_roi_dir', quad_roi_dir, ...
        'overlap_fig', overlap_fig, ...
        'overlap_png', overlap_png, ...
        'overlap_roi_dir', overlap_roi_dir, ...
        'correlation_fig', correlation_fig, ...
        'correlation_png', correlation_png, ...
        'mat_file', mat_file, ...
        'created_at', datetime("now")));
end

function correlation_stats = compute_dual_correlation_statistics(voltage_integral, calcium_normalized)
[paired_corr, shuffled_corr] = compute_paired_vs_shuffled_correlation(voltage_integral, calcium_normalized, 'Pearson');
[paired_spearman, shuffled_spearman] = compute_paired_vs_shuffled_correlation(voltage_integral, calcium_normalized, 'Spearman');

correlation_stats = struct( ...
    'pearson', struct( ...
        'paired', paired_corr, ...
        'shuffled', shuffled_corr, ...
        'p_value', compare_correlation_groups(paired_corr, shuffled_corr)), ...
    'spearman', struct( ...
        'paired', paired_spearman, ...
        'shuffled', shuffled_spearman, ...
        'p_value', compare_correlation_groups(paired_spearman, shuffled_spearman)));
end

function [paired_values, shuffled_values] = compute_paired_vs_shuffled_correlation(voltage_integral, calcium_normalized, corr_type)
nrois = size(voltage_integral, 2);
paired_values = NaN(nrois, 1);

for roi_idx = 1:nrois
    paired_values(roi_idx) = corr(voltage_integral(:, roi_idx), calcium_normalized(:, roi_idx), ...
        'Type', corr_type, 'Rows', 'complete');
end

if nrois > 1
    shuffled_values = NaN(nrois * (nrois - 1), 1);
    write_idx = 1;
    for roi_idx = 1:nrois
        for shuffled_idx = 1:nrois
            if shuffled_idx == roi_idx
                continue;
            end
            shuffled_values(write_idx) = corr(voltage_integral(:, roi_idx), calcium_normalized(:, shuffled_idx), ...
                'Type', corr_type, 'Rows', 'complete');
            write_idx = write_idx + 1;
        end
    end
    shuffled_values = shuffled_values(1:write_idx-1);
else
    nshuffles = 100;
    shuffled_values = NaN(nshuffles, 1);
    for shuffle_idx = 1:nshuffles
        shift_amount = randi([1, max(1, size(calcium_normalized, 1) - 1)], 1, 1);
        shuffled_values(shuffle_idx) = corr( ...
            voltage_integral(:, 1), ...
            circshift(calcium_normalized(:, 1), shift_amount), ...
            'Type', corr_type, 'Rows', 'complete');
    end
end

paired_values = paired_values(isfinite(paired_values));
shuffled_values = shuffled_values(isfinite(shuffled_values));
end

function p_value = compare_correlation_groups(group_a, group_b)
group_a = group_a(isfinite(group_a));
group_b = group_b(isfinite(group_b));

if isempty(group_a) || isempty(group_b)
    p_value = NaN;
    return;
end

try
    [~, p_value] = ttest2(group_a, group_b, 'Vartype', 'unequal');
catch
    p_value = NaN;
end
end

function plot_dual_correlation_summary(correlation_stats, metric_name)
figure('Color', 'w', 'Position', [120, 120, 1100, 500]);

subplot(1, 2, 1);
plot_correlation_distribution(gca, correlation_stats.pearson.paired, correlation_stats.pearson.shuffled, ...
    'Pearson r', correlation_stats.pearson.p_value, [0.85 0.25 0.25]);

subplot(1, 2, 2);
plot_correlation_distribution(gca, correlation_stats.spearman.paired, correlation_stats.spearman.shuffled, ...
    'Spearman \rho', correlation_stats.spearman.p_value, [0.2 0.45 0.85]);

sgtitle(sprintf('Dual %s Correlation: Paired ROI vs Shuffled ROI', upper(metric_name)));
end

function plot_correlation_distribution(ax, paired_values, shuffled_values, metric_label, p_value, point_color)
group_labels = [repmat({'Paired'}, numel(paired_values), 1); repmat({'Shuffled'}, numel(shuffled_values), 1)];
group_values = [paired_values(:); shuffled_values(:)];

boxplot(ax, group_values, group_labels, 'Colors', [0.2 0.2 0.2], 'Symbol', '');
hold(ax, 'on');

x_paired = 1 + 0.10 * (rand(size(paired_values)) - 0.5);
x_shuffled = 2 + 0.10 * (rand(size(shuffled_values)) - 0.5);
scatter(ax, x_paired, paired_values, 28, point_color, 'filled', 'MarkerFaceAlpha', 0.75);
scatter(ax, x_shuffled, shuffled_values, 18, [0.55 0.55 0.55], 'filled', 'MarkerFaceAlpha', 0.35);

y_limits = ylim(ax);
y_text = y_limits(2) - 0.06 * diff(y_limits);
text(ax, 1.5, y_text, sprintf('p = %s', format_p_value(p_value)), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
    'FontSize', 9, 'FontWeight', 'bold');

ylabel(ax, metric_label);
title(ax, sprintf('%s Distribution', metric_label));
grid(ax, 'on');
box(ax, 'off');
end

function [calcium_deconv, calcium_spikes, info] = build_calcium_deconvolution(calcium_traces)
% Build a calcium event-like representation that can be compared more
% directly with voltage after calcium kinetics blur the original response.
nrois = size(calcium_traces, 2);
calcium_deconv = zeros(size(calcium_traces));
calcium_spikes = zeros(size(calcium_traces));
used_foopsi = false;
failed_messages = strings(0, 1);

for roi_idx = 1:nrois
    current_trace = double(calcium_traces(:, roi_idx));
    if exist('deconvolveCa', 'file') == 2
        try
            [c, s] = deconvolveCa(current_trace, 'method', 'foopsi', ...
                'sn', 0.05, 'pars', 0.992, 'optimize_pars', 1, 'optimize_smin', 1);
            calcium_deconv(:, roi_idx) = normalize_columns_to_unit_range(c);
            calcium_spikes(:, roi_idx) = normalize_columns_to_unit_range(s);
            used_foopsi = true;
            continue;
        catch ME
            failed_messages(end+1, 1) = string(ME.message); %#ok<AGROW>
        end
    end

    calcium_deconv(:, roi_idx) = normalize_columns_to_unit_range(current_trace);
    spike_proxy = [0; max(0, diff(current_trace))];
    calcium_spikes(:, roi_idx) = normalize_columns_to_unit_range(spike_proxy);
end

info = struct( ...
    'method', ternary(used_foopsi, 'deconvolveCa_foopsi_with_fallback', 'fallback_positive_diff_proxy'), ...
    'foopsi_available', exist('deconvolveCa', 'file') == 2, ...
    'failed_messages', failed_messages);
end

function plot_dual_quad_summary(t_calcium, t_voltage, calcium_trace, deconv_trace, integral_trace, voltage_trace, nrois, metric_name, stim_windows)
linemaxroi = 3;
plotlines = ceil(nrois / linemaxroi);
axisOpts = {'XTickLabel', [], 'YTickLabel', [], 'TickLength', [0 0], ...
    'XColor', 'none', 'YColor', 'none', 'Box', 'off', 'Color', 'none'};
row_labels = {'Calcium', 'Deconv Ca', 'Accum V', 'Voltage'};
row_colors = {[0.1 0.65 0.2], [0.1 0.35 0.9], [0.75 0.15 0.75], [0.85 0.15 0.15]};
calcium_ylim = compute_shared_ylim(calcium_trace);
deconv_ylim = compute_shared_ylim(deconv_trace);
integral_ylim = compute_shared_ylim(integral_trace);
voltage_ylim = compute_shared_ylim(voltage_trace);

fig = figure('Name', sprintf('Dual %s Quad Summary', upper(metric_name)), ...
    'Color', 'w', 'Units', 'normalized', 'Position', [0.06 0.05 0.88 0.86]);

for i = 1:nrois
    col = mod(i-1, linemaxroi);
    row = floor((i-1) / linemaxroi);
    base_idx = row * 4 * linemaxroi + col + 1;

    ax1 = subplot(plotlines * 4, linemaxroi, base_idx);
    plot(t_calcium, calcium_trace(:, i), 'g', 'LineWidth', 1);
    if nargin >= 9 && isstruct(stim_windows) && isfield(stim_windows, 'calcium') && isfield(stim_windows.calcium, 'stim_time_ranges')
        add_stim_shading(ax1, stim_windows.calcium, stim_windows.condition_index, stim_windows.condition_colors, 0.12, stim_windows.trial_labels, stim_windows.block_labels);
    end
    title(sprintf('ROI %d', i), 'FontSize', 8, 'FontWeight', 'normal');
    xlim(ax1, [min(t_calcium), max(t_calcium)]);
    ylim(ax1, calcium_ylim);
    set(ax1, axisOpts{:});
    add_trace_badge(ax1, row_labels{1}, row_colors{1});
    add_axis_scalebar(ax1, t_calcium, calcium_trace(:, i), row_colors{1}, '');

    ax2 = subplot(plotlines * 4, linemaxroi, base_idx + linemaxroi);
    plot(t_calcium, deconv_trace(:, i), 'b', 'LineWidth', 1);
    if nargin >= 9 && isstruct(stim_windows) && isfield(stim_windows, 'calcium') && isfield(stim_windows.calcium, 'stim_time_ranges')
        add_stim_shading(ax2, stim_windows.calcium, stim_windows.condition_index, stim_windows.condition_colors, 0.12, stim_windows.trial_labels, stim_windows.block_labels);
    end
    xlim(ax2, [min(t_calcium), max(t_calcium)]);
    ylim(ax2, deconv_ylim);
    set(ax2, axisOpts{:});
    add_trace_badge(ax2, row_labels{2}, row_colors{2});
    add_axis_scalebar(ax2, t_calcium, deconv_trace(:, i), row_colors{2}, '');

    ax3 = subplot(plotlines * 4, linemaxroi, base_idx + 2 * linemaxroi);
    plot(t_calcium, integral_trace(:, i), 'm', 'LineWidth', 1);
    if nargin >= 9 && isstruct(stim_windows) && isfield(stim_windows, 'calcium') && isfield(stim_windows.calcium, 'stim_time_ranges')
        add_stim_shading(ax3, stim_windows.calcium, stim_windows.condition_index, stim_windows.condition_colors, 0.12, stim_windows.trial_labels, stim_windows.block_labels);
    end
    xlim(ax3, [min(t_calcium), max(t_calcium)]);
    ylim(ax3, integral_ylim);
    set(ax3, axisOpts{:});
    add_trace_badge(ax3, row_labels{3}, row_colors{3});
    add_axis_scalebar(ax3, t_calcium, integral_trace(:, i), row_colors{3}, '');

    ax4 = subplot(plotlines * 4, linemaxroi, base_idx + 3 * linemaxroi);
    plot(t_voltage, voltage_trace(:, i), 'r', 'LineWidth', 1);
    if nargin >= 9 && isstruct(stim_windows) && isfield(stim_windows, 'voltage') && isfield(stim_windows.voltage, 'stim_time_ranges')
        add_stim_shading(ax4, stim_windows.voltage, stim_windows.condition_index, stim_windows.condition_colors, 0.12, stim_windows.trial_labels, stim_windows.block_labels);
    end
    xlim(ax4, [min(t_voltage), max(t_voltage)]);
    ylim(ax4, voltage_ylim);
    set(ax4, axisOpts{:});
    add_trace_badge(ax4, row_labels{4}, row_colors{4});
    add_axis_scalebar(ax4, t_voltage, voltage_trace(:, i), row_colors{4}, '');
end

sgtitle(sprintf('Dual %s Summary: Calcium / Deconvolution / Accumulated Voltage / Voltage', upper(metric_name)));
set(fig, 'Position', get(0, 'Screensize'));
end

function save_dual_quad_per_roi(t_calcium, t_voltage, calcium_trace, deconv_trace, integral_trace, voltage_trace, nrois, metric_name, stim_windows, roi_save_dir)
% For large ROI sets, save one editable quad figure per ROI so each paired
% calcium/voltage comparison remains easy to inspect during debugging.
if nrois <= 5
    return;
end

if ~exist(roi_save_dir, 'dir')
    mkdir(roi_save_dir);
end

row_labels = {'Calcium', 'Deconv Ca', 'Accum V', 'Voltage'};
row_colors = {[0.1 0.65 0.2], [0.1 0.35 0.9], [0.75 0.15 0.75], [0.85 0.15 0.15]};
calcium_ylim = compute_shared_ylim(calcium_trace);
deconv_ylim = compute_shared_ylim(deconv_trace);
integral_ylim = compute_shared_ylim(integral_trace);
voltage_ylim = compute_shared_ylim(voltage_trace);

for roi_idx = 1:nrois
    fig = figure('Name', sprintf('Dual %s Quad ROI %03d', upper(metric_name), roi_idx), ...
        'Color', 'w', 'Position', [120, 80, 1200, 900]);
    layout = tiledlayout(fig, 4, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

    ax1 = nexttile(layout);
    plot(t_calcium, calcium_trace(:, roi_idx), 'Color', row_colors{1}, 'LineWidth', 1.1);
    apply_dual_quad_roi_axis(ax1, t_calcium, calcium_trace(:, roi_idx), calcium_ylim, row_labels{1}, row_colors{1}, stim_windows, 'calcium');
    title(ax1, sprintf('ROI %d | Dual %s Quad Pair', roi_idx, upper(metric_name)), 'FontSize', 10, 'FontWeight', 'normal');

    ax2 = nexttile(layout);
    plot(t_calcium, deconv_trace(:, roi_idx), 'Color', row_colors{2}, 'LineWidth', 1.1);
    apply_dual_quad_roi_axis(ax2, t_calcium, deconv_trace(:, roi_idx), deconv_ylim, row_labels{2}, row_colors{2}, stim_windows, 'calcium');

    ax3 = nexttile(layout);
    plot(t_calcium, integral_trace(:, roi_idx), 'Color', row_colors{3}, 'LineWidth', 1.1);
    apply_dual_quad_roi_axis(ax3, t_calcium, integral_trace(:, roi_idx), integral_ylim, row_labels{3}, row_colors{3}, stim_windows, 'calcium');

    ax4 = nexttile(layout);
    plot(t_voltage, voltage_trace(:, roi_idx), 'Color', row_colors{4}, 'LineWidth', 1.1);
    apply_dual_quad_roi_axis(ax4, t_voltage, voltage_trace(:, roi_idx), voltage_ylim, row_labels{4}, row_colors{4}, stim_windows, 'voltage');
    xlabel(ax4, 'Time (s)');

    fig_file = fullfile(roi_save_dir, sprintf('ROI_%03d.fig', roi_idx));
    png_file = fullfile(roi_save_dir, sprintf('ROI_%03d.png', roi_idx));
    save_figure_bundle_preserve_layout(fig, fig_file, png_file);
    close(fig);
end
fprintf('Dual %s per-ROI quad figures saved to: %s\n', metric_name, roi_save_dir);
end

function apply_dual_quad_roi_axis(ax, t_axis, y_values, y_limits, trace_label, trace_color, stim_windows, channel_name)
hold(ax, 'on');
if isstruct(stim_windows) && isfield(stim_windows, channel_name) && isfield(stim_windows.(channel_name), 'stim_time_ranges')
    add_stim_shading(ax, stim_windows.(channel_name), stim_windows.condition_index, stim_windows.condition_colors, 0.12, stim_windows.trial_labels, stim_windows.block_labels);
end
plot(ax, t_axis, y_values, 'Color', trace_color, 'LineWidth', 1.1);
xlim(ax, [min(t_axis), max(t_axis)]);
ylim(ax, y_limits);
ylabel(ax, trace_label);
grid(ax, 'on');
box(ax, 'off');
add_trace_badge(ax, trace_label, trace_color);
add_axis_scalebar(ax, t_axis, y_values, trace_color, '');
end

function plot_dual_overlap_summary(t_calcium, t_voltage, calcium_trace, voltage_trace, nrois, metric_name, stim_windows)
xlimit = [0, max([t_calcium(:); t_voltage(:)])];
roi_per_row = 3;
ncols = min(roi_per_row, max(1, nrois));
nrows = max(1, ceil(nrois / roi_per_row));

fig = figure('Name', sprintf('Dual %s Overlap', upper(metric_name)), ...
    'Color', 'w', 'Position', [100, 100, max(1200, 420 * ncols), max(420, 260 * nrows)]);

for i = 1:nrois
    ax = subplot(nrows, ncols, i);
    plot_single_dual_overlap_axis( ...
        ax, t_calcium, t_voltage, calcium_trace(:, i), voltage_trace(:, i), i, xlimit, stim_windows);

    if i <= (nrows - 1) * ncols
        ax.XTickLabel = [];
    else
        xlabel(ax, 'Time (s)');
    end
end

sgtitle(sprintf('Dual %s Overlap: Dual yyaxis, separately normalized', upper(metric_name)));
end

function save_dual_overlap_per_roi(t_calcium, t_voltage, calcium_trace, voltage_trace, nrois, metric_name, stim_windows, output_dir)
if ~isfolder(output_dir)
    mkdir(output_dir);
end
xlimit = [0, max([t_calcium(:); t_voltage(:)])];

for roi_idx = 1:nrois
    fig = figure('Name', sprintf('Dual %s Overlap ROI %d', upper(metric_name), roi_idx), ...
        'Color', 'w', 'Position', [100, 100, 1400, 360]);
    ax = axes(fig);
    plot_single_dual_overlap_axis( ...
        ax, t_calcium, t_voltage, calcium_trace(:, roi_idx), voltage_trace(:, roi_idx), roi_idx, xlimit, stim_windows);
    xlabel(ax, 'Time (s)');
    fig_file = fullfile(output_dir, sprintf('ROI_%03d.fig', roi_idx));
    png_file = fullfile(output_dir, sprintf('ROI_%03d.png', roi_idx));
    save_figure_bundle_preserve_layout(fig, fig_file, png_file);
    close(fig);
end
end

function plot_single_dual_overlap_axis(ax, t_calcium, t_voltage, calcium_trace_one, voltage_trace_one, roi_idx, xlimit, stim_windows)
hold(ax, 'on');
if isstruct(stim_windows) && isfield(stim_windows, 'voltage') && isfield(stim_windows.voltage, 'stim_time_ranges')
    add_stim_shading(ax, stim_windows.voltage, stim_windows.condition_index, stim_windows.condition_colors, 0.12, stim_windows.trial_labels, stim_windows.block_labels);
end

voltage_norm = normalize_columns_to_unit_range(voltage_trace_one);
calcium_norm = normalize_columns_to_unit_range(calcium_trace_one);

yyaxis(ax, 'left');
plot(ax, t_voltage, voltage_norm, 'r', 'LineWidth', 1);
ylim(ax, [0 1]);
ax.YColor = [0.85 0.15 0.15];
ylabel(ax, 'Voltage (norm)');

yyaxis(ax, 'right');
plot(ax, t_calcium, calcium_norm, 'g', 'LineWidth', 1);
ylim(ax, [0 1]);
ax.YColor = [0.1 0.65 0.2];
ylabel(ax, 'Calcium (norm)');

xlim(ax, xlimit);
ax.Box = 'off';
ax.Color = 'none';
ax.TickLength = [0.01 0.01];
add_overlap_roi_corner_label(ax, roi_idx);
end

function add_overlap_roi_corner_label(ax, roi_idx)
text(ax, 0.03, 0.95, sprintf('ROI %d', roi_idx), ...
    'Units', 'normalized', ...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'top', ...
    'FontSize', 9, ...
    'FontWeight', 'bold', ...
    'Color', [0.1 0.1 0.1], ...
    'BackgroundColor', [1 1 1], ...
    'Margin', 2, ...
    'Interpreter', 'none');
end

function trace_out = normalize_columns_to_unit_range(trace_in)
trace_in = double(trace_in);
if isvector(trace_in)
    trace_in = trace_in(:);
end

trace_out = zeros(size(trace_in));
for i = 1:size(trace_in, 2)
    current = trace_in(:, i);
    current_min = min(current, [], 'omitnan');
    current_max = max(current, [], 'omitnan');
    if ~isfinite(current_min) || ~isfinite(current_max) || abs(current_max - current_min) < eps
        trace_out(:, i) = zeros(size(current));
    else
        trace_out(:, i) = (current - current_min) ./ (current_max - current_min);
    end
end
end

function plot_dual_stacked(voltage_traces, calcium_traces, t_voltage, t_calcium, nrois, stim_windows)
linemaxroi = 3;
plotlines = ceil(nrois / linemaxroi);
xlimit = [0, max([t_voltage(:); t_calcium(:)])];
ylimitv = [min(voltage_traces, [], 'all'), max(voltage_traces, [], 'all')];
ylimitc = [min(calcium_traces, [], 'all'), max(calcium_traces, [], 'all')];

figure('Color', 'w');
axisOpts = {'XTickLabel', [], 'YTickLabel', [], 'TickLength', [0 0], ...
            'XColor', 'none', 'YColor', 'none', 'Box', 'off', 'Color', 'none'};

for i = 0:nrois-1
    row = floor(i / linemaxroi);
    col = mod(i, linemaxroi);
    calcium_idx = row * 2 * linemaxroi + col + 1;
    voltage_idx = calcium_idx + linemaxroi;

    ax_voltage = subplot(plotlines * 2, linemaxroi, voltage_idx);
    plot(t_voltage, voltage_traces(:, i+1), 'r');
    xlim(xlimit);
    ylim(ylimitv);
    if nargin >= 6 && isstruct(stim_windows) && isfield(stim_windows, 'voltage')
        add_stim_shading(ax_voltage, stim_windows.voltage, stim_windows.condition_index, stim_windows.condition_colors, 0.12, stim_windows.trial_labels, stim_windows.block_labels);
    end
    set(ax_voltage, 'YDir', 'reverse', axisOpts{:});
    title(ax_voltage, sprintf('ROI %d', i+1), 'FontSize', 7, 'FontWeight', 'normal');

    ax_calcium = subplot(plotlines * 2, linemaxroi, calcium_idx);
    plot(t_calcium, calcium_traces(:, i+1), 'g');
    xlim(xlimit);
    ylim(ylimitc);
    if nargin >= 6 && isstruct(stim_windows) && isfield(stim_windows, 'calcium')
        add_stim_shading(ax_calcium, stim_windows.calcium, stim_windows.condition_index, stim_windows.condition_colors, 0.12, stim_windows.trial_labels, stim_windows.block_labels);
    end
    set(ax_calcium, axisOpts{:});
end
end

function [fig_filename, png_filename] = plot_dual_bleach_overview( ...
    traces_voltage_bleach_removed, traces_voltage_input, baseline_voltage, t_voltage, bleach_mode_voltage, ...
    traces_calcium_bleach_removed, traces_calcium_input, baseline_calcium, t_calcium, bleach_mode_calcium, ...
    save_path)

fig = figure('Name', 'Dual Bleaching Correction Overview', 'Color', 'w');
set(fig, 'Position', get(0, 'Screensize'));

ax11 = subplot(2, 2, 1);
plot_bleach_panel_fit( ...
    traces_voltage_input, baseline_voltage, t_voltage, ...
    sprintf('Voltage Original And Fitted Baseline (%s)', bleach_mode_voltage));

ax12 = subplot(2, 2, 2);
plot_bleach_panel_residual( ...
    traces_voltage_bleach_removed, t_voltage, ...
    sprintf('Voltage Bleach-Removed Residuals (%s)', bleach_mode_voltage));

ax21 = subplot(2, 2, 3);
plot_bleach_panel_fit( ...
    traces_calcium_input, baseline_calcium, t_calcium, ...
    sprintf('Calcium Original And Fitted Baseline (%s)', bleach_mode_calcium));

ax22 = subplot(2, 2, 4);
plot_bleach_panel_residual( ...
    traces_calcium_bleach_removed, t_calcium, ...
    sprintf('Calcium Bleach-Removed Residuals (%s)', bleach_mode_calcium));

linkaxes([ax11, ax12], 'x');
linkaxes([ax21, ax22], 'x');

fig_filename = fullfile(save_path, '2_dual_bleach_correction_stacked.fig');
png_filename = fullfile(save_path, '2_dual_bleach_correction_stacked.png');
save_figure_bundle(fig, fig_filename, png_filename);
close(fig);
end

function plot_bleach_panel_fit(traces_input, baseline, t_axis, panel_title)
nrois = size(traces_input, 2);
hold on;
spacing_raw = mean(std(traces_input, 0, 1)) * 5;
if ~isfinite(spacing_raw) || spacing_raw <= 0
    spacing_raw = 1;
end

for r = 1:nrois
    offset = (nrois - r) * spacing_raw;
    plot(t_axis, traces_input(:, r) + offset, 'LineWidth', 0.5);
    plot(t_axis, baseline(:, r) + offset, 'r', 'LineWidth', 1.2);

    if mod(r, 5) == 0 || r == 1 || r == nrois
        text(t_axis(1), offset, [' ROI ', num2str(r)], 'FontSize', 8, 'FontWeight', 'bold');
    end
end

title(panel_title);
xlabel('Time (s)');
ylabel('Stacked Magnitude');
grid on;
axis tight;
end

function plot_bleach_panel_residual(traces_bleach_removed, t_axis, panel_title)
nrois = size(traces_bleach_removed, 2);
hold on;
spacing_corr = mean(std(traces_bleach_removed, 0, 1)) * 8;
if ~isfinite(spacing_corr) || spacing_corr <= 0
    spacing_corr = 1;
end

for r = 1:nrois
    offset = (nrois - r) * spacing_corr;
    plot(t_axis, traces_bleach_removed(:, r) + offset, 'k', 'LineWidth', 0.5);
    line([t_axis(1) t_axis(end)], [offset offset], 'Color', [0.3 0.7 1], 'LineStyle', '--');
end

title(panel_title);
xlabel('Time (s)');
ylabel('Stacked Magnitude');
grid on;
axis tight;
end

function [fig_filename, png_filename] = plot_dual_noise_reference_comparison( ...
    traces_voltage_bleach, voltage_noise_reference, ...
    traces_calcium_bleach, calcium_noise_reference, ...
    t_voltage, t_calcium, save_path)

fig = figure('Name', 'Dual Noise Reference Comparison', 'Color', 'w', ...
    'Position', [150, 150, 1200, 900]);

subplot(2, 1, 1);
plot_trace_reference_comparison(traces_voltage_bleach, voltage_noise_reference, t_voltage, 'Voltage: Bleach Removed vs Noise Reference');

subplot(2, 1, 2);
plot_trace_reference_comparison(traces_calcium_bleach, calcium_noise_reference, t_calcium, 'Calcium: Bleach Removed vs Noise Reference');

fig_filename = fullfile(save_path, '3_dual_noise_reference_comparison.fig');
png_filename = fullfile(save_path, '3_dual_noise_reference_comparison.png');
save_figure_bundle(fig, fig_filename, png_filename);
close(fig);
end

function plot_trace_reference_comparison(signal_traces, reference_traces, t_axis, panel_title)
nrois = size(signal_traces, 2);
max_display = min(10, nrois);
roi_to_show = round(linspace(1, nrois, max_display));
spacing = max(signal_traces(:)) * 0.8;
if ~isfinite(spacing) || spacing <= 0
    spacing = 1;
end

hold on;
for i = 1:length(roi_to_show)
    r_idx = roi_to_show(i);
    offset = (i - 1) * spacing;
    plot(t_axis, signal_traces(:, r_idx) + offset, 'Color', [0.7 0.7 0.7], 'HandleVisibility', 'off');
    plot(t_axis, reference_traces(:, r_idx) + offset, 'LineWidth', 1);
end

set(gca, 'YTick', (0:max_display-1) * spacing, 'YTickLabel', string(roi_to_show));
title(panel_title);
xlabel('Time (s)');
ylabel('ROI Index');
legend('Noise Reference');
grid on;
end

function label_out = capitalize_label(label_in)
label_in = char(string(label_in));
if isempty(label_in)
    label_out = label_in;
else
    label_out = [upper(label_in(1)), label_in(2:end)];
end
end

function [voltage_integral, calcium_normalized, paired_corr, paired_spearman] = build_integral_comparison(voltage_traces, calcium_traces, voltage_polarity, calcium_polarity)
% Convert voltage into a leaky accumulated signal so it can be compared
% more fairly with slower calcium dynamics. This is a comparison bridge,
% not a claim of exact biophysical reconstruction.
nrois = size(voltage_traces, 2);
calcium_normalized = zeros(size(calcium_traces));
voltage_integral = zeros(size(calcium_traces));
paired_corr = NaN(nrois, 1);
paired_spearman = NaN(nrois, 1);

g = 20;
k = 0.8;
gthr = 0.5;
kth = 0.21;
dt = 1 / 200;
ftest = @(x, v) x + (g * max(v - gthr, 0) - k * max(x - kth, 0)) * dt;

for i = 1:nrois
    v = normalize(voltage_polarity * voltage_traces(:, i), 'range');
    c = normalize(calcium_polarity * calcium_traces(:, i), 'range');
    accum = zeros(length(v), 1);

    for j = 2:length(v)
        accum(j) = ftest(accum(j-1), v(j));
    end

    window_size = max(1, round(length(accum) / length(c)));
    accum_smoothed = movmean(accum, window_size);
    accum_resampled = interp1(linspace(1, length(accum), length(accum)), accum_smoothed, ...
        linspace(1, length(accum), length(c)))';

    voltage_integral(:, i) = normalize(accum_resampled, 'range');
    calcium_normalized(:, i) = c;
    paired_corr(i) = corr(voltage_integral(:, i), calcium_normalized(:, i), 'Rows', 'complete');
    paired_spearman(i) = corr(voltage_integral(:, i), calcium_normalized(:, i), 'Type', 'Spearman', 'Rows', 'complete');
end
end

function [dual_info, voltage_results, calcium_results, dual_results, stim_results, stim_context, stim_windows] = load_saved_analysis_only_context(reuse_results_path)
% Reload one finished Dual_analysis3 output folder and reuse everything up
% through ROI selection so the later trace-analysis sections can run
% without repeating movie loading, motion correction, or ROI drawing.
if strlength(string(reuse_results_path)) == 0
    error('analysis_mode=analysis_only requires reuse_results_path.');
end

fprintf('Reusing saved results from: %s\n', reuse_results_path);

dual_info_path = fullfile(reuse_results_path, 'dual_info.mat');
voltage_results_path = fullfile(reuse_results_path, 'voltage_results.mat');
calcium_results_path = fullfile(reuse_results_path, 'calcium_results.mat');
dual_results_path = fullfile(reuse_results_path, 'dual_results.mat');
stim_results_path = fullfile(reuse_results_path, 'stim_results.mat');

dual_info = load_required_struct(dual_info_path, 'dual_info');
voltage_results = load_required_struct(voltage_results_path, 'voltage_results');
calcium_results = load_required_struct(calcium_results_path, 'calcium_results');
if isfile(dual_results_path)
    dual_results = load_required_struct(dual_results_path, 'dual_results');
else
    dual_results = struct();
end
if isfile(stim_results_path)
    stim_results = load_required_struct(stim_results_path, 'stim_results');
else
    stim_results = struct();
end

stim_context = dual_info.stim_context;
stim_windows = struct();
if isstruct(stim_results) && isfield(stim_results, 'windows')
    stim_windows = stim_results.windows;
end
fprintf('Analysis-only reuse loaded: saved ROI/channel/stim results are available for downstream sections.\n');
end

function [dual_info, voltage_results, calcium_results, dual_results, stim_results] = seed_analysis_only_output_folder(save_path, reuse_results_path, dual_info, voltage_results, calcium_results, dual_results, stim_results)
% Create an analysis-only fork: load traces from the old result folder, but
% write all rerun outputs into a new folder. The old folder stays read-only
% from this point onward unless the caller explicitly uses it as save_path.
if ~isfolder(save_path)
    mkdir(save_path);
end

source_roi_file = "";
if isstruct(dual_results) && isfield(dual_results, 'registration') ...
        && isfield(dual_results.registration, 'info') ...
        && isfield(dual_results.registration.info, 'roi_file')
    source_roi_file = string(dual_results.registration.info.roi_file);
end

fork_roi_file = "";
if strlength(source_roi_file) > 0 && isfile(char(source_roi_file))
    fork_roi_file = string(fullfile(save_path, '1_dual_roi_results.mat'));
    if ~strcmpi(char(source_roi_file), char(fork_roi_file))
        copyfile(char(source_roi_file), char(fork_roi_file));
    end
    dual_results.registration.info.source_roi_file = source_roi_file;
    dual_results.registration.info.roi_file = char(fork_roi_file);
end

if isstruct(dual_info)
    if isfield(dual_info, 'save_path')
        dual_info.source_save_path = string(dual_info.save_path);
    else
        dual_info.source_save_path = string(reuse_results_path);
    end
    dual_info.save_path = char(string(save_path));
    dual_info.analysis_only_source_path = string(reuse_results_path);
    dual_info.analysis_only_fork_created_at = datetime("now");
end

save(fullfile(save_path, 'dual_info.mat'), 'dual_info');
save(fullfile(save_path, 'voltage_results.mat'), 'voltage_results', '-v7.3');
save(fullfile(save_path, 'calcium_results.mat'), 'calcium_results', '-v7.3');
save(fullfile(save_path, 'dual_results.mat'), 'dual_results', '-v7.3');
save(fullfile(save_path, 'stim_results.mat'), 'stim_results', '-v7.3');

fprintf('Analysis-only fork initialized.\n');
fprintf('  source results: %s\n', reuse_results_path);
fprintf('  output folder : %s\n', save_path);
if strlength(fork_roi_file) > 0
    fprintf('  ROI context copied to: %s\n', fork_roi_file);
end
end

function save_figure_bundle(fig_handle, fig_file, png_file)
% Save the editable FIG plus a large PNG snapshot. PNG export is done
% after maximizing the figure so saved images match the on-screen layout
% more closely than the MATLAB default small window capture.
prepare_figure_for_png_export(fig_handle);
save_figure_outputs_resilient(fig_handle, fig_file, png_file, true);
end

function save_figure_bundle_preserve_layout(fig_handle, fig_file, png_file)
if ~isgraphics(fig_handle, 'figure')
    return;
end
save_figure_outputs_resilient(fig_handle, fig_file, png_file, false);
end

function save_figure_outputs_resilient(fig_handle, fig_file, png_file, prepared_for_export)
if nargin < 4
    prepared_for_export = false;
end
if ~isgraphics(fig_handle, 'figure')
    warning('Dual_analysis3:InvalidFigureHandle', ...
        'Figure output skipped because the figure handle is invalid: %s', png_file);
    return;
end

try
    saveas(fig_handle, fig_file, 'fig');
catch ME
    warning('Dual_analysis3:SaveFigFailed', ...
        'Failed to save FIG file %s: %s', fig_file, ME.message);
end

if ~isgraphics(fig_handle, 'figure')
    warning('Dual_analysis3:FigureClosedBeforePng', ...
        'PNG export skipped because the figure closed before export: %s', png_file);
    return;
end

try
    drawnow;
    exportgraphics(fig_handle, png_file, 'Resolution', 150);
catch ME
    warning('Dual_analysis3:ExportGraphicsFailed', ...
        'exportgraphics failed for %s: %s. Trying saveas PNG fallback.', png_file, ME.message);
    if isgraphics(fig_handle, 'figure')
        try
            if ~prepared_for_export
                drawnow;
            end
            saveas(fig_handle, png_file, 'png');
        catch ME2
            warning('Dual_analysis3:SavePngFallbackFailed', ...
                'PNG fallback also failed for %s: %s', png_file, ME2.message);
        end
    end
end
end

function prepare_figure_for_png_export(fig_handle)
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

function time_frequency_results = run_dual_time_frequency_section( ...
    voltage_results, calcium_results, ...
    t_voltage_default, t_calcium_default, ...
    freq_voltage_default, freq_calcium_default, ...
    stim_windows, save_path, ...
    voltage_polarity, calcium_polarity)
% Build the final spectral summary directly from saved processed traces.
% This keeps the last section independent from the earlier movie-based
% steps and makes iterative figure changes much cheaper to rerun.
[voltage_trace, voltage_stage, t_voltage, freq_voltage] = resolve_time_frequency_channel_input( ...
    voltage_results, {'snr', 'sensitivity', 'bleach_removed'}, t_voltage_default, freq_voltage_default);
[calcium_trace, calcium_stage, t_calcium, freq_calcium] = resolve_time_frequency_channel_input( ...
    calcium_results, {'snr', 'sensitivity', 'bleach_removed'}, t_calcium_default, freq_calcium_default);

params = struct();
params.voltage_stage = string(voltage_stage);
params.calcium_stage = string(calcium_stage);
params.min_freq_hz = 0.5;
params.max_freq_hz = min([80, freq_voltage / 2 - eps, freq_calcium / 2 - eps]);
params.wavelet_voices_per_octave = 12;
params.wavelet_name = "amor";

fprintf('Time-frequency input stages | voltage=%s | calcium=%s\n', ...
    string(voltage_stage), string(calcium_stage));
fprintf('Time-frequency frequency range | min=%.2f Hz | max=%.2f Hz\n', ...
    params.min_freq_hz, params.max_freq_hz);
fprintf('Time-frequency settings | spectrum=direct FFT | wavelet=%s\n', ...
    params.wavelet_name);

cleanup_legacy_oscillation_outputs(save_path);

voltage_display = double(voltage_polarity) * double(voltage_trace);
calcium_display = double(calcium_polarity) * double(calcium_trace);

voltage_tf = analyze_population_time_frequency( ...
    voltage_display, freq_voltage, t_voltage, params, 'Voltage');
calcium_tf = analyze_population_time_frequency( ...
    calcium_display, freq_calcium, t_calcium, params, 'Calcium');

[fourier_fig, fourier_png] = plot_population_fourier_summary( ...
    voltage_tf, calcium_tf, stim_windows, params, save_path);
[wavelet_fig, wavelet_png, wavelet_roi_dir] = plot_population_wavelet_summary( ...
    voltage_tf, calcium_tf, stim_windows, save_path);

time_frequency_results = struct( ...
    'parameters', params, ...
    'voltage', voltage_tf, ...
    'calcium', calcium_tf, ...
    'visualizations', struct( ...
        'fourier_summary', struct('fig_file', fourier_fig, 'png_file', fourier_png), ...
        'wavelet_summary', struct('fig_file', wavelet_fig, 'png_file', wavelet_png, 'roi_pair_dir', wavelet_roi_dir)));
end

function result = analyze_population_time_frequency(trace_matrix, frame_rate, t_axis, params, channel_name)
trace_matrix = double(trace_matrix);
t_axis = double(t_axis(:));
nframes = size(trace_matrix, 1);
nrois = size(trace_matrix, 2);

spectrum_frequency = [];
spectrum_amplitude = [];
peak_frequency_hz = NaN(nrois, 1);
peak_amplitude = NaN(nrois, 1);
signal_rms = NaN(nrois, 1);

for roi_idx = 1:nrois
    x = sanitize_trace_for_spectrum(trace_matrix(:, roi_idx));
    signal_rms(roi_idx) = rms(x);
    [frequency_i, amplitude_i] = compute_trace_fft_spectrum(x, frame_rate);
    if isempty(spectrum_frequency)
        spectrum_frequency = frequency_i;
        spectrum_amplitude = NaN(numel(frequency_i), nrois);
    end
    if isempty(amplitude_i)
        continue;
    end
    spectrum_amplitude(:, roi_idx) = amplitude_i;

    valid_mask = frequency_i >= params.min_freq_hz & frequency_i <= params.max_freq_hz;
    valid_frequency = frequency_i(valid_mask);
    valid_amplitude = amplitude_i(valid_mask);
    if isempty(valid_frequency)
        continue;
    end

    [peak_value, max_idx] = max(valid_amplitude);
    peak_frequency_hz(roi_idx) = valid_frequency(max_idx);
    peak_amplitude(roi_idx) = peak_value;
end

representative_roi = find(peak_amplitude == max(peak_amplitude, [], 'omitnan'), 1, 'first');
if isempty(representative_roi)
    representative_roi = 1;
end
representative_trace = sanitize_trace_for_spectrum(trace_matrix(:, representative_roi));
representative_wavelet = compute_wavelet_scalogram(representative_trace, frame_rate, t_axis, params);

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
        'trace', representative_trace, ...
        'peak_frequency_hz', peak_frequency_hz(representative_roi), ...
        'wavelet', representative_wavelet), ...
    'summary', struct( ...
        'nrois', nrois, ...
        'median_peak_frequency_hz', median(peak_frequency_hz, 'omitnan'), ...
        'median_peak_amplitude', median(peak_amplitude, 'omitnan')));
end

function x = sanitize_trace_for_spectrum(x)
x = double(x(:));
if isempty(x) || all(~isfinite(x))
    x = zeros(size(x));
    return;
end
x = fillmissing(x, 'linear', 'EndValues', 'nearest');
x = x - mean(x, 'omitnan');
end

function [frequency, amplitude] = compute_trace_fft_spectrum(x, frame_rate)
x = sanitize_trace_for_spectrum(x);
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

function wavelet_data = compute_wavelet_scalogram(x, frame_rate, t_axis, params)
x = sanitize_trace_for_spectrum(x);
if numel(x) < 8
    wavelet_data = struct('time', double(t_axis(:)), 'frequency', zeros(0, 1), 'power', zeros(0, 0));
    return;
end

[wt, frequency] = cwt(x, frame_rate);
valid_mask = frequency >= params.min_freq_hz & frequency <= params.max_freq_hz;
if numel(t_axis) ~= numel(x)
    t_axis = (0:numel(x) - 1)' / frame_rate;
end
wavelet_data = struct( ...
    'time', double(t_axis(:)), ...
    'frequency', frequency(valid_mask), ...
    'power', abs(wt(valid_mask, :)).^2);
end

function [fig_file, png_file] = plot_population_fourier_summary(voltage_tf, calcium_tf, stim_windows, params, save_path)
fig = figure('Color', 'w', 'Name', 'Population FFT Summary');
tiledlayout(fig, 4, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
plot_channel_fft_summary(voltage_tf, 'r', 'Voltage FFT Amplitude', [0, params.max_freq_hz], false);

nexttile;
plot_channel_fft_summary(calcium_tf, 'g', 'Calcium FFT Amplitude', [0, params.max_freq_hz], false);

nexttile;
plot_channel_fft_summary(voltage_tf, 'r', 'Voltage FFT Amplitude', [0, 10], true);

nexttile;
plot_channel_fft_summary(calcium_tf, 'g', 'Calcium FFT Amplitude', [0, 10], true);

nexttile;
plot_flash_fft_comparison( ...
    voltage_tf, resolve_channel_stim_windows(stim_windows, 'voltage'), params, 'r', 'Voltage Flash FFT');

nexttile;
plot_flash_fft_comparison( ...
    calcium_tf, resolve_channel_stim_windows(stim_windows, 'calcium'), params, 'g', 'Calcium Flash FFT');

nexttile;
plot_roi_peak_summary(voltage_tf, [0, 10], 'Voltage ROI Peak Frequency (0-10 Hz)');

nexttile;
plot_roi_peak_summary(calcium_tf, [0, 10], 'Calcium ROI Peak Frequency (0-10 Hz)');

fig_file = fullfile(save_path, '8_fourier_summary.fig');
png_file = fullfile(save_path, '8_fourier_summary.png');
save_figure_bundle(fig, fig_file, png_file);
close(fig);
end

function [fig_file, png_file, roi_pair_dir] = plot_population_wavelet_summary(voltage_tf, calcium_tf, stim_windows, save_path)
fig = figure('Color', 'w', 'Name', 'Population Wavelet Summary');
nrois = size(voltage_tf.trace_matrix, 2);
voltage_windows = resolve_channel_stim_windows(stim_windows, 'voltage');
calcium_windows = resolve_channel_stim_windows(stim_windows, 'calcium');
tiledlayout(fig, nrois, 4, 'Padding', 'compact', 'TileSpacing', 'compact');

for roi_idx = 1:nrois
    nexttile;
    plot_roi_trace_with_stim( ...
        gca, voltage_tf, roi_idx, voltage_windows, 'r', 'Voltage Trace');

    nexttile;
    plot_roi_wavelet_scalogram( ...
        gca, voltage_tf, roi_idx, voltage_windows, 'Voltage Wavelet');

    nexttile;
    plot_roi_trace_with_stim( ...
        gca, calcium_tf, roi_idx, calcium_windows, 'g', 'Calcium Trace');

    nexttile;
    plot_roi_wavelet_scalogram( ...
        gca, calcium_tf, roi_idx, calcium_windows, 'Calcium Wavelet');
end

fig_file = fullfile(save_path, '8_wavelet_summary.fig');
png_file = fullfile(save_path, '8_wavelet_summary.png');
save_figure_bundle(fig, fig_file, png_file);
close(fig);

roi_pair_dir = "";
if nrois > 5
    roi_pair_dir = fullfile(save_path, '8_wavelet_roi_pairs');
    save_wavelet_pair_per_roi(voltage_tf, calcium_tf, voltage_windows, calcium_windows, roi_pair_dir);
end
end

function save_wavelet_pair_per_roi(voltage_tf, calcium_tf, voltage_windows, calcium_windows, roi_pair_dir)
% Save each ROI as its own voltage/calcium trace-plus-wavelet pair. This
% keeps the original population wavelet overview while making large ROI
% sets easier to inspect one matched ROI at a time.
if ~exist(roi_pair_dir, 'dir')
    mkdir(roi_pair_dir);
end

nrois = size(voltage_tf.trace_matrix, 2);
for roi_idx = 1:nrois
    fig = figure('Color', 'w', 'Name', sprintf('Wavelet ROI Pair %03d', roi_idx), ...
        'Position', [120, 80, 1400, 760]);
    layout = tiledlayout(fig, 2, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

    nexttile(layout);
    plot_roi_trace_with_stim(gca, voltage_tf, roi_idx, voltage_windows, 'r', 'Voltage Trace');
    title(sprintf('ROI %d | Voltage Trace', roi_idx), 'FontSize', 10, 'FontWeight', 'normal');

    nexttile(layout);
    plot_roi_wavelet_scalogram(gca, voltage_tf, roi_idx, voltage_windows, 'Voltage Wavelet');

    nexttile(layout);
    plot_roi_trace_with_stim(gca, calcium_tf, roi_idx, calcium_windows, 'g', 'Calcium Trace');
    title(sprintf('ROI %d | Calcium Trace', roi_idx), 'FontSize', 10, 'FontWeight', 'normal');

    nexttile(layout);
    plot_roi_wavelet_scalogram(gca, calcium_tf, roi_idx, calcium_windows, 'Calcium Wavelet');

    fig_file = fullfile(roi_pair_dir, sprintf('ROI_%03d.fig', roi_idx));
    png_file = fullfile(roi_pair_dir, sprintf('ROI_%03d.png', roi_idx));
    save_figure_bundle_preserve_layout(fig, fig_file, png_file);
    close(fig);
end
fprintf('Wavelet per-ROI pair figures saved to: %s\n', roi_pair_dir);
end

function plot_channel_fft_summary(channel_result, trace_color, title_text, display_band_hz, annotate_peak)
frequency = channel_result.fft.frequency;
amplitude = channel_result.fft.amplitude;
if isempty(frequency) || isempty(amplitude)
    text(0.5, 0.5, 'FFT unavailable', 'Units', 'normalized', 'HorizontalAlignment', 'center');
    axis off;
    return;
end

display_mask = frequency >= display_band_hz(1) & frequency <= display_band_hz(2);
if ~any(display_mask)
    text(0.5, 0.5, 'No spectral points in display band', 'Units', 'normalized', 'HorizontalAlignment', 'center');
    axis off;
    return;
end

plot(frequency(display_mask), amplitude(display_mask, :), 'Color', [0.82 0.82 0.82], 'LineWidth', 0.8);
hold on;
mean_amplitude = mean(amplitude, 2, 'omitnan');
plot(frequency(display_mask), mean_amplitude(display_mask), trace_color, 'LineWidth', 2);
if annotate_peak
    [display_peak_hz, display_peak_amplitude] = find_band_peak_from_amplitude( ...
        frequency, mean_amplitude, display_band_hz);
    if isfinite(display_peak_hz)
        xline(display_peak_hz, '--', sprintf('%.2f Hz', display_peak_hz), ...
            'Color', trace_color, 'LineWidth', 1.2, 'LabelVerticalAlignment', 'middle');
        scatter(display_peak_hz, display_peak_amplitude, 42, trace_color, 'filled');
        text(display_peak_hz, display_peak_amplitude, sprintf('  %.2f Hz', display_peak_hz), ...
        'Color', trace_color, 'FontSize', 9, 'VerticalAlignment', 'bottom');
    end
end
xlabel('Frequency (Hz)');
ylabel('Single-Sided Amplitude');
if annotate_peak
    title(sprintf('%s | 0-10 Hz zoom', title_text));
else
    title(sprintf('%s | full range', title_text));
end
xlim(display_band_hz);
grid on;
end

function plot_roi_peak_summary(channel_result, display_band_hz, title_text)
peak_frequency_hz = NaN(size(channel_result.roi.peak_frequency_hz(:)));
peak_amplitude = NaN(size(channel_result.roi.peak_amplitude(:)));
roi_index = (1:numel(peak_frequency_hz))';
for roi_idx = 1:numel(peak_frequency_hz)
    [peak_frequency_hz(roi_idx), peak_amplitude(roi_idx)] = find_band_peak_from_amplitude( ...
        channel_result.fft.frequency, channel_result.fft.amplitude(:, roi_idx), display_band_hz);
end
valid = isfinite(peak_frequency_hz) & isfinite(peak_amplitude);
if ~any(valid)
    text(0.5, 0.5, 'ROI peak summary unavailable', 'Units', 'normalized', 'HorizontalAlignment', 'center');
    axis off;
    return;
end

scatter(roi_index(valid), peak_frequency_hz(valid), 42, peak_amplitude(valid), 'filled', 'MarkerFaceAlpha', 0.8);
hold on;
plot(roi_index(valid), peak_frequency_hz(valid), '-', 'Color', [0.75 0.75 0.75]);
xlabel('ROI');
ylabel('Peak Frequency (Hz)');
ylim(display_band_hz);
title(title_text);
cb = colorbar;
cb.Label.String = 'Peak Amplitude';
grid on;
end

function plot_representative_trace_with_stim(ax, channel_result, channel_windows, trace_color, title_text)
t = channel_result.time(:);
x = channel_result.representative.trace(:);
[y_label, scale_suffix, stage_label] = resolve_time_frequency_trace_label(channel_result);
plot(ax, t, x, 'Color', trace_color, 'LineWidth', 1.2);
apply_trace_axis_limits(ax, x);
if isstruct(channel_windows) && ~isempty(fieldnames(channel_windows))
    nshades = count_channel_shading_ranges(channel_windows);
    add_stim_shading( ...
        ax, channel_windows, ...
        ones(max(1, nshades), 1), ...
        repmat([0.75 0.75 0.75], max(1, nshades), 1), ...
        0.16);
    hold(ax, 'on');
end
plot(ax, t, x, 'Color', trace_color, 'LineWidth', 1.2);
xlabel(ax, 'Time (s)');
ylabel(ax, y_label);
title(ax, sprintf('%s (%s) | ROI %d | peak %.2f Hz', ...
    title_text, stage_label, channel_result.representative.roi_index, channel_result.representative.peak_frequency_hz));
grid(ax, 'on');
add_axis_scalebar(ax, t, x, trace_color, scale_suffix);
end

function plot_roi_trace_with_stim(ax, channel_result, roi_idx, channel_windows, trace_color, title_text)
t = channel_result.time(:);
x = double(channel_result.trace_matrix(:, roi_idx));
[y_label, scale_suffix, stage_label] = resolve_time_frequency_trace_label(channel_result);
plot(ax, t, x, 'Color', trace_color, 'LineWidth', 1.0);
apply_trace_axis_limits(ax, x);
if isstruct(channel_windows) && ~isempty(fieldnames(channel_windows))
    nshades = count_channel_shading_ranges(channel_windows);
    add_stim_shading( ...
        ax, channel_windows, ...
        ones(max(1, nshades), 1), ...
        repmat([0.75 0.75 0.75], max(1, nshades), 1), ...
        0.16);
    hold(ax, 'on');
end
plot(ax, t, x, 'Color', trace_color, 'LineWidth', 1.0);
xlabel(ax, 'Time (s)');
ylabel(ax, y_label);
title(ax, sprintf('%s (%s) | ROI %d | peak %.2f Hz', ...
    title_text, stage_label, roi_idx, channel_result.roi.peak_frequency_hz(roi_idx)));
grid(ax, 'on');
add_axis_scalebar(ax, t, x, trace_color, scale_suffix);
end

function [y_label, scale_suffix, stage_label] = resolve_time_frequency_trace_label(channel_result)
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
[y_label, scale_suffix, stage_label] = describe_trace_stage_for_display(stage_name);
end

function [y_label, scale_suffix, stage_label] = describe_trace_stage_for_display(stage_name)
stage_name = lower(string(stage_name));
switch stage_name
    case "snr"
        y_label = 'SNR (signed)';
        scale_suffix = 'SNR';
        stage_label = 'SNR';
    case "sensitivity"
        y_label = 'Sensitivity (signed)';
        scale_suffix = 'Sensitivity';
        stage_label = 'Sensitivity';
    case "bleach_removed"
        y_label = 'Bleach-Removed Signal';
        scale_suffix = 'a.u.';
        stage_label = 'bleach removed';
    case "bg_removed"
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

function apply_trace_axis_limits(ax, x)
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

function [peak_frequency_hz, peak_amplitude] = find_band_peak_from_amplitude(frequency, amplitude, display_band_hz)
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
[peak_value, max_idx] = max(band_amplitude);
peak_frequency_hz = band_frequency(max_idx);
peak_amplitude = peak_value;
end

function plot_flash_fft_comparison(channel_result, channel_windows, ~, trace_color, title_text)
[baseline_trace, stim_trace, t_local, supported] = extract_flash_average_segments(channel_result, channel_windows);
if ~supported
    text(0.5, 0.5, 'Flash pre/post FFT unavailable', 'Units', 'normalized', 'HorizontalAlignment', 'center');
    axis off;
    return;
end

[f_base, a_base] = compute_trace_fft_spectrum(baseline_trace, channel_result.frame_rate);
[f_stim, a_stim] = compute_trace_fft_spectrum(stim_trace, channel_result.frame_rate);
base_mask = f_base >= 0 & f_base <= 10;
stim_mask = f_stim >= 0 & f_stim <= 10;
plot(f_base(base_mask), a_base(base_mask), 'Color', [0.55 0.55 0.55], 'LineWidth', 1.5);
hold on;
plot(f_stim(stim_mask), a_stim(stim_mask), 'Color', trace_color, 'LineWidth', 2);
[peak_base_hz, peak_base_amp] = find_band_peak_from_amplitude(f_base, a_base, [0, 10]);
[peak_stim_hz, peak_stim_amp] = find_band_peak_from_amplitude(f_stim, a_stim, [0, 10]);
if isfinite(peak_base_hz)
    scatter(peak_base_hz, peak_base_amp, 36, [0.35 0.35 0.35], 'filled');
end
if isfinite(peak_stim_hz)
    scatter(peak_stim_hz, peak_stim_amp, 36, trace_color, 'filled');
end
xlabel('Frequency (Hz)');
ylabel('Single-Sided Amplitude');
title(sprintf('%s | baseline vs response', title_text));
legend({'Baseline', 'Response'}, 'Location', 'best');
xlim([0, 10]);
grid on;
end

function [baseline_trace, stim_trace, t_local, supported] = extract_flash_average_segments(channel_result, channel_windows)
baseline_trace = [];
stim_trace = [];
t_local = [];
supported = isstruct(channel_windows) ...
    && isfield(channel_windows, 'baseline_frames') ...
    && isfield(channel_windows, 'stim_frames') ...
    && ~isempty(channel_windows.baseline_frames) ...
    && ~isempty(channel_windows.stim_frames);
if ~supported
    return;
end

trace = channel_result.representative.trace(:);
baseline_trials = extract_aligned_trials(trace, channel_windows.baseline_frames);
stim_trials = extract_aligned_trials(trace, channel_windows.stim_frames);
if isempty(baseline_trials) || isempty(stim_trials)
    supported = false;
    return;
end
baseline_trace = mean(baseline_trials, 2, 'omitnan');
stim_trace = mean(stim_trials, 2, 'omitnan');
t_local = (0:numel(baseline_trace) - 1)' / channel_result.frame_rate;
end

function trials = extract_aligned_trials(trace, frame_ranges)
trace = double(trace(:));
frame_ranges = round(frame_ranges);
ntrials = size(frame_ranges, 1);
trial_lengths = frame_ranges(:, 2) - frame_ranges(:, 1) + 1;
valid_lengths = trial_lengths(isfinite(trial_lengths) & trial_lengths > 1);
if isempty(valid_lengths)
    trials = [];
    return;
end
target_length = min(valid_lengths);
trials = NaN(target_length, ntrials);
count = 0;
for idx = 1:ntrials
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

function nshades = count_channel_shading_ranges(channel_windows)
nshades = 0;
if ~isstruct(channel_windows)
    return;
end
if isfield(channel_windows, 'shading_time_ranges') && ~isempty(channel_windows.shading_time_ranges)
    nshades = size(channel_windows.shading_time_ranges, 1);
elseif isfield(channel_windows, 'block_time_ranges') && ~isempty(channel_windows.block_time_ranges)
    nshades = size(channel_windows.block_time_ranges, 1);
elseif isfield(channel_windows, 'stim_time_ranges') && ~isempty(channel_windows.stim_time_ranges)
    nshades = size(channel_windows.stim_time_ranges, 1);
end
end

function plot_wavelet_scalogram(ax, channel_result, channel_windows, title_text, frequency_band_hz)
if nargin < 5 || isempty(frequency_band_hz)
    frequency_band_hz = [];
end
spec = channel_result.representative.wavelet;
if isempty(spec.frequency) || isempty(spec.power)
    text(ax, 0.5, 0.5, 'Wavelet map unavailable', 'Units', 'normalized', 'HorizontalAlignment', 'center');
    axis(ax, 'off');
    return;
end

frequency = spec.frequency;
power = spec.power;
if ~isempty(frequency_band_hz)
    band_mask = frequency >= frequency_band_hz(1) & frequency <= frequency_band_hz(2);
    if any(band_mask)
        frequency = frequency(band_mask);
        power = power(band_mask, :);
    end
end

% CWT frequencies are not linearly spaced. Using imagesc would stretch the
% rows onto a linear y-axis and can make the same power ridge appear at the
% wrong frequency after zooming. Draw the scalogram on the true frequency
% coordinates instead.
[frequency, sort_idx] = sort(frequency(:), 'ascend');
power = power(sort_idx, :);
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
if isempty(frequency_band_hz)
    title(ax, sprintf('%s | ROI %d', title_text, channel_result.representative.roi_index));
else
    title(ax, sprintf('%s | ROI %d | %.0f-%.0f Hz', ...
        title_text, channel_result.representative.roi_index, frequency_band_hz(1), frequency_band_hz(2)));
end
cb = colorbar(ax);
cb.Label.String = 'Wavelet Power';
hold(ax, 'on');
if isfinite(channel_result.representative.peak_frequency_hz)
    yline(ax, channel_result.representative.peak_frequency_hz, 'w--', 'LineWidth', 1.0);
end
if ~isempty(frequency_band_hz)
    ymin = max(min(frequency(frequency > 0)), max(frequency_band_hz(1), 0.5));
    ymax = min(max(frequency), frequency_band_hz(2));
    ylim(ax, [ymin, ymax]);
    yticks(ax, [1 3 10 30]);
else
    ymin = max(min(frequency(frequency > 0)), 0.5);
    ymax = max(frequency);
    ylim(ax, [ymin, ymax]);
    yticks(ax, [1 3 10 30 80]);
end
overlay_stim_boundaries(ax, channel_windows);
end

function plot_roi_wavelet_scalogram(ax, channel_result, roi_idx, channel_windows, title_text)
trace = double(channel_result.trace_matrix(:, roi_idx));
params = channel_result.analysis_parameters;
spec = compute_wavelet_scalogram(trace, channel_result.frame_rate, channel_result.time, params);
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
ymin = max(min(frequency(frequency > 0)), 0.5);
ymax = max(frequency);
ylim(ax, [ymin, ymax]);
yticks(ax, [1 3 10 30 80]);
overlay_stim_boundaries(ax, channel_windows);
end

function overlay_stim_boundaries(ax, channel_windows)
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

function channel_windows = resolve_channel_stim_windows(stim_windows, channel_name)
channel_windows = struct();
if isstruct(stim_windows) && isfield(stim_windows, channel_name)
    channel_windows = stim_windows.(channel_name);
end
end

function cleanup_legacy_oscillation_outputs(save_path)
legacy_files = { ...
    '8_oscillation_waveform_summary.fig', ...
    '8_oscillation_waveform_summary.png', ...
    '8_oscillation_spectral_summary.fig', ...
    '8_oscillation_spectral_summary.png', ...
    '8_oscillation_dual_summary.fig', ...
    '8_oscillation_dual_summary.png', ...
    '8_oscillation_condition_summary.fig', ...
    '8_oscillation_condition_summary.png', ...
    '8_oscillation_results.mat', ...
    '8_fourier_zoom_summary.fig', ...
    '8_fourier_zoom_summary.png', ...
    '8_wavelet_flash_comparison.fig', ...
    '8_wavelet_flash_comparison.png'};
for i = 1:numel(legacy_files)
    legacy_path = fullfile(save_path, legacy_files{i});
    if isfile(legacy_path)
        delete(legacy_path);
    end
end
end

function [trace_data, stage_name, t_axis, frame_rate] = resolve_time_frequency_channel_input(results, preferred_stages, t_axis_default, frame_rate_default)
[trace_data, stage_name] = resolve_preferred_trace_stage(results, preferred_stages);
t_axis = [];
frame_rate = frame_rate_default;
if isstruct(results) && isfield(results, 'trace_results') && isfield(results.trace_results, stage_name)
    stage_struct = results.trace_results.(stage_name);
    if isfield(stage_struct, 'time')
        t_axis = double(stage_struct.time(:));
    end
    if isfield(stage_struct, 'frame_rate')
        frame_rate = double(stage_struct.frame_rate);
    end
end

if isempty(t_axis)
    if ~isempty(t_axis_default)
        t_axis = double(t_axis_default(:));
    else
        t_axis = (0:size(trace_data, 1) - 1)' / frame_rate;
    end
end
end

function frame_rate = resolve_saved_channel_frame_rate(dual_info, role_name)
frame_rate = NaN;
if isfield(dual_info, 'camera_cfg')
    camera_cfg = dual_info.camera_cfg;
    for idx = 1:numel(camera_cfg)
        if isfield(camera_cfg(idx), 'role') && strcmpi(string(camera_cfg(idx).role), role_name)
            if isfield(camera_cfg(idx), 'frame_rate') && ~isempty(camera_cfg(idx).frame_rate)
                frame_rate = double(camera_cfg(idx).frame_rate);
                return;
            end
        end
    end
end
if ~isfinite(frame_rate)
    error('Cannot resolve frame_rate for %s from dual_info.', role_name);
end
end

function frame_count = resolve_saved_channel_frame_count(movie_info)
frame_count = NaN;
if isfield(movie_info, 'frame_count') && ~isempty(movie_info.frame_count)
    frame_count = double(movie_info.frame_count);
elseif isfield(movie_info, 'movie_size') && numel(movie_info.movie_size) >= 3
    frame_count = double(movie_info.movie_size(3));
end
if ~isfinite(frame_count)
    error('Cannot resolve frame_count from saved movie_info.');
end
end

function [ncols, nrows] = resolve_saved_analysis_frame_size(voltage_movie_info, calcium_movie_info)
ncols = NaN;
nrows = NaN;
if isfield(voltage_movie_info, 'analysis_frame_size') && numel(voltage_movie_info.analysis_frame_size) >= 2
    ncols = double(voltage_movie_info.analysis_frame_size(1));
    nrows = double(voltage_movie_info.analysis_frame_size(2));
elseif isfield(calcium_movie_info, 'analysis_frame_size') && numel(calcium_movie_info.analysis_frame_size) >= 2
    ncols = double(calcium_movie_info.analysis_frame_size(1));
    nrows = double(calcium_movie_info.analysis_frame_size(2));
end
if ~isfinite(ncols) || ~isfinite(nrows)
    error('Cannot resolve saved analysis frame size from movie_info.');
end
end

function value = load_required_struct(mat_path, variable_name)
if ~isfile(mat_path)
    error('Required saved result file is missing: %s', mat_path);
end
tmp = load(mat_path, variable_name);
if ~isfield(tmp, variable_name)
    error('Variable %s is missing from %s', variable_name, mat_path);
end
value = tmp.(variable_name);
end

function save_explicit_dual_results_summary( ...
    save_path, ...
    dual_info, dual_info_path, ...
    voltage_results, voltage_results_path, ...
    calcium_results, calcium_results_path, ...
    dual_results, dual_results_path, ...
    stim_results, stim_results_path, ...
    time_frequency_results, time_frequency_results_path)
explicit_dual_results = struct( ...
    'dual_info', dual_info, ...
    'result_files', struct( ...
        'dual_info', dual_info_path, ...
        'voltage_results', voltage_results_path, ...
        'calcium_results', calcium_results_path, ...
        'dual_results', dual_results_path, ...
        'stim_results', stim_results_path, ...
        'time_frequency_results', time_frequency_results_path), ...
    'available_voltage_trace_stages', string(fieldnames(voltage_results.trace_results)), ...
    'available_calcium_trace_stages', string(fieldnames(calcium_results.trace_results)), ...
    'dual_result_fields', string(fieldnames(dual_results)), ...
    'stim_result_fields', string(fieldnames(stim_results)), ...
    'time_frequency_result_fields', string(fieldnames(time_frequency_results)));
save(fullfile(save_path, '-1_explicit_dual_results.mat'), 'explicit_dual_results', '-v7.3');
end

function print_section(section_name)
fprintf('\n============================================================\n');
fprintf('[Section] %s\n', section_name);
fprintf('============================================================\n');
end

function formula_text = describe_bleach_formula(bleach_mode)
switch lower(bleach_mode)
    case 'highpass'
        formula_text = 'baseline = highpass estimate; bleach_removed = input - baseline';
    case 'linear'
        formula_text = 'bleach_removed = detrend(input, 1); baseline = input - bleach_removed';
    case 'exp2'
        formula_text = 'baseline = fit_exp2(input); bleach_removed = input - baseline';
    otherwise
        formula_text = 'unknown bleaching formula';
end
end

function [voltage_results, calcium_results] = load_channel_results(voltage_results_path, calcium_results_path, voltage_results, calcium_results)
% Refresh channel result containers from disk so later sections can run
% independently of the current workspace state.
if (isempty(voltage_results) || ~isstruct(voltage_results) || ~isfield(voltage_results, 'trace_results')) && isfile(voltage_results_path)
    tmp = load(voltage_results_path, 'voltage_results');
    if isfield(tmp, 'voltage_results')
        voltage_results = tmp.voltage_results;
    end
end

if (isempty(calcium_results) || ~isstruct(calcium_results) || ~isfield(calcium_results, 'trace_results')) && isfile(calcium_results_path)
    tmp = load(calcium_results_path, 'calcium_results');
    if isfield(tmp, 'calcium_results')
        calcium_results = tmp.calcium_results;
    end
end
end

function [dual_results, rois, traces_voltage_raw, traces_calcium_raw, roi_results_file, nrois] = load_dual_roi_context(dual_results_path, dual_results)
% Reload the shared ROI context that ties both channels together. This is
% the minimum common state needed by most later trace-processing sections.
rois = struct();
traces_voltage_raw = [];
traces_calcium_raw = [];
roi_results_file = '';
nrois = [];

if (isempty(dual_results) || ~isstruct(dual_results) || ~isfield(dual_results, 'registration')) && isfile(dual_results_path)
    tmp = load(dual_results_path, 'dual_results');
    if isfield(tmp, 'dual_results')
        dual_results = tmp.dual_results;
    end
end

if isstruct(dual_results) && isfield(dual_results, 'registration') ...
        && isfield(dual_results.registration, 'info') ...
        && isfield(dual_results.registration.info, 'roi_file')
    roi_results_file = dual_results.registration.info.roi_file;
end

if ~isempty(roi_results_file) && isfile(roi_results_file)
    tmp = load(roi_results_file);
    if isfield(tmp, 'rois')
        rois = tmp.rois;
    end
    if isfield(tmp, 'traces_voltage_raw')
        traces_voltage_raw = tmp.traces_voltage_raw;
    end
    if isfield(tmp, 'traces_calcium_raw')
        traces_calcium_raw = tmp.traces_calcium_raw;
    end
    if isfield(tmp, 'nrois')
        nrois = tmp.nrois;
    elseif ~isempty(traces_voltage_raw)
        nrois = size(traces_voltage_raw, 2);
    end
end
end

function bootstrap_volpy_voltage_reanalysis(cycle_path, save_path, source_results_path, auto_run, force_rerun, use_existing, flip_signal, voltage_polarity, calcium_polarity, calcium_smoothing_window)
% Build a standalone voltage-reanalysis bundle that reuses the existing
% dual ROI/calcium/stim context while replacing only the voltage trace
% stages with VolPy-derived outputs saved under a separate folder.
source_results_path = resolve_volpy_source_results_path(cycle_path, source_results_path);
[source_dual_info, source_voltage_results, source_calcium_results, source_dual_results, source_stim_results, source_stim_context, source_stim_windows] = ...
    load_saved_analysis_only_context(source_results_path);
source_dual_results_path = fullfile(source_results_path, 'dual_results.mat');
[source_dual_results, rois, ~, ~, roi_results_file] = load_dual_roi_context(source_dual_results_path, source_dual_results);
if isempty(fieldnames(rois)) || ~isfield(rois, 'bwmask') || isempty(rois.bwmask)
    error('VolPy voltage re-analysis requires an existing dual ROI mask in the source results.');
end

voltage_movie_path = resolve_channel_movie_path_from_dual_info(source_dual_info, "voltage");
if ~isfile(voltage_movie_path)
    error('Voltage movie for VolPy re-analysis was not found: %s', voltage_movie_path);
end

voltage_frame_rate = double(source_voltage_results.movie_info.frame_rate);
volpy_output_dir = fullfile(save_path, '0_volpy_backend');
if ~isfolder(volpy_output_dir)
    mkdir(volpy_output_dir);
end
roi_mask_file = fullfile(volpy_output_dir, 'volpy_roi_mask.mat');
mask = uint16(rois.bwmask);
save(roi_mask_file, 'mask');

result_mat_path = fullfile(volpy_output_dir, 'volpy_results.mat');
volpy_backend_info = run_volpy_backend_pipeline_dual( ...
    voltage_movie_path, voltage_frame_rate, flip_signal, ...
    auto_run, force_rerun, use_existing, ...
    volpy_output_dir, result_mat_path, roi_mask_file);
volpy_data = load_volpy_backend_results_dual(char(volpy_backend_info.result_mat));

voltage_results = source_voltage_results;
voltage_results = preserve_standard_voltage_stage(voltage_results, 'raw');
voltage_results = preserve_standard_voltage_stage(voltage_results, 'bleach_removed');
voltage_results = preserve_standard_voltage_stage(voltage_results, 'baseline');
voltage_results = preserve_standard_voltage_stage(voltage_results, 'noise_reference');
voltage_results = preserve_standard_voltage_stage(voltage_results, 'noise');
voltage_results = preserve_standard_voltage_stage(voltage_results, 'sensitivity');
voltage_results = preserve_standard_voltage_stage(voltage_results, 'snr');

voltage_results.movie_info.motion.applied = true;
voltage_results.movie_info.motion.method = 'VolPy internal MotionCorrect';
voltage_results.movie_info.motion.shift_file = char(volpy_backend_info.motion_corrected_file);
voltage_results.movie_info.motion.parameter_file = char(volpy_backend_info.result_mat);
voltage_results.movie_info.updated_at = datetime("now");

voltage_results = store_trace_stage( ...
    voltage_results, 'raw', volpy_data.t, {}, roi_results_file, ...
    voltage_results.movie_info, 'volpy_trace', struct('source', char(volpy_backend_info.result_mat)));
voltage_results = store_trace_stage( ...
    voltage_results, 'bleach_removed', volpy_data.t, {'raw'}, roi_results_file, ...
    voltage_results.movie_info, 'volpy_import_compat', struct('source_stage', 'volpy_t'));
voltage_results = store_trace_stage( ...
    voltage_results, 'baseline', volpy_data.f0, {'volpy_t'}, roi_results_file, ...
    voltage_results.movie_info, 'volpy_F0', struct());
voltage_results = store_trace_stage( ...
    voltage_results, 'noise_reference', volpy_data.t_rec, {'volpy_t'}, roi_results_file, ...
    voltage_results.movie_info, 'volpy_reconstructed_spike_trace', struct());
voltage_results = store_trace_stage( ...
    voltage_results, 'noise', volpy_data.noise, {'volpy_t', 'noise_reference'}, roi_results_file, ...
    voltage_results.movie_info, 'volpy_residual', struct('expression', 'volpy_t - noise_reference'));
voltage_results = store_trace_stage( ...
    voltage_results, 'sensitivity', volpy_data.dff, {'volpy_t', 'baseline'}, roi_results_file, ...
    voltage_results.movie_info, 'volpy_dff_import', struct());
voltage_results = store_trace_stage( ...
    voltage_results, 'snr', volpy_data.snr_trace, {'volpy_t', 'noise'}, roi_results_file, ...
    voltage_results.movie_info, 'volpy_trace_divided_by_residual_std', struct());
voltage_results = store_trace_stage( ...
    voltage_results, 'volpy_t', volpy_data.t, {'raw'}, roi_results_file, ...
    voltage_results.movie_info, 'volpy_trace', struct('source', char(volpy_backend_info.result_mat)));
voltage_results = store_trace_stage( ...
    voltage_results, 'volpy_ts', volpy_data.ts, {'volpy_t'}, roi_results_file, ...
    voltage_results.movie_info, 'volpy_matched_filter_trace', struct());
voltage_results = store_trace_stage( ...
    voltage_results, 'volpy_t_rec', volpy_data.t_rec, {'volpy_t'}, roi_results_file, ...
    voltage_results.movie_info, 'volpy_reconstructed_spike_trace', struct());
voltage_results = store_trace_stage( ...
    voltage_results, 'volpy_subthreshold', volpy_data.t_sub, {'volpy_t'}, roi_results_file, ...
    voltage_results.movie_info, 'volpy_subthreshold', struct());
voltage_results = store_trace_stage( ...
    voltage_results, 'volpy_dff', volpy_data.dff, {'volpy_t', 'baseline'}, roi_results_file, ...
    voltage_results.movie_info, 'volpy_dff', struct());
voltage_results = store_trace_stage( ...
    voltage_results, 'volpy_noise', volpy_data.noise, {'volpy_t', 'volpy_t_rec'}, roi_results_file, ...
    voltage_results.movie_info, 'volpy_residual', struct('expression', 'volpy_t - volpy_t_rec'));

calcium_results = source_calcium_results;
dual_results = source_dual_results;
stim_results = source_stim_results;
dual_info = source_dual_info;
dual_info.analysis_name = 'Dual_analysis3';
dual_info.analysis_backend = 'volpy_voltage_reanalysis';
dual_info.save_path = save_path;
dual_info.created_at = datetime("now");
dual_info.volpy_source_results_path = source_results_path;
dual_info.volpy_backend = struct( ...
    'output_dir', string(volpy_backend_info.output_dir), ...
    'result_mat', string(volpy_backend_info.result_mat), ...
    'motion_corrected_file', string(volpy_backend_info.motion_corrected_file), ...
    'roi_mask_file', string(roi_mask_file), ...
    'flip_signal', logical(flip_signal), ...
    'voltage_movie_path', string(voltage_movie_path), ...
    'calcium_smoothing_window', calcium_smoothing_window, ...
    'voltage_polarity', voltage_polarity, ...
    'calcium_polarity', calcium_polarity);
if ~isstruct(dual_results)
    dual_results = struct();
end
dual_results.volpy_voltage_reanalysis = struct( ...
    'data', struct(), ...
    'info', struct( ...
        'source_results_path', string(source_results_path), ...
        'result_mat', string(volpy_backend_info.result_mat), ...
        'motion_corrected_file', string(volpy_backend_info.motion_corrected_file), ...
        'roi_mask_file', string(roi_mask_file), ...
        'created_at', datetime("now")));
stim_results.info = rmfield_if_exists(source_stim_context, {'logs', 'method_manifest'});
stim_results.windows = source_stim_windows;

save(fullfile(save_path, 'dual_info.mat'), 'dual_info');
save(fullfile(save_path, 'voltage_results.mat'), 'voltage_results', '-v7.3');
save(fullfile(save_path, 'calcium_results.mat'), 'calcium_results', '-v7.3');
save(fullfile(save_path, 'dual_results.mat'), 'dual_results', '-v7.3');
save(fullfile(save_path, 'stim_results.mat'), 'stim_results', '-v7.3');
fprintf('VolPy voltage re-analysis bundle prepared in: %s\n', save_path);
fprintf('VolPy backend result: %s\n', volpy_backend_info.result_mat);
end

function source_results_path = resolve_volpy_source_results_path(cycle_path, source_results_path)
source_results_path = string(source_results_path);
if strlength(source_results_path) > 0
    if ~isfolder(source_results_path)
        error('Provided volpy_source_results_path does not exist: %s', source_results_path);
    end
    return;
end
explicit_result = find_latest_standard_explicit_result(cycle_path);
if strlength(explicit_result) == 0
    error('No existing standard Dual_analysis3 result was found under this cycle for VolPy re-analysis.');
end
source_results_path = string(fileparts(explicit_result));
end

function explicit_result = find_latest_standard_explicit_result(cycle_path)
explicit_result = "";
listing = dir(fullfile(cycle_path, 'Dual_analysis3', '**', '-1_explicit_dual_results.mat'));
if isempty(listing)
    return;
end
[~, newest_idx] = max([listing.datenum]);
explicit_result = string(fullfile(listing(newest_idx).folder, listing(newest_idx).name));
end

function movie_path = resolve_channel_movie_path_from_dual_info(dual_info, role_name)
movie_path = "";
if ~isstruct(dual_info) || ~isfield(dual_info, 'camera_source')
    return;
end
camera_source = dual_info.camera_source;
camera_cfg = dual_info.camera_cfg;
for idx = 1:min(numel(camera_source), numel(camera_cfg))
    if strcmpi(string(camera_cfg(idx).role), string(role_name))
        movie_path = string(camera_source(idx).path);
        movie_path = resolve_movie_file_from_camera_source_path(movie_path);
        return;
    end
end
end

function movie_file = resolve_movie_file_from_camera_source_path(camera_source_path)
camera_source_path = string(camera_source_path);
movie_file = camera_source_path;
if strlength(movie_file) == 0
    return;
end
if isfile(movie_file)
    return;
end
if ~isfolder(movie_file)
    movie_file = "";
    return;
end
listing = dir(fullfile(char(movie_file), '*.tif'));
if isempty(listing)
    listing = dir(fullfile(char(movie_file), '*.tiff'));
end
if isempty(listing)
    movie_file = "";
    return;
end
if numel(listing) == 1
    movie_file = string(fullfile(listing(1).folder, listing(1).name));
    return;
end
[~, newest_idx] = max([listing.datenum]);
movie_file = string(fullfile(listing(newest_idx).folder, listing(newest_idx).name));
end

function info = run_volpy_backend_pipeline_dual(input_movie_path, frame_rate, flip_signal, auto_run, force_rerun, use_existing, output_dir, result_mat_path, roi_mask_file)
python_exe = 'C:\Users\DELL\anaconda3\envs\caiman\python.exe';
script_path = fullfile(fileparts(mfilename('fullpath')), 'python_seg', 'run_volpy_backend.py');

if ~isfile(python_exe)
    error('VolPy backend python interpreter not found: %s', python_exe);
end
if ~isfile(script_path)
    error('VolPy backend script not found: %s', script_path);
end
if ~isfolder(output_dir)
    mkdir(output_dir);
end

result_exists = isfile(result_mat_path);
should_run = force_rerun || ~(result_exists && use_existing);
cmdout = "existing result reused";
if should_run
    if ~auto_run && ~force_rerun
        error('VolPy backend result not found and auto-run is disabled: %s', result_mat_path);
    end
    flip_text = ternary(flip_signal, 'true', 'false');
    cmd = sprintf('"%s" "%s" "%s" --output-dir "%s" --frame-rate %.12g --flip-signal %s --roi-mask "%s"', ...
        python_exe, script_path, input_movie_path, output_dir, frame_rate, flip_text, roi_mask_file);
    [status, raw_cmdout] = system(cmd);
    cmdout = string(raw_cmdout);
    if status ~= 0
        error('VolPy backend command failed:\n%s', raw_cmdout);
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
    'command_output', string(cmdout));
end

function data = load_volpy_backend_results_dual(result_mat_path)
s = load(result_mat_path);
data = struct();
data.mask = uint16(s.mask);
data.t = cell_columns_to_matrix_dual(s.t);
data.ts = cell_columns_to_matrix_dual(s.ts);
data.t_rec = cell_columns_to_matrix_dual(s.t_rec);
data.t_sub = cell_columns_to_matrix_dual(s.t_sub);
if isfield(s, 'F0')
    data.f0 = cell_columns_to_matrix_dual(s.F0);
else
    data.f0 = ones(size(data.t));
end
if isfield(s, 'dFF')
    data.dff = cell_columns_to_matrix_dual(s.dFF);
else
    baseline_safe = data.f0;
    baseline_safe(abs(baseline_safe) < eps) = 1;
    data.dff = data.t ./ baseline_safe;
end

nframes = size(data.t, 1);
nrois = size(data.t, 2);
data.noise = data.t - data.t_rec;
data.snr_trace = zeros(size(data.t));
for i = 1:nrois
    noise_std = std(data.noise(:, i), 0, 1);
    if ~isfinite(noise_std) || noise_std <= eps
        noise_std = 1;
    end
    data.snr_trace(:, i) = data.t(:, i) ./ noise_std;
end

data.spikes = cell(1, nrois);
data.peak_amplitude = cell(1, nrois);
data.peaks_polarity = cell(1, nrois);
for i = 1:nrois
    idx = round(double(s.spikes{i}(:)));
    idx = idx(idx >= 1 & idx <= nframes);
    data.spikes{i} = idx;
    data.peak_amplitude{i} = data.t(idx, i);
    polarity = 1;
    if isfield(s, 'templates') && numel(s.templates) >= i && ~isempty(s.templates{i})
        template_i = double(s.templates{i}(:));
        template_i = template_i(~isnan(template_i));
        if ~isempty(template_i)
            polarity = sign(sum(template_i));
        end
    end
    if ~isfinite(polarity) || polarity == 0
        polarity = 1;
    end
    data.peaks_polarity{i} = polarity;
end

if isfield(s, 'snr')
    data.snr_scalar = double(s.snr(:));
else
    data.snr_scalar = zeros(nrois, 1);
end
if isfield(s, 'locality')
    data.locality = logical(s.locality(:));
else
    data.locality = true(nrois, 1);
end
if isfield(s, 'low_spikes')
    data.low_spikes = logical(s.low_spikes(:));
else
    data.low_spikes = false(nrois, 1);
end
end

function matrix = cell_columns_to_matrix_dual(cell_values)
if ~iscell(cell_values)
    matrix = double(cell_values);
    return;
end
if isempty(cell_values)
    matrix = [];
    return;
end
column_cells = cellfun(@(v) double(v(:)), cell_values, 'UniformOutput', false);
matrix = cell2mat(column_cells);
end

function results = preserve_standard_voltage_stage(results, stage_name)
backup_name = "standard_" + string(stage_name);
if has_trace_stage(results, stage_name) && ~has_trace_stage(results, char(backup_name))
    results.trace_results.(char(backup_name)) = results.trace_results.(stage_name);
end
end

function nframes = infer_frame_count_from_channel_results(results)
nframes = [];
if isstruct(results) && isfield(results, 'movie_info') && isfield(results.movie_info, 'frame_count') ...
        && ~isempty(results.movie_info.frame_count)
    nframes = double(results.movie_info.frame_count);
    return;
end
preferred_fields = fieldnames(results.trace_results);
for idx = 1:numel(preferred_fields)
    stage = preferred_fields{idx};
    if isfield(results.trace_results.(stage), 'data') && ~isempty(results.trace_results.(stage).data)
        nframes = size(results.trace_results.(stage).data, 1);
        return;
    end
end
error('Could not infer frame count from saved channel results.');
end

function run_dual_post_trace_sections( ...
    save_path, ...
    dual_info, dual_info_path, ...
    voltage_results, voltage_results_path, ...
    calcium_results, calcium_results_path, ...
    dual_results, dual_results_path, ...
    stim_results, stim_results_path, ...
    stim_context, stim_windows, ...
    t_voltage, t_calcium, ...
    freq_voltage, freq_calcium, ...
    nframes_voltage, nframes_calcium, ...
    voltage_polarity, calcium_polarity, calcium_smoothing_window, reuse_saved_stim_windows)
print_section('Channel Summary Plots');
fprintf('Saving channel summary plots...\n');
[voltage_results, calcium_results] = load_channel_results( ...
    voltage_results_path, calcium_results_path, voltage_results, calcium_results);
traces_voltage_raw = fetch_trace_stage(voltage_results, 'raw');
has_calcium_raw = has_trace_stage(calcium_results, 'raw_smoothed') || has_trace_stage(calcium_results, 'raw');
has_voltage_sensitivity = has_trace_stage(voltage_results, 'sensitivity');
has_voltage_snr = has_trace_stage(voltage_results, 'snr');
has_calcium_sensitivity = has_trace_stage(calcium_results, 'sensitivity_smoothed') || has_trace_stage(calcium_results, 'sensitivity');
has_calcium_snr = has_trace_stage(calcium_results, 'snr_smoothed') || has_trace_stage(calcium_results, 'snr');
if has_voltage_sensitivity
    voltage_sensitivity = fetch_trace_stage(voltage_results, 'sensitivity');
else
    voltage_sensitivity = [];
end
if has_voltage_snr
    voltage_snr = fetch_trace_stage(voltage_results, 'snr');
else
    voltage_snr = [];
end
if has_calcium_raw
    [traces_calcium_raw, calcium_raw_stage] = resolve_preferred_trace_stage(calcium_results, {'raw_smoothed', 'raw'});
else
    traces_calcium_raw = [];
    calcium_raw_stage = '';
end
if has_calcium_sensitivity
    [calcium_sensitivity, calcium_sensitivity_stage] = resolve_preferred_trace_stage(calcium_results, {'sensitivity_smoothed', 'sensitivity'});
else
    calcium_sensitivity = [];
    calcium_sensitivity_stage = '';
end
if has_calcium_snr
    [calcium_snr, calcium_snr_stage] = resolve_preferred_trace_stage(calcium_results, {'snr_smoothed', 'snr'});
else
    calcium_snr = [];
    calcium_snr_stage = '';
end
summary_fig = figure('Color', 'w');
subplot(2, 3, 1); hold on; title('Voltage Raw'); plot_offset_stage(traces_voltage_raw, t_voltage);
subplot(2, 3, 2); hold on; title('Voltage Sensitivity'); plot_optional_stage(voltage_polarity * voltage_sensitivity, t_voltage, has_voltage_sensitivity, 'Sensitivity stage unavailable');
subplot(2, 3, 3); hold on; title('Voltage SNR'); plot_optional_stage(voltage_polarity * voltage_snr, t_voltage, has_voltage_snr, 'SNR stage unavailable');
subplot(2, 3, 4); hold on; title(sprintf('Calcium Raw (%s, window=%d)', strrep(char(calcium_raw_stage), '_', '\_'), calcium_smoothing_window)); plot_optional_stage(traces_calcium_raw, t_calcium, has_calcium_raw, 'Raw stage unavailable');
subplot(2, 3, 5); hold on; title(sprintf('Calcium Sensitivity (%s, window=%d)', strrep(char(calcium_sensitivity_stage), '_', '\_'), calcium_smoothing_window)); plot_optional_stage(calcium_polarity * calcium_sensitivity, t_calcium, has_calcium_sensitivity, 'Sensitivity stage unavailable');
subplot(2, 3, 6); hold on; title(sprintf('Calcium SNR (%s, window=%d)', strrep(char(calcium_snr_stage), '_', '\_'), calcium_smoothing_window)); plot_optional_stage(calcium_polarity * calcium_snr, t_calcium, has_calcium_snr, 'SNR stage unavailable');
save_figure_bundle(summary_fig, fullfile(save_path, '4_dual_trace_summary.fig'), fullfile(save_path, '4_dual_trace_summary.png'));
close(summary_fig);

% ROI-aligned population view for analysis-only / post-trace reruns.
% This is intentionally kept inline, matching the main-script section, so
% the figure can be debugged by setting breakpoints directly in the rerun
% path.
%
% Display contract:
%   - one y-row represents one paired ROI;
%   - calcium sensitivity is the background heatmap;
%   - voltage sensitivity is overlaid as a red trace in the same ROI row;
%   - voltage keeps its saved high-rate time axis and is not resampled;
%   - the x-limits follow the calcium heatmap time axis.
%
% Scaling contract:
%   - find the voltage ROI with the largest displayed max-min range;
%   - use that range as the global scale;
%   - center every ROI by its own midpoint before plotting;
%   - the reference ROI fills one full row, and other ROIs shrink relative
%     to that same scale.
%
% Calcium heatmap baseline/color contract:
%   - prefer sensitivity_smoothed, falling back to sensitivity;
%   - subtract each ROI's 1st percentile as displayed sensitivity zero;
%   - clip values below that percentile-zero baseline to 0 for display;
%   - use [0 max] on the percentile-zeroed display matrix.
print_section('ROI Calcium Heatmap With Voltage Trace');
fprintf('Saving ROI-aligned calcium heatmap with overlaid voltage sensitivity traces...\n');
[voltage_results, calcium_results] = load_channel_results( ...
    voltage_results_path, calcium_results_path, voltage_results, calcium_results);

% Initialize a result record before validation so this section leaves a
% clear status in dual_results even when it skips on old or incomplete
% analysis-only result folders.
roi_heatmap_trace_result = struct( ...
    'status', "skipped", ...
    'reason', "", ...
    'fig_file', "", ...
    'png_file', "", ...
    'mat_file', "", ...
    'created_at', datetime("now"));

% Required inputs are one voltage sensitivity stage and one calcium
% sensitivity stage. Calcium prefers sensitivity_smoothed because that is
% the final calcium display stage created after metric computation.
has_voltage_sensitivity = has_trace_stage(voltage_results, 'sensitivity');
has_calcium_sensitivity = has_trace_stage(calcium_results, 'sensitivity_smoothed') ...
    || has_trace_stage(calcium_results, 'sensitivity');
if ~has_voltage_sensitivity || ~has_calcium_sensitivity
    roi_heatmap_trace_result.reason = sprintf('Missing required sensitivity stage: voltage=%d calcium=%d.', ...
        has_voltage_sensitivity, has_calcium_sensitivity);
    fprintf('Skipping ROI calcium heatmap/voltage trace plot: %s\n', roi_heatmap_trace_result.reason);
else
    % Fetch trace stages and apply display polarity up front. Everything
    % below uses these display matrices, including the saved diagnostic MAT.
    voltage_heatmap_stage = 'sensitivity';
    voltage_heatmap_display = double(voltage_polarity) * double(fetch_trace_stage(voltage_results, voltage_heatmap_stage));
    [calcium_heatmap_display, calcium_heatmap_stage] = resolve_preferred_trace_stage(calcium_results, {'sensitivity_smoothed', 'sensitivity'});
    calcium_heatmap_display = double(calcium_polarity) * double(calcium_heatmap_display);

    % Repair time-axis length mismatches defensively without resampling the
    % trace data. Too-long axes are truncated; too-short axes are extended
    % with the median frame interval. The rule is stored in time_axis_info.
    t_voltage_heatmap = double(t_voltage(:));
    voltage_time_info = struct('channel', "voltage", 'input_time_points', numel(t_voltage_heatmap), ...
        'trace_frames', size(voltage_heatmap_display, 1), 'rule', "unchanged");
    if numel(t_voltage_heatmap) > size(voltage_heatmap_display, 1)
        t_voltage_heatmap = t_voltage_heatmap(1:size(voltage_heatmap_display, 1));
        voltage_time_info.rule = "truncated time axis to trace frame count";
    elseif numel(t_voltage_heatmap) < size(voltage_heatmap_display, 1)
        if numel(t_voltage_heatmap) >= 2
            dt_heatmap = median(diff(t_voltage_heatmap), 'omitnan');
            if ~isfinite(dt_heatmap) || dt_heatmap <= 0
                dt_heatmap = 1;
            end
            first_t_heatmap = t_voltage_heatmap(1);
        else
            dt_heatmap = 1;
            first_t_heatmap = 0;
        end
        t_voltage_heatmap = first_t_heatmap + (0:size(voltage_heatmap_display, 1)-1)' * dt_heatmap;
        voltage_time_info.rule = "extended time axis using median dt";
    end

    t_calcium_heatmap = double(t_calcium(:));
    calcium_time_info = struct('channel', "calcium", 'input_time_points', numel(t_calcium_heatmap), ...
        'trace_frames', size(calcium_heatmap_display, 1), 'rule', "unchanged");
    if numel(t_calcium_heatmap) > size(calcium_heatmap_display, 1)
        t_calcium_heatmap = t_calcium_heatmap(1:size(calcium_heatmap_display, 1));
        calcium_time_info.rule = "truncated time axis to trace frame count";
    elseif numel(t_calcium_heatmap) < size(calcium_heatmap_display, 1)
        if numel(t_calcium_heatmap) >= 2
            dt_heatmap = median(diff(t_calcium_heatmap), 'omitnan');
            if ~isfinite(dt_heatmap) || dt_heatmap <= 0
                dt_heatmap = 1;
            end
            first_t_heatmap = t_calcium_heatmap(1);
        else
            dt_heatmap = 1;
            first_t_heatmap = 0;
        end
        t_calcium_heatmap = first_t_heatmap + (0:size(calcium_heatmap_display, 1)-1)' * dt_heatmap;
        calcium_time_info.rule = "extended time axis using median dt";
    end

    % Pair ROIs by saved column order. If an older result folder has unequal
    % voltage/calcium ROI counts, plot only the shared prefix and record the
    % original counts in the result metadata.
    nrois_voltage_heatmap = size(voltage_heatmap_display, 2);
    nrois_calcium_heatmap = size(calcium_heatmap_display, 2);
    nrois_heatmap = min(nrois_voltage_heatmap, nrois_calcium_heatmap);
    if nrois_heatmap == 0
        roi_heatmap_trace_result.reason = 'No shared ROI columns are available in the selected sensitivity stages.';
        fprintf('Skipping ROI calcium heatmap/voltage trace plot: %s\n', roi_heatmap_trace_result.reason);
    else
        if nrois_voltage_heatmap ~= nrois_calcium_heatmap
            fprintf('ROI heatmap/trace ROI mismatch: voltage=%d calcium=%d; plotting first %d paired ROIs.\n', ...
                nrois_voltage_heatmap, nrois_calcium_heatmap, nrois_heatmap);
        end
        voltage_heatmap_display = voltage_heatmap_display(:, 1:nrois_heatmap);
        calcium_heatmap_display_raw = calcium_heatmap_display(:, 1:nrois_heatmap);

        % Voltage overlay scaling:
        %   1. compute each ROI's displayed max-min range;
        %   2. use the largest ROI range as the single global scale;
        %   3. center each ROI by its own midpoint before scaling.
        %
        % This makes the reference ROI fill exactly one heatmap row from
        % bottom to top, while all other ROIs are comparable to that same
        % reference amplitude.
        voltage_roi_min = min(voltage_heatmap_display, [], 1, 'omitnan');
        voltage_roi_max = max(voltage_heatmap_display, [], 1, 'omitnan');
        voltage_roi_range = voltage_roi_max - voltage_roi_min;
        [voltage_global_range, voltage_reference_roi] = max(voltage_roi_range);
        if ~isfinite(voltage_global_range) || voltage_global_range <= 0
            voltage_global_range = 1;
            voltage_reference_roi = 1;
        end
        voltage_roi_mid = (voltage_roi_min + voltage_roi_max) / 2;
        voltage_roi_mid(~isfinite(voltage_roi_mid)) = 0;

        % Calcium heatmap baseline rule:
        %   each ROI uses its own 1st percentile as displayed sensitivity 0.
        % Values below that baseline are clipped to 0 for display, while
        % the raw polarity-adjusted calcium matrix is still saved below.
        calcium_heatmap_zero_percentile = 1;
        calcium_heatmap_baseline = prctile(calcium_heatmap_display_raw, calcium_heatmap_zero_percentile, 1);
        calcium_heatmap_baseline(~isfinite(calcium_heatmap_baseline)) = 0;
        calcium_heatmap_display_zeroed = calcium_heatmap_display_raw - calcium_heatmap_baseline;
        calcium_heatmap_display_zeroed(calcium_heatmap_display_zeroed < 0) = 0;

        % After percentile-zeroing, color limits are intentionally [0 max].
        % A flat or empty heatmap falls back to [0 1] so clim remains valid.
        finite_calcium = calcium_heatmap_display_zeroed(isfinite(calcium_heatmap_display_zeroed));
        if isempty(finite_calcium)
            calcium_heatmap_clim = [0, 1];
            calcium_heatmap_clim_rule = "fallback_empty_to_0_1";
        else
            calcium_heatmap_max = max(finite_calcium);
            if calcium_heatmap_max > 0
                calcium_heatmap_clim = [0, calcium_heatmap_max];
                calcium_heatmap_clim_rule = "roi_1st_percentile_zero_to_max";
            else
                calcium_heatmap_clim = [0, 1];
                calcium_heatmap_clim_rule = "fallback_flat_to_0_1";
            end
        end
        if calcium_heatmap_clim(1) == calcium_heatmap_clim(2)
            calcium_heatmap_clim = calcium_heatmap_clim + [-0.5, 0.5];
            calcium_heatmap_clim_rule = calcium_heatmap_clim_rule + "_expanded_flat_range";
        end

        % Inline black-blue-white colormap from the MLX reference. Keeping
        % it inline here makes this section easy to tune interactively.
        gamma_val_calcium = 0.8;
        color_pivot_calcium = 0.1;
        cmap_n = 256;
        base_ice = zeros(cmap_n, 3);
        base_ice(:, 3) = linspace(0, 1, cmap_n);
        cyan_start_node = max(1, floor(cmap_n * 0.1));
        base_ice(cyan_start_node:end, 2) = linspace(0, 1, cmap_n - cyan_start_node + 1);
        white_start_node = max(1, floor(cmap_n * 0.9));
        base_ice(white_start_node:end, 1) = linspace(0, 1, cmap_n - white_start_node + 1);
        x_old = linspace(0, 1, cmap_n);
        x_new = linspace(0, 1, cmap_n) .^ gamma_val_calcium;
        warped = interp1([0, color_pivot_calcium, 1], [0, 0.5, 1], x_new, 'linear', 'extrap');
        warped = min(max(warped, 0), 1);
        calcium_heatmap_cmap = interp1(x_old, base_ice, warped);

        % Diagnostic transfer-curve plot for tuning the calcium heatmap.
        % The dashed curve shows x^gamma and the solid curve shows the
        % final colormap index after the pivot warp. This makes it easier
        % to see why changing gamma or pivot changes background contrast.
        gamma_curve_input = linspace(0, 1, cmap_n);
        gamma_curve_after_gamma = gamma_curve_input .^ gamma_val_calcium;
        gamma_curve_warped = interp1([0, color_pivot_calcium, 1], [0, 0.5, 1], ...
            gamma_curve_after_gamma, 'linear', 'extrap');
        gamma_curve_warped = min(max(gamma_curve_warped, 0), 1);
        if gamma_val_calcium > 0
            gamma_curve_pivot_input = color_pivot_calcium .^ (1 / gamma_val_calcium);
        else
            gamma_curve_pivot_input = NaN;
        end

        colormap_curve_fig = figure( ...
            'Name', 'Calcium Heatmap Colormap Transfer Curve', ...
            'Color', 'w', ...
            'Units', 'pixels', ...
            'Position', [120, 120, 920, 620]);
        curve_layout = tiledlayout(colormap_curve_fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
        ax_curve = nexttile(curve_layout, 1);
        plot(ax_curve, gamma_curve_input, gamma_curve_after_gamma, '--', 'Color', [0.25 0.25 0.25], 'LineWidth', 1.3);
        hold(ax_curve, 'on');
        plot(ax_curve, gamma_curve_input, gamma_curve_warped, 'b-', 'LineWidth', 1.8);
        if isfinite(gamma_curve_pivot_input)
            xline(ax_curve, gamma_curve_pivot_input, ':', sprintf('pivot x=%.3g', gamma_curve_pivot_input), ...
                'Color', [0.2 0.2 0.2], 'LabelVerticalAlignment', 'bottom');
        end
        yline(ax_curve, 0.5, ':', 'colormap midpoint', 'Color', [0.45 0.45 0.45]);
        xlabel(ax_curve, 'Normalized calcium value after clim');
        ylabel(ax_curve, 'Colormap lookup index');
        title(ax_curve, sprintf('Calcium colormap transfer | gamma=%.3g, pivot=%.3g', ...
            gamma_val_calcium, color_pivot_calcium));
        legend(ax_curve, {'x^{gamma}', 'warped final index'}, 'Location', 'southeast');
        grid(ax_curve, 'on');
        ylim(ax_curve, [0 1]);

        ax_strip = nexttile(curve_layout, 2);
        image(ax_strip, gamma_curve_input, 1, reshape(calcium_heatmap_cmap, [1, cmap_n, 3]));
        set(ax_strip, 'YTick', [], 'TickDir', 'out');
        xlabel(ax_strip, 'Normalized calcium value after clim');
        title(ax_strip, 'Resulting black-blue-white color strip');
        xlim(ax_strip, [0 1]);

        colormap_curve_fig_file = fullfile(save_path, '4_dual_roi_calcium_heatmap_colormap_curve.fig');
        colormap_curve_png_file = fullfile(save_path, '4_dual_roi_calcium_heatmap_colormap_curve.png');
        save_figure_bundle_preserve_layout(colormap_curve_fig, colormap_curve_fig_file, colormap_curve_png_file);
        close(colormap_curve_fig);

        % Build a tall single-axis figure: the heatmap supplies row
        % backgrounds and every voltage trace is drawn in the same ROI row
        % coordinate system.
        fig_height = min(2200, max(650, 24 * nrois_heatmap + 180));
        roi_heatmap_fig = figure( ...
            'Name', 'ROI Calcium Heatmap With Voltage Sensitivity Trace', ...
            'Color', 'w', ...
            'Units', 'pixels', ...
            'Position', [80, 80, 1500, fig_height]);
        ax = axes(roi_heatmap_fig);
        imagesc(ax, t_calcium_heatmap, 1:nrois_heatmap, calcium_heatmap_display_zeroed');
        set(ax, 'YDir', 'normal', 'TickDir', 'out', 'Layer', 'top');
        colormap(ax, calcium_heatmap_cmap);
        clim(ax, calcium_heatmap_clim);
        hold(ax, 'on');
        for roi_idx = 1:nrois_heatmap
            % y=roi_idx is the row center. After subtracting each ROI's own
            % midpoint, division by voltage_global_range maps the largest
            % ROI swing to exactly +/-0.5 row units.
            y_trace = roi_idx + (voltage_heatmap_display(:, roi_idx) - voltage_roi_mid(roi_idx)) / voltage_global_range;
            plot(ax, t_voltage_heatmap, y_trace, 'Color', [1.0, 0.12, 0.02], 'LineWidth', 0.45);
        end
        xlim(ax, [min(t_calcium_heatmap), max(t_calcium_heatmap)]);
        ylim(ax, [0.5, nrois_heatmap + 0.5]);
        if nrois_heatmap <= 20
            yticks(ax, 1:nrois_heatmap);
        else
            roi_tick_step = max(1, ceil(nrois_heatmap / 20));
            yticks(ax, unique([1:roi_tick_step:nrois_heatmap, nrois_heatmap]));
        end
        xlabel(ax, 'Time (s)');
        ylabel(ax, 'ROI');
        title(ax, sprintf('Calcium Sensitivity Heatmap (%s) + Voltage Sensitivity Trace', ...
            strrep(char(calcium_heatmap_stage), '_', '\_')));
        cb = colorbar(ax);
        ylabel(cb, 'Calcium sensitivity display value');
        grid(ax, 'on');
        box(ax, 'off');

        roi_heatmap_fig_file = fullfile(save_path, '4_dual_roi_calcium_heatmap_voltage_trace.fig');
        roi_heatmap_png_file = fullfile(save_path, '4_dual_roi_calcium_heatmap_voltage_trace.png');
        roi_heatmap_mat_file = fullfile(save_path, '4_dual_roi_calcium_heatmap_voltage_trace.mat');
        save_figure_bundle_preserve_layout(roi_heatmap_fig, roi_heatmap_fig_file, roi_heatmap_png_file);
        close(roi_heatmap_fig);

        % Persist a compact provenance record so analysis_only reruns can be
        % audited later: selected stages, polarity, scaling reference ROI,
        % color limits, ROI counts, and time-axis repair rules.
        roi_heatmap_trace_result = struct( ...
            'status', "completed", ...
            'reason', "", ...
            'fig_file', string(roi_heatmap_fig_file), ...
            'png_file', string(roi_heatmap_png_file), ...
            'mat_file', string(roi_heatmap_mat_file), ...
            'colormap_curve_fig_file', string(colormap_curve_fig_file), ...
            'colormap_curve_png_file', string(colormap_curve_png_file), ...
            'input_stages', struct('voltage', string(voltage_heatmap_stage), 'calcium', string(calcium_heatmap_stage)), ...
            'parameters', struct( ...
                'voltage_polarity', voltage_polarity, ...
                'calcium_polarity', calcium_polarity, ...
                'calcium_smoothing_window', calcium_smoothing_window, ...
                'calcium_colormap', "custom_ice_adjust_inline", ...
                'calcium_colormap_gamma', gamma_val_calcium, ...
                'calcium_colormap_pivot', color_pivot_calcium, ...
                'calcium_zero_rule', "per-ROI 1st percentile subtracted, then values below zero clipped", ...
                'calcium_zero_percentile', calcium_heatmap_zero_percentile, ...
                'calcium_clim', calcium_heatmap_clim, ...
                'calcium_clim_rule', string(calcium_heatmap_clim_rule), ...
                'voltage_scale_rule', "global max-min range; each ROI centered by its own midpoint", ...
                'voltage_global_range', voltage_global_range, ...
                'voltage_reference_roi', voltage_reference_roi, ...
                'roi_count_plotted', nrois_heatmap, ...
                'roi_count_voltage', nrois_voltage_heatmap, ...
                'roi_count_calcium', nrois_calcium_heatmap, ...
                'x_axis_source', "calcium time axis", ...
                'voltage_trace_time_rule', "overlay voltage trace using voltage seconds without resampling"), ...
            'time_axis_info', struct('voltage', voltage_time_info, 'calcium', calcium_time_info), ...
            'created_at', datetime("now"));

        % Save the actual display matrices used by the figure. These are
        % polarity-adjusted and truncated to paired ROI columns, so they are
        % intentionally not identical to the raw trace stage matrices.
        save(roi_heatmap_mat_file, ...
            'roi_heatmap_trace_result', ...
            'voltage_heatmap_display', ...
            'calcium_heatmap_display_raw', ...
            'calcium_heatmap_display_zeroed', ...
            'calcium_heatmap_baseline', ...
            't_voltage_heatmap', ...
            't_calcium_heatmap', ...
            'voltage_roi_min', ...
            'voltage_roi_max', ...
            'voltage_roi_range', ...
            'voltage_reference_roi', ...
            'gamma_curve_input', ...
            'gamma_curve_after_gamma', ...
            'gamma_curve_warped', ...
            'gamma_curve_pivot_input', ...
            '-v7.3');
    end
end
if ~isfield(dual_results, 'visualizations') || ~isstruct(dual_results.visualizations)
    dual_results.visualizations = struct();
end
dual_results.visualizations.roi_calcium_heatmap_voltage_trace = roi_heatmap_trace_result;
save(dual_results_path, 'dual_results', '-v7.3');

print_section('Dual Comparison');
fprintf('Building dual-channel comparison results...\n');
if stim_context.supported
    if reuse_saved_stim_windows && isstruct(stim_results) && isfield(stim_results, 'windows') && ~isempty(stim_results.windows)
        stim_windows = stim_results.windows;
        fprintf('VolPy re-analysis: reusing saved stim windows from the standard dual results.\n');
    elseif reuse_saved_stim_windows && ~(isfield(stim_context, 'logs') && isstruct(stim_context.logs))
        fprintf(['VolPy re-analysis: saved stim windows are unavailable and saved stim_context has no logs. ' ...
            'Rebuilding stim metadata from the current cycle manifests/logs.\n']);
        voltage_idx_local = find(strcmpi(string({dual_info.camera_cfg.role}), "voltage"), 1, 'first');
        calcium_idx_local = find(strcmpi(string({dual_info.camera_cfg.role}), "calcium"), 1, 'first');
        raw_dual_input_cfg_local = struct();
        if isfield(dual_info, 'manual_options') && isstruct(dual_info.manual_options) ...
                && isfield(dual_info.manual_options, 'raw_dual_input_cfg')
            raw_dual_input_cfg_local = dual_info.manual_options.raw_dual_input_cfg;
        end
        [cycle_manifest_reload, record_manifest_reload, ~, input_layout_reload] = resolve_dual_camera_sources( ...
            cycle_path, dual_info.camera_cfg, [], raw_dual_input_cfg_local);
        stim_context = resolve_stim_context( ...
            cycle_path, cycle_manifest_reload, record_manifest_reload, ...
            dual_info.camera_cfg(voltage_idx_local).camera_index, dual_info.camera_cfg(calcium_idx_local).camera_index, ...
            [], input_layout_reload);
        stim_windows = build_visual_stim_windows(stim_context, freq_voltage, freq_calcium, nframes_voltage, nframes_calcium);
    else
        stim_windows = build_visual_stim_windows(stim_context, freq_voltage, freq_calcium, nframes_voltage, nframes_calcium);
    end
    print_stim_context_summary(stim_context, stim_windows, 'Dual Comparison');
else
    print_stim_context_summary(stim_context, struct(), 'Dual Comparison');
end
[voltage_results, calcium_results] = load_channel_results( ...
    voltage_results_path, calcium_results_path, voltage_results, calcium_results);
has_sensitivity_pair = has_trace_stage(voltage_results, 'sensitivity') ...
    && (has_trace_stage(calcium_results, 'sensitivity_smoothed') || has_trace_stage(calcium_results, 'sensitivity'));
has_snr_pair = has_trace_stage(voltage_results, 'snr') ...
    && (has_trace_stage(calcium_results, 'snr_smoothed') || has_trace_stage(calcium_results, 'snr'));
if has_sensitivity_pair
    voltage_sensitivity = fetch_trace_stage(voltage_results, 'sensitivity');
    [calcium_sensitivity, calcium_sensitivity_stage] = resolve_preferred_trace_stage(calcium_results, {'sensitivity_smoothed', 'sensitivity'});
else
    voltage_sensitivity = [];
    calcium_sensitivity = [];
    calcium_sensitivity_stage = '';
end
if has_snr_pair
    voltage_snr = fetch_trace_stage(voltage_results, 'snr');
    [calcium_snr, calcium_snr_stage] = resolve_preferred_trace_stage(calcium_results, {'snr_smoothed', 'snr'});
else
    voltage_snr = [];
    calcium_snr = [];
    calcium_snr_stage = '';
end
comparison_data = struct();
comparison_info = struct();
if has_sensitivity_pair
    comparison_sensitivity = build_dual_metric_comparison( ...
        'sensitivity', voltage_sensitivity, calcium_sensitivity, ...
        t_voltage, t_calcium, calcium_smoothing_window, save_path, stim_windows, voltage_polarity, calcium_polarity);
    comparison_data.sensitivity = comparison_sensitivity.data;
    comparison_info.sensitivity = comparison_sensitivity.info;
end
if has_snr_pair
    comparison_snr = build_dual_metric_comparison( ...
        'snr', voltage_snr, calcium_snr, ...
        t_voltage, t_calcium, calcium_smoothing_window, save_path, stim_windows, voltage_polarity, calcium_polarity);
    comparison_data.snr = comparison_snr.data;
    comparison_info.snr = comparison_snr.info;
end
if ~has_sensitivity_pair && ~has_snr_pair
    error('VolPy dual comparison cannot run because neither sensitivity nor snr stage is available in both channels.');
end
dual_results.comparison = struct( ...
    'data', comparison_data, ...
    'info', struct( ...
        'method', 'dual_metric_comparison_with_quad_and_overlap_plots', ...
        'parameters', struct('calcium_smoothing_window', calcium_smoothing_window), ...
        'input_stages', struct( ...
            'voltage_sensitivity', 'sensitivity', ...
            'voltage_snr', 'snr', ...
            'calcium_sensitivity', string(calcium_sensitivity_stage), ...
            'calcium_snr', string(calcium_snr_stage)), ...
        'available_stage_pairs', struct('sensitivity', has_sensitivity_pair, 'snr', has_snr_pair), ...
        'comparison_files', comparison_info, ...
        'created_at', datetime("now")));
save(dual_results_path, 'dual_results', '-v7.3');

print_section('Stim-Aware Analysis');
fprintf('Running stimulus-specific analysis after dual comparison...\n');
[voltage_results, calcium_results] = load_channel_results( ...
    voltage_results_path, calcium_results_path, voltage_results, calcium_results);
has_voltage_sensitivity = has_trace_stage(voltage_results, 'sensitivity');
has_voltage_snr = has_trace_stage(voltage_results, 'snr');
has_calcium_sensitivity = has_trace_stage(calcium_results, 'sensitivity_smoothed') || has_trace_stage(calcium_results, 'sensitivity');
has_calcium_snr = has_trace_stage(calcium_results, 'snr_smoothed') || has_trace_stage(calcium_results, 'snr');
if has_voltage_sensitivity
    voltage_sensitivity = fetch_trace_stage(voltage_results, 'sensitivity');
else
    voltage_sensitivity = [];
end
if has_voltage_snr
    voltage_snr = fetch_trace_stage(voltage_results, 'snr');
else
    voltage_snr = [];
end
if has_calcium_sensitivity
    [calcium_sensitivity, calcium_sensitivity_stage] = resolve_preferred_trace_stage(calcium_results, {'sensitivity_smoothed', 'sensitivity'});
else
    calcium_sensitivity = [];
    calcium_sensitivity_stage = '';
end
if has_calcium_snr
    [calcium_snr, calcium_snr_stage] = resolve_preferred_trace_stage(calcium_results, {'snr_smoothed', 'snr'});
else
    calcium_snr = [];
    calcium_snr_stage = '';
end
stim_results.info = rmfield_if_exists(stim_context, {'logs', 'method_manifest'});
stim_results.windows = stim_windows;
if stim_context.supported
    print_stim_context_summary(stim_context, stim_windows, 'Stim-Aware Analysis');
    stim_results.status = ternary(isstruct(stim_windows) && isfield(stim_windows, 'supported') && stim_windows.supported, ...
        'supported', 'unsupported_window_layout');
    if stim_windows.supported
        if ~(has_voltage_snr && has_calcium_snr)
            error('Stim-aware analysis requires snr stage in both voltage_results and calcium_results.');
        end
        stim_trace_fig = figure('Color', 'w');
        subplot(2, 2, 1); hold on; title('Voltage Sensitivity'); plot_optional_stage(voltage_polarity * voltage_sensitivity, t_voltage, has_voltage_sensitivity, 'Sensitivity stage unavailable'); add_stim_shading(gca, stim_windows.voltage, stim_windows.condition_index, stim_windows.condition_colors, 0.14, stim_windows.trial_labels, stim_windows.block_labels);
        subplot(2, 2, 2); hold on; title('Voltage SNR'); plot_optional_stage(voltage_polarity * voltage_snr, t_voltage, has_voltage_snr, 'SNR stage unavailable'); add_stim_shading(gca, stim_windows.voltage, stim_windows.condition_index, stim_windows.condition_colors, 0.14, stim_windows.trial_labels, stim_windows.block_labels);
        subplot(2, 2, 3); hold on; title(sprintf('Calcium Sensitivity (%s)', strrep(char(calcium_sensitivity_stage), '_', '\_'))); plot_optional_stage(calcium_polarity * calcium_sensitivity, t_calcium, has_calcium_sensitivity, 'Sensitivity stage unavailable'); add_stim_shading(gca, stim_windows.calcium, stim_windows.condition_index, stim_windows.condition_colors, 0.14, stim_windows.trial_labels, stim_windows.block_labels);
        subplot(2, 2, 4); hold on; title(sprintf('Calcium SNR (%s)', strrep(char(calcium_snr_stage), '_', '\_'))); plot_optional_stage(calcium_polarity * calcium_snr, t_calcium, has_calcium_snr, 'SNR stage unavailable'); add_stim_shading(gca, stim_windows.calcium, stim_windows.condition_index, stim_windows.condition_colors, 0.14, stim_windows.trial_labels, stim_windows.block_labels);
        save_figure_bundle(stim_trace_fig, fullfile(save_path, '4_dual_trace_summary_with_stim.fig'), fullfile(save_path, '4_dual_trace_summary_with_stim.png'));
        close(stim_trace_fig);

        voltage_stim_metrics_snr = compute_stim_trial_metrics(voltage_snr, stim_windows.voltage, voltage_polarity, freq_voltage);
        calcium_stim_metrics_snr = compute_stim_trial_metrics(calcium_snr, stim_windows.calcium, calcium_polarity, freq_calcium);
        stim_results.response = struct( ...
            'analysis_trace_stage', struct( ...
                'snr', struct('voltage', "snr", 'calcium', string(calcium_snr_stage)), ...
                'sensitivity', struct('voltage', "sensitivity", 'calcium', string(calcium_sensitivity_stage))), ...
            'plot_trace_stage', struct('voltage_sensitivity', "sensitivity", 'voltage_snr', "snr", ...
                'calcium_sensitivity', string(calcium_sensitivity_stage), 'calcium_snr', string(calcium_snr_stage)), ...
            'calcium_smoothing_window', calcium_smoothing_window, ...
            'voltage_polarity', voltage_polarity, ...
            'calcium_polarity', calcium_polarity, ...
            'snr', struct('voltage', voltage_stim_metrics_snr, 'calcium', calcium_stim_metrics_snr));

        if has_voltage_sensitivity && has_calcium_sensitivity
            voltage_stim_metrics_sensitivity = compute_stim_trial_metrics(voltage_sensitivity, stim_windows.voltage, voltage_polarity, freq_voltage);
            calcium_stim_metrics_sensitivity = compute_stim_trial_metrics(calcium_sensitivity, stim_windows.calcium, calcium_polarity, freq_calcium);
            stim_results.response.sensitivity = struct('voltage', voltage_stim_metrics_sensitivity, 'calcium', calcium_stim_metrics_sensitivity);
        else
            voltage_stim_metrics_sensitivity = [];
            calcium_stim_metrics_sensitivity = [];
        end

        if stim_windows.is_grating
            [voltage_tuning, calcium_tuning, per_roi_tuning] = compute_dual_grating_tuning_by_roi( ...
                voltage_stim_metrics_snr.delta_mean, calcium_stim_metrics_snr.delta_mean, stim_windows.orientations);
            stim_results.analysis_kind = 'grating_tuning';
            stim_results.tuning = struct('voltage', voltage_tuning, 'calcium', calcium_tuning, 'per_roi', per_roi_tuning);
            per_roi_tuning = save_grating_tuning_by_roi(per_roi_tuning, save_path, '7_stim_tuning_by_roi');
            stim_results.tuning.per_roi = per_roi_tuning;
            save(fullfile(save_path, '7_stim_tuning_summary.mat'), 'voltage_tuning', 'calcium_tuning', 'per_roi_tuning');
            stim_tuning_fig = figure('Color', 'w');
            subplot(2, 2, 1); plot_tuning_population(voltage_tuning.unique_orientations, voltage_tuning.response_by_condition, 'r', 'Voltage');
            subplot(2, 2, 2); plot_tuning_population(calcium_tuning.unique_orientations, calcium_tuning.response_by_condition, 'g', 'Calcium');
            subplot(2, 2, 3); scatter(voltage_tuning.gosi, calcium_tuning.gosi, 28, 'filled'); xlabel('Voltage gOSI'); ylabel('Calcium gOSI'); title('Orientation Selectivity'); grid on;
            subplot(2, 2, 4); scatter(voltage_tuning.pref_dir, calcium_tuning.pref_dir, 28, 'filled'); hold on; plot([0 360], [0 360], 'k--'); xlim([0 360]); ylim([0 360]); xlabel('Voltage Pref. Dir (deg)'); ylabel('Calcium Pref. Dir (deg)'); title('Preferred Direction'); grid on;
            save_figure_bundle(stim_tuning_fig, fullfile(save_path, '7_stim_tuning_summary.fig'), fullfile(save_path, '7_stim_tuning_summary.png'));
            close(stim_tuning_fig);
        else
            voltage_block_summary_snr = summarize_block_by_condition(voltage_snr, stim_windows.voltage, stim_windows.block_labels, voltage_polarity, stim_windows.trial_labels);
            calcium_block_summary_snr = summarize_block_by_condition(calcium_snr, stim_windows.calcium, stim_windows.block_labels, calcium_polarity, stim_windows.trial_labels);
            voltage_delta_summary_snr = summarize_stim_by_condition(voltage_stim_metrics_snr, stim_windows);
            calcium_delta_summary_snr = summarize_stim_by_condition(calcium_stim_metrics_snr, stim_windows);
            stim_results.analysis_kind = 'condition_response';
            stim_results.condition_summary = struct( ...
                'snr', struct( ...
                    'block', struct('voltage', voltage_block_summary_snr, 'calcium', calcium_block_summary_snr), ...
                    'delta', struct('voltage', voltage_delta_summary_snr, 'calcium', calcium_delta_summary_snr)));
            plot_stim_condition_summary(voltage_block_summary_snr, calcium_block_summary_snr, stim_context, save_path, '7_stim_block_summary_snr', 'SNR Block Response', 'SNR');
            plot_stim_condition_summary(voltage_delta_summary_snr, calcium_delta_summary_snr, stim_context, save_path, '7_stim_delta_summary_snr', 'SNR Delta Response', '\Delta SNR (stim - baseline)');
            if ~isempty(voltage_stim_metrics_sensitivity) && ~isempty(calcium_stim_metrics_sensitivity)
                voltage_block_summary_sensitivity = summarize_block_by_condition(voltage_sensitivity, stim_windows.voltage, stim_windows.block_labels, voltage_polarity, stim_windows.trial_labels);
                calcium_block_summary_sensitivity = summarize_block_by_condition(calcium_sensitivity, stim_windows.calcium, stim_windows.block_labels, calcium_polarity, stim_windows.trial_labels);
                voltage_delta_summary_sensitivity = summarize_stim_by_condition(voltage_stim_metrics_sensitivity, stim_windows);
                calcium_delta_summary_sensitivity = summarize_stim_by_condition(calcium_stim_metrics_sensitivity, stim_windows);
                stim_results.condition_summary.sensitivity = struct( ...
                    'block', struct('voltage', voltage_block_summary_sensitivity, 'calcium', calcium_block_summary_sensitivity), ...
                    'delta', struct('voltage', voltage_delta_summary_sensitivity, 'calcium', calcium_delta_summary_sensitivity));
                plot_stim_condition_summary(voltage_block_summary_sensitivity, calcium_block_summary_sensitivity, stim_context, save_path, '7_stim_block_summary_sensitivity', 'Sensitivity Block Response', 'Sensitivity');
                plot_stim_condition_summary(voltage_delta_summary_sensitivity, calcium_delta_summary_sensitivity, stim_context, save_path, '7_stim_delta_summary_sensitivity', 'Sensitivity Delta Response', '\Delta Sensitivity (stim - baseline)');
            end
        end
    end
else
    stim_results.status = 'not_applicable_or_missing_inputs';
end
save(stim_results_path, 'stim_results', '-v7.3');

print_section('Population Time-Frequency Analysis');
fprintf('Running Fourier and wavelet analysis after dual/stim summaries...\n');
[voltage_results, calcium_results] = load_channel_results( ...
    voltage_results_path, calcium_results_path, voltage_results, calcium_results);
time_frequency_results = run_dual_time_frequency_section( ...
    voltage_results, calcium_results, ...
    t_voltage, t_calcium, ...
    freq_voltage, freq_calcium, ...
    stim_windows, ...
    save_path, ...
    voltage_polarity, calcium_polarity);
time_frequency_record = build_section_record( ...
    'Population Time-Frequency Analysis', ...
    'Summarize processed dual-channel traces with direct FFT spectra and wavelet time-frequency maps, while reusing the saved processed trace stages as section inputs.', ...
    struct( ...
        'voltage_results_file', string(voltage_results_path), ...
        'calcium_results_file', string(calcium_results_path), ...
        'stim_results_file', string(stim_results_path), ...
        'voltage_stage', string(time_frequency_results.parameters.voltage_stage), ...
        'calcium_stage', string(time_frequency_results.parameters.calcium_stage)), ...
    time_frequency_results.parameters, ...
    struct( ...
        'fourier_summary_png', string(time_frequency_results.visualizations.fourier_summary.png_file), ...
        'wavelet_summary_png', string(time_frequency_results.visualizations.wavelet_summary.png_file)), ...
    struct(), ...
    "To rerun the VolPy voltage re-analysis pipeline, rerun Dual_analysis3 with analysis_backend = 'volpy_voltage_reanalysis'.");
time_frequency_results.record = time_frequency_record;
time_frequency_results_path = fullfile(save_path, '8_time_frequency_results.mat');
save(time_frequency_results_path, 'time_frequency_results', '-v7.3');

print_section('Save Explicit Results');
save_explicit_dual_results_summary( ...
    save_path, ...
    dual_info, dual_info_path, ...
    voltage_results, voltage_results_path, ...
    calcium_results, calcium_results_path, ...
    dual_results, dual_results_path, ...
    stim_results, stim_results_path, ...
    time_frequency_results, time_frequency_results_path);

fprintf('Dual_analysis3 VolPy voltage re-analysis finished.\n');
fprintf('Results saved to: %s\n', save_path);
end

function data = fetch_trace_stage(results, stage_name)
if ~isstruct(results) || ~isfield(results, 'trace_results')
    error('Channel results are not loaded.');
end

if isfield(results.trace_results, stage_name) && isfield(results.trace_results.(stage_name), 'data')
    data = results.trace_results.(stage_name).data;
    return;
end

error('Trace stage "%s" is not available in the saved channel results.', stage_name);
end

function record = build_section_record(section_name, description, input_struct, parameter_struct, output_struct, excluded_inputs_struct, rerun_note)
% Build a compact "interface description" for one analysis section.
%
% Natural-language meaning:
%   each record is a self-contained summary of what the section did,
%   what went in, what came out, and what would still be needed to rerun it.
%
% Interface:
%   record.input           -> saved non-movie inputs used by this section
%   record.parameters      -> algorithm choices and tunable settings
%   record.output          -> direct outputs of this section
%   record.excluded_inputs -> required inputs that were intentionally not saved
%   record.rerun.note      -> plain-language rerun instruction
%
% Example:
%   record.output = struct('traces_voltage_bg', traces_voltage_bg)
%   record.excluded_inputs.movie_voltage = 'not saved'
record = struct( ...
    'section_name', string(section_name), ...
    'description', string(description), ...
    'input', input_struct, ...
    'parameters', parameter_struct, ...
    'output', output_struct, ...
    'excluded_inputs', excluded_inputs_struct, ...
    'rerun', struct( ...
        'note', string(rerun_note), ...
        'movie_saved', false), ...
    'created_at', datetime("now"));
end

function [data, stage_name] = resolve_preferred_trace_stage(results, preferred_stages)
% Choose the best available saved stage from a priority list. This lets
% newer analyses prefer refined stages while staying compatible with older
% result files that only contain earlier stage names.
for i = 1:numel(preferred_stages)
    stage_name = preferred_stages{i};
    if has_trace_stage(results, stage_name)
        data = fetch_trace_stage(results, stage_name);
        return;
    end
end

error('None of the requested trace stages are available: %s', strjoin(preferred_stages, ', '));
end

function tf = has_trace_stage(results, stage_name)
tf = isstruct(results) ...
    && isfield(results, 'trace_results') ...
    && isfield(results.trace_results, stage_name) ...
    && isfield(results.trace_results.(stage_name), 'data');
end

function s = rmfield_if_exists(s, field_names)
if isempty(s) || ~isstruct(s)
    return;
end
for i = 1:numel(field_names)
    if isfield(s, field_names{i})
        s = rmfield(s, field_names{i});
    end
end
end

function out = ternary(condition, true_value, false_value)
if condition
    out = true_value;
else
    out = false_value;
end
end
