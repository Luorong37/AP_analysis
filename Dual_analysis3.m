% Dual_analysis3 - Dual-camera voltage/calcium companion analysis
% Dual_analysis3 notes:
%   - This script is intended for dual-camera acquisitions saved by the
%     rebuilt Hamamatsu workflow under Rec*/Cycle*/Cam*_label paths.
%   - It prefers cycle_manifest.mat / record_manifest.mat when available,
%     then falls back to Cam1_* / Cam2_* folders inside the selected cycle.
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
%   correct_offset = false;
%   run('Dual_analysis3.m')
%
% Important:
%   this script no longer clears the workspace automatically.
%   - If you run it directly by itself, CLEAR old variables manually first.
%   - If you run it from an outer batch script, do not clear the workspace,
%     because the outer script may be intentionally passing config values in.
clc;

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
%   correct_offset = false;
%   run('Dual_analysis3.m');

%% Input Setup
% This block collects the main user-tunable analysis settings in one place
% so most reruns only require edits here.
nowtime = string(datetime('now'));
nowtime = strrep(nowtime, ':', '-');
fprintf('Initializing dual analysis...\n');

if ~exist('cycle_path', 'var') || isempty(cycle_path)
    cycle_path = 'V:\Luorong\Invivo\26.04.15_dual-color_P195\Methods5_default\Rec1_2026-04-15_20-58-57\Cycle1';
end
[record_path_for_name, cycle_name_for_name] = fileparts(cycle_path);
[~, record_name_for_name] = fileparts(record_path_for_name);
if ~exist('analysis_run_name', 'var') || isempty(analysis_run_name)
    analysis_run_name = sprintf('%s_%s_%s', record_name_for_name, cycle_name_for_name, char(nowtime));
end

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

if ~exist('map_bin', 'var') || isempty(map_bin)
    map_bin = 4;
end
if ~exist('calcium_smoothing_window', 'var') || isempty(calcium_smoothing_window)
    calcium_smoothing_window = 40;
end
if ~exist('downsample_window', 'var') || isempty(downsample_window)
    downsample_window = calcium_smoothing_window; % legacy label used by dual comparison outputs
end
if ~exist('bleach_mode_voltage', 'var') || isempty(bleach_mode_voltage)
    bleach_mode_voltage = 'linear';   % 'linear' | 'highpass' | 'exp2'
end
if ~exist('bleach_mode_calcium', 'var') || isempty(bleach_mode_calcium)
    bleach_mode_calcium = 'exp2';     % 'linear' | 'highpass' | 'exp2'
end
if ~exist('voltage_polarity', 'var') || isempty(voltage_polarity)
    voltage_polarity = -1;            % default display/analysis polarity for voltage traces
end
if ~exist('calcium_polarity', 'var') || isempty(calcium_polarity)
    calcium_polarity = 1;             % default display/analysis polarity for calcium traces
end
if ~exist('correct_offset', 'var') || isempty(correct_offset)
    correct_offset = false;           % true when manually estimating inter-camera ROI offset
end
if ~exist('reuse_roi_file', 'var') || isempty(reuse_roi_file)
    reuse_roi_file = '';              % e.g. fullfile(save_path, '1_dual_roi_results.mat')
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
        'highpass', true);
end

if numel(camera_cfg) ~= 2
    error('Dual_analysis3 currently expects exactly two cameras.');
end

if ~exist('save_path', 'var') || isempty(save_path)
    save_path = fullfile(cycle_path, 'Dual_analysis3', analysis_run_name);
end
mkdir(save_path);

if gpu && isempty(gcp('nocreate'))
    gcp;
end

%% Resolve Rebuilt Inputs
% Prefer rebuilt manifests so the analysis follows the acquisition layout
% recorded during acquisition instead of guessing from folder names alone.
print_section('Resolve Rebuilt Inputs');
[cycle_manifest, record_manifest, camera_source] = resolve_dual_camera_sources(cycle_path, camera_cfg);

voltage_idx = find(strcmpi(string({camera_cfg.role}), "voltage"), 1, 'first');
calcium_idx = find(strcmpi(string({camera_cfg.role}), "calcium"), 1, 'first');
if isempty(voltage_idx) || isempty(calcium_idx)
    error('camera_cfg must contain one voltage camera and one calcium camera.');
end

stim_context = resolve_stim_context( ...
    cycle_path, cycle_manifest, record_manifest, ...
    camera_cfg(voltage_idx).camera_index, camera_cfg(calcium_idx).camera_index);

dual_info = struct();
dual_info.analysis_name = 'Dual_analysis3';
dual_info.cycle_path = cycle_path;
dual_info.save_path = save_path;
dual_info.created_at = datetime("now");
dual_info.map_bin = map_bin;
dual_info.downsample_window = downsample_window;
dual_info.calcium_smoothing_window = calcium_smoothing_window;
dual_info.bleach_mode = struct( ...
    'voltage', bleach_mode_voltage, ...
    'calcium', bleach_mode_calcium);
dual_info.polarity = struct( ...
    'voltage', voltage_polarity, ...
    'calcium', calcium_polarity);
dual_info.correct_offset = correct_offset;
dual_info.camera_cfg = camera_cfg;
dual_info.camera_source = camera_source;
dual_info.role_order = struct( ...
    'voltage_camera_index', camera_cfg(voltage_idx).camera_index, ...
    'calcium_camera_index', camera_cfg(calcium_idx).camera_index);
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
currentScript = mfilename('fullpath');
if isempty(currentScript)
    currentScript = which('Dual_analysis3.m');
end
copy_analysis_code(currentScript, code_path);

%% Motion Correction
% Motion is estimated once on voltage and then applied to calcium.
% The processing intent is to keep both channels spatially locked. If the
% two channels were corrected independently, one biological ROI could end
% up drifting to different places in voltage and calcium.
print_section('Motion Correction');
fprintf('Applying shared motion correction from voltage to calcium...\n');
fprintf('Motion model: rigid NoRMCorre estimated on voltage, then apply_shifts to calcium.\n');
[movie_voltage_3d, movie_calcium_3d, voltage_motion_info, calcium_motion_info] = run_shared_motion_correction( ...
    movie_voltage_3d, movie_calcium_3d, save_path, motion_cfg);

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
fprintf('Motion complete | shared shift file: %s\n', voltage_motion_info.shift_file);

%% Create Sensitivity Maps
% These maps are quick summary images used mainly to guide ROI selection.
% They are visual aids rather than final quantitative outputs.
print_section('Create Sensitivity Maps');
fprintf('Creating voltage/calcium maps...\n');
fprintf('Map method: create_map(movie, nrows, ncols, map_bin), map_bin=%d\n', map_bin);
map_voltage = create_map(movie_voltage, nrows, ncols, map_bin);
map_calcium = create_map(movie_calcium, nrows, ncols, map_bin, 'calcium');

figure('Color', 'w');
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
saveas(gcf, fullfile(save_path, '0_dual_sensitivity_map.fig'), 'fig');
saveas(gcf, fullfile(save_path, '0_dual_sensitivity_map.png'), 'png');

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
    if isfield(roi_cache, 'rois')
        mask_voltage = roi_cache.rois.bwmask;
        mask_calcium = roi_cache.rois.bwmask_ca;
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

[rois, traces_voltage_raw, traces_calcium_raw, offset] = select_ROI_dual( ...
    movie_voltage, movie_calcium, nrows, ncols, correct_offset, ...
    map_voltage, map_calcium, mask_voltage, mask_calcium, offset);
% This is the main transition from movie space to trace space. After this
% point, most analysis steps work on ROI-by-time matrices rather than raw
% image stacks.

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
        'correct_offset', correct_offset, ...
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
        'offset_mode', ternary(correct_offset, "manual", "reused_or_zero"), ...
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
plot_dual_raw_traces( ...
    traces_voltage_raw, traces_calcium_raw, ...
    t_voltage, t_calcium);
saveas(gcf, raw_trace_fig, 'fig');
saveas(gcf, raw_trace_png, 'png');
dual_results.visualizations.raw = struct( ...
    'data', struct(), ...
    'info', struct( ...
        'fig_file', raw_trace_fig, ...
        'png_file', raw_trace_png, ...
        'created_at', datetime("now")));
save(dual_results_path, 'dual_results', '-v7.3');
fprintf('Raw trace overview saved to: %s\n', raw_trace_png);
fprintf('ROI complete | nrois=%d | offset=[%.3f %.3f] | map source=%s | roi file=%s\n', ...
    nrois, offset(1), offset(2), roi_map_source, roi_results_file);

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
fprintf('Removing background for both channels...\n');
[voltage_results, calcium_results] = load_channel_results( ...
    voltage_results_path, calcium_results_path, voltage_results, calcium_results);
[dual_results, rois, traces_voltage_raw, traces_calcium_raw, roi_results_file, nrois] = ...
    load_dual_roi_context(dual_results_path, dual_results);

rois_voltage = struct('bwmask', rois.bwmask, 'boundary', {rois.boundary}, 'position', {rois.position});
rois_calcium = struct('bwmask', rois.bwmask_ca, 'boundary', {rois.boundary_ca}, 'position', {rois.position_ca});

if ~exist('movie_voltage', 'var') || ~exist('movie_calcium', 'var')
    if has_trace_stage(voltage_results, 'bg_removed') && has_trace_stage(calcium_results, 'bg_removed')
        fprintf('Background stage already exists. Reusing saved bg_removed traces without recomputation.\n');
        fprintf('Current bg_removed parent stage: voltage=%s | calcium=%s\n', ...
            strjoin(string(voltage_results.trace_results.bg_removed.info.parent_results), ','), ...
            strjoin(string(calcium_results.trace_results.bg_removed.info.parent_results), ','));
    else
        error(['Background removal needs movie data in memory to recompute. ' ...
            'Run the loading/motion/ROI sections first, or skip this section and let downstream sections use raw traces.']);
    end
else
    [background_voltage, background_fit_voltage, ~, traces_voltage_bg, background_mask_voltage] = ...
        remove_background(movie_voltage, ncols, nrows, rois_voltage, freq_voltage, 1);
    [background_calcium, background_fit_calcium, ~, traces_calcium_bg, background_mask_calcium] = ...
        remove_background(movie_calcium, ncols, nrows, rois_calcium, freq_calcium, 1);

    background_results_file = fullfile(save_path, '1_dual_background_results.mat');
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

background_results_file = fullfile(save_path, '1_dual_background_results.mat');
if isfile(background_results_file)
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
    fprintf('Background summary plot skipped because background result file is missing.\n');
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
figure('Color', 'w');
subplot(2, 3, 1); hold on; title('Voltage Raw'); plot_offset_stage(traces_voltage_raw, t_voltage);
subplot(2, 3, 2); hold on; title('Voltage Sensitivity'); plot_optional_stage(voltage_polarity * voltage_sensitivity, t_voltage, has_voltage_sensitivity, 'Sensitivity stage unavailable');
subplot(2, 3, 3); hold on; title('Voltage SNR'); plot_optional_stage(voltage_polarity * voltage_snr, t_voltage, has_voltage_snr, 'SNR stage unavailable');
subplot(2, 3, 4); hold on; title(sprintf('Calcium Raw (%s, window=%d)', strrep(char(calcium_raw_stage), '_', '\_'), calcium_smoothing_window)); plot_optional_stage(traces_calcium_raw, t_calcium, has_calcium_raw, 'Raw stage unavailable');
subplot(2, 3, 5); hold on; title(sprintf('Calcium Sensitivity (%s, window=%d)', strrep(char(calcium_sensitivity_stage), '_', '\_'), calcium_smoothing_window)); plot_optional_stage(calcium_polarity * calcium_sensitivity, t_calcium, has_calcium_sensitivity, 'Sensitivity stage unavailable');
subplot(2, 3, 6); hold on; title(sprintf('Calcium SNR (%s, window=%d)', strrep(char(calcium_snr_stage), '_', '\_'), calcium_smoothing_window)); plot_optional_stage(calcium_polarity * calcium_snr, t_calcium, has_calcium_snr, 'SNR stage unavailable');
saveas(gcf, fullfile(save_path, '4_dual_trace_summary.fig'), 'fig');
saveas(gcf, fullfile(save_path, '4_dual_trace_summary.png'), 'png');
fprintf('Summary plot stage availability | voltage sensitivity=%d snr=%d | calcium sensitivity=%d snr=%d\n', ...
    has_voltage_sensitivity, has_voltage_snr, has_calcium_sensitivity, has_calcium_snr);
fprintf('Channel summary polarity | voltage=%d | calcium=%d\n', voltage_polarity, calcium_polarity);
fprintf('Channel summary calcium stages | raw=%s | sensitivity=%s | snr=%s\n', ...
    string(calcium_raw_stage), string(calcium_sensitivity_stage), string(calcium_snr_stage));
fprintf('Channel summary calcium final smoothing window=%d\n', calcium_smoothing_window);

%% Dual Comparison
% This section asks whether the processed voltage and calcium traces tell a
% consistent ROI-by-ROI story after each channel finishes its own
% preprocessing.
print_section('Dual Comparison');
fprintf('Building dual-channel comparison results...\n');
stim_windows = struct();
if stim_context.supported
    stim_windows = build_visual_stim_windows(stim_context, freq_voltage, freq_calcium, nframes_voltage, nframes_calcium);
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
    t_voltage, t_calcium, downsample_window, save_path, stim_windows, voltage_polarity, calcium_polarity);
    comparison_data.sensitivity = comparison_sensitivity.data;
    comparison_info.sensitivity = comparison_sensitivity.info;
else
    fprintf('Skipping dual sensitivity comparison: sensitivity stage missing in saved channel results\n');
end

if has_snr_pair
    comparison_snr = build_dual_metric_comparison( ...
        'snr', voltage_snr, calcium_snr, ...
        t_voltage, t_calcium, downsample_window, save_path, stim_windows, voltage_polarity, calcium_polarity);
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
        'parameters', struct('downsample_window', downsample_window, 'calcium_smoothing_window', calcium_smoothing_window), ...
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

        figure('Color', 'w');
        subplot(2, 2, 1); hold on; title('Voltage Sensitivity'); plot_optional_stage(voltage_polarity * voltage_sensitivity, t_voltage, has_voltage_sensitivity, 'Sensitivity stage unavailable'); add_stim_shading(gca, stim_windows.voltage, stim_windows.condition_index, stim_windows.condition_colors, 0.14, stim_windows.trial_labels, stim_windows.block_labels);
        subplot(2, 2, 2); hold on; title('Voltage SNR'); plot_optional_stage(voltage_polarity * voltage_snr, t_voltage, has_voltage_snr, 'SNR stage unavailable'); add_stim_shading(gca, stim_windows.voltage, stim_windows.condition_index, stim_windows.condition_colors, 0.14, stim_windows.trial_labels, stim_windows.block_labels);
        subplot(2, 2, 3); hold on; title(sprintf('Calcium Sensitivity (%s)', strrep(char(calcium_sensitivity_stage), '_', '\_'))); plot_optional_stage(calcium_polarity * calcium_sensitivity, t_calcium, has_calcium_sensitivity, 'Sensitivity stage unavailable'); add_stim_shading(gca, stim_windows.calcium, stim_windows.condition_index, stim_windows.condition_colors, 0.14, stim_windows.trial_labels, stim_windows.block_labels);
        subplot(2, 2, 4); hold on; title(sprintf('Calcium SNR (%s)', strrep(char(calcium_snr_stage), '_', '\_'))); plot_optional_stage(calcium_polarity * calcium_snr, t_calcium, has_calcium_snr, 'SNR stage unavailable'); add_stim_shading(gca, stim_windows.calcium, stim_windows.condition_index, stim_windows.condition_colors, 0.14, stim_windows.trial_labels, stim_windows.block_labels);
        saveas(gcf, fullfile(save_path, '4_dual_trace_summary_with_stim.fig'), 'fig');
        saveas(gcf, fullfile(save_path, '4_dual_trace_summary_with_stim.png'), 'png');

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
            voltage_tuning = compute_grating_tuning(voltage_stim_metrics_snr.delta_mean, stim_windows.orientations);
            calcium_tuning = compute_grating_tuning(calcium_stim_metrics_snr.delta_mean, stim_windows.orientations);
            stim_results.analysis_kind = 'grating_tuning';
            stim_results.tuning = struct( ...
                'voltage', voltage_tuning, ...
                'calcium', calcium_tuning);

            save(fullfile(save_path, '7_stim_tuning_summary.mat'), 'voltage_tuning', 'calcium_tuning');

            figure('Color', 'w');
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
            saveas(gcf, fullfile(save_path, '7_stim_tuning_summary.fig'), 'fig');
            saveas(gcf, fullfile(save_path, '7_stim_tuning_summary.png'), 'png');
        else
            fprintf('Stim analysis branch: condition response\n');
            fprintf('  OSI/DSI: skipped\n');
            fprintf('  reason: non-grating stimulus type = %s\n', stim_context.stim_type);
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

%% Save Explicit Results
% Save one bundled snapshot for manual inspection. The modular section
% files remain the main working outputs, while this file is a convenient
% "open everything at once" archive.
print_section('Save Explicit Results');
explicit_dual_results = struct( ...
    'dual_info', dual_info, ...
    'voltage_results', voltage_results, ...
    'calcium_results', calcium_results, ...
    'dual_results', dual_results, ...
    'stim_results', stim_results);
save(fullfile(save_path, '-1_explicit_dual_results.mat'), 'explicit_dual_results', '-v7.3');

fprintf('Dual_analysis3 finished.\n');
fprintf('Results saved to: %s\n', save_path);

function [cycle_manifest, record_manifest, camera_source] = resolve_dual_camera_sources(cycle_path, camera_cfg)
% Resolve where each camera movie came from and which label belongs to it.
% The intent is to keep analysis inputs tied to acquisition metadata when
% available, while still supporting older folder-only datasets.
cycle_manifest = [];
record_manifest = [];
camera_source = repmat(struct('path', '', 'label', "", 'original_frame_size', [NaN, NaN]), 1, numel(camera_cfg));

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
if isempty(currentScript)
    return;
end
[requiredFiles, ~] = matlab.codetools.requiredFilesAndProducts(currentScript);
for k = 1:length(requiredFiles)
    [~, name, ext] = fileparts(requiredFiles{k});
    copyfile(requiredFiles{k}, fullfile(code_path, [name, ext]));
end
fprintf('All required code files were copied to %s\n', code_path);
end

function [voltage_corrected, calcium_corrected, voltage_motion_info, calcium_motion_info] = run_shared_motion_correction(voltage_movie, calcium_movie, save_path, cfg)
% Estimate one shared motion model on voltage, then apply the same shifts
% to calcium so both channels remain spatially coupled after correction.
voltage_corrected = voltage_movie;
calcium_corrected = calcium_movie;

voltage_motion_info = struct( ...
    'applied', false, ...
    'method', '', ...
    'shift_file', '', ...
    'parameter_file', '', ...
    'highpass', cfg.highpass, ...
    'use_saved_shift', cfg.use_saved_shift, ...
    'shared_with_role', 'calcium');
calcium_motion_info = struct( ...
    'applied', false, ...
    'method', '', ...
    'shift_file', '', ...
    'parameter_file', '', ...
    'highpass', cfg.highpass, ...
    'use_saved_shift', cfg.use_saved_shift, ...
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

if cfg.use_saved_shift && ~isempty(cfg.saved_shift_file)
    % Reusing a previously estimated shift field is useful when iterating
    % on later sections without recomputing motion every time.
    shift_data = load(cfg.saved_shift_file);
    if isfield(shift_data, 'options_r')
        options_r = shift_data.options_r;
    end
    shifts_r = shift_data.shifts_r;
    copyfile(cfg.saved_shift_file, shift_res_path);
else
    if cfg.highpass
        % High-pass preprocessing helps the registration focus more on
        % structure than on slow brightness drift.
        movie_for_estimation = create_temp_highpass(voltage_single);
    else
        movie_for_estimation = voltage_single;
    end
    [~, shifts_r, ~] = normcorre_batch(movie_for_estimation, options_r);
    save(shift_res_path, 'shifts_r', 'options_r', '-v7.3');
end

voltage_corrected = apply_shifts(voltage_single, shifts_r, options_r);
calcium_corrected = apply_shifts(calcium_single, shifts_r, options_r);
save(params_save_path, 'options_r', 'cfg');

voltage_motion_info.applied = true;
voltage_motion_info.method = 'NoRMCorre_rigid_shared_voltage_reference';
voltage_motion_info.shift_file = shift_res_path;
voltage_motion_info.parameter_file = params_save_path;

calcium_motion_info.applied = true;
calcium_motion_info.method = 'reuse_voltage_motion_shifts';
calcium_motion_info.shift_file = shift_res_path;
calcium_motion_info.parameter_file = params_save_path;
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

function stim_context = resolve_stim_context(cycle_path, cycle_manifest, record_manifest, voltage_camera_index, calcium_camera_index)
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
    'block_labels', strings(0, 1));
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
    if iscell(raw_labels)
        labels = string(raw_labels(:));
    else
        labels = string(raw_labels(:));
    end
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

if strcmp(selected_program, 'gray_blue_gray') || contains(selected_program, 'blue') || contains(selected_label, 'blue')
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
fprintf('Stim context summary [%s]:\n', caller_name);
fprintf('  record mode: %s\n', string(stim_context.recordmode));
fprintf('  selected program: %s\n', string(stim_context.selected_program));
fprintf('  selected label: %s\n', string(stim_context.selected_label));
fprintf('  classified stim type: %s\n', string(stim_context.stim_type));
fprintf('  logs available: %d | stimSpec available: %d | sync available: %d\n', ...
    stim_context.has_logs, stim_context.has_stimSpec, ...
    stim_context.has_logs && isfield(stim_context.logs, 'sync') && istable(stim_context.logs.sync));
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

if isfield(channel_windows, 'block_time_ranges') && ~isempty(channel_windows.block_time_ranges) && numel(block_labels) == size(channel_windows.block_time_ranges, 1)
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

if contains(label, 'black')
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
saveas(fig, fig_filename, 'fig');
saveas(fig, png_filename, 'png');
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
saveas(fig, fullfile(save_path, [file_stem '.fig']), 'fig');
saveas(fig, fullfile(save_path, [file_stem '.png']), 'png');
end

function comparison = build_dual_metric_comparison(metric_name, voltage_traces, calcium_traces, t_voltage, t_calcium, downsample_window, save_path, stim_windows, voltage_polarity, calcium_polarity)
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
overlap_fig = fullfile(save_path, sprintf('5_dual_%s_overlap.fig', metric_name));
overlap_png = fullfile(save_path, sprintf('5_dual_%s_overlap.png', metric_name));
correlation_fig = fullfile(save_path, sprintf('6_dual_%s_correlation.fig', metric_name));
correlation_png = fullfile(save_path, sprintf('6_dual_%s_correlation.png', metric_name));
mat_file = fullfile(save_path, sprintf('6_dual_%s_comparison.mat', metric_name));

plot_dual_quad_summary( ...
    t_calcium, t_voltage, calcium_normalized, calcium_deconv_display, voltage_integral, voltage_display, ...
    size(calcium_metric, 2), metric_name, stim_windows);
saveas(gcf, quad_fig, 'fig');
saveas(gcf, quad_png, 'png');

plot_dual_overlap_summary( ...
    t_calcium, t_voltage, calcium_display, voltage_display, ...
    size(calcium_metric, 2), metric_name, stim_windows);
saveas(gcf, overlap_fig, 'fig');
saveas(gcf, overlap_png, 'png');

plot_dual_correlation_summary(correlation_stats, metric_name);
saveas(gcf, correlation_fig, 'fig');
saveas(gcf, correlation_png, 'png');

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
        'legacy_downsample_window_argument', downsample_window, ...
        'quad_fig', quad_fig, ...
        'quad_png', quad_png, ...
        'overlap_fig', overlap_fig, ...
        'overlap_png', overlap_png, ...
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

function plot_dual_overlap_summary(t_calcium, t_voltage, calcium_trace, voltage_trace, nrois, metric_name, stim_windows)
xlimit = [0, max([t_calcium(:); t_voltage(:)])];

fig = figure('Name', sprintf('Dual %s Overlap', upper(metric_name)), ...
    'Color', 'w', 'Position', [100, 100, 1200, max(400, 220 * nrois)]);

for i = 1:nrois
    ax = subplot(nrois, 1, i);
    hold(ax, 'on');
    if nargin >= 7 && isstruct(stim_windows) && isfield(stim_windows, 'voltage') && isfield(stim_windows.voltage, 'stim_time_ranges')
        add_stim_shading(ax, stim_windows.voltage, stim_windows.condition_index, stim_windows.condition_colors, 0.12, stim_windows.trial_labels, stim_windows.block_labels);
    end

    voltage_norm = normalize_columns_to_unit_range(voltage_trace(:, i));
    calcium_norm = normalize_columns_to_unit_range(calcium_trace(:, i));

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
    title(ax, sprintf('ROI %d', i), 'FontSize', 9, 'FontWeight', 'bold');

    if i < nrois
        ax.XTickLabel = [];
    else
        xlabel(ax, 'Time (s)');
    end
end

sgtitle(sprintf('Dual %s Overlap: Dual yyaxis, separately normalized', upper(metric_name)));
set(fig, 'Position', get(0, 'Screensize'));
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
saveas(fig, fig_filename, 'fig');
saveas(fig, png_filename, 'png');
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
saveas(fig, fig_filename, 'fig');
saveas(fig, png_filename, 'png');
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
