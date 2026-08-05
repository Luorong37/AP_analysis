% run_Dual_analysis3_rec - Run Dual_analysis3 across all cycles in one record
%
% Workflow:
% 1. Choose one reference cycle.
% 2. Create or reuse one ROI file from that reference cycle.
% 3. Reuse the same ROI file for every other Cycle* under the same record.
%
% Two main analysis workflows:
%   1. dual3_rec_control.preset='full': create one new timestamped result folder per
%      cycle and run the complete analysis.
%   2. dual3_rec_control.preset='analysis_only': find one previous complete result
%      folder per cycle, reuse its ROI/traces/accepted peaks/stim context,
%      and write only the downstream rerun into a new timestamped folder.
%
% Analysis-only fork details:
%   Set dual3_rec_control.preset='analysis_only'. The
%   script will find each Cycle*/Dual_analysis3 source result, reuse its
%   saved traces/ROI/info, and write the new analysis outputs under
%   Cycle*/Dual_analysis3/<run_name>. The source folder is
%   read as input and is not used as save_path.
%
% Motion-only workflow:
%   Set dual3_rec_control.preset='motion_only'. dual3_rec_control.rec_path may
%   point either to one Rec*
%   folder or to a Methods* folder containing multiple Rec* folders. The
%   script runs Dual_analysis3 motion correction for every discovered
%   Rec*/Cycle* and stops each cycle immediately after motion outputs are
%   saved.
%
% This script keeps Dual_analysis3 as a normal runnable script.
% It passes parameters in from an outer scope by calling run(...).

%% Input Setup
% Recommended workflow control:
%
%   dual3_rec_control = struct();
%   dual3_rec_control.preset = "full";
%   dual3_rec_control.rec_path = '...\Rec1_...';
%
% Or reuse each Cycle's latest complete result:
%
%   dual3_rec_control = struct();
%   dual3_rec_control.preset = "analysis_only";
%   dual3_rec_control.rec_path = '...\Rec1_...';
%
% dual3_rec_control is the only batch workflow-control input. Missing fields
% use the preset below; unrelated workspace variables are ignored.
% Signal-processing parameters below remain separate and are passed to
% Dual_analysis3 only when the selected preset needs them.

% -------------------------------------------------------------------------
% A. Workflow control presets
% -------------------------------------------------------------------------
% dual3_rec_control fields:
%
%   rec_path
%       Batch scope. For "full" and "analysis_only", use one Rec* folder
%       containing Cycle* folders. "motion_only" also accepts a Methods*
%       folder and discovers Rec*/Cycle* folders below it.
%
%   preset
%       Selects the per-cycle workflow:
%       "full"          calls Dual_analysis3 preset="full" for each Cycle.
%       "reuse_motion"  loads movies and applies saved motion shifts.
%       "reuse_trace"   reuses saved traces and reruns peaks/stim analysis.
%       "skip_motion"   loads movies but skips motion and registration.
%       "custom"        requires all section actions explicitly.
%                       These direct presets forward sections/source controls.
%                       When reuse_reference_roi=true, the reference Cycle ROI
%                       is passed separately through roi_path to every Cycle.
%       "analysis_only" finds a reusable result in each Cycle, then calls
%                       Dual_analysis3 preset="analysis_only".
%       "motion_only"   calls Dual_analysis3 preset="motion_only" and stops
%                       each Cycle after motion outputs.
%
%   sections
%       Struct forwarded to Dual_analysis3 as dual3_control.sections in
%       direct named/custom presets. Empty means use the selected preset unchanged.
%
%   source_results_path
%       Forwarded to Dual_analysis3 as dual3_control.source_results_path in
%       direct named/custom presets. Empty means no reused source result is supplied.
%
%   roi_path
%       Forwarded to Dual_analysis3 as dual3_control.roi_path in direct preset
%       modes. When shared-reference ROI is enabled, reference_roi_path
%       takes priority, then roi_path, then automatic reference-Cycle lookup.
%
%   output_root
%       Destination root folder name used by analysis_only under each Cycle.
%       The default "Dual_analysis3" writes beside normal Dual3 results.
%       It does not change where full or motion_only writes.
%
%   run_tag
%       Timestamp/name suffix shared by the analysis_only batch. Empty creates
%       one timestamp when run_Dual_analysis3_rec starts. The complete per-cycle
%       folder remains <Rec>_<Cycle>_<run_tag>.
%
%   reference_cycle_name
%       Cycle used to create or locate the shared ROI in full and direct modes.
%       It is not used to choose traces in analysis_only mode.
%
%   reference_roi_path
%       Folder containing the shared ROI for full/direct mode. Dual3 and AP3 ROI
%       files are supported. A valid value bypasses automatic ROI discovery
%       unless redraw_reference_roi=true. This folder is used only as
%       Dual_analysis3 roi_path, not as the source for motion, trace, peak,
%       or other reused sections.
%
%   reuse_reference_roi
%       Full/direct mode. When true and reference_roi_path is empty, load the
%       newest supported ROI already present in the reference Cycle. When
%       false, create a new reference ROI if no explicit folder was supplied.
%
%   redraw_reference_roi
%       Full mode only. When true, rerun the reference Cycle and create a new
%       ROI even if an old or explicit reference ROI is available.
%
%   reference_offset_mode
%       Registration method used while creating the reference ROI:
%       "none", "manual_points", or "matlab_register".
%
%   include_reference_cycle
%       Full mode only. false excludes the reference Cycle from the later
%       batch loop because it already supplied the ROI; true runs it again as
%       a normal selected Cycle and creates another timestamped result.
%
%   cycles
%       Optional string list such as ["Cycle1","Cycle3"]. Empty processes all
%       discovered Cycle* folders. Nonmatching cycles are excluded.
%
%   skip_existing
%       true skips cycles with a completed matching output; false always creates
%       another timestamped run. The marker is a final explicit bundle, or a
%       completed motion output for motion_only. Empty defaults
%       to true only for motion_only and false for all other presets.
%
%   record_average_only
%       true skips all per-cycle Dual3 execution, finds existing cycle result
%       folders, and only generates the Rec-level stimulus-specific average.
%
%   record_average
%       true generates the Rec-level average after per-cycle processing.
%       The output type follows the recording stimulus type; flash and
%       grating results are not mixed. record_average_only also forces this
%       option true.
%
%   stop_on_error
%       true stops the batch at the first failed Cycle. false records the
%       error, saves a per-cycle error log, and continues with later cycles.
%
%   reuse_saved_peaks
%       analysis_only only. true validates and reuses accepted voltage peaks.
%       false sets peak="run", so peak finding/manual editing runs separately
%       for every Cycle before stimulus analysis.
%
%   write_diary
%       true saves the batch command-window output as a timestamped log under
%       the Rec folder. false prints only to the current command window.
%
%   run_motion
%       Full mode request passed to Dual3: true runs motion correction, false
%       skips it. analysis_only forces false; motion_only forces true.
%
%   cycle_offset_mode
%       Requested inter-camera registration mode for non-reference full-mode
%       Cycle runs: "none", "manual_points", or "matlab_register". Dual3
%       forces it to "none" whenever that cycle's registration section skips.
%
%   reuse_offset
%       Existing [x y] channel offset passed to reference and per-cycle Dual3
%       runs when no new registration offset is estimated.
%
% Common examples:
%
%   1. Full analysis of all Cycles using Cycle1 as the ROI reference:
%      dual3_rec_control = struct( ...
%          'rec_path', 'D:\Data\Rec1_2026-01-01_12-00-00', ...
%          'preset', "full", ...
%          'reference_cycle_name', "Cycle1", ...
%          'record_average', true);
%
%   2. Full analysis using an explicitly selected ROI folder:
%      dual3_rec_control = struct( ...
%          'rec_path', 'D:\Data\Rec1_2026-01-01_12-00-00', ...
%          'preset', "full", ...
%          'reference_roi_path', ...
%              'D:\Data\Rec1_2026-01-01_12-00-00\Cycle1\Dual_analysis3\previous_run');
%
%   3. Analysis-only for every Cycle, reusing saved accepted peaks:
%      dual3_rec_control = struct( ...
%          'rec_path', 'D:\Data\Rec1_2026-01-01_12-00-00', ...
%          'preset', "analysis_only", ...
%          'reuse_saved_peaks', true, ...
%          'record_average', true);
%
%   4. Analysis-only with peak finding/manual editing repeated per Cycle:
%      dual3_rec_control = struct( ...
%          'rec_path', 'D:\Data\Rec1_2026-01-01_12-00-00', ...
%          'preset', "analysis_only", ...
%          'reuse_saved_peaks', false);
%
%   5. Process selected Cycles only and continue after a failed Cycle:
%      dual3_rec_control = struct( ...
%          'rec_path', 'D:\Data\Rec1_2026-01-01_12-00-00', ...
%          'preset', "full", ...
%          'cycles', ["Cycle2","Cycle4"], ...
%          'stop_on_error', false);
%
%   6. Motion correction for all Rec*/Cycle* folders under one Methods folder:
%      dual3_rec_control = struct( ...
%          'rec_path', 'D:\Data\Methods1', ...
%          'preset', "motion_only", ...
%          'skip_existing', true);
%
%   7. Generate only the Rec-level average from existing cycle results:
%      dual3_rec_control = struct( ...
%          'rec_path', 'D:\Data\Rec1_2026-01-01_12-00-00', ...
%          'preset', "analysis_only", ...
%          'record_average_only', true);
%
% Run after defining the control:
%      run('run_Dual_analysis3_rec.m');
if ~exist('dual3_rec_control', 'var') || isempty(dual3_rec_control)
    dual3_rec_control = struct();
end
dual3_rec_control_defaults = struct( ...
    'rec_path', '', ...
    'preset', "", ...
    'mode', "", ...
    'sections', struct(), ...
    'source_results_path', "", ...
    'roi_path', "", ...
    'output_root', "Dual_analysis3", ...
    'run_tag', "", ...
    'reference_cycle_name', "Cycle1", ...
    'reference_roi_path', "", ...
    'reuse_reference_roi', true, ...
    'redraw_reference_roi', false, ...
    'reference_offset_mode', "none", ...
    'include_reference_cycle', false, ...
    'cycles', strings(0, 1), ...
    'skip_existing', [], ...
    'record_average_only', false, ...
    'record_average', false, ...
    'stop_on_error', true, ...
    'reuse_saved_peaks', true, ...
    'write_diary', true, ...
    'run_motion', true, ...
    'cycle_offset_mode', "none", ...
    'reuse_offset', [0, 0]);
dual3_rec_control_effective = apply_control_defaults_rec(dual3_rec_control, dual3_rec_control_defaults);

rec_path = char(string(dual3_rec_control_effective.rec_path));
requested_rec_preset = strtrim(string(dual3_rec_control_effective.preset));
legacy_rec_mode = strtrim(string(dual3_rec_control_effective.mode));
if strlength(requested_rec_preset) == 0 && strlength(legacy_rec_mode) > 0
    requested_rec_preset = legacy_rec_mode;
    warning('run_Dual_analysis3_rec:LegacyModeField', ...
        'dual3_rec_control.mode is deprecated; use dual3_rec_control.preset.');
elseif strlength(requested_rec_preset) > 0 && strlength(legacy_rec_mode) > 0 ...
        && ~strcmpi(requested_rec_preset, legacy_rec_mode)
    error(['dual3_rec_control.preset and legacy dual3_rec_control.mode disagree. ' ...
        'Set only preset, or make both values identical.']);
elseif strlength(requested_rec_preset) == 0
    requested_rec_preset = "full";
end
[preset_name, preset_sections, preset_overlay_sections, resolved_sections, preset_resolution] = ...
    resolve_dual3_sections(requested_rec_preset, dual3_rec_control_effective.sections);
requested_preset_before_canonicalization = preset_resolution.requested_preset;
if preset_name == "full" && ~isfield(preset_overlay_sections, 'motion')
    if logical(dual3_rec_control_effective.run_motion)
        preset_overlay_sections.motion = "run";
    else
        preset_overlay_sections.motion = "skip";
    end
elseif preset_name == "analysis_only" && ~isfield(preset_overlay_sections, 'peak')
    if logical(dual3_rec_control_effective.reuse_saved_peaks)
        preset_overlay_sections.peak = "reuse";
    else
        preset_overlay_sections.peak = "run";
    end
end
[preset_name, preset_sections, preset_overlay_sections, resolved_sections] = ...
    resolve_dual3_sections(preset_name, preset_overlay_sections);
analysis_mode = char(preset_name); % Saved for compatibility with older batch summaries.
dual3_rec_control_effective.requested_preset = requested_preset_before_canonicalization;
dual3_rec_control_effective.preset = preset_name;
dual3_rec_control_effective.mode = preset_name;
dual3_rec_control_effective.preset_sections = preset_sections;
dual3_rec_control_effective.overlay_sections = preset_overlay_sections;
dual3_rec_control_effective.resolved_sections = resolved_sections;
analysis_only_rec_mode = preset_name == "analysis_only";
motion_only_rec_mode = preset_name == "motion_only";
direct_preset_rec_mode = any(preset_name == ["reuse_motion", "reuse_trace", "skip_motion", "custom"]);
analysis_only_output_root_dir = char(string(dual3_rec_control_effective.output_root));
analysis_only_run_tag = string(dual3_rec_control_effective.run_tag);
reference_cycle_name = char(string(dual3_rec_control_effective.reference_cycle_name));
reference_roi_path = char(string(dual3_rec_control_effective.reference_roi_path));
direct_source_results_path = string(dual3_rec_control_effective.source_results_path);
direct_roi_path = string(dual3_rec_control_effective.roi_path);
direct_execution_plan = resolved_sections;
reference_roi_file = "";
if ~motion_only_rec_mode && ~analysis_only_rec_mode ...
        && strlength(string(reference_roi_path)) > 0
    if ~isfolder(reference_roi_path)
        error('reference_roi_path folder does not exist: %s', reference_roi_path);
    end
    reference_roi_file = resolve_roi_file_from_folder_rec(reference_roi_path);
end
reuse_existing_reference_roi = logical(dual3_rec_control_effective.reuse_reference_roi);
rerun_reference_cycle_for_roi = logical(dual3_rec_control_effective.redraw_reference_roi);
direct_reference_roi_enabled = direct_preset_rec_mode && ( ...
    reuse_existing_reference_roi ...
    || rerun_reference_cycle_for_roi ...
    || strlength(string(reference_roi_path)) > 0 ...
    || strlength(direct_roi_path) > 0);
reference_correct_offset_mode = dual3_rec_control_effective.reference_offset_mode;
run_reference_cycle_in_batch = logical(dual3_rec_control_effective.include_reference_cycle);
cycle_name_filter = string(dual3_rec_control_effective.cycles);
skip_cycles_with_existing_results = dual3_rec_control_effective.skip_existing;
run_record_average_only = logical(dual3_rec_control_effective.record_average_only);
run_record_average = logical(dual3_rec_control_effective.record_average);
stop_on_cycle_error = logical(dual3_rec_control_effective.stop_on_error);
analysis_only_reuse_saved_peak_results = resolved_sections.peak == "reuse";
enable_batch_diary = logical(dual3_rec_control_effective.write_diary);
run_motion_correction = ismember(resolved_sections.motion, ["run", "reuse"]);
correct_offset_mode = dual3_rec_control_effective.cycle_offset_mode;
reuse_offset = dual3_rec_control_effective.reuse_offset;
source_result_root_dir = 'Dual_analysis3';
result_root_dir = source_result_root_dir;
average_output_dir_name = 'Dual_analysis3_rec_average';
batch_summary_file_name = 'Dual_analysis3_rec_batch_summary.mat';
if analysis_only_rec_mode && strlength(analysis_only_run_tag) == 0
    analysis_only_run_tag = string(datetime('now', 'Format', 'yyyy-MM-dd HH-mm-ss'));
end
dual3_rec_control_effective.run_tag = analysis_only_run_tag;
if analysis_only_rec_mode
    run_motion_correction = false;
    correct_offset_mode = "none";
elseif motion_only_rec_mode
    run_motion_correction = true;
    correct_offset_mode = "none";
end
dual3_rec_control_effective.run_motion = run_motion_correction;
dual3_rec_control_effective.reuse_saved_peaks = analysis_only_reuse_saved_peak_results;

% Whether the reference cycle should enter manual inter-camera ROI offset
% correction before ROI drawing.
% Modes:
%   'none'            -> keep reuse_offset / [0 0]
%   'manual_points'   -> click one matching point on the two average images
%   'matlab_register' -> use MATLAB built-in translation registration
% Legacy true/false is still accepted and maps to manual_points/none.
reference_correct_offset_mode = normalize_correct_offset_mode_rec(reference_correct_offset_mode);
correct_offset_mode = normalize_correct_offset_mode_rec(correct_offset_mode);
dual3_rec_control_effective.cycle_offset_mode = string(correct_offset_mode);

% Skip a cycle if it already contains previous matching results. In
% motion-only mode this checks completed motion-correction outputs; in the
% normal workflow it checks completed Dual_analysis3 ROI results.
if isempty(skip_cycles_with_existing_results)
    skip_cycles_with_existing_results = motion_only_rec_mode;
end
dual3_rec_control_effective.skip_existing = logical(skip_cycles_with_existing_results);
run_record_average = logical(run_record_average) || logical(run_record_average_only);
fprintf('[Batch] Effective workflow control:\n');
fprintf('  preset=%s\n', string(preset_name));
fprintf('  rec_path=%s\n', string(rec_path));
fprintf('  cycle_filter_count=%d\n', numel(cycle_name_filter));
fprintf('  skip_existing=%d\n', logical(skip_cycles_with_existing_results));
fprintf('  record_average=%d\n', run_record_average);
fprintf('  run_motion=%d\n', run_motion_correction);
fprintf('  cycle_offset_mode=%s\n', string(correct_offset_mode));
if analysis_only_rec_mode
    fprintf('  analysis_only_output_root=%s\n', string(analysis_only_output_root_dir));
    fprintf('  reuse_saved_peaks=%d\n', analysis_only_reuse_saved_peak_results);
elseif direct_preset_rec_mode
    fprintf('  dual3_preset=%s\n', string(analysis_mode));
    fprintf('  dual3_source_results_path=%s\n', string(direct_source_results_path));
    fprintf('  dual3_roi_path=%s\n', string(direct_roi_path));
    fprintf('  shared_reference_roi=%d | reference_cycle=%s\n', ...
        direct_reference_roi_enabled, string(reference_cycle_name));
end
print_dual3_rec_execution_plan( ...
    analysis_mode, reuse_existing_reference_roi, rerun_reference_cycle_for_roi, ...
    strlength(string(reference_roi_file)) > 0, analysis_only_reuse_saved_peak_results, ...
    skip_cycles_with_existing_results, run_record_average, ...
    preset_sections, preset_overlay_sections, resolved_sections, direct_reference_roi_enabled);

% -------------------------------------------------------------------------
% B. Algorithm Parameters: common channel settings
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
% -------------------------------------------------------------------------
% C. Algorithm Parameters: shared Dual_analysis3 overrides
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
    run_background_removal =true;
end
if ~exist('voltage_polarity', 'var') || isempty(voltage_polarity)
    voltage_polarity = -1;
end
if ~exist('calcium_polarity', 'var') || isempty(calcium_polarity)
    calcium_polarity = 1;
end
if ~exist('voltage_peak_polarity_mode', 'var') || isempty(voltage_peak_polarity_mode)
    voltage_peak_polarity_mode = "global";
end
if ~exist('voltage_peak_min_prominence', 'var') || isempty(voltage_peak_min_prominence)
    voltage_peak_min_prominence = 0.5;
end
if ~exist('voltage_peak_min_prominence_mode', 'var') || isempty(voltage_peak_min_prominence_mode)
    voltage_peak_min_prominence_mode = "relative_factor";
end
if ~exist('voltage_peak_min_distance_frames', 'var') || isempty(voltage_peak_min_distance_frames)
    voltage_peak_min_distance_frames = max(2, round(0.01 * voltage_frame_rate));
end
if ~exist('voltage_peak_min_height', 'var') || isempty(voltage_peak_min_height)
    voltage_peak_min_height = 0;
end
if ~exist('run_manual_voltage_peak_edit', 'var') || isempty(run_manual_voltage_peak_edit)
    run_manual_voltage_peak_edit = true;
end
voltage_peak_params = struct( ...
    'polarity_mode', string(voltage_peak_polarity_mode), ...
    'min_peak_prominence', voltage_peak_min_prominence, ...
    'min_peak_prominence_mode', string(voltage_peak_min_prominence_mode), ...
    'min_peak_distance_frames', voltage_peak_min_distance_frames, ...
    'min_peak_height', voltage_peak_min_height, ...
    'run_manual_edit', logical(run_manual_voltage_peak_edit));
if ~exist('gpu', 'var') || isempty(gpu)
    gpu = true;
end
if ~exist('motion_cfg', 'var') || isempty(motion_cfg)
    motion_cfg = struct( ...
        'enabled', true, ...
        'use_saved_shift', false, ...
        'saved_shift_file', '', ...
        'saved_calcium_shift_file', '', ...
        'highpass', true, ...
        'auto_reuse_previous_shift', false);
end
motion_cfg.enabled = logical(run_motion_correction);

%% Resolve Cycles
if exist('dual_script_path', 'var') && ~isempty(dual_script_path)
    dual_script_path = char(string(dual_script_path));
else
    script_dir = fileparts(mfilename('fullpath'));
    dual_script_path = fullfile(script_dir, 'Dual_analysis3.m');
end
if ~isfile(dual_script_path)
    path_match = which('Dual_analysis3.m');
    if ~isempty(path_match)
        dual_script_path = path_match;
    end
end
if ~isfile(dual_script_path)
    error(['Cannot find Dual_analysis3.m. Put Dual_analysis3.m next to run_Dual_analysis3_rec.m, ' ...
        'add the AP_analysis code folder to the MATLAB path, or set dual_script_path explicitly.']);
end
fprintf('[Batch] Dual_analysis3 script:\n  %s\n', dual_script_path);

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

if motion_only_rec_mode
    %% Run Motion-Only Across Rec/Cycle Folders
    fprintf('\n============================================================\n');
    fprintf('[Batch] Motion-Only Across Cycles\n');
    fprintf('============================================================\n');
    motion_cycle_entries = resolve_motion_only_cycle_entries_rec(rec_path, cycle_name_filter);
    fprintf('Motion-only scope root: %s\n', rec_path);
    fprintf('Discovered cycles: %d\n', numel(motion_cycle_entries));

    batch_results = repmat(struct( ...
        'cycle_name', "", ...
        'cycle_path', "", ...
        'status', "", ...
        'save_path', "", ...
        'roi_file', "", ...
        'message', ""), 0, 1);

    for i = 1:numel(motion_cycle_entries)
        current_entry = motion_cycle_entries(i);
        current_cycle_label = current_entry.label;
        current_cycle_path = char(current_entry.cycle_path);

        if skip_cycles_with_existing_results
            existing_motion = find_latest_motion_correction_result(current_cycle_path);
            if strlength(existing_motion) > 0
                batch_results(end+1, 1) = struct( ...
                    'cycle_name', current_cycle_label, ...
                    'cycle_path', string(current_cycle_path), ...
                    'status', "skipped_existing_motion_correction", ...
                    'save_path', string(fileparts(existing_motion)), ...
                    'roi_file', "", ...
                    'message', existing_motion);
                fprintf('Skipping %s because motion-correction output already exists:\n  %s\n', current_cycle_label, existing_motion);
                continue;
            end
        end

        fprintf('\n[Batch] Motion-only for %s ...\n', current_cycle_label);
        try
            run_result = run_dual_cycle_motion_only( ...
                dual_script_path, current_cycle_path, preset_overlay_sections, ...
                camera_cfg, map_bin, calcium_smoothing_window, ...
                bleach_mode_voltage, bleach_mode_calcium, run_background_removal, ...
                voltage_polarity, calcium_polarity, ...
                reuse_offset, gpu, motion_cfg);

            batch_results(end+1, 1) = struct( ...
                'cycle_name', current_cycle_label, ...
                'cycle_path', string(current_cycle_path), ...
                'status', "completed", ...
                'save_path', string(run_result.save_path), ...
                'roi_file', "", ...
                'message', string(run_result.motion_result_file));
        catch ME
            error_report = getReport(ME, 'extended', 'hyperlinks', 'off');
            error_log_file = save_cycle_error_report_rec(rec_path, sanitize_path_component_rec(current_cycle_label), error_report);
            batch_results(end+1, 1) = struct( ...
                'cycle_name', current_cycle_label, ...
                'cycle_path', string(current_cycle_path), ...
                'status', "failed", ...
                'save_path', "", ...
                'roi_file', "", ...
                'message', string(error_report));
            fprintf(2, '[Batch] %s motion-only failed:\n%s\n', current_cycle_label, error_report);
            fprintf(2, '[Batch] Error report saved to:\n  %s\n', error_log_file);
            if stop_on_cycle_error
                rethrow(ME);
            end
        end
    end

    reference_roi_file = "";
    save(batch_summary_file, 'batch_results', 'reference_roi_file', 'analysis_mode', ...
        'motion_cycle_entries', 'dual3_rec_control_effective');
    fprintf('\nMotion-only batch summary saved to:\n  %s\n', batch_summary_file);
    return;
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
        existing_result = find_latest_explicit_result(current_cycle_path);
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
        reference_roi_file = string(find_latest_roi_file_rec(reference_cycle_path));
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
        existing_run_i = find_latest_explicit_result(cycle_path_i);
        if strlength(existing_run_i) > 0
            fprintf('%s | existing result marker detected:\n  %s\n', cycle_dirs(i).name, existing_run_i);
        elseif direct_preset_rec_mode && strlength(direct_source_results_path) == 0
            reusable_source_i = find_latest_reuse_source_rec( ...
                cycle_path_i, direct_execution_plan, direct_roi_path);
            if strlength(reusable_source_i) > 0
                fprintf('%s | reusable direct-preset source detected:\n  %s\n', ...
                    cycle_dirs(i).name, reusable_source_i);
            else
                fprintf('%s | no processed %s run detected.\n', cycle_dirs(i).name, source_result_root_dir);
            end
        else
            fprintf('%s | no processed %s run detected.\n', cycle_dirs(i).name, source_result_root_dir);
        end
    end

    if direct_preset_rec_mode
        %% Run Direct Dual3 Presets
        fprintf('\n============================================================\n');
        fprintf('[Batch] Run Dual3 Direct Preset\n');
        fprintf('============================================================\n');
        fprintf('Dual3 preset: %s\n', string(analysis_mode));

        batch_results = repmat(struct( ...
            'cycle_name', "", ...
            'cycle_path', "", ...
            'status', "", ...
            'save_path', "", ...
            'roi_file', "", ...
            'message', ""), 0, 1);

        direct_reference_result = struct('save_path', "", 'roi_file', "");
        direct_shared_roi_path = direct_roi_path;
        if direct_reference_roi_enabled
            reference_cycle_path = fullfile(rec_path, reference_cycle_name);
            if ~isfolder(reference_cycle_path)
                error('Reference cycle folder does not exist: %s', reference_cycle_path);
            end

            fprintf('\n============================================================\n');
            fprintf('[Batch] Resolve Direct-Preset Reference ROI\n');
            fprintf('============================================================\n');
            fprintf('Reference cycle: %s\n', reference_cycle_name);

            if strlength(string(reference_roi_file)) == 0 && strlength(direct_roi_path) > 0
                if ~isfolder(direct_roi_path)
                    error('dual3_rec_control.roi_path folder does not exist: %s', direct_roi_path);
                end
                reference_roi_file = resolve_roi_file_from_folder_rec(char(direct_roi_path));
                fprintf('Using explicit roi_path as the shared reference ROI:\n  %s\n', ...
                    reference_roi_file);
            end

            if strlength(string(reference_roi_file)) == 0 ...
                    && reuse_existing_reference_roi && ~rerun_reference_cycle_for_roi
                reference_roi_file = string(find_latest_roi_file_rec(reference_cycle_path));
                if strlength(reference_roi_file) > 0
                    fprintf('Reusing latest ROI file from %s:\n  %s\n', ...
                        reference_cycle_name, reference_roi_file);
                end
            end

            if strlength(string(reference_roi_file)) == 0 || rerun_reference_cycle_for_roi
                if direct_execution_plan.roi ~= "run"
                    error(['No reusable ROI was found in %s, but preset %s has roi="%s" and ' ...
                        'cannot create a reference ROI. Provide reference_roi_path/roi_path, ' ...
                        'or use a custom preset with roi="run".'], ...
                        reference_cycle_name, string(analysis_mode), direct_execution_plan.roi);
                end

                reference_source_results_path = direct_source_results_path;
                if strlength(reference_source_results_path) == 0
                    reference_source_results_path = find_latest_reuse_source_rec( ...
                        reference_cycle_path, direct_execution_plan, "");
                    if plan_requires_source_rec(direct_execution_plan, "")
                        if strlength(reference_source_results_path) == 0
                            error(['No result folder satisfying the %s reuse sections was found ' ...
                                'under reference cycle %s.'], string(analysis_mode), reference_cycle_path);
                        end
                        fprintf('[Batch] Auto-selected reference-Cycle reuse source:\n  %s\n', ...
                            reference_source_results_path);
                    end
                end

                fprintf('Running %s to create/recreate the reference ROI...\n', reference_cycle_name);
                direct_reference_result = run_dual_cycle_direct_preset( ...
                    dual_script_path, reference_cycle_path, string(analysis_mode), ...
                    reference_source_results_path, "", preset_overlay_sections, ...
                    reference_correct_offset_mode, ...
                    camera_cfg, map_bin, calcium_smoothing_window, ...
                    bleach_mode_voltage, bleach_mode_calcium, run_background_removal, ...
                    voltage_polarity, calcium_polarity, ...
                    reuse_offset, gpu, motion_cfg, voltage_peak_params);
                reference_roi_file = string(direct_reference_result.roi_file);
                fprintf('Reference ROI file created:\n  %s\n', reference_roi_file);
            end

            if strlength(string(reference_roi_file)) == 0 || ~isfile(reference_roi_file)
                error('Direct-preset reference ROI file could not be resolved: %s', reference_roi_file);
            end
            direct_shared_roi_path = string(fileparts(char(reference_roi_file)));
            fprintf('Shared ROI folder for non-reference Cycles:\n  %s\n', direct_shared_roi_path);
        end

        for i = 1:numel(cycle_dirs)
            current_cycle_name = string(cycle_dirs(i).name);
            current_cycle_path = fullfile(cycle_dirs(i).folder, cycle_dirs(i).name);

            if direct_reference_roi_enabled ...
                    && strcmpi(current_cycle_name, string(reference_cycle_name)) ...
                    && ~run_reference_cycle_in_batch
                batch_results(end+1, 1) = struct( ...
                    'cycle_name', current_cycle_name, ...
                    'cycle_path', string(current_cycle_path), ...
                    'status', "skipped_reference_already_used_for_roi", ...
                    'save_path', string(direct_reference_result.save_path), ...
                    'roi_file', reference_roi_file, ...
                    'message', "Reference cycle already used to create/reuse the shared ROI.");
                fprintf('Skipping %s in the direct-preset batch loop because it supplied the shared ROI.\n', ...
                    current_cycle_name);
                continue;
            end

            if skip_cycles_with_existing_results
                existing_result = find_latest_explicit_result(current_cycle_path);
                if strlength(existing_result) > 0
                    batch_results(end+1, 1) = struct( ...
                        'cycle_name', current_cycle_name, ...
                        'cycle_path', string(current_cycle_path), ...
                        'status', "skipped_existing_result", ...
                        'save_path', fileparts(existing_result), ...
                        'roi_file', ternary_rec(direct_reference_roi_enabled, reference_roi_file, ""), ...
                        'message', existing_result);
                    fprintf('Skipping %s because the final output already exists:\n  %s\n', current_cycle_name, existing_result);
                    continue;
                end
            end

            fprintf('\n[Batch] Dual3 preset %s for %s ...\n', string(analysis_mode), current_cycle_name);
            try
                cycle_execution_plan = direct_execution_plan;
                cycle_sections = preset_overlay_sections;
                cycle_roi_path = direct_roi_path;
                if direct_reference_roi_enabled
                    cycle_execution_plan.roi = "reuse";
                    cycle_sections.roi = "reuse";
                    cycle_roi_path = direct_shared_roi_path;
                    fprintf('[Batch] %s shared ROI source:\n  %s\n', ...
                        current_cycle_name, reference_roi_file);
                end

                cycle_source_results_path = direct_source_results_path;
                if strlength(cycle_source_results_path) == 0
                    cycle_source_results_path = find_latest_reuse_source_rec( ...
                        current_cycle_path, cycle_execution_plan, cycle_roi_path);
                    if plan_requires_source_rec(cycle_execution_plan, cycle_roi_path)
                        if strlength(cycle_source_results_path) == 0
                            error(['No result folder satisfying the %s reuse sections was found under %s. ' ...
                                'Set dual3_rec_control.source_results_path explicitly or provide a valid prior result.'], ...
                                string(analysis_mode), current_cycle_path);
                        end
                        fprintf('[Batch] Auto-selected reuse source for %s:\n  %s\n', ...
                            current_cycle_name, cycle_source_results_path);
                    end
                end
                run_result = run_dual_cycle_direct_preset( ...
                    dual_script_path, current_cycle_path, string(analysis_mode), ...
                    cycle_source_results_path, cycle_roi_path, cycle_sections, ...
                    correct_offset_mode, ...
                    camera_cfg, map_bin, calcium_smoothing_window, ...
                    bleach_mode_voltage, bleach_mode_calcium, run_background_removal, ...
                    voltage_polarity, calcium_polarity, ...
                    reuse_offset, gpu, motion_cfg, voltage_peak_params);

                batch_results(end+1, 1) = struct( ...
                    'cycle_name', current_cycle_name, ...
                    'cycle_path', string(current_cycle_path), ...
                    'status', "completed", ...
                    'save_path', string(run_result.save_path), ...
                    'roi_file', string(run_result.roi_file), ...
                    'message', "Dual3 preset " + string(analysis_mode));
            catch ME
                error_report = getReport(ME, 'extended', 'hyperlinks', 'off');
                error_log_file = save_cycle_error_report_rec(rec_path, current_cycle_name, error_report);
                batch_results(end+1, 1) = struct( ...
                    'cycle_name', current_cycle_name, ...
                    'cycle_path', string(current_cycle_path), ...
                    'status', "failed", ...
                    'save_path', "", ...
                    'roi_file', ternary_rec(direct_reference_roi_enabled, reference_roi_file, ""), ...
                    'message', string(error_report));
                fprintf(2, '[Batch] %s direct preset failed:\n%s\n', current_cycle_name, error_report);
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
    elseif analysis_only_rec_mode
        %% Run Analysis-Only Forks From Existing Cycle Results
        fprintf('\n============================================================\n');
        fprintf('[Batch] Run Analysis-Only Forks\n');
        fprintf('============================================================\n');
        fprintf('Source root: Dual_analysis3\n');
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
            source_result = find_latest_reusable_analysis_result(current_cycle_path);

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
                    analysis_only_output_root_dir, analysis_only_run_tag, preset_overlay_sections, ...
                    camera_cfg, map_bin, calcium_smoothing_window, ...
                    bleach_mode_voltage, bleach_mode_calcium, run_background_removal, ...
                    voltage_polarity, calcium_polarity, ...
                    reuse_offset, gpu, motion_cfg, voltage_peak_params);

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
        reference_existing_run = find_latest_explicit_result(reference_cycle_path);
        if strlength(reference_existing_run) > 0
            fprintf('Reference cycle already has a result marker:\n  %s\n', reference_existing_run);
        end

        if strlength(string(reference_roi_file)) == 0
            reference_roi_file = "";
        else
            reference_roi_file = string(reference_roi_file);
        end

        if strlength(reference_roi_file) == 0 && reuse_existing_reference_roi && ~rerun_reference_cycle_for_roi
            reference_roi_file = string(find_latest_roi_file_rec(reference_cycle_path));
            if strlength(reference_roi_file) > 0
                fprintf('Reusing latest ROI file already present in reference cycle:\n  %s\n', reference_roi_file);
            end
        end

        if strlength(reference_roi_file) == 0 || rerun_reference_cycle_for_roi
            fprintf('Running reference cycle to create/recreate ROI...\n');
            reference_result = run_dual_cycle_with_overrides( ...
                dual_script_path, reference_cycle_path, "", preset_overlay_sections, reference_correct_offset_mode, ...
                camera_cfg, map_bin, calcium_smoothing_window, ...
                bleach_mode_voltage, bleach_mode_calcium, run_background_removal, ...
                voltage_polarity, calcium_polarity, ...
                reuse_offset, gpu, motion_cfg, voltage_peak_params);
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
                    dual_script_path, current_cycle_path, fileparts(char(reference_roi_file)), ...
                    preset_overlay_sections, correct_offset_mode, ...
                    camera_cfg, map_bin, calcium_smoothing_window, ...
                    bleach_mode_voltage, bleach_mode_calcium, run_background_removal, ...
                    voltage_polarity, calcium_polarity, ...
                    reuse_offset, gpu, motion_cfg, voltage_peak_params);

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

    save(batch_summary_file, 'batch_results', 'reference_roi_file', 'analysis_mode', ...
        'analysis_only_output_root_dir', 'analysis_only_run_tag', ...
        'voltage_polarity', 'calcium_polarity', 'calcium_smoothing_window', ...
        'dual3_rec_control_effective');
    fprintf('\nBatch summary saved to:\n  %s\n', batch_summary_file);
end

%% Optional Record-Level Analysis
% Rec-level averaging is implemented in its own script. The saved batch
% summary is the interface, so the same exact Cycle result set can be
% analyzed again later without rerunning per-Cycle Dual_analysis3.
if run_record_average
    fprintf('\n============================================================\n');
    fprintf('[Batch] Calling Dual_rec_analysis3\n');
    fprintf('============================================================\n');
    dual_rec_script_path = fullfile(fileparts(mfilename('fullpath')), 'Dual_rec_analysis3.m');
    if ~isfile(dual_rec_script_path)
        error('Cannot find Dual_rec_analysis3.m next to run_Dual_analysis3_rec.m.');
    end
    dual_rec_analysis3_control = struct( ...
        'batch_summary_path', string(batch_summary_file), ...
        'output_dir_name', string(average_output_dir_name));
    run(dual_rec_script_path);
else
    fprintf('\n[Batch] Rec-level analysis skipped. Set record_average=true to call Dual_rec_analysis3.\n');
end

%% Local Functions
function control = apply_control_defaults_rec(control, defaults)
if ~isstruct(control) || ~isscalar(control)
    error('dual3_rec_control must be one scalar struct.');
end
unknown_fields = setdiff(fieldnames(control), fieldnames(defaults));
if ~isempty(unknown_fields)
    error('Unknown dual3_rec_control field(s): %s', strjoin(unknown_fields, ', '));
end
default_fields = fieldnames(defaults);
for idx = 1:numel(default_fields)
    field_name = default_fields{idx};
    if ~isfield(control, field_name) || isempty(control.(field_name))
        control.(field_name) = defaults.(field_name);
    end
end
end

function print_dual3_rec_execution_plan( ...
    preset_name, reuse_reference_roi, redraw_reference_roi, has_reference_roi_file, ...
    reuse_saved_peaks, skip_existing, run_average, ...
    preset_sections, overlay_sections, resolved_sections, direct_reference_roi_enabled)
preset_name = lower(string(preset_name));
rows = strings(0, 3);
rows(end+1, :) = ["Discover Cycle folders", "RUN", "Apply control.cycles filter"];
switch preset_name
    case "full"
        if redraw_reference_roi
            rows(end+1, :) = ["Reference ROI", "RUN", "Force a new ROI on the reference cycle"];
        elseif reuse_reference_roi || has_reference_roi_file
            rows(end+1, :) = ["Reference ROI", "REUSE", "Use the latest/reference configured ROI"];
        else
            rows(end+1, :) = ["Reference ROI", "RUN", "Create ROI when no explicit reference ROI is supplied"];
        end
        rows(end+1, :) = ["Per-cycle movie processing", "RUN", "Run complete Dual_analysis3 for selected cycles"];
        rows(end+1, :) = ["Previous traces and peaks", "SKIP", "Full mode computes new results"];
    case {"reuse_motion", "reuse_trace", "skip_motion", "custom"}
        if direct_reference_roi_enabled
            rows(end+1, :) = ["Reference ROI orchestration", "RUN/REUSE", "Resolve one shared ROI and overlay roi=reuse per Cycle"];
        else
            rows(end+1, :) = ["Reference ROI orchestration", "SKIP", "Pass the requested Dual3 preset directly per Cycle"];
        end
        rows(end+1, :) = ["Per-cycle Dual3 preset", "RUN", "Call Dual_analysis3 with preset=" + preset_name];
        rows(end+1, :) = ["Dual3 sections overrides", "PASS", "Forward dual3_rec_control.sections to each Cycle"];
    case "analysis_only"
        rows(end+1, :) = ["Find previous result per Cycle", "RUN", "Select latest complete reusable result"];
        rows(end+1, :) = ["Per-cycle movie/ROI/trace processing", "SKIP", "Read saved ROI and traces"];
        if reuse_saved_peaks
            rows(end+1, :) = ["Per-cycle accepted peaks", "REUSE", "Validate accepted peaks from each source result"];
        else
            rows(end+1, :) = ["Per-cycle accepted peaks", "RUN", "Recompute peaks from saved voltage sensitivity traces"];
        end
        rows(end+1, :) = ["Per-cycle plots / DSI / OSI", "RUN", "Write new analysis_only result folder"];
        rows(end+1, :) = ["Previous result folders", "READ ONLY", "Never overwrite source folders"];
    case "motion_only"
        rows(end+1, :) = ["Per-cycle motion correction", "RUN", "Run motion stage only"];
        rows(end+1, :) = ["ROI / traces / peaks / tuning", "STOP", "Do not run later sections"];
end
if skip_existing
    rows(end+1, :) = ["Existing matching output", "SKIP", "Do not rerun a completed matching output"];
else
    rows(end+1, :) = ["Existing matching output", "RUN AGAIN", "Create another timestamped output"];
end
if run_average
    rows(end+1, :) = ["Record-level average", "RUN", "Generate stimulus-type-specific Rec summary"];
else
    rows(end+1, :) = ["Record-level average", "SKIP", "Keep per-cycle outputs only"];
end

fprintf('[Batch] Execution plan:\n');
fprintf('  %-38s %-10s %s\n', 'Stage', 'Action', 'Meaning');
for idx = 1:size(rows, 1)
    fprintf('  %-38s %-10s %s\n', rows(idx, 1), rows(idx, 2), rows(idx, 3));
end
fprintf('[Batch] Base section resolution before per-Cycle ROI orchestration:\n');
fprintf('  %-18s %-10s %-10s %-10s\n', 'Section', 'Preset', 'Overlay', 'Resolved');
section_names = fieldnames(resolved_sections);
for idx = 1:numel(section_names)
    section_name = section_names{idx};
    fprintf('  %-18s %-10s %-10s %-10s\n', section_name, ...
        get_section_action_rec(preset_sections, section_name, "CUSTOM"), ...
        get_section_action_rec(overlay_sections, section_name, "-"), ...
        upper(string(resolved_sections.(section_name))));
end
end

function action = get_section_action_rec(sections, section_name, missing_value)
if isstruct(sections) && isfield(sections, section_name)
    action = upper(string(sections.(section_name)));
else
    action = string(missing_value);
end
end

function result = run_dual_cycle_with_overrides( ...
    dual_script_path, cycle_path, reuse_roi_path, section_overlays, correct_offset_mode, ...
    camera_cfg, map_bin, calcium_smoothing_window, ...
    bleach_mode_voltage, bleach_mode_calcium, run_background_removal, ...
    voltage_polarity, calcium_polarity, ...
    reuse_offset, gpu, motion_cfg, voltage_peak_params)

% All variables defined in this function are visible to Dual_analysis3.m
% when it is executed via run(...). This keeps Dual_analysis3 as a script
% while still letting an outer batch driver provide overrides.
[record_path_for_name, cycle_name_for_name] = fileparts(cycle_path);
[~, record_name_for_name] = fileparts(record_path_for_name);
time_tag = string(datetime('now', 'Format', 'yyyy-MM-dd HH-mm-ss-SSS'));
analysis_run_name = sprintf('%s_%s_%s', record_name_for_name, cycle_name_for_name, char(time_tag));
if ~isstruct(section_overlays) || ~isscalar(section_overlays)
    error('dual3_rec_control.sections must be one scalar struct.');
end
cycle_sections = section_overlays;
if ~isfield(cycle_sections, 'motion')
    cycle_sections.motion = ternary_rec(logical(motion_cfg.enabled), "run", "skip");
end
if strlength(string(reuse_roi_path)) > 0
    cycle_sections.registration = "skip";
    cycle_sections.roi = "reuse";
else
    if ~isfield(cycle_sections, 'registration')
        cycle_sections.registration = "run";
    end
    if ~isfield(cycle_sections, 'roi')
        cycle_sections.roi = "run";
    end
end
dual3_control = struct( ...
    'cycle_path', string(cycle_path), ...
    'preset', "full", ...
    'run_name', string(analysis_run_name), ...
    'source_results_path', "", ...
    'roi_path', string(reuse_roi_path), ...
    'offset_mode', string(correct_offset_mode), ...
    'reuse_offset', reuse_offset, ...
    'sections', cycle_sections);
voltage_peak_polarity_mode = voltage_peak_params.polarity_mode;
voltage_peak_min_prominence = voltage_peak_params.min_peak_prominence;
voltage_peak_min_prominence_mode = voltage_peak_params.min_peak_prominence_mode;
voltage_peak_min_distance_frames = voltage_peak_params.min_peak_distance_frames;
voltage_peak_min_height = voltage_peak_params.min_peak_height;
run_manual_voltage_peak_edit = voltage_peak_params.run_manual_edit;

run(dual_script_path);

result = struct( ...
    'save_path', save_path, ...
    'roi_file', fullfile(save_path, '1_dual_roi_results.mat'));
end

function result = run_dual_cycle_direct_preset( ...
    dual_script_path, cycle_path, dual3_preset, source_results_path, roi_path, sections, correct_offset_mode, ...
    camera_cfg, map_bin, calcium_smoothing_window, ...
    bleach_mode_voltage, bleach_mode_calcium, run_background_removal, ...
    voltage_polarity, calcium_polarity, ...
    reuse_offset, gpu, motion_cfg, voltage_peak_params)

% Direct preset mode: the Rec runner only iterates Cycle folders and passes
% the requested Dual_analysis3 preset/control fields through.
[record_path_for_name, cycle_name_for_name] = fileparts(cycle_path);
[~, record_name_for_name] = fileparts(record_path_for_name);
time_tag = string(datetime('now', 'Format', 'yyyy-MM-dd HH-mm-ss-SSS'));
analysis_run_name = sprintf('%s_%s_%s', record_name_for_name, cycle_name_for_name, char(time_tag));
analysis_run_name = sanitize_path_component_rec(analysis_run_name);
if ~isstruct(sections)
    error('dual3_rec_control.sections must be a struct when provided.');
end
dual3_control = struct( ...
    'cycle_path', string(cycle_path), ...
    'preset', string(dual3_preset), ...
    'run_name', string(analysis_run_name), ...
    'source_results_path', string(source_results_path), ...
    'roi_path', string(roi_path), ...
    'offset_mode', string(correct_offset_mode), ...
    'reuse_offset', reuse_offset, ...
    'sections', sections);
voltage_peak_polarity_mode = voltage_peak_params.polarity_mode;
voltage_peak_min_prominence = voltage_peak_params.min_peak_prominence;
voltage_peak_min_prominence_mode = voltage_peak_params.min_peak_prominence_mode;
voltage_peak_min_distance_frames = voltage_peak_params.min_peak_distance_frames;
voltage_peak_min_height = voltage_peak_params.min_peak_height;
run_manual_voltage_peak_edit = voltage_peak_params.run_manual_edit;

run(dual_script_path);

result = struct( ...
    'save_path', save_path, ...
    'roi_file', fullfile(save_path, '1_dual_roi_results.mat'));
end

function result = run_dual_cycle_motion_only( ...
    dual_script_path, cycle_path, section_overlays, ...
    camera_cfg, map_bin, calcium_smoothing_window, ...
    bleach_mode_voltage, bleach_mode_calcium, run_background_removal, ...
    voltage_polarity, calcium_polarity, ...
    reuse_offset, gpu, motion_cfg)

% Run Dual_analysis3 only through shared motion correction. The script
% returns before channel registration, ROI selection, trace extraction, and
% comparison sections.
[record_path_for_name, cycle_name_for_name] = fileparts(cycle_path);
[~, record_name_for_name] = fileparts(record_path_for_name);
time_tag = string(datetime('now', 'Format', 'yyyy-MM-dd HH-mm-ss-SSS'));
analysis_run_name = sprintf('%s_%s_%s', ...
    record_name_for_name, cycle_name_for_name, char(time_tag));
analysis_run_name = sanitize_path_component_rec(analysis_run_name);
dual3_control = struct( ...
    'cycle_path', string(cycle_path), ...
    'preset', "motion_only", ...
    'run_name', string(analysis_run_name), ...
    'offset_mode', "none", ...
    'reuse_offset', reuse_offset, ...
    'sections', section_overlays);
motion_cfg.enabled = true;

run(dual_script_path);

result = struct( ...
    'save_path', save_path, ...
    'motion_result_file', fullfile(save_path, 'shared_motion_shifts_result.mat'));
end

function result = run_dual_cycle_analysis_only_fork( ...
    dual_script_path, cycle_path, source_results_path, output_root_dir, run_tag, section_overlays, ...
    camera_cfg, map_bin, calcium_smoothing_window, ...
    bleach_mode_voltage, bleach_mode_calcium, run_background_removal, ...
    voltage_polarity, calcium_polarity, ...
    reuse_offset, gpu, motion_cfg, voltage_peak_params)

% Reuse saved ROI/trace/stim context from an existing Dual_analysis3 folder,
% but write the rerun products to a new folder under output_root_dir.
[record_path_for_name, cycle_name_for_name] = fileparts(cycle_path);
[~, record_name_for_name] = fileparts(record_path_for_name);
analysis_run_name = sprintf('%s_%s_%s', ...
    record_name_for_name, cycle_name_for_name, char(string(run_tag)));
analysis_run_name = sanitize_path_component_rec(analysis_run_name);
save_path = fullfile(cycle_path, char(string(output_root_dir)), analysis_run_name);
dual3_control = struct( ...
    'cycle_path', string(cycle_path), ...
    'preset', "analysis_only", ...
    'run_name', string(analysis_run_name), ...
    'output_path', string(save_path), ...
    'source_results_path', string(source_results_path), ...
    'offset_mode', "none", ...
    'reuse_offset', reuse_offset, ...
    'sections', section_overlays);
voltage_peak_polarity_mode = voltage_peak_params.polarity_mode;
voltage_peak_min_prominence = voltage_peak_params.min_peak_prominence;
voltage_peak_min_prominence_mode = voltage_peak_params.min_peak_prominence_mode;
voltage_peak_min_distance_frames = voltage_peak_params.min_peak_distance_frames;
voltage_peak_min_height = voltage_peak_params.min_peak_height;
run_manual_voltage_peak_edit = voltage_peak_params.run_manual_edit;

run(dual_script_path);

result = struct( ...
    'save_path', save_path, ...
    'roi_file', fullfile(save_path, '1_dual_roi_results.mat'), ...
    'source_save_path', string(source_results_path));
end

function entries = resolve_motion_only_cycle_entries_rec(scope_path, cycle_name_filter)
% Accept either one Rec* folder containing Cycle* children or one Methods*
% folder containing Rec*/Cycle* children.
scope_path = char(string(scope_path));
cycle_name_filter = string(cycle_name_filter(:));
entries = repmat(struct('rec_name', "", 'cycle_name', "", 'label', "", 'cycle_path', ""), 0, 1);

direct_cycles = dir(fullfile(scope_path, 'Cycle*'));
direct_cycles = direct_cycles([direct_cycles.isdir]);
direct_cycles = sort_cycle_dirs(direct_cycles);
if ~isempty(direct_cycles)
    [~, rec_name] = fileparts(scope_path);
    for idx = 1:numel(direct_cycles)
        cycle_name = string(direct_cycles(idx).name);
        label = cycle_name;
        if should_keep_motion_cycle_rec(cycle_name, label, cycle_name_filter)
            entries(end+1, 1) = struct( ...
                'rec_name', string(rec_name), ...
                'cycle_name', cycle_name, ...
                'label', label, ...
                'cycle_path', string(fullfile(direct_cycles(idx).folder, direct_cycles(idx).name))); %#ok<AGROW>
        end
    end
    if isempty(entries)
        error('No direct Cycle* folders remain after applying cycle_name_filter.');
    end
    return;
end

rec_dirs = dir(fullfile(scope_path, 'Rec*'));
rec_dirs = rec_dirs([rec_dirs.isdir]);
if isempty(rec_dirs)
    error('No Cycle* folders or Rec*/Cycle* folders found under: %s', scope_path);
end
[~, rec_order] = sort(string({rec_dirs.name}));
rec_dirs = rec_dirs(rec_order);
for rec_idx = 1:numel(rec_dirs)
    rec_name = string(rec_dirs(rec_idx).name);
    cycle_dirs = dir(fullfile(rec_dirs(rec_idx).folder, rec_dirs(rec_idx).name, 'Cycle*'));
    cycle_dirs = cycle_dirs([cycle_dirs.isdir]);
    cycle_dirs = sort_cycle_dirs(cycle_dirs);
    for cycle_idx = 1:numel(cycle_dirs)
        cycle_name = string(cycle_dirs(cycle_idx).name);
        label = rec_name + "/" + cycle_name;
        if should_keep_motion_cycle_rec(cycle_name, label, cycle_name_filter)
            entries(end+1, 1) = struct( ...
                'rec_name', rec_name, ...
                'cycle_name', cycle_name, ...
                'label', label, ...
                'cycle_path', string(fullfile(cycle_dirs(cycle_idx).folder, cycle_dirs(cycle_idx).name))); %#ok<AGROW>
        end
    end
end

if isempty(entries)
    error('No Rec*/Cycle* folders remain after applying cycle_name_filter.');
end
end

function tf = should_keep_motion_cycle_rec(cycle_name, label, cycle_name_filter)
if isempty(cycle_name_filter)
    tf = true;
else
    tf = ismember(string(cycle_name), cycle_name_filter) || ismember(string(label), cycle_name_filter);
end
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

function roi_file = find_latest_roi_file_rec(cycle_path)
roi_file = "";
search_patterns = { ...
    fullfile(cycle_path, '**', '1_dual_roi_results.mat'), ...
    fullfile(cycle_path, '**', '1_raw_ROI.mat'), ...
    fullfile(cycle_path, '**', 'roi.mat'), ...
    fullfile(cycle_path, '**', '*ROI*.mat')};
listing = struct([]);
for pattern_idx = 1:numel(search_patterns)
    matches = dir(search_patterns{pattern_idx});
    listing = [listing; matches(~[matches.isdir])]; %#ok<AGROW>
end
if isempty(listing)
    return;
end
candidate_paths = string(arrayfun( ...
    @(entry) fullfile(entry.folder, entry.name), listing, 'UniformOutput', false));
[~, unique_idx] = unique(lower(candidate_paths), 'stable');
listing = listing(unique_idx);
[~, newest_order] = sort([listing.datenum], 'descend');
for candidate_idx = newest_order(:)'
    candidate = fullfile(listing(candidate_idx).folder, listing(candidate_idx).name);
    if validate_roi_mask_file_rec(candidate)
        roi_file = string(candidate);
        return;
    end
end
end

function roi_file = resolve_roi_file_from_folder_rec(roi_folder)
roi_folder = char(string(roi_folder));
listing = dir(fullfile(roi_folder, '*.mat'));
if isempty(listing)
    error('reference_roi_path contains no MAT files: %s', roi_folder);
end

priority_names = ["1_dual_roi_results.mat", "1_raw_ROI.mat", "roi.mat"];
listing_names = string({listing.name});
for priority_idx = 1:numel(priority_names)
    match_idx = find(strcmpi(listing_names, priority_names(priority_idx)), 1, 'first');
    if isempty(match_idx)
        continue;
    end
    candidate = fullfile(listing(match_idx).folder, listing(match_idx).name);
    [valid_roi, reason] = validate_roi_mask_file_rec(candidate);
    if ~valid_roi
        error('ROI file %s is invalid: %s', candidate, reason);
    end
    roi_file = string(candidate);
    return;
end

generic_idx = find(contains(lower(listing_names), 'roi'));
valid_files = strings(0, 1);
for candidate_idx = generic_idx(:)'
    candidate = fullfile(listing(candidate_idx).folder, listing(candidate_idx).name);
    if validate_roi_mask_file_rec(candidate)
        valid_files(end+1, 1) = string(candidate); %#ok<AGROW>
    end
end
if isempty(valid_files)
    error(['No valid ROI MAT file was found directly inside reference_roi_path: %s. ' ...
        'Expected rois.bwmask, bwmask, or mask.'], roi_folder);
end
if numel(valid_files) > 1
    error('Multiple valid ROI MAT files were found in reference_roi_path:\n  %s', ...
        strjoin(valid_files, newline + "  "));
end
roi_file = valid_files(1);
end

function [is_valid, reason] = validate_roi_mask_file_rec(roi_file)
is_valid = false;
reason = "";
try
    roi_data = load(roi_file);
catch ME
    reason = "MAT file could not be loaded: " + string(ME.message);
    return;
end
if isfield(roi_data, 'rois') && isstruct(roi_data.rois) ...
        && isfield(roi_data.rois, 'bwmask') && ~isempty(roi_data.rois.bwmask)
    is_valid = true;
elseif isfield(roi_data, 'bwmask') && ~isempty(roi_data.bwmask)
    is_valid = true;
elseif isfield(roi_data, 'mask') && ~isempty(roi_data.mask)
    is_valid = true;
else
    reason = "missing nonempty rois.bwmask, bwmask, or mask";
end
end

function explicit_result = find_latest_explicit_result(cycle_path)
explicit_result = "";
result_root_dir = 'Dual_analysis3';

% A completed result is identified only by the final explicit bundle. ROI-only
% folders remain available to reference-ROI discovery but do not trigger
% skip_existing for a full/direct analysis run.
explicit_listing = dir(fullfile(cycle_path, result_root_dir, '**', '-1_explicit_dual_results.mat'));
explicit_result = select_latest_file_from_listing_rec(explicit_listing);
end

function explicit_result = find_latest_reusable_analysis_result(cycle_path)
explicit_result = "";
result_root_dir = 'Dual_analysis3';
listing = dir(fullfile(cycle_path, result_root_dir, '**', 'voltage_results.mat'));
if isempty(listing)
    return;
end

[~, order] = sort([listing.datenum], 'descend');
listing = listing(order);
required_files = {'dual_info.mat', 'calcium_results.mat', '1_dual_roi_results.mat'};
for idx = 1:numel(listing)
    candidate_dir = listing(idx).folder;
    if all(cellfun(@(name) isfile(fullfile(candidate_dir, name)), required_files))
        explicit_result = string(fullfile(candidate_dir, listing(idx).name));
        return;
    end
end
end

function requires_source = plan_requires_source_rec(execution_plan, roi_path)
source_plan = execution_plan;
if source_plan.roi == "reuse" && strlength(strtrim(string(roi_path))) > 0
    source_plan.roi = "skip";
end
requires_source = any(structfun(@(value) string(value) == "reuse", source_plan));
end

function source_path = find_latest_reuse_source_rec(cycle_path, execution_plan, roi_path)
source_path = "";
source_plan = execution_plan;
if source_plan.roi == "reuse" && strlength(strtrim(string(roi_path))) > 0
    source_plan.roi = "skip";
end
if ~any(structfun(@(value) string(value) == "reuse", source_plan))
    return;
end

result_root = fullfile(cycle_path, 'Dual_analysis3');
if ~isfolder(result_root)
    return;
end

if source_plan.motion == "reuse"
    anchor_listing = [ ...
        dir(fullfile(result_root, '**', 'shared_motion_shifts_result.mat')); ...
        dir(fullfile(result_root, '**', 'motion_shifts_result.mat'))];
elseif any([source_plan.input, source_plan.registration, source_plan.trace, source_plan.peak] == "reuse")
    anchor_listing = dir(fullfile(result_root, '**', 'voltage_results.mat'));
elseif source_plan.roi == "reuse"
    anchor_listing = dir(fullfile(result_root, '**', '1_dual_roi_results.mat'));
elseif any([source_plan.visualization, source_plan.comparison] == "reuse")
    anchor_listing = dir(fullfile(result_root, '**', 'dual_results.mat'));
elseif source_plan.stim == "reuse"
    anchor_listing = dir(fullfile(result_root, '**', 'stim_results.mat'));
elseif source_plan.frequency == "reuse"
    anchor_listing = dir(fullfile(result_root, '**', '8_time_frequency_results.mat'));
else
    return;
end
if isempty(anchor_listing)
    return;
end

[~, order] = sort([anchor_listing.datenum], 'descend');
anchor_listing = anchor_listing(order);
visited_folders = strings(0, 1);
for idx = 1:numel(anchor_listing)
    candidate = string(anchor_listing(idx).folder);
    if any(visited_folders == candidate)
        continue;
    end
    visited_folders(end+1, 1) = candidate; %#ok<AGROW>
    if reuse_source_has_required_files_rec(char(candidate), source_plan)
        source_path = candidate;
        return;
    end
end
end

function is_valid = reuse_source_has_required_files_rec(candidate, source_plan)
required_files = strings(0, 1);
if source_plan.input == "reuse"
    required_files(end+1:end+3, 1) = ["dual_info.mat"; "voltage_results.mat"; "calcium_results.mat"];
end
if any([source_plan.registration, source_plan.trace, source_plan.peak] == "reuse")
    required_files(end+1, 1) = "voltage_results.mat";
end
if source_plan.trace == "reuse"
    required_files(end+1, 1) = "calcium_results.mat";
end
if source_plan.roi == "reuse"
    required_files(end+1:end+2, 1) = ["1_dual_roi_results.mat"; "dual_results.mat"];
end
if any([source_plan.visualization, source_plan.comparison] == "reuse")
    required_files(end+1, 1) = "dual_results.mat";
end
if source_plan.stim == "reuse"
    required_files(end+1, 1) = "stim_results.mat";
end
if source_plan.frequency == "reuse"
    required_files(end+1, 1) = "8_time_frequency_results.mat";
end
required_files = unique(required_files, 'stable');
is_valid = all(arrayfun(@(name) isfile(fullfile(candidate, name)), required_files));
if is_valid && source_plan.motion == "reuse"
    is_valid = isfile(fullfile(candidate, 'shared_motion_shifts_result.mat')) || ...
        isfile(fullfile(candidate, 'motion_shifts_result.mat'));
end
end

function motion_result = find_latest_motion_correction_result(cycle_path)
motion_result = "";
motion_roots = "Dual_analysis3";
candidate_listing = [];
for root_idx = 1:numel(motion_roots)
    current_listing = dir(fullfile(cycle_path, char(motion_roots(root_idx)), '**', 'shared_motion_shifts_result.mat'));
    candidate_listing = [candidate_listing; current_listing(:)]; %#ok<AGROW>
end
if isempty(candidate_listing)
    return;
end
[~, newest_order] = sort([candidate_listing.datenum], 'descend');
candidate_listing = candidate_listing(newest_order);
for idx = 1:numel(candidate_listing)
    candidate_file = fullfile(candidate_listing(idx).folder, candidate_listing(idx).name);
    has_voltage_tif = ~isempty(dir(fullfile(candidate_listing(idx).folder, 'voltage_motion_corrected_ds*.tif')));
    has_calcium_tif = ~isempty(dir(fullfile(candidate_listing(idx).folder, 'calcium_motion_corrected_ds*.tif')));
    if has_voltage_tif && has_calcium_tif
        motion_result = string(candidate_file);
        return;
    end
end
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
explicit_listing = filter_result_listing_by_preset_rec(explicit_listing, "analysis_only");
explicit_result = select_latest_file_from_listing_rec(explicit_listing);
end

function listing = filter_result_listing_by_preset_rec(listing, expected_preset)
keep = false(size(listing));
for idx = 1:numel(listing)
    dual_info_file = fullfile(listing(idx).folder, 'dual_info.mat');
    if ~isfile(dual_info_file)
        continue;
    end
    try
        saved = load(dual_info_file, 'dual_info');
        keep(idx) = isfield(saved, 'dual_info') ...
            && isstruct(saved.dual_info) ...
            && isfield(saved.dual_info, 'workflow_control') ...
            && isstruct(saved.dual_info.workflow_control) ...
            && isfield(saved.dual_info.workflow_control, 'preset') ...
            && strcmpi(string(saved.dual_info.workflow_control.preset), string(expected_preset));
    catch
        keep(idx) = false;
    end
end
listing = listing(keep);
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


function out = ternary_rec(condition, true_value, false_value)
if condition
    out = true_value;
else
    out = false_value;
end
end
