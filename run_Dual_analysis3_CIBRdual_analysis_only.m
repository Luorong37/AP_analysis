%% run_Dual_analysis3_CIBRdual_analysis_only
% Re-run saved ROI-after Dual_analysis3 sections for the CIBR dual dataset.
% This keeps the original saved movie/ROI/trace results and regenerates
% downstream analysis figures, including the ROI calcium heatmap with
% overlaid voltage sensitivity traces.

clearvars;
clc;

%% 1. Input / Output Paths
cycle_path = 'I:\1_Data\3b. Dual-color imaging in SCN\2024.09.05_P2A-G8s\20240905-153142POA';

analysis_mode = 'analysis_only';
reuse_results_path = 'I:\1_Data\3b. Dual-color imaging in SCN\2024.09.05_P2A-G8s\Dual_analysis3\2024.09_20240905-153142POA_2026-05-26 16-50-50';
save_path = '';
analysis_run_name = '';
analysis_backend = 'default';

%% 2. Camera Role / Orientation / Frame Rate
voltage_frame_rate = 400;
calcium_frame_rate = 10;

voltage_transpose_movie = false;
calcium_transpose_movie = false;

camera_cfg(1) = struct( ...
    'camera_index', 1, ...
    'role', "calcium", ...
    'transpose_before_analysis', calcium_transpose_movie, ...
    'frame_rate', calcium_frame_rate);

camera_cfg(2) = struct( ...
    'camera_index', 2, ...
    'role', "voltage", ...
    'transpose_before_analysis', voltage_transpose_movie, ...
    'frame_rate', voltage_frame_rate);

%% 3. Raw Two-Folder Compatibility
raw_dual_green_suffix = '_Green';
raw_dual_primary_role = "voltage";
raw_dual_green_role = "calcium";
camera_source_override = [];

%% 4. Motion / ROI Reuse
run_motion_correction = false;
motion_cfg = struct( ...
    'enabled', run_motion_correction, ...
    'use_saved_shift', false, ...
    'saved_shift_file', '', ...
    'highpass', true, ...
    'auto_reuse_previous_shift', true);

reuse_roi_file = '';
correct_offset_mode = 'none';
reuse_offset = [-1.22674328942136,	2.39504549476061];

%% 5. Trace Processing / Display Settings
map_bin = 4;
calcium_smoothing_window = 1;
bleach_mode_voltage = 'linear';
bleach_mode_calcium = 'exp2';
run_background_removal = true;
voltage_polarity = -1;
calcium_polarity = 1;

%% 6. VolPy / Stim / Runtime
volpy_source_results_path = '';
volpy_auto_run = false;
volpy_force_rerun = false;
volpy_use_existing = false;
volpy_flip_signal = false;
stim_context_override = [];
gpu = false;

%% 7. Run Dual_analysis3
script_dir = fileparts(mfilename('fullpath'));
run(fullfile(script_dir, 'Dual_analysis3.m'));
