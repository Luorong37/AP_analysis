%% run_Dual_analysis3_template
% Template for calling Dual_analysis3.m from an outer script.
%
% Usage:
%   1. Save this file in the same folder as Dual_analysis3.m.
%   2. Edit the parameter blocks below.
%   3. Run this template script.
%
% Note:
% Dual_analysis3.m does not clear the workspace automatically. This wrapper
% clears local variables first so old settings are not reused by accident.

clearvars;
clc;

%% 1. Input / Output Paths
% Point to one rebuilt Rec*/Cycle* folder, or to one folder from a raw
% two-folder pair if using the raw dual compatibility block below.
cycle_path = 'V:\Luorong\Invivo\26.06.02_dual_color_WT\Methods2_drifting_grating\Rec1_2026-06-02_21-05-24\Cycle1';

% Optional name for this run. Leave empty to let Dual_analysis3 build a
% name from record, cycle, and timestamp.
analysis_run_name = '';

% Optional output folder.
%   - Leave empty in analysis_only mode to write back into reuse_results_path.
%   - Set this to a new folder to reuse old traces while saving all rerun
%     figures/results into the new folder instead of overwriting the source.
%   - In full mode, leave empty for:
%     <cycle or raw-pair base>\Dual_analysis3\<analysis_run_name>
save_path = '';

%% 2. Analysis Mode
% 'full'          -> load movies, run ROI / processing / summaries
% 'analysis_only' -> reuse an existing Dual_analysis3 output folder and
%                    rerun ROI-after analysis sections from saved results
analysis_mode = 'analysis_only';

% Existing Dual_analysis3 output folder for analysis_only reruns. This must
% be a Dual_analysis3 result folder containing dual_info.mat,
% voltage_results.mat, calcium_results.mat, and dual_results.mat.
reuse_results_path = 'V:\Luorong\Invivo\26.06.02_dual_color_WT\Methods2_drifting_grating\Rec1_2026-06-02_21-05-24\Cycle1\Dual_analysis3\Rec1_2026-06-02_21-05-24_Cycle1_2026-06-03 00-30-07-384';

% Backend:
% 'default'                  -> standard dual analysis
% 'volpy_voltage_reanalysis' -> reuse a standard Dual_analysis3 run as the
%                               calcium/stim/ROI source, then reanalyze
%                               voltage with VolPy-derived events
analysis_backend = 'default';

%% 3. Camera Role / Orientation / Frame Rate
% Role-level settings are mapped onto camera_cfg below by Dual_analysis3.
% Default legacy mapping:
%   Cam1 -> calcium, transposed
%   Cam2 -> voltage, not transposed
voltage_frame_rate = 400;
calcium_frame_rate = 400;

voltage_transpose_movie = false;
calcium_transpose_movie = true;

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

%% 4. Raw Two-Folder Compatibility
% Use this when cycle_path is not a rebuilt Rec*/Cycle* folder.
% Default meaning:
%   primary folder       -> voltage
%   primary folder_Green -> calcium
raw_dual_green_suffix = '_Green';
raw_dual_primary_role = "voltage";
raw_dual_green_role = "calcium";

%% 5. Explicit Movie Source Override
% Leave [] for normal manifest / Cam* auto-detection.
% Use this only when you want to pass exact movie/folder paths manually.
%
% Example:
% camera_source_override(1) = struct( ...
%     'camera_index', 1, ...
%     'path', 'E:\1_Data\...\Cam1_Ca', ...
%     'label', "Cam1_Ca", ...
%     'original_frame_size', [NaN, NaN]);
% camera_source_override(2) = struct( ...
%     'camera_index', 2, ...
%     'path', 'E:\1_Data\...\Cam2_Voltage', ...
%     'label', "Cam2_Voltage", ...
%     'original_frame_size', [NaN, NaN]);
camera_source_override = [];

%% 6. Motion Correction
% Motion is estimated on voltage and applied to calcium to keep channels
% spatially locked.
run_motion_correction = false;

motion_cfg = struct( ...
    'enabled', run_motion_correction, ...
    'use_saved_shift', false, ...
    'saved_shift_file', '', ...
    'highpass', true, ...
    'auto_reuse_previous_shift', true);

%% 7. ROI Reuse / Channel Offset
% Reuse a previously saved ROI file. If reuse_results_path is provided,
% Dual_analysis3 gives it priority over reuse_roi_file.
reuse_roi_file = '';

% Inter-camera ROI offset mode:
% 'none'            -> use reuse_offset or [0 0]
% 'manual_points'   -> manually click one matching voltage/calcium point
% 'matlab_register' -> estimate translation with MATLAB registration
correct_offset_mode = 'none';

% Offset convention:
%   voltage_position = calcium_position + offset
% Leave [] unless manually forcing a known [x y] offset.
reuse_offset = [];

%% 8. Trace Processing
map_bin = 4;

% Calcium-only smoothing window after main metric computation.
calcium_smoothing_window = 40;

% Bleach correction mode:
% 'linear' | 'highpass' | 'exp2'
bleach_mode_voltage = 'linear';
bleach_mode_calcium = 'exp2';

% Background removal can be slow and interactive depending on ROI state.
run_background_removal = true;

% Display / analysis polarity.
voltage_polarity = -1;
calcium_polarity = 1;

%% 9. VolPy Voltage Reanalysis Options
% Used only when analysis_backend = 'volpy_voltage_reanalysis'.
volpy_source_results_path = '';
volpy_auto_run = false;
volpy_force_rerun = false;
volpy_use_existing = false;
volpy_flip_signal = false;

%% 10. Stimulus Metadata Override
% Leave [] for normal manifest/log discovery.
% Provide a struct here only if stimulus metadata/logs need to be supplied
% by the caller.
stim_context_override = [];

%% 11. Runtime Options
% gpu=true lets Dual_analysis3 try to open a parallel pool for full movie
% loading / motion stages. It continues without an explicit pool if that
% startup fails.
gpu = true;

%% 12. Run Dual_analysis3
script_dir = fileparts(mfilename('fullpath'));
run(fullfile(script_dir, 'Dual_analysis3.m'));
