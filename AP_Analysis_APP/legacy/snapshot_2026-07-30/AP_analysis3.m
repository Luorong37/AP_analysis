% AP_analysis3 - Single-channel voltage imaging analysis with dual-style result tracking
%
% Helper / readme:
%   - movie_info records only the movie source, geometry, and motion state.
%     The movie array itself is intentionally not duplicated in results.
%   - trace_results stores each trace-processing stage under one unified
%     structure: results.trace_results.<stage>.data + .info
%   - peak_results stores peak-detection and peak-gating stages in the
%     same style.
%   - ap_results stores downstream AP event, summary, sequence, and
%     density outputs without saving the whole workspace.
%   - each major section also saves a compact "section record" that says:
%       what inputs were used
%       what parameters were used
%       what outputs were produced
%       what was intentionally not saved
%       how to rerun the section later
%
% Main saved trace stages:
%   raw -> bg_removed -> bleach_removed -> baseline -> denoised
%   -> noise_reference -> noise -> sensitivity -> snr
%
% Standalone usage note:
%   - this script no longer clears the workspace automatically.
%   - if you want a clean standalone run, clear manually before running.
%   - if this script is called by an outer batch script, keeping the
%     workspace allows that script to pass override parameters in.
clc ;

%% Loading Raw Data
% Configure the input movie path, acquisition rate, and output analysis folder.
nowtime = string(datetime( 'now'));
% Replace colons with hyphens to get the desired output format
nowtime = strrep(nowtime , ':', '-');
fprintf('Loading...\n')


% support for folder, .tif, .tiff, .bin.
if ~exist('folder_path', 'var') || isempty(folder_path)
    folder_path = 'I:\1_Data\4e. Duplex imaging\25.10.29_DUPLEX_SCN\';
end
if ~exist('file', 'var') || isempty(file)
    file = 'R5';  % must add format.do not add '\' at last
end
if ~exist('bin', 'var') || isempty(bin)
    bin = 1;
end

if ~exist('freq', 'var') || isempty(freq)
    freq = 400; % Hz
end
if ~exist('gpu', 'var') || isempty(gpu)
    gpu = true; % defined gpu open
end

if ~exist('analysis_mode', 'var') || isempty(analysis_mode)
    analysis_mode = "full"; % "full" | "analysis_only"
else
    analysis_mode = string(analysis_mode);
end

if ~exist('transpose_before_analysis', 'var') || isempty(transpose_before_analysis)
    transpose_before_analysis = false;
end
if ~exist('analysis_run_name', 'var') || isempty(analysis_run_name)
    analysis_run_name = char(nowtime);
end

if ~exist('analysis_backend', 'var') || isempty(analysis_backend)
    analysis_backend = "classic"; %'volpy'"classic"
else
    analysis_backend = string(analysis_backend);
end
analysis_backend = lower(strtrim(string(analysis_backend)));
if ~exist('reuse_results_path', 'var') || isempty(reuse_results_path)
    reuse_results_path = "";
else
    reuse_results_path = string(reuse_results_path);
end


%%%%%%%%% motion part %%%%%%%%%%%%
if ~exist('run_motion_correction', 'var') || isempty(run_motion_correction)
    run_motion_correction = false;
end
if ~exist('motion_cfg', 'var') || isempty(motion_cfg)
    motion_cfg = struct( ...
        'enabled', true, ...
        'use_saved_shift', false, ...
        'saved_shift_file', '', ...
        'highpass', true, ...
        'nonrigid', true, ...
        'plot_metrics', true, ...
        'save_downsampled_tif', true);
end
motion_cfg.enabled = logical(run_motion_correction);
if ~isfield(motion_cfg, 'use_saved_shift') || isempty(motion_cfg.use_saved_shift)
    motion_cfg.use_saved_shift = false;
end
if ~isfield(motion_cfg, 'saved_shift_file') || isempty(motion_cfg.saved_shift_file)
    motion_cfg.saved_shift_file = '';
end
if ~isfield(motion_cfg, 'highpass') || isempty(motion_cfg.highpass)
    motion_cfg.highpass = true;
end
if ~isfield(motion_cfg, 'nonrigid') || isempty(motion_cfg.nonrigid)
    motion_cfg.nonrigid = true;
end
if ~isfield(motion_cfg, 'plot_metrics') || isempty(motion_cfg.plot_metrics)
    motion_cfg.plot_metrics = true;
end
if ~isfield(motion_cfg, 'save_downsampled_tif') || isempty(motion_cfg.save_downsampled_tif)
    motion_cfg.save_downsampled_tif = true;
end

%%%%%%%%% volpy part %%%%%%%%%%%%
if ~exist('volpy_auto_run', 'var') || isempty(volpy_auto_run)
    volpy_auto_run = true;
end
if ~exist('volpy_force_rerun', 'var') || isempty(volpy_force_rerun)
    volpy_force_rerun = false;
end
if ~exist('volpy_use_existing', 'var') || isempty(volpy_use_existing)
    volpy_use_existing = true;
end
if ~exist('volpy_input_file', 'var')
    volpy_input_file = "";
end
if ~exist('volpy_work_dir', 'var')
    volpy_work_dir = "";
end
if ~exist('volpy_frame_rate', 'var')
    volpy_frame_rate = [];
end
if ~exist('volpy_flip_signal', 'var') || isempty(volpy_flip_signal)
    volpy_flip_signal = true;
end
if ~exist('volpy_skip_motion_correction', 'var') || isempty(volpy_skip_motion_correction)
    volpy_skip_motion_correction = false;
end
if ~exist('volpy_roi_mask_file', 'var')
    volpy_roi_mask_file = "";
end
if ~exist('volpy_temp_dir', 'var')
    volpy_temp_dir = "";
end
if ~exist('volpy_result_mat', 'var')
    volpy_result_mat = "";
end
if ~exist('volpy_size_min', 'var') || isempty(volpy_size_min)
    volpy_size_min = 5;
end
if ~exist('volpy_size_max', 'var') || isempty(volpy_size_max)
    volpy_size_max = 22;
end
if ~exist('volpy_corr_window_seconds', 'var') || isempty(volpy_corr_window_seconds)
    volpy_corr_window_seconds = 4;
end
if ~exist('volpy_corr_stride_seconds', 'var') || isempty(volpy_corr_stride_seconds)
    volpy_corr_stride_seconds = 4;
end
if ~exist('volpy_corr_baseline_seconds', 'var') || isempty(volpy_corr_baseline_seconds)
    volpy_corr_baseline_seconds = 1;
end
if ~exist('volpy_corr_remove_baseline', 'var') || isempty(volpy_corr_remove_baseline)
    volpy_corr_remove_baseline = true;
end
if ~exist('volpy_corr_gaussian_blur', 'var') || isempty(volpy_corr_gaussian_blur)
    volpy_corr_gaussian_blur = false;
end
if ~exist('volpy_summary_mode', 'var') || isempty(volpy_summary_mode)
    volpy_summary_mode = "mean_mean_corr";
else
    volpy_summary_mode = string(volpy_summary_mode);
end
if ~exist('volpy_mrcnn_input_mode', 'var') || isempty(volpy_mrcnn_input_mode)
    volpy_mrcnn_input_mode = "raw_zscore";
else
    volpy_mrcnn_input_mode = string(volpy_mrcnn_input_mode);
end
if ~exist('volpy_mrcnn_confidence_threshold', 'var') || isempty(volpy_mrcnn_confidence_threshold)
    volpy_mrcnn_confidence_threshold = 0.7;
end
if ~exist('run_manual_peak_gating', 'var') || isempty(run_manual_peak_gating)
    run_manual_peak_gating = analysis_backend ~= "volpy";
end
if ~exist('run_manual_peak_refinement', 'var') || isempty(run_manual_peak_refinement)
    run_manual_peak_refinement = analysis_backend ~= "volpy";
end
if ~exist('show_roi_saved_dialog', 'var') || isempty(show_roi_saved_dialog)
    show_roi_saved_dialog = analysis_backend ~= "volpy";
end
if ~exist('save_failed_peak_debug', 'var') || isempty(save_failed_peak_debug)
    save_failed_peak_debug = false;
end
if ~exist('enable_peak_redirection', 'var') || isempty(enable_peak_redirection)
    enable_peak_redirection = true;
end


if exist('v','var')
    matim =true;
    [folder_path,file,exten] = fileparts(save_file);
    file = [file,exten];
else
    matim = false; % imaging via matlab
end
% loadtype = "mat"; % default "tif" = loading tif file, or "mat" = loading .mat file
% mat_name = 'LEDcyan_rec1_cycles1.mat';     % if loading .mat file

% -----------------------------------------------------------

% Create an analysis folder

file_path = fullfile(folder_path, file);
if isfolder(file_path)
    file_name = file;
    %file_dir = dir(file_path);
    %[~, ~, file_extension] = fileparts(file_dir(3).name);
else
    % [~, file_name, file_extension] = fileparts(file_path)
    % Comment removed after encoding repair.
    [folder_path, file_name, fext] = fileparts(file_path);
    % [folder_path, file_name, ~] = fileparts(file_path);
end

% create a folder for analysis
if ~exist('save_path', 'var') || isempty(save_path)
    if strcmpi(string(analysis_mode), "analysis_only") && strlength(reuse_results_path) > 0
        save_path = char(reuse_results_path);
    else
        save_path = fullfile(folder_path, strcat(file_name, '_Analysis'), analysis_run_name);
    end
end
mkdir(save_path);
analysis_only_mode = strcmpi(string(analysis_mode), "analysis_only");
if analysis_only_mode && analysis_backend == "volpy"
    warning('AP_analysis3:AnalysisOnlyClassicPipeline', ...
        'analysis_only resumes after bg_removed and uses the classic downstream trace pipeline. Switching analysis_backend from volpy to classic.');
    analysis_backend = "classic";
end

% Parallel pool setup happens once, before any movie loading. Downstream
% helpers only reuse an existing pool and fall back to serial work when no
% pool is available.
parallel_pool_info = struct( ...
    'requested', logical(gpu && analysis_backend ~= "volpy" && ~analysis_only_mode), ...
    'ready', false, ...
    'num_workers', 0, ...
    'message', "", ...
    'initialized_before_movie_load', true);
if parallel_pool_info.requested
    fprintf('\n===== Parallel Pool Setup =====\n');
    [~, parallel_pool_ready, parallel_pool_info] = initialize_ap_parallel_pool(parallel_pool_info);
    parallel_pool_info.ready = parallel_pool_ready;
    gpu = parallel_pool_ready;
else
    parallel_pool_info.message = "Parallel pool not requested for this mode.";
end

% Load image file

if analysis_only_mode
    fprintf('Loading saved AP_analysis3 context for analysis_only mode...\n');
    analysis_only_source_path = save_path;
    if strlength(reuse_results_path) > 0
        analysis_only_source_path = char(reuse_results_path);
    end
    [movie_info, trace_results, ~, ~, ap_context] = ...
        load_saved_ap_analysis_only_context(analysis_only_source_path);
    if ~strcmpi(char(string(save_path)), char(string(analysis_only_source_path)))
        copied_roi_results_file = fullfile(save_path, '1_raw_ROI.mat');
        copyfile(ap_context.roi_results_file, copied_roi_results_file);
        ap_context.roi_results_file = copied_roi_results_file;
    end
    file_path = char(movie_info.source_path);
    freq = double(movie_info.frame_rate);
    ncols = double(movie_info.analysis_frame_size(1));
    nrows = double(movie_info.analysis_frame_size(2));
    movie = [];
    traces_raw = ap_context.traces_raw;
    traces_background_removed = ap_context.traces_background_removed;
    nframes = size(traces_raw, 1);
    avg_image = ap_context.avg_image;
    map = ap_context.map;
    mask = ap_context.mask;
    rois = ap_context.rois;
    roi_results_file = ap_context.roi_results_file;
    trace_results = keep_ap_trace_results(trace_results, ["raw", "bg_removed"]);
    peak_results = struct();
    ap_results = struct();
    nrois = size(traces_raw, 2);
    matim = true;
elseif analysis_backend == "volpy"
    movie = [];
    tif_info = imfinfo(file_path);
    nframes = numel(tif_info);
    ncols = tif_info(1).Height;
    nrows = tif_info(1).Width;
elseif ~matim
    fprintf("Start loading movie, please wait...\n");
    tload = tic;
    % if loadtype == "tif"
    [movie, ncols, nrows, nframes] = load_movie_with_fallback(file_path);
    % end
    % if loadtype == "mat"
    %     mat_path = fullfile(file_path, mat_name);
    %     load(mat_path);
    %     [ncols, nrows, nframes] = size(movie);
    % end
    fprintf("Finished loading movie after %.2d s\n", toc(tload));
else
    [ncols, nrows, nframes] = size(movie);
    movie = reshape(movie,ncols*nrows, []);
end

% Presetting

% Define parameters
dt = 1 / freq; % Calculate time axis
% colors = [lines(7);hsv(5);spring(3);winter(3);gray(3)];
colors = lines(256);
t = (1:nframes) * dt;
options.colors = colors;
if ~analysis_only_mode || ~exist('map', 'var')
    map = [];
end
if ~analysis_only_mode || ~exist('mask', 'var')
    mask = [];
end
x = (1:nframes)' * dt;
ap_qc_trace = [];
ap_qc_sensitivity = [];
ap_qc_SNR = [];
ap_qc_baseline_subtract = false;

if analysis_only_mode
    if ~exist('movie_vol_2D', 'var') || isempty(movie_vol_2D)
        movie_vol_2D = avg_image;
    end
elseif analysis_backend == "volpy"
    movie_vol_2D = zeros(ncols, nrows, 'single');
    avg_image = zeros(ncols, nrows, 'single');
else
    try
        movie_vol_2D = reshape(mean(movie,2), ncols, nrows, []);
    catch ME
        movie_vol_2D = squeeze(mean(movie,3));
    end

    avg_image  = (movie_vol_2D - min(movie_vol_2D(:))) / (max(movie_vol_2D(:)) - min(movie_vol_2D(:)));
end

% movie_info follows the same layout as Dual_analysis3 so later scripts can
% inspect single-channel and dual-channel outputs in a similar way.
if ~analysis_only_mode
movie_info = struct( ...
    'role', "voltage", ...
    'label', "single_channel_voltage", ...
    'source_path', file_path, ...
    'frame_rate', freq, ...
    'frame_count', nframes, ...
    'original_frame_size', [ncols, nrows], ...
    'analysis_frame_size', [ncols, nrows], ...
    'transpose_before_analysis', transpose_before_analysis, ...
    'parallel_pool_info', parallel_pool_info, ...
    'motion', struct( ...
        'applied', false, ...
        'method', '', ...
        'shift_file', '', ...
        'parameter_file', ''), ...
    'created_at', datetime("now"), ...
    'updated_at', datetime("now"));
end
movie_info_path = fullfile(save_path, 'movie_info.mat');
save(movie_info_path, 'movie_info');

if ~analysis_only_mode
analysis_info = struct( ...
    'analysis_name', 'AP_analysis3', ...
    'save_path', save_path, ...
    'source_path', file_path, ...
    'frame_rate', freq, ...
    'bin', bin, ...
    'gpu', gpu, ...
    'analysis_mode', string(analysis_mode), ...
    'reuse_results_path', string(reuse_results_path), ...
    'run_motion_correction', logical(run_motion_correction), ...
    'motion_cfg', motion_cfg, ...
    'parallel_pool_info', parallel_pool_info, ...
    'transpose_before_analysis', transpose_before_analysis, ...
    'created_at', datetime("now"));
else
analysis_info = struct( ...
    'analysis_name', 'AP_analysis3', ...
    'save_path', save_path, ...
    'source_path', file_path, ...
    'frame_rate', freq, ...
    'bin', bin, ...
    'gpu', gpu, ...
    'parallel_pool_info', parallel_pool_info, ...
    'transpose_before_analysis', movie_info.transpose_before_analysis, ...
    'analysis_mode', string(analysis_mode), ...
    'reuse_results_path', string(save_path), ...
    'run_motion_correction', logical(run_motion_correction), ...
    'motion_cfg', motion_cfg, ...
    'created_at', datetime("now"));
end

% trace_results is the main trace provenance container. Each stage stores
% a data matrix plus metadata describing its parents and processing method.
if ~analysis_only_mode
trace_results = struct();
end
trace_results_path = fullfile(save_path, 'trace_results.mat');
save(trace_results_path, 'trace_results', '-v7.3');

if ~analysis_only_mode
peak_results = struct();
end
peak_results_path = fullfile(save_path, 'peak_results.mat');
if ~analysis_only_mode
save(peak_results_path, 'peak_results', '-v7.3');
end

if ~analysis_only_mode
ap_results = struct();
end
ap_results_path = fullfile(save_path, 'ap_results.mat');
if ~analysis_only_mode
save(ap_results_path, 'ap_results', '-v7.3');
end




% Save code
code_path = fullfile(save_path,'Code');
mkdir(code_path);
currentScript = which("AP_analysis3.m");
% Encoding-corrupted comment removed.
[requiredFiles, ~] = matlab.codetools.requiredFilesAndProducts(currentScript);
% Encoding-corrupted comment removed.
for k = 1:length(requiredFiles)
    [~, name, ext] = fileparts(requiredFiles{k});
    copyfile(requiredFiles{k}, fullfile(code_path, [name, ext]));
end
% Comment removed after encoding repair.
fprintf('All codes have been copied to %s\n', code_path);

if string(analysis_backend) == "classic"
if ~analysis_only_mode
%% Motion Correction
% Estimate or reuse motion shifts, update movie_info, and keep movie in memory only.
fprintf('Initializing Motion Correction...\n');
% Encoding-corrupted comment removed.
if ismatrix(movie)
    movie = reshape(movie, ncols, nrows, []);
end
 
% Comment removed after encoding repair.
if motion_cfg.enabled

apply_only = logical(motion_cfg.use_saved_shift);
loadMC     = 0;
Norigid    = logical(motion_cfg.nonrigid);      % use nonrigid correction
hp         = logical(motion_cfg.highpass);       % use high-pass filtering for shift estimation
template   = [];                                % optional manual correction template
plotmetric = logical(motion_cfg.plot_metrics);
dssave     = logical(motion_cfg.save_downsampled_tif);
% Comment removed after encoding repair.
init_batch = 100; % can be modified manually

options_r = NoRMCorreSetParms('d1',ncols,'d2',nrows,'bin_width',200,'max_shift',50,'us_fac',30,'iter',1,'correct_bidir',false);
options_nr = NoRMCorreSetParms('d1',ncols,'d2',nrows,'bin_width',200,'max_shift',30,'us_fac',30, ...
    'grid_size',[128,128],'overlap_pre',[32,32],'mot_uf',4,'max_dev', [5,5],'iter',1,'correct_bidir',false);


% Comment removed after encoding repair.
[folder_path, file_name, fext] = fileparts(file_path);
shift_res_path = fullfile(save_path, 'motion_shifts_result.mat');
params_save_path = fullfile(save_path, 'motion_correction_para.mat');

% Comment removed after encoding repair.
if loadMC
    [Mr, ncols, nrows, ~] = load_movie(mc_path);
else
    t1 = tic;

    % Comment removed after encoding repair.
    movie = single(movie);
    % Comment removed after encoding repair.

    if apply_only
        % Encoding-corrupted comment removed.
        if strlength(string(motion_cfg.saved_shift_file)) > 0
            saved_shift_path = char(string(motion_cfg.saved_shift_file));
        else
            [shift_filename, shift_foldername] = uigetfile(save_path, 'Select motion shift file');
            if isequal(shift_filename,0), return; end
            saved_shift_path = fullfile(shift_foldername, shift_filename);
        end
        S = load(saved_shift_path);
        plotmetric = false;

        fprintf(' -> Mode: Apply existing shifts...\n');
        Mr = apply_shifts(movie, S.shifts_r, options_r);
        if Norigid && isfield(S, 'shifts_nr')
            Mpr = apply_shifts(Mr, S.shifts_nr, S.options_nr);
        else
            Norigid = false;
        end
    else
        % Comment removed after encoding repair.
        if hp
            % Encoding-corrupted comment removed.
            fprintf(' -> High-pass mode: Filtering for better estimation...\n');
            Y = create_temp_highpass(movie);
        else
            Y = movie;
        end

        fprintf(' -> Mode: Estimating shifts with rigid motion...\n');
        % Encoding-corrupted comment removed.
        if plotmetric || Norigid
            [M1, shifts_r, template1] = normcorre_batch(Y, options_r);
        else
            [~, shifts_r, template1] = normcorre_batch(Y, options_r);
            clear Y
        end
        Mr = apply_shifts(movie, shifts_r, options_r);

        % Encoding-corrupted comment removed.

        % Encoding-corrupted comment removed.
        shifts_nr = [];
        if Norigid
            fprintf(' -> Mode: Estimating shifts with Norigid motion...\n');
            [M2, shifts_nr, ~] = normcorre_batch(M1, options_nr, template1);
            Mpr = apply_shifts(Mr, shifts_nr, options_nr);
        end

        % Comment removed after encoding repair.
        fprintf(' -> Saving shifts and params...\n');
        save(shift_res_path, 'shifts_r', 'shifts_nr', 'template1', '-v7.3');
        save(params_save_path, 'options_r', 'options_nr', 'hp', 'Norigid');
    end

    tsub = 40;
    % Encoding-corrupted comment removed.
    fprintf(' -> Downsampling corrected movie (tsub = 40) for saving...\n');

    % Comment removed after encoding repair.
    
    if dssave
        % Encoding-corrupted comment removed.
        % Encoding-corrupted comment removed.
        if Norigid && exist('Mpr', 'var')
            Mr_ds = downsample_data(Mpr, 'time', tsub);
        else
            Mr_ds = downsample_data(Mr, 'time', tsub);
        end

        % Comment removed after encoding repair.
        nn_ds = quantile(Mr_ds(:), 0.0005);
        mm_ds = quantile(Mr_ds(:), 0.99995);

        % Encoding-corrupted comment removed.
        % Encoding-corrupted comment removed.
        Mr_ds = (Mr_ds - nn_ds) / (mm_ds - nn_ds) * 65535;
        Mr_ds(Mr_ds < 0) = 0;
        Mr_ds(Mr_ds > 65535) = 65535;

        % Comment removed after encoding repair.
        save_name_ds = fullfile(save_path, ['motion_corrected_ds' mat2str(tsub) '.tif']);
        array2tif(uint16(Mr_ds), save_name_ds);

        % Encoding-corrupted comment removed.
        movie_vol_2D = mean(Mr_ds, 3);
        fprintf(' -> Downsampled movie saved. \n');
    else
        % Comment removed after encoding repair.
        save_name = fullfile(folder_path, [file_name, '_motion_correction.tif']);
        array2tif(uint16(Mr), save_name);
        movie_vol_2D = mean(Mr, 3);
        fprintf('Finished in %d s\n', round(toc(t1)));
    end

end









% Comment removed after encoding repair.
avg_image = (movie_vol_2D - min(movie_vol_2D(:))) ./ (max(movie_vol_2D(:)) - min(movie_vol_2D(:)));



if plotmetric
    fprintf(' -> Computing metrics. \n');
    % compute metrics

    [cY,mY,vY] = motion_metrics(Y,options_r.max_shift);
    [cYf,mYf,vYf] = motion_metrics(reshape(movie, ncols,nrows, []),options_r.max_shift);

    [cM1,mM1,vM1] = motion_metrics(M1,options_r.max_shift);
    [cM1f,mM1f,vM1f] = motion_metrics(Mr,options_r.max_shift);

    % plot rigid shifts and metrics
    shifts_rplot = squeeze(cat(3,shifts_r(:).shifts));
    figure;
    subplot(311); plot(shifts_rplot);
    title('Rigid shifts','fontsize',14,'fontweight','bold');
    legend('y-shifts','x-shifts');
    subplot(312); plot(t,cY,t,cM1);
    title('Correlation coefficients on filtered movie','fontsize',14,'fontweight','bold');
    legend('raw','rigid');
     set(gca,'Xtick',[],'Ylim',[0.8,1])
    subplot(313); plot(t,cYf,t,cM1f);
     set(gca,'Xtick',[],'Ylim',[0.8,1])
    title('Correlation coefficients on full movie','fontsize',14,'fontweight','bold');
    legend('raw','rigid');
    saveas(gcf,fullfile(save_path,'rigidshifts.fig'))
    saveas(gcf,fullfile(save_path,'rigidshifts.png'))
    if Norigid
        % plot shifts   compute metrics

        [cM2,mM2,vM2] = motion_metrics(M2,options_nr.max_shift);
        [cM2f,mM2f,vM2f] = motion_metrics(Mpr,options_nr.max_shift);
        shifts_nrplot = cat(ndims(shifts_nr(1).shifts)+1,shifts_nr(:).shifts);
        shifts_nrplot = reshape(shifts_nrplot,[],ndims(Y)-1,nframes);
        shifts_x = squeeze(shifts_nrplot(:,2,:))';
        shifts_y = squeeze(shifts_nrplot(:,1,:))';

        patch_id = 1:size(shifts_x,2);
        str = strtrim(cellstr(int2str(patch_id.')));
        str = cellfun(@(x) ['patch # ',x],str,'un',0);

        figure;
        ax1 = subplot(311); plot(t,cY,t,cM1,t,cM2); legend('raw data','rigid','non-rigid'); title('correlation coefficients for filtered data','fontsize',14,'fontweight','bold')
        set(gca,'Xtick',[],'XLim',[0,nframes-3])
        ax2 = subplot(312); plot(shifts_x); hold on; plot( shifts_rplot(:,2),'--k','linewidth',2); title('displacements along x','fontsize',14,'fontweight','bold')
        set(gca,'Xtick',[],'Ylim',[0,1])
 
        ax3 = subplot(313); plot(shifts_y); hold on; plot( shifts_rplot(:,1),'--k','linewidth',2); title('displacements along y','fontsize',14,'fontweight','bold')
        set(gca,'Xtick',[],'Ylim',[0,1])
        xlabel('timestep','fontsize',14,'fontweight','bold')
        linkaxes([ax1,ax2,ax3],'x')
        saveas(gcf,fullfile(save_path,'nonrigidshifts.fig'))
        saveas(gcf,fullfile(save_path,'nonrigidshifts.png'))
    end
end

if ~Norigid
    movie = reshape(uint16(Mr), ncols*nrows, []);
else
    movie = reshape(uint16(Mpr), ncols*nrows, []);
end

% Encoding-corrupted comment removed.
movie_info.motion.applied = true;
movie_info.motion.method = 'NoRMCorre';
movie_info.motion.shift_file = shift_res_path;
movie_info.motion.parameter_file = params_save_path;
movie_info.updated_at = datetime("now");
save(movie_info_path, 'movie_info');

else
    fprintf('Motion correction skipped by run_motion_correction=false.\n');
    movie_vol_2D = squeeze(mean(movie, 3));
    avg_range = max(movie_vol_2D(:)) - min(movie_vol_2D(:));
    if avg_range > 0
        avg_image = (movie_vol_2D - min(movie_vol_2D(:))) ./ avg_range;
    else
        avg_image = zeros(size(movie_vol_2D), 'like', movie_vol_2D);
    end
    movie = reshape(uint16(movie), ncols*nrows, []);
    movie_info.motion.applied = false;
    movie_info.motion.method = 'skipped_by_run_motion_correction';
    movie_info.motion.shift_file = '';
    movie_info.motion.parameter_file = '';
    movie_info.updated_at = datetime("now");
    save(movie_info_path, 'movie_info');
end

% Encoding-corrupted comment removed.
%% Create Sensitivity Map
% The active temporary high-pass helper is defined at the end of this file.
% Generate a quick activity map from the motion-aligned movie for ROI selection.
t1 = tic; % Start a timer
fprintf('Creating a map...\n')
% if the map cannot figure out active cells, please large the bin.
mapbin = 4; % defined bin = 4

[quick_map] = create_map(movie, nrows, ncols, mapbin);

% Encoding-corrupted comment removed.
quick_map(1,:) = quick_map(5,:);
quick_map(2,:) = quick_map(5,:);
quick_map(3,:) = quick_map(5,:);
quick_map(4,:) = quick_map(5,:);

map = quick_map;

% Visualize correlation coefficients as heatmap
figure()
imagesc(quick_map);
colormap turbo;
title('Sensitivity Map');
axis image;
colorbar;
fig_filename = fullfile(save_path, '0_Sensitivity_Map.fig');
png_filename = fullfile(save_path, '0_Sensitivity_Map.png');
mat_filename = fullfile(save_path, '0_Sensitivity_Map.mat');

saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');
save(mat_filename, 'map');


t2 = toc(t1); % Get the elapsed time
fprintf('Finished mask creating after %d s\n',round(t2))
%% Load Or Create Mask
% Reuse saved ROIs or create a cellpose mask before manual ROI extraction.
if exist('mask_method_override', 'var') && ~isempty(mask_method_override)
    methods = char(string(mask_method_override));
    fprintf('Mask method override: %s\n', methods);
else
    methods = questdlg('Load previous data?','load data','previous ROIs','cellpose','volpy ROI');
end
switch methods
    case 'previous ROIs'
        if exist('mask_file_override', 'var') && ~isempty(mask_file_override)
            rois_data = load(char(string(mask_file_override)));
        else
            [roi_filename,roi_foldername] = uigetfile(save_path);
            rois_data = load(fullfile(roi_foldername,roi_filename));
        end
        try
            rois = rois_data.rois;
            mask = rois.bwmask;
        catch ME

            if isfield(rois_data, 'bwmask')
                rois.bwmask = rois_data.bwmask;
            elseif isfield(rois_data, 'mask')
                rois.bwmask = rois_data.mask;
            else
                rethrow(ME);
            end
            mask = rois.bwmask;
            boundaries = bwboundaries(mask, 'noholes');
    
        % Comment removed after encoding repair.
        % Encoding-corrupted comment removed.
        current_boundary = boundaries{1};
        rois.boundary = [current_boundary(:,2), current_boundary(:,1)]; % [X, Y]
        
        % Comment removed after encoding repair.
        % Encoding-corrupted comment removed.
        % Encoding-corrupted comment removed.
        stats = regionprops(mask, 'BoundingBox');
        rois.Position = stats(1).BoundingBox;

        end
    case 'cellpose'
        fprintf('cellpose running......\n');
        cp = cellpose(ExecutionEnvironment="gpu");
        avgdia = 25;
        % gamma_image = imadjust(map,[],[],1.2); % recommend raise the gamma factor from 1 to 4.
        gamma_image = imadjust(avg_image,[],[],1.2); % recommend raise the gamma factor from 1 to 4.
        mask = segmentCells2D(cp, gamma_image , ImageCellDiameter = avgdia,FlowErrorThreshold = 3,CellThreshold = -6);% CellThreshold = -2, ,  FlowErrorThreshold = 2
        fprintf('%d cells are found.\n',max(mask(:)));
        figure()
        % Comment removed after encoding repair.
        overlayImage = labeloverlay(gamma_image, mask,'Transparency', 0.6);
        % Comment removed after encoding repair.
        imshow(overlayImage);
        title('cellpose mask');
        fig_filename = fullfile(save_path, '0_cellpose_mask.fig');
        png_filename = fullfile(save_path, '0_cellpose_mask.png');
        mat_filename = fullfile(save_path, '0_cellpose_mask.mat');

        saveas(gcf, fig_filename, 'fig');
        saveas(gcf, png_filename, 'png');
        save(mat_filename, 'mask');
    case 'volpy ROI'
        fprintf('volpy ROI running......\n');
        if exist('save_name_ds', 'var') && isfile(save_name_ds)
            volpy_input_path = save_name_ds;
            volpy_frame_rate = freq / tsub;
        elseif exist('save_name', 'var') && isfile(save_name)
            volpy_input_path = save_name;
            volpy_frame_rate = freq;
        else
            volpy_input_path = file_path;
            volpy_frame_rate = freq;
        end
        [mask, volpy_roi_info] = run_volpy_roi_segmentation(volpy_input_path, save_path, volpy_frame_rate);
        fprintf('%d cells are found.\n', max(mask(:)));
        figure()
        overlayImage = labeloverlay(avg_image, mask, 'Transparency', 0.6);
        imshow(overlayImage);
        title('volpy ROI mask');
        fig_filename = fullfile(save_path, '0_volpy_mask.fig');
        png_filename = fullfile(save_path, '0_volpy_mask.png');
        mat_filename = fullfile(save_path, '0_volpy_mask.mat');
        saveas(gcf, fig_filename, 'fig');
        saveas(gcf, png_filename, 'png');
        save(mat_filename, 'mask', 'volpy_roi_info');
    case 'manual'
        fprintf('manual ROI mode selected. Skipping automatic mask creation.\n');
end
%% Select ROI
% Extract raw ROI traces and register them as trace_results.raw.
t1 = tic; % Start a timer
figure()
% with or without Mask and Map
if exist('methods','var')
    if any(strcmp(methods,{'No','','manual'} ))
        [rois, traces] = select_ROI(movie, ncols, nrows, mask, map);
        nrois = max(rois.bwmask,[],'all');
    else
        [~, traces] = select_ROI(movie, ncols, nrows, mask, map);
        nrois = max(mask(:));
        rois.bwmask = mask;
    end                               
else
    [rois, traces] = select_ROI(movie, ncols, nrows, mask, map);
    nrois = max(rois.bwmask,[],'all');
end
% try
%     bwmask = rois.bwmask;
% catch ME
% 
% end
traces_raw = traces;

% if do not need a map, run the following code:
% map = [];
% [bwmask, traces] = select_ROI(movie, nrows, ncols, t, mask, map);

fig_filename = fullfile(save_path, '1_raw_trace.fig');
png_filename = fullfile(save_path, '1_raw_trace.png');
roi_results_file = fullfile(save_path, '1_raw_ROI.mat');

saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');
roi_record = build_section_record( ...
    'Select ROI', ...
    'Select voltage ROIs and extract the raw trace matrix that becomes the starting point of all later trace stages.', ...
    struct( ...
        'nrows', nrows, ...
        'ncols', ncols, ...
        'map_used', ~isempty(map), ...
        'mask_used', ~isempty(mask)), ...
    struct( ...
        'selection_function', 'select_ROI', ...
        'transpose_before_analysis', transpose_before_analysis), ...
    struct( ...
        'rois', rois, ...
        'avg_image', avg_image, ...
        'traces_raw', traces_raw, ...
        'nrois', nrois), ...
    struct( ...
        'movie', 'not saved'), ...
    'To rerun ROI selection, reload the movie externally and call select_ROI with the saved geometry / mask / map settings from this record.');
save(roi_results_file, 'rois', 'avg_image', 'traces', 'roi_record');

trace_results = store_trace_stage( ...
    trace_results, 'raw', traces_raw, {}, roi_results_file, ...
    movie_info, 'roi_selection', struct());
save(trace_results_path, 'trace_results', '-v7.3');

if show_roi_saved_dialog
    h = msgbox('All rois saved.', 'Done', 'help');
    uiwait(h);
end


tif_filename = fullfile(save_path, '1_Averaged.tif');

tiffile = Tiff(tif_filename, 'w');

tiffile.setTag('ImageLength', nrows);
tiffile.setTag('ImageWidth', ncols);
tiffile.setTag('Photometric', Tiff.Photometric.MinIsBlack);
tiffile.setTag('BitsPerSample', 16);
tiffile.setTag('SamplesPerPixel', 1);
tiffile.setTag('RowsPerStrip', 16);
tiffile.setTag('PlanarConfiguration', Tiff.PlanarConfiguration.Chunky);
tiffile.setTag('Compression', Tiff.Compression.None);

tiffile.setTag('Software', 'MATLAB');
tiffile.write(uint16(movie_vol_2D));
tiffile.writeDirectory();
tiffile.close();


%close(gcf);
end
%% Background Removal
% Background removal is the first persisted trace-processing stage.
fprintf('Removing Background...\n')
background_results_file = fullfile(save_path, '1_background_results.mat');
background_computed = exist('movie', 'var') == 1 && ~isempty(movie);
if analysis_only_mode
fprintf('Analysis-only starts after background removal. Reusing saved bg_removed traces.\n');
elseif background_computed
[background, background_fitted, traces_bgcorr, traces_bgfitcorr, background_mask]...
    = remove_background(movie, ncols, nrows, rois, freq, bin);
traces_bg_removed = traces_bgfitcorr;
background_record = build_section_record( ...
    'Background Removal', ...
    'Estimate an ROI-matched background signal and subtract it from the raw traces.', ...
    struct( ...
        'parent_stage', "raw", ...
        'traces_raw', traces_raw, ...
        'rois', rois, ...
        'nrows', nrows, ...
        'ncols', ncols), ...
    struct( ...
        'function_name', 'remove_background', ...
        'freq', freq, ...
        'bin', bin), ...
    struct( ...
        'background', background, ...
        'background_fitted', background_fitted, ...
        'background_mask', background_mask, ...
        'traces_bgcorr', traces_bgcorr, ...
        'traces_bg_removed', traces_bg_removed), ...
    struct( ...
        'movie', 'not saved'), ...
    'To rerun this section, reload the movie externally and call remove_background with the saved ROI definition and parameters from this record.');
save(background_results_file, 'background', 'background_fitted', 'background_mask', 'traces_bgcorr', 'traces_bgfitcorr', 'background_record');
trace_results = store_trace_stage( ...
    trace_results, 'bg_removed', traces_bg_removed, {'raw'}, roi_results_file, ...
    movie_info, 'remove_background', struct('bin', bin, 'freq', freq));
traces_background_removed = traces_bg_removed;
else
fprintf('Background removal skipped because movie data are not in memory. Downstream sections will use raw traces.\n');
traces_background_removed = traces_raw;
background_record = build_section_record( ...
    'Background Removal', ...
    'Skip background removal because movie data are not in memory.', ...
    struct( ...
        'parent_stage', "raw", ...
        'traces_raw', traces_raw, ...
        'rois', rois, ...
        'nrows', nrows, ...
        'ncols', ncols), ...
    struct( ...
        'function_name', 'remove_background', ...
        'skipped', true, ...
        'reason', 'movie_not_in_memory'), ...
    struct( ...
        'traces_bg_removed', traces_background_removed), ...
    struct( ...
        'movie', 'not available'), ...
    'To recompute background removal, rerun full analysis or reload the movie before this section.');
save(background_results_file, 'traces_background_removed', 'background_record');
trace_results = store_trace_stage( ...
    trace_results, 'bg_removed', traces_background_removed, {'raw'}, roi_results_file, ...
    movie_info, 'skipped_background_removal', struct('reason', 'movie_not_in_memory'));
end
save(trace_results_path, 'trace_results', '-v7.3');
fprintf('Finished\n')

% Comment removed after encoding repair.
if background_computed
fprintf('Generating summary plots...\n')

% Encoding-corrupted comment removed.
num_rois = max(rois.bwmask(:));
colors = lines(num_rois);

figure('Name', 'Background Correction Summary', 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.7]);

% Encoding-corrupted comment removed.
subplot(1, 2, 1);
movie2D_mean = mean(reshape(movie, ncols, nrows, []), 3);
imshow(movie2D_mean, [], 'InitialMagnification', 'fit');
hold on;

for i = 1:num_rois
    % Encoding-corrupted comment removed.
    current_color = colors(i, :);

    % Comment removed after encoding repair.
    roi_bw = (rois.bwmask == i);
    roi_boundaries = bwboundaries(roi_bw);
    for k = 1:length(roi_boundaries)
        boundary = roi_boundaries{k};
        plot(boundary(:, 2), boundary(:, 1), 'Color', current_color, 'LineWidth', 1, 'DisplayName', ['ROI ', num2str(i)]);
    end

    % Comment removed after encoding repair.
    bg_bw = (background_mask == i);
    bg_boundaries = bwboundaries(bg_bw);
    for k = 1:length(bg_boundaries)
        boundary = bg_boundaries{k};
        plot(boundary(:, 2), boundary(:, 1), ':', 'Color', current_color, 'LineWidth', 1);
    end

    % Encoding-corrupted comment removed.
    stats = regionprops(roi_bw, 'Centroid');
    if ~isempty(stats)
        text(stats(1).Centroid(1), stats(1).Centroid(2), num2str(i), ...
            'Color', current_color, 'FontWeight', 'bold', 'FontSize', 12, 'HorizontalAlignment', 'center');
    end
end
title('Image: ROI (Solid) & Background (Dotted)');

% Encoding-corrupted comment removed.
subplot(1, 2, 2);
hold on;

% Encoding-corrupted comment removed.
p1 = plot(nan, nan, 'Color', [0.7 0.7 0.7]);
p2 = plot(nan, nan, 'k--', 'LineWidth', 1);
p3 = plot(nan, nan, 'k', 'LineWidth', 1);

for i = 1:num_rois
    current_color = colors(i, :);

    % Comment removed after encoding repair.
    raw_signal = traces_bgfitcorr(:, i) + background_fitted(:, i);

    % Encoding-corrupted comment removed.
    plot(t, raw_signal, 'Color', [current_color, 0.3], 'LineWidth', 0.5);

    % Comment removed after encoding repair.
    plot(t, background_fitted(:, i), '--', 'Color', current_color, 'LineWidth', 1);

    % Encoding-corrupted comment removed.
    plot(t, traces_bgfitcorr(:, i), 'Color', current_color, 'LineWidth', 1);
end

xlabel('Times(s)');
ylabel('Intensity');
title('Traces: Raw (Faded), BG (Dash), Corrected (Solid)');
grid on;

% Comment removed after encoding repair.
legend([p1, p2, p3], {'Raw Signal', 'Fitted Background', 'Corrected Signal'}, 'Location', 'northeast');

% Comment removed after encoding repair.
saveas(gcf, fullfile(save_path, 'background_correction_summary.png'));
fprintf('Summary plot saved to: %s\n', save_path);
end
%% Bleaching Removal
% bleaching_removed is defined as the residual after baseline removal.

if ~exist('bleachmode', 'var') || isempty(bleachmode)
    bleachmode = 'linear';% 'linear' 'highpass' 'exp2'
end

fprintf('Correcting Bleaching (Mode: %s)...\n', bleachmode);

switch bleachmode
    case 'highpass'
        fc = 0.5/t(end);
        [traces_bleaching_removed, baseline] = highpass_bleach_remove(traces_background_removed, freq, fc);

    case 'linear'
        % Encoding-corrupted comment removed.
        traces_bleaching_removed = detrend(traces_background_removed, 1);
        baseline = traces_background_removed - traces_bleaching_removed;

    case 'exp2'
        % Encoding-corrupted comment removed.
        [~, baseline] = fit_exp2(traces_background_removed);
        traces_bleaching_removed = traces_background_removed - baseline;
end

% Comment removed after encoding repair.
plot_corrected(traces_bleaching_removed, traces_background_removed, baseline, t, save_path);

switch bleachmode
    case 'highpass'
        bleaching_parameters = struct('fc', fc);
    otherwise
        bleaching_parameters = struct();
end

bleach_results_file = fullfile(save_path, '2_bleach_results.mat');
bleach_record = build_section_record( ...
    'Bleaching Removal', ...
    'Remove the slow bleaching baseline from the background-corrected traces and save both the residual trace and fitted baseline.', ...
    struct( ...
        'parent_stage', "bg_removed", ...
        'traces_bg_removed', traces_background_removed, ...
        'time_axis', t, ...
        'frame_rate', freq), ...
    struct( ...
        'mode', bleachmode, ...
        'parameters', bleaching_parameters, ...
        'formula', describe_bleach_formula(bleachmode)), ...
    struct( ...
        'traces_bleach_removed', traces_bleaching_removed, ...
        'baseline', baseline), ...
    struct(), ...
    'To rerun this section, load the saved bg_removed traces and parameters from this record, then apply the same bleach-removal mode.');
save(bleach_results_file, 'traces_bleaching_removed', 'baseline', 'bleachmode', 'bleaching_parameters', 'bleach_record');

trace_results = store_trace_stage( ...
    trace_results, 'bleach_removed', traces_bleaching_removed, {'bg_removed'}, roi_results_file, ...
    movie_info, bleachmode, bleaching_parameters);
trace_results = store_trace_stage( ...
    trace_results, 'baseline', baseline, {'bg_removed'}, roi_results_file, ...
    movie_info, bleachmode, bleaching_parameters);
save(trace_results_path, 'trace_results', '-v7.3');

fprintf('Finished Bleaching Correction.\n');
%{
function plot_corrected(traces_corrected, traces, baseline, t, save_path)
% Comment removed after encoding repair.
nrois = size(traces_corrected, 2);

% Comment removed after encoding repair.
fig = figure('Name', 'Bleaching Correction Overview', 'Color', 'w');
set(fig, 'Position', get(0, 'Screensize'));

% Encoding-corrupted comment removed.
ax_fit = subplot(1, 2, 1); hold on;
% Comment removed after encoding repair.
spacing_raw = mean(std(traces, 0, 1)) * 5;

for r = 1:nrois
    offset = (nrois - r) * spacing_raw;
    % Comment removed after encoding repair.
    plot(t, traces(:, r) + offset, 'LineWidth', 0.5, 'DisplayName', 'Original');
    % Comment removed after encoding repair.
    plot(t, baseline(:, r) + offset, 'r', 'LineWidth', 1.2, 'DisplayName', 'Baseline');

    if mod(r, 5) == 0 || r == 1 || r == nrois
        text(t(1), offset, [' ROI ', num2str(r)], 'FontSize', 8, 'FontWeight', 'bold');
    end
end
title(['Original Traces & Fitted Baselines (N=', num2str(nrois), ')']);
xlabel('Time (s)'); ylabel('Stacked Magnitude');
grid on; axis tight;

% Encoding-corrupted comment removed.
ax_corr = subplot(1, 2, 2); hold on;
spacing_corr = mean(std(traces_corrected, 0, 1)) * 8;

for r = 1:nrois
    offset = (nrois - r) * spacing_corr;
    % Comment removed after encoding repair.
    plot(t, traces_corrected(:, r) + offset, 'k', 'LineWidth', 0.5);
    % Encoding-corrupted comment removed.
    line([t(1) t(end)], [offset offset], 'Color', [0.3 0.7 1, 0.5], 'LineStyle', '--');
end
title('Corrected Traces (Residuals)');
xlabel('Time (s)'); ylabel('Stacked Magnitude');
grid on; axis tight;

% Comment removed after encoding repair.
linkaxes([ax_fit, ax_corr], 'x');

% Comment removed after encoding repair.
fig_filename = fullfile(save_path, '2_bleach_correction_stacked.fig');
png_filename = fullfile(save_path, '2_bleach_correction_stacked.png');
saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');
end
%}
% %% Bleaching Correction
%
% bleachmode = 'linear';% 'linear' 'highpass' 'exp2'
%
%      plot_corrected1(traces_corrected, traces, baseline, save_path)
%         % plot
%         fig = figure();
%         set(fig,'Position',get(0,'Screensize'));
%         fit_axe = subplot(1,3,1);
%         fited_axe = subplot(1,3,2);
%         baseline_axe = subplot(1,3,3);
%
%         % plot
%         for i = 1: size(traces_corrected,2)
%             plot(traces(:,i),'Parent',fit_axe);
%             hold(fit_axe, 'on');
%             plot(traces_corrected(:,i),'Parent',fited_axe);
%             hold(fited_axe, 'on');
%
%             plot( baseline(:,i),'LineWidth',1,'Parent',baseline_axe);
%             hold(baseline_axe, 'on');
%
%
%         end
%         hold off;
%
%         % note
%         title(fit_axe, 'Original and Fitted Curves');
%         title(fited_axe, 'Corrected Traces');
%         ylim(fited_axe,[-200,+200])
%         title(baseline_axe, 'Baseline');
%         legend(fit_axe, 'Original Trace', 'Fitted Curve');
%
%         fig_filename = fullfile(save_path, '2_fitted_trace.fig');
%         png_filename = fullfile(save_path, '2_fitted_trace.png');
%
%         saveas(gcf, fig_filename, 'fig');
%         saveas(gcf, png_filename, 'png');
%     end
%
% fprintf('Correcting Bleaching...\n')
%
% %highpass bleach
% switch bleachmode
%     case 'highpass'
%         fc = 15/t(end);%15
%         [traces_corrected,baseline] = highpass_bleach_remove(traces_bgfitcorr,freq,fc );
%     case 'linear'
%         traces_corrected = detrend(traces_bgfitcorr,1 );
%         baseline = traces_bgfitcorr-traces_corrected;
%     case 'exp2'
%         % [traces_corrected, fitted_curves, params, gof] = fit_exp2(traces_bgfitcorr);
%         [traces_corrected,fit_y_fit] = fit_bleach(traces_bgfitcorr,x);
%         baseline = traces_bgfitcorr-traces_corrected;
% end
%
%
% plot_corrected(traces_corrected, traces, baseline, save_path);
% fprintf('Finished\n');
%% Wavelet Denoising
% Denoising is applied on top of bleaching_removed traces.

% Comment removed after encoding repair.
fprintf('Wavelet Denoising...\n')
if ~exist('Dnmethods', 'var') || isempty(Dnmethods)
    Dnmethods = 'FDR';
end
if ~exist('Dnlevel', 'var') || isempty(Dnlevel)
    Dnlevel = 8;
end
if ~exist('Wavename', 'var') || isempty(Wavename)
    Wavename = 'bior6.8';
end
traces_denoised = wdenoise(traces_bleaching_removed, Dnlevel ,DenoisingMethod=Dnmethods,Wavelet = Wavename);

traces_filename = fullfile(save_path, '2_processed_traces.mat');
save(traces_filename,"traces_denoised", "traces_raw", 'traces_bleaching_removed', 'baseline', 'Dnmethods','Dnlevel')
denoise_parameters = struct( ...
    'level', Dnlevel, ...
    'denoising_method', Dnmethods, ...
    'wavelet_name', Wavename);
denoise_record = build_section_record( ...
    'Wavelet Denoising', ...
    'Apply wavelet denoising to the bleach-corrected traces. This denoised trace is also saved as noise_reference for later SNR calculation.', ...
    struct( ...
        'parent_stage', "bleach_removed", ...
        'source_file', bleach_results_file), ...
    denoise_parameters, ...
    struct( ...
        'saved_variables', ["traces_denoised"]), ...
    struct(), ...
    'To rerun this section, load the saved bleach_removed traces and apply wdenoise with the saved parameters.');
save(fullfile(save_path, '2_denoise_results.mat'), ...
    'traces_denoised', 'denoise_parameters', 'denoise_record');
trace_results = store_trace_stage( ...
    trace_results, 'denoised', traces_denoised, {'bleach_removed'}, roi_results_file, ...
    movie_info, 'wdenoise', denoise_parameters);
trace_results = store_trace_stage( ...
    trace_results, 'noise_reference', traces_denoised, {'bleach_removed'}, roi_results_file, ...
    movie_info, 'wdenoise', denoise_parameters);
save(trace_results_path, 'trace_results', '-v7.3');
fprintf('Finished\n');

% Encoding-corrupted comment removed.
fprintf('Plotting all ROIs comparison...\n');

nrois = size(traces_denoised, 2);
t = (1:size(traces_denoised, 1))';


% Encoding-corrupted comment removed.

figure('Name', 'Stacked ROI Comparison', 'Color', 'w', 'Position', [150, 150, 1000, 900]);
max_display = min(10, nrois);
roi_to_show = round(linspace(1, nrois, max_display));

spacing = max(traces_bleaching_removed(:)) * 0.8;

hold on;
for i = 1:length(roi_to_show)
    r_idx = roi_to_show(i);
    offset = (i-1) * spacing;

    % Comment removed after encoding repair.
    plot(t, traces_bleaching_removed(:, r_idx) + offset, 'Color', [0.7 0.7 0.7], 'HandleVisibility', 'off');
    % Comment removed after encoding repair.
    plot(t, traces_denoised(:, r_idx) + offset, 'LineWidth', 1);
end

set(gca, 'YTick', (0:max_display-1)*spacing, 'YTickLabel', string(roi_to_show));
title('Comparison of Selected ROIs (Stacked View)');
xlabel('Time Points'); ylabel('ROI Index');
legend('Denoised Signal');
grid on; hold off;

fprintf('Finished plotting all requested ROI comparisons.\n');
%% Sensitivity And SNR
% Derived metrics are saved into trace_results with explicit lineage.

traces_sensitivity = traces_bleaching_removed./baseline;
noise = traces_bleaching_removed-traces_denoised;
traces_SNR = traces_bleaching_removed./std(noise);

traces_filename = fullfile(save_path, '3_calculated_traces.mat');
save(traces_filename,"traces_sensitivity", 'traces_bleaching_removed', 'baseline', 'noise','traces_SNR')
metric_record = build_section_record( ...
    'Sensitivity And SNR', ...
    'Use the bleach-corrected trace, fitted baseline, and wavelet noise reference to compute noise, sensitivity, and SNR traces.', ...
    struct( ...
        'parent_stage', "bleach_removed", ...
        'noise_reference_stage', "noise_reference", ...
        'source_files', ["2_bleach_results.mat", "2_denoise_results.mat"]), ...
    struct( ...
        'noise_formula', 'bleach_removed - noise_reference', ...
        'sensitivity_formula', 'bleach_removed ./ baseline', ...
        'snr_formula', 'bleach_removed ./ std(noise)'), ...
    struct( ...
        'saved_variables', ["noise", "traces_sensitivity", "traces_SNR"]), ...
    struct(), ...
    'To rerun this section, load the saved bleach_removed / baseline / noise_reference traces and apply the formulas recorded here.');
save(fullfile(save_path, '3_metric_results.mat'), ...
    'traces_sensitivity', 'noise', 'traces_SNR', 'metric_record');
trace_results = store_trace_stage( ...
    trace_results, 'noise', noise, {'bleach_removed', 'noise_reference'}, roi_results_file, ...
    movie_info, 'residual', struct('expression', 'bleach_removed - noise_reference'));
trace_results = store_trace_stage( ...
    trace_results, 'sensitivity', traces_sensitivity, {'bleach_removed', 'baseline'}, roi_results_file, ...
    movie_info, 'ratio', struct('expression', 'bleach_removed ./ baseline'));
trace_results = store_trace_stage( ...
    trace_results, 'snr', traces_SNR, {'bleach_removed', 'noise'}, roi_results_file, ...
    movie_info, 'signal_divided_by_noise_std', struct('noise_reference_stage', 'noise_reference'));
save(trace_results_path, 'trace_results', '-v7.3');
figure()
%plot raw
subplot(1,3,1);
title('raw');
hold on;
[~] = offset_plot(traces_raw,t);


% plot sensitivity
subplot(1,3,2);
title('Sensitivity');
hold on;
[~] = offset_plot(traces_sensitivity,t);

% plot SNR
subplot(1,3,3);
title('SNR');
hold on;
[~] = offset_plot(traces_SNR,t);

fig_filename = fullfile(save_path, '4_SNR.fig');
png_filename = fullfile(save_path, '4_SNR.png');
trace_filename = fullfile(save_path, '4_SNR.mat');

saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');
%% Backend-Specific Peak Detection
% The unified sensitivity peak workflow below runs after both backends have
% produced trace_results.sensitivity.

elseif ~analysis_only_mode && string(analysis_backend) == "volpy"
fprintf('Initializing VolPy backend...\n');
if strlength(string(volpy_input_file)) == 0
    volpy_input_file = string(file_path);
end
if strlength(string(volpy_work_dir)) == 0
    volpy_work_dir = string(fullfile(save_path, '0_volpy_backend'));
end
if isempty(volpy_frame_rate)
    volpy_frame_rate = freq;
end
if strlength(string(volpy_result_mat)) == 0
    volpy_result_mat = string(fullfile(volpy_work_dir, 'volpy_results.mat'));
end

volpy_backend_info = run_volpy_backend_pipeline( ...
    char(volpy_input_file), save_path, volpy_frame_rate, logical(volpy_flip_signal), ...
    logical(volpy_skip_motion_correction), char(volpy_roi_mask_file), ...
    logical(volpy_auto_run), logical(volpy_force_rerun), logical(volpy_use_existing), ...
    char(volpy_work_dir), char(volpy_result_mat), char(volpy_temp_dir), volpy_size_min, volpy_size_max, ...
    volpy_corr_window_seconds, volpy_corr_stride_seconds, volpy_corr_baseline_seconds, ...
    logical(volpy_corr_remove_baseline), logical(volpy_corr_gaussian_blur), ...
    char(volpy_summary_mode), char(volpy_mrcnn_input_mode), volpy_mrcnn_confidence_threshold);
volpy_data = load_volpy_backend_results(char(volpy_backend_info.result_mat));

mask = volpy_data.mask;
rois = struct('bwmask', mask);
avg_image = volpy_data.avg_image;
traces_raw = volpy_data.t;
traces = traces_raw;
nrois = size(traces_raw, 2);

movie_info.motion.applied = true;
movie_info.motion.method = 'VolPy frontend';
movie_info.motion.shift_file = char(volpy_backend_info.motion_corrected_file);
movie_info.motion.parameter_file = char(volpy_backend_info.result_mat);
movie_info.updated_at = datetime("now");
save(movie_info_path, 'movie_info');

roi_results_file = fullfile(save_path, '1_raw_ROI.mat');
roi_record = build_section_record( ...
    'VolPy Frontend Import', ...
    'Import VolPy-generated ROI masks and trace outputs so AP_analysis3 can reuse the downstream AP event and summary stages.', ...
    struct( ...
        'backend_result_file', string(volpy_backend_info.result_mat), ...
        'backend_work_dir', string(volpy_backend_info.output_dir), ...
        'input_movie', string(volpy_backend_info.input_movie_path)), ...
    struct( ...
        'frame_rate', volpy_frame_rate, ...
        'flip_signal', logical(volpy_flip_signal)), ...
    struct( ...
        'nrois', nrois, ...
        'mask', mask), ...
    struct( ...
        'movie', 'not saved'), ...
    'To rerun this section, rerun the VolPy backend launcher or point AP_analysis3 to an existing volpy_results.mat file.');
save(roi_results_file, 'rois', 'avg_image', 'traces', 'roi_record');

trace_results = store_trace_stage( ...
    trace_results, 'raw', traces_raw, {}, roi_results_file, ...
    movie_info, 'volpy_import', struct('source', char(volpy_backend_info.result_mat)));
trace_results = store_trace_stage( ...
    trace_results, 'bleach_removed', volpy_data.event_trace, {'raw'}, roi_results_file, ...
    movie_info, 'volpy_event_trace_compat', struct('source_stage', 'volpy_t_native'));
trace_results = store_trace_stage( ...
    trace_results, 'volpy_t', traces_raw, {'raw'}, roi_results_file, ...
    movie_info, 'volpy_trace_native_orientation', struct('source', char(volpy_backend_info.result_mat)));
trace_results = store_trace_stage( ...
    trace_results, 'volpy_ts', volpy_data.ts, {'volpy_t'}, roi_results_file, ...
    movie_info, 'volpy_matched_filter_trace_native_orientation', struct());
trace_results = store_trace_stage( ...
    trace_results, 'volpy_t_rec', volpy_data.t_rec, {'volpy_t'}, roi_results_file, ...
    movie_info, 'volpy_reconstructed_spike_trace_native_orientation', struct());
trace_results = store_trace_stage( ...
    trace_results, 'volpy_subthreshold', volpy_data.t_sub, {'volpy_t'}, roi_results_file, ...
    movie_info, 'volpy_subthreshold_native_orientation', struct());
trace_results = store_trace_stage( ...
    trace_results, 'volpy_dff', volpy_data.dff, {'volpy_t'}, roi_results_file, ...
    movie_info, 'volpy_dff_native_orientation', struct());
trace_results = store_trace_stage( ...
    trace_results, 'volpy_event_trace', volpy_data.event_trace, {'volpy_t'}, roi_results_file, ...
    movie_info, 'volpy_native_trace', struct());
trace_results = store_trace_stage( ...
    trace_results, 'volpy_event_sensitivity', volpy_data.event_sensitivity, {'volpy_t'}, roi_results_file, ...
    movie_info, 'volpy_native_dff', struct());
trace_results = store_trace_stage( ...
    trace_results, 'volpy_qc_event_trace', volpy_data.qc_event_trace, {'volpy_t_rec'}, roi_results_file, ...
    movie_info, 'volpy_spike_isolated_trace_for_ap_qc', struct());
trace_results = store_trace_stage( ...
    trace_results, 'volpy_qc_event_sensitivity', volpy_data.qc_event_sensitivity, {'volpy_t_rec'}, roi_results_file, ...
    movie_info, 'volpy_spike_isolated_sensitivity_for_ap_qc', struct());
trace_results = store_trace_stage( ...
    trace_results, 'volpy_qc_snr', volpy_data.qc_snr_trace, {'volpy_qc_event_trace', 'volpy_noise'}, roi_results_file, ...
    movie_info, 'volpy_spike_isolated_trace_divided_by_residual_std', struct());
trace_results = store_trace_stage( ...
    trace_results, 'volpy_noise', volpy_data.noise, {'volpy_t', 'volpy_t_rec'}, roi_results_file, ...
    movie_info, 'volpy_residual_native_orientation', struct('expression', 'volpy_t - volpy_t_rec'));
trace_results = store_trace_stage( ...
    trace_results, 'sensitivity', volpy_data.event_sensitivity, {'volpy_event_trace'}, roi_results_file, ...
    movie_info, 'volpy_native_dff_import', struct());
trace_results = store_trace_stage( ...
    trace_results, 'snr', volpy_data.snr_trace, {'volpy_event_trace', 'volpy_noise'}, roi_results_file, ...
    movie_info, 'volpy_native_trace_divided_by_residual_std', struct());
save(trace_results_path, 'trace_results', '-v7.3');

traces_bleaching_removed = volpy_data.event_trace;
traces_sensitivity = volpy_data.event_sensitivity;
traces_SNR = volpy_data.snr_trace;
noise = volpy_data.noise;
ap_qc_trace = volpy_data.qc_event_trace;
ap_qc_sensitivity = volpy_data.qc_event_sensitivity;
ap_qc_SNR = volpy_data.qc_snr_trace;
ap_qc_baseline_subtract = true;
peak_trace_result = 'volpy_event_trace';

peaks_index = volpy_data.spikes;
peaks_amplitude = volpy_data.peak_amplitude;
peaks_polarity = volpy_data.peaks_polarity;

peak_results.volpy_detected.data = struct( ...
    'index', {peaks_index}, ...
    'amplitude', {peaks_amplitude}, ...
    'polarity', {peaks_polarity}, ...
    'snr', volpy_data.snr_scalar, ...
    'locality', volpy_data.locality, ...
    'low_spikes', volpy_data.low_spikes);
peak_results.volpy_detected.info = struct( ...
    'result_name', 'volpy_detected', ...
    'parent_results', {{}}, ...
    'trace_result', peak_trace_result, ...
    'method', 'volpy_spike_extraction', ...
    'parameters', struct( ...
        'backend_result_file', string(volpy_backend_info.result_mat), ...
        'flip_signal', logical(volpy_flip_signal)), ...
    'created_at', datetime("now"));
peak_results.current_result = 'volpy_detected';
save(fullfile(save_path, '4_volpy_peak_import_results.mat'), 'peaks_index', 'peaks_amplitude', 'peaks_polarity', 'volpy_backend_info');
save(peak_results_path, 'peak_results', '-v7.3');
save_volpy_backend_compat_outputs( ...
    save_path, avg_image, mask, traces_raw, volpy_data.t_sub, volpy_data.t_rec, ...
    traces_sensitivity, traces_SNR, t, peaks_index);
end

%% Sensitivity Peak Detection And Manual Editing
% Adjustable detection/edit parameters. Prominence uses the AP2-style
% relative trace-amplitude factor by default.
if string(analysis_backend) ~= "volpy"
peak_results = struct();
if ~exist('peak_polarity_mode', 'var') || isempty(peak_polarity_mode)
    peak_polarity_mode = "auto"; % "auto", "positive", or "negative"
else
    peak_polarity_mode = string(peak_polarity_mode);
end
if ~exist('peak_min_prominence', 'var') || isempty(peak_min_prominence)
    peak_min_prominence = 0.3;
end
if ~exist('peak_min_prominence_mode', 'var') || isempty(peak_min_prominence_mode)
    peak_min_prominence_mode = "relative_factor"; % "relative_factor" or "absolute"
else
    peak_min_prominence_mode = string(peak_min_prominence_mode);
end
if ~exist('peak_min_distance_frames', 'var') || isempty(peak_min_distance_frames)
    peak_min_distance_frames = max(2, round(0.003 * freq));
end
if ~exist('peak_min_height', 'var') || isempty(peak_min_height)
    peak_min_height = 0;
end
if ~exist('run_manual_peak_edit', 'var') || isempty(run_manual_peak_edit)
    run_manual_peak_edit = run_manual_peak_refinement;
end
peak_trace_result = 'sensitivity';
peak_detect_params = struct( ...
    'polarity_mode', peak_polarity_mode, ...
    'min_peak_prominence', peak_min_prominence, ...
    'min_peak_prominence_mode', peak_min_prominence_mode, ...
    'min_peak_distance_frames', peak_min_distance_frames, ...
    'min_peak_height', peak_min_height);
[peak_table_detected, peaks_index, peaks_amplitude, peaks_polarity] = ...
    detect_sensitivity_peaks(traces_sensitivity, freq, peak_detect_params);

peak_results.sensitivity_detected.data = struct( ...
    'peak_table', peak_table_detected, ...
    'index', {peaks_index}, ...
    'amplitude', {peaks_amplitude}, ...
    'polarity', {peaks_polarity});
peak_results.sensitivity_detected.info = struct( ...
    'result_name', 'sensitivity_detected', ...
    'parent_results', {{}}, ...
    'trace_result', peak_trace_result, ...
    'method', 'detect_sensitivity_peaks', ...
    'parameters', peak_detect_params, ...
    'created_at', datetime("now"));
peak_results.current_result = 'sensitivity_detected';
peak_detect_record = build_section_record( ...
    'Sensitivity Peak Detection', ...
    'Run automatic peak finding directly on trace_results.sensitivity and save an event-level peak table.', ...
    struct( ...
        'trace_stage', "sensitivity", ...
        'input_trace', traces_sensitivity), ...
    peak_detect_params, ...
    struct( ...
        'peak_table', peak_table_detected, ...
        'peaks_index', {peaks_index}, ...
        'peaks_amplitude', {peaks_amplitude}, ...
        'peaks_polarity', {peaks_polarity}), ...
    struct(), ...
    'To rerun this section, load trace_results.sensitivity and call detect_sensitivity_peaks with the saved parameters.');
save(fullfile(save_path, '4_sensitivity_peak_detect_results.mat'), ...
    'peak_table_detected', 'peaks_index', 'peaks_amplitude', 'peaks_polarity', ...
    'peak_detect_params', 'peak_detect_record');
save(peak_results_path, 'peak_results', '-v7.3');

edit_history = struct('action', {}, 'roi', {}, 'peak_id', {}, 'index_before', {}, ...
    'index_after', {}, 'mode', {}, 'created_at', {});
manual_edit_params = struct( ...
    'enabled', logical(run_manual_peak_edit), ...
    'min_peak_distance_frames', peak_min_distance_frames, ...
    'frame_rate_hz', freq, ...
    'default_mode', "delete_box", ...
    'keys', "A add box, D delete box, N/space next ROI, R reset ROI, Q finish");
if run_manual_peak_edit
    [peak_table_final, edit_history] = edit_sensitivity_peaks( ...
        traces_sensitivity, peak_table_detected, peaks_polarity, manual_edit_params);
    manual_edit_method = 'manual_add_delete_peak_edit';
else
    peak_table_final = peak_table_detected;
    manual_edit_method = 'manual_add_delete_peak_edit_skipped';
end
[peaks_index, peaks_amplitude, peaks_polarity] = ...
    peak_table_to_peak_cells(peak_table_final, nrois);

peak_results.sensitivity_manually_edited.data = struct( ...
    'peak_table', peak_table_final, ...
    'edit_history', edit_history, ...
    'index', {peaks_index}, ...
    'amplitude', {peaks_amplitude}, ...
    'polarity', {peaks_polarity});
peak_results.sensitivity_manually_edited.info = struct( ...
    'result_name', 'sensitivity_manually_edited', ...
    'parent_results', {{'sensitivity_detected'}}, ...
    'trace_result', peak_trace_result, ...
    'method', manual_edit_method, ...
    'parameters', manual_edit_params, ...
    'created_at', datetime("now"));
peak_results.current_result = 'sensitivity_manually_edited';
if ~isfield(peak_results, 'sensitivity_detected') || ...
        ~isfield(peak_results, 'sensitivity_manually_edited') || ...
        ~isfield(peak_results, 'current_result')
    error('AP_analysis3:PeakResultGenerationFailed', ...
        'Sensitivity peak detection ran but failed to generate peak_results.');
end
peak_edit_record = build_section_record( ...
    'Sensitivity Manual Peak Editing', ...
    ternary(run_manual_peak_edit, ...
        'Interactively add and delete peaks on sensitivity traces while preserving the automatic peak table and edit history.', ...
        'Skip interactive editing and accept the automatically detected sensitivity peaks.'), ...
    struct( ...
        'parent_peak_stage', "sensitivity_detected", ...
        'trace_stage', "sensitivity", ...
        'peak_table_detected', peak_table_detected), ...
    manual_edit_params, ...
    struct( ...
        'peak_table_final', peak_table_final, ...
        'edit_history', edit_history, ...
        'peaks_index', {peaks_index}, ...
        'peaks_amplitude', {peaks_amplitude}, ...
        'peaks_polarity', {peaks_polarity}), ...
    struct(), ...
    'To rerun this section, load the detected peak table and repeat manual add/delete editing on trace_results.sensitivity.');
save(fullfile(save_path, '4_sensitivity_peak_edit_results.mat'), ...
    'peak_table_final', 'edit_history', 'peaks_index', 'peaks_amplitude', ...
    'peaks_polarity', 'manual_edit_params', 'peak_edit_record');
save(peak_results_path, 'peak_results', '-v7.3');
else
    peak_trace_result = 'volpy_event_trace';
    run_manual_peak_edit = false;
end

if isempty(ap_qc_trace)
    ap_qc_trace = traces_bleaching_removed;
end
if isempty(ap_qc_sensitivity)
    ap_qc_sensitivity = traces_sensitivity;
end
if isempty(ap_qc_SNR)
    ap_qc_SNR = traces_SNR;
end


%% AP Events From Accepted Peaks
% Build event-level AP measurements. FWHM is measured, but does not gate or
% redirect accepted peaks.
if ~exist('AP_window_width', 'var') || isempty(AP_window_width)
    AP_window_width = 15; % frames on each side of the peak
end
offset_width = 0;
enable_peak_redirection = false;
% AP_list = AP_statistic(nrois, peaks_index, peaks_amplitude, traces_corrected, traces_sensitivity, traces_SNR, AP_window_width, nframes, dt, peaks_polarity, save_path);

% function [AP_list,peaks_index_corrected]  = AP_statistic(nrois, peaks_index, peaks_amplitude, traces_corrected, traces_sensitivity, traces_SNR, AP_window_width, nframes, dt, peaks_polarity, save_path)
AP_list = cell(1, nrois);

% each trace
for i = 1:nrois % i for trace
    peaks_num = length(peaks_index{i});
    each_trace_amp = ap_qc_trace(:,i);
    each_trace_sensitivity = ap_qc_sensitivity(:,i);
    each_trace_SNR = ap_qc_SNR(:,i);
    AP_list{i} = cell(1, length(peaks_index{i}));
    fprintf('ROI %d processing\n',i);
    j = 1;
    % each peak
    while j <= peaks_num % j for peak

        peak_index_ij = peaks_index{i}(j);
        peak_amp_ij = peaks_amplitude{i}(j);

        % keep in board
        AP_start_index = max(1, peak_index_ij - AP_window_width);
        AP_end_index = min(nframes, peak_index_ij + AP_window_width);
        AP_index = AP_start_index : AP_end_index;

        % search
        AP_amp = each_trace_amp(AP_start_index:AP_end_index)';
        AP_sensitivity = each_trace_sensitivity(AP_start_index:AP_end_index)';
        AP_SNR = each_trace_SNR(AP_start_index:AP_end_index)';

        % fill NaN
        if 1 > peak_index_ij - AP_window_width
            AP_amp = [NaN(1,0 - (peak_index_ij - AP_window_width)+1), AP_amp];
            AP_sensitivity = [NaN(1,0 - (peak_index_ij - AP_window_width)+1),AP_sensitivity];
            AP_SNR = [NaN(1,0 - (peak_index_ij - AP_window_width)+1),AP_SNR];
        elseif nframes < peak_index_ij + AP_window_width
            AP_amp = [AP_amp, NaN(1,peak_index_ij + AP_window_width - nframes)];
            AP_sensitivity = [AP_sensitivity, NaN(1,peak_index_ij + AP_window_width - nframes)];
            AP_SNR = [AP_SNR, NaN(1,peak_index_ij + AP_window_width - nframes)];
        end

        if ap_qc_baseline_subtract
            AP_amp = baseline_subtract_event_window(AP_amp, AP_window_width + 1);
            AP_sensitivity = baseline_subtract_event_window(AP_sensitivity, AP_window_width + 1);
            AP_SNR = baseline_subtract_event_window(AP_SNR, AP_window_width + 1);
        end

        % Calculate;
        Amplitude = abs(peak_amp_ij);
        Sensitivity = AP_sensitivity(AP_window_width+1)*100 ;
        SNR = abs(AP_SNR(AP_window_width+1));
        FWHM = calculate_FWHM(AP_amp, dt, peaks_polarity{i}, false);
        % sprintf('roi % d peaks %d FWHM:%d',i,j,FWHM);

        % save AP data
        each_AP = struct('Trace', i, 'AP_number', j, 'AP_index',AP_index, ...
            'AP_amp',AP_amp,'Amplitude', Amplitude,'FWHM',FWHM, ...
            'AP_sensitivity',AP_sensitivity,'Sensitivity',Sensitivity, ...
            'AP_SNR', AP_SNR, 'SNR', SNR);
        AP_list{i}{j} = each_AP;
        j = j + 1;
    end


end

save(fullfile(save_path,'5_accepted_peaks.mat'),'peaks_index','peaks_polarity','peaks_amplitude')
if ~isfield(peak_results, 'current_result') || isempty(peak_results.current_result)
    error('AP_analysis3:MissingPeakCurrentResult', ...
        'No peak result stage was produced before AP event extraction. Check the sensitivity peak detection section.');
end
event_parent_peak_result = string(peak_results.current_result);
peak_results.accepted_for_events.data = struct( ...
    'index', {peaks_index}, ...
    'amplitude', {peaks_amplitude}, ...
    'polarity', {peaks_polarity});
peak_results.accepted_for_events.info = struct( ...
    'result_name', 'accepted_for_events', ...
    'parent_results', {{char(event_parent_peak_result)}}, ...
    'trace_result', peak_trace_result, ...
    'method', 'accepted_peaks_no_fwhm_gate', ...
    'parameters', struct( ...
        'ap_window_width', AP_window_width, ...
        'fwhm_used_for_gating', false, ...
        'peak_redirection_enabled', false), ...
    'created_at', datetime("now"));
peak_results.current_result = 'accepted_for_events';
save(peak_results_path, 'peak_results', '-v7.3');

ap_results.events.data = AP_list;
ap_results.events.info = struct( ...
    'result_name', 'events', ...
    'parent_peak_result', 'accepted_for_events', ...
    'parent_trace_results', {{'bleach_removed', 'sensitivity', 'snr'}}, ...
    'method', 'ap_window_statistics_without_fwhm_gate', ...
    'parameters', struct( ...
        'ap_window_width', AP_window_width, ...
        'qc_trace_source', ternary(ap_qc_baseline_subtract, 'volpy_qc_event_trace', 'bleach_removed'), ...
        'baseline_subtracted', logical(ap_qc_baseline_subtract), ...
        'fwhm_used_for_gating', false, ...
        'peak_redirection_enabled', false), ...
    'created_at', datetime("now"));
ap_results.current_result = 'events';
ap_event_record = build_section_record( ...
    'AP Events From Accepted Peaks', ...
    'Build event-level AP windows from accepted peaks without FWHM gating or peak redirection.', ...
    struct( ...
        'parent_peak_stage', event_parent_peak_result, ...
        'trace_stage', string(peak_trace_result), ...
        'traces_bleach_removed', traces_bleaching_removed, ...
        'traces_sensitivity', traces_sensitivity, ...
        'traces_snr', traces_SNR), ...
    struct( ...
        'ap_window_width', AP_window_width, ...
        'frame_dt', dt, ...
        'fwhm_used_for_gating', false, ...
        'peak_redirection_enabled', false), ...
    struct( ...
        'peaks_index', {peaks_index}, ...
        'peaks_amplitude', {peaks_amplitude}, ...
        'peaks_polarity', {peaks_polarity}, ...
        'AP_list', AP_list), ...
    struct(), ...
    'To rerun this section, load the accepted peak stage and trace stages, then rerun AP window extraction.');
save(fullfile(save_path, '5_sensitivity_ap_event_results.mat'), 'peaks_index', 'peaks_amplitude', 'peaks_polarity', 'AP_list', 'ap_event_record');
save(ap_results_path, 'ap_results', '-v7.3');

%% Manual Peak Gate Compatibility
% Manual add/delete editing above is the active curation step. Keep a full
% acceptance mask here so legacy AP summary code can run unchanged.

run_manual_peak_gating = exist('legacy_manual_peak_gating', 'var') == 1 && logical(legacy_manual_peak_gating);
if run_manual_peak_gating
peaks_index_manually_gated = peaks_index;

for i = 1:nrois
    figure()
    set(gcf,'Position',[0,0,2500,1800])
    hold on;

    % Encoding-corrupted comment removed.
    plot(traces_SNR(:,i).*peaks_polarity{i});
    peaks_x = peaks_index{i};
    peaks_y = traces_SNR(peaks_x,i).*peaks_polarity{i};
    plot(peaks_x, peaks_y,'v','MarkerFaceColor','r');

    title(sprintf('ROI %d: draw polygon to delete peaks | press R to redraw', i));

    % Comment removed after encoding repair.
    redraw = true;
    while redraw
        % Encoding-corrupted comment removed.
        poly = drawpolygon(gca);
        wait(poly);

        % Encoding-corrupted comment removed.
        xv = poly.Position(:,1);
        yv = poly.Position(:,2);

        % Comment removed after encoding repair.
        [in, ~] = inpolygon(peaks_x, peaks_y, xv, yv);

        % Encoding-corrupted comment removed.
        selected_plot = plot(peaks_x(in), peaks_y(in), ...
            'v', 'MarkerFaceColor', 'g', 'MarkerSize', 10);

        % Comment removed after encoding repair.
        choice = questdlg(sprintf('Selected %d peaks. Confirm deletion?', sum(in)), ...
            'Confirm selection', 'Confirm', 'Redraw(R)', 'Cancel', 'Redraw(R)');

        switch choice
            case 'Confirm'
                manual_gated_index = in;
                redraw = false;
                fprintf('ROI %d: confirmed deletion of %d peaks\n', i, sum(in));

            case 'Redraw(R)'
                % Comment removed after encoding repair.
                delete(poly);
                delete(selected_plot);
                fprintf('ROI %d: redraw selection...\n', i);

            case 'Cancel'
                % Comment removed after encoding repair.
                delete(poly);
                delete(selected_plot);
                fprintf('ROI %d: selection canceled\n', i);
                manual_gated_index = false(size(peaks_x));
                redraw = false;
        end
    end

    % Comment removed after encoding repair.
    peaks_index_manually_gated{i}(manual_gated_index) = NaN;
    peaks_index_manually_gated{i}(~manual_gated_index) = 1;

    % Comment removed after encoding repair.
    mkdir(fullfile(save_path,'Manually gated peaks'));
    saveas(gcf, fullfile(save_path,'Manually gated peaks', ...
        sprintf('noi %d, peaks %d.png', i, length(peaks_x))))
    saveas(gcf, fullfile(save_path,'Manually gated peaks', ...
        sprintf('noi %d, peaks %d.fig', i, length(peaks_x))))
    close(gcf)
end
else
peaks_index_manually_gated = cell(1, nrois);
for i = 1:nrois
    peaks_index_manually_gated{i} = ones(size(peaks_index{i}));
end
end

save(fullfile(save_path,'Manually gated peaks.mat'),'peaks_index_manually_gated')
peaks_index_gated = peaks_index;
peaks_amplitude_gated = peaks_amplitude;
for i = 1:nrois
    peaks_index_gated{i} = peaks_index{i} .* peaks_index_manually_gated{i};
    peaks_amplitude_gated{i} = peaks_amplitude{i} .* peaks_index_manually_gated{i};
end
peak_results.manually_gated.data = struct( ...
    'index_mask', {peaks_index_manually_gated}, ...
    'index', {peaks_index_gated}, ...
    'amplitude', {peaks_amplitude_gated}, ...
    'polarity', {peaks_polarity});
peak_results.manually_gated.info = struct( ...
    'result_name', 'manually_gated', ...
    'parent_results', {{'accepted_for_events'}}, ...
    'trace_result', 'sensitivity', ...
    'method', 'legacy_full_acceptance_mask', ...
    'parameters', struct('skipped', true), ...
    'created_at', datetime("now"));
peak_results.current_result = 'manually_gated';
peak_gate_record = build_section_record( ...
    'Manual Peak Gate Compatibility', ...
    'Skip the legacy polygon gate and accept all peaks from the add/delete editing workflow.', ...
    struct( ...
        'parent_peak_stage', "accepted_for_events", ...
        'trace_stage', "sensitivity", ...
        'traces_sensitivity', traces_sensitivity, ...
        'peaks_index_before_gate', {peaks_index}, ...
        'peaks_polarity', {peaks_polarity}), ...
    struct('skipped', true), ...
    struct( ...
        'peaks_index_manually_gated', {peaks_index_manually_gated}, ...
        'peaks_index_gated', {peaks_index_gated}, ...
        'peaks_amplitude_gated', {peaks_amplitude_gated}), ...
    struct(), ...
    'No rerun is needed; manual peak curation is recorded in sensitivity_manually_edited.');
save(fullfile(save_path, '5_manual_peak_gate_results.mat'), 'peaks_index_manually_gated', 'peaks_index_gated', 'peaks_amplitude_gated', 'peak_gate_record');
save(peak_results_path, 'peak_results', '-v7.3');
%% AP Summary Statistics
% Summarize event-level AP measurements into per-ROI tables and figures.


% Encoding-corrupted comment removed.
avg_FWHM = zeros(length(AP_list), 1);
avg_sensitivity = zeros(length(AP_list), 1);
avg_SNR = zeros(length(AP_list), 1);
AP_number = zeros(length(AP_list), 1);
ROI_number = zeros(length(AP_list), 1);

AP_data.amp = cell(1, nrois);
AP_data.FWHM = cell(1, nrois);
AP_data.sensitivity = cell(1, nrois);
AP_data.SNR = cell(1, nrois);
AP_data.index = cell(1, nrois);
for i = 1:nrois
    AP_data.amp{i} = zeros(0, 2*AP_window_width+1);
    AP_data.FWHM{i} = [];
    AP_data.sensitivity{i} = [];
    AP_data.SNR{i} = [];
    AP_data.index{i} = [];
    ROI_number(i) = i;
    avg_FWHM(i) = NaN;
    avg_sensitivity(i) = NaN;
    avg_SNR(i) = NaN;
end


% Encoding-corrupted comment removed.
table_name = fullfile(save_path,'AP_data.xlsx');
for i = 1:length(AP_list)

    if ~isempty(AP_list{i}) && any(~cellfun('isempty', AP_list{i}))
        AP_i = AP_list{i};

        % Encoding-corrupted comment removed.
        number_i = zeros(length(AP_i), 1);
        amp_i = zeros(length(AP_i), 2*AP_window_width+1);
        FWHM_i = zeros(length(AP_i), 1);
        sensitivity_i = zeros(length(AP_i), 1);
        SNR_i = zeros(length(AP_i), 1);
        index_i = zeros(length(AP_i), 1);

        for j = 1:length(AP_i)
            if ~isempty(peaks_index_manually_gated)
                each_AP = AP_i{j};
                number_i(j) = each_AP.AP_number;
                amp_i(j,:)  = each_AP.AP_amp .* peaks_index_manually_gated{i}(j);
                FWHM_i(j)  = each_AP.FWHM*peaks_index_manually_gated{i}(j);
                sensitivity_i(j)  = each_AP.Sensitivity*peaks_index_manually_gated{i}(j);
                SNR_i(j)  = each_AP.SNR*peaks_index_manually_gated{i}(j);
                index_i(j) = peaks_index{i}(j)*peaks_index_manually_gated{i}(j);
            else
                each_AP = AP_i{j};
                number_i(j) = each_AP.AP_number;
                amp_i(j,:)  = each_AP.AP_amp;
                FWHM_i(j)  = each_AP.FWHM;
                sensitivity_i(j)  = each_AP.Sensitivity;
                SNR_i(j)  = each_AP.SNR;
                index_i(j) = peaks_index{i}(j);
            end
        end

        % Encoding-corrupted comment removed.
        T = table(number_i, amp_i(:,2*AP_window_width+1), FWHM_i, sensitivity_i, SNR_i, index_i, ...
            'VariableNames', {'Number', 'Amplitude', 'FWHM (ms)', 'Sensitivity', 'SNR', 'Index'});

        % Encoding-corrupted comment removed.
        sheet_name = string(['ROI ' num2str(i)]);
        writetable(T,table_name, 'Sheet', sheet_name);

        % save average value
        % avg_amp(i) = mean(amp_i,'omitmissing');
        avg_FWHM(i) = mean(FWHM_i,'omitmissing');
        avg_sensitivity(i) = mean(sensitivity_i,'omitmissing');
        avg_SNR(i) = mean(SNR_i,'omitmissing');
        AP_number(i) = number_i(end) - sum(isnan(peaks_index_manually_gated{i}));
        ROI_number(i) = i;
        AP_data.amp{i} = amp_i;
        AP_data.FWHM{i} = FWHM_i;
        AP_data.sensitivity{i} = sensitivity_i;
        AP_data.SNR{i} = SNR_i;
        AP_data.index{i} = index_i;
    end
end
fprintf('Finished statistic AP\n')
%% Save AP Summary Figures

figure()
FWHM_axe = subplot(1,3,3);hold on;xlim([0,nrois+1]);
sensitivity_axe = subplot(1,3,2);hold on;xlim([0,nrois+1]);
SNR_axe = subplot(1,3,1);hold on;xlim([0,nrois+1]);
sgtitle('AP statistic');

% Comment removed after encoding repair.
allFWHM = [];
alldff = [];
allSNR = [];
Labels = [];

% Comment removed after encoding repair.
for i = 1:nrois
    % Encoding-corrupted comment removed.
    currentFWHM = AP_data.FWHM{i};
    currentdff = AP_data.sensitivity{i};
    currentSNR = AP_data.SNR{i};
    % Comment removed after encoding repair.
    allFWHM  = [allFWHM; currentFWHM(:)];
    alldff  = [alldff; currentdff(:)];
    allSNR  = [allSNR; currentSNR(:)];
    % Encoding-corrupted comment removed.
    Labels = [Labels; i * ones(length(currentFWHM), 1)];
end
boxchart(Labels, allFWHM, 'Parent',FWHM_axe,'MarkerStyle','x','JitterOutliers','on');

xlabel('ROI number','Parent',FWHM_axe);
ylabel('FWHM (ms)','Parent',FWHM_axe);

boxchart(Labels, alldff*-1, 'Parent',sensitivity_axe,'MarkerStyle','x','JitterOutliers','on');hold on;
xlabel('ROI number','Parent',sensitivity_axe);
ylabel('Sensitiviy','Parent',sensitivity_axe);

boxchart(Labels, allSNR,'Parent',SNR_axe,'MarkerStyle','x','JitterOutliers','on');hold on;
xlabel('ROI number','Parent',SNR_axe);
ylabel('SNR','Parent',SNR_axe);

fig_filename = fullfile(save_path, '5_AP statistic.fig');
png_filename = fullfile(save_path, '5_AP statistic.png');

saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');

T_ave = table(ROI_number, AP_number,  avg_FWHM, avg_sensitivity, avg_SNR, ...
    'VariableNames', {'ROI Number','AP Number',  'Average FWHM (ms)', 'Average Sensitivity', 'Average SNR'});
writetable(T_ave, table_name, 'Sheet', 'Average');
AP_data_filename = fullfile(save_path, 'AP_data.mat');
save(AP_data_filename, "AP_data",'AP_list')
ap_results.summary.data = struct( ...
    'AP_data', AP_data, ...
    'AP_list', AP_list, ...
    'ROI_number', ROI_number, ...
    'AP_number', AP_number, ...
    'avg_FWHM', avg_FWHM, ...
    'avg_sensitivity', avg_sensitivity, ...
    'avg_SNR', avg_SNR);
ap_results.summary.info = struct( ...
    'result_name', 'summary', ...
    'parent_peak_result', 'manually_gated', ...
    'parent_trace_results', {{'bleach_removed', 'sensitivity', 'snr'}}, ...
    'method', 'ap_summary_statistics', ...
    'parameters', struct( ...
        'table_name', table_name), ...
    'created_at', datetime("now"));
ap_results.current_result = 'summary';
ap_summary_record = build_section_record( ...
    'AP Summary Statistics', ...
    'Aggregate event-level AP measurements into per-ROI tables and summary statistics.', ...
    struct( ...
        'parent_ap_stage', "events", ...
        'parent_peak_stage', "manually_gated", ...
        'AP_list', AP_list, ...
        'peaks_index_manually_gated', {peaks_index_manually_gated}), ...
    struct( ...
        'table_name', table_name, ...
        'ap_window_width', AP_window_width), ...
    struct( ...
        'AP_data', AP_data, ...
        'ROI_number', ROI_number, ...
        'AP_number', AP_number, ...
        'avg_FWHM', avg_FWHM, ...
        'avg_sensitivity', avg_sensitivity, ...
        'avg_SNR', avg_SNR), ...
    struct(), ...
    'To rerun this section, load the saved AP_list and peak gate mask, then repeat the ROI-wise aggregation and table export.');
save(fullfile(save_path, '6_ap_summary_results.mat'), 'AP_data', 'AP_list', 'ROI_number', 'AP_number', 'avg_FWHM', 'avg_sensitivity', 'avg_SNR', 'ap_summary_record');
save(ap_results_path, 'ap_results', '-v7.3');
fprintf('Finished statistic AP\n')
% end

fprintf('AP_data.xlsx saved.\n');
%% Trend Analysis
nframe = size(traces_SNR,1);
trendlength = floor(nframe/freq);
trendbin = 5;
trendpart = floor(trendlength/trendbin);
trend_avgSNR = NaN(1,trendpart);
trend_stdSNR = NaN(1,trendpart);
all_avgSNR = NaN(nrois,trendpart);

trend_avgdff = NaN(1,trendpart);
trend_stddff = NaN(1,trendpart);
all_avgdff = NaN(nrois,trendpart);

trend_avgFWHM = NaN(1,trendpart);
trend_stdFWHM = NaN(1,trendpart);
all_avgFWHM = NaN(nrois,trendpart);

trend_avgFR = NaN(1,trendpart);
trend_stdFR = NaN(1,trendpart);
all_avgFR = NaN(nrois,trendpart);

trend_ISI = cell(nrois,trendpart);
all_ISI = cell(nrois,1);

for t = 1:trendpart
    startindex = (t-1) * trendbin * freq + 1;
    endindex = min(nframe, t * trendbin * freq);
    current_avgSNR = NaN(1,nrois);
    current_avgdff = NaN(1,nrois);
    current_avgFWHM = NaN(1,nrois);
    current_avgFR = NaN(1,nrois);

    % current_stdSNR = zeros(1,24);
    for i = 1:nrois
        indice = find((AP_data.index{i} >= startindex) & (AP_data.index{i} <= endindex));
        current_avgSNR(i) = mean(AP_data.SNR{i}(indice),'omitmissing');
        % current_stdSNR(i) = std(AP_data.SNR{i}(indice),'omitmissing');
        all_avgSNR(i,t) = current_avgSNR(i);

        current_avgdff(i) = mean(AP_data.sensitivity{i}(indice),'omitmissing');
        all_avgdff(i,t) = current_avgdff(i);

        current_avgFWHM(i) = mean(AP_data.FWHM{i}(indice),'omitmissing');
        all_avgFWHM(i,t) = current_avgFWHM(i);

        current_avgFR(i) = numel(indice) / trendbin;
        all_avgFR(i,t) = current_avgFR(i);


        all_ISI{i} = diff(AP_data.index{i})/freq;
        trend_ISI{i,t} = all_ISI{i}(indice(indice<=length(all_ISI{i})));
    end

    trend_avgSNR(t) = mean(current_avgSNR,'omitmissing');
    trend_stdSNR(t) = std(current_avgSNR,'omitmissing');

    trend_avgdff(t) = mean(current_avgdff,'omitmissing');
    trend_stddff(t) = std(current_avgdff,'omitmissing');

    trend_avgFWHM(t) = mean(current_avgFWHM,'omitmissing');
    trend_stdFWHM(t) = std(current_avgFWHM,'omitmissing');

    trend_avgFR(t) = mean(current_avgFR,'omitmissing');
    trend_stdFR(t) = std(current_avgFR,'omitmissing');


end

figure('Position', [100, 100, 900, 1100])
trendx = (1:trendpart) * trendbin;
subplot(3,1,1)
title('SNR');hold on;
% for i = 1:nrois
%     plot(trendx, all_avgSNR(i,:),'Color',[0.8,0.8,0.8])
% end
plot_trend_band(trendx, trend_avgSNR, trend_stdSNR, [0.70 0.78 0.92]);

xlabel('Time (s)')
ylabel('SNR')

subplot(3,1,2)
title('Sensitivity');hold on;
% for i = 1:nrois
%     plot(trendx, all_avgdff(i,:)*-1,'Color',[0.8,0.8,0.8])
% end
plot_trend_band(trendx, trend_avgdff*-1, abs(trend_stddff), [0.92 0.74 0.74]);

xlabel('Time (s)')
ylabel('Sensitivity (-%)')


subplot(3,1,3)
title('FWHM');hold on;
% for i = 1:nrois
%     plot(trendx, all_avgFWHM(i,:)*1000,'Color',[0.8,0.8,0.8])
% end
plot_trend_band(trendx, trend_avgFWHM, trend_stdFWHM, [0.74 0.86 0.74]);

xlabel('Time (s)')
ylabel('FWHM (ms)')

saveas(gcf,fullfile(save_path,'Trend Analysis.png'))
saveas(gcf,fullfile(save_path,'Trend Analysis.fig'))

save(fullfile(save_path,'Trend Analysis.mat'),'all_avgdff', 'all_avgFR', 'all_avgFWHM', 'all_avgSNR', ...
    'trend_avgdff', 'trend_avgFR', 'trend_avgFWHM', 'trend_avgSNR', 'trend_avgSNR', 'trend_stddff', 'trend_stdFR', 'trend_stdFWHM', 'trend_stdSNR', ...
    'trendx', 'trendbin', 'trendpart')

%% Average AP Plots
% Plot average AP waveforms using sensitivity and SNR aligned to detected events.
figure();
% Comment removed after encoding repair.
plot_cols = sum(cellfun('isempty',AP_list)==0)+1;
plot_col = 0;
% Initialize arrays to store all traces for final average calculation
trace_AP_mean = [];

for i = 1:nrois % i for trace
    % Comment removed after encoding repair.
    peaks_num = length(peaks_index{i});
    if cellfun(['isempt' ...
            'y'],AP_list{i}) == 0
        plot_col = plot_col + 1;
        subplot(2,ceil(plot_cols/2),plot_col);
        set(gca,'color','none');

        % set y axis direction
        if peaks_polarity{i} == -1
            set(gca,'YDir','reverse')
            hold on;
        end

        % get each AP
        AP_i = zeros(peaks_num, AP_window_width*2+1);
        for j = 1:peaks_num
            each_AP = AP_list{i}{j};
            AP_i(j,:) = each_AP.AP_sensitivity;
            % plot each AP
            % plot((1:AP_window_width*2+1)*dt, each_AP.AP_sensitivity','Color',[0.8 0.8 0.8]);
            hold on;
        end

        AP_mean = mean(AP_i, 1, 'omitnan');
        AP_sd = std(AP_i, 0, 1, 'omitnan');

        % plot average AP for each trace
        subplot(2,ceil(plot_cols/2),plot_col);
        current_color = colors(mod(i - 1, size(colors, 1)) + 1, :);
        plot((1:AP_window_width*2+1)*dt, AP_mean,'Color',current_color,'LineWidth',2);
        hold on;
        title(sprintf('ROI %d\n',i));
        fill([(1:AP_window_width*2+1)*dt, fliplr((1:AP_window_width*2+1)*dt)], ...
            [AP_mean + AP_sd, fliplr(AP_mean - AP_sd)], ...
            current_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        hold on;

        % plot average AP for average trace

        % plot((1:AP_window_width*2+1)*dt, mean(AP_i,1,'omitnan'),'Color',colors(i,:),'LineWidth',1);
        % hold on;
        trace_AP_mean = [trace_AP_mean; AP_mean*peaks_polarity{i}];
    end
end

if ~isempty(trace_AP_mean)
    overall_mean = mean(trace_AP_mean, 1, 'omitnan');
    overall_sem = std(trace_AP_mean, 0, 1, 'omitnan') / sum(~cellfun('isempty', AP_list));
    subplot(2, ceil(plot_cols / 2), plot_cols);set(gca,'color','none');hold on;
    fill([(1:AP_window_width*2+1)*dt, fliplr((1:AP_window_width*2+1)*dt)], ...
        [overall_mean + overall_sem, fliplr(overall_mean - overall_sem)], ...
        [0.8 0.8 0.8], 'EdgeColor', 'none');

    title('Averaged of All');
    plot((1:AP_window_width*2+1) * dt, overall_mean, 'Color', 'k', 'LineWidth', 1);
    hold on;
end

sgtitle('Averaged Sensitivity');
hold on;
fig_filename = fullfile(save_path, '6_average_AP_sensitivity.fig');
png_filename = fullfile(save_path, '6_average_AP_sensitivity.png');
mat_filename = fullfile(save_path, '6_average_AP_sensitivity.mat');
saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');
save(mat_filename, 'AP_list', 'peaks_index', 'nrois', 'AP_window_width', 'dt', ...
    'peaks_polarity', 'colors', 'trace_AP_mean', 'overall_mean', 'overall_sem');


% Plot average AP SNR with SD
figure();
plot_cols = sum(cellfun('isempty', AP_list) == 0) + 1;
plot_col = 0;
% Initialize arrays to store all traces for final average calculation
trace_AP_SNR_mean = [];

for i = 1:nrois
    peaks_num = length(peaks_index{i});
    if ~isempty(AP_list{i})
        plot_col = plot_col + 1;
        subplot(2, ceil(plot_cols / 2), plot_col);
        set(gca,'color','none');

        % Set y axis direction if necessary
        if peaks_polarity{i} == -1
            set(gca, 'YDir', 'reverse');
            hold on;
        end

        % Get each AP and compute mean & SEM
        AP_i = zeros(peaks_num, AP_window_width*2 + 1);
        for j = 1:peaks_num
            each_AP = AP_list{i}{j};
            AP_i(j,:) = each_AP.AP_SNR;
            %plot((1:AP_window_width*2+1) * dt, each_AP.AP_SNR', 'Color', [0.8 0.8 0.8]);
            hold on;
        end

        AP_mean = mean(AP_i, 1, 'omitnan');
        AP_sd = std(AP_i, 0, 1, 'omitnan');

        % Plot mean with Sd
        current_color = colors(mod(i - 1, size(colors, 1)) + 1, :);
        plot((1:AP_window_width*2+1) * dt, AP_mean, 'Color', current_color, 'LineWidth', 2);
        fill([(1:AP_window_width*2+1)*dt, fliplr((1:AP_window_width*2+1)*dt)], ...
            [AP_mean + AP_sd, fliplr(AP_mean - AP_sd)], ...
            current_color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        title(sprintf('ROI %d', i));
        hold on;

        % % Plot average AP for all traces
        % subplot(2, ceil(plot_cols / 2), plot_cols);
        % plot((1:AP_window_width*2+1) * dt, AP_mean*peaks_polarity(i), 'Color', [0.8 0.8 0.8]);
        % hold on;

        % Collect traces for final average calculation
        trace_AP_SNR_mean = [trace_AP_SNR_mean; AP_mean*peaks_polarity{i}];
    end

end

% Plot overall average and SEM in the last subplot
if ~isempty(trace_AP_SNR_mean)
    overall_mean = mean(trace_AP_SNR_mean, 1, 'omitnan');
    overall_sem = std(trace_AP_SNR_mean, 0, 1, 'omitnan') / sum(~cellfun('isempty', AP_list));
    subplot(2, ceil(plot_cols / 2), plot_cols);hold on;
    set(gca,'color','none');
    fill([(1:AP_window_width*2+1)*dt, fliplr((1:AP_window_width*2+1)*dt)], ...
        [overall_mean + overall_sem, fliplr(overall_mean - overall_sem)], ...
        [0.8 0.8 0.8], 'EdgeColor', 'none');

    title('Averaged of All');
    plot((1:AP_window_width*2+1) * dt, overall_mean, 'Color', 'k', 'LineWidth', 1);
    hold on;
end
sgtitle('Average SNR with SD');
fig_filename = fullfile(save_path, '6_average_AP_SNR_with_SD.fig');
png_filename = fullfile(save_path, '6_average_AP_SNR_with_SD.png');

saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');


% %% drafting (optional)
% % plot SNR
% figure;
% title('SNR');
% hold on;
% [~] = offset_plot(traces_SNR,t);
% 
% cycle_gd = t(end)/32;
% 
% cycles = 32;
% starts = 0;
% for i = 1:cycles/2
%     fill([starts,starts+cycle_gd,starts+cycle_gd,starts],[0,0,sum(max(traces_SNR))*3,sum(max(traces_SNR))*3],'k','FaceAlpha',0.2)
%     starts = starts+cycle_gd*2 ;
% end
% 
% fig_filename = fullfile(save_path, '4_grafting_SNR.fig');
% png_filename = fullfile(save_path, '4_grafting_SNR.png');
% trace_filename = fullfile(save_path, '4_grafting_SNR.mat');
% 
% saveas(gcf, fig_filename, 'fig');
% saveas(gcf, png_filename, 'png');

%% AP Sequence And Density
% Build AP sequence and density summaries from the final AP statistics.
% Comment removed after encoding repair.
valid_rois = [];
roi_ap_counts = [];

for roi = 1:size(traces_SNR, 2)
    if ~isempty(AP_data.index{roi}) && ~all(isnan(AP_data.index{roi}))
        % Comment removed after encoding repair.
        ap_indices = AP_data.index{roi};
        valid_ap_count = sum(~isnan(ap_indices));
        
        if valid_ap_count > 0
            valid_rois = [valid_rois, roi];
            roi_ap_counts = [roi_ap_counts, valid_ap_count];
        end
    end
end

fprintf('AP sequence summary\n');
fprintf('Total ROI count: %d\n', size(traces_SNR, 2));
fprintf('ROIs with APs: %d\n', length(valid_rois));
fprintf('ROIs without APs: %d\n', size(traces_SNR, 2) - length(valid_rois));
fprintf('ROI IDs with APs: %s\n', mat2str(valid_rois));

% Comment removed after encoding repair.
if isempty(valid_rois)
    fprintf('No ROIs with APs were found. Skipping AP sequence plot.\n');
    return;
end

% Comment removed after encoding repair.
figure('Position', [100, 100, 1400, 800]);

% Comment removed after encoding repair.
peaks_img = zeros(size(traces_SNR, 1), length(valid_rois));

% Comment removed after encoding repair.
hold on;
for idx = 1:length(valid_rois)
    roi = valid_rois(idx);
    ap_indices = AP_data.index{roi};
    
    % Encoding-corrupted comment removed.
    valid_ap_indices = ap_indices(~isnan(ap_indices));
    
    if ~isempty(valid_ap_indices)
        % Comment removed after encoding repair.
        h = plot(valid_ap_indices, idx, '|k', 'LineWidth', 1, 'MarkerSize', 8);
        
        % Comment removed after encoding repair.
        for i = 1:length(valid_ap_indices)
            text(valid_ap_indices(i), idx, sprintf('ROI%d', roi), ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
                'FontSize', 8, 'Color', 'blue', 'Visible', 'off');
        end
        
        % Comment removed after encoding repair.
        for i = 1:length(valid_ap_indices)
            index = valid_ap_indices(i);
            if index >= 1 && index <= size(peaks_img, 1)
                peaks_img(index, idx) = 1;
            end
        end
    end
end

% Encoding-corrupted comment removed.
xlabel('Frame', 'FontSize', 12);
ylabel('ROI number', 'FontSize', 12);
title(sprintf('AP sequence (%d ROIs with APs, total AP count %d)', length(valid_rois), sum(roi_ap_counts)), 'FontSize', 14);

% Comment removed after encoding repair.
set(gca, 'YTick', 1:length(valid_rois));
set(gca, 'YTickLabel', arrayfun(@num2str, valid_rois, 'UniformOutput', false));
grid on;
box on;

% Encoding-corrupted comment removed.
xlim([1, size(traces_SNR, 1)]);
ylim([0.5, length(valid_rois)+0.5]);


% Comment removed after encoding repair.
fig_filename = fullfile(save_path, '8_AP_Sequence.fig');
png_filename = fullfile(save_path, '8_AP_Sequence.png');
saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');

% Comment removed after encoding repair.
print(fullfile(save_path, '8_AP_Sequence_highres.png'), '-dpng', '-r300');

% Comment removed after encoding repair.
peaks_img_valid = peaks_img;
save(fullfile(save_path, 'ap_sequence_data.mat'), ...
    'peaks_img_valid', 'valid_rois', 'roi_ap_counts');
ap_results.sequence.data = struct( ...
    'peaks_img_valid', peaks_img_valid, ...
    'valid_rois', valid_rois, ...
    'roi_ap_counts', roi_ap_counts);
ap_results.sequence.info = struct( ...
    'result_name', 'sequence', ...
    'parent_ap_result', 'summary', ...
    'method', 'ap_sequence_visualization', ...
    'parameters', struct(), ...
    'created_at', datetime("now"));
ap_results.current_result = 'sequence';
save(ap_results_path, 'ap_results', '-v7.3');

% Encoding-corrupted comment removed.
figure('Position', [100, 100, 1400, 400]);
ap_density = sum(peaks_img, 2);
plot(1:length(ap_density), ap_density, 'b-', 'LineWidth', 1);
xlabel('Frame', 'FontSize', 12);
ylabel('AP count', 'FontSize', 12);
title('AP density over time', 'FontSize', 14);
grid on;
box on;

% Encoding-corrupted comment removed.
fig_filename = fullfile(save_path, '9_AP_Density.fig');
png_filename = fullfile(save_path, '9_AP_Density.png');
saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');
ap_results.density.data = struct( ...
    'ap_density', ap_density);
ap_results.density.info = struct( ...
    'result_name', 'density', ...
    'parent_ap_result', 'sequence', ...
    'method', 'ap_density_visualization', ...
    'parameters', struct(), ...
    'created_at', datetime("now"));
ap_results.current_result = 'density';
save(ap_results_path, 'ap_results', '-v7.3');

% Encoding-corrupted comment removed.
fprintf('\nAP sequence plot details\n');
fprintf('ROIs shown in figure: %d\n', length(valid_rois));
fprintf('Total AP count: %d\n', sum(roi_ap_counts));
fprintf('Mean AP count per ROI with APs: %.2f\n', mean(roi_ap_counts));
fprintf('AP density (AP/frame): %.4f\n', sum(roi_ap_counts)/size(traces_SNR, 1));
fprintf('Figures saved to:\n');
fprintf('  - %s\n', fullfile(save_path, '8_AP_Sequence.png'));
fprintf('  - %s\n', fullfile(save_path, '9_AP_Density.png'));

%% ISI Statistic Phenotype QC
% Integrated revised ISI / burst / firing phenotype QC.
% This AP3-native section uses variables already present in the workspace.
% It does not reload movie_info, trace_results, or peak_results from disk.
if ~exist('isi_stat_params_override', 'var') || isempty(isi_stat_params_override)
    isi_stat_params_override = struct();
end
params_override = isi_stat_params_override;
params = default_isi_stat_params_revised();
params = merge_struct_fields(params, params_override);
params = finalize_isi_stat_params(params);

% Prepare current AP3 context
recording_duration_s = nframes / freq;
frame_dt_ms = 1000 / freq;

if exist('peaks_index_gated', 'var') == 1 && ~isempty(peaks_index_gated)
    peak_indices = peaks_index_gated;
    peak_source = 'workspace.peaks_index_gated';
elseif exist('peaks_index', 'var') == 1 && ~isempty(peaks_index)
    peak_indices = peaks_index;
    peak_source = 'workspace.peaks_index';
else
    error('AP_analysis3:MissingIsiPeakIndices', ...
        'ISI statistic QC requires peaks_index_gated or peaks_index in the current workspace.');
end

if exist('traces_sensitivity', 'var') == 1 && ~isempty(traces_sensitivity)
    trace_for_qc = traces_sensitivity;
    trace_source = 'workspace.traces_sensitivity';
elseif exist('traces_SNR', 'var') == 1 && ~isempty(traces_SNR)
    trace_for_qc = traces_SNR;
    trace_source = 'workspace.traces_SNR';
else
    trace_for_qc = [];
    trace_source = 'none';
end

nrois = max([numel(peak_indices), size(trace_for_qc, 2)]);
if nrois == 0
    error('AP_analysis3:NoIsiRois', ...
        'No ROI traces or peak indices are available for ISI statistic QC.');
end
if numel(peak_indices) < nrois
    peak_indices(numel(peak_indices)+1:nrois) = {[]};
end

% Convert accepted peak frames to spike timestamps
spike_times_s = cell(1, nrois);
for roi_idx = 1:nrois
    idx = peak_indices{roi_idx};
    idx = idx(isfinite(idx) & idx > 0);
    idx = unique(idx(:)');
    if params.frame_index_origin_is_zero_time
        % MATLAB frame 1 corresponds to t = 0. This does not affect ISI, but
        % it gives correct timing for FR windows and plot overlays.
        t = (idx - 1) ./ freq;
    else
        t = idx ./ freq;
    end
    t = t(t >= 0 & t <= recording_duration_s);
    spike_times_s{roi_idx} = sort(t);
end

output_dir = fullfile(save_path, '10_isi_stat_phenotype_qc_revised');
if ~isfolder(output_dir)
    mkdir(output_dir);
end

% Analyze all neurons
neuron_rows = repmat(empty_neuron_row(), nrois, 1);
burst_event_rows = empty_burst_event_row();
burst_event_rows(1) = [];
per_neuron_metrics = cell(nrois, 1);

for roi_idx = 1:nrois
    t_spikes = sort(spike_times_s{roi_idx}(:)');
    metrics = analyze_one_spontaneous_neuron(t_spikes, recording_duration_s, params, roi_idx, frame_dt_ms);
    per_neuron_metrics{roi_idx} = metrics;
    neuron_rows(roi_idx) = metrics_to_neuron_row(metrics, roi_idx, peak_source, trace_source);

    if ~isempty(metrics.bursts)
        new_rows = bursts_to_event_rows(metrics.bursts, roi_idx);
        burst_event_rows = [burst_event_rows; new_rows(:)]; %#ok<AGROW>
    end

    if params.make_per_neuron_figures
        plot_one_neuron_qc(metrics, trace_for_qc, trace_source, freq, output_dir, params);
    end
end

neuron_table = struct2table(neuron_rows);
if isempty(burst_event_rows)
    burst_event_table = struct2table(empty_burst_event_row());
    burst_event_table(1,:) = [];
else
    burst_event_table = struct2table(burst_event_rows);
end

% Optional threshold sensitivity
sensitivity_table = table();
if params.run_threshold_sensitivity
    sensitivity_table = run_burst_threshold_sensitivity(spike_times_s, recording_duration_s, params, frame_dt_ms);
end

% Save outputs
isi_stat_qc_summary = struct();
isi_stat_qc_summary.params = params;
isi_stat_qc_summary.results_path = save_path;
isi_stat_qc_summary.peak_source = peak_source;
isi_stat_qc_summary.trace_source = trace_source;
isi_stat_qc_summary.frame_rate_hz = freq;
isi_stat_qc_summary.frame_dt_ms = frame_dt_ms;
isi_stat_qc_summary.recording_duration_s = recording_duration_s;
isi_stat_qc_summary.neuron_table = neuron_table;
isi_stat_qc_summary.burst_event_table = burst_event_table;
isi_stat_qc_summary.sensitivity_table = sensitivity_table;
isi_stat_qc_summary.per_neuron_metrics = per_neuron_metrics;

save(fullfile(output_dir, 'isi_stat_qc_summary_revised.mat'), 'isi_stat_qc_summary');
writetable(neuron_table, fullfile(output_dir, 'neuron_summary_revised.csv'));
writetable(burst_event_table, fullfile(output_dir, 'burst_event_table_revised.csv'));
if ~isempty(sensitivity_table)
    writetable(sensitivity_table, fullfile(output_dir, 'threshold_sensitivity_revised.csv'));
end

plot_population_summary(neuron_table, output_dir, params);
plot_roi_isi_violin(per_neuron_metrics, output_dir, params);

fprintf('Revised ISI statistic QC saved to %s\n', output_dir);
fprintf('Using peaks from %s; using trace from %s.\n', peak_source, trace_source);

%% Save Explicit Results
% Comment removed after encoding repair.
save_filename = fullfile(save_path, '-1_explicit_results.mat');
% Save the explicit result bundle alongside the stage-specific result files.

% Comment removed after encoding repair.
% Persist the lightweight result objects instead of serializing the entire workspace.
results_summary = struct();
results_summary.analysis_info = analysis_info;
results_summary.result_files = struct( ...
    'movie_info', movie_info_path, ...
    'trace_results', trace_results_path, ...
    'peak_results', peak_results_path, ...
    'ap_results', ap_results_path);
results_summary.available_trace_stages = string(fieldnames(trace_results));
results_summary.available_peak_stages = string(setdiff(fieldnames(peak_results), {'current_result'}, 'stable'));
results_summary.available_ap_stages = string(setdiff(fieldnames(ap_results), {'current_result'}, 'stable'));
results_summary.current_peak_result = ternary(isfield(peak_results, 'current_result'), string(peak_results.current_result), "");
results_summary.current_ap_result = ternary(isfield(ap_results, 'current_result'), string(ap_results.current_result), "");
findmode_summary = "";
if exist('findmode', 'var') == 1
    findmode_summary = string(findmode);
end
bleachmode_summary = "";
if exist('bleachmode', 'var') == 1
    bleachmode_summary = string(bleachmode);
end
parts_summary = NaN;
if exist('parts', 'var') == 1
    parts_summary = parts;
end
min_peak_prominence_factor_summary = NaN;
if exist('MinPeakProminence_factor', 'var') == 1
    min_peak_prominence_factor_summary = MinPeakProminence_factor;
end
min_peak_distance_factor_summary = NaN;
if exist('MinPeakDistance_factor', 'var') == 1
    min_peak_distance_factor_summary = MinPeakDistance_factor;
end
min_peak_height_summary = NaN;
if exist('MinPeakHeight', 'var') == 1
    min_peak_height_summary = MinPeakHeight;
end
peak_polarity_mode_summary = "";
if exist('peak_polarity_mode', 'var') == 1
    peak_polarity_mode_summary = string(peak_polarity_mode);
end
peak_min_prominence_summary = NaN;
if exist('peak_min_prominence', 'var') == 1
    peak_min_prominence_summary = peak_min_prominence;
end
peak_min_distance_frames_summary = NaN;
if exist('peak_min_distance_frames', 'var') == 1
    peak_min_distance_frames_summary = peak_min_distance_frames;
end
peak_min_height_sensitivity_summary = NaN;
if exist('peak_min_height', 'var') == 1
    peak_min_height_sensitivity_summary = peak_min_height;
end
results_summary.parameters = struct( ...
    'findmode', findmode_summary, ...
    'bleachmode', bleachmode_summary, ...
    'parts', parts_summary, ...
    'min_peak_prominence_factor', min_peak_prominence_factor_summary, ...
    'min_peak_distance_factor', min_peak_distance_factor_summary, ...
    'min_peak_height', min_peak_height_summary, ...
    'sensitivity_peak_polarity_mode', peak_polarity_mode_summary, ...
    'sensitivity_peak_min_prominence', peak_min_prominence_summary, ...
    'sensitivity_peak_min_distance_frames', peak_min_distance_frames_summary, ...
    'sensitivity_peak_min_height', peak_min_height_sensitivity_summary);

save(movie_info_path, 'movie_info');
save(trace_results_path, 'trace_results', '-v7.3');
save(peak_results_path, 'peak_results', '-v7.3');
save(ap_results_path, 'ap_results', '-v7.3');
save(save_filename, 'results_summary', '-v7.3');



%%

%% Per-neuron analysis
function metrics = analyze_one_spontaneous_neuron(t_spikes, recording_duration_s, params, neuron_id, frame_dt_ms)
spike_count = numel(t_spikes);
ISI_s = diff(t_spikes);
ISI_ms = ISI_s * 1000;

metrics = struct();
metrics.neuron_id = neuron_id;
metrics.spike_times_s = t_spikes;
metrics.recording_duration_s = recording_duration_s;
metrics.spike_count = spike_count;
metrics.mean_firing_rate_hz = spike_count / max(recording_duration_s, eps);
metrics.ISI_s = ISI_s;
metrics.ISI_ms = ISI_ms;
metrics.median_ISI_ms = safe_median(ISI_ms);   % Q50
metrics.min_ISI_ms = safe_min(ISI_ms);
metrics.max_ISI_ms = safe_max(ISI_ms);
metrics.frame_dt_ms = frame_dt_ms;

% Basic spike-detection sanity checks.
metrics.suspicious_short_ISI_count = sum(ISI_ms < params.suspicious_short_isi_ms);
metrics.suspicious_short_ISI_fraction = safe_div(metrics.suspicious_short_ISI_count, max(numel(ISI_ms), 0));
metrics.temporal_resolution_warning = frame_dt_ms > params.burst_start_isi_threshold_ms / params.min_frames_per_burst_interval;

% Firing-rate stationarity.
[centers_s, fr_hz, fit_slope, first_FR, second_FR] = compute_sliding_fr(t_spikes, recording_duration_s, params.window_s, params.window_step_s);
metrics.window_centers_s = centers_s;
metrics.window_FR_hz = fr_hz;
metrics.FR_slope_hz_per_s = fit_slope;
metrics.FR_slope_hz_per_min = fit_slope * 60;
metrics.first_half_FR_hz = first_FR;
metrics.second_half_FR_hz = second_FR;
metrics.second_minus_first_FR_hz = second_FR - first_FR;
metrics.half_diff_fraction = abs(metrics.second_minus_first_FR_hz) / max(metrics.mean_firing_rate_hz, eps);
metrics.window_FR_CV = safe_div(nanstd_compat(fr_hz), nanmean_compat(fr_hz));
metrics.is_rate_drifting = abs(metrics.FR_slope_hz_per_min) >= params.rate_drift_slope_threshold_hz_per_min || ...
    metrics.half_diff_fraction >= params.rate_drift_half_diff_fraction;

% ISI-based burst detection with start/end thresholds and boundary quiescence.
[bursts, burst_spike_mask] = detect_bursts_from_isi_revised(t_spikes, params);
metrics.bursts = bursts;
metrics.burst_spike_mask = burst_spike_mask;
metrics.burst_count = numel(bursts);
metrics.burst_rate_per_min = metrics.burst_count / max(recording_duration_s, eps) * 60;
metrics.burst_spike_count = sum(burst_spike_mask);
metrics.burst_spike_fraction = safe_div(metrics.burst_spike_count, spike_count);

if isempty(bursts)
    metrics.median_spikes_per_burst = NaN;
    metrics.median_burst_duration_ms = NaN;
    metrics.median_intra_burst_ISI_ms = NaN;
    metrics.median_intra_burst_frequency_hz = NaN;
    metrics.burst_onsets_s = [];
    metrics.IBI_s = [];
    metrics.median_IBI_s = NaN;
    metrics.CV_IBI = NaN;
    metrics.edge_truncated_burst_count = 0;
else
    metrics.median_spikes_per_burst = safe_median([bursts.spikes_per_burst]);
    metrics.median_burst_duration_ms = safe_median([bursts.burst_duration_ms]);
    all_intra_ms = [bursts.intra_burst_ISI_ms];
    metrics.median_intra_burst_ISI_ms = safe_median(all_intra_ms);
    metrics.median_intra_burst_frequency_hz = 1000 / metrics.median_intra_burst_ISI_ms;
    metrics.burst_onsets_s = [bursts.burst_onset_s];
    metrics.IBI_s = diff(metrics.burst_onsets_s);
    metrics.median_IBI_s = safe_median(metrics.IBI_s);
    if numel(metrics.IBI_s) >= 2 && nanmean_compat(metrics.IBI_s) > 0
        metrics.CV_IBI = nanstd_compat(metrics.IBI_s) / nanmean_compat(metrics.IBI_s);
    else
        metrics.CV_IBI = NaN;
    end
    metrics.edge_truncated_burst_count = sum([bursts.edge_truncated]);
end

% Non-burst tonic regularity.
reg = compute_nonburst_regularity(t_spikes, burst_spike_mask, params);
metrics.nonburst_ISI_s = reg.nonburst_ISI_s;
metrics.nonburst_segment_count = reg.nonburst_segment_count;
metrics.nonburst_ISI_count = numel(reg.nonburst_ISI_s);
metrics.CV_ISI_nonburst = reg.CV_ISI_nonburst;
metrics.LV_nonburst = reg.LV_nonburst;
metrics.LvR_nonburst = reg.LvR_nonburst;
metrics.LvR_negative_term_fraction = reg.LvR_negative_term_fraction;

% Whole-train metrics for reference only. Do not use whole-train LvR to decide burst phenotype.
metrics.CV_ISI_all = NaN;
metrics.LV_all = NaN;
metrics.LvR_all = NaN;
if numel(ISI_s) >= 2 && nanmean_compat(ISI_s) > 0
    metrics.CV_ISI_all = nanstd_compat(ISI_s) / nanmean_compat(ISI_s);
    whole_reg = compute_local_variation_from_adjacent(ISI_s(1:end-1), ISI_s(2:end), params.LvR_refractory_R_s);
    metrics.LV_all = whole_reg.LV;
    metrics.LvR_all = whole_reg.LvR;
end

% High-frequency tonic guard: many short ISIs but no quiescence-defined bursts.
metrics.short_ISI_fraction = safe_div(sum(ISI_ms <= params.burst_start_isi_threshold_ms), numel(ISI_ms));
metrics.sustained_high_frequency_candidate = metrics.short_ISI_fraction >= params.high_short_isi_fraction_threshold && ...
    metrics.burst_spike_fraction < params.high_burst_fraction_threshold;

metrics.phenotype_label = assign_firing_phenotype_revised(metrics, params);
end

function [bursts, burst_spike_mask] = detect_bursts_from_isi_revised(t_spikes, params)
% ISI burst detector using start/end thresholds.
% A burst starts when an ISI <= burst_start_isi_threshold_ms and is bounded
% by quiescence: previous and following ISIs should be >= burst_end_isi_threshold_ms
% unless params.allow_edge_bursts is true. Once a burst starts, it continues
% until the first ISI >= burst_end_isi_threshold_ms. Candidate bursts must
% contain at least min_spikes_per_burst spikes and enough fast ISIs.

bursts = struct('burst_onset_s', {}, 'burst_offset_s', {}, 'burst_duration_ms', {}, ...
    'spike_indices_in_burst', {}, 'spikes_per_burst', {}, 'intra_burst_ISI_ms', {}, ...
    'pre_burst_ISI_ms', {}, 'post_burst_ISI_ms', {}, 'fast_ISI_count', {}, ...
    'edge_truncated', {});

spike_count = numel(t_spikes);
burst_spike_mask = false(1, spike_count);
if spike_count < params.min_spikes_per_burst
    return;
end

ISI_ms = diff(t_spikes) * 1000;
start_thr = params.burst_start_isi_threshold_ms;
end_thr = params.burst_end_isi_threshold_ms;
min_spikes = params.min_spikes_per_burst;
min_fast_isi_count = max(params.min_fast_isi_count_per_burst, min_spikes - 1);

k = 1;
burst_id = 0;
while k <= numel(ISI_ms)
    if ISI_ms(k) > start_thr
        k = k + 1;
        continue;
    end

    % Candidate starts at spike k. Require quiescence before the burst unless
    % the burst starts at the recording edge and edge bursts are allowed.
    if k == 1
        pre_isi = NaN;
        pre_ok = params.allow_edge_bursts;
        edge_start = true;
    else
        pre_isi = ISI_ms(k - 1);
        pre_ok = pre_isi >= end_thr;
        edge_start = false;
    end
    if params.require_quiescence_boundaries && ~pre_ok
        k = k + 1;
        continue;
    end

    run_start = k;
    run_end = k;
    % Continue burst until a clear long interval appears. The long interval
    % itself is not included inside the burst.
    while run_end < numel(ISI_ms) && ISI_ms(run_end + 1) < end_thr
        run_end = run_end + 1;
    end

    spike_indices = run_start:(run_end + 1);
    intra_ms = ISI_ms(run_start:run_end);

    if run_end < numel(ISI_ms)
        post_isi = ISI_ms(run_end + 1);
        post_ok = post_isi >= end_thr;
        edge_end = false;
    else
        post_isi = NaN;
        post_ok = params.allow_edge_bursts;
        edge_end = true;
    end

    fast_isi_count = sum(intra_ms <= start_thr);
    duration_ms = (t_spikes(spike_indices(end)) - t_spikes(spike_indices(1))) * 1000;

    candidate_ok = numel(spike_indices) >= min_spikes && ...
        fast_isi_count >= min_fast_isi_count && ...
        duration_ms <= params.max_burst_duration_ms;
    if params.require_quiescence_boundaries
        candidate_ok = candidate_ok && post_ok;
    end

    if candidate_ok
        burst_id = burst_id + 1;
        bursts(burst_id).burst_onset_s = t_spikes(spike_indices(1)); %#ok<AGROW>
        bursts(burst_id).burst_offset_s = t_spikes(spike_indices(end));
        bursts(burst_id).burst_duration_ms = duration_ms;
        bursts(burst_id).spike_indices_in_burst = spike_indices;
        bursts(burst_id).spikes_per_burst = numel(spike_indices);
        bursts(burst_id).intra_burst_ISI_ms = intra_ms;
        bursts(burst_id).pre_burst_ISI_ms = pre_isi;
        bursts(burst_id).post_burst_ISI_ms = post_isi;
        bursts(burst_id).fast_ISI_count = fast_isi_count;
        bursts(burst_id).edge_truncated = edge_start || edge_end;
        burst_spike_mask(spike_indices) = true;
    end

    k = run_end + 1;
end
end

function reg = compute_nonburst_regularity(t_spikes, burst_spike_mask, params)
reg = struct();
reg.nonburst_ISI_s = [];
reg.nonburst_segment_count = 0;
reg.CV_ISI_nonburst = NaN;
reg.LV_nonburst = NaN;
reg.LvR_nonburst = NaN;
reg.LvR_negative_term_fraction = NaN;

nonburst_indices = find(~burst_spike_mask);
if numel(nonburst_indices) < 3
    return;
end

adjacent_isi_1 = [];
adjacent_isi_2 = [];
segment_breaks = [0 find(diff(nonburst_indices) > 1) numel(nonburst_indices)];
for seg_idx = 1:(numel(segment_breaks) - 1)
    segment_indices = nonburst_indices(segment_breaks(seg_idx)+1:segment_breaks(seg_idx+1));
    if numel(segment_indices) < 2
        continue;
    end
    reg.nonburst_segment_count = reg.nonburst_segment_count + 1;
    segment_times = t_spikes(segment_indices);
    segment_isi = diff(segment_times);
    reg.nonburst_ISI_s = [reg.nonburst_ISI_s segment_isi]; %#ok<AGROW>
    if numel(segment_isi) >= 2
        adjacent_isi_1 = [adjacent_isi_1 segment_isi(1:end-1)]; %#ok<AGROW>
        adjacent_isi_2 = [adjacent_isi_2 segment_isi(2:end)]; %#ok<AGROW>
    end
end

if numel(reg.nonburst_ISI_s) >= 2 && nanmean_compat(reg.nonburst_ISI_s) > 0
    reg.CV_ISI_nonburst = nanstd_compat(reg.nonburst_ISI_s) / nanmean_compat(reg.nonburst_ISI_s);
end
if ~isempty(adjacent_isi_1)
    lv = compute_local_variation_from_adjacent(adjacent_isi_1, adjacent_isi_2, params.LvR_refractory_R_s);
    reg.LV_nonburst = lv.LV;
    reg.LvR_nonburst = lv.LvR;
    reg.LvR_negative_term_fraction = lv.negative_term_fraction;
end
end

function lv = compute_local_variation_from_adjacent(isi1_s, isi2_s, R_s)
isi_sum = isi1_s + isi2_s;
valid = isfinite(isi_sum) & isi_sum > 0 & isfinite(isi1_s) & isfinite(isi2_s);
isi1_s = isi1_s(valid);
isi2_s = isi2_s(valid);
isi_sum = isi_sum(valid);

lv = struct('LV', NaN, 'LvR', NaN, 'negative_term_fraction', NaN);
if isempty(isi_sum)
    return;
end
ratio_sq = ((isi2_s - isi1_s) ./ isi_sum).^2;
LV_terms = 3 * ratio_sq;
lv.LV = nanmean_compat(LV_terms);

% Revised local variation with refractory correction. Shinomoto et al.
% compensate for a refractory period R by subtracting R from adjacent ISIs;
% the first-order expansion gives a multiplicative factor (1 + 4R/(ISI_i+ISI_{i+1})).
% R should not be interpreted as a burst threshold; it is only a correction
% term for local ISI irregularity metrics.
refractory_factor = 1 + 4 * R_s ./ isi_sum;
LvR_terms = 3 * refractory_factor .* ratio_sq;
lv.LvR = nanmean_compat(LvR_terms);
% This field is kept for compatibility. With the correct first-order LvR
% factor it should be zero unless invalid inputs were already filtered out.
lv.LvR_negative_term_fraction = 0;
lv.negative_term_fraction = 0;
end

function label = assign_firing_phenotype_revised(metrics, params)
if metrics.spike_count < params.min_spike_count
    label = "insufficient_data";
    return;
end
if metrics.suspicious_short_ISI_fraction >= params.suspicious_short_isi_fraction_threshold
    label = "check_peak_detection_short_ISI";
    return;
end

is_burst_like = metrics.burst_spike_fraction >= params.high_burst_fraction_threshold;
is_mixed = metrics.burst_spike_fraction >= params.low_burst_fraction_threshold && ...
    metrics.burst_spike_fraction < params.high_burst_fraction_threshold;

if metrics.is_rate_drifting && is_burst_like
    label = "rate_drifting_burst_like";
elseif metrics.is_rate_drifting
    label = "rate_drifting_firing";
elseif metrics.sustained_high_frequency_candidate
    label = "sustained_high_frequency_tonic_candidate";
elseif is_burst_like
    if metrics.burst_count >= params.min_burst_count_for_periodicity && isfinite(metrics.CV_IBI)
        if metrics.CV_IBI <= params.low_CV_IBI_threshold
            label = "periodic_bursting";
        elseif metrics.CV_IBI >= params.high_CV_IBI_threshold
            label = "irregular_bursting";
        else
            label = "bursting_intermediate_IBI_variability";
        end
    else
        label = "burst_like_unclassified_periodicity";
    end
elseif is_mixed
    label = "mixed_tonic_burst";
elseif isfinite(metrics.LvR_nonburst) && metrics.LvR_nonburst <= params.low_LvR_threshold
    label = "locally_regular_tonic";
elseif isfinite(metrics.LvR_nonburst) && metrics.LvR_nonburst >= params.high_LvR_threshold
    label = "locally_irregular_tonic";
else
    label = "tonic_unclassified_regularity";
end
end

%% FR and sensitivity helpers
function [centers_s, fr_hz, fit_slope, first_FR, second_FR] = compute_sliding_fr(t_spikes, duration_s, window_s, step_s)
if duration_s <= 0
    centers_s = [];
    fr_hz = [];
    fit_slope = NaN;
    first_FR = NaN;
    second_FR = NaN;
    return;
end
window_s = min(window_s, duration_s);
step_s = min(step_s, window_s);
starts = 0:step_s:max(duration_s - window_s, 0);
if isempty(starts)
    starts = 0;
end
centers_s = starts + window_s / 2;
fr_hz = NaN(size(centers_s));
for i = 1:numel(starts)
    t0 = starts(i);
    t1 = min(t0 + window_s, duration_s);
    in_window = t_spikes >= t0 & t_spikes < t1;
    fr_hz(i) = sum(in_window) / max(t1 - t0, eps);
end

if numel(centers_s) >= 2 && sum(isfinite(fr_hz)) >= 2
    idx = isfinite(fr_hz);
    p = polyfit(centers_s(idx), fr_hz(idx), 1);
    fit_slope = p(1);
else
    fit_slope = NaN;
end
first_FR = sum(t_spikes <= duration_s / 2) / max(duration_s / 2, eps);
second_FR = sum(t_spikes > duration_s / 2) / max(duration_s / 2, eps);
end

function sensitivity_table = run_burst_threshold_sensitivity(spike_times_s, recording_duration_s, params, frame_dt_ms)
rows = [];
row_idx = 0;
thresholds = params.sensitivity_burst_start_thresholds_ms(:)';
for neuron_id = 1:numel(spike_times_s)
    for thr = thresholds
        params_i = params;
        params_i.burst_start_isi_threshold_ms = thr;
        if params.sensitivity_scale_end_threshold
            params_i.burst_end_isi_threshold_ms = max(params.min_burst_end_isi_threshold_ms, params.burst_end_to_start_ratio * thr);
        end
        metrics_i = analyze_one_spontaneous_neuron(spike_times_s{neuron_id}, recording_duration_s, params_i, neuron_id, frame_dt_ms);
        row_idx = row_idx + 1;
        rows(row_idx).neuron_id = neuron_id; %#ok<AGROW>
        rows(row_idx).burst_start_isi_threshold_ms = thr;
        rows(row_idx).burst_end_isi_threshold_ms = params_i.burst_end_isi_threshold_ms;
        rows(row_idx).burst_count = metrics_i.burst_count;
        rows(row_idx).burst_spike_fraction = metrics_i.burst_spike_fraction;
        rows(row_idx).CV_IBI = metrics_i.CV_IBI;
        rows(row_idx).phenotype_label = string(metrics_i.phenotype_label);
    end
end
if isempty(rows)
    sensitivity_table = table();
else
    sensitivity_table = struct2table(rows);
end
end

%% Row conversion
function row = empty_neuron_row()
row = struct();
row.neuron_id = NaN;
row.recording_duration_s = NaN;
row.spike_count = NaN;
row.mean_firing_rate_hz = NaN;
row.median_ISI_ms = NaN;
row.min_ISI_ms = NaN;
row.max_ISI_ms = NaN;
row.frame_dt_ms = NaN;
row.temporal_resolution_warning = false;
row.suspicious_short_ISI_count = NaN;
row.suspicious_short_ISI_fraction = NaN;
row.FR_slope_hz_per_min = NaN;
row.first_half_FR_hz = NaN;
row.second_half_FR_hz = NaN;
row.second_minus_first_FR_hz = NaN;
row.half_diff_fraction = NaN;
row.window_FR_CV = NaN;
row.is_rate_drifting = false;
row.short_ISI_fraction = NaN;
row.burst_count = NaN;
row.burst_rate_per_min = NaN;
row.burst_spike_count = NaN;
row.burst_spike_fraction = NaN;
row.edge_truncated_burst_count = NaN;
row.median_spikes_per_burst = NaN;
row.median_burst_duration_ms = NaN;
row.median_intra_burst_ISI_ms = NaN;
row.median_intra_burst_frequency_hz = NaN;
row.median_IBI_s = NaN;
row.CV_IBI = NaN;
row.nonburst_segment_count = NaN;
row.nonburst_ISI_count = NaN;
row.CV_ISI_nonburst = NaN;
row.LV_nonburst = NaN;
row.LvR_nonburst = NaN;
row.LvR_negative_term_fraction = NaN;
row.CV_ISI_all = NaN;
row.LV_all = NaN;
row.LvR_all = NaN;
row.sustained_high_frequency_candidate = false;
row.phenotype_label = "";
row.peak_source = "";
row.trace_source = "";
end

function row = metrics_to_neuron_row(m, neuron_id, peak_source, trace_source)
row = empty_neuron_row();
row.neuron_id = neuron_id;
row.recording_duration_s = m.recording_duration_s;
row.spike_count = m.spike_count;
row.mean_firing_rate_hz = m.mean_firing_rate_hz;
row.median_ISI_ms = m.median_ISI_ms;
row.min_ISI_ms = m.min_ISI_ms;
row.max_ISI_ms = m.max_ISI_ms;
row.frame_dt_ms = m.frame_dt_ms;
row.temporal_resolution_warning = m.temporal_resolution_warning;
row.suspicious_short_ISI_count = m.suspicious_short_ISI_count;
row.suspicious_short_ISI_fraction = m.suspicious_short_ISI_fraction;
row.FR_slope_hz_per_min = m.FR_slope_hz_per_min;
row.first_half_FR_hz = m.first_half_FR_hz;
row.second_half_FR_hz = m.second_half_FR_hz;
row.second_minus_first_FR_hz = m.second_minus_first_FR_hz;
row.half_diff_fraction = m.half_diff_fraction;
row.window_FR_CV = m.window_FR_CV;
row.is_rate_drifting = m.is_rate_drifting;
row.short_ISI_fraction = m.short_ISI_fraction;
row.burst_count = m.burst_count;
row.burst_rate_per_min = m.burst_rate_per_min;
row.burst_spike_count = m.burst_spike_count;
row.burst_spike_fraction = m.burst_spike_fraction;
row.edge_truncated_burst_count = m.edge_truncated_burst_count;
row.median_spikes_per_burst = m.median_spikes_per_burst;
row.median_burst_duration_ms = m.median_burst_duration_ms;
row.median_intra_burst_ISI_ms = m.median_intra_burst_ISI_ms;
row.median_intra_burst_frequency_hz = m.median_intra_burst_frequency_hz;
row.median_IBI_s = m.median_IBI_s;
row.CV_IBI = m.CV_IBI;
row.nonburst_segment_count = m.nonburst_segment_count;
row.nonburst_ISI_count = m.nonburst_ISI_count;
row.CV_ISI_nonburst = m.CV_ISI_nonburst;
row.LV_nonburst = m.LV_nonburst;
row.LvR_nonburst = m.LvR_nonburst;
row.LvR_negative_term_fraction = m.LvR_negative_term_fraction;
row.CV_ISI_all = m.CV_ISI_all;
row.LV_all = m.LV_all;
row.LvR_all = m.LvR_all;
row.sustained_high_frequency_candidate = m.sustained_high_frequency_candidate;
row.phenotype_label = string(m.phenotype_label);
row.peak_source = string(peak_source);
row.trace_source = string(trace_source);
end

function row = empty_burst_event_row()
row = struct();
row.neuron_id = NaN;
row.burst_id = NaN;
row.burst_onset_s = NaN;
row.burst_offset_s = NaN;
row.burst_duration_ms = NaN;
row.spikes_per_burst = NaN;
row.median_intra_burst_ISI_ms = NaN;
row.pre_burst_ISI_ms = NaN;
row.post_burst_ISI_ms = NaN;
row.fast_ISI_count = NaN;
row.edge_truncated = false;
end

function rows = bursts_to_event_rows(bursts, neuron_id)
rows = repmat(empty_burst_event_row(), numel(bursts), 1);
for i = 1:numel(bursts)
    rows(i).neuron_id = neuron_id;
    rows(i).burst_id = i;
    rows(i).burst_onset_s = bursts(i).burst_onset_s;
    rows(i).burst_offset_s = bursts(i).burst_offset_s;
    rows(i).burst_duration_ms = bursts(i).burst_duration_ms;
    rows(i).spikes_per_burst = bursts(i).spikes_per_burst;
    rows(i).median_intra_burst_ISI_ms = safe_median(bursts(i).intra_burst_ISI_ms);
    rows(i).pre_burst_ISI_ms = bursts(i).pre_burst_ISI_ms;
    rows(i).post_burst_ISI_ms = bursts(i).post_burst_ISI_ms;
    rows(i).fast_ISI_count = bursts(i).fast_ISI_count;
    rows(i).edge_truncated = bursts(i).edge_truncated;
end
end

%% Plotting
function plot_one_neuron_qc(m, trace_for_qc, trace_source, freq, output_dir, params)
if params.figures_visible
    visible_mode = 'on';
else
    visible_mode = 'off';
end
fig = figure('Visible', visible_mode, 'Position', [100, 100, 1150, 1250], ...
    'Color', 'w', 'Name', sprintf('Neuron %03d revised ISI QC', m.neuron_id));

subplot(6,1,1); hold on;
if ~isempty(trace_for_qc) && size(trace_for_qc, 2) >= m.neuron_id
    trace_i = trace_for_qc(:, m.neuron_id);
    if params.frame_index_origin_is_zero_time
        t_trace = (0:(numel(trace_i)-1)) ./ freq;
    else
        t_trace = (1:numel(trace_i)) ./ freq;
    end
    plot(t_trace, trace_i, 'k');
    spike_frames = round(m.spike_times_s * freq) + double(params.frame_index_origin_is_zero_time);
    spike_frames = max(1, min(numel(trace_i), spike_frames));
    if ~isempty(spike_frames)
        plot(m.spike_times_s, trace_i(spike_frames), 'rv', 'MarkerFaceColor', 'r', 'MarkerSize', 4);
    end
    y_lim = ylim;
    for b = 1:numel(m.bursts)
        patch([m.bursts(b).burst_onset_s, m.bursts(b).burst_offset_s, ...
            m.bursts(b).burst_offset_s, m.bursts(b).burst_onset_s], ...
            [y_lim(1), y_lim(1), y_lim(2), y_lim(2)], [1.0, 0.80, 0.45], ...
            'FaceAlpha', 0.35, 'EdgeColor', 'none');
    end
    uistack(findobj(gca, 'Type', 'line'), 'top');
end
xlim([0, m.recording_duration_s]);
title(sprintf('Neuron %03d: %s', m.neuron_id, char(m.phenotype_label)), 'Interpreter', 'none');
ylabel(trace_source, 'Interpreter', 'none');
box off;

subplot(6,1,2); hold on;
plot(m.window_centers_s, m.window_FR_hz, '-o', 'Color', [0.10, 0.35, 0.75], ...
    'MarkerFaceColor', [0.10, 0.35, 0.75]);
if numel(m.window_centers_s) >= 2 && isfinite(m.FR_slope_hz_per_s)
    intercept = nanmean_compat(m.window_FR_hz) - m.FR_slope_hz_per_s * nanmean_compat(m.window_centers_s);
    plot(m.window_centers_s, m.FR_slope_hz_per_s * m.window_centers_s + intercept, '--', ...
        'Color', [0.85, 0.20, 0.20], 'LineWidth', 1.2);
end
xlim([0, m.recording_duration_s]);
ylabel('FR (Hz)');
title(sprintf('Sliding FR: slope %.3g Hz/min, half diff %.3g', m.FR_slope_hz_per_min, m.half_diff_fraction));
box off;

subplot(6,1,3); hold on;
valid_isi_ms = m.ISI_ms(isfinite(m.ISI_ms) & m.ISI_ms > 0);
if ~isempty(valid_isi_ms)
    histogram(log10(valid_isi_ms), params.isi_histogram_num_bins, ...
        'FaceColor', [0.30, 0.30, 0.30], 'EdgeColor', 'none');
end
xline(log10(params.burst_start_isi_threshold_ms), 'r-', 'LineWidth', 1.5);
xline(log10(params.burst_end_isi_threshold_ms), 'b--', 'LineWidth', 1.2);
if ~isempty(params.isi_histogram_xlim_log10_ms)
    xlim(params.isi_histogram_xlim_log10_ms);
end
ylabel('Count');
title(sprintf('log10(ISI ms): start <= %.3g ms, end/quiescence >= %.3g ms', ...
    params.burst_start_isi_threshold_ms, params.burst_end_isi_threshold_ms));
box off;

subplot(6,1,4); hold on;
if numel(m.ISI_ms) >= 2
    plot(log10(m.ISI_ms(1:end-1)), log10(m.ISI_ms(2:end)), 'k.', 'MarkerSize', 8);
    xline(log10(params.burst_start_isi_threshold_ms), 'r-');
    yline(log10(params.burst_start_isi_threshold_ms), 'r-');
end
xlabel('log10(ISI_i ms)');
ylabel('log10(ISI_{i+1} ms)');
title('ISI return map');
box off;

subplot(6,1,5); hold on;
if ~isempty(m.IBI_s)
    histogram(m.IBI_s, params.ibi_histogram_num_bins, 'FaceColor', [0.35, 0.55, 0.80], 'EdgeColor', 'none');
end
if ~isempty(params.ibi_histogram_xlim_s)
    xlim(params.ibi_histogram_xlim_s);
end
xlabel('IBI (s)');
ylabel('Count');
title(sprintf('Burst IBI: burst count %d, CV_IBI %.3g', m.burst_count, m.CV_IBI));
box off;

subplot(6,1,6); hold on;
if ~isempty(m.spike_times_s)
    plot(m.spike_times_s(~m.burst_spike_mask), zeros(1, sum(~m.burst_spike_mask)), 'k|', 'MarkerSize', 10);
    plot(m.spike_times_s(m.burst_spike_mask), zeros(1, sum(m.burst_spike_mask)), 'r|', 'MarkerSize', 12);
end
xlim([0, m.recording_duration_s]);
ylim([-1, 1]);
yticks([]);
xlabel('Time (s)');
title(sprintf('Burst fraction %.3g, non-burst LvR %.3g, non-burst LV %.3g', ...
    m.burst_spike_fraction, m.LvR_nonburst, m.LV_nonburst));
box off;

saveas(fig, fullfile(output_dir, sprintf('neuron_%03d_isi_stat_qc_revised.fig', m.neuron_id)), 'fig');
saveas(fig, fullfile(output_dir, sprintf('neuron_%03d_isi_stat_qc_revised.png', m.neuron_id)), 'png');
if ~params.figures_visible
    close(fig);
end
end

function plot_population_summary(neuron_table, output_dir, params)
if isempty(neuron_table)
    return;
end
if params.figures_visible
    visible_mode = 'on';
else
    visible_mode = 'off';
end
fig = figure('Visible', visible_mode, 'Position', [100, 100, 950, 900], ...
    'Color', 'w', 'Name', 'Population revised ISI summary');

subplot(2,1,1); hold on;
scatter(neuron_table.burst_spike_fraction, neuron_table.LvR_nonburst, 50, 'filled', ...
    'MarkerFaceColor', [0.20, 0.45, 0.70], 'MarkerFaceAlpha', 0.75);
xline(params.low_burst_fraction_threshold, 'Color', [0.7 0.7 0.7], 'LineStyle', '--');
xline(params.high_burst_fraction_threshold, 'r--');
yline(params.low_LvR_threshold, 'b--');
yline(params.high_LvR_threshold, 'b--');
xlabel('Burst spike fraction');
ylabel('Non-burst LvR');
title('Burst burden vs non-burst local regularity');
box off;

subplot(2,1,2); hold on;
burst_like = neuron_table.burst_spike_fraction >= params.high_burst_fraction_threshold;
scatter(neuron_table.burst_spike_fraction(burst_like), neuron_table.CV_IBI(burst_like), 55, 'filled', ...
    'MarkerFaceColor', [0.80, 0.25, 0.20], 'MarkerFaceAlpha', 0.75);
yline(params.low_CV_IBI_threshold, 'k--');
yline(params.high_CV_IBI_threshold, 'k--');
xlabel('Burst spike fraction');
ylabel('CV IBI');
title('Burst-like neurons: burst burden vs IBI variability');
box off;

saveas(fig, fullfile(output_dir, 'population_isi_stat_summary_revised.fig'), 'fig');
saveas(fig, fullfile(output_dir, 'population_isi_stat_summary_revised.png'), 'png');
if ~params.figures_visible
    close(fig);
end
% --- Max ISI vs Q50 ISI figure ---
fig2 = figure('Visible', visible_mode, 'Position', [150, 150, 700, 600], ...
    'Color', 'w', 'Name', 'Max ISI vs Q50 ISI');

valid = isfinite(neuron_table.max_ISI_ms) & isfinite(neuron_table.median_ISI_ms) & ...
        neuron_table.max_ISI_ms > 0 & neuron_table.median_ISI_ms > 0;

scatter(neuron_table.max_ISI_ms(valid), neuron_table.median_ISI_ms(valid), ...
    60, 'filled', 'MarkerFaceColor', [0.20, 0.45, 0.70], 'MarkerFaceAlpha', 0.8);
hold on;

% Add neuron ID labels to each point.
for i = find(valid)'
    text(neuron_table.max_ISI_ms(i), neuron_table.median_ISI_ms(i), ...
        sprintf('  %d', neuron_table.neuron_id(i)), ...
        'FontSize', 9, 'Color', [0.1 0.1 0.1]);
end

set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('Max ISI (ms)');
ylabel('Q50 ISI / Median ISI (ms)');
title('Max ISI vs Q50 ISI');
grid on;
box off;

saveas(fig2, fullfile(output_dir, 'population_maxISI_vs_Q50ISI.fig'), 'fig');
saveas(fig2, fullfile(output_dir, 'population_maxISI_vs_Q50ISI.png'), 'png');

if ~params.figures_visible
    close(fig2);
end
end
function plot_roi_isi_violin(per_neuron_metrics, output_dir, params)
% Plot ROI-wise ISI distributions as violin plots.
% y-axis uses raw ISI_ms.

if isempty(per_neuron_metrics)
    return;
end

if params.figures_visible
    visible_mode = 'on';
else
    visible_mode = 'off';
end

nrois = numel(per_neuron_metrics);

fig = figure('Visible', visible_mode, 'Position', [120, 120, 1100, 650], ...
    'Color', 'w', 'Name', 'ROI-wise ISI violin plot');
hold on;

max_width = 0.38;
min_points_for_violin = 3;

for roi_idx = 1:nrois
    m = per_neuron_metrics{roi_idx};
    if isempty(m) || ~isfield(m, 'ISI_ms')
        continue;
    end

    isi_ms = m.ISI_ms(:);
    isi_ms = isi_ms(isfinite(isi_ms) & isi_ms > 0);

    if isempty(isi_ms)
        continue;
    end

    y = isi_ms;

    % Draw violin if enough ISIs exist.
    if numel(y) >= min_points_for_violin
        if exist('ksdensity', 'file') == 2
            y_min = min(y);
            y_max = max(y);

            if y_min == y_max
                y_grid = y_min + linspace(-0.05, 0.05, 50);
                dens = ones(size(y_grid));
            else
                y_grid = linspace(y_min, y_max, 100);
                dens = ksdensity(y, y_grid);
            end
        else
            % Fallback if Statistics Toolbox is unavailable.
            edges = linspace(min(y), max(y), 30);
            if numel(unique(edges)) < 2
                y_grid = y(1) + linspace(-0.05, 0.05, 50);
                dens = ones(size(y_grid));
            else
                counts = histcounts(y, edges, 'Normalization', 'pdf');
                y_grid = edges(1:end-1) + diff(edges)/2;
                dens = counts;
            end
        end

        if max(dens) > 0
            dens = dens ./ max(dens) * max_width;
        else
            dens = zeros(size(dens));
        end

        x_left = roi_idx - dens;
        x_right = roi_idx + dens;

        patch([x_left, fliplr(x_right)], [y_grid, fliplr(y_grid)], ...
            [0.55, 0.70, 0.85], ...
            'FaceAlpha', 0.55, ...
            'EdgeColor', [0.25, 0.35, 0.45], ...
            'LineWidth', 0.8);
    end

    % Overlay individual ISIs with small horizontal jitter.
    jitter_width = 0.08;
    x_jitter = roi_idx + (rand(size(y)) - 0.5) * 2 * jitter_width;
    plot(x_jitter, y, '.', ...
        'Color', [0.15, 0.15, 0.15], ...
        'MarkerSize', 5);

    % Overlay Q25, median, Q75.
    q25 = prctile(y, 25);
    q50 = prctile(y, 50);
    q75 = prctile(y, 75);

    plot([roi_idx - 0.25, roi_idx + 0.25], [q50, q50], ...
        'k-', 'LineWidth', 2.0);

    plot([roi_idx, roi_idx], [q25, q75], ...
        'k-', 'LineWidth', 1.2);

    plot([roi_idx - 0.15, roi_idx + 0.15], [q25, q25], ...
        'k-', 'LineWidth', 1.0);
    plot([roi_idx - 0.15, roi_idx + 0.15], [q75, q75], ...
        'k-', 'LineWidth', 1.0);
end

xlim([0.5, nrois + 0.5]);
xticks(1:nrois);
xlabel('ROI / neuron ID');

ylabel('ISI (ms)');
title('ROI-wise ISI distribution');

grid on;
box off;

saveas(fig, fullfile(output_dir, 'roi_isi_violin_ISI_ms.fig'), 'fig');
saveas(fig, fullfile(output_dir, 'roi_isi_violin_ISI_ms.png'), 'png');

if ~params.figures_visible
    close(fig);
end
end
%% Parameters
function params = default_isi_stat_params_revised()
params = struct();

% Timestamp convention.
params.frame_index_origin_is_zero_time = true;  % frame 1 -> t = 0 s.

% Firing-rate stationarity.
params.window_s = 10;
params.window_step_s = 5;
params.rate_drift_slope_threshold_hz_per_min = 0.05;
params.rate_drift_half_diff_fraction = 0.5;

% Data sufficiency and spike-detection QC.
params.min_spike_count = 10;
params.suspicious_short_isi_ms = 1.0;          % suspicious duplicate/false split interval.
params.suspicious_short_isi_fraction_threshold = 0.05;
params.min_frames_per_burst_interval = 3;      % warn if burst threshold spans fewer frames.

% Burst detection. These are deliberately exposed for sensitivity analysis.
params.burst_start_isi_threshold_ms = 20;      % high-frequency interval required to start/qualify a burst.
params.burst_end_isi_threshold_ms = NaN;       % if NaN, set to max(50 ms, 3*start threshold).
params.min_burst_end_isi_threshold_ms = 50;
params.burst_end_to_start_ratio = 3;
params.min_spikes_per_burst = 3;
params.min_fast_isi_count_per_burst = NaN;     % if NaN, set to min_spikes_per_burst - 1.
params.max_burst_duration_ms = Inf;
params.require_quiescence_boundaries = true;   % key guard against sustained high-frequency tonic misclassification.
params.allow_edge_bursts = false;              % edge-truncated candidates are excluded by default.
params.high_short_isi_fraction_threshold = 0.50;

% Burst interpretation.
params.min_burst_count_for_periodicity = 5;
params.low_burst_fraction_threshold = 0.10;
params.high_burst_fraction_threshold = 0.30;
params.low_CV_IBI_threshold = 0.50;
params.high_CV_IBI_threshold = 1.00;

% Non-burst local regularity.
params.low_LvR_threshold = 0.80;
params.high_LvR_threshold = 1.20;
params.LvR_refractory_R_s = 0.005;

% Plots and output.
params.make_per_neuron_figures = true;
params.figures_visible = false;
params.isi_histogram_num_bins = 100;
params.ibi_histogram_num_bins = 20;
params.isi_histogram_xlim_log10_ms = [];
params.ibi_histogram_xlim_s = [];

% Sensitivity analysis.
params.run_threshold_sensitivity = true;
params.sensitivity_burst_start_thresholds_ms = [10 20 30];
params.sensitivity_scale_end_threshold = true;
end

function params = finalize_isi_stat_params(params)
if ~isfield(params, 'burst_start_isi_threshold_ms') || isempty(params.burst_start_isi_threshold_ms)
    params.burst_start_isi_threshold_ms = 20;
end
if ~isfield(params, 'burst_end_isi_threshold_ms') || isempty(params.burst_end_isi_threshold_ms) || isnan(params.burst_end_isi_threshold_ms)
    params.burst_end_isi_threshold_ms = max(params.min_burst_end_isi_threshold_ms, ...
        params.burst_end_to_start_ratio * params.burst_start_isi_threshold_ms);
end
if ~isfield(params, 'min_fast_isi_count_per_burst') || isempty(params.min_fast_isi_count_per_burst) || isnan(params.min_fast_isi_count_per_burst)
    params.min_fast_isi_count_per_burst = max(1, params.min_spikes_per_burst - 1);
end
params.LvR_refractory_R_s = double(params.LvR_refractory_R_s);
end


%% Small compatibility utilities
function base = merge_struct_fields(base, override)
if ~isstruct(override)
    return;
end
names = fieldnames(override);
for i = 1:numel(names)
    base.(names{i}) = override.(names{i});
end
end

function y = safe_div(a, b)
if isempty(b) || b == 0 || ~isfinite(b)
    y = NaN;
else
    y = a ./ b;
end
end

function y = safe_median(x)
x = x(isfinite(x));
if isempty(x)
    y = NaN;
else
    y = median(x);
end
end

function y = safe_min(x)
x = x(isfinite(x));
if isempty(x)
    y = NaN;
else
    y = min(x);
end
end

function y = nanmean_compat(x)
x = x(isfinite(x));
if isempty(x)
    y = NaN;
else
    y = mean(x);
end
end

function y = nanstd_compat(x)
x = x(isfinite(x));
if numel(x) <= 1
    y = NaN;
else
    y = std(x);
end
end
function y = safe_max(x)
x = x(isfinite(x));
if isempty(x)
    y = NaN;
else
    y = max(x);
end
end


function [peak_table, peaks_index, peaks_amplitude, peaks_polarity] = detect_sensitivity_peaks(traces_sensitivity, freq, params)
nrois = size(traces_sensitivity, 2);
peaks_polarity = cell(1, nrois);
peak_table = empty_peak_table();
next_peak_id = 1;
min_distance = max(1, round(params.min_peak_distance_frames));
created_at = datetime("now");

for roi_idx = 1:nrois
    trace_i = traces_sensitivity(:, roi_idx);
    polarity = choose_peak_polarity(trace_i, params.polarity_mode);
    peaks_polarity{roi_idx} = polarity;
    plot_trace = trace_i * polarity;

    min_prominence = resolve_peak_min_prominence(plot_trace, params);
    [~, peak_x] = findpeaks(plot_trace, ...
        'MinPeakProminence', min_prominence, ...
        'MinPeakDistance', min_distance, ...
        'MinPeakHeight', params.min_peak_height);

    if isempty(peak_x)
        continue;
    end

    peak_count = numel(peak_x);
    new_rows = table( ...
        (next_peak_id:next_peak_id+peak_count-1)', ...
        repmat(roi_idx, peak_count, 1), ...
        peak_x(:), ...
        peak_x(:) ./ freq, ...
        repmat(polarity, peak_count, 1), ...
        trace_i(peak_x(:)), ...
        repmat("accepted", peak_count, 1), ...
        repmat("auto", peak_count, 1), ...
        NaN(peak_count, 1), ...
        repmat("sensitivity_detected", peak_count, 1), ...
        repmat(created_at, peak_count, 1), ...
        'VariableNames', peak_table.Properties.VariableNames);
    peak_table = [peak_table; new_rows];
    next_peak_id = next_peak_id + peak_count;
end

[peaks_index, peaks_amplitude, peaks_polarity] = peak_table_to_peak_cells(peak_table, nrois, peaks_polarity);
end

function min_prominence = resolve_peak_min_prominence(plot_trace, params)
prominence_value = params.min_peak_prominence;
prominence_mode = "relative_factor";
if isfield(params, 'min_peak_prominence_mode') && ~isempty(params.min_peak_prominence_mode)
    prominence_mode = lower(string(params.min_peak_prominence_mode));
end

switch prominence_mode
    case "absolute"
        min_prominence = prominence_value;
    otherwise
        trace_mean = mean(plot_trace, 'omitnan');
        trace_max = max(plot_trace, [], 'omitnan');
        trace_scale = trace_max - trace_mean;
        if ~isfinite(trace_scale) || trace_scale < 0
            trace_scale = 0;
        end
        min_prominence = prominence_value * trace_scale;
end

if ~isfinite(min_prominence) || min_prominence < 0
    min_prominence = 0;
end
end

function [peak_table, edit_history] = edit_sensitivity_peaks(traces_sensitivity, peak_table, peaks_polarity, params)
nrois = size(traces_sensitivity, 2);
edit_history = struct('action', {}, 'roi', {}, 'peak_id', {}, 'index_before', {}, ...
    'index_after', {}, 'mode', {}, 'created_at', {});
max_peak_id = 0;
if height(peak_table) > 0
    max_peak_id = max(peak_table.peak_id);
end
next_peak_id = max_peak_id + 1;
quit_editor = false;

for roi_idx = 1:nrois
    if quit_editor
        break;
    end

    mode = string(params.default_mode);
    reset_table = peak_table;
    fig = figure('Name', sprintf('ROI %d sensitivity peak editor', roi_idx));
    set(fig, 'Position', get(0, 'Screensize'));

    while ishandle(fig)
        ax = draw_peak_editor(fig, traces_sensitivity(:, roi_idx), peak_table, ...
            peaks_polarity{roi_idx}, roi_idx, mode);
        key = wait_for_peak_editor_key(fig);
        if ~ishandle(fig)
            break;
        end

        switch key
            case "a"
                mode = "add_box";
                [box_position, box_ok] = get_peak_editor_box(ax);
                if ~box_ok || ~ishandle(fig)
                    continue;
                end
                [peak_table, edit_history, next_peak_id, added_count] = add_peaks_in_box( ...
                    traces_sensitivity(:, roi_idx), peak_table, edit_history, next_peak_id, ...
                    roi_idx, peaks_polarity{roi_idx}, params, box_position, mode);
                fprintf('ROI %d: added %d peaks from rectangle.\n', roi_idx, added_count);
            case "d"
                mode = "delete_box";
                [box_position, box_ok] = get_peak_editor_box(ax);
                if ~box_ok || ~ishandle(fig)
                    continue;
                end
                [peak_table, edit_history, deleted_count] = delete_peaks_in_box( ...
                    traces_sensitivity(:, roi_idx), peak_table, edit_history, ...
                    roi_idx, peaks_polarity{roi_idx}, box_position, mode);
                fprintf('ROI %d: deleted %d peaks from rectangle.\n', roi_idx, deleted_count);
            case "r"
                peak_table = reset_roi_peak_table(peak_table, reset_table, roi_idx);
                edit_history(end+1) = make_edit_history("reset", roi_idx, ...
                    NaN, NaN, NaN, mode); %#ok<SAGROW>
            case "q"
                quit_editor = true;
                close(fig);
                break;
            case {"n", "space", "return"}
                close(fig);
                break;
        end
    end
end
end

function ax = draw_peak_editor(fig, trace_i, peak_table, polarity, roi_idx, mode)
figure(fig);
clf(fig);
set(fig, 'KeyPressFcn', @(src, event) set_peak_editor_key(event, fig));
ax = axes('Parent', fig, 'Tag', 'sensitivityPeakEditorAxes', ...
    'Units', 'normalized', 'Position', [0.06 0.08 0.91 0.86]);
plot(ax, trace_i * polarity, 'k'); hold(ax, 'on');
accepted = peak_table.roi == roi_idx & peak_table.status == "accepted";
deleted = peak_table.roi == roi_idx & peak_table.status == "deleted";
accepted_idx = peak_table.index(accepted);
deleted_idx = peak_table.index(deleted);
if ~isempty(accepted_idx)
    plot(ax, accepted_idx, trace_i(accepted_idx) * polarity, 'rv', 'MarkerFaceColor', 'r');
end
if ~isempty(deleted_idx)
    plot(ax, deleted_idx, trace_i(deleted_idx) * polarity, 'x', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
end
title(ax, sprintf('ROI %d | Mode: %s | A add, D delete, N/Space next, R reset, Q quit', roi_idx, upper(mode)));
xlabel(ax, 'Frame');
ylabel(ax, 'Sensitivity x polarity');
grid(ax, 'on');
hold(ax, 'off');
end

function set_peak_editor_key(event, fig)
if any(strcmp(event.Key, {'a', 'd', 'n', 'space', 'return', 'r', 'q'}))
    fig.UserData.key = event.Key;
end
end

function key = wait_for_peak_editor_key(fig)
fig.UserData.key = [];
waitfor(fig, 'UserData');
if ishandle(fig) && isstruct(fig.UserData) && isfield(fig.UserData, 'key')
    key = string(fig.UserData.key);
else
    key = "";
end
end

function [box_position, box_ok] = get_peak_editor_box(ax)
box_position = [NaN NaN NaN NaN];
box_ok = false;
if isempty(ax)
    return;
end
try
    set(get(ax, 'Parent'), 'CurrentAxes', ax);
    rect = drawrectangle(ax);
    wait(rect);
    if isvalid(rect)
        box_position = rect.Position;
        delete(rect);
        box_ok = all(isfinite(box_position)) && box_position(3) > 0 && box_position(4) > 0;
    end
catch ME
    warning('AP_analysis3:PeakEditRectangleFailed', ...
        'Rectangle selection failed: %s', ME.message);
    return;
end
end

function [peak_table, edit_history, next_peak_id, added_count] = add_peaks_in_box( ...
    trace_i, peak_table, edit_history, next_peak_id, roi_idx, polarity, params, box_position, mode)
added_count = 0;
candidate_indices = find_box_peak_indices(trace_i, polarity, box_position, params);
if isempty(candidate_indices)
    return;
end

accepted = peak_table.roi == roi_idx & peak_table.status == "accepted";
existing_indices = peak_table.index(accepted);
min_distance = max(1, round(params.min_peak_distance_frames));

for idx = reshape(candidate_indices, 1, [])
    if ~isempty(existing_indices) && any(abs(existing_indices - idx) < min_distance)
        continue;
    end
    peak_table = add_manual_peak_row(peak_table, next_peak_id, roi_idx, idx, ...
        trace_i(idx), polarity, params.frame_rate_hz);
    edit_history(end+1) = make_edit_history("add_box", roi_idx, ... %#ok<SAGROW>
        next_peak_id, NaN, idx, mode);
    existing_indices = [existing_indices; idx]; %#ok<AGROW>
    next_peak_id = next_peak_id + 1;
    added_count = added_count + 1;
end
end

function [peak_table, edit_history, deleted_count] = delete_peaks_in_box( ...
    trace_i, peak_table, edit_history, roi_idx, polarity, box_position, mode)
deleted_count = 0;
accepted_rows = find(peak_table.roi == roi_idx & peak_table.status == "accepted");
if isempty(accepted_rows)
    return;
end

peak_x = peak_table.index(accepted_rows);
peak_y = trace_i(peak_x) * polarity;
inside_box = points_in_box(peak_x, peak_y, box_position);
rows_to_delete = accepted_rows(inside_box);
for row_idx = reshape(rows_to_delete, 1, [])
    peak_table.status(row_idx) = "deleted";
    edit_history(end+1) = make_edit_history("delete_box", roi_idx, ... %#ok<SAGROW>
        peak_table.peak_id(row_idx), peak_table.index(row_idx), NaN, mode);
    deleted_count = deleted_count + 1;
end
end

function peak_indices = find_box_peak_indices(trace_i, polarity, box_position, params)
plot_trace = trace_i * polarity;
x_min = box_position(1);
x_max = box_position(1) + box_position(3);
y_min = box_position(2);
y_max = box_position(2) + box_position(4);
frame_range = max(1, ceil(x_min)):min(numel(trace_i), floor(x_max));
if isempty(frame_range)
    peak_indices = [];
    return;
end

segment_trace = plot_trace(frame_range);
min_distance = max(1, round(params.min_peak_distance_frames));
[peak_y, peak_x_rel] = findpeaks(segment_trace, 'MinPeakDistance', min_distance);
peak_indices = frame_range(peak_x_rel);
inside_y = peak_y >= y_min & peak_y <= y_max;
peak_indices = peak_indices(inside_y);
peak_indices = unique(peak_indices(:));
end

function inside_box = points_in_box(x_values, y_values, box_position)
x_min = box_position(1);
x_max = box_position(1) + box_position(3);
y_min = box_position(2);
y_max = box_position(2) + box_position(4);
inside_box = x_values >= x_min & x_values <= x_max & y_values >= y_min & y_values <= y_max;
end

function plot_trend_band(x_values, mean_values, std_values, band_color)
x_values = x_values(:)';
mean_values = mean_values(:)';
std_values = std_values(:)';
valid = isfinite(x_values) & isfinite(mean_values) & isfinite(std_values);
if any(valid)
    x_valid = x_values(valid);
    mean_valid = mean_values(valid);
    std_valid = abs(std_values(valid));
    fill([x_valid fliplr(x_valid)], ...
        [mean_valid + std_valid fliplr(mean_valid - std_valid)], ...
        band_color, 'EdgeColor', 'none', 'FaceAlpha', 0.35);
end
plot(x_values, mean_values, 'k', 'LineWidth', 2);
end

function peak_table = empty_peak_table()
peak_table = table( ...
    zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), zeros(0,1), ...
    strings(0,1), strings(0,1), zeros(0,1), strings(0,1), NaT(0,1), ...
    'VariableNames', {'peak_id', 'roi', 'index', 'time_s', 'polarity', ...
    'amplitude_sensitivity', 'status', 'source', 'parent_peak_id', ...
    'created_stage', 'created_at'});
end

function polarity = choose_peak_polarity(trace_i, mode)
mode = lower(string(mode));
switch mode
    case "positive"
        polarity = 1;
    case "negative"
        polarity = -1;
    otherwise
        if abs(max(trace_i) - mean(trace_i, 'omitnan')) >= abs(min(trace_i) - mean(trace_i, 'omitnan'))
            polarity = 1;
        else
            polarity = -1;
        end
end
end

function [peaks_index, peaks_amplitude, peaks_polarity] = peak_table_to_peak_cells(peak_table, nrois, fallback_polarity)
if nargin < 3
    fallback_polarity = cell(1, nrois);
end
peaks_index = cell(1, nrois);
peaks_amplitude = cell(1, nrois);
peaks_polarity = cell(1, nrois);

for roi_idx = 1:nrois
    accepted = peak_table.roi == roi_idx & peak_table.status == "accepted";
    rows_i = peak_table(accepted, :);
    if height(rows_i) > 0
        [~, order] = sort(rows_i.index);
        rows_i = rows_i(order, :);
        peaks_index{roi_idx} = rows_i.index;
        peaks_amplitude{roi_idx} = rows_i.amplitude_sensitivity;
        peaks_polarity{roi_idx} = rows_i.polarity(1);
    else
        peaks_index{roi_idx} = [];
        peaks_amplitude{roi_idx} = [];
        if numel(fallback_polarity) >= roi_idx && ~isempty(fallback_polarity{roi_idx})
            peaks_polarity{roi_idx} = fallback_polarity{roi_idx};
        else
            peaks_polarity{roi_idx} = 1;
        end
    end
end
end

function peak_table = add_manual_peak_row(peak_table, peak_id, roi_idx, index, amplitude, polarity, freq)
new_row = table(peak_id, roi_idx, index, index / freq, polarity, amplitude, "accepted", ...
    "manual_add", NaN, "sensitivity_manually_edited", datetime("now"), ...
    'VariableNames', peak_table.Properties.VariableNames);
peak_table = [peak_table; new_row];
end

function peak_table = reset_roi_peak_table(peak_table, reset_table, roi_idx)
peak_table(peak_table.roi == roi_idx, :) = [];
peak_table = [peak_table; reset_table(reset_table.roi == roi_idx, :)];
end

function entry = make_edit_history(action, roi_idx, peak_id, index_before, index_after, mode)
entry = struct( ...
    'action', string(action), ...
    'roi', roi_idx, ...
    'peak_id', peak_id, ...
    'index_before', index_before, ...
    'index_after', index_after, ...
    'mode', string(mode), ...
    'created_at', datetime("now"));
end

function [movie_info, trace_results, peak_results, ap_results, ctx] = load_saved_ap_analysis_only_context(reuse_results_path)
reuse_results_path = char(string(reuse_results_path));
if isempty(reuse_results_path) || ~isfolder(reuse_results_path)
    error('AP_analysis3:AnalysisOnlyPathMissing', ...
        'analysis_mode=analysis_only requires reuse_results_path or save_path to point to an existing AP_analysis3 result folder.');
end

movie_info_file = fullfile(reuse_results_path, 'movie_info.mat');
trace_results_file = fullfile(reuse_results_path, 'trace_results.mat');
if ~isfile(movie_info_file) || ~isfile(trace_results_file)
    error('AP_analysis3:AnalysisOnlyMissingResults', ...
        'analysis_only requires movie_info.mat and trace_results.mat in %s.', reuse_results_path);
end

s = load(movie_info_file, 'movie_info');
movie_info = s.movie_info;
s = load(trace_results_file, 'trace_results');
trace_results = s.trace_results;

peak_results = struct();
peak_results_file = fullfile(reuse_results_path, 'peak_results.mat');
if isfile(peak_results_file)
    s = load(peak_results_file, 'peak_results');
    peak_results = s.peak_results;
end

ap_results = struct();
ap_results_file = fullfile(reuse_results_path, 'ap_results.mat');
if isfile(ap_results_file)
    s = load(ap_results_file, 'ap_results');
    ap_results = s.ap_results;
end

ctx = struct();
ctx.traces_raw = get_ap_trace_stage_data(trace_results, 'raw', []);
ctx.traces_background_removed = get_ap_trace_stage_data(trace_results, 'bg_removed', []);

roi_file = fullfile(reuse_results_path, '1_raw_ROI.mat');
ctx.roi_results_file = roi_file;
ctx.avg_image = [];
ctx.rois = struct();
if ~isfile(roi_file)
    error('AP_analysis3:AnalysisOnlyMissingROIResults', ...
        'analysis_only starts after ROI selection and requires %s.', roi_file);
end
roi_data = load(roi_file);
if isfield(roi_data, 'avg_image'), ctx.avg_image = roi_data.avg_image; end
if isfield(roi_data, 'rois'), ctx.rois = roi_data.rois; end
if isempty(ctx.traces_raw) && isfield(roi_data, 'traces')
    ctx.traces_raw = roi_data.traces;
end

map_file = fullfile(reuse_results_path, '0_Sensitivity_Map.mat');
ctx.map = [];
if isfile(map_file)
    map_data = load(map_file);
    if isfield(map_data, 'map'), ctx.map = map_data.map; end
end

mask_file = fullfile(reuse_results_path, '0_cellpose_mask.mat');
ctx.mask = [];
if isfile(mask_file)
    mask_data = load(mask_file);
    if isfield(mask_data, 'mask'), ctx.mask = mask_data.mask; end
end

required_fields = {'traces_raw', 'traces_background_removed'};
missing = required_fields(cellfun(@(name) isempty(ctx.(name)), required_fields));
if ~isempty(missing)
    error('AP_analysis3:AnalysisOnlyIncompleteTraceResults', ...
        'analysis_only starts after background removal and requires saved trace stages: %s.', strjoin(missing, ', '));
end
end

function trace_results = keep_ap_trace_results(trace_results, stage_names)
stage_names = string(stage_names);
if ~isstruct(trace_results)
    trace_results = struct();
    return;
end
source_results = trace_results;
trace_results = struct();
for i = 1:numel(stage_names)
    stage = char(stage_names(i));
    if isfield(source_results, stage)
        trace_results.(stage) = source_results.(stage);
    elseif isfield(source_results, 'trace_results') && isfield(source_results.trace_results, stage)
        trace_results.(stage) = source_results.trace_results.(stage);
    end
end
end

function data = get_ap_trace_stage_data(trace_results, stage_name, default_value)
data = default_value;
if isstruct(trace_results) && isfield(trace_results, 'trace_results') && ...
        isfield(trace_results.trace_results, stage_name) && ...
        isfield(trace_results.trace_results.(stage_name), 'data')
    data = trace_results.trace_results.(stage_name).data;
elseif isstruct(trace_results) && isfield(trace_results, stage_name) && ...
        isfield(trace_results.(stage_name), 'data')
    data = trace_results.(stage_name).data;
end
end

function [poolObj, poolReady, poolInfo] = initialize_ap_parallel_pool(poolInfo)
poolObj = [];
poolReady = false;

if exist('gcp', 'file') ~= 2
    poolInfo.message = "Parallel Computing Toolbox is not available.";
    warning('AP_analysis3:ParallelUnavailable', '%s Continuing in serial mode.', poolInfo.message);
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
    warning('AP_analysis3:ParpoolUnavailable', ...
        ['Parallel pool could not be started (%s). ' ...
        'Continuing in serial mode.'], ME.message);
end
end

function results = store_trace_stage(results, stage_name, data, parent_results, roi_file, movie_info, method, parameters)
% Save one named trace stage together with the minimum provenance needed
% to understand how that stage was produced.
%
% Natural-language meaning:
%   this is the common bookkeeping step used after a new trace
%   representation has been computed.
%
% Interface example:
%   results.trace_results.snr.data
%   results.trace_results.snr.info.parent_results -> {'bleach_removed','noise'}
%   results.trace_results.snr.info.method         -> 'signal_divided_by_noise_std'
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

function record = build_section_record(section_name, description, input_struct, parameter_struct, output_struct, excluded_inputs_struct, rerun_note)
% Build a compact "interface description" for one analysis section.
%
% Natural-language meaning:
%   each record is a small summary of what the section did, what went in,
%   what came out, and what would still be needed to rerun it.
%
% Interface:
%   record.input           -> saved non-movie inputs used by this section
%   record.parameters      -> algorithm choices and tunable settings
%   record.output          -> direct outputs of this section
%   record.excluded_inputs -> required inputs that were intentionally not saved
%   record.rerun.note      -> plain-language rerun instruction
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

function formula_text = describe_bleach_formula(bleach_mode)
switch lower(string(bleach_mode))
    case "linear"
        formula_text = "bleach_removed = detrend(bg_removed, 1); baseline = bg_removed - bleach_removed";
    case "highpass"
        formula_text = "bleach_removed = highpass_bleach_remove(bg_removed, freq, fc); baseline returned by highpass helper";
    case "exp2"
        formula_text = "baseline = fit_exp2(bg_removed); bleach_removed = bg_removed - baseline";
    otherwise
        formula_text = "unknown bleach formula";
end
end

function out = ternary(condition, true_value, false_value)
if condition
    out = true_value;
else
    out = false_value;
end
end

function [mask, info] = run_volpy_roi_segmentation(input_movie_path, save_path, frame_rate)
python_exe = 'C:\Users\DELL\anaconda3\envs\caiman\python.exe';
script_path = fullfile(fileparts(mfilename('fullpath')), 'python_seg', 'segment_voltage_imaging_volpy.py');
output_dir = fullfile(save_path, '0_volpy_roi');

if ~isfile(python_exe)
    error('VolPy ROI python interpreter not found: %s', python_exe);
end
if ~isfile(script_path)
    error('VolPy ROI script not found: %s', script_path);
end
if ~isfolder(output_dir)
    mkdir(output_dir);
end

cmd = sprintf('"%s" "%s" "%s" --output-dir "%s" --frame-rate %.12g', ...
    python_exe, script_path, input_movie_path, output_dir, frame_rate);
[status, cmdout] = system(cmd);
if status ~= 0
    error('VolPy ROI command failed:\n%s', cmdout);
end

mask_file = fullfile(output_dir, 'mask.mat');
if ~isfile(mask_file)
    error('VolPy ROI did not create mask.mat in %s', output_dir);
end

mask_data = load(mask_file, 'mask');
mask = mask_data.mask;
info = struct( ...
    'input_movie_path', string(input_movie_path), ...
    'output_dir', string(output_dir), ...
    'frame_rate', frame_rate, ...
    'script_path', string(script_path), ...
    'command_output', string(cmdout));
end

function info = run_volpy_backend_pipeline(input_movie_path, save_path, frame_rate, flip_signal, skip_motion_correction, roi_mask_file, auto_run, force_rerun, use_existing, output_dir, result_mat_path, temp_dir_override, size_min, size_max, corr_window_seconds, corr_stride_seconds, corr_baseline_seconds, corr_remove_baseline, corr_gaussian_blur, summary_mode, mrcnn_input_mode, mrcnn_confidence_threshold)
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
    corr_remove_baseline_text = ternary(corr_remove_baseline, 'true', 'false');
    corr_gaussian_blur_text = ternary(corr_gaussian_blur, 'true', 'false');
    skip_motion_text = ternary(skip_motion_correction, 'true', 'false');
    cmd = sprintf(['"%s" "%s" "%s" --output-dir "%s" --frame-rate %.12g --flip-signal %s ' ...
        '--size-min %d --size-max %d --corr-window-sec %.12g --corr-stride-sec %.12g ' ...
        '--corr-baseline-sec %.12g --corr-remove-baseline %s --corr-gaussian-blur %s ' ...
        '--summary-mode "%s" --mrcnn-input-mode "%s" --mrcnn-confidence-threshold %.12g ' ...
        '--skip-motion-correction %s'], ...
        python_exe, script_path, input_movie_path, output_dir, frame_rate, flip_text, ...
        round(size_min), round(size_max), corr_window_seconds, corr_stride_seconds, ...
        corr_baseline_seconds, corr_remove_baseline_text, corr_gaussian_blur_text, ...
        summary_mode, mrcnn_input_mode, mrcnn_confidence_threshold, skip_motion_text);
    if strlength(string(roi_mask_file)) > 0
        cmd = sprintf('%s --roi-mask "%s"', cmd, char(roi_mask_file));
    end
    if strlength(string(temp_dir_override)) > 0
        cmd = sprintf('%s --caiman-temp-dir "%s"', cmd, char(temp_dir_override));
    end
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
    'temp_dir_override', string(temp_dir_override), ...
    'manifest_path', string(manifest_path), ...
    'motion_corrected_file', motion_corrected_file, ...
    'memmap_file', memmap_file, ...
    'frame_rate', frame_rate, ...
    'flip_signal', flip_signal, ...
    'skip_motion_correction', skip_motion_correction, ...
    'roi_mask_file', string(roi_mask_file), ...
    'size_min', size_min, ...
    'size_max', size_max, ...
    'corr_window_seconds', corr_window_seconds, ...
    'corr_stride_seconds', corr_stride_seconds, ...
    'corr_baseline_seconds', corr_baseline_seconds, ...
    'corr_remove_baseline', corr_remove_baseline, ...
    'corr_gaussian_blur', corr_gaussian_blur, ...
    'summary_mode', string(summary_mode), ...
    'mrcnn_input_mode', string(mrcnn_input_mode), ...
    'mrcnn_confidence_threshold', mrcnn_confidence_threshold, ...
    'script_path', string(script_path), ...
    'command_output', string(cmdout));
end

function data = load_volpy_backend_results(result_mat_path)
s = load(result_mat_path);

data = struct();
data.mask = uint16(s.mask);
if isfield(s, 'flip_signal')
    data.flip_signal = logical(s.flip_signal(1));
else
    data.flip_signal = false;
end
data.native_sign = ternary(data.flip_signal, -1, 1);

data.t = data.native_sign * cell_columns_to_matrix(s.t);
data.ts = data.native_sign * cell_columns_to_matrix(s.ts);
data.t_rec = data.native_sign * cell_columns_to_matrix(s.t_rec);
data.t_sub = data.native_sign * cell_columns_to_matrix(s.t_sub);
if isfield(s, 'dFF')
    data.dff = data.native_sign * cell_columns_to_matrix(s.dFF);
else
    data.dff = data.t;
end
if isfield(s, 'F0')
    data.f0 = cell_columns_to_matrix(s.F0);
else
    data.f0 = ones(size(data.t));
end

nframes = size(data.t, 1);
nrois = size(data.t, 2);
data.event_trace = data.t;
data.event_sensitivity = data.dff;
valid_f0 = abs(data.f0) > eps;
data.noise = data.t - data.t_rec;
data.snr_trace = zeros(size(data.t));
data.qc_event_trace = data.t_rec;
data.qc_event_sensitivity = zeros(size(data.t_rec));
data.qc_event_sensitivity(valid_f0) = data.t_rec(valid_f0) ./ data.f0(valid_f0);
data.qc_snr_trace = zeros(size(data.t_rec));
for i = 1:nrois
    noise_std = std(data.noise(:, i), 0, 1);
    if ~isfinite(noise_std) || noise_std <= eps
        noise_std = 1;
    end
    data.snr_trace(:, i) = data.event_trace(:, i) ./ noise_std;
    data.qc_snr_trace(:, i) = data.qc_event_trace(:, i) ./ noise_std;
end

data.spikes = cell(1, nrois);
data.peak_amplitude = cell(1, nrois);
data.peaks_polarity = cell(1, nrois);
for i = 1:nrois
    idx = round(double(s.spikes{i}(:)));
    idx = idx(idx >= 1 & idx <= nframes);
    data.spikes{i} = idx;
    data.peak_amplitude{i} = data.t(idx, i);

    polarity = data.native_sign;
    if isfield(s, 'templates') && numel(s.templates) >= i && ~isempty(s.templates{i})
        template_i = double(s.templates{i}(:));
        template_i = template_i(~isnan(template_i));
        if ~isempty(template_i)
            polarity = sign(sum(template_i)) * data.native_sign;
        end
    end
    if ~isfinite(polarity) || polarity == 0
        polarity = data.native_sign;
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

mean_image_path = fullfile(fileparts(result_mat_path), 'mean_image.tif');
if isfile(mean_image_path)
    avg_image = double(imread(mean_image_path));
    avg_min = min(avg_image(:));
    avg_max = max(avg_image(:));
    if avg_max > avg_min
        avg_image = (avg_image - avg_min) ./ (avg_max - avg_min);
    else
        avg_image = zeros(size(avg_image));
    end
    data.avg_image = avg_image;
else
    data.avg_image = double(data.mask > 0);
end
end

function matrix = cell_columns_to_matrix(cell_values)
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

function save_volpy_backend_compat_outputs(save_path, avg_image, mask, traces_raw, traces_sub, traces_rec, traces_sensitivity, traces_snr, t, peaks_index)
%SAVE_VOLPY_BACKEND_COMPAT_OUTPUTS Recreate the most useful classic-stage figures for volpy backend runs.

% 0_volpy_mask
fig = figure('Name', 'VolPy ROI Mask', 'Color', 'w');
overlay_image = labeloverlay(avg_image, mask, 'Transparency', 0.6);
imshow(overlay_image);
title('volpy ROI mask');
saveas(fig, fullfile(save_path, '0_volpy_mask.fig'), 'fig');
saveas(fig, fullfile(save_path, '0_volpy_mask.png'), 'png');
save(fullfile(save_path, '0_volpy_mask.mat'), 'mask');
close(fig);

% 0_average_roi_overlay
fig = figure('Name', 'Average With ROI Labels', 'Color', 'w');
imshow(avg_image, [], 'InitialMagnification', 'fit');
hold on;
stats = regionprops(mask, 'Centroid');
num_rois = max(mask(:));
if num_rois < 1
    num_rois = 1;
end
colors = lines(num_rois);
for roi_idx = 1:max(mask(:))
    roi_mask = mask == roi_idx;
    roi_boundaries = bwboundaries(roi_mask, 'noholes');
    current_color = colors(mod(roi_idx - 1, size(colors, 1)) + 1, :);
    for b_idx = 1:numel(roi_boundaries)
        boundary = roi_boundaries{b_idx};
        plot(boundary(:, 2), boundary(:, 1), 'Color', current_color, 'LineWidth', 1.2);
    end
    if roi_idx <= numel(stats) && ~isempty(stats(roi_idx).Centroid)
        centroid = stats(roi_idx).Centroid;
        text(centroid(1), centroid(2), sprintf('%d', roi_idx), ...
            'Color', 'y', 'FontWeight', 'bold', 'FontSize', 10, ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
end
title('Average image with ROI labels');
saveas(fig, fullfile(save_path, '0_average_roi_overlay.fig'), 'fig');
saveas(fig, fullfile(save_path, '0_average_roi_overlay.png'), 'png');
close(fig);

% 1_raw_trace
fig = figure('Name', 'Raw Trace (VolPy)', 'Color', 'w');
[~] = offset_plot(traces_raw, t);
title('Raw Trace (VolPy backend import)');
saveas(fig, fullfile(save_path, '1_raw_trace.fig'), 'fig');
saveas(fig, fullfile(save_path, '1_raw_trace.png'), 'png');
close(fig);

% 1_Averaged.tif
avg_uint16 = uint16(round(mat2gray(avg_image) * 65535));
tifffile = Tiff(fullfile(save_path, '1_Averaged.tif'), 'w');
tifffile.setTag('ImageLength', size(avg_uint16, 1));
tifffile.setTag('ImageWidth', size(avg_uint16, 2));
tifffile.setTag('Photometric', Tiff.Photometric.MinIsBlack);
tifffile.setTag('BitsPerSample', 16);
tifffile.setTag('SamplesPerPixel', 1);
tifffile.setTag('RowsPerStrip', 16);
tifffile.setTag('PlanarConfiguration', Tiff.PlanarConfiguration.Chunky);
tifffile.setTag('Compression', Tiff.Compression.None);
tifffile.setTag('Software', 'MATLAB');
tifffile.write(avg_uint16);
tifffile.writeDirectory();
tifffile.close();

% background_correction_summary.png (compatibility summary)
fig = figure('Name', 'VolPy Frontend Summary', 'Color', 'w', 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.7]);
subplot(1, 2, 1);
imshow(overlay_image);
title('VolPy ROI overlay');
subplot(1, 2, 2);
hold on;
plot(t, traces_raw, 'Color', [0.7 0.7 0.7]);
plot(t, traces_sub, '--', 'LineWidth', 1);
plot(t, traces_rec, 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Signal');
title('VolPy trace decomposition');
grid on;
saveas(fig, fullfile(save_path, 'background_correction_summary.png'));
close(fig);

% 2_bleach_correction_stacked (compatibility figure)
fig = figure('Name', 'VolPy Trace Decomposition', 'Color', 'w');
ax_left = subplot(1, 2, 1);
hold(ax_left, 'on');
spacing = max(std(traces_raw, 0, 1), [], 'all');
if ~isfinite(spacing) || spacing <= 0
    spacing = 1;
end
for r = 1:size(traces_raw, 2)
    offset = (size(traces_raw, 2) - r) * spacing * 8;
    plot(ax_left, t, traces_raw(:, r) + offset, 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
    plot(ax_left, t, traces_sub(:, r) + offset, '--', 'Color', [0.2 0.45 0.9], 'LineWidth', 0.8);
end
title(ax_left, 'VolPy raw and subthreshold');
xlabel(ax_left, 'Time (s)');
ylabel(ax_left, 'Stacked magnitude');
grid(ax_left, 'on');
axis(ax_left, 'tight');

ax_right = subplot(1, 2, 2);
hold(ax_right, 'on');
for r = 1:size(traces_raw, 2)
    offset = (size(traces_raw, 2) - r) * spacing * 8;
    plot(ax_right, t, traces_raw(:, r) + offset, 'Color', [0.75 0.75 0.75], 'LineWidth', 0.5);
    plot(ax_right, t, traces_rec(:, r) + offset, 'Color', [0.9 0.15 0.15], 'LineWidth', 0.8);
end
title(ax_right, 'VolPy raw and reconstructed spikes');
xlabel(ax_right, 'Time (s)');
ylabel(ax_right, 'Stacked magnitude');
grid(ax_right, 'on');
axis(ax_right, 'tight');
saveas(fig, fullfile(save_path, '2_bleach_correction_stacked.fig'), 'fig');
saveas(fig, fullfile(save_path, '2_bleach_correction_stacked.png'), 'png');
close(fig);

% 4_SNR
fig = figure('Name', 'VolPy SNR Summary', 'Color', 'w');
subplot(1, 3, 1);
title('raw');
hold on;
[~] = offset_plot(traces_raw, t);
subplot(1, 3, 2);
title('Sensitivity');
hold on;
[~] = offset_plot(traces_sensitivity, t);
subplot(1, 3, 3);
title('SNR');
hold on;
[~] = offset_plot(traces_snr, t);
saveas(fig, fullfile(save_path, '4_SNR.fig'), 'fig');
saveas(fig, fullfile(save_path, '4_SNR.png'), 'png');
save(fullfile(save_path, '4_SNR.mat'), 'traces_raw', 'traces_sensitivity', 'traces_snr');
close(fig);

% 4_Sensitivity_stack
fig = figure('Name', 'VolPy Sensitivity Stack', 'Color', 'w');
[~] = offset_plot(traces_sensitivity, t);
title('Sensitivity');
saveas(fig, fullfile(save_path, '4_Sensitivity_stack.fig'), 'fig');
saveas(fig, fullfile(save_path, '4_Sensitivity_stack.png'), 'png');
save(fullfile(save_path, '4_Sensitivity_stack.mat'), 'traces_sensitivity');
close(fig);

% 4_Sensitivity_stack_with_AP_ticks
fig = figure('Name', 'VolPy Sensitivity Stack with AP ticks', 'Color', 'w');
ax = axes(fig);
hold(ax, 'on');
offset_array = offset_plot(traces_sensitivity, t);
add_peak_tick_marks(ax, traces_sensitivity, t, peaks_index, offset_array);
title(ax, 'Sensitivity with detected AP ticks');
xlabel(ax, 'Time (s)');
ylabel(ax, 'Stacked Magnitude');
saveas(fig, fullfile(save_path, '4_Sensitivity_stack_with_AP_ticks.fig'), 'fig');
saveas(fig, fullfile(save_path, '4_Sensitivity_stack_with_AP_ticks.png'), 'png');
save(fullfile(save_path, '4_Sensitivity_stack_with_AP_ticks.mat'), 'traces_sensitivity', 'peaks_index');
close(fig);

% 4_SNR_stack
fig = figure('Name', 'VolPy SNR Stack', 'Color', 'w');
[~] = offset_plot(traces_snr, t);
title('SNR');
saveas(fig, fullfile(save_path, '4_SNR_stack.fig'), 'fig');
saveas(fig, fullfile(save_path, '4_SNR_stack.png'), 'png');
save(fullfile(save_path, '4_SNR_stack.mat'), 'traces_snr');
close(fig);

% 4_SNR_stack_with_AP_ticks
fig = figure('Name', 'VolPy SNR Stack with AP ticks', 'Color', 'w');
ax = axes(fig);
hold(ax, 'on');
offset_array = offset_plot(traces_snr, t);
add_peak_tick_marks(ax, traces_snr, t, peaks_index, offset_array);
title(ax, 'SNR with detected AP ticks');
xlabel(ax, 'Time (s)');
ylabel(ax, 'Stacked Magnitude');
saveas(fig, fullfile(save_path, '4_SNR_stack_with_AP_ticks.fig'), 'fig');
saveas(fig, fullfile(save_path, '4_SNR_stack_with_AP_ticks.png'), 'png');
save(fullfile(save_path, '4_SNR_stack_with_AP_ticks.mat'), 'traces_snr', 'peaks_index');
close(fig);

% Peak overview image for easier QC in volpy mode
fig = figure('Name', 'VolPy Peak Overview', 'Color', 'w');
nrois = size(traces_raw, 2);
for i = 1:nrois
    subplot(nrois, 1, i);
    plot(t, traces_raw(:, i), 'k', 'LineWidth', 0.8);
    hold on;
    idx = peaks_index{i};
    idx = idx(idx >= 1 & idx <= size(traces_raw, 1));
    if ~isempty(idx)
        scatter(t(idx), traces_raw(idx, i), 12, 'r', 'filled');
    end
    title(sprintf('ROI %d peaks (VolPy backend)', i));
    xlabel('Time (s)');
    ylabel('Signal');
    grid on;
end
saveas(fig, fullfile(save_path, '4_peak_detect_overview.png'), 'png');
close(fig);
end

function add_peak_tick_marks(ax, traces, t, peaks_index, offset_array)
%ADD_PEAK_TICK_MARKS Draw short black lines at detected AP peak positions on a stacked trace plot.
hold(ax, 'on');
if isempty(traces)
    return;
end
trace_ranges = max(traces, [], 1) - min(traces, [], 1);
global_range = median(trace_ranges(isfinite(trace_ranges)));
if ~isfinite(global_range) || global_range <= 0
    global_range = 1;
end

for roi_idx = 1:size(traces, 2)
    idx = peaks_index{roi_idx};
    idx = idx(idx >= 1 & idx <= size(traces, 1));
    if isempty(idx)
        continue;
    end
    trace_i = traces(:, roi_idx);
    local_range = max(trace_i) - min(trace_i);
    if ~isfinite(local_range) || local_range <= 0
        local_range = global_range;
    end
    tick_half_height = max(global_range * 0.05, local_range * 0.12);
    y_center = trace_i(idx) + offset_array(roi_idx);
    x_values = t(idx);
    for peak_idx = 1:numel(idx)
        line(ax, [x_values(peak_idx) x_values(peak_idx)], ...
            [y_center(peak_idx) - tick_half_height, y_center(peak_idx) + tick_half_height], ...
            'Color', 'k', 'LineWidth', 0.6, 'Clipping', 'on');
    end
end
end

function segment_out = baseline_subtract_event_window(segment_in, center_idx)
segment_out = segment_in;
if isempty(segment_in)
    return;
end
center_idx = max(1, min(numel(segment_in), round(center_idx)));
pre_end = max(1, center_idx - 1);
if pre_end < 1
    return;
end
pre_segment = segment_in(1:pre_end);
pre_segment = pre_segment(isfinite(pre_segment));
if isempty(pre_segment)
    return;
end
baseline = median(pre_segment, 'omitnan');
if ~isfinite(baseline)
    return;
end
segment_out = segment_in - baseline;
end


function [movie, ncols, nrows, nframes] = load_movie_with_fallback(file_path)
%LOAD_MOVIE_WITH_FALLBACK Prefer the original load_movie path and only fall back
% to serial TIFF loading when parallel pool startup or parallel TIFF access fails.
try
    [movie, ncols, nrows, nframes] = load_movie(file_path);
catch ME
    warning('AP_analysis3:loadMovieFallback', ...
        'load_movie failed (%s). Falling back to serial TIFF loading.', ME.message);
    info = imfinfo(file_path);
    nframes = numel(info);
    ncols = info(1).Height;
    nrows = info(1).Width;
    movie = zeros(ncols * nrows, nframes, 'uint16');
    for frame_idx = 1:nframes
        frame = imread(file_path, frame_idx, 'Info', info);
        movie(:, frame_idx) = reshape(frame, ncols * nrows, 1);
    end
end
end


function Y_hp = create_temp_highpass(movie)
%CREATE_TEMP_HIGHPASS Build a temporary high-pass filtered movie.
% Used only during motion estimation to improve shift estimation quality.
gSig = 7;
gSiz = 17;
psf = fspecial('gaussian', round(2 * gSiz), gSig);
ind = (psf >= max(psf(:,1)));
psf = psf - mean(psf(ind));
psf(~ind) = 0;
Y_hp = imfilter(movie, psf, 'symmetric');
end

function plot_corrected(traces_bleaching_removed, traces_background_removed, baseline, t, save_path)
%PLOT_CORRECTED Save a two-panel summary of baseline fitting and bleaching removal.
nrois = size(traces_bleaching_removed, 2);

fig = figure('Name', 'Bleaching Correction Overview', 'Color', 'w');
set(fig, 'Position', get(0, 'Screensize'));

ax_fit = subplot(1, 2, 1);
hold(ax_fit, 'on');
spacing_raw = mean(std(traces_background_removed, 0, 1)) * 5;

for r = 1:nrois
    offset = (nrois - r) * spacing_raw;
    plot(t, traces_background_removed(:, r) + offset, 'LineWidth', 0.5, 'DisplayName', 'Original');
    plot(t, baseline(:, r) + offset, 'r', 'LineWidth', 1.2, 'DisplayName', 'Baseline');

    if mod(r, 5) == 0 || r == 1 || r == nrois
        text(t(1), offset, [' ROI ', num2str(r)], 'FontSize', 8, 'FontWeight', 'bold');
    end
end
title(['Original Traces & Fitted Baselines (N=', num2str(nrois), ')']);
xlabel('Time (s)');
ylabel('Stacked Magnitude');
grid on;
axis tight;

ax_corr = subplot(1, 2, 2);
hold(ax_corr, 'on');
spacing_corr = mean(std(traces_bleaching_removed, 0, 1)) * 8;

for r = 1:nrois
    offset = (nrois - r) * spacing_corr;
    plot(t, traces_bleaching_removed(:, r) + offset, 'k', 'LineWidth', 0.5);
    line([t(1) t(end)], [offset offset], 'Color', [0.3 0.7 1, 0.5], 'LineStyle', '--');
end
title('Bleaching-Removed Traces (Residuals)');
xlabel('Time (s)');
ylabel('Stacked Magnitude');
grid on;
axis tight;

linkaxes([ax_fit, ax_corr], 'x');

fig_filename = fullfile(save_path, '2_bleach_correction_stacked.fig');
png_filename = fullfile(save_path, '2_bleach_correction_stacked.png');
saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');
end
