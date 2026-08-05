function [ctx, outcome] = motion(ctx, action)
%MOTION Coordinate shared single-channel motion algorithms.
%   One-channel workflows estimate/apply one model. In a voltage/calcium
%   workflow, voltage is the reference: equal frame counts share its shift
%   field; unequal frame counts estimate or reuse a separate calcium model.
%   No shift field is truncated, padded, or temporally interpolated.

action = normalize_action(action);
outcome = struct( ...
    'section', "motion", ...
    'action', action, ...
    'state', "completed", ...
    'channels', strings(0, 1), ...
    'message', "");

if action == "skip"
    outcome.state = "skipped";
    outcome.message = "Motion correction skipped.";
    return;
end
[channels, roles] = validate_context(ctx);
outcome.channels = roles(:);

% Load inputs: select the reference/secondary movies and bound profiles.
reference_idx = find(roles == "voltage", 1, 'first');
if isempty(reference_idx)
    reference_idx = 1;
end
secondary_idx = setdiff(1:numel(channels), reference_idx, 'stable');
if numel(channels) == 2 && roles(secondary_idx) ~= "calcium"
    error('AAA:Sections:InvalidDualMotionRoles', ...
        'Two-channel motion requires one voltage and one calcium profile.');
end

reference_movie = channels(reference_idx).data.movie_3d;
reference_profile = channels(reference_idx).profile;
reference_frames = size(reference_movie, 3);
reference_role=roles(reference_idx);

% Run algorithm: resolve models first, then apply them without alignment.
aaa.sections.emit_phase(ctx,"motion","resolve_model","running", ...
    "Estimating or loading motion model.",reference_role);
reference_model = resolve_model( ...
    action, reference_movie, reference_profile, ctx, true);

secondary_model = struct();
shared_model = false;
if ~isempty(secondary_idx)
    secondary_movie = channels(secondary_idx).data.movie_3d;
    if ~isequal(spatial_size(reference_movie), spatial_size(secondary_movie))
        error('AAA:Sections:MotionGeometryMismatch', ...
            'Channel movies must share spatial geometry before motion correction.');
    end
    secondary_frames = size(secondary_movie, 3);
    shared_model = reference_frames == secondary_frames;
    if shared_model
        secondary_model = reference_model;
        aaa.sections.emit_phase(ctx,"motion","share_model","completed", ...
            "Frame counts match; sharing the reference shift field.",roles(secondary_idx));
    else
        aaa.sections.emit_phase(ctx,"motion","resolve_model","running", ...
            "Estimating or loading independent motion model.",roles(secondary_idx));
        secondary_model = resolve_model( ...
            action, secondary_movie, channels(secondary_idx).profile, ctx, false);
    end
end

% Resolve and validate every required model before applying any correction.
aaa.sections.emit_phase(ctx,"motion","apply_motion","running", ...
    "Applying motion correction.",reference_role);
[reference_corrected, reference_application] = ...
    apply_motion(reference_movie, reference_model);
if ~isempty(secondary_idx)
    aaa.sections.emit_phase(ctx,"motion","apply_motion","running", ...
        "Applying motion correction.",roles(secondary_idx));
    [secondary_corrected, secondary_application] = ...
        apply_motion(secondary_movie, secondary_model);
else
    secondary_corrected = [];
    secondary_application = struct();
end

% Return results: write motion auxiliaries and update channel state.
output_path = resolve_output_path(ctx);
aaa.sections.emit_phase(ctx,"motion","persist","running", ...
    "Saving motion models, metrics, and corrected outputs.");
[reference_paths, secondary_paths] = persist_motion_outputs( ...
    output_path, action, reference_model, secondary_model, ...
    shared_model, reference_profile, ...
    conditional_profile(channels, secondary_idx), ...
    reference_movie, reference_corrected, ...
    conditional_movie(channels, secondary_idx), secondary_corrected);

reference_info = build_reference_info( ...
    reference_model, reference_paths, action, ...
    ~isempty(secondary_idx) && shared_model);
channels(reference_idx) = store_motion_result( ...
    channels(reference_idx), reference_corrected, reference_model, ...
    reference_application, reference_info);

if ~isempty(secondary_idx)
    secondary_info = build_secondary_info( ...
        secondary_model, secondary_paths, reference_paths, action, shared_model);
    channels(secondary_idx) = store_motion_result( ...
        channels(secondary_idx), secondary_corrected, secondary_model, ...
        secondary_application, secondary_info);
end

ctx.channels = channels;
outcome.shared_shift = shared_model;
if isempty(secondary_idx)
    outcome.message = "Estimated/applied one channel motion model.";
elseif shared_model
    outcome.message = "Applied one voltage-estimated shift field to both channels.";
else
    outcome.message = [ ...
        "Frame counts differ; applied independent voltage and calcium " ...
        "motion models without temporal alignment."];
end
end

function action = normalize_action(action)
action = lower(strtrim(string(action)));
if ~isscalar(action) || ~ismember(action, ["run","reuse","skip"])
    error('AAA:Sections:InvalidMotionAction', ...
        'Motion action must be run, reuse, or skip.');
end
end

function [channels, roles] = validate_context(ctx)
if ~isstruct(ctx) || ~isscalar(ctx) ...
        || ~isfield(ctx, 'channels') || isempty(ctx.channels) ...
        || ~isstruct(ctx.channels) || numel(ctx.channels) > 2
    error('AAA:Sections:InvalidMotionContext', ...
        'Motion context must contain one or two channels.');
end
channels = ctx.channels;
roles = strings(1, numel(channels));
for idx = 1:numel(channels)
    channel = channels(idx);
    if ~isfield(channel, 'profile') || ~isstruct(channel.profile) ...
            || ~isscalar(channel.profile) ...
            || ~isfield(channel.profile, 'role') ...
            || strlength(string(channel.profile.role)) == 0
        error('AAA:Sections:InvalidMotionProfile', ...
            'Channel %d must contain one bound profile with role.', idx);
    end
    roles(idx) = lower(string(channel.profile.role));
    if ~isfield(channel, 'data') || ~isstruct(channel.data) ...
            || ~isfield(channel.data, 'movie_3d') ...
            || isempty(channel.data.movie_3d) ...
            || ~isfield(channel.data, 'movie_info')
        error('AAA:Sections:MotionInputMissing', ...
            'Channel %s has no loaded movie input.', roles(idx));
    end
    if ~isfield(channel, 'output_files') || ~isstruct(channel.output_files)
        error('AAA:Sections:InvalidMotionOutputFiles', ...
            'Channel %s output_files must be a struct.', roles(idx));
    end
end
if numel(unique(roles)) ~= numel(roles)
    error('AAA:Sections:DuplicateMotionRole', ...
        'Motion channel roles must be unique.');
end
end

function model = resolve_model(action, movie_3d, profile, ctx, is_reference)
if action == "run"
    model = estimate_motion(movie_3d, profile);
    return;
end
source_file = resolve_shift_source(profile, ctx, is_reference);
model = load_saved_model(source_file, movie_3d, profile, is_reference);
end

function source_file = resolve_shift_source(profile, ctx, is_reference)
source_file = "";
motion_profile = require_motion_profile(profile);
if ~is_reference && isfield(motion_profile, 'saved_calcium_shift_file') ...
        && strlength(string(motion_profile.saved_calcium_shift_file)) > 0
    source_file = string(motion_profile.saved_calcium_shift_file);
elseif isfield(motion_profile, 'saved_shift_file') ...
        && strlength(string(motion_profile.saved_shift_file)) > 0
    source_file = string(motion_profile.saved_shift_file);
end

if strlength(source_file) == 0
    source_root = resolve_source_results_path(ctx);
    if strlength(source_root) > 0
        if is_reference
            candidates = ["shared_motion_shifts_result.mat", ...
                "motion_shifts_result.mat"];
        else
            candidates = "calcium_motion_shifts_result.mat";
        end
        for idx = 1:numel(candidates)
            candidate = string(fullfile(source_root, candidates(idx)));
            if isfile(candidate)
                source_file = candidate;
                break;
            end
        end
    end
end
if strlength(source_file) == 0 || ~isfile(source_file)
    if is_reference
        label = 'reference/shared';
    else
        label = 'independent calcium';
    end
    error('AAA:Sections:MissingMotionShiftFile', ...
        'motion="reuse" requires a valid %s shift file.', label);
end
end

function source_root = resolve_source_results_path(ctx)
source_root = "";
if isfield(ctx, 'source_results_path') ...
        && strlength(string(ctx.source_results_path)) > 0
    source_root = string(ctx.source_results_path);
elseif isfield(ctx, 'request') && isstruct(ctx.request) ...
        && isfield(ctx.request, 'workflow') ...
        && isfield(ctx.request.workflow, 'source_results_path') ...
        && strlength(string(ctx.request.workflow.source_results_path)) > 0
    source_root = string(ctx.request.workflow.source_results_path);
end
end

function model = load_saved_model(source_file, movie_3d, profile, is_reference)
saved = load(char(source_file));
if is_reference
    if ~isfield(saved, 'shifts_r')
        error('AAA:Sections:InvalidReferenceShiftFile', ...
            'Saved reference motion file does not contain shifts_r: %s', ...
            source_file);
    end
    shifts = saved.shifts_r;
    if isfield(saved, 'options_r')
        options = saved.options_r;
    else
        options = build_options(movie_3d, profile);
    end
else
    if isfield(saved, 'shifts_c')
        shifts = saved.shifts_c;
    elseif isfield(saved, 'shifts_r')
        shifts = saved.shifts_r;
    else
        error('AAA:Sections:InvalidCalciumShiftFile', ...
            ['Saved calcium motion file contains neither shifts_c nor ' ...
             'shifts_r: %s'], source_file);
    end
    if isfield(saved, 'options_c')
        options = saved.options_c;
    elseif isfield(saved, 'options_r')
        options = saved.options_r;
    else
        options = build_options(movie_3d, profile);
    end
end

frame_count = size(movie_3d, 3);
if numel(shifts) ~= frame_count
    error('AAA:Sections:SavedShiftFrameMismatch', ...
        ['Saved shifts contain %d frames but the current movie has %d. ' ...
         'The shift field will not be truncated or interpolated.'], ...
        numel(shifts), frame_count);
end
model = struct( ...
    'method', "NoRMCorre_rigid", ...
    'role', string(profile.role), ...
    'shifts', shifts, ...
    'options', options, ...
    'frame_count', frame_count, ...
    'frame_size', spatial_size(movie_3d), ...
    'parameters', require_motion_profile(profile), ...
    'source_shift_file', string(source_file), ...
    'reused', true, ...
    'created_at', datetime('now'));
end

function options = build_options(movie_3d, profile)
params = require_motion_profile(profile);
[ncols, nrows, ~] = size(movie_3d);
required = {'bin_width','max_shift','upsample_factor', ...
    'iterations','correct_bidir'};
for idx = 1:numel(required)
    if ~isfield(params, required{idx}) || isempty(params.(required{idx}))
        error('AAA:Sections:IncompleteMotionProfile', ...
            'profile.motion is missing %s.', required{idx});
    end
end
options = NoRMCorreSetParms( ...
    'd1', ncols, 'd2', nrows, ...
    'bin_width', params.bin_width, ...
    'max_shift', params.max_shift, ...
    'us_fac', params.upsample_factor, ...
    'iter', params.iterations, ...
    'correct_bidir', params.correct_bidir);
end

function params = require_motion_profile(profile)
if ~isstruct(profile) || ~isscalar(profile) ...
        || ~isfield(profile, 'motion') || ~isstruct(profile.motion) ...
        || ~isscalar(profile.motion)
    error('AAA:Sections:IncompleteMotionProfile', ...
        'Bound channel profile must contain one scalar motion struct.');
end
params = profile.motion;
end

function output_path = resolve_output_path(ctx)
output_path = "";
if isfield(ctx, 'output_path') && strlength(string(ctx.output_path)) > 0
    output_path = string(ctx.output_path);
elseif isfield(ctx, 'paths') && isstruct(ctx.paths) ...
        && isfield(ctx.paths, 'output_path') ...
        && strlength(string(ctx.paths.output_path)) > 0
    output_path = string(ctx.paths.output_path);
elseif isfield(ctx, 'request') && isstruct(ctx.request) ...
        && isfield(ctx.request, 'workflow') ...
        && isfield(ctx.request.workflow, 'output_path') ...
        && strlength(string(ctx.request.workflow.output_path)) > 0
    output_path = string(ctx.request.workflow.output_path);
end
if strlength(output_path) > 0 && ~isfolder(output_path)
    mkdir(output_path);
end
end

function [reference_paths, secondary_paths] = persist_motion_outputs( ...
        output_path, action, reference_model, secondary_model, ...
        shared_model, reference_profile, secondary_profile, ...
        reference_raw, reference_corrected, secondary_raw, secondary_corrected)
reference_paths = empty_paths();
secondary_paths = empty_paths();
if strlength(output_path) == 0
    return;
end

shifts_r = reference_model.shifts;
options_r = reference_model.options;
cfg = build_saved_cfg(action, reference_profile, secondary_profile);
reference_paths.shift_file = string(fullfile(output_path, ...
    'shared_motion_shifts_result.mat'));
reference_paths.parameter_file = string(fullfile(output_path, ...
    'shared_motion_correction_para.mat'));
save(reference_paths.shift_file, 'shifts_r', 'options_r', '-v7.3');
save(reference_paths.parameter_file, 'options_r', 'cfg');

if ~isempty(fieldnames(secondary_model)) && ~shared_model
    shifts_c = secondary_model.shifts;
    options_c = secondary_model.options;
    secondary_paths.shift_file = string(fullfile(output_path, ...
        'calcium_motion_shifts_result.mat'));
    secondary_paths.parameter_file = string(fullfile(output_path, ...
        'calcium_motion_correction_para.mat'));
    calcium_motion_context = struct( ...
        'reason', "voltage_calcium_frame_count_mismatch", ...
        'voltage_frame_count', size(reference_raw, 3), ...
        'calcium_frame_count', size(secondary_raw, 3), ...
        'source_shift_file', secondary_model.source_shift_file, ...
        'auto_reused_previous_shift', false, ...
        'created_at', datetime('now'));
    save(secondary_paths.shift_file, 'shifts_c', 'options_c', '-v7.3');
    save(secondary_paths.parameter_file, ...
        'options_c', 'cfg', 'calcium_motion_context');
end

params = require_motion_profile(reference_profile);
save_tif = require_output_flag(params, 'save_downsampled_tif');
plot_metrics = require_output_flag(params, 'plot_metrics');
if save_tif || plot_metrics
    [factor, q_low, q_high] = require_qc_parameters(params);
else
    factor = [];
    q_low = [];
    q_high = [];
end
if save_tif
    reference_paths.downsampled_tif_file = string(fullfile(output_path, ...
        sprintf('voltage_motion_corrected_ds%d.tif', factor)));
    save_downsampled_motion_tif(reference_corrected, factor, ...
        reference_paths.downsampled_tif_file, q_low, q_high);
    if ~isempty(secondary_corrected)
        secondary_paths.downsampled_tif_file = string(fullfile(output_path, ...
            sprintf('calcium_motion_corrected_ds%d.tif', factor)));
        save_downsampled_motion_tif(secondary_corrected, factor, ...
            secondary_paths.downsampled_tif_file, q_low, q_high);
    end
end
if plot_metrics
    reference_paths.metrics_file = string(fullfile(output_path, ...
        'shared_motion_metrics.mat'));
    reference_paths.metrics_fig = string(fullfile(output_path, ...
        'shared_motion_metrics.fig'));
    reference_paths.metrics_png = string(fullfile(output_path, ...
        'shared_motion_metrics.png'));
    save_motion_metrics(reference_raw, reference_corrected, ...
        reference_model, reference_paths, factor, "Shared");
    if ~isempty(secondary_corrected) && ~shared_model
        secondary_paths.metrics_file = string(fullfile(output_path, ...
            'calcium_motion_metrics.mat'));
        secondary_paths.metrics_fig = string(fullfile(output_path, ...
            'calcium_motion_metrics.fig'));
        secondary_paths.metrics_png = string(fullfile(output_path, ...
            'calcium_motion_metrics.png'));
        save_motion_metrics(secondary_raw, secondary_corrected, ...
            secondary_model, secondary_paths, factor, "Calcium");
    end
end
end

function cfg = build_saved_cfg(action, reference_profile, secondary_profile)
params = require_motion_profile(reference_profile);
cfg = struct( ...
    'enabled', true, ...
    'use_saved_shift', action == "reuse", ...
    'saved_shift_file', char(resolve_optional_field(params, 'saved_shift_file', "")), ...
    'saved_calcium_shift_file', '', ...
    'highpass', logical(params.highpass), ...
    'auto_reuse_previous_shift', false);
if ~isempty(fieldnames(secondary_profile))
    secondary_params = require_motion_profile(secondary_profile);
    cfg.saved_calcium_shift_file = char(resolve_optional_field( ...
        secondary_params, 'saved_shift_file', ...
        resolve_optional_field(secondary_params, ...
        'saved_calcium_shift_file', "")));
end
end

function value = resolve_optional_field(target, name, fallback)
if isfield(target, name) && ~isempty(target.(name))
    value = string(target.(name));
else
    value = string(fallback);
end
end

function value = require_output_flag(params, name)
if ~isfield(params, name) || isempty(params.(name)) ...
        || ~(islogical(params.(name)) && isscalar(params.(name)))
    error('AAA:Sections:IncompleteMotionProfile', ...
        'profile.motion.%s must be one public logical value.', name);
end
value = params.(name);
end

function [factor, q_low, q_high] = require_qc_parameters(params)
required = {'downsample_factor','display_quantile_low', ...
    'display_quantile_high'};
for idx = 1:numel(required)
    if ~isfield(params, required{idx}) || isempty(params.(required{idx}))
        error('AAA:Sections:IncompleteMotionProfile', ...
            'profile.motion is missing %s.', required{idx});
    end
end
factor = double(params.downsample_factor);
q_low = double(params.display_quantile_low);
q_high = double(params.display_quantile_high);
if ~isscalar(factor) || ~isfinite(factor) || factor < 1 ...
        || factor ~= round(factor)
    error('AAA:Sections:InvalidMotionQcParameter', ...
        'Motion QC downsample_factor must be a positive integer.');
end
if ~isscalar(q_low) || ~isscalar(q_high) ...
        || ~isfinite(q_low) || ~isfinite(q_high) ...
        || q_low < 0 || q_high > 1 || q_low >= q_high
    error('AAA:Sections:InvalidMotionQcParameter', ...
        'Motion QC quantiles must satisfy 0 <= low < high <= 1.');
end
end

function paths = empty_paths()
paths = struct( ...
    'shift_file', "", ...
    'parameter_file', "", ...
    'downsampled_tif_file', "", ...
    'metrics_file', "", ...
    'metrics_fig', "", ...
    'metrics_png', "");
end

function info = build_reference_info(model, paths, action, shared_with_secondary)
info = common_motion_info(model, paths, action);
info.method = 'NoRMCorre_rigid_shared_voltage_reference';
if shared_with_secondary
    info.shared_with_role = 'calcium';
else
    info.shared_with_role = '';
end
end

function info = build_secondary_info(model, paths, reference_paths, action, shared)
if shared
    info = common_motion_info(model, reference_paths, action);
    info.method = 'reuse_voltage_motion_shifts';
    info.source_role = 'voltage';
    info.frame_count_fallback = false;
else
    info = common_motion_info(model, paths, action);
    info.method = 'NoRMCorre_rigid_independent_calcium';
    info.source_role = 'calcium';
    info.frame_count_fallback = true;
    info.parameters = model.options;
end
end

function info = common_motion_info(model, paths, action)
info = struct( ...
    'applied', true, ...
    'method', '', ...
    'shift_file', char(paths.shift_file), ...
    'source_shift_file', char(model.source_shift_file), ...
    'parameter_file', char(paths.parameter_file), ...
    'downsample_factor', NaN, ...
    'downsampled_tif_file', char(paths.downsampled_tif_file), ...
    'metrics_file', char(paths.metrics_file), ...
    'metrics_fig', char(paths.metrics_fig), ...
    'metrics_png', char(paths.metrics_png), ...
    'highpass', logical(model.parameters.highpass), ...
    'use_saved_shift', action == "reuse", ...
    'auto_reused_previous_shift', false);
if isfield(model.parameters, 'downsample_factor')
    info.downsample_factor = double(model.parameters.downsample_factor);
end
end

function channel = store_motion_result( ...
        channel, corrected, model, application, motion_info)
channel.data.movie_3d = corrected;
channel.data.movie_info.motion = motion_info;
channel.data.movie_info.analysis_frame_size = spatial_size(corrected);
channel.data.movie_info.frame_count = size(corrected, 3);
channel.data.movie_info.updated_at = datetime('now');
channel.results.movie_info = channel.data.movie_info;
channel.data.time = (1:size(corrected, 3))' ...
    / channel.data.movie_info.frame_rate;
channel.profile.roi.frame_size = spatial_size(corrected);
channel.results.motion = struct('model',model,'application',application, ...
    'info',motion_info);
channel.output_files.motion = motion_output_files(motion_info);
end

function files=motion_output_files(info)
files=strings(0,1);
for name=["shift_file","parameter_file","downsampled_tif_file", ...
        "metrics_file","metrics_fig","metrics_png"]
    if isfield(info,char(name)) && strlength(string(info.(char(name))))>0
        files(end+1,1)=string(info.(char(name))); %#ok<AGROW>
    end
end
files=unique(files,'stable');
end

function profile = conditional_profile(channels, idx)
if isempty(idx)
    profile = struct();
else
    profile = channels(idx).profile;
end
end

function movie_3d = conditional_movie(channels, idx)
if isempty(idx)
    movie_3d = [];
else
    movie_3d = channels(idx).data.movie_3d;
end
end

function save_downsampled_motion_tif( ...
        movie_3d, factor, output_file, q_low, q_high)
movie_ds = downsample_movie_time_mean(movie_3d, factor);
movie_scaled = scale_movie_to_uint16(movie_ds, q_low, q_high);
aaa.io.write_tiff_stack(uint16(movie_scaled), output_file);
end

function movie_ds = downsample_movie_time_mean(movie_3d, factor)
[ncols, nrows, nframes] = size(movie_3d);
nblocks = floor(nframes / factor);
if nblocks < 1
    movie_ds = single(movie_3d);
    return;
end
movie_ds = zeros(ncols, nrows, nblocks, 'single');
for idx = 1:nblocks
    first_frame = (idx - 1) * factor + 1;
    last_frame = idx * factor;
    movie_ds(:, :, idx) = mean( ...
        single(movie_3d(:, :, first_frame:last_frame)), 3, 'omitnan');
end
end

function movie_scaled = scale_movie_to_uint16(movie_3d, q_low, q_high)
low_value = quantile(movie_3d(:), q_low);
high_value = quantile(movie_3d(:), q_high);
if ~isfinite(low_value) || ~isfinite(high_value) ...
        || high_value <= low_value
    movie_scaled = zeros(size(movie_3d), 'single');
    return;
end
movie_scaled = (single(movie_3d) - low_value) ...
    / (high_value - low_value) * 65535;
movie_scaled(movie_scaled < 0) = 0;
movie_scaled(movie_scaled > 65535) = 65535;
end

function save_motion_metrics(raw, corrected, model, paths, factor, label)
if exist('motion_metrics', 'file') ~= 2
    error('AAA:Sections:MissingMotionMetrics', ...
        'NoRMCorre motion_metrics is required for motion QC.');
end
raw_ds = downsample_movie_time_mean(raw, factor);
[correlation_raw, mean_raw, variance_raw] = ...
    motion_metrics(raw_ds, model.options.max_shift);
corrected_ds = downsample_movie_time_mean(corrected, factor);
[correlation_corrected, mean_corrected, variance_corrected] = ...
    motion_metrics(corrected_ds, model.options.max_shift);
shifts_plot = squeeze(cat(3, model.shifts(:).shifts));
metrics_data = struct( ...
    'rigid_shifts', shifts_plot, ...
    'metric_downsample_factor', factor, ...
    'filtered', struct( ...
        'status', 'skipped', ...
        'reason', ['High-pass filtered metrics are skipped to avoid an ' ...
        'extra full-movie apply_shifts and motion_metrics memory peak.']), ...
    'full', struct( ...
        'correlation_raw', correlation_raw, ...
        'correlation_corrected', correlation_corrected, ...
        'mean_raw', mean_raw, ...
        'mean_corrected', mean_corrected, ...
        'variance_raw', variance_raw, ...
        'variance_corrected', variance_corrected));
save(paths.metrics_file, 'metrics_data', '-v7.3');

fig = figure('Color', 'w', 'Name', sprintf('%s Motion Metrics', label));
subplot(3, 1, 1);
plot(shifts_plot, 'LineWidth', 1);
if label == "Shared"
    title('Shared Rigid Shifts');
    metric_channel = "voltage";
else
    title(sprintf('%s Rigid Shifts', label));
    metric_channel = lower(label);
end
legend('y-shifts', 'x-shifts');
grid on;
subplot(3, 1, 2);
text(0.5, 0.5, sprintf( ...
    'High-pass filtered metrics skipped\nFull-%s metrics use time downsample x%d', ...
    metric_channel, factor), ...
    'Units', 'normalized', 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle');
axis off;
subplot(3, 1, 3);
plot(correlation_raw, 'Color', [0.45 0.45 0.45], 'LineWidth', 1);
hold on;
plot(correlation_corrected, 'r', 'LineWidth', 1.2);
if label == "Shared"
    title('Correlation Coefficients On Full Voltage Movie');
else
    title(sprintf('Correlation Coefficients On Full %s Movie', label));
end
legend('raw', 'corrected');
ylim([0.8, 1]);
grid on;
xlabel(sprintf('Downsampled frame (x%d)', factor));
savefig(fig, paths.metrics_fig);
exportgraphics(fig, paths.metrics_png, 'Resolution', 200);
close(fig);
end

function value = spatial_size(movie_3d)
value = [size(movie_3d, 1), size(movie_3d, 2)];
end
