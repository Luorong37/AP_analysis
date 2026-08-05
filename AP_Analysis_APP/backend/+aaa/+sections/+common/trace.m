function [ctx, outcome] = trace(ctx, action)
%TRACE Run or strictly reuse the shared channel trace pipeline.
% This handler coordinates files and stage provenance. Scientific formulas
% are implemented by ordinary functions in AP_Analysis_APP/functions.

action = normalize_action(action);
outcome = struct( ...
    'section', "trace", ...
    'action', action, ...
    'state', "completed", ...
    'channels', channel_roles(ctx), ...
    'message', "");
if action == "skip"
    outcome.state = "skipped";
    outcome.message = "Trace section skipped.";
    return;
end

validate_context(ctx, action);
require_output_path(ctx);

% Load inputs: reuse saved stages, or obtain current raw traces/profiles.
if action == "reuse"
    aaa.sections.emit_phase(ctx,"trace","reuse_load","running", ...
        "Loading complete saved trace stages without movie pixels.");
    [ctx,source_files,frame_rates,compatibility_report] = reuse_trace_results(ctx);
    outcome.source_files = source_files;
    outcome.saved_frame_rates = frame_rates;
    outcome.compatibility_report = compatibility_report;
    outcome.movie_loaded = false;
    outcome.message = "Reused complete saved trace stages without loading a movie. " + ...
        "Saved frame rates are authoritative: " + format_frame_rates(frame_rates) + ".";
    return;
end

run_background = logical(ctx.request.params.common.run_background_removal);
background_by_role = struct();

% Run algorithm: process each role with the same bound-profile functions.
for channel_idx = 1:numel(ctx.channels)
    channel = ctx.channels(channel_idx);
    role = char(lower(string(channel.profile.role)));
    role_text=string(role);
    roi_file = raw_roi_file(channel.results);
    raw = fetch_trace_stage(channel.results, 'raw');

    parent_stage = 'raw';
    trace_input = raw;
    if run_background
        aaa.sections.emit_phase(ctx,"trace","background","running", ...
            "Removing fitted background.",role_text);
        require_background_inputs(channel);
        [background, background_info] = remove_background( ...
            channel.data.movie_3d, channel.data.roi_mask, channel.profile);
        trace_input = background.roi_minus_background_fitted;
        parent_stage = 'bg_removed';
        channel.results = store_trace_stage( ...
            channel.results, 'bg_removed', trace_input, {'raw'}, roi_file, ...
            channel.results.movie_info, 'remove_background', ...
            channel.profile.roi.background);
        background.info = background_info;
        background_by_role.(role) = background;
    end

    bleach_spec = channel.profile.trace.bleach;
    aaa.sections.emit_phase(ctx,"trace","bleach","running", ...
        "Removing bleach trend and calculating baseline.",role_text);
    bleach_spec.time_axis = (1:size(trace_input, 1))' ...
        / channel.profile.frame_rate;
    [bleach_removed, baseline, bleach_parameters] = ...
        remove_bleach(trace_input, bleach_spec);
    channel.results = store_trace_stage( ...
        channel.results, 'bleach_removed', bleach_removed, ...
        {parent_stage}, roi_file, channel.results.movie_info, ...
        char(string(bleach_spec.mode)), bleach_parameters);
    channel.results = store_trace_stage( ...
        channel.results, 'baseline', baseline, {parent_stage}, roi_file, ...
        channel.results.movie_info, char(string(bleach_spec.mode)), ...
        bleach_parameters);

    aaa.sections.emit_phase(ctx,"trace","metrics","running", ...
        "Computing role-specific trace metrics.",role_text);
    [metrics, metric_info] = compute_trace_metrics( ...
        bleach_removed, baseline, channel.profile);
    channel.results = store_metric_stages( ...
        channel.results, metrics, metric_info, roi_file);

    smoothing = channel.profile.trace.smoothing;
    if lower(string(smoothing.method)) ~= "none"
        aaa.sections.emit_phase(ctx,"trace","smoothing","running", ...
            "Computing configured smoothed stages.",role_text);
        channel.results = store_smoothed_stages( ...
            channel.results, raw, metrics, smoothing, roi_file);
    end

    channel.output_files.trace = strings(0,1);
    ctx.channels(channel_idx) = channel;
end

% Return results: register compact trace provenance for unified saving.
aaa.sections.emit_phase(ctx,"trace","persist","running", ...
    "Preparing the unified trace result.");
ctx.shared.results.trace_details=struct( ...
    'background',background_qc_only(background_by_role));
outcome.message = sprintf('Processed trace stages for %d channel(s).', ...
    numel(ctx.channels));
end

function value=background_qc_only(value)
roles=fieldnames(value);
for idx=1:numel(roles)
    if isfield(value.(roles{idx}),'roi_minus_background_fitted')
        value.(roles{idx})=rmfield( ...
            value.(roles{idx}),'roi_minus_background_fitted');
    end
end
end

function action = normalize_action(action)
action = lower(strtrim(string(action)));
if ~isscalar(action) || ~ismember(action, ["run","reuse","skip"])
    error('AAA:Sections:InvalidTraceAction', ...
        'Trace action must be run, reuse, or skip.');
end
end

function validate_context(ctx, action)
if ~isstruct(ctx) || ~isscalar(ctx) || ~isfield(ctx, 'channels') ...
        || isempty(ctx.channels) || ~isfield(ctx, 'request')
    error('AAA:Sections:InvalidTraceContext', ...
        'Trace context must contain one or more channels.');
end
for idx = 1:numel(ctx.channels)
    channel = ctx.channels(idx);
    incomplete = ~isfield(channel,'profile') || ~isfield(channel.profile,'role') ...
        || ~isfield(channel.profile,'trace') || ~isfield(channel,'results');
    if action == "run"
        incomplete = incomplete || ~isfield(channel.results,'movie_info') ...
            || ~isfield(channel.results,'trace_results');
    end
    if incomplete
        error('AAA:Sections:IncompleteTraceChannel', ...
            'Channel %d is missing profile, movie_info, or trace_results.', idx);
    end
end
end

function roles = channel_roles(ctx)
roles = strings(0,1);
if ~isstruct(ctx) || ~isfield(ctx, 'channels') || isempty(ctx.channels)
    return;
end
profiles = [ctx.channels.profile];
roles = string({profiles.role})';
end

function output_path = require_output_path(ctx)
output_path = string(ctx.output_path);
if strlength(output_path) == 0
    error('AAA:Sections:MissingOutputPath', ...
        'Trace section requires a nonempty output path.');
end
if ~isfolder(output_path)
    mkdir(output_path);
end
end

function require_background_inputs(channel)
if ~isfield(channel, 'data') || ~isfield(channel.data, 'movie_3d') ...
        || isempty(channel.data.movie_3d) ...
        || ~isfield(channel.data, 'roi_mask') ...
        || isempty(channel.data.roi_mask)
    error('AAA:Sections:BackgroundInputMissing', ...
        ['Background removal requires the current movie and ROI mask. ' ...
         'Disable background removal or run input/ROI first.']);
end
end

function roi_file = raw_roi_file(results)
roi_file = "";
if isfield(results.trace_results, 'raw') ...
        && isfield(results.trace_results.raw, 'info') ...
        && isfield(results.trace_results.raw.info, 'roi_file')
    roi_file = string(results.trace_results.raw.info.roi_file);
end
end

function results = store_metric_stages(results, metrics, info, roi_file)
movie_info = results.movie_info;
results = store_trace_stage(results, 'noise_reference', ...
    metrics.noise_reference, {'bleach_removed'}, roi_file, movie_info, ...
    info.noise_reference_method, info.noise_reference_parameters);
results = store_trace_stage(results, 'noise', metrics.noise, ...
    {'bleach_removed','noise_reference'}, roi_file, movie_info, ...
    'residual', struct('expression','bleach_removed - noise_reference'));
results = store_trace_stage(results, 'sensitivity', metrics.sensitivity, ...
    {'bleach_removed','baseline'}, roi_file, movie_info, ...
    'ratio', struct('expression','bleach_removed ./ baseline'));
results = store_trace_stage(results, 'snr', metrics.snr, ...
    {'bleach_removed','noise'}, roi_file, movie_info, ...
    info.snr_method, info.snr_parameters);
end

function results = store_smoothed_stages(results, raw, metrics, smoothing, roi_file)
names = {'raw','sensitivity','snr'};
inputs = {raw, metrics.sensitivity, metrics.snr};
for idx = 1:numel(names)
    [data, smooth_info] = smooth_trace(inputs{idx}, smoothing);
    stage_name = [names{idx}, '_smoothed'];
    results = store_trace_stage(results, stage_name, data, {names{idx}}, ...
        roi_file, results.movie_info, smooth_info.method, ...
        smooth_info.parameters);
end
end

function [ctx,source_files,frame_rates,report] = reuse_trace_results(ctx)
source_path = string(ctx.request.workflow.source_results_path);
if strlength(source_path) == 0 || ~isfolder(source_path)
    error('AAA:Sections:TraceReuseSourceMissing', ...
        'trace="reuse" requires workflow.source_results_path.');
end
if ~aaa.io.manifest_allows_reuse(source_path)
    error('AAA:Sections:TraceReuseManifestRejected', ...
        'Source manifest is not completed and cannot be reused: %s',source_path);
end
spec=struct('mode',ctx.mode,'scope',ctx.scope,'runtime',ctx.runtime, ...
    'max_mat_files',discovery_max_files(ctx.request));
[loaded_results,report]=load_analysis_results(source_path,spec);
source_files=report.absolute_source_files;
frame_rates=struct();
if isfield(loaded_results,'roi'),ctx.shared.results.roi=loaded_results.roi;end
if isfield(loaded_results,'maps'),ctx.shared.results.maps=loaded_results.maps;end
if isfield(loaded_results,'registration')
    ctx.shared.results.registration=loaded_results.registration;
end
if isfield(loaded_results,'trace_details')
    ctx.shared.results.trace_details=loaded_results.trace_details;
end
for channel_idx = 1:numel(ctx.channels)
    role = lower(string(ctx.channels(channel_idx).profile.role));
    if ~isfield(loaded_results,char(role))
        error('AAA:Sections:TraceReuseRoleMissing', ...
            'No compatible saved trace was found for role %s.',role);
    end
    loaded=loaded_results.(char(role));
    validate_complete_trace_results(loaded, ctx.channels(channel_idx).profile);
    saved_frame_rate = saved_frame_rate_from_results(loaded,role);
    ctx.channels(channel_idx).profile.frame_rate = saved_frame_rate;
    if isfield(ctx.channels(channel_idx).profile,'trace') ...
            && isfield(ctx.channels(channel_idx).profile.trace,'bleach')
        ctx.channels(channel_idx).profile.trace.bleach.frame_rate = saved_frame_rate;
    end
    [nframes,movie_info] = reused_trace_metadata(loaded,saved_frame_rate, ...
        ctx.channels(channel_idx).profile.role,source_path);
    loaded.movie_info = movie_info;
    ctx.channels(channel_idx).results = loaded;
    ctx.channels(channel_idx).data.movie_info = movie_info;
    ctx.channels(channel_idx).data.time = ...
        (1:nframes)' / double(ctx.channels(channel_idx).profile.frame_rate);
    ctx.channels(channel_idx).output_files.trace = strings(0,1);
    frame_rates.(char(role))=saved_frame_rate;
end
ctx.input_geometry=struct('action',"not_required",'common_size',[]);
ctx.input_layout=struct('action',"reuse_trace", ...
    'source_results_path',source_path,'source_files',source_files, ...
    'saved_frame_rates',frame_rates,'movie_loaded',false);
end

function frame_rate=saved_frame_rate_from_results(results,role)
values=[];
if isfield(results,'movie_info')&&isfield(results.movie_info,'frame_rate')
    values(end+1)=double(results.movie_info.frame_rate); %#ok<AGROW>
end
if isfield(results,'trace_results')&&isfield(results.trace_results,'raw') ...
        &&isfield(results.trace_results.raw,'frame_rate')
    values(end+1)=double(results.trace_results.raw.frame_rate); %#ok<AGROW>
end
values=values(isfinite(values)&values>0);
if isempty(values)
    error('AAA:Sections:TraceReuseFrameRateMissing', ...
        'Saved trace role %s has no authoritative frame rate.',role);
end
if max(values)-min(values)>max(1e-9,1e-9*max(values))
    error('AAA:Sections:TraceReuseFrameRateMismatch', ...
        'Saved frame rates disagree for role %s.',role);
end
frame_rate=values(1);
end

function [nframes,movie_info] = reused_trace_metadata(results,frame_rate,role,source_path)
raw=fetch_trace_stage(results,'raw');
nframes=size(raw,1);
if nframes<1
    error('AAA:Sections:TraceReuseEmpty','Saved raw trace has no frames.');
end
if isfield(results,'movie_info') && isstruct(results.movie_info)
    movie_info=results.movie_info;
else
    movie_info=struct();
end
movie_info.frame_count=nframes;
movie_info.frame_rate=double(frame_rate);
movie_info.role=string(role);
movie_info.reuse_source_path=string(source_path);
movie_info.movie_loaded=false;
movie_info=enrich_saved_source_metadata(movie_info,role);
end

function info=enrich_saved_source_metadata(info,role)
if isfield(info,'source_files') && ~isempty(info.source_files),return;end
if ~isfield(info,'source_path') || strlength(string(info.source_path))==0,return;end
spec=struct('camera_index',1,'role',string(role));
try
    [~,report]=aaa.helpers.common.resolve_movie_sources(info.source_path,spec);
    if ~report.valid||isempty(report.sources),return;end
    source=report.sources(1);
    info.logical_source_path=source.logical_path;
    info.source_path=source.path;
    info.source_files=source.files;
    info.source_part_count=source.part_count;
    info.source_is_split=source.is_split;
    info.source_resolution_method="reuse_saved_source_"+source.resolution_method;
catch
    % Reuse remains valid when the historical movie has moved or is offline.
end
end

function value=discovery_max_files(request)
value=500;
if isfield(request,'workflow') ...
        && isfield(request.workflow,'discovery_max_mat_files')
    value=request.workflow.discovery_max_mat_files;
end
end

function message=format_frame_rates(values)
roles=string(fieldnames(values)); parts=strings(1,numel(roles));
for idx=1:numel(roles)
    parts(idx)=roles(idx)+"="+string(values.(char(roles(idx))))+" Hz";
end
message=strjoin(parts,", ");
end

function validate_complete_trace_results(results, profile)
required = {'raw','bleach_removed','baseline','noise_reference', ...
    'noise','sensitivity','snr'};
if lower(string(profile.trace.smoothing.method)) ~= "none"
    required = [required, {'raw_smoothed','sensitivity_smoothed','snr_smoothed'}];
end
for idx = 1:numel(required)
    fetch_trace_stage(results, required{idx});
end
nframes = size(fetch_trace_stage(results, 'raw'), 1);
nrois = size(fetch_trace_stage(results, 'raw'), 2);
for idx = 2:numel(required)
    data = fetch_trace_stage(results, required{idx});
    if ~isequal(size(data), [nframes, nrois])
        error('AAA:Sections:TraceReuseShapeMismatch', ...
            'Saved trace stage %s does not match raw trace shape.', required{idx});
    end
end
end
