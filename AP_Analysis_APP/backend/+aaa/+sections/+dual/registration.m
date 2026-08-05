function [ctx, outcome] = registration(ctx, action)
%REGISTRATION Coordinate dual-channel registration and persistence.
% register_channels.m owns the scientific/interactive registration unit;
% this handler owns request binding, reuse, and result files.

action = normalize_action(action);
outcome = struct('section',"registration",'action',action, ...
    'state',"completed",'channels',["voltage";"calcium"],'message',"");
if action == "skip"
    outcome.state = "skipped";
    outcome.message = "Channel registration skipped.";
    return;
end
validate_context(ctx);
output_path = require_output_path(ctx);

% Load inputs: obtain current movies or a reusable registration result.
if action == "reuse"
    aaa.sections.emit_phase(ctx,"registration","reuse_load","running", ...
        "Loading saved channel registration.");
    [info, source_file, compatibility_report] = load_registration(ctx);
    preview = struct();
    message = "Reused saved dual-channel registration.";
else
    % Run algorithm: register current movies using bound role profiles.
    voltage_idx = role_index(ctx, "voltage");
    calcium_idx = role_index(ctx, "calcium");
    spec = struct( ...
        'mode',string(ctx.request.workflow.offset_mode), ...
        'reuse_offset',ctx.request.workflow.reuse_offset);
    phase_state="running";
    if contains(lower(string(spec.mode)),"manual"),phase_state="waiting";end
    aaa.sections.emit_phase(ctx,"registration","register_channels",phase_state, ...
        "Computing channel registration; manual mode may wait for user input.");
    [info, preview] = register_channels( ...
        ctx.channels(voltage_idx).data.movie_3d, ...
        ctx.channels(calcium_idx).data.movie_3d, ...
        ctx.channels(voltage_idx).profile, ...
        ctx.channels(calcium_idx).profile, spec);
    info.preview_files=render_registration_preview(preview,info.mode,output_path);
    source_file = "";
    compatibility_report = struct();
    message = "Computed dual-channel coordinate registration.";
end

% Return results: register offsets, preview references, and channel metadata.
info.source_file = string(source_file);
ctx.shared.results.registration = info;
ctx.shared.data.registration_preview = preview;
for idx = 1:numel(ctx.channels)
    if ~isfield(ctx.channels(idx).results, 'movie_info')
        ctx.channels(idx).results.movie_info = struct();
    end
    ctx.channels(idx).results.movie_info.channel_registration = info;
    ctx.channels(idx).output_files.registration = strings(0,1);
end
aaa.sections.emit_phase(ctx,"registration","persist","running", ...
    "Preparing the unified registration result.");
if isfield(info,'preview_files')
    ctx.shared.output_files.registration=info.preview_files;
else
    ctx.shared.output_files.registration=strings(0,1);
end
outcome.source_files=string(source_file);
outcome.compatibility_report=compatibility_report;
outcome.message = message;
end

function action = normalize_action(action)
action = lower(strtrim(string(action)));
if ~isscalar(action) || ~ismember(action,["run","reuse","skip"])
    error('AAA:Sections:InvalidRegistrationAction', ...
        'Registration action must be run, reuse, or skip.');
end
end

function validate_context(ctx)
if ~isstruct(ctx) || ~isscalar(ctx) || string(ctx.mode) ~= "dual" ...
        || numel(ctx.channels) ~= 2
    error('AAA:Sections:InvalidRegistrationContext', ...
        'Registration requires one dual context with two channels.');
end
for idx = 1:numel(ctx.channels)
    channel = ctx.channels(idx);
    if ~isfield(channel.profile,'role') || ~isfield(channel.data,'movie_3d') ...
            || isempty(channel.data.movie_3d)
        error('AAA:Sections:RegistrationMovieMissing', ...
            'Registration requires loaded voltage and calcium movies.');
    end
end
end

function [info, source_file, report] = load_registration(ctx)
source = string(ctx.request.workflow.source_results_path);
if strlength(source) == 0
    error('AAA:Sections:RegistrationReuseSourceMissing', ...
        'registration="reuse" requires workflow.source_results_path.');
end
[results,report]=load_analysis_results(source,struct( ...
    'mode',ctx.mode,'scope',ctx.scope,'runtime',ctx.runtime));
if isfield(results,'registration')
    info=results.registration;
elseif isfield(results,'voltage')&&isfield(results.voltage,'movie_info') ...
        &&isfield(results.voltage.movie_info,'channel_registration')
    info=results.voltage.movie_info.channel_registration;
else
    error('AAA:Sections:RegistrationReuseMissing', ...
        'No complete registration result was found at %s.',source);
end
source_file="";if ~isempty(report.absolute_source_files),source_file=report.absolute_source_files(1);end
validate_registration(info,source_file);
end

function validate_registration(info, source)
if ~isstruct(info) || ~isscalar(info) || ~isfield(info,'offset_xy') ...
        || ~isnumeric(info.offset_xy) || numel(info.offset_xy) ~= 2 ...
        || any(~isfinite(info.offset_xy))
    error('AAA:Sections:InvalidRegistrationReuse', ...
        'Registration result %s lacks a finite offset_xy.', source);
end
end

function output_path = require_output_path(ctx)
output_path = string(ctx.output_path);
if strlength(output_path) == 0
    error('AAA:Sections:MissingOutputPath', ...
        'Registration requires a nonempty output path.');
end
if ~isfolder(output_path), mkdir(output_path); end
end

function idx = role_index(ctx, role)
profiles = [ctx.channels.profile];
idx = find(string({profiles.role}) == string(role),1);
if isempty(idx)
    error('AAA:Sections:MissingRegistrationRole', ...
        'Dual context is missing role %s.', role);
end
end
