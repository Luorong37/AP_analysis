function [ctx, receipt] = run_plan(ctx, plan, entries)
%RUN_PLAN Execute an ordered run/reuse/skip section plan.
% Handlers use the common signature [ctx,outcome] = handler(ctx,action).
% Handler failures stop the plan and are returned in ctx/receipt. The upper
% workflow owns final manifest/log handling and production rethrow behavior.

validate_context(ctx);
if nargin < 2 || ~isstruct(plan) || ~isscalar(plan)
    error('AAA:Sections:InvalidPlan', 'plan must be one scalar struct.');
end
if nargin < 3 || isempty(entries)
    entries = aaa.sections.registry(ctx.mode);
end
validate_registry(entries, ctx.mode);
[ordered_names, plan] = normalize_plan(plan, entries);

ctx = prepare_execution(ctx, plan);
for idx = 1:numel(ordered_names)
    section = ordered_names(idx);
    action = plan.(char(section));
    started_at = datetime('now');
    ctx.execution.current_section = section;
    ctx.execution.current_action = action;
    ctx.execution.current_channel = "";

    if action == "skip"
        message = "Section skipped by the execution plan.";
        record = make_record(section, action, "skipped", ...
            started_at, datetime('now'), message, strings(0,1), "", empty_error());
        ctx.execution.history(end + 1, 1) = record;
        emit_section_event(ctx, section, action, "skipped", message, ...
            strings(0,1), "");
        continue;
    end

    start_state = "running";
    start_message = "Running section.";
    if action == "reuse"
        start_state = "reusing";
        start_message = "Reusing section results.";
    end
    emit_section_event(ctx, section, action, start_state, start_message, ...
        strings(0,1), "");

    entry_index = find([entries.name] == section, 1, 'first');
    handler = entries(entry_index).handler;
    if isempty(handler)
        ME = MException('AAA:Sections:MissingHandler', ...
            'No handler is registered for enabled section "%s".', section);
        ctx = record_failure(ctx, section, action, started_at, ...
            strings(0,1), "", ME);
        break;
    end

    try
        [candidate_ctx, raw_outcome] = handler(ctx, action);
        validate_context(candidate_ctx);
        outcome = normalize_outcome(candidate_ctx, raw_outcome, action);
        if outcome.state == "failed"
            ctx = candidate_ctx;
            ME = outcome_exception(outcome, section);
            ctx = record_failure(ctx, section, action, started_at, ...
                outcome.channels, outcome.failed_channel, ME, outcome.message);
            break;
        end


        [candidate_ctx, result_record] = save_analysis_results( ...
            candidate_ctx, section, outcome);
        outcome.result_record = result_record;
        ctx = candidate_ctx;

        final_state = "completed";
        default_message = "Section completed.";
        if action == "reuse"
            final_state = "reused";
            default_message = "Section results reused.";
        end
        message = outcome.message;
        if strlength(message) == 0
            message = default_message;
        end
        completed_at = datetime('now');
        record = make_record(section, action, final_state, ...
            started_at, completed_at, message, outcome.channels, "", empty_error());
        ctx.execution.history(end + 1, 1) = record;
        ctx = update_channel_status(ctx, outcome.channels, final_state, ...
            section, message, empty_error());
        emit_section_event(ctx, section, action, final_state, message, ...
            outcome.channels, "");
    catch ME
        ctx = record_failure(ctx, section, action, started_at, ...
            strings(0,1), "", ME);
        break;
    end
end

if ctx.execution.status ~= "failed"
    ctx.execution.status = "completed";
    ctx.execution.completed_at = datetime('now');
    ctx.execution.current_section = "";
    ctx.execution.current_action = "";
    ctx.execution.current_channel = "";
end
receipt = aaa.sections.build_receipt(ctx);
end

function validate_context(ctx)
if ~isstruct(ctx) || ~isscalar(ctx)
    error('AAA:Sections:InvalidContext', ...
        'ctx must be one scalar struct created by aaa.sections.create_context.');
end
required = {'mode','scope','request','input_path','output_path', ...
    'runtime','channels','execution'};
missing = setdiff(required, fieldnames(ctx));
if ~isempty(missing)
    error('AAA:Sections:InvalidContext', ...
        'ctx is missing field(s): %s.', strjoin(missing, ', '));
end
if ~ismember(string(ctx.mode), ["ap","dual"])
    error('AAA:Sections:InvalidContext', 'ctx.mode must be ap or dual.');
end
if ~isstruct(ctx.channels)
    error('AAA:Sections:InvalidContext', 'ctx.channels must be a struct array.');
end
for idx = 1:numel(ctx.channels)
    if ~isfield(ctx.channels(idx), 'profile') ...
            || ~isstruct(ctx.channels(idx).profile) ...
            || ~isfield(ctx.channels(idx).profile, 'role')
        error('AAA:Sections:InvalidContext', ...
            'Every channel must contain profile.role.');
    end
end
end

function validate_registry(entries, mode)
if ~isstruct(entries) || isempty(entries) ...
        || ~all(isfield(entries, {'name','mode','handler_name','handler'}))
    error('AAA:Sections:InvalidRegistry', ...
        'entries must be a nonempty registry struct array.');
end
names = string({entries.name});
if any(strlength(names) == 0) || numel(unique(names)) ~= numel(names)
    error('AAA:Sections:InvalidRegistry', ...
        'Registry section names must be nonempty and unique.');
end
if any(string({entries.mode}) ~= string(mode))
    error('AAA:Sections:RegistryModeMismatch', ...
        'Every registry entry must match ctx.mode.');
end
for idx = 1:numel(entries)
    if ~(isempty(entries(idx).handler) ...
            || isa(entries(idx).handler, 'function_handle'))
        error('AAA:Sections:InvalidRegistry', ...
            'Registry handler for %s is invalid.', entries(idx).name);
    end
end
end

function [ordered_names, plan] = normalize_plan(plan, entries)
plan_names = string(fieldnames(plan));
registry_names = [entries.name]';
unknown = setdiff(plan_names, registry_names);
if ~isempty(unknown)
    error('AAA:Sections:UnknownPlanSection', ...
        'Plan contains unregistered section(s): %s.', strjoin(unknown, ', '));
end
for idx = 1:numel(plan_names)
    name = plan_names(idx);
    action = lower(strtrim(string(plan.(char(name)))));
    if ~isscalar(action) || ~ismember(action, ["run","reuse","skip"])
        error('AAA:Sections:InvalidAction', ...
            'Section %s action must be run, reuse, or skip.', name);
    end
    plan.(char(name)) = action;
end
ordered_names = registry_names(ismember(registry_names, plan_names));
end

function ctx = prepare_execution(ctx, plan)
if ~isfield(ctx.execution, 'history') || ~isstruct(ctx.execution.history)
    ctx.execution.history = empty_history();
end
if ~isfield(ctx.execution, 'started_at') || isnat(ctx.execution.started_at)
    ctx.execution.started_at = datetime('now');
end
ctx.execution.status = "running";
ctx.execution.plan = merge_struct(ctx.execution.plan, plan);
ctx.execution.completed_at = NaT;
ctx.execution.failed_section = "";
ctx.execution.failed_channel = "";
ctx.execution.error = empty_error();
ctx.execution.exception = [];
end

function outcome = normalize_outcome(ctx, outcome, action)
if nargin < 2 || isempty(outcome)
    outcome = struct();
end
if ~isstruct(outcome) || ~isscalar(outcome)
    error('AAA:Sections:InvalidOutcome', ...
        'Section handler outcome must be one scalar struct.');
end
defaults = struct( ...
    'state', "", ...
    'message', "", ...
    'channels', strings(0,1), ...
    'failed_channel', "", ...
    'error', empty_error(), ...
    'exception', []);
outcome = merge_struct(defaults, outcome);
outcome.state = lower(strtrim(string(outcome.state)));
outcome.message = string(outcome.message);
outcome.channels = unique(string(outcome.channels(:)), 'stable');
outcome.channels = outcome.channels(strlength(outcome.channels) > 0);
outcome.failed_channel = strtrim(string(outcome.failed_channel));
if ~isscalar(outcome.state) || ~ismember(outcome.state, ...
        ["","completed","reused","failed"])
    error('AAA:Sections:InvalidOutcomeState', ...
        'Handler outcome state is invalid for action %s.', action);
end
if ~isscalar(outcome.message) || ~isscalar(outcome.failed_channel)
    error('AAA:Sections:InvalidOutcome', ...
        'Handler message and failed_channel must be scalar text.');
end
known_roles = channel_roles(ctx);
unknown_roles = setdiff(outcome.channels, known_roles);
if ~isempty(unknown_roles)
    error('AAA:Sections:InvalidOutcomeChannel', ...
        'Handler reported unknown channel(s): %s.', strjoin(unknown_roles, ', '));
end
if strlength(outcome.failed_channel) > 0 ...
        && ~ismember(outcome.failed_channel, known_roles)
    error('AAA:Sections:InvalidOutcomeChannel', ...
        'Handler reported unknown failed channel: %s.', outcome.failed_channel);
end
if outcome.state ~= "failed" && strlength(outcome.failed_channel) > 0
    error('AAA:Sections:InvalidOutcome', ...
        'failed_channel is legal only when outcome.state="failed".');
end
if outcome.state == "failed" && strlength(outcome.failed_channel) > 0 ...
        && ~ismember(outcome.failed_channel, outcome.channels)
    outcome.channels(end + 1, 1) = outcome.failed_channel;
end
end

function roles = channel_roles(ctx)
profiles = [ctx.channels.profile];
roles = unique(string({profiles.role})', 'stable');
end

function ctx = record_failure(ctx, section, action, started_at, ...
        channels, failed_channel, ME, message)
if nargin < 9 || strlength(string(message)) == 0
    message = string(ME.message);
else
    message = string(message);
end
failure = exception_struct(ME);
completed_at = datetime('now');
record = make_record(section, action, "failed", started_at, ...
    completed_at, message, channels, failed_channel, failure);
ctx.execution.history(end + 1, 1) = record;
ctx.execution.status = "failed";
ctx.execution.completed_at = completed_at;
ctx.execution.failed_section = section;
ctx.execution.failed_channel = string(failed_channel);
ctx.execution.error = failure;
ctx.execution.exception = ME;
ctx.execution.current_channel = string(failed_channel);
ctx = update_channel_status(ctx, string(failed_channel), "failed", ...
    section, message, failure);
emit_section_event(ctx, section, action, "failed", message, ...
    channels, failed_channel);
end

function ctx = update_channel_status(ctx, roles, state, section, message, error_value)
roles = string(roles(:));
for role_idx = 1:numel(roles)
    role = roles(role_idx);
    for channel_idx = 1:numel(ctx.channels)
        if string(ctx.channels(channel_idx).profile.role) == role
            ctx.channels(channel_idx).status = struct( ...
                'state', string(state), ...
                'section', string(section), ...
                'message', string(message), ...
                'error', error_value);
            break;
        end
    end
end
end

function emit_section_event(ctx, section, action, state, message, channels, failed_channel)
event = struct( ...
    'type', "section", ...
    'mode', string(ctx.mode), ...
    'scope', string(ctx.scope), ...
    'input_path', string(ctx.input_path), ...
    'output_path', string(ctx.output_path), ...
    'section', string(section), ...
    'state', string(state), ...
    'message', string(message), ...
    'timestamp', datetime('now'), ...
    'details', struct( ...
        'action', string(action), ...
        'channels', string(channels(:)), ...
        'failed_channel', string(failed_channel)));
aaa.io.emit_event(ctx.runtime, event);
end

function record = make_record(section, action, state, started_at, ...
        completed_at, message, channels, failed_channel, error_value)
record = struct( ...
    'section', string(section), ...
    'action', string(action), ...
    'state', string(state), ...
    'started_at', started_at, ...
    'completed_at', completed_at, ...
    'message', string(message), ...
    'channels', string(channels(:)), ...
    'failed_channel', string(failed_channel), ...
    'error', error_value);
end

function ME = outcome_exception(outcome, section)
if isa(outcome.exception, 'MException')
    ME = outcome.exception;
    return;
end
identifier = "AAA:Sections:HandlerReportedFailure";
message = "Section handler reported failure.";
if isstruct(outcome.error) && isscalar(outcome.error)
    if isfield(outcome.error, 'identifier') ...
            && strlength(string(outcome.error.identifier)) > 0
        identifier = string(outcome.error.identifier);
    end
    if isfield(outcome.error, 'message') ...
            && strlength(string(outcome.error.message)) > 0
        message = string(outcome.error.message);
    end
end
if strlength(outcome.message) > 0
    message = outcome.message;
end
ME = MException(char(identifier), '%s', ...
    char("Section " + section + ": " + message));
end

function value = exception_struct(ME)
value = struct( ...
    'identifier', string(ME.identifier), ...
    'message', string(ME.message));
end

function output = merge_struct(base, patch)
output = base;
names = fieldnames(patch);
for idx = 1:numel(names)
    output.(names{idx}) = patch.(names{idx});
end
end

function history = empty_history()
history = repmat(make_record("", "", "", NaT, NaT, "", ...
    strings(0,1), "", empty_error()), 0, 1);
end

function value = empty_error()
value = struct('identifier',"",'message',"");
end
