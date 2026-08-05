function event = emit_event(runtime, event)
%EMIT_EVENT Normalize, log, and deliver one AAA status event.

if nargin < 1 || isempty(runtime)
    runtime = struct();
end
if nargin < 2 || isempty(event)
    event = struct();
end
if ~isstruct(runtime) || ~isscalar(runtime)
    error('AAA:IO:InvalidRuntime', 'runtime must be one scalar struct.');
end
if ~isstruct(event) || ~isscalar(event)
    error('AAA:IO:InvalidEvent', 'event must be one scalar struct.');
end

defaults = struct( ...
    'type', "status", ...
    'mode', "", ...
    'scope', "", ...
    'input_path', "", ...
    'cycle', "", ...
    'section', "", ...
    'state', "", ...
    'message', "", ...
    'output_path', "", ...
    'timestamp', datetime('now'), ...
    'details', struct());
event = merge_struct(defaults, event);
event.type = string(event.type);
event.mode = string(event.mode);
event.scope = string(event.scope);
event.input_path = string(event.input_path);
event.cycle = string(event.cycle);
event.section = string(event.section);
event.state = lower(string(event.state));
event.message = string(event.message);
event.output_path = string(event.output_path);
if isempty(event.timestamp)
    event.timestamp = datetime('now');
end

if isfield(runtime, 'logger') && ~isempty(runtime.logger)
    log_level = "INFO";
    if event.state == "failed"
        log_level = "ERROR";
    elseif event.type == "warning"
        log_level = "WARNING";
    end
    log_message = format_event_message(event);
    try
        if isa(runtime.logger, 'aaa.io.AnalysisLogger')
            runtime.logger.write(log_level, log_message);
        elseif isa(runtime.logger, 'function_handle')
            runtime.logger(log_level, log_message);
        end
    catch ME
        warning('AAA:IO:LoggerCallbackFailed', ...
            'AAA logger callback failed: %s', ME.message);
    end
end

if isfield(runtime, 'notify') && ~isempty(runtime.notify)
    if ~isa(runtime.notify, 'function_handle')
        error('AAA:IO:InvalidNotifyCallback', ...
            'runtime.notify must be a function handle.');
    end
    try
        runtime.notify(event);
    catch ME
        warning('AAA:IO:NotifyCallbackFailed', ...
            'AAA notify callback failed: %s', ME.message);
    end
end
end

function message = format_event_message(event)
parts = strings(0,1);
if strlength(event.mode) > 0, parts(end+1,1) = "mode=" + event.mode; end
if strlength(event.scope) > 0, parts(end+1,1) = "scope=" + event.scope; end
if strlength(event.cycle) > 0, parts(end+1,1) = "cycle=" + event.cycle; end
if strlength(event.section) > 0, parts(end+1,1) = "section=" + event.section; end
if strlength(event.state) > 0, parts(end+1,1) = "state=" + event.state; end
if strlength(event.message) > 0, parts(end+1,1) = event.message; end
message = strjoin(parts, ' | ');
end

function output = merge_struct(base, patch)
output = base;
names = fieldnames(patch);
for idx = 1:numel(names)
    output.(names{idx}) = patch.(names{idx});
end
end
