function [runtime, logger, log_file, owns_logger] = ...
    attach_analysis_logger(runtime, request, output_path)
%ATTACH_ANALYSIS_LOGGER Add the standard result-folder logger to runtime.
% Existing caller-supplied loggers are preserved.  OWNS_LOGGER tells the
% workflow whether it should close/delete the returned logger at completion.

if nargin < 1 || isempty(runtime)
    runtime = struct();
end
if nargin < 2 || ~isstruct(request) || ~isscalar(request)
    error('AAA:IO:InvalidLoggerRequest', ...
        'request must be one scalar struct.');
end
if nargin < 3 || strlength(strtrim(string(output_path))) == 0
    error('AAA:IO:MissingLoggerOutputPath', ...
        'A nonempty output path is required for analysis logging.');
end
if ~isstruct(runtime) || ~isscalar(runtime)
    error('AAA:IO:InvalidRuntime', 'runtime must be one scalar struct.');
end

write_log = true;
if isfield(request, 'workflow') && isfield(request.workflow, 'write_analysis_log')
    write_log = logical(request.workflow.write_analysis_log);
end
if ~write_log
    logger = [];
    log_file = "";
    owns_logger = false;
    return;
end

log_file = string(fullfile(string(output_path), 'analysis.log'));
existing_logger = [];
if isfield(runtime, 'logger') && ~isempty(runtime.logger)
    existing_logger = runtime.logger;
end
if isa(existing_logger, 'aaa.io.AnalysisLogger') ...
        && strcmpi(existing_logger.FilePath, log_file)
    logger = existing_logger;
    owns_logger = false;
    return;
end

% A nested Cycle logger still writes its own result-folder log, but the
% parent Record logger owns command-window echo so each event is shown once.
logger = aaa.io.AnalysisLogger(log_file, isempty(existing_logger));
if isempty(existing_logger)
    runtime.logger = logger;
else
    runtime.logger = @(level, message) ...
        write_both(existing_logger, logger, level, message);
end
owns_logger = true;
end

function write_both(first_logger, result_logger, level, message)
write_one(result_logger, level, message);
write_one(first_logger, level, message);
end

function write_one(logger, level, message)
if isa(logger, 'aaa.io.AnalysisLogger')
    logger.write(level, message);
elseif isa(logger, 'function_handle')
    logger(level, message);
else
    error('AAA:IO:InvalidLogger', ...
        'runtime.logger must be an AnalysisLogger or function handle.');
end
end
