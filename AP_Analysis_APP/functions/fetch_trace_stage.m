function data = fetch_trace_stage(results, stage_name)
%FETCH_TRACE_STAGE Return the data stored under one exact trace stage name.

if ~isstruct(results) || ~isfield(results, 'trace_results')
    error('Channel results are not loaded.');
end

if isfield(results.trace_results, stage_name) ...
        && isfield(results.trace_results.(stage_name), 'data')
    data = results.trace_results.(stage_name).data;
    return;
end

error('Trace stage "%s" is not available in the saved channel results.', stage_name);
end
