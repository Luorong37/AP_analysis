function [data, stage_name] = resolve_trace_stage(results, preferred_stages)
%RESOLVE_TRACE_STAGE Select the first available stage in caller priority order.
%
% No fallback is added beyond preferred_stages. The returned stage_name has
% the same scalar text type as the corresponding preferred_stages element.

for idx = 1:numel(preferred_stages)
    stage_name = preferred_stages{idx};
    if has_trace_stage(results, stage_name)
        data = fetch_trace_stage(results, stage_name);
        return;
    end
end

error('None of the requested trace stages are available: %s', ...
    strjoin(preferred_stages, ', '));
end

function tf = has_trace_stage(results, stage_name)
tf = isstruct(results) ...
    && isfield(results, 'trace_results') ...
    && isfield(results.trace_results, stage_name) ...
    && isfield(results.trace_results.(stage_name), 'data');
end
