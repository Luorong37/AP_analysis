function [request, source_name] = direct_run_request(mode, scope)
%DIRECT_RUN_REQUEST Resolve a zero-input Editor Run from the base workspace.
% Explicit function arguments remain authoritative.  This helper is used
% only when a workflow is launched by pressing Run on its file.

mode = lower(strtrim(string(mode)));
scope = lower(strtrim(string(scope)));
request = aaa.schema.default_request(mode, scope);
source_name = "defaults";

if scope == "record"
    mode_name = mode + "3_rec_request";
else
    mode_name = mode + "3_request";
end
% A mode/scope-specific request is more intentional than the generic
% fallback when both variables happen to exist in the base workspace.
candidates = [mode_name, "aaa_request"];
for idx = 1:numel(candidates)
    name = candidates(idx);
    exists = evalin('base', sprintf("exist('%s','var')", name));
    if exists ~= 1
        continue;
    end
    candidate = evalin('base', char(name));
    if ~isstruct(candidate) || ~isscalar(candidate)
        error('AAA:DirectRun:InvalidWorkspaceRequest', ...
            'Base-workspace variable %s must be one scalar request struct.', name);
    end
    request = candidate;
    source_name = name;
    return;
end
end
