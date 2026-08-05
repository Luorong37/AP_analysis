function [record,match_kind,cycle] = runtime_preflight(runtime,request)
%RUNTIME_PREFLIGHT Return a matching authoritative App preflight record.

record=struct();match_kind="";cycle=struct();
if ~isstruct(runtime)||~isfield(runtime,'preflight')|| ...
        ~isstruct(runtime.preflight)||~isscalar(runtime.preflight)
    return;
end
candidate=runtime.preflight;
if ~isfield(candidate,'valid')||~logical(candidate.valid),return;end
record=candidate;
key=request_key(request);
if isfield(candidate,'resolved_request_key')&& ...
        string(candidate.resolved_request_key)==key
    match_kind="exact";
    return;
end
if ~isfield(request,'scope')||string(request.scope)~="single"|| ...
        ~isfield(candidate,'inspection')|| ...
        ~isfield(candidate.inspection,'cycles')
    return;
end
cycles=candidate.inspection.cycles;
input_path=normalized_path(string(request.input_path));
for idx=1:numel(cycles)
    if normalized_path(string(cycles(idx).path))==input_path
        cycle=cycles(idx);
        match_kind="cycle";
        return;
    end
end
end

function key=request_key(request)
try,key=string(jsonencode(request));catch,key="";end
end

function value=normalized_path(value)
value=replace(strtrim(string(value)),"\","/");
while strlength(value)>3&&endsWith(value,"/"),value=extractBefore(value,strlength(value));end
value=lower(value);
end
