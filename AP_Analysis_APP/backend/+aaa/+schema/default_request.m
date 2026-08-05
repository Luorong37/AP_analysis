function request = default_request(mode, scope)
%DEFAULT_REQUEST Build one complete App-visible AAA request.

if nargin < 1 || isempty(mode)
    mode = "dual";
end
if nargin < 2 || isempty(scope)
    scope = "single";
end
mode = normalize_choice(mode, ["ap","dual"], 'mode');
scope = normalize_choice(scope, ["single","record"], 'scope');

request = struct( ...
    'schema_version', "1.0.0", ...
    'mode', mode, ...
    'scope', scope, ...
    'input_path', "", ...
    'workflow', struct('preset',"full",'sections',struct(), ...
        'source_results_path',"",'roi_path',"",'output_path',"",'run_name',""), ...
    'params', struct('common',struct(),'ap',struct(),'dual',struct()), ...
    'record', struct());

catalog = aaa.schema.parameter_catalog(mode, scope);
for idx = 1:numel(catalog)
    request = set_path(request, catalog(idx).key, catalog(idx).default);
end

% Sections remains an overlay. The complete plan is always produced by
% resolve_sections so a preset change never leaves stale actions behind.
request.workflow.sections = struct();
end

function value = normalize_choice(value, allowed, label)
value = lower(strtrim(string(value)));
if ~isscalar(value) || ~ismember(value, allowed)
    error('AAA:Schema:InvalidChoice', '%s must be one of: %s.', ...
        label, strjoin(allowed, ', '));
end
end

function target = set_path(target, path, value)
parts = split(string(path), '.');
target = set_recursive(target, cellstr(parts), value);
end

function target = set_recursive(target, parts, value)
name = parts{1};
if isscalar(parts)
    target.(name) = value;
    return;
end
if ~isfield(target, name) || ~isstruct(target.(name)) || ~isscalar(target.(name))
    target.(name) = struct();
end
target.(name) = set_recursive(target.(name), parts(2:end), value);
end
