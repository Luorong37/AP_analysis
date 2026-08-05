function entries = registry(mode, overrides)
%REGISTRY Return the ordered section-to-handler mapping for one mode.
% Missing handlers remain explicit empty entries during staged extraction.
% A run/reuse action against an empty entry fails in run_plan; skip never
% invokes a handler. Tests and callers may replace handlers with overrides.

if nargin < 1 || isempty(mode)
    mode = "dual";
end
if nargin < 2 || isempty(overrides)
    overrides = struct();
end
mode = lower(strtrim(string(mode)));
if ~isscalar(mode) || ~ismember(mode, ["ap", "dual"])
    error('AAA:Sections:InvalidMode', ...
        'mode must be "ap" or "dual".');
end
if ~isstruct(overrides) || ~isscalar(overrides)
    error('AAA:Sections:InvalidRegistryOverrides', ...
        'Registry overrides must be one scalar struct.');
end

catalog = aaa.schema.section_catalog(mode);
entries = repmat(entry_template(mode), numel(catalog), 1);
for idx = 1:numel(catalog)
    name = catalog(idx).name;
    qualified_name = handler_name(mode, name);
    entries(idx).name = name;
    entries(idx).handler_name = qualified_name;
    entries(idx).handler = resolve_handler(qualified_name);
end

override_names = fieldnames(overrides);
known_names = cellstr([entries.name]);
unknown_names = setdiff(override_names, known_names);
if ~isempty(unknown_names)
    error('AAA:Sections:UnknownRegistryOverride', ...
        'Unknown %s section override(s): %s.', ...
        mode, strjoin(unknown_names, ', '));
end
for idx = 1:numel(override_names)
    name = override_names{idx};
    value = overrides.(name);
    if ~(isempty(value) || isa(value, 'function_handle'))
        error('AAA:Sections:InvalidHandlerOverride', ...
            'Handler override for %s must be a function handle or empty.', name);
    end
    entry_index = find([entries.name] == string(name), 1, 'first');
    entries(entry_index).handler = value;
    if isa(value, 'function_handle')
        entries(entry_index).handler_name = string(func2str(value));
    else
        entries(entry_index).handler_name = "";
    end
end
end

function entry = entry_template(mode)
entry = struct( ...
    'name', "", ...
    'mode', mode, ...
    'handler_name', "", ...
    'handler', []);
end

function qualified_name = handler_name(mode, name)
common_names = ["input","motion","roi","trace","peak", ...
    "stim","visualization","export"];
if ismember(name, common_names)
    qualified_name = "aaa.sections.common." + name;
else
    qualified_name = "aaa.sections." + mode + "." + name;
end
end

function handler = resolve_handler(qualified_name)
handler = [];
if strlength(qualified_name) == 0
    return;
end
if ~isempty(which(char(qualified_name)))
    handler = str2func(char(qualified_name));
end
end
