function [request, report] = validate_request(request)
%VALIDATE_REQUEST Apply defaults and validate an AAA request without I/O.
% The function returns errors instead of throwing so the App can display
% all actionable problems before enabling Run.

report = struct('valid',false,'errors',strings(0,1),'warnings',strings(0,1), ...
    'execution_plan',struct(),'canonical_preset',"",'resolution_info',struct());

if nargin < 1 || ~isstruct(request) || ~isscalar(request)
    request = struct();
    report.errors(end+1,1) = "Request must be one scalar struct.";
    return;
end

[mode, mode_error] = normalize_choice_field(request, 'mode', ["ap","dual"]);
[scope, scope_error] = normalize_choice_field(request, 'scope', ["single","record"]);
if strlength(mode_error) > 0
    report.errors(end+1,1) = mode_error;
end
if strlength(scope_error) > 0
    report.errors(end+1,1) = scope_error;
end
if ~isempty(report.errors)
    return;
end

defaults = aaa.schema.default_request(mode, scope);
request = merge_struct(defaults, request);
request.mode = mode;
request.scope = scope;

unknown_errors = validate_unknown_fields(request, mode, scope);
report.errors = [report.errors; unknown_errors];

catalog = aaa.schema.parameter_catalog(mode, scope);
for idx = 1:numel(catalog)
    entry = catalog(idx);
    [exists, value] = get_path(request, entry.key);
    if ~exists
        report.errors(end+1,1) = "Missing parameter: " + entry.key;
        continue;
    end
    [ok, normalized, message] = validate_value(value, entry);
    if ok
        if isfield(entry, 'editable') && ~entry.editable ...
                && ~isequaln(normalized, entry.default)
            report.errors(end+1,1) = entry.key + ...
                ": is read-only in the current backend and must remain " + ...
                value_summary(entry.default) + "."; %#ok<AGROW>
        else
            request = set_path(request, entry.key, normalized);
        end
    else
        report.errors(end+1,1) = entry.key + ": " + message;
    end
end

try
    [canonical, ~, overlay, plan, info] = aaa.schema.resolve_sections( ...
        mode, request.workflow.preset, request.workflow.sections);
    request.workflow.preset = canonical;
    request.workflow.sections = overlay;
    report.execution_plan = plan;
    report.canonical_preset = canonical;
    report.resolution_info = info;
catch ME
    report.errors(end+1,1) = string(ME.message);
end

input_path = strtrim(string(request.input_path));
if ~isscalar(input_path) || strlength(input_path) == 0
    report.errors(end+1,1) = "input_path is required.";
elseif scope == "record" && ~isfolder(input_path)
    report.errors(end+1,1) = "Record input_path must be an existing folder: " + input_path;
elseif scope == "single" && ~(isfolder(input_path) || isfile(input_path))
    report.errors(end+1,1) = "Single input_path does not exist: " + input_path;
end

if ~isempty(fieldnames(report.execution_plan))
    [request, report] = validate_reuse_paths(request, report);
end
if scope == "record"
    [request, report] = validate_record_controls(request, report);
end
if request.params.common.background.outer_distance ...
        <= request.params.common.background.inner_distance
    report.errors(end+1,1) = ...
        "Background outer_distance must be greater than inner_distance.";
end
if mode == "dual"
    if string(request.params.dual.camera1_role) == string(request.params.dual.camera2_role)
        report.errors(end+1,1) = "Dual camera1_role and camera2_role must be different.";
    end
    if request.params.common.motion.display_quantile_low >= request.params.common.motion.display_quantile_high
        report.errors(end+1,1) = "Motion display low quantile must be lower than the high quantile.";
    end
end
report.errors = unique(report.errors, 'stable');
report.warnings = unique(report.warnings, 'stable');
report.valid = isempty(report.errors);
end

function text = value_summary(value)
if isempty(value)
    text = "[]";
elseif isstring(value) && isscalar(value)
    text = "\"" + value + "\"";
elseif ischar(value)
    text = "\"" + string(value) + "\"";
elseif isnumeric(value) || islogical(value)
    text = string(mat2str(value));
else
    text = "the catalog default";
end
end

function [request, report] = validate_reuse_paths(request, report)
plan = report.execution_plan;
roi_path = strtrim(string(request.workflow.roi_path));
source_path = strtrim(string(request.workflow.source_results_path));
if strlength(roi_path) > 0 && plan.roi ~= "reuse"
    report.errors(end+1,1) = ...
        "workflow.roi_path is set, but the ROI section is not reuse.";
end
if strlength(roi_path) > 0 && ~(isfolder(roi_path) || isfile(roi_path))
    report.errors(end+1,1) = "ROI source does not exist: " + roi_path;
end

source_plan = plan;
if isfield(source_plan, 'roi') && source_plan.roi == "reuse" && strlength(roi_path) > 0
    source_plan.roi = "skip";
end
needs_source = any(structfun(@(value) string(value) == "reuse", source_plan));
if needs_source && strlength(source_path) == 0
    if request.scope == "record"
        report.warnings(end+1,1) = ...
            "Record reuse source is unresolved; the Record scheduler must resolve and validate it per Cycle.";
    else
        report.errors(end+1,1) = ...
            "No compatible reuse source was resolved. Set workflow.source_results_path or enable auto-detection under a Cycle/Record input path.";
    end
elseif strlength(source_path) > 0 && ~isfolder(source_path)
    report.errors(end+1,1) = "Source result folder does not exist: " + source_path;
end
if plan.roi == "reuse" && strlength(roi_path) == 0 && strlength(source_path) == 0
    if request.scope == "record"
        report.warnings(end+1,1) = ...
            "Record ROI source is unresolved; the Record scheduler must resolve a shared/per-Cycle ROI.";
    else
        report.errors(end+1,1) = "roi=reuse requires workflow.roi_path or workflow.source_results_path.";
    end
end
request.workflow.roi_path = roi_path;
request.workflow.source_results_path = source_path;
end

function [request, report] = validate_record_controls(request, report)
record = request.record;
if record.reuse_reference_roi && record.redraw_reference_roi
    report.errors(end+1,1) = ...
        "reuse_reference_roi and redraw_reference_roi cannot both be selected.";
end
if record.record_average_only && ~record.record_average
    report.errors(end+1,1) = ...
        "record_average_only requires the Record average checkbox to be selected.";
end
if record.reference_roi_only && strlength(string(record.reference_cycle_name)) == 0
    report.errors(end+1,1) = "reference_roi_only requires a reference_cycle_name.";
end
if request.mode == "dual" && isfield(record, 'run_motion') ...
        && isfield(report.execution_plan, 'motion')
    plan_runs_motion = ismember(report.execution_plan.motion, ["run","reuse"]);
    if logical(record.run_motion) ~= plan_runs_motion
        report.warnings(end+1,1) = ...
            "record.run_motion differs from the resolved section table; the section table is authoritative.";
    end
end
end

function errors = validate_unknown_fields(request, mode, scope)
errors = strings(0,1);
allowed_top = {'schema_version','mode','scope','input_path','workflow','params','record'};
errors = append_unknown(errors, request, allowed_top, 'request');

catalog = aaa.schema.parameter_catalog(mode, scope);
allowed_paths = [catalog.key];
leaf_paths = list_leaf_paths(request, "");
    fixed_paths = ["schema_version","mode","scope","workflow.sections", ...
        "params.common","params.ap","params.dual","record"];
section_names = [aaa.schema.section_catalog(mode).name];
for idx = 1:numel(leaf_paths)
    path = leaf_paths(idx);
    if startsWith(path, "workflow.sections.")
        section_name = extractAfter(path, "workflow.sections.");
        if ~ismember(section_name, section_names)
            errors(end+1,1) = "Unknown section field: " + section_name; %#ok<AGROW>
        end
    elseif ~ismember(path, allowed_paths) && ~ismember(path, fixed_paths)
        errors(end+1,1) = "Unknown request field: " + path; %#ok<AGROW>
    end
end
end

function errors = append_unknown(errors, value, allowed, label)
names = fieldnames(value);
unknown = setdiff(names, allowed);
for idx = 1:numel(unknown)
    errors(end+1,1) = "Unknown " + string(label) + " field: " + string(unknown{idx}); %#ok<AGROW>
end
end

function paths = list_leaf_paths(value, prefix)
paths = strings(0,1);
if ~isstruct(value) || ~isscalar(value) || isempty(fieldnames(value))
    if strlength(prefix) > 0
        paths(end+1,1) = prefix;
    end
    return;
end
names = fieldnames(value);
for idx = 1:numel(names)
    name = string(names{idx});
    child_path = name;
    if strlength(prefix) > 0
        child_path = prefix + "." + name;
    end
    child = value.(names{idx});
    if isstruct(child) && isscalar(child) && ~isempty(fieldnames(child))
        paths = [paths; list_leaf_paths(child, child_path)]; %#ok<AGROW>
    else
        paths(end+1,1) = child_path; %#ok<AGROW>
    end
end
end

function [ok, normalized, message] = validate_value(value, entry)
ok = true;
normalized = value;
message = "";
type = entry.value_type;
if isempty(value) && isempty(entry.default)
    return;
end
switch type
    case {"string","path"}
        normalized = string(value);
        ok = isscalar(normalized);
        message = "must be one text value.";
    case "string_vector"
        try
            normalized = string(value(:));
        catch
            ok = false;
        end
        message = "must be a text vector.";
    case "logical"
        ok = (islogical(value) || isnumeric(value)) && isscalar(value) && ...
            (islogical(value) || ismember(double(value), [0 1]));
        if ok, normalized = logical(value); end
        message = "must be one checkbox/logical value.";
    case "double"
        ok = isnumeric(value) && isreal(value) && isscalar(value);
        if ok, normalized = double(value); end
        message = "must be one numeric value.";
    case "integer"
        ok = isnumeric(value) && isreal(value) && isscalar(value) && ...
            isfinite(value) && value == round(value);
        if ok, normalized = double(value); end
        message = "must be one finite integer.";
    case "numeric_vector"
        ok = isnumeric(value) && isreal(value) && (isempty(value) || isvector(value));
        if ok, normalized = double(value); end
        message = "must be a numeric vector.";
    case "enum"
        choices = entry.choices;
        if isnumeric(value) && isscalar(value)
            normalized = double(value);
            ok = ismember(string(value), choices);
        else
            normalized = lower(strtrim(string(value)));
            ok = isscalar(normalized) && any(strcmpi(normalized, choices));
            if ok
                matched = find(strcmpi(normalized, choices), 1, 'first');
                normalized = choices(matched);
            end
        end
        message = "must be one of: " + strjoin(choices, ', ') + ".";
    case "struct"
        ok = isstruct(value);
        message = "must be a struct.";
    otherwise
        ok = false;
        message = "has unsupported catalog type " + type + ".";
end

if ok && isnumeric(normalized) && ~isempty(normalized)
    finite_values = normalized(isfinite(normalized));
    if ~isempty(entry.minimum) && any(finite_values < entry.minimum)
        ok = false;
        message = "must be >= " + string(entry.minimum) + ".";
    elseif ~isempty(entry.maximum) && any(finite_values > entry.maximum)
        ok = false;
        message = "must be <= " + string(entry.maximum) + ".";
    end
end
end

function [value, message] = normalize_choice_field(request, name, allowed)
message = "";
if ~isfield(request, name)
    value = "";
    message = "Missing request." + name + ".";
    return;
end
value = lower(strtrim(string(request.(name))));
if ~isscalar(value) || ~ismember(value, allowed)
    message = "request." + name + " must be one of: " + strjoin(allowed, ', ') + ".";
end
end

function output = merge_struct(base, override)
output = base;
names = fieldnames(override);
for idx = 1:numel(names)
    name = names{idx};
    if isfield(output, name) && isstruct(output.(name)) && isscalar(output.(name)) ...
            && isstruct(override.(name)) && isscalar(override.(name))
        output.(name) = merge_struct(output.(name), override.(name));
    else
        output.(name) = override.(name);
    end
end
end

function [exists, value] = get_path(source, path)
parts = cellstr(split(string(path), '.'));
value = source;
exists = true;
for idx = 1:numel(parts)
    if ~isstruct(value) || ~isscalar(value) || ~isfield(value, parts{idx})
        exists = false;
        value = [];
        return;
    end
    value = value.(parts{idx});
end
end

function target = set_path(target, path, value)
parts = cellstr(split(string(path), '.'));
target = set_recursive(target, parts, value);
end

function target = set_recursive(target, parts, value)
name = parts{1};
if isscalar(parts)
    target.(name) = value;
else
    if ~isfield(target, name) || ~isstruct(target.(name)) || ~isscalar(target.(name))
        target.(name) = struct();
    end
    target.(name) = set_recursive(target.(name), parts(2:end), value);
end
end
