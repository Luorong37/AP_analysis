function [request, resolution] = resolve_request_reuse_source(request, target)
%RESOLVE_REQUEST_REUSE_SOURCE Resolve one single-Cycle request before validation.

if nargin < 2 || isempty(target), target = infer_target(request); end
resolution = struct('attempted',false,'source_path',"",'report',struct(), ...
    'frame_rate_overrides',repmat(struct('role',"",'ui_value',NaN, ...
    'saved_value',NaN),0,1));
if ~isstruct(request) || ~isscalar(request) ...
        || ~isfield(request,'workflow') || ~isstruct(request.workflow)
    return;
end
mode = lower(string(request.mode));
preset = string(request.workflow.preset);
sections = struct();
if isfield(request.workflow,'sections') && isstruct(request.workflow.sections)
    sections = request.workflow.sections;
end
try
    [~,~,~,plan] = aaa.schema.resolve_sections(mode, preset, sections);
catch
    return;
end
if ~plan_requires_source(plan, request.workflow)
    return;
end
auto_detect = true;
if isfield(request.workflow,'auto_detect_reuse_source')
    auto_detect = logical(request.workflow.auto_detect_reuse_source);
end
explicit = strtrim(string(request.workflow.source_results_path));
if strlength(explicit) > 0
    search_root = explicit;
elseif ~auto_detect
    return;
else
    search_root = string(request.input_path);
end
resolution.attempted = true;
[source_path, report] = aaa.helpers.common.resolve_reuse_source( ...
    mode, search_root, plan, target, request.workflow);
resolution.source_path = source_path;
resolution.report = report;
if strlength(source_path) > 0
    request.workflow.source_results_path = source_path;
    [request,resolution.frame_rate_overrides] = ...
        apply_saved_frame_rates(request,report);
else
    request.workflow.source_results_path = "";
end

function [request,overrides] = apply_saved_frame_rates(request,report)
overrides=repmat(struct('role',"",'ui_value',NaN,'saved_value',NaN),0,1);
if ~isfield(report,'chosen_inspection') ...
        || ~isstruct(report.chosen_inspection) ...
        || ~isfield(report.chosen_inspection,'saved_frame_rates')
    return;
end
rates=report.chosen_inspection.saved_frame_rates;
roles=string(fieldnames(rates));
for idx=1:numel(roles)
    role=roles(idx); saved=double(rates.(char(role)));
    if ~isscalar(saved)||~isfinite(saved)||saved<=0,continue;end
    if string(request.mode)=="ap"
        ui=double(request.params.ap.freq);
        request.params.ap.freq=saved;
    else
        key=char(role+"_frame_rate");
        if ~isfield(request.params.dual,key),continue;end
        ui=double(request.params.dual.(key));
        request.params.dual.(key)=saved;
    end
    overrides(end+1,1)=struct('role',role,'ui_value',ui, ...
        'saved_value',saved); %#ok<AGROW>
end
end
end

function tf = plan_requires_source(plan, workflow)
source_plan = plan;
roi_path = "";
if isfield(workflow,'roi_path'), roi_path = string(workflow.roi_path); end
if isfield(source_plan,'roi') && source_plan.roi == "reuse" && strlength(roi_path) > 0
    source_plan.roi = "skip";
end
tf = any(structfun(@(value) string(value) == "reuse", source_plan));
end

function target = infer_target(request)
target = struct('cycle_name',"",'label',"");
if ~isfield(request,'input_path'), return; end
[~,name] = fileparts(char(string(request.input_path)));
if startsWith(string(name), "Cycle", 'IgnoreCase', true)
    target.cycle_name = string(name);
    target.label = string(name);
end
end
