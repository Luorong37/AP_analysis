function receipt = run_analysis3(request, runtime)
%RUN_ANALYSIS3 Validate and dispatch one AP/Dual single/record request.

if nargin < 1 || isempty(request)
    request = aaa.schema.default_request("dual", "single");
end
if nargin < 2 || isempty(runtime)
    runtime = struct();
end
[preflight,preflight_match]=aaa.helpers.common.runtime_preflight(runtime,request);
if preflight_match=="exact"&&isfield(preflight,'resolved_request')
    request=preflight.resolved_request;
elseif isstruct(request) && isfield(request,'scope') ...
        && lower(string(request.scope)) == "single"
    [request,~] = aaa.helpers.common.resolve_request_reuse_source(request);
end

[request, validation] = aaa.schema.validate_request(request);
if ~validation.valid
    error('AAA:InvalidRequest', 'AAA request validation failed:\n%s', ...
        strjoin(validation.errors, newline));
end

aaa.io.emit_event(runtime, struct( ...
    'type',"workflow",'mode',request.mode,'scope',request.scope, ...
    'input_path',request.input_path,'state',"running", ...
    'message',"AAA workflow started."));

workflow_name = select_workflow(request.mode, request.scope);
try
    receipt = feval(workflow_name, request, runtime);
    receipt = normalize_receipt(receipt, request);
    aaa.io.emit_event(runtime, struct( ...
        'type',"workflow",'mode',request.mode,'scope',request.scope, ...
        'input_path',request.input_path,'state',receipt.status, ...
        'message',"AAA workflow finished.",'output_path',receipt.output_path));
catch ME
    aaa.io.emit_event(runtime, struct( ...
        'type',"workflow",'mode',request.mode,'scope',request.scope, ...
        'input_path',request.input_path,'state',"failed", ...
        'message',string(ME.message),'details',exception_details(ME)));
    rethrow(ME);
end
end

function workflow_name = select_workflow(mode, scope)
key = string(mode) + ":" + string(scope);
switch key
    case "ap:single"
        workflow_name = 'aaa.workflows.AP_analysis3';
    case "ap:record"
        workflow_name = 'aaa.workflows.run_AP_analysis3_rec';
    case "dual:single"
        workflow_name = 'aaa.workflows.Dual_analysis3';
    case "dual:record"
        workflow_name = 'aaa.workflows.run_Dual_analysis3_rec';
    otherwise
        error('AAA:UnknownWorkflow', 'Unsupported AAA workflow: %s', key);
end
end

function receipt = normalize_receipt(receipt, request)
if nargin < 1 || isempty(receipt)
    receipt = struct();
end
if ~isstruct(receipt) || ~isscalar(receipt)
    error('AAA:InvalidReceipt', 'Workflow receipt must be one scalar struct.');
end
defaults = struct( ...
    'status',"completed", ...
    'mode',string(request.mode), ...
    'scope',string(request.scope), ...
    'input_path',string(request.input_path), ...
    'output_path',string(request.workflow.output_path), ...
    'manifest_file',"", ...
    'completed_sections',strings(0,1), ...
    'failed_section',"", ...
    'error',struct('identifier',"",'message',"",'stack',struct([])));
names = fieldnames(defaults);
for idx = 1:numel(names)
    name = names{idx};
    if ~isfield(receipt, name) || isempty(receipt.(name))
        receipt.(name) = defaults.(name);
    end
end
receipt.status = lower(string(receipt.status));
receipt.mode = string(receipt.mode);
receipt.scope = string(receipt.scope);
receipt.input_path = string(receipt.input_path);
receipt.output_path = string(receipt.output_path);
receipt.manifest_file = string(receipt.manifest_file);
receipt.completed_sections = string(receipt.completed_sections(:));
receipt.failed_section = string(receipt.failed_section);
end

function details = exception_details(ME)
details = struct('identifier',string(ME.identifier), ...
    'message',string(ME.message),'stack',{ME.stack});
end
