function receipt = run_single(request,runtime)
%RUN_SINGLE Execute one AP/Dual request through the section runtime.
% This is the sole single-run backend used by the App and record runners.

if nargin < 2 || isempty(runtime), runtime = struct(); end
[preflight,preflight_match]=aaa.helpers.common.runtime_preflight(runtime,request);
reuse_resolution=struct('attempted',false,'source_path',"",'report',struct(), ...
    'frame_rate_overrides',repmat(struct('role',"",'ui_value',NaN,'saved_value',NaN),0,1));
if preflight_match=="exact"&&isfield(preflight,'reuse_resolution') ...
        &&isstruct(preflight.reuse_resolution) ...
        &&isfield(preflight.reuse_resolution,'attempted')
    reuse_resolution=preflight.reuse_resolution;
else
    [request,reuse_resolution]=aaa.helpers.common.resolve_request_reuse_source(request);
end
[request,report] = aaa.schema.validate_request(request);
if ~report.valid
    error('AAA:Workflow:InvalidRequest','%s',strjoin(report.errors,newline));
end
if string(request.scope) ~= "single"
    error('AAA:Workflow:WrongScope','run_single requires scope="single".');
end
if strlength(preflight_match)==0
    input_inspection=aaa.helpers.common.inspect_input_structure( ...
        request,report.execution_plan);
    if ~input_inspection.valid
        error('AAA:Workflow:InvalidInputStructure','%s', ...
            strjoin(input_inspection.errors,newline));
    end
else
    aaa.io.emit_event(runtime,struct('type',"preflight",'mode',request.mode, ...
        'scope',request.scope,'input_path',request.input_path,'state',"reused", ...
        'message',"Matching authoritative preflight reused; backend input scan skipped."));
end
request = bind_output_path(request);
[request,report] = aaa.schema.validate_request(request);
plan = report.execution_plan;

% These AP stages were explicitly excluded from this migration phase.
excluded = strings(0,1);
if string(request.mode) == "ap"
    for name = ["ap_events","ap_statistics"]
        if plan.(char(name)) ~= "skip", excluded(end+1,1) = name; end %#ok<AGROW>
        plan.(char(name)) = "skip";
    end
end

output_path = string(request.workflow.output_path);
if ~isfolder(output_path), mkdir(output_path); end
[runtime,logger,log_file,owns_logger] = ...
    aaa.io.attach_analysis_logger(runtime,request,output_path);
cleanup = onCleanup(@() close_logger(logger,owns_logger)); %#ok<NASGU>
if reuse_resolution.attempted
    if strlength(reuse_resolution.source_path) > 0
        logger.write("INFO", "Reuse source resolved: " + reuse_resolution.source_path);
    else
        logger.write("ERROR", "Reuse source auto-detection found no compatible result under: " + ...
            string(reuse_resolution.report.search_root));
    end
    for idx=1:numel(reuse_resolution.frame_rate_overrides)
        item=reuse_resolution.frame_rate_overrides(idx);
        logger.write("INFO","Reuse frame rate authority: role="+item.role+ ...
            ", saved="+item.saved_value+" Hz, UI/request before reuse="+ ...
            item.ui_value+" Hz. Saved value is used downstream.");
    end
end
[~,manifest_file] = aaa.io.create_manifest(request,plan,output_path);
update_manifest(manifest_file,struct('status',"running", ...
    'started_at',datetime("now"),'log_file',string(log_file), ...
    'excluded_sections',excluded));

ctx = aaa.sections.create_context(request,runtime);
entries = aaa.sections.registry(request.mode);
[ctx,receipt] = aaa.sections.run_plan(ctx,plan,entries);
receipt.manifest_file = string(manifest_file);
receipt.roi_file = expected_roi_file(request,output_path);
receipt.output_files = list_output_files(output_path);
update_manifest(manifest_file,struct( ...
    'status',receipt.status, ...
    'completed_sections',receipt.completed_sections, ...
    'failed_section',receipt.failed_section, ...
    'section_status',receipt.section_status, ...
    'completed_at',receipt.completed_at, ...
    'output_files',receipt.output_files, ...
    'input_sources',manifest_input_sources(ctx), ...
    'excluded_sections',excluded, ...
    'error',receipt.error));
if receipt.status == "failed"
    if isfield(ctx.execution,'exception') && isa(ctx.execution.exception,'MException')
        rethrow(ctx.execution.exception);
    end
    error('AAA:Workflow:SectionFailed','Section %s failed: %s', ...
        receipt.failed_section,receipt.error.message);
end
end

function sources=manifest_input_sources(ctx)
template=struct('role',"",'camera_index',NaN,'logical_path',"", ...
    'resolved_path',"",'physical_files',strings(0,1), ...
    'is_split',false,'part_count',0,'resolution_method',"", ...
    'frame_rate',NaN,'frame_count',0,'reuse_source_path',"", ...
    'reuse_result_file',"", ...
    'movie_loaded',false);
sources=repmat(template,numel(ctx.channels),1);
for idx=1:numel(ctx.channels)
    channel=ctx.channels(idx); item=template;
    item.role=string(channel.profile.role);
    item.camera_index=channel.camera_index;
    if isfield(channel,'data') && isfield(channel.data,'movie_info')
        info=channel.data.movie_info;
        item=copy_if_present(item,info,'logical_path','logical_source_path');
        item=copy_if_present(item,info,'resolved_path','source_path');
        item=copy_if_present(item,info,'physical_files','source_files');
        item=copy_if_present(item,info,'is_split','source_is_split');
        item=copy_if_present(item,info,'part_count','source_part_count');
        item=copy_if_present(item,info,'resolution_method','source_resolution_method');
        item=copy_if_present(item,info,'frame_rate','frame_rate');
        item=copy_if_present(item,info,'frame_count','frame_count');
        item=copy_if_present(item,info,'reuse_source_path','reuse_source_path');
        if isfield(info,'movie_loaded'),item.movie_loaded=logical(info.movie_loaded);
        else,item.movie_loaded=isfield(channel.data,'movie_3d')&&~isempty(channel.data.movie_3d);end
    end
    if isfield(ctx,'input_layout')&&isstruct(ctx.input_layout) ...
            &&isfield(ctx.input_layout,'source_files') ...
            &&~isempty(ctx.input_layout.source_files)
        item.reuse_result_file=string(ctx.input_layout.source_files(1));
    end
    sources(idx)=item;
end
end

function target=copy_if_present(target,source,target_name,source_name)
if isstruct(source)&&isfield(source,source_name)
    target.(target_name)=source.(source_name);
end
end

function request = bind_output_path(request)
if strlength(string(request.workflow.run_name)) == 0
    request.workflow.run_name = "AAA_" + string(datetime('now','Format','yyyyMMdd_HHmmss_SSS'));
end
if strlength(string(request.workflow.output_path)) > 0, return; end
input_path = string(request.input_path);
if isfile(input_path), base = string(fileparts(input_path)); else, base = input_path; end
folder = upper(string(request.mode)) + "_analysis3";
request.workflow.output_path = fullfile(base,folder,request.workflow.run_name);
end

function close_logger(logger,owns_logger)
if owns_logger && ~isempty(logger) && isa(logger,'aaa.io.AnalysisLogger')
    logger.close();
end
end

function update_manifest(file_path,patch)
if strlength(string(file_path)) > 0 && isfile(file_path)
    aaa.io.update_manifest(file_path,patch);
end
end

function file_path = expected_roi_file(request,output_path) %#ok<INUSD>
file_path = string(fullfile(output_path,'results','roi.mat'));
if ~isfile(file_path), file_path = ""; end
end

function files = list_output_files(output_path)
listing = dir(fullfile(output_path,'**','*'));
listing = listing(~[listing.isdir]);
files = string(fullfile({listing.folder},{listing.name}))';
end
