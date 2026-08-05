function [ctx, record] = save_analysis_results(ctx, section, outcome)
%SAVE_ANALYSIS_RESULTS Write one section into the unified AAA result schema.
% Every successful/reused section enters through this function. Scientific
% movie pixels are never copied into a result file.

if nargin < 3 || isempty(outcome), outcome = struct(); end
if ~isstruct(ctx) || ~isscalar(ctx) || ~isstruct(outcome) || ~isscalar(outcome)
    error('AAA:Results:InvalidSaveInput', ...
        'ctx and outcome must be scalar structs.');
end
section = lower(strtrim(string(section)));
if ~isscalar(section) || strlength(section) == 0
    error('AAA:Results:InvalidSection','section must be nonempty scalar text.');
end
output_path = absolute_path(string(ctx.output_path));
if strlength(output_path) == 0
    error('AAA:Results:MissingOutputPath', ...
        'A unified result requires ctx.output_path.');
end
results_path = string(fullfile(output_path,'results'));
if ~isfolder(results_path), mkdir(results_path); end

% A trace-only reuse still needs a standalone raw-trace dependency. Store
% that dependency once as roi.mat instead of duplicating raw in trace.mat.
if section=="trace" && ~isfile(fullfile(results_path,'roi.mat')) ...
        && context_has_raw(ctx)
    dependency_outcome=struct('state',"reused",'action',"reuse", ...
        'source_files',outcome_strings(outcome,'source_files'));
    if isfield(outcome,'compatibility_report')
        dependency_outcome.compatibility_report=outcome.compatibility_report;
    end
    [ctx,~]=save_analysis_results(ctx,"roi",dependency_outcome);
end

roles = channel_roles(ctx);
saved_at = datetime('now');
state = outcome_field(outcome,'state',"completed");
action = outcome_field(outcome,'action',"run");
if action == "reuse" && state ~= "failed", state = "reused"; end
source_results = "";
if isfield(ctx,'request') && isfield(ctx.request,'workflow') ...
        && isfield(ctx.request.workflow,'source_results_path')
    value = string(ctx.request.workflow.source_results_path);
    if strlength(value) > 0, source_results = absolute_path(value); end
end
source_files = outcome_strings(outcome,'source_files');
compatibility_report = empty_compatibility_report();
if isfield(outcome,'compatibility_report') ...
        && isstruct(outcome.compatibility_report) ...
        && ~isempty(outcome.compatibility_report)
    compatibility_report = normalize_compatibility_report( ...
        outcome.compatibility_report);
end

result = struct( ...
    'schema',"AAA.result.v1", ...
    'section',section, ...
    'mode',string(ctx.mode), ...
    'scope',string(ctx.scope), ...
    'status',state, ...
    'roles',roles, ...
    'data',section_payload(ctx,section), ...
    'provenance',struct( ...
        'action',action, ...
        'input_path',absolute_path(string(ctx.input_path)), ...
        'output_path',output_path, ...
        'source_results',source_results, ...
        'source_files',absolute_paths(source_files,source_results), ...
        'compatibility_report',compatibility_report), ...
    'output_files',section_output_files(ctx,section,outcome,output_path), ...
    'saved_at',saved_at);

result_file = absolute_path(string(fullfile(results_path,section + ".mat")));
atomic_save(result_file,result);
record = result_record(section,result_file,output_path,result,source_results);
ctx = register_result(ctx,section,record,result.output_files);
update_manifest(ctx,record,result);
emit_save_event(ctx,record);
end

function data = section_payload(ctx,section)
data = struct();
switch section
    case "input"
        data = role_payload(ctx,@input_payload);
        if isfield(ctx,'input_layout'), data.input_layout = ctx.input_layout; end
        if isfield(ctx,'input_geometry'), data.input_geometry = ctx.input_geometry; end
    case "motion"
        data = role_payload(ctx,@motion_payload);
    case "registration"
        data = shared_result(ctx,'registration');
    case "roi"
        data.roi = shared_result(ctx,'roi');
        maps = shared_result(ctx,'maps');
        if ~isempty(fieldnames(maps)), data.maps = maps; end
        data = merge_struct(data,role_payload(ctx,@roi_payload));
    case "trace"
        data = role_payload(ctx,@trace_payload);
        details=shared_result(ctx,'trace_details');
        if ~isempty(fieldnames(details)),data.details=details;end
    case "peak"
        data = role_payload(ctx,@peak_payload);
    case "stim"
        data = shared_result(ctx,'stim');
    case "comparison"
        data = shared_result(ctx,'comparison');
    case "frequency"
        data = shared_result(ctx,'frequency');
    case "visualization"
        data = shared_output(ctx,'visualization');
    case "export"
        execution=ctx.execution;
        execution.status="completed";
        execution.completed_at=datetime('now');
        execution.current_section="";
        execution.current_action="";
        execution.current_channel="";
        data = struct('request',ctx.request,'execution',execution, ...
            'index',shared_result(ctx,'export'), ...
            'result_index',shared_result(ctx,'result_index'));
    case "record"
        data=shared_result(ctx,'record');
    otherwise
        error('AAA:Results:UnknownSection', ...
            'No unified save rule exists for section %s.',section);
end
end

function tf=context_has_raw(ctx)
tf=false;
for idx=1:numel(ctx.channels)
    if isfield(ctx.channels(idx).results,'trace_results') ...
            && isfield(ctx.channels(idx).results.trace_results,'raw')
        tf=true;return;
    end
end
end

function value = input_payload(channel)
value = struct('source',channel.source,'movie_info',struct());
if isfield(channel,'data') && isfield(channel.data,'movie_info')
    value.movie_info = channel.data.movie_info;
elseif isfield(channel,'results') && isfield(channel.results,'movie_info')
    value.movie_info = channel.results.movie_info;
end
end

function value = motion_payload(channel)
value = struct('movie_info',struct(),'motion',struct());
if isfield(channel,'results') && isfield(channel.results,'movie_info')
    value.movie_info = channel.results.movie_info;
end
if isfield(channel,'results') && isfield(channel.results,'motion')
    value.motion = channel.results.motion;
elseif isfield(value.movie_info,'motion')
    value.motion = struct('info',value.movie_info.motion);
end
end

function value = roi_payload(channel)
value = struct('movie_info',struct(),'trace_results',struct());
if isfield(channel.results,'movie_info'), value.movie_info=channel.results.movie_info; end
if isfield(channel.results,'trace_results') ...
        && isfield(channel.results.trace_results,'raw')
    value.trace_results.raw = channel.results.trace_results.raw;
end
end

function value = trace_payload(channel)
value = struct('movie_info',struct(),'trace_results',struct());
if isfield(channel.results,'movie_info'), value.movie_info=channel.results.movie_info; end
if isfield(channel.results,'trace_results')
    value.trace_results = channel.results.trace_results;
    if isfield(value.trace_results,'raw')
        value.trace_results = rmfield(value.trace_results,'raw');
    end
end
end

function value = peak_payload(channel)
value = struct();
if isfield(channel.results,'peak_results')
    value.peak_results = channel.results.peak_results;
end
end

function output = role_payload(ctx,builder)
output = struct();
for idx = 1:numel(ctx.channels)
    role = char(lower(string(ctx.channels(idx).profile.role)));
    output.(role) = builder(ctx.channels(idx));
end
end

function value = shared_result(ctx,name)
value = struct();
if isfield(ctx,'shared') && isfield(ctx.shared,'results') ...
        && isfield(ctx.shared.results,name)
    value = ctx.shared.results.(name);
end
end

function value = shared_output(ctx,name)
value = struct();
if isfield(ctx,'shared') && isfield(ctx.shared,'output_files') ...
        && isfield(ctx.shared.output_files,name)
    value = ctx.shared.output_files.(name);
end
end

function files = section_output_files(ctx,section,outcome,output_path)
values = {};
if isfield(outcome,'output_files'), values{end+1}=outcome.output_files; end %#ok<AGROW>
if isfield(ctx,'shared') && isfield(ctx.shared,'output_files') ...
        && isfield(ctx.shared.output_files,char(section))
    values{end+1}=ctx.shared.output_files.(char(section)); %#ok<AGROW>
end
for idx=1:numel(ctx.channels)
    if isfield(ctx.channels(idx),'output_files') ...
            && isfield(ctx.channels(idx).output_files,char(section))
        values{end+1}=ctx.channels(idx).output_files.(char(section)); %#ok<AGROW>
    end
end
files=strings(0,1);
for idx=1:numel(values)
    files=[files;collect_paths(values{idx})]; %#ok<AGROW>
end
files=absolute_paths(files,output_path);
inside=false(size(files));
root=lower(char(output_path+filesep));
for idx=1:numel(files)
    inside(idx)=startsWith(lower(char(files(idx))),root);
end
files=unique(files(inside & arrayfun(@isfile,files)),'stable');
end

function files=collect_paths(value)
files=strings(0,1);
if ischar(value) || isstring(value)
    files=string(value(:));
elseif iscell(value)
    for idx=1:numel(value),files=[files;collect_paths(value{idx})];end %#ok<AGROW>
elseif isstruct(value)
    for item=1:numel(value)
        names=fieldnames(value(item));
        for idx=1:numel(names)
            files=[files;collect_paths(value(item).(names{idx}))]; %#ok<AGROW>
        end
    end
end
files=files(strlength(files)>0);
end

function ctx=register_result(ctx,section,record,output_files)
if ~isfield(ctx.shared,'output_files') || ~isstruct(ctx.shared.output_files)
    ctx.shared.output_files=struct();
end
ctx.shared.output_files.(char(section))=output_files;
if ~isfield(ctx.shared.results,'result_index') ...
        || ~isstruct(ctx.shared.results.result_index)
    ctx.shared.results.result_index=struct();
end
ctx.shared.results.result_index.(char(section))=record;
end

function update_manifest(ctx,record,result)
manifest_file=string(fullfile(ctx.output_path,'analysis_manifest.mat'));
if ~isfile(manifest_file),return;end
manifest=aaa.io.load_manifest(manifest_file);
if ~isfield(manifest,'results') || ~isstruct(manifest.results)
    manifest.results=repmat(record,0,1);
end
if ~isempty(manifest.results)
    keep=string({manifest.results.section})~=record.section;
    manifest.results=manifest.results(keep);
end
manifest.results(end+1,1)=record;
manifest.results_schema="AAA.results.v1";
manifest.output_files=unique([string(manifest.output_files(:)); ...
    record.file;string(result.output_files(:))],'stable');
report=result.provenance.compatibility_report;
if strlength(report.adapter_used)>0
    if ~isfield(manifest,'compatibility_history') ...
            || ~isstruct(manifest.compatibility_history)
        manifest.compatibility_history=repmat(report,0,1);
    end
    append=true;
    if ~isempty(manifest.compatibility_history)
        same_source=string({manifest.compatibility_history.absolute_source_path})' ...
            ==string(report.absolute_source_path);
        same_adapter=string({manifest.compatibility_history.adapter_used})' ...
            ==string(report.adapter_used);
        append=~any(same_source&same_adapter);
    end
    if append,manifest.compatibility_history(end+1,1)=report;end
end
aaa.io.save_manifest(manifest,manifest_file);
end

function record=result_record(section,file,root,result,source_results)
record=struct('section',section,'file',file, ...
    'relative_file',relative_path(file,root),'schema',result.schema, ...
    'roles',result.roles,'status',result.status,'saved_at',result.saved_at, ...
    'source_results',source_results);
end

function atomic_save(file,result)
folder=string(fileparts(file));
temporary=string([tempname(folder),'.mat']);
cleanup=onCleanup(@()delete_if_present(temporary));
save(temporary,'result','-v7.3');
[ok,message]=movefile(temporary,file,'f');
if ~ok
    error('AAA:Results:SaveFailed','Could not replace %s: %s',file,message);
end
clear cleanup;
end

function emit_save_event(ctx,record)
if ~isfield(ctx,'runtime'),return;end
event=struct('type',"result",'mode',string(ctx.mode), ...
    'scope',string(ctx.scope),'input_path',string(ctx.input_path), ...
    'output_path',string(ctx.output_path),'section',record.section, ...
    'state',"saved",'message',"Unified result saved: "+record.file, ...
    'timestamp',datetime('now'),'details',record);
aaa.io.emit_event(ctx.runtime,event);
end

function roles=channel_roles(ctx)
roles=strings(0,1);
if ~isfield(ctx,'channels') || isempty(ctx.channels),return;end
profiles=[ctx.channels.profile];
roles=string({profiles.role})';
end

function value=outcome_field(outcome,name,default)
value=default;
if isfield(outcome,name) && ~isempty(outcome.(name))
    value=string(outcome.(name));
end
end

function values=outcome_strings(outcome,name)
values=strings(0,1);
if isfield(outcome,name),values=string(outcome.(name));values=values(:);end
values=values(strlength(values)>0);
end

function paths=absolute_paths(paths,base)
paths=string(paths(:));
for idx=1:numel(paths)
    value=paths(idx);
    if ~is_absolute(value) && strlength(string(base))>0
        value=fullfile(base,value);
    end
    paths(idx)=absolute_path(value);
end
paths=unique(paths(strlength(paths)>0),'stable');
end

function value=absolute_path(value)
value=strtrim(string(value));
if strlength(value)==0,return;end
if ~is_absolute(value),value=fullfile(pwd,value);end
try,value=string(char(java.io.File(char(value)).getCanonicalPath()));
catch,value=string(char(java.io.File(char(value)).getAbsolutePath()));end
end

function tf=is_absolute(value)
text=char(string(value));
tf=~isempty(regexp(text,'^[A-Za-z]:[\\/]','once')) ...
    || startsWith(text,'\\') || startsWith(text,'/');
end

function value=relative_path(file,root)
file=string(file);root=string(root);
prefix=root+filesep;
if startsWith(lower(file),lower(prefix))
    value=extractAfter(file,strlength(prefix));
else
    value=file;
end
end

function report=normalize_compatibility_report(value)
report=empty_compatibility_report();
names=fieldnames(report);
for idx=1:numel(names)
    if isfield(value,names{idx}),report.(names{idx})=value.(names{idx});end
end
end

function report=empty_compatibility_report()
report=struct('source_path',"",'absolute_source_path',"", ...
    'detected_schema',"",'adapter_used',"", ...
    'source_files',strings(0,1),'absolute_source_files',strings(0,1), ...
    'relative_source_files',strings(0,1),'frame_rate_source',struct(), ...
    'mapped_fields',strings(0,1),'missing_fields',strings(0,1), ...
    'warnings',strings(0,1),'loaded_at',NaT);
end

function output=merge_struct(base,patch)
output=base;names=fieldnames(patch);
for idx=1:numel(names),output.(names{idx})=patch.(names{idx});end
end

function delete_if_present(path)
if isfile(path),delete(path);end
end
