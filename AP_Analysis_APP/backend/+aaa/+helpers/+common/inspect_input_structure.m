function inspection = inspect_input_structure(request,plan,options)
%INSPECT_INPUT_STRUCTURE Preflight Cycle/Record inputs without loading movies.
% Manual reuse roots are authoritative. Automatic candidates are still
% reported for comparison but never replace a nonempty manual path.

if nargin<3||isempty(options),options=struct();end
if ~isfield(options,'scan_reuse'),options.scan_reuse=true;end
if ~isfield(options,'progress'),options.progress=[];end
scan_reuse=logical(options.scan_reuse);

inspection=struct('valid',true,'mode',string(request.mode), ...
    'requested_scope',string(request.scope),'inferred_scope',"", ...
    'input_path',string(request.input_path),'manual_reuse_root', ...
    string(request.workflow.source_results_path),'errors',strings(0,1), ...
    'warnings',strings(0,1),'cycles',repmat(empty_cycle(),0,1));
input_path=strtrim(string(request.input_path));
if ~(isfile(input_path)||isfolder(input_path))
    inspection.errors(end+1,1)="Input path does not exist: "+input_path;
    inspection.valid=false;
    return;
end

[items,inferred,discovery_error]=inspection_items(request,input_path);
inspection.inferred_scope=inferred;
if strlength(discovery_error)>0
    inspection.errors(end+1,1)=discovery_error;
    inspection.valid=false;
    return;
end

requires_movie=isfield(plan,'input') && string(plan.input)=="run";
requires_reuse=plan_requires_source(plan,request.workflow);
manual_root=strtrim(string(request.workflow.source_results_path));
auto_detect=true;
if isfield(request.workflow,'auto_detect_reuse_source')
    auto_detect=logical(request.workflow.auto_detect_reuse_source);
end
[~,~,~,inventory_plan]=aaa.schema.resolve_sections( ...
    request.mode,"reuse_trace",struct());

for idx=1:numel(items)
    item=items(idx);
    report_progress(options,sprintf('%s (%d/%d): checking movie inputs.', ...
        string(item.label),idx,numel(items)));
    cycle=empty_cycle();
    cycle.label=string(item.label);
    cycle.path=string(item.cycle_path);
    cycle.empty=is_empty_path(cycle.path);
    [movie_sources,movie_report]=aaa.helpers.common.resolve_movie_sources( ...
        cycle.path,channel_specs(request));
    cycle.standard_input=movie_report.valid;
    cycle.movie_sources=movie_sources;
    cycle.input_evidence=movie_evidence(movie_report);
    if ~movie_report.valid
        cycle.input_resolution_errors=movie_report.errors;
    end

    if string(request.scope)=="single"
        target=struct('cycle_name',"",'label',"");
    else
        target=struct('cycle_name',string(item.cycle_name), ...
            'label',string(item.label));
    end
    if scan_reuse
        report_progress(options,sprintf('%s (%d/%d): searching saved results.', ...
            string(item.label),idx,numel(items)));
        [auto_source,auto_report]=aaa.helpers.common.resolve_reuse_source( ...
            request.mode,cycle.path,choose_plan(requires_reuse,plan,inventory_plan), ...
            target,request.workflow);
        cycle.auto_reuse_candidate=auto_source;
        cycle.auto_candidate_status=string(auto_report.status);
        cycle.valid_output_count=sum([auto_report.candidates.valid]);
        cycle.output_candidate_count=numel(auto_report.candidates);
        cycle.auto_selection_rule=string(auto_report.selection_rule);
        cycle.auto_frame_rates=chosen_frame_rates(auto_report);
        inspection.warnings=[inspection.warnings;string(auto_report.warnings(:))];

        if strlength(manual_root)>0
            report_progress(options,sprintf('%s (%d/%d): checking manual Reuse source.', ...
                string(item.label),idx,numel(items)));
            [manual_source,manual_report]=aaa.helpers.common.resolve_reuse_source( ...
                request.mode,manual_root,choose_plan(requires_reuse,plan,inventory_plan), ...
                target,request.workflow);
            cycle.manual_reuse_candidate=manual_source;
            cycle.manual_candidate_status=string(manual_report.status);
            cycle.reuse_selection="manual";
            cycle.selected_reuse_source=manual_source;
            cycle.manual_overrode_auto=strlength(manual_source)>0 ...
                && strlength(auto_source)>0 && manual_source~=auto_source;
            cycle.manual_selection_rule=string(manual_report.selection_rule);
            cycle.manual_frame_rates=chosen_frame_rates(manual_report);
            inspection.warnings=[inspection.warnings;string(manual_report.warnings(:))];
        elseif requires_reuse && auto_detect
            cycle.reuse_selection="automatic";
            cycle.selected_reuse_source=auto_source;
        elseif requires_reuse
            cycle.reuse_selection="manual_required";
        else
            cycle.reuse_selection="inventory_only";
        end
    else
        cycle.reuse_selection="pending";
    end

    if cycle.empty && requires_movie
        inspection.errors(end+1,1)=cycle.label+": Cycle/input is empty."; %#ok<AGROW>
    elseif requires_movie && ~cycle.standard_input
        detail="";
        if ~isempty(cycle.input_resolution_errors)
            detail=" Details: "+strjoin(cycle.input_resolution_errors," | ");
        end
        inspection.errors(end+1,1)=cycle.label+ ...
            ": no standard movie input was found. This workflow cannot start. "+ ...
            "Choose reuse_trace/analysis_only and enter a manual Reuse source, "+ ...
            "or correct the acquisition folder."+detail; %#ok<AGROW>
    end
    if scan_reuse&&requires_reuse&&strlength(cycle.selected_reuse_source)==0
        if strlength(manual_root)>0
            inspection.errors(end+1,1)=cycle.label+ ...
                ": manual Reuse source has no compatible completed result: "+manual_root; %#ok<AGROW>
        elseif ~auto_detect
            inspection.errors(end+1,1)=cycle.label+ ...
                ": automatic reuse detection is disabled; enter Reuse source manually."; %#ok<AGROW>
        else
            inspection.errors(end+1,1)=cycle.label+ ...
                ": no compatible completed result was auto-detected. "+ ...
                "Enter Reuse source manually."; %#ok<AGROW>
        end
    end
    inspection.cycles(end+1,1)=cycle; %#ok<AGROW>
end
inspection.errors=unique(inspection.errors,'stable');
inspection.warnings=unique(inspection.warnings,'stable');
inspection.valid=isempty(inspection.errors);
end

function report_progress(options,message)
if ~isfield(options,'progress')||~isa(options.progress,'function_handle'),return;end
try
    options.progress(string(message));
catch
end
end

function [items,inferred,message]=inspection_items(request,input_path)
items=repmat(struct('label',"",'cycle_name',"",'cycle_path',""),0,1);
message="";
if isfile(input_path)
    inferred="single";
    items=struct('label',last_part(input_path),'cycle_name',"",'cycle_path',input_path);
    return;
end
[~,name]=fileparts(input_path);
if startsWith(string(name),"Cycle",'IgnoreCase',true)
    inferred="single";
    items=struct('label',string(name),'cycle_name',string(name),'cycle_path',input_path);
    return;
end
try
    found=aaa.helpers.common.discover_record_cycles(input_path,[],true);
    inferred="record";
    items=repmat(struct('label',"",'cycle_name',"",'cycle_path',""),numel(found),1);
    for idx=1:numel(found)
        items(idx).label=found(idx).label;
        items(idx).cycle_name=found(idx).cycle_name;
        items(idx).cycle_path=found(idx).cycle_path;
    end
catch
    inferred="single";
    if string(request.scope)=="record"
        message="Record input contains no direct Cycle* or nested Rec*/Cycle* folders: "+input_path;
    else
        items=struct('label',string(name),'cycle_name',"",'cycle_path',input_path);
    end
end
end

function specs=channel_specs(request)
if string(request.mode)=="ap"
    specs=struct('camera_index',double(request.params.ap.camera_index), ...
        'role',"voltage");
    return;
end
camera_roles=[string(request.params.dual.camera1_role), ...
    string(request.params.dual.camera2_role)];
roles=["voltage","calcium"];
specs=repmat(struct('camera_index',NaN,'role',""),1,2);
for idx=1:2
    specs(idx).role=roles(idx);
    specs(idx).camera_index=find(camera_roles==roles(idx),1,'first');
end
end

function evidence=movie_evidence(report)
evidence=strings(0,1);
if ~report.valid,return;end
for idx=1:numel(report.sources)
    source=report.sources(idx);
    evidence(end+1,1)=source.role+" camera "+source.camera_index+ ...
        ": "+source.path+"; parts="+source.part_count+ ...
        "; method="+source.resolution_method; %#ok<AGROW>
end
end

function values=chosen_frame_rates(report)
values=struct();
if ~isstruct(report)||~isfield(report,'chosen_inspection') ...
        || ~isstruct(report.chosen_inspection) ...
        || ~isfield(report.chosen_inspection,'saved_frame_rates')
    return;
end
values=report.chosen_inspection.saved_frame_rates;
end

function tf=is_empty_path(path_value)
if isfile(path_value),tf=false;return;end
listing=dir(path_value);
names=string({listing.name});
tf=~any(names~="." & names~=".." & ~startsWith(names,"."));
end

function tf=plan_requires_source(plan,workflow)
source_plan=plan;
if isfield(source_plan,'roi') && source_plan.roi=="reuse" ...
        && strlength(string(workflow.roi_path))>0
    source_plan.roi="skip";
end
tf=any(structfun(@(value)string(value)=="reuse",source_plan));
end

function plan=choose_plan(use_requested,requested,inventory)
if use_requested,plan=requested;else,plan=inventory;end
end

function value=last_part(path_value)
[~,name,ext]=fileparts(path_value);value=string([name,ext]);
end

function cycle=empty_cycle()
cycle=struct('label',"",'path',"",'empty',false,'standard_input',false, ...
    'input_evidence',strings(0,1),'input_resolution_errors',strings(0,1), ...
    'movie_sources',{{}},'output_candidate_count',0, ...
    'valid_output_count',0,'auto_reuse_candidate',"", ...
    'auto_candidate_status',"",'manual_reuse_candidate',"", ...
    'manual_candidate_status',"",'selected_reuse_source',"", ...
    'reuse_selection',"",'manual_overrode_auto',false, ...
    'auto_selection_rule',"",'manual_selection_rule',"", ...
    'auto_frame_rates',struct(),'manual_frame_rates',struct());
end
