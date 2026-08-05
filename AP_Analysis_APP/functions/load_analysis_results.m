function [results, report] = load_analysis_results(source, spec)
%LOAD_ANALYSIS_RESULTS Read current AAA results or adapt historical outputs.
% This file is the only compatibility boundary for legacy AP/Dual result
% containers. File names are discovery hints; payload structure decides use.

if nargin < 2 || isempty(spec), spec=struct(); end
if ~isstruct(spec) || ~isscalar(spec)
    error('AAA:Results:InvalidLoadSpec','spec must be one scalar struct.');
end
source=strtrim(string(source));
if ~isscalar(source) || strlength(source)==0 || ~(isfolder(source)||isfile(source))
    error('AAA:Results:MissingSource','Result source does not exist: %s',source);
end
absolute_source=absolute_path(source);
report=empty_report(source,absolute_source);

[files,is_current]=current_result_files(absolute_source,spec);
if is_current
    [results,report]=load_current_results(files,report);
else
    [results,report]=adapt_legacy_results(absolute_source,spec,report);
end
report.loaded_at=datetime('now');
report.source_files=report.absolute_source_files;
report.relative_source_files=relative_paths( ...
    report.absolute_source_files,folder_root(absolute_source));
emit_report(spec,report);
end

function [files,found]=current_result_files(source,spec)
files=strings(0,1);found=false;
if isfile(source)
    try
        names=string({whos('-file',char(source)).name});
        if ismember("result",names),files=source;found=true;end
    catch
    end
    return;
end
results_dir=string(fullfile(source,'results'));
if ~isfolder(results_dir),return;end
listing=dir(fullfile(results_dir,'*.mat'));
if isempty(listing),return;end
files=string(fullfile({listing.folder},{listing.name}))';
if isfield(spec,'sections')&&~isempty(spec.sections)
    wanted=lower(string(spec.sections(:)));
    names=strings(size(files));
    for idx=1:numel(files),[~,name]=fileparts(files(idx));names(idx)=lower(string(name));end
    files=files(ismember(names,wanted));
    if isempty(files),return;end
end
order=current_section_order(files);
files=files(order);found=true;
end

function order=current_section_order(files)
canonical=["input","motion","registration","roi","trace","peak", ...
    "stim","comparison","frequency","visualization","export"];
rank=zeros(numel(files),1)+numel(canonical)+1;
for idx=1:numel(files)
    [~,name]=fileparts(files(idx));
    found=find(canonical==lower(string(name)),1);
    if ~isempty(found),rank(idx)=found;end
end
[~,order]=sortrows([rank,(1:numel(files))']);
end

function [results,report]=load_current_results(files,report)
results=struct();
for idx=1:numel(files)
    loaded=load(files(idx),'result');
    if ~isfield(loaded,'result') || ~isstruct(loaded.result) ...
            || ~isscalar(loaded.result) || ~isfield(loaded.result,'schema') ...
            || string(loaded.result.schema)~="AAA.result.v1"
        error('AAA:Results:InvalidCurrentResult', ...
            'Invalid unified result file: %s',files(idx));
    end
    results=merge_current_section(results,loaded.result);
    report.absolute_source_files(end+1,1)=absolute_path(files(idx)); %#ok<AGROW>
    report.mapped_fields(end+1,1)=string(loaded.result.section); %#ok<AGROW>
end
report.detected_schema="AAA.results.v1";
report.adapter_used="load_current_results";
report.absolute_source_files=unique(report.absolute_source_files,'stable');
report.mapped_fields=unique(report.mapped_fields,'stable');
report.frame_rate_source=frame_rate_sources(results,report.absolute_source_files);
end

function results=merge_current_section(results,result)
section=lower(string(result.section));data=result.data;
switch section
    case {"input","motion","peak"}
        results=merge_roles(results,data);
    case "trace"
        results=merge_roles(results,data);
        if isfield(data,'details'),results.trace_details=data.details;end
    case "roi"
        results=merge_roles(results,data);
        if isfield(data,'roi'),results.roi=data.roi;end
        if isfield(data,'maps'),results.maps=data.maps;end
    case "stim"
        results.stim=data;
    case "registration"
        results.registration=data;
    case "comparison"
        results.comparison=data;
    case "frequency"
        results.frequency=data;
    case "visualization"
        results.visualization=data;
    case "export"
        results.export=data;
    otherwise
        results.(char(section))=data;
end
end

function results=merge_roles(results,data)
for role=["voltage","calcium"]
    if ~isfield(data,char(role)),continue;end
    if ~isfield(results,char(role)),results.(char(role))=struct();end
    results.(char(role))=merge_struct_recursive( ...
        results.(char(role)),data.(char(role)));
end
end

function [results,report]=adapt_legacy_results(source,spec,report)
files=legacy_candidate_files(source,spec);
if isempty(files)
    error('AAA:Results:NoLegacyResults', ...
        'No MAT result payloads were found at %s.',source);
end
results=struct(); channel_candidates=struct('voltage',[], 'calcium',[]);
channel_meta=struct('voltage',empty_candidate(), 'calcium',empty_candidate());
recognized=false;
for idx=1:numel(files)
    file=files(idx);
    try
        header=whos('-file',char(file));
    catch ME
        report.warnings(end+1,1)="Unreadable MAT header: "+file+ ...
            " ("+string(ME.message)+")"; %#ok<AGROW>
        continue;
    end
    is_struct=strcmp({header.class},'struct');
    is_scalar=arrayfun(@(x)isequal(x.size,[1 1]),header);
    is_scalar=reshape(is_scalar,size(is_struct));
    scalar_structs=header(is_struct&is_scalar);
    file_loaded=struct();
    for var_idx=1:numel(scalar_structs)
        name=string(scalar_structs(var_idx).name);
        try
            loaded=load(file,char(name));value=loaded.(char(name));
        catch ME
            report.warnings(end+1,1)="Could not read "+file+"::"+name+ ...
                " ("+string(ME.message)+")"; %#ok<AGROW>
            continue;
        end
        [role,channel,score]=legacy_channel(value,name,spec);
        if strlength(role)>0
            candidate=struct('value',channel,'file',file,'variable',name, ...
                'score',score,'timestamp',file_time(file));
            [channel_candidates,channel_meta,accepted,warning]= ...
                choose_channel_candidate(channel_candidates,channel_meta,role,candidate);
            if strlength(warning)>0,report.warnings(end+1,1)=warning;end %#ok<AGROW>
            if accepted,recognized=true;end
            continue;
        end
        [results,mapped]=map_legacy_value(results,value,name);
        if strlength(mapped)>0
            recognized=true;report.mapped_fields(end+1,1)=mapped; %#ok<AGROW>
            report.absolute_source_files(end+1,1)=absolute_path(file); %#ok<AGROW>
        end
    end
    if any(string({header.name})=="rois")
        try,file_loaded=load(file);catch,file_loaded=struct();end
        [results,mapped]=adapt_legacy_roi_file(results,file_loaded);
        if strlength(mapped)>0
            recognized=true;report.mapped_fields(end+1,1)=mapped; %#ok<AGROW>
            report.absolute_source_files(end+1,1)=absolute_path(file); %#ok<AGROW>
        end
    end
end
for role=["voltage","calcium"]
    candidate=channel_candidates.(char(role));
    if ~isempty(candidate)
        if isfield(results,char(role))
            results.(char(role))=merge_struct_recursive( ...
                results.(char(role)),candidate);
        else
            results.(char(role))=candidate;
        end
        meta=channel_meta.(char(role));
        report.absolute_source_files(end+1,1)=absolute_path(meta.file); %#ok<AGROW>
        report.mapped_fields(end+1,1)=role; %#ok<AGROW>
    end
end
results=apply_manifest_sources(results,source,report);
if ~recognized || isempty(fieldnames(results))
    error('AAA:Results:UnrecognizedLegacyResults', ...
        'MAT files were found but no supported AP/Dual result schema was recognized at %s.',source);
end
report.detected_schema="legacy AP/Dual results";
report.adapter_used="adapt_legacy_results";
report.absolute_source_files=unique(report.absolute_source_files,'stable');
report.mapped_fields=unique(report.mapped_fields,'stable');
report.warnings=unique(report.warnings,'stable');
report.frame_rate_source=frame_rate_sources(results,report.absolute_source_files);
expected=expected_roles(spec);
for role=expected
    if ~isfield(results,char(role))
        report.missing_fields(end+1,1)=role; %#ok<AGROW>
    end
end
end

function files=legacy_candidate_files(source,spec)
max_files=500;
if isfield(spec,'max_mat_files'),max_files=double(spec.max_mat_files);end
if isfile(source),files=source;return;end
files=strings(0,1);
manifest_file=fullfile(source,'analysis_manifest.mat');
if isfile(manifest_file)
    try
        saved=load(manifest_file,'manifest');
        if isfield(saved,'manifest') && isfield(saved.manifest,'output_files')
            listed=string(saved.manifest.output_files(:));
            for idx=1:numel(listed)
                file=listed(idx);if ~isfile(file),file=fullfile(source,file);end
                if isfile(file) && lower(string(fileparts_extension(file)))==".mat"
                    files(end+1,1)=absolute_path(file); %#ok<AGROW>
                end
            end
        end
    catch
    end
end
listing=dir(fullfile(source,'*.mat'));
if ~isempty(listing)
    [~,order]=sort([listing.datenum],'descend');listing=listing(order);
    files=[files;string(fullfile({listing.folder},{listing.name}))']; %#ok<AGROW>
end
files=unique(files,'stable');
if numel(files)>max_files,files=files(1:max_files);end
end

function [role,channel,score]=legacy_channel(value,name,spec)
role="";channel=struct();score=0;
if ~isstruct(value)||~isscalar(value),return;end
if isfield(value,'trace_results')&&isstruct(value.trace_results)
    channel=value;score=100+numel(fieldnames(value.trace_results));
    if isfield(value,'peak_results'),score=score+10;end
    if isfield(value,'movie_info'),score=score+5;end
elseif has_trace_stages(value)
    channel=struct('trace_results',value);score=50+numel(fieldnames(value));
else
    return;
end
role=saved_role(channel,name,spec);
if strlength(role)==0,channel=struct();score=0;end
end

function role=saved_role(value,name,spec)
role="";
if isfield(value,'movie_info')&&isstruct(value.movie_info) ...
        && isfield(value.movie_info,'role')
    role=lower(strtrim(string(value.movie_info.role)));
end
if ~ismember(role,["voltage","calcium"])
    text=lower(string(name));
    if contains(text,"voltage"),role="voltage";
    elseif contains(text,"calcium"),role="calcium";end
end
if strlength(role)==0 && isfield(spec,'mode') ...
        && lower(string(spec.mode))=="ap"
    role="voltage";
end
if ~ismember(role,["voltage","calcium"]),role="";end
end

function tf=has_trace_stages(value)
tf=isstruct(value)&&isscalar(value)&&isfield(value,'raw') ...
    && isstruct(value.raw)&&isfield(value.raw,'data');
end

function [store,meta,accepted,warning]=choose_channel_candidate(store,meta,role,candidate)
accepted=false;warning="";name=char(role);current=meta.(name);
if strlength(current.file)==0
    store.(name)=candidate.value;meta.(name)=rmfield(candidate,'value');accepted=true;return;
end
replace=candidate.score>current.score ...
    || (candidate.score==current.score && candidate.timestamp>current.timestamp);
if replace
    warning="Multiple semantic "+role+" containers found; selected "+ ...
        candidate.file+"::"+candidate.variable+" by completeness, then newest time.";
    store.(name)=candidate.value;meta.(name)=rmfield(candidate,'value');accepted=true;
else
    warning="Multiple semantic "+role+" containers found; retained "+ ...
        current.file+"::"+current.variable+" by completeness, then newest time.";
end
end

function [results,mapped]=map_legacy_value(results,value,name)
mapped="";key=lower(string(name));
if key=="stim_results"
    if isfield(value,'artifacts') && ~isfield(value,'output_files')
        value.output_files=value.artifacts;value=rmfield(value,'artifacts');
    end
    results.stim=value;mapped="stim";
elseif key=="dual_results"
    if isfield(value,'visualizations') && ~isfield(value,'visualization')
        value.visualization=value.visualizations;value=rmfield(value,'visualizations');
    end
    fields=fieldnames(value);
    for idx=1:numel(fields),results.(fields{idx})=value.(fields{idx});end
    mapped="dual shared results";
elseif key=="time_frequency_results"
    results.frequency=value;mapped="frequency";
elseif key=="registration_info"
    results.registration=value;mapped="registration";
elseif key=="record_average"
    results.record=value;mapped="record";
elseif key=="record_grating_tuning"
    if ~isfield(results,'record'),results.record=struct();end
    results.record.grating=value;mapped="record.grating";
elseif key=="peak_results"
    if ~isfield(results,'voltage'),results.voltage=struct();end
    results.voltage.peak_results=value;mapped="voltage.peak_results";
end
end

function results=apply_manifest_sources(results,source,report) %#ok<INUSD>
if ~isfolder(source),return;end
file=fullfile(source,'analysis_manifest.mat');if ~isfile(file),return;end
try,saved=load(file,'manifest');catch,return;end
if ~isfield(saved,'manifest')||~isfield(saved.manifest,'input_sources') ...
        ||~isstruct(saved.manifest.input_sources),return;end
for idx=1:numel(saved.manifest.input_sources)
    item=saved.manifest.input_sources(idx);
    if ~isfield(item,'role'),continue;end
    role=lower(string(item.role));if ~ismember(role,["voltage","calcium"]) ...
            ||~isfield(results,char(role)),continue;end
    channel=results.(char(role));
    if ~isfield(channel,'movie_info'),channel.movie_info=struct();end
    if isfield(item,'frame_rate')&&isfinite(double(item.frame_rate)) ...
            &&double(item.frame_rate)>0
        if isfield(channel.movie_info,'frame_rate') ...
                &&isfinite(double(channel.movie_info.frame_rate))
            tolerance=max(1e-9,1e-9*max(abs([double(item.frame_rate), ...
                double(channel.movie_info.frame_rate)])));
            if abs(double(item.frame_rate)-double(channel.movie_info.frame_rate))>tolerance
                error('AAA:Results:ManifestFrameRateMismatch', ...
                    'Saved trace and manifest frame rates disagree for role %s.',role);
            end
        else
            channel.movie_info.frame_rate=double(item.frame_rate);
        end
    end
    channel.movie_info=copy_manifest(channel.movie_info,item, ...
        'logical_source_path','logical_path');
    channel.movie_info=copy_manifest(channel.movie_info,item, ...
        'source_path','resolved_path');
    channel.movie_info=copy_manifest(channel.movie_info,item, ...
        'source_files','physical_files');
    channel.movie_info=copy_manifest(channel.movie_info,item, ...
        'source_is_split','is_split');
    channel.movie_info=copy_manifest(channel.movie_info,item, ...
        'source_part_count','part_count');
    channel.movie_info=copy_manifest(channel.movie_info,item, ...
        'source_resolution_method','resolution_method');
    channel.movie_info=copy_manifest(channel.movie_info,item, ...
        'reuse_source_path','reuse_source_path');
    results.(char(role))=channel;
end
end

function target=copy_manifest(target,source,target_name,source_name)
if isfield(source,source_name)&&~isempty(source.(source_name))
    target.(target_name)=source.(source_name);
end
end

function [results,mapped]=adapt_legacy_roi_file(results,saved)
mapped="";
if ~isfield(saved,'rois')||~isstruct(saved.rois),return;end
results.roi=struct('rois',saved.rois,'info',struct('legacy',true));
if isfield(saved,'offset'),results.roi.info.offset_xy=saved.offset;end
if isfield(saved,'traces_voltage_raw')
    results=store_legacy_raw(results,"voltage",saved.traces_voltage_raw,saved);
end
if isfield(saved,'traces_calcium_raw')
    results=store_legacy_raw(results,"calcium",saved.traces_calcium_raw,saved);
end
if isfield(saved,'traces')
    results=store_legacy_raw(results,"voltage",saved.traces,saved);
end
mapped="roi";
end

function results=store_legacy_raw(results,role,data,saved)
if ~isfield(results,char(role)),results.(char(role))=struct();end
if ~isfield(results.(char(role)),'trace_results')
    results.(char(role)).trace_results=struct();
end
info=struct('data',data);
if isfield(saved,'frame_rate'),info.frame_rate=saved.frame_rate;end
results.(char(role)).trace_results.raw=info;
end

function sources=frame_rate_sources(results,files)
sources=struct();
for role=["voltage","calcium"]
    if ~isfield(results,char(role)),continue;end
    value=results.(char(role));rate=NaN;field="";
    if isfield(value,'movie_info')&&isfield(value.movie_info,'frame_rate')
        rate=double(value.movie_info.frame_rate);field="movie_info.frame_rate";
    elseif isfield(value,'trace_results')&&isfield(value.trace_results,'raw') ...
            && isfield(value.trace_results.raw,'frame_rate')
        rate=double(value.trace_results.raw.frame_rate);field="trace_results.raw.frame_rate";
    end
    if isscalar(rate)&&isfinite(rate)&&rate>0
        sources.(char(role))=struct('value',rate,'field',field,'files',files);
    end
end
end

function roles=expected_roles(spec)
roles=strings(0,1);
if isfield(spec,'roles'),roles=lower(string(spec.roles(:)))';
elseif isfield(spec,'mode')
    if lower(string(spec.mode))=="dual",roles=["voltage","calcium"];
    elseif lower(string(spec.mode))=="ap",roles="voltage";end
end
end

function emit_report(spec,report)
message="Result load: schema="+report.detected_schema+ ...
    ", adapter="+report.adapter_used+", source="+report.absolute_source_path+".";
if isfield(spec,'progress')&&isa(spec.progress,'function_handle')
    try,spec.progress(message);catch,end
end
if isfield(spec,'runtime')
    event=struct('type',"compatibility",'mode',get_spec(spec,'mode',""), ...
        'scope',get_spec(spec,'scope',""),'input_path',report.absolute_source_path, ...
        'output_path',"",'section',"load_results",'state',"completed", ...
        'message',message,'timestamp',datetime('now'),'details',report);
    aaa.io.emit_event(spec.runtime,event);
end
end

function value=get_spec(spec,name,default)
value=default;if isfield(spec,name),value=string(spec.(name));end
end

function report=empty_report(source,absolute_source)
report=struct('source_path',string(source), ...
    'absolute_source_path',string(absolute_source), ...
    'detected_schema',"",'adapter_used',"", ...
    'source_files',strings(0,1),'absolute_source_files',strings(0,1), ...
    'relative_source_files',strings(0,1),'frame_rate_source',struct(), ...
    'mapped_fields',strings(0,1),'missing_fields',strings(0,1), ...
    'warnings',strings(0,1),'loaded_at',NaT);
end

function output=merge_struct_recursive(base,patch)
output=base;names=fieldnames(patch);
for idx=1:numel(names)
    name=names{idx};
    if isfield(output,name)&&isstruct(output.(name))&&isscalar(output.(name)) ...
            && isstruct(patch.(name))&&isscalar(patch.(name))
        output.(name)=merge_struct_recursive(output.(name),patch.(name));
    else
        output.(name)=patch.(name);
    end
end
end

function values=relative_paths(files,root)
values=string(files(:));prefix=string(root)+filesep;
for idx=1:numel(values)
    if startsWith(lower(values(idx)),lower(prefix))
        values(idx)=extractAfter(values(idx),strlength(prefix));
    end
end
end

function root=folder_root(source)
if isfolder(source),root=source;else,root=string(fileparts(source));end
end

function ext=fileparts_extension(path)
[~,~,ext]=fileparts(path);
end

function value=file_time(path)
info=dir(path);if isempty(info),value=NaT;else,value=datetime(info.datenum,'ConvertFrom','datenum');end
end

function value=absolute_path(value)
value=strtrim(string(value));if strlength(value)==0,return;end
if ~is_absolute(value),value=fullfile(pwd,value);end
try,value=string(char(java.io.File(char(value)).getCanonicalPath()));
catch,value=string(char(java.io.File(char(value)).getAbsolutePath()));end
end

function tf=is_absolute(value)
text=char(string(value));
tf=~isempty(regexp(text,'^[A-Za-z]:[\\/]','once')) ...
    || startsWith(text,'\\')||startsWith(text,'/');
end

function value=empty_candidate()
value=struct('file',"",'variable',"",'score',-Inf,'timestamp',NaT);
end
