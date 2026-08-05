function report = inspect_result_folder(folder, mode, plan, max_mat_files)
%INSPECT_RESULT_FOLDER Validate saved results through the one compatibility loader.

if nargin<4||isempty(max_mat_files),max_mat_files=500;end
folder=strtrim(string(folder));mode=lower(strtrim(string(mode)));
report=struct('folder',folder,'valid',false,'reason',"", ...
    'result_file_count',0,'inspected_file_count',0, ...
    'trace_bindings',repmat(empty_binding(),0,1), ...
    'saved_frame_rates',struct(),'warnings',strings(0,1), ...
    'compatibility_report',struct());
if ~isfolder(folder),report.reason="not a folder";return;end
if ~ismember(mode,["ap","dual"])
    error('AAA:Reuse:InvalidMode','mode must be ap or dual.');
end
try
    [results,compatibility]=load_analysis_results(folder,struct( ...
        'mode',mode,'scope',"single",'max_mat_files',max_mat_files, ...
        'sections',["input","roi","trace"]));
catch ME
    report.reason=string(ME.message);return;
end
report.compatibility_report=compatibility;
report.result_file_count=numel(compatibility.absolute_source_files);
report.inspected_file_count=report.result_file_count;
report.warnings=compatibility.warnings;
needs_trace=reused(plan,'input')||reused(plan,'trace');
if ~needs_trace
    report.valid=~isempty(fieldnames(results));
    if report.valid,report.reason="recognized result content";else,report.reason="no recognized result content";end
    return;
end
if mode=="ap",roles="voltage";else,roles=["voltage","calcium"];end
for role=roles
    if ~isfield(results,char(role))
        report.reason="missing semantic trace role "+role;return;
    end
    [binding,reason]=validate_role(results.(char(role)),role,compatibility);
    if strlength(reason)>0,report.reason=reason;return;end
    report.trace_bindings(end+1,1)=binding; %#ok<AGROW>
    report.saved_frame_rates.(char(role))=binding.frame_rate;
end
report.valid=true;report.reason="complete semantic trace result";
end

function [binding,reason]=validate_role(value,role,compatibility)
binding=empty_binding();reason="";
if ~isstruct(value)||~isfield(value,'trace_results')
    reason="saved trace role "+role+" has no trace_results";return;
end
required=["raw","bleach_removed","baseline","noise_reference", ...
    "noise","sensitivity","snr"];
if role=="calcium"
    required=[required,"raw_smoothed","sensitivity_smoothed","snr_smoothed"];
end
stages=value.trace_results;expected=[];
for name=required
    if ~isfield(stages,char(name))||~isstruct(stages.(char(name))) ...
            ||~isfield(stages.(char(name)),'data') ...
            ||isempty(stages.(char(name)).data)||~ismatrix(stages.(char(name)).data)
        reason="saved trace role "+role+" is missing valid stage "+name;return;
    end
    current=size(stages.(char(name)).data);
    if isempty(expected),expected=current;
    elseif ~isequal(current,expected),reason="trace stage size mismatch for "+role;return;end
end
rate=saved_frame_rate(value);
if ~isfinite(rate)||rate<=0
    reason="saved trace role "+role+" has no authoritative frame rate";return;
end
file=choose_role_file(compatibility.absolute_source_files);
binding=struct('role',role,'file_path',file,'variable_name',"result", ...
    'container_kind',"unified_results",'nested_field',"", ...
    'frame_rate',rate,'frame_count',expected(1),'roi_count',expected(2), ...
    'stage_names',string(fieldnames(stages)),'reason',"complete shared trace schema");
end

function rate=saved_frame_rate(value)
values=[];
if isfield(value,'movie_info')&&isfield(value.movie_info,'frame_rate')
    values(end+1)=double(value.movie_info.frame_rate); %#ok<AGROW>
end
if isfield(value,'trace_results')&&isfield(value.trace_results,'raw') ...
        &&isfield(value.trace_results.raw,'frame_rate')
    values(end+1)=double(value.trace_results.raw.frame_rate); %#ok<AGROW>
end
values=values(isfinite(values)&values>0);
if isempty(values),rate=NaN;return;end
if max(values)-min(values)>max(1e-9,1e-9*max(values))
    error('AAA:Reuse:SavedFrameRateMismatch', ...
        'Saved movie_info and raw trace frame rates disagree.');
end
rate=values(1);
end

function file=choose_role_file(files)
files=string(files(:));file="";if isempty(files),return;end
idx=find(endsWith(lower(files),string(filesep)+"trace.mat") ...
    |endsWith(lower(files),"/trace.mat"),1);
if isempty(idx),idx=1;end
file=files(idx);
end

function tf=reused(plan,name)
tf=isstruct(plan)&&isfield(plan,name)&&string(plan.(name))=="reuse";
end

function binding=empty_binding()
binding=struct('role',"",'file_path',"",'variable_name',"", ...
    'container_kind',"",'nested_field',"",'frame_rate',NaN, ...
    'frame_count',0,'roi_count',0,'stage_names',strings(0,1),'reason',"");
end
