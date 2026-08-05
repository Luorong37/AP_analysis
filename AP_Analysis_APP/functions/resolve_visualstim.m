function [context, windows] = resolve_visualstim(source_path, channels, override, params)
%RESOLVE_VISUALSTIM Discover visual-stim metadata and build role-keyed windows.
% Frame conversion follows the frozen Dual rule: sync row 2 is the anchor,
% camera indices are rounded, and traces are never interpolated.

if nargin < 3 || isempty(override), override = struct(); end
if nargin < 4 || isempty(params), params = struct(); end
validateattributes(channels, {'struct'}, {'nonempty'});
context = empty_context();
windows = struct();

[context, windows, resolved] = apply_override(context, override);
if resolved, return; end

root = resolve_root(source_path);
[cycle_manifest, cycle_file] = load_named_struct(fullfile(root,'cycle_manifest.mat'),'manifest');
record_file = absolute_from(root,resolve_reference(cycle_manifest, 'record_manifest', fullfile(fileparts(root),'record_manifest.mat')));
[record_manifest, record_file] = load_named_struct(record_file,'manifest');
method_file = resolve_reference(cycle_manifest, 'method_manifest', "");
if strlength(method_file)>0,method_file=absolute_from(root,method_file);end
if strlength(method_file)==0, method_file=absolute_from(fileparts(record_file),resolve_reference(record_manifest,'method_manifest',"")); end
if strlength(method_file)==0
    candidates = [fullfile(root,'method_manifest.mat'), ...
        fullfile(fileparts(root),'method_manifest.mat'), ...
        fullfile(fileparts(fileparts(root)),'method_manifest.mat')];
    method_file = first_file(candidates);
end
[method_manifest, method_file] = load_named_struct(method_file,'manifest');

context.cycle_manifest_path = string(cycle_file);
context.record_manifest_path = string(record_file);
context.method_manifest_path = string(method_file);
context.recordmode = first_nested_string({cycle_manifest,record_manifest,method_manifest}, {'spec','recordmode'}, "unknown");
if has_nested(method_manifest,{'spec','stimSpec'})
    context.stimSpec = method_manifest.spec.stimSpec;
end
if has_nested(method_manifest,{'actual','stimRuntime'})
    context.stimRuntime = method_manifest.actual.stimRuntime;
end
if isfield(context.stimSpec,'selectedProgram'), context.selected_program=string(context.stimSpec.selectedProgram); end
if isfield(context.stimSpec,'selectedLabel'), context.selected_label=string(context.stimSpec.selectedLabel); end

logs_file = absolute_from(root,resolve_manifest_file(cycle_manifest,'logs_mat',fullfile(root,'logs.mat')));
[logs_value, logs_file] = load_named_struct(logs_file,'logs');
context.logs_path = string(logs_file);
context.logs = logs_value;
context=merge_override(context,override);
if ~strcmpi(context.recordmode,'visualstim')
    context.reason = "recordmode is not visualstim";
    return;
end
if isempty(fieldnames(context.stimSpec))
    context.reason = "method manifest has no stimSpec";
    return;
end
if isempty(fieldnames(context.logs)) || ~isfield(context.logs,'sync') || ~istable(context.logs.sync)
    context.reason = "visualstim source has no logs.sync table";
    return;
end

roles = strings(1,numel(channels));
for idx=1:numel(channels)
    role=lower(string(channels(idx).profile.role)); roles(idx)=role;
    if isfield(context.sync_vars,char(role)),sync_name=string(context.sync_vars.(char(role)));
    elseif isfield(context,char(role+"_sync_var")),sync_name=string(context.(char(role+"_sync_var")));
    else,sync_name=find_sync_variable(context.logs.sync,channels(idx).camera_index);end
    if strlength(sync_name)==0
        context.reason="No sync variable for camera "+channels(idx).camera_index;
        return;
    end
    context.sync_vars.(char(role))=sync_name;
end
context.stim_type=classify_stim(context.stimSpec);
windows=build_windows(context,channels,params);
context.supported=isfield(windows,'supported') && windows.supported;
if ~context.supported && strlength(context.reason)==0
    context.reason="No complete stimulus trials fit every requested channel.";
end
context.logs=struct(); % Do not duplicate potentially large logs in results.
end

function value=empty_context()
value=struct('recordmode',"unknown",'supported',false,'reason',"", ...
    'stim_type',"none",'selected_program',"",'selected_label',"", ...
    'stimSpec',struct(),'stimRuntime',struct(),'logs',struct(), ...
    'logs_path',"",'cycle_manifest_path',"",'record_manifest_path',"", ...
    'method_manifest_path',"",'sync_vars',struct());
end

function [context,windows,resolved]=apply_override(context,value)
windows=struct(); resolved=false;
if ~isstruct(value)||isempty(fieldnames(value)),return;end
if isfield(value,'stim_windows'),windows=value.stim_windows;resolved=true;end
if isfield(value,'windows'),windows=value.windows;resolved=true;end
fields=setdiff(fieldnames(value),{'stim_windows','windows'});
for idx=1:numel(fields),context.(fields{idx})=value.(fields{idx});end
if resolved
    context.supported=isstruct(windows)&&isfield(windows,'supported')&&windows.supported;
    if ~context.supported && isstruct(windows)&&~isempty(fieldnames(windows)),context.supported=true;end
    return;
end
if isfield(context,'logs')&&isfield(context,'stimSpec')&&isfield(context.logs,'sync')&&istable(context.logs.sync)
    resolved=false;
end
end
function context=merge_override(context,value)
if ~isstruct(value),return;end;fields=setdiff(fieldnames(value),{'stim_windows','windows'});for idx=1:numel(fields),context.(fields{idx})=value.(fields{idx});end
end

function root=resolve_root(source)
root=string(source);
if isfile(root),root=string(fileparts(root));end
end

function [value,file]=load_named_struct(file,varname)
value=struct(); file=string(file);
if strlength(file)==0||~isfile(file),return;end
s=load(file);
if isfield(s,varname)&&isstruct(s.(varname)),value=s.(varname);
elseif isfield(s,'manifest')&&isstruct(s.manifest),value=s.manifest;
elseif isfield(s,'logs')&&isstruct(s.logs),value=s.logs;
elseif isscalar(s)&&isstruct(s),value=s;
end
end

function file=resolve_reference(manifest,name,fallback)
file=string(fallback);
if isstruct(manifest)&&isfield(manifest,'refs')&&isstruct(manifest.refs)&&isfield(manifest.refs,name)
    file=string(manifest.refs.(name));
end
end
function file=resolve_manifest_file(manifest,name,fallback)
file=string(fallback);
if isstruct(manifest)&&isfield(manifest,'artifacts')&&isstruct(manifest.artifacts)&&isfield(manifest.artifacts,name)
    file=string(manifest.artifacts.(name));
end
end
function file=first_file(values)
file=""; for idx=1:numel(values),if isfile(values(idx)),file=string(values(idx));return;end,end
end
function file=absolute_from(base,file)
file=string(file);if strlength(file)==0||isfile(file)||isfolder(file),return;end
if ~isempty(regexp(char(file),'^[A-Za-z]:[\\/]|^\\\\','once')),return;end
file=string(fullfile(base,file));
end
function tf=has_nested(s,path)
tf=isstruct(s); for idx=1:numel(path),if ~tf||~isfield(s,path{idx}),tf=false;return;end;s=s.(path{idx});end
end
function value=first_nested_string(values,path,fallback)
value=string(fallback); for idx=1:numel(values),s=values{idx};if has_nested(s,path),for j=1:numel(path),s=s.(path{j});end;value=string(s);return;end,end
end

function name=find_sync_variable(sync,camera_index)
names=string(sync.Properties.VariableNames); low=lower(names);
mask=contains(low,'camera') & contains(low,string(camera_index));
hit=find(mask,1); if isempty(hit),fallback=3+double(camera_index);if fallback<=numel(names),hit=fallback;end,end
if isempty(hit),name="";else,name=names(hit);end
end

function kind=classify_stim(spec)
text=lower(join([string(fieldnames(spec)); string(get_field(spec,'selectedProgram',"")); string(get_field(spec,'selectedLabel',""))]," "));
if contains(text,'random')&&contains(text,'grating'),kind="visualstim_random_grating";
elseif contains(text,'grating')||isfield(spec,'orientations'),kind="visualstim_grating";
elseif contains(text,'flash'),kind="visualstim_flash";
elseif contains(text,'blue'),kind="visualstim_blue";
elseif contains(text,'luminance')||contains(text,'gray')||contains(text,'white')||contains(text,'black'),kind="visualstim_luminance";
elseif contains(text,'flicker')||contains(text,'contrast')||isfield(spec,'flicker')||isfield(spec,'contrastReverse'),kind="visualstim_flicker";
else,kind="visualstim_generic";end
end

function windows=build_windows(context,channels,params)
switch context.stim_type
    case "visualstim_random_grating", windows=grating_windows(context,channels,true);
    case "visualstim_grating", windows=grating_windows(context,channels,false);
    case "visualstim_flash", windows=flash_windows(context,channels,params);
    case {"visualstim_blue","visualstim_luminance"}, windows=block_windows(context,channels);
    case "visualstim_flicker", windows=flicker_windows(context,channels);
    otherwise
        if isfield(context.stimSpec,'blockSequence'),windows=block_windows(context,channels);else,windows=grating_windows(context,channels,false);end
end
end

function windows=grating_windows(context,channels,is_random)
sync=context.logs.sync; ifi=resolve_ifi(context,sync);
base_n=max(1,round(double(get_field(context.stimSpec,'isi',0))/ifi));
stim_n=max(1,round(double(get_field(context.stimSpec,'duration',ifi))/ifi));
if is_random
    if ~isfield(context.logs,'stimulus')||~isfield(context.logs.stimulus,'angleSequence'),error('AAA:Stim:MissingAngleSequence','Random grating requires logs.stimulus.angleSequence.');end
    orientation=double(context.logs.stimulus.angleSequence(:));
else
    orientation=double(get_field(context.stimSpec,'orientations',0)); orientation=orientation(:);
end
available=floor((height(sync)-1)/(base_n+stim_n)); n=min(numel(orientation),available);
if is_random && n<numel(orientation),error('AAA:Stim:IncompleteRandomGrating','Sync table has fewer complete trials than angleSequence.');end
orientation=orientation(1:n); rows=zeros(n,2); stim=zeros(n,2);
for i=1:n,b=2+(i-1)*(base_n+stim_n);rows(i,:)=[b,b+base_n-1];stim(i,:)=[b+base_n,b+base_n+stim_n-1];end
[labels,index]=condition_labels(string(orientation));
windows=finalize(context,channels,rows,stim,orientation,index,labels);
windows.is_grating=true; windows.angle_sequence=orientation';
if is_random
    counts=groupcounts(categorical(orientation));
    if numel(unique(counts))>1,error('AAA:Stim:UnbalancedRandomGrating','Random grating repeats per direction are unequal.');end
end
end

function windows=flash_windows(context,channels,params)
[entries,labels,colors,durations]=sequence(context.stimSpec);
if isempty(entries),windows=empty_windows(context);return;end
ifi=resolve_ifi(context,context.logs.sync); starts=2+cumsum([0;round(durations(1:end-1)/ifi)]);
flash=find(~baseline_like(labels)); base_n=max(1,round(double(get_field(params,'flash_baseline_window_s',2))/ifi));
resp_n=max(1,round(double(get_field(params,'flash_response_window_s',2))/ifi));
b=zeros(numel(flash),2);s=b;trial_labels=strings(numel(flash),1);
for i=1:numel(flash),on=starts(flash(i));b(i,:)=[max(2,on-base_n),on-1];s(i,:)=[on,on+resp_n-1];trial_labels(i)=labels(flash(i));end
[condition,index]=condition_labels(trial_labels);
windows=finalize(context,channels,b,s,zeros(numel(flash),1),index,condition);
windows.condition_colors=condition_colors(condition,labels,colors); windows.shading_labels=trial_labels;
end

function windows=block_windows(context,channels)
[~,labels,colors,durations]=sequence(context.stimSpec); if isempty(labels),windows=empty_windows(context);return;end
ifi=resolve_ifi(context,context.logs.sync); starts=2+cumsum([0;round(durations(1:end-1)/ifi)]); ends=starts+round(durations/ifi)-1;
events=find(~baseline_like(labels));b=zeros(numel(events),2);s=b;trial_labels=strings(numel(events),1);
for i=1:numel(events),e=events(i);base=max(1,e-1);b(i,:)=[starts(base),ends(base)];s(i,:)=[starts(e),ends(e)];trial_labels(i)=labels(e);end
[condition,index]=condition_labels(trial_labels);windows=finalize(context,channels,b,s,zeros(numel(events),1),index,condition);
windows.block_labels=labels;windows.condition_colors=condition_colors(condition,labels,colors);
end

function windows=flicker_windows(context,channels)
ifi=resolve_ifi(context,context.logs.sync); f=get_field(context.stimSpec,'flicker',get_field(context.stimSpec,'contrastReverse',struct()));
hz=double(get_field(f,'frequencyHz',1));duration=double(get_field(f,'duration',0));half=max(1,round((1/(2*hz))/ifi));count=max(2,round(duration*hz*2));
starts=2+(0:count-1)'*half; ranges=[starts,starts+half-1];b=ranges(1:end-1,:);s=ranges(2:end,:);
labels=repmat(["phase_1";"phase_2"],ceil(size(s,1)/2),1);labels=labels(1:size(s,1));[condition,index]=condition_labels(labels);
windows=finalize(context,channels,b,s,zeros(size(s,1),1),index,condition);windows.block_labels=labels;
end

function windows=finalize(context,channels,base_rows,stim_rows,orientations,condition_index,condition_labels_value)
windows=empty_windows(context); n=size(stim_rows,1); valid=true(n,1);
for idx=1:numel(channels)
    role=char(lower(string(channels(idx).profile.role))); series=double(context.logs.sync.(char(context.sync_vars.(role))));
    offset=series(2)-1; nf=frame_count(channels(idx)); bf=round([series(base_rows(:,1))-offset,series(base_rows(:,2))-offset]);sf=round([series(stim_rows(:,1))-offset,series(stim_rows(:,2))-offset]);
    bf(:,1)=max(1,bf(:,1));sf(:,2)=min(nf,sf(:,2));valid=valid&all(isfinite([bf,sf]),2)&bf(:,1)<=bf(:,2)&sf(:,1)<=sf(:,2)&bf(:,2)<=nf&sf(:,1)>=1;
    windows.(role)=struct('baseline_frames',bf,'stim_frames',sf);
end
roles=fieldnames(context.sync_vars);for idx=1:numel(roles),windows.(roles{idx}).baseline_frames=windows.(roles{idx}).baseline_frames(valid,:);windows.(roles{idx}).stim_frames=windows.(roles{idx}).stim_frames(valid,:);end
windows.supported=any(valid);windows.trial_count=sum(valid);windows.condition_index=condition_index(valid);windows.condition_labels=condition_labels_value;windows.trial_labels=condition_labels_value(condition_index(valid));windows.orientations=orientations(valid);windows.unique_orientations=unique(orientations(valid));windows.is_grating=false;
for idx=1:numel(channels),role=char(lower(string(channels(idx).profile.role)));fps=double(channels(idx).profile.frame_rate);windows.(role).baseline_time_ranges=windows.(role).baseline_frames/fps;windows.(role).stim_time_ranges=windows.(role).stim_frames/fps;end
end

function n=frame_count(channel)
if isfield(channel,'data')&&isfield(channel.data,'movie_info')&&isfield(channel.data.movie_info,'frame_count'),n=double(channel.data.movie_info.frame_count);
elseif isfield(channel,'results')&&isfield(channel.results,'movie_info')&&isfield(channel.results.movie_info,'frame_count'),n=double(channel.results.movie_info.frame_count);
else,[x,~]=resolve_trace_stage(channel.results,{'sensitivity','snr','bleach_removed','raw'});n=size(x,1);end
end
function value=empty_windows(context)
value=struct('supported',false,'stim_type',context.stim_type,'selected_program',context.selected_program,'trial_count',0,'condition_index',zeros(0,1),'condition_labels',strings(0,1),'trial_labels',strings(0,1),'condition_colors',zeros(0,3),'orientations',zeros(0,1),'unique_orientations',zeros(0,1),'is_grating',false);
end
function ifi=resolve_ifi(context,sync)
ifi=double(get_field(context.stimRuntime,'ifi',NaN));if ~isfinite(ifi)||ifi<=0
    names=string(sync.Properties.VariableNames);hit=find(strcmpi(names,'PTB_VBL_Time'),1);if isempty(hit),error('AAA:Stim:MissingIFI','Cannot resolve stimulus IFI.');end
    ifi=median(diff(double(sync.(char(names(hit))))),'omitnan');
end
end
function [entries,labels,colors,durations]=sequence(spec)
entries=[];labels=strings(0,1);colors=zeros(0,3);durations=zeros(0,1);if ~isfield(spec,'blockSequence'),return;end
entries=spec.blockSequence;if iscell(entries),entries=[entries{:}];end
n=numel(entries);labels=strings(n,1);colors=nan(n,3);durations=zeros(n,1);
for i=1:n,labels(i)=string(get_field(entries(i),'label',get_field(entries(i),'name',"block_"+i)));durations(i)=double(get_field(entries(i),'duration',get_field(entries(i),'durationSeconds',0)));c=get_field(entries(i),'color',[NaN NaN NaN]);if isnumeric(c)&&numel(c)>=3,colors(i,:)=double(c(1:3));end,end
repeat=max(1,round(double(get_field(spec,'repeatCount',1))));labels=repmat(labels,repeat,1);colors=repmat(colors,repeat,1);durations=repmat(durations,repeat,1);
end
function tf=baseline_like(labels)
low=lower(string(labels));tf=contains(low,'gray')|contains(low,'baseline')|contains(low,'rest')|contains(low,'isi');
end
function [labels,index]=condition_labels(values)
[labels,~,index]=unique(string(values(:)),'stable');
end
function colors=condition_colors(condition,labels,allcolors)
colors=lines(max(1,numel(condition)));for i=1:numel(condition),hit=find(labels==condition(i)&all(isfinite(allcolors),2),1);if ~isempty(hit),colors(i,:)=allcolors(hit,:);end,end
end
function value=get_field(s,name,fallback)
if isstruct(s)&&isfield(s,name)&&~isempty(s.(name)),value=s.(name);else,value=fallback;end
end
