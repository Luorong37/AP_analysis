function output_files = render_record_grating_trials(entries,grating,output_dir,spec,progress)
%RENDER_RECORD_GRATING_TRIALS Render saved random-grating Trial/Cycle views.
% Trial traces come directly from each saved sensitivity stage. Baseline and
% stimulus frame windows are concatenated without interpolation and aligned
% to stimulus onset. This renderer does not recalculate traces or peaks.

if nargin<4||~isstruct(spec),spec=struct();end
if nargin<5,progress=[];end
output_dir=string(output_dir);
output_files=struct('status',"not_applicable",'reason',"", ...
    'stack',struct(),'heatmap',struct(),'cycle_trace',struct(), ...
    'cycle_metrics',struct(),'files',strings(0,1));
if isempty(entries)
    output_files.reason="No Cycle entries were supplied.";return;
end
random_flags=arrayfun(@is_random_grating_entry,entries);
if ~any(random_flags)
    output_files.reason="Record is not identified as random grating.";return;
end
if ~all(random_flags)
    error('AAA:Record:MixedGratingOrder', ...
        'Random and non-random grating Cycles cannot share Trial/Cycle plots.');
end
roles=string(get_field(spec,'roles',fieldnames(entries(1).channels)));
roles=lower(roles(:));
angles=resolve_angles(entries);
roi_counts=validate_entries(entries,roles);
if numel(unique(struct2array(roi_counts)))~=1
    error('AAA:Record:RoleROIMismatch', ...
        'Combined Cycle trace figures require equal ROI counts across roles.');
end

stack_root=fullfile(output_dir,'9_rg_stack');
heatmap_root=fullfile(output_dir,'10_rg_hm');
cycle_root=fullfile(output_dir,'11_rg_trace');
metric_root=fullfile(output_dir,'12_rg_cycle_metrics');
make_folder(stack_root);make_folder(heatmap_root);make_folder(cycle_root);

for role=roles'
    role_name=char(role);nrois=roi_counts.(role_name);
    report(progress,"Rendering "+role+" direction-grouped Trial stacks.");
    stack_files=render_role_stacks(entries,role,angles,nrois,stack_root,spec,progress);
    output_files.stack.(role_name)=stack_files;
    output_files.files=[output_files.files;stack_files.files]; %#ok<AGROW>

    report(progress,"Rendering "+role+" direction-grouped Trial heatmaps.");
    heatmap_files=render_role_heatmaps(entries,role,angles,nrois,heatmap_root,spec,progress);
    output_files.heatmap.(role_name)=heatmap_files;
    output_files.files=[output_files.files;heatmap_files.files]; %#ok<AGROW>
end

report(progress,"Rendering per-Cycle continuous traces with direction labels.");
cycle_files=render_cycle_traces(entries,roles,angles,roi_counts.(char(roles(1))), ...
    cycle_root,spec,progress);
output_files.cycle_trace=cycle_files;
output_files.files=[output_files.files;cycle_files.files];

if isstruct(grating)&&isfield(grating,'status')&&string(grating.status)=="completed"
    make_folder(metric_root);
    report(progress,"Rendering per-Cycle DSI/OSI comparisons.");
    metric_files=render_cycle_metrics(grating,roles,metric_root,progress);
    output_files.cycle_metrics=metric_files;
    output_files.files=[output_files.files;metric_files.files];
end
output_files.status="completed";
output_files.reason="";
output_files.angles=angles;
output_files.roles=roles;
output_files.files=unique(output_files.files(strlength(output_files.files)>0),'stable');
end

function tf=is_random_grating_entry(entry)
tf=false;
if isfield(entry,'stim_windows')&&isstruct(entry.stim_windows) ...
        &&isfield(entry.stim_windows,'stim_type')
    kind=lower(string(entry.stim_windows.stim_type));
    tf=contains(kind,"random")&&contains(kind,"grating");
end
if ~tf&&isfield(entry,'stim_analysis_kind')
    tf=lower(string(entry.stim_analysis_kind))=="random_grating_tuning";
end
end

function angles=resolve_angles(entries)
values=cell(numel(entries),1);
for idx=1:numel(entries)
    windows=entries(idx).stim_windows;
    if ~isfield(windows,'orientations')||isempty(windows.orientations)
        error('AAA:Record:MissingOrientations', ...
            'Cycle %s has no saved grating orientations.',entries(idx).cycle_name);
    end
    values{idx}=mod(double(windows.orientations(:)),360);
end
angles=sort(unique(vertcat(values{:})))';
end

function counts=validate_entries(entries,roles)
counts=struct();
for role=roles'
    role_name=char(role);count=NaN;
    for idx=1:numel(entries)
        if ~isfield(entries(idx).channels,role_name) ...
                ||~isfield(entries(idx).channels.(role_name),'sensitivity')
            error('AAA:Record:MissingSensitivity', ...
                'Cycle %s is missing %s sensitivity.',entries(idx).cycle_name,role);
        end
        stage=entries(idx).channels.(role_name).sensitivity;
        current=size(stage.data,2);
        if ~isfinite(count),count=current;
        elseif current~=count
            error('AAA:Record:ROIMismatch', ...
                '%s ROI count differs in Cycle %s.',role,entries(idx).cycle_name);
        end
        validate_trial_windows(entries(idx),role,size(stage.data,1));
    end
    counts.(role_name)=count;
end
end

function validate_trial_windows(entry,role,nframes)
role_name=char(role);windows=entry.stim_windows;
if ~isfield(windows,role_name)||~isfield(windows.(role_name),'baseline_frames') ...
        ||~isfield(windows.(role_name),'stim_frames')
    error('AAA:Record:MissingTrialWindows', ...
        'Cycle %s is missing %s baseline/stimulus frames.',entry.cycle_name,role);
end
orientations=windows.orientations(:);
baseline=double(windows.(role_name).baseline_frames);
stimulus=double(windows.(role_name).stim_frames);
if size(baseline,2)~=2||~isequal(size(baseline),size(stimulus)) ...
        ||size(baseline,1)~=numel(orientations) ...
        ||any(~isfinite([baseline(:);stimulus(:)])) ...
        ||any([baseline(:);stimulus(:)]~=round([baseline(:);stimulus(:)])) ...
        ||any([baseline(:);stimulus(:)]<1) ...
        ||any([baseline(:);stimulus(:)]>nframes) ...
        ||any(baseline(:,1)>baseline(:,2))||any(stimulus(:,1)>stimulus(:,2)) ...
        ||any(stimulus(:,1)<baseline(:,1))
    error('AAA:Record:InvalidTrialWindows', ...
        'Cycle %s has invalid %s random-grating frame windows.',entry.cycle_name,role);
end
end

function result=render_role_stacks(entries,role,angles,nrois,root,spec,progress)
role_name=char(role);result=struct('overview',struct(),'roi_dir',"",'files',strings(0,1));
colors=hsv(max(1,numel(angles)));
fig=figure('Color','w','Visible','off','Position',grid_position(numel(angles),nrois));
layout=tiledlayout(fig,nrois,numel(angles),'Padding','compact','TileSpacing','compact');
for roi=1:nrois
    trials=collect_trials(entries,role,roi,angles,spec);
    for angle_idx=1:numel(angles)
        ax=nexttile(layout);subset=trials(abs([trials.angle]-angles(angle_idx))<1e-9);
        title_text="";if roi==1,title_text=sprintf('%g deg',angles(angle_idx));end
        plot_trial_stack(ax,subset,colors(angle_idx,:),title_text);
        if angle_idx==1,ylabel(ax,sprintf('ROI %03d',roi));end
    end
end
sgtitle(fig,upper(role)+" sensitivity | random-grating Trials | columns: direction, rows: ROI");
overview=save_bundle(fig,root,role+"_stack");close(fig);
result.overview=overview;result.files=[result.files;struct_files(overview)];

roi_dir=fullfile(root,'roi',role_name);make_folder(roi_dir);result.roi_dir=string(roi_dir);
for roi=1:nrois
    report(progress,sprintf('Trial stacks: %s ROI %d/%d.',role,roi,nrois));
    trials=collect_trials(entries,role,roi,angles,spec);
    fig=figure('Color','w','Visible','off','Position',grid_position(numel(angles),1));
    layout=tiledlayout(fig,1,numel(angles),'Padding','compact','TileSpacing','compact');
    for angle_idx=1:numel(angles)
        ax=nexttile(layout);subset=trials(abs([trials.angle]-angles(angle_idx))<1e-9);
        plot_trial_stack(ax,subset,colors(angle_idx,:),sprintf('%g deg',angles(angle_idx)));
        if angle_idx==1,ylabel(ax,'Trials');end
    end
    sgtitle(fig,sprintf('%s sensitivity | ROI %03d | columns: direction',upper(role),roi));
    bundle=save_bundle(fig,roi_dir,string(sprintf('r%03d_stack',roi)));close(fig);
    result.files=[result.files;struct_files(bundle)]; %#ok<AGROW>
end
end

function result=render_role_heatmaps(entries,role,angles,nrois,root,spec,progress)
role_name=char(role);result=struct('overview',struct(),'roi_dir',"", ...
    'color_limit',[],'files',strings(0,1));
limit=resolve_color_limit(entries,role,angles,nrois,spec);
fig=figure('Color','w','Visible','off','Position',grid_position(numel(angles),nrois));
layout=tiledlayout(fig,nrois,numel(angles),'Padding','compact','TileSpacing','compact');last_ax=[];
for roi=1:nrois
    trials=collect_trials(entries,role,roi,angles,spec);
    for angle_idx=1:numel(angles)
        last_ax=nexttile(layout);subset=trials(abs([trials.angle]-angles(angle_idx))<1e-9);
        title_text="";if roi==1,title_text=sprintf('%g deg',angles(angle_idx));end
        plot_trial_heatmap(last_ax,subset,limit,title_text);
        if angle_idx==1,ylabel(last_ax,sprintf('ROI %03d',roi));end
    end
end
colormap(fig,redblue_map(256));add_shared_colorbar(last_ax,role);
sgtitle(fig,upper(role)+" sensitivity heatmap | columns: direction, rows: ROI");
overview=save_bundle(fig,root,role+"_heatmap");close(fig);
result.overview=overview;result.color_limit=limit;
result.files=[result.files;struct_files(overview)];

roi_dir=fullfile(root,'roi',role_name);make_folder(roi_dir);result.roi_dir=string(roi_dir);
for roi=1:nrois
    report(progress,sprintf('Trial heatmaps: %s ROI %d/%d.',role,roi,nrois));
    trials=collect_trials(entries,role,roi,angles,spec);
    fig=figure('Color','w','Visible','off','Position',grid_position(numel(angles),1));
    layout=tiledlayout(fig,1,numel(angles),'Padding','compact','TileSpacing','compact');last_ax=[];
    for angle_idx=1:numel(angles)
        last_ax=nexttile(layout);subset=trials(abs([trials.angle]-angles(angle_idx))<1e-9);
        plot_trial_heatmap(last_ax,subset,limit,sprintf('%g deg',angles(angle_idx)));
        if angle_idx==1,ylabel(last_ax,'Trials');end
    end
    colormap(fig,redblue_map(256));add_shared_colorbar(last_ax,role);
    sgtitle(fig,sprintf('%s sensitivity heatmap | ROI %03d',upper(role),roi));
    bundle=save_bundle(fig,roi_dir,string(sprintf('r%03d_heatmap',roi)));close(fig);
    result.files=[result.files;struct_files(bundle)]; %#ok<AGROW>
end
end

function result=render_cycle_traces(entries,roles,angles,nrois,root,spec,progress)
result=struct('root',string(root),'files',strings(0,1));
colors=hsv(max(1,numel(angles)));
for roi=1:nrois
    roi_dir=fullfile(root,sprintf('r%03d',roi));make_folder(roi_dir);
    for cycle_idx=1:numel(entries)
        report(progress,sprintf('Continuous traces: ROI %d/%d, Cycle %d/%d.', ...
            roi,nrois,cycle_idx,numel(entries)));
        fig=figure('Color','w','Visible','off','Position',[80 80 1500 max(520,330*numel(roles))]);
        layout=tiledlayout(fig,numel(roles),1,'Padding','compact','TileSpacing','compact');
        for role_idx=1:numel(roles)
            role=roles(role_idx);role_name=char(role);entry=entries(cycle_idx);
            stage=entry.channels.(role_name).sensitivity;
            polarity=resolve_polarity(entry,role,size(stage.data,2),spec);
            trace=double(stage.data(:,roi))*polarity(roi);
            ax=nexttile(layout);plot_continuous_trace(ax,stage.time,trace, ...
                entry.stim_windows,role,angles,colors);
            title(ax,sprintf('%s | %s | ROI %03d',upper(role),entry.cycle_name,roi), ...
                'Interpreter','none');ylabel(ax,'Sensitivity');
            if role_idx==numel(roles),xlabel(ax,'Time in Cycle (s)');end
        end
        sgtitle(fig,sprintf('Random grating | ROI %03d | %s',roi,entries(cycle_idx).cycle_name), ...
            'Interpreter','none');
        stem=safe_name(entries(cycle_idx).cycle_name);
        bundle=save_bundle(fig,roi_dir,stem);close(fig);
        result.files=[result.files;struct_files(bundle)]; %#ok<AGROW>
    end
end
end

function result=render_cycle_metrics(grating,roles,root,progress)
result=struct('variants',struct(),'files',strings(0,1));
variants=metric_variants(grating,roles);
for idx=1:numel(variants)
    field=variants(idx);value=grating.(char(field));target=fullfile(root,field);make_folder(target);
    report(progress,"Cycle metrics: "+field+".");
    files=strings(0,1);cycles=string(value.included_cycles(:));
    dsi=double(value.metrics_per_cycle.dsi);osi=double(value.metrics_per_cycle.osi);
    fig=figure('Color','w','Visible','off','Position',[100 80 1450 720]);
    layout=tiledlayout(fig,1,2,'Padding','compact','TileSpacing','compact');
    metric_heatmap(nexttile(layout),dsi,cycles,'DSI');
    metric_heatmap(nexttile(layout),osi,cycles,'OSI');
    sgtitle(fig,field+" | per-Cycle metrics",'Interpreter','none');
    bundle=save_bundle(fig,target,"cycle_metrics");close(fig);files=[files;struct_files(bundle)];
    roi_dir=fullfile(target,'roi');make_folder(roi_dir);
    for roi=1:size(dsi,1)
        fig=figure('Color','w','Visible','off','Position',[120 100 1200 650]);
        layout=tiledlayout(fig,2,1,'Padding','compact','TileSpacing','compact');
        metric_line(nexttile(layout),dsi(roi,:),cycles,'DSI');
        metric_line(nexttile(layout),osi(roi,:),cycles,'OSI');
        sgtitle(fig,sprintf('%s | ROI %03d',field,roi),'Interpreter','none');
        roi_bundle=save_bundle(fig,roi_dir,string(sprintf('r%03d_metrics',roi)));close(fig);
        files=[files;struct_files(roi_bundle)]; %#ok<AGROW>
    end
    result.variants.(char(field))=struct('root',string(target),'files',files);
    result.files=[result.files;files]; %#ok<AGROW>
end
end

function variants=metric_variants(grating,roles)
variants=strings(0,1);
if ismember("voltage",roles)
    if completed(grating,'voltage_peak_count'),variants(end+1,1)="voltage_peak_count";end
    if completed(grating,'voltage_sensitivity'),variants(end+1,1)="voltage_sensitivity";end
    if isempty(variants)&&completed(grating,'voltage'),variants(end+1,1)="voltage";end
end
if ismember("calcium",roles)&&completed(grating,'calcium'),variants(end+1,1)="calcium";end
end

function trials=collect_trials(entries,role,roi,angles,spec)
template=struct('trace',[],'time',[],'frame_rate',NaN,'angle',NaN, ...
    'cycle_name',"",'trial_index',NaN);
trials=repmat(template,0,1);role_name=char(role);
for angle_idx=1:numel(angles)
    target=angles(angle_idx);
    for cycle_idx=1:numel(entries)
        entry=entries(cycle_idx);windows=entry.stim_windows;
        orientations=mod(double(windows.orientations(:)),360);
        baseline=round(double(windows.(role_name).baseline_frames));
        stimulus=round(double(windows.(role_name).stim_frames));
        stage=entry.channels.(role_name).sensitivity;
        polarity=resolve_polarity(entry,role,size(stage.data,2),spec);
        matching=find(abs(orientations-target)<1e-9);
        for trial_idx=matching(:)'
            frame_index=(baseline(trial_idx,1):stimulus(trial_idx,2))';
            item=template;item.trace=double(stage.data(frame_index,roi))*polarity(roi);
            item.time=double(frame_index-stimulus(trial_idx,1))/double(stage.frame_rate);
            item.frame_rate=double(stage.frame_rate);item.angle=target;
            item.cycle_name=entry.cycle_name;item.trial_index=trial_idx;
            trials(end+1,1)=item; %#ok<AGROW>
        end
    end
end
end

function polarity=resolve_polarity(entry,role,nrois,spec)
role_name=char(role);polarity=[];
if isfield(entry,'stim_metrics')&&isstruct(entry.stim_metrics) ...
        &&isfield(entry.stim_metrics,'sensitivity') ...
        &&isfield(entry.stim_metrics.sensitivity,role_name) ...
        &&isfield(entry.stim_metrics.sensitivity.(role_name),'polarity')
    polarity=entry.stim_metrics.sensitivity.(role_name).polarity;
end
if isempty(polarity)
    name=role_name+"_polarity";
    if isfield(spec,name),polarity=spec.(name);end
end
polarity=double(polarity(:)');
if isscalar(polarity),polarity=repmat(polarity,1,nrois);end
if numel(polarity)~=nrois||any(~isfinite(polarity))||any(polarity==0)
    error('AAA:Record:MissingPolarity', ...
        'Saved or Record-control polarity is invalid for %s in Cycle %s.', ...
        role,entry.cycle_name);
end
polarity(polarity>0)=1;polarity(polarity<0)=-1;
end

function plot_trial_stack(ax,trials,color,title_text)
if isempty(trials),unavailable(ax);title(ax,title_text);return;end
hold(ax,'on');spacing=stack_spacing({trials.trace});n=numel(trials);
offsets=(n-(1:n))*spacing;
labels=strings(n,1);
for idx=1:n
    plot(ax,trials(idx).time,trials(idx).trace+offsets(idx),'Color',color,'LineWidth',.8);
    labels(idx)=sprintf('%s | T%d',trials(idx).cycle_name,trials(idx).trial_index);
end
xline(ax,0,'--k','Stim onset','LineWidth',1);
all_time=vertcat(trials.time);all_trace=vertcat(trials.trace);
xlim(ax,[min(all_time) max(all_time)]);ylim(ax,[min(all_trace)-.5*spacing max(all_trace)+offsets(1)+.5*spacing]);
[ticks,order]=sort(offsets,'ascend');set(ax,'YTick',ticks,'YTickLabel',labels(order), ...
    'TickLabelInterpreter','none','TickDir','out');xlabel(ax,'Time from onset (s)');title(ax,title_text);box(ax,'off');
end

function plot_trial_heatmap(ax,trials,limit,title_text)
if isempty(trials),unavailable(ax);title(ax,title_text);return;end
[matrix,time,labels]=trial_matrix(trials);image=imagesc(ax,time,1:size(matrix,1),matrix);
set(image,'AlphaData',isfinite(matrix));set(ax,'Color',[.88 .88 .88], ...
    'YDir','normal','YTick',1:numel(labels),'YTickLabel',labels, ...
    'TickLabelInterpreter','none','TickDir','out');clim(ax,limit);hold(ax,'on');
xline(ax,0,'--k','Stim onset','LineWidth',1);xlabel(ax,'Time from onset (s)');title(ax,title_text);box(ax,'off');
end

function [matrix,time,labels]=trial_matrix(trials)
rates=[trials.frame_rate];reference=rates(1);tolerance=max(1e-9,abs(reference)*1e-9);
if any(abs(rates-reference)>tolerance)
    error('AAA:Record:TrialFrameRateMismatch', ...
        'Random-grating Trials have inconsistent saved frame rates within one role.');
end
relative=arrayfun(@(x)round(double(x.time(:))*reference),trials,'UniformOutput',false);
frame_min=min(cellfun(@min,relative));frame_max=max(cellfun(@max,relative));frames=frame_min:frame_max;
matrix=NaN(numel(trials),numel(frames));labels=strings(numel(trials),1);
for idx=1:numel(trials)
    columns=relative{idx}-frame_min+1;matrix(idx,columns)=double(trials(idx).trace(:));
    labels(idx)=sprintf('%s | T%d',trials(idx).cycle_name,trials(idx).trial_index);
end
time=double(frames)/reference;
end

function limit=resolve_color_limit(entries,role,angles,nrois,spec)
maximum=0;
for roi=1:nrois
    trials=collect_trials(entries,role,roi,angles,spec);
    for idx=1:numel(trials)
        value=max(abs(double(trials(idx).trace)),[],'omitnan');
        if isfinite(value),maximum=max(maximum,value);end
    end
end
if maximum<=0||~isfinite(maximum),maximum=1;end
limit=[-maximum maximum];
end

function plot_continuous_trace(ax,time,trace,windows,role,angles,colors)
time=double(time(:));trace=double(trace(:));hold(ax,'on');
finite=trace(isfinite(trace));if isempty(finite),limits=[-1 1];else,limits=[min(finite) max(finite)];if limits(2)<=limits(1),limits=limits+[-.5 .5];end,end
ylim(ax,limits);if ~isempty(time),xlim(ax,[min(time) max(time)]);end
role_windows=windows.(char(role));stim=round(double(role_windows.stim_frames));orientation=mod(double(windows.orientations(:)),360);
for idx=1:size(stim,1)
    color_idx=find(abs(angles-orientation(idx))<1e-9,1);if isempty(color_idx),color_idx=1;end
    start=time(stim(idx,1));stop=time(stim(idx,2));shade=1-.30*(1-colors(color_idx,:));
    patch(ax,[start stop stop start],[limits(1) limits(1) limits(2) limits(2)],shade, ...
        'FaceAlpha',.18,'EdgeColor','none');
    text(ax,mean([start stop]),limits(2),sprintf('%g deg',orientation(idx)), ...
        'HorizontalAlignment','center','VerticalAlignment','top','FontSize',7,'Clipping','on');
end
plot(ax,time,trace,'Color',[.12 .12 .12],'LineWidth',.9);set(ax,'TickDir','out');box(ax,'off');
end

function metric_heatmap(ax,values,cycles,label)
imagesc(ax,1:numel(cycles),1:size(values,1),values);set(ax,'YDir','normal', ...
    'XTick',1:numel(cycles),'XTickLabel',cycles,'TickLabelInterpreter','none');
xtickangle(ax,35);xlabel(ax,'Cycle');ylabel(ax,'ROI');title(ax,label);colorbar(ax);
end
function metric_line(ax,values,cycles,label)
plot(ax,1:numel(cycles),values,'-o','LineWidth',1.4);set(ax,'XTick',1:numel(cycles), ...
    'XTickLabel',cycles,'TickLabelInterpreter','none');xtickangle(ax,35);ylabel(ax,label);box(ax,'off');
end
function tf=completed(s,name)
tf=isfield(s,char(name))&&isstruct(s.(char(name))) ...
    &&isfield(s.(char(name)),'status')&&string(s.(char(name)).status)=="completed";
end
function spacing=stack_spacing(cells)
values=cellfun(@(x)std(double(x),0,'omitnan'),cells);spacing=6*median(values,'omitnan');
if ~isfinite(spacing)||spacing<=0
    ranges=cellfun(@(x)range(double(x)),cells);spacing=max(1,median(ranges,'omitnan'));
end
end
function position=grid_position(nangles,nrows)
position=[40 40 min(2600,max(1400,380*nangles)) min(2800,max(850,300*nrows))];
end
function add_shared_colorbar(ax,role)
if isempty(ax)||~isvalid(ax),return;end;cb=colorbar(ax);cb.Layout.Tile='east';cb.Label.String=upper(role)+" sensitivity";
end
function map=redblue_map(n)
half=floor(n/2);map=[linspace(0,1,half)' linspace(0,1,half)' ones(half,1); ...
    ones(n-half,1) linspace(1,0,n-half)' linspace(1,0,n-half)'];
end
function unavailable(ax)
text(ax,.5,.5,'No valid Trials','Units','normalized','HorizontalAlignment','center');axis(ax,'off');
end
function files=save_bundle(fig,folder,stem)
stem=string(stem);files=struct('fig_file',string(fullfile(folder,stem+".fig")), ...
    'png_file',string(fullfile(folder,stem+".png")));
savefig(fig,files.fig_file);exportgraphics(fig,files.png_file,'Resolution',150);
end
function files=struct_files(value)
names=fieldnames(value);files=strings(0,1);
for idx=1:numel(names),candidate=string(value.(names{idx}));if isscalar(candidate)&&strlength(candidate)>0,files(end+1,1)=candidate;end,end
end
function make_folder(folder),if ~isfolder(folder),mkdir(folder);end,end
function name=safe_name(value)
name=string(regexprep(char(string(value)),'[<>:"/\\|?*]','_'));name=regexprep(name,'\s+','_');if strlength(name)==0,name="unnamed";end
end
function report(progress,message)
if ~isa(progress,'function_handle'),return;end;try,progress(string(message));catch,end
end
function value=get_field(s,name,fallback),if isfield(s,name),value=s.(name);else,value=fallback;end,end
