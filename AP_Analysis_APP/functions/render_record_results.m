function output_files = render_record_results(record_average,grating,entries,output_dir,spec,progress)
%RENDER_RECORD_RESULTS Render role-aware Record summaries for AP or Dual.
% Dual-compatible names are retained where one plot has a direct analogue;
% AP uses voltage-only layouts and never fabricates a calcium channel.

if nargin<5||~isstruct(spec),spec=struct();end
if nargin<6,progress=[];end
output_dir=string(output_dir);if ~isfolder(output_dir),mkdir(output_dir);end
roles=string(record_average.info.roles(:));mode=lower(string(get_field(spec,'mode',record_average.info.mode)));
spec.roles=roles;
output_files=struct('trace',struct(),'grating',struct(), ...
    'grating_trials',struct(),'files',strings(0,1));
for role=roles'
    for stage=["sensitivity","snr"]
        data=record_average.channels.(char(role)).(char(stage));
        stem=trace_stem(mode,role,stage);
        files=plot_trace_summary(data,role,stage,output_dir,stem,record_average.stim_windows);
        output_files.trace.(char(role)).(char(stage))=files;
        output_files.files=[output_files.files;struct_files(files)]; %#ok<AGROW>
    end
end
if isstruct(grating)&&isfield(grating,'status')&&grating.status=="completed"
    output_files.grating=render_grating(grating,roles,mode,output_dir);
    output_files.files=[output_files.files;string(output_files.grating.files(:))];
    output_files.grating_trials=render_record_grating_trials( ...
        entries,grating,output_dir,spec,progress);
    output_files.files=[output_files.files;string(output_files.grating_trials.files(:))];
end
output_files.files=unique(output_files.files(strlength(output_files.files)>0),'stable');
end

function files=plot_trace_summary(data,role,stage,folder,stem,windows)
fig=figure('Color','w','Visible','off','Name',char(role+" "+stage+" Record summary"), ...
    'Position',[80 80 1450 820]);layout=tiledlayout(fig,2,2,'Padding','compact','TileSpacing','compact');
ax=nexttile(layout);imagesc(ax,data.time,1:data.ncycles,data.cycle_roi_mean');axis(ax,'xy');colorbar(ax);xlabel(ax,'Relative time (s)');ylabel(ax,'Cycle');title(ax,role+" "+stage+" | Cycle x time");shade_boundaries(ax,windows,role);
ax=nexttile(layout);hold(ax,'on');colors=parula(max(2,data.ncycles));for idx=1:data.ncycles,plot(ax,data.time,data.cycle_roi_mean(:,idx),'Color',colors(idx,:));end;xlabel(ax,'Relative time (s)');ylabel(ax,'ROI mean');title(ax,'Per-cycle ROI mean');shade_ranges(ax,windows,role);
ax=nexttile(layout);population=squeeze(mean(data.average,2,'omitnan'));plot(ax,data.time,population,'k','LineWidth',1.5);xlabel(ax,'Relative time (s)');ylabel(ax,'Mean');title(ax,'Record average');shade_ranges(ax,windows,role);
ax=nexttile(layout);imagesc(ax,data.time,1:data.nrois,data.average');axis(ax,'xy');colorbar(ax);xlabel(ax,'Relative time (s)');ylabel(ax,'ROI');title(ax,'Cycle-average ROI traces');shade_boundaries(ax,windows,role);
sgtitle(fig,sprintf('%s %s Record summary | %d cycles',upper(char(role)),char(stage),data.ncycles));
files=save_bundle(fig,folder,stem);close(fig);
end

function output_files=render_grating(summary,roles,mode,folder)
output_files=struct('files',strings(0,1),'variants',struct());
variants=repmat(struct('field',"",'label',"",'units',"",'role',""),0,1);
if ismember("voltage",roles)
    if completed(summary.voltage_peak_count),variants(end+1)=struct('field',"voltage_peak_count",'label',"Voltage peak count",'units',"Peaks / trial",'role',"voltage");end %#ok<AGROW>
    if completed(summary.voltage_sensitivity),variants(end+1)=struct('field',"voltage_sensitivity",'label',"Voltage mean sensitivity",'units',"Mean sensitivity",'role',"voltage");end %#ok<AGROW>
    if isempty(variants)&&completed(summary.voltage),variants(end+1)=struct('field',"voltage",'label',"Voltage response",'units',"Response",'role',"voltage");end %#ok<AGROW>
end
if ismember("calcium",roles)&&completed(summary.calcium)
    variants(end+1)=struct('field',"calcium",'label',"Calcium mean sensitivity",'units',"Mean sensitivity",'role',"calcium"); %#ok<AGROW>
end
for idx=1:numel(variants)
    variant=variants(idx);value=summary.(char(variant.field));
    stem=grating_stem(mode,variant.field);
    files=plot_grating_variant(summary.orientations,value,variant,folder,stem);
    output_files.variants.(char(variant.field))=files;
    output_files.files=[output_files.files;struct_files(files)]; %#ok<AGROW>
    roi_files=plot_grating_rois(summary.orientations,value,variant,folder,stem+"_roi");
    output_files.files=[output_files.files;roi_files]; %#ok<AGROW>
    csv_file=write_metric_table(value,folder,stem+"_metrics.csv");output_files.files(end+1,1)=csv_file; %#ok<AGROW>
    cycle_files=plot_cycle_tuning(summary.orientations,value,variant,folder,cycle_stem(mode,variant.field));
    output_files.files=[output_files.files;cycle_files]; %#ok<AGROW>
end
output_files.files=unique(output_files.files,'stable');
end

function files=plot_grating_variant(angles,value,variant,folder,stem)
fig=figure('Color','w','Visible','off','Position',[100 90 1300 760]);layout=tiledlayout(fig,2,2,'Padding','compact','TileSpacing','compact');
ax=nexttile(layout);population_curve(ax,angles,value,variant);
ax=nexttile(layout);metric_distribution(ax,value,'dsi','DSI');
ax=nexttile(layout);metric_distribution(ax,value,'osi','OSI');
ax=nexttile(layout);radar_plot(ax,angles,mean(value.stim_mean,1,'omitnan'),variant.label);
sgtitle(fig,variant.label+" Record tuning | "+numel(value.included_cycles)+" cycles");files=save_bundle(fig,folder,stem);close(fig);
end

function files=plot_grating_rois(angles,value,variant,folder,dir_name)
target=fullfile(folder,dir_name);if ~isfolder(target),mkdir(target);end
nrois=size(value.stim_mean,1);files=strings(0,1);
for roi=1:nrois
    fig=figure('Color','w','Visible','off','Position',[120 100 1000 650]);layout=tiledlayout(fig,1,2,'Padding','compact','TileSpacing','compact');
    ax=nexttile(layout);errorbar(ax,angles,value.stim_mean(roi,:),value.stim_sem(roi,:),'-o','LineWidth',1.4);hold(ax,'on');errorbar(ax,angles,value.nonstim_mean(roi,:),value.nonstim_sem(roi,:),'--o');xlabel(ax,'Direction (deg)');ylabel(ax,variant.units);title(ax,variant.label+" ROI "+roi);legend(ax,{'Stim','Baseline'});
    ax=nexttile(layout);bar(ax,[value.metrics_mean_across_cycles.dsi.mean(roi),value.metrics_from_mean_curve.dsi(roi);value.metrics_mean_across_cycles.osi.mean(roi),value.metrics_from_mean_curve.osi(roi)]);set(ax,'XTickLabel',{'DSI','OSI'});legend(ax,{'Mean cycle metric','Metric from mean curve'});ylim(ax,[0 1]);
    stem=string(sprintf('r%03d_tuning',roi));bundle=save_bundle(fig,target,stem);close(fig);files=[files;struct_files(bundle)]; %#ok<AGROW>
end
end

function files=plot_cycle_tuning(angles,value,variant,folder,dir_name)
target=fullfile(folder,dir_name);if ~isfolder(target),mkdir(target);end
nrois=size(value.stim_per_cycle,1);ncycles=size(value.stim_per_cycle,3);files=strings(0,1);
for roi=1:nrois
    fig=figure('Color','w','Visible','off','Position',[100 80 1200 720]);ax=axes(fig);hold(ax,'on');colors=parula(max(2,ncycles));
    for cycle=1:ncycles,plot(ax,angles,value.stim_per_cycle(roi,:,cycle),'-o','Color',colors(cycle,:));end
    plot(ax,angles,value.stim_mean(roi,:),'k-o','LineWidth',2.2,'MarkerFaceColor','k');xlabel(ax,'Direction (deg)');ylabel(ax,variant.units);title(ax,variant.label+" ROI "+roi+" per-cycle tuning");legend(ax,[cellstr(value.included_cycles);{'Mean'}],'Location','eastoutside');
    stem=string(sprintf('r%03d_cycle_tuning',roi));bundle=save_bundle(fig,target,stem);close(fig);files=[files;struct_files(bundle)]; %#ok<AGROW>
end
end

function file=write_metric_table(value,folder,name)
nrois=size(value.stim_mean,1);ROI=(1:nrois)';
DSI_MeanAcrossCycles=value.metrics_mean_across_cycles.dsi.mean;DSI_SEM=value.metrics_mean_across_cycles.dsi.sem;DSI_FromMeanCurve=value.metrics_from_mean_curve.dsi;
OSI_MeanAcrossCycles=value.metrics_mean_across_cycles.osi.mean;OSI_SEM=value.metrics_mean_across_cycles.osi.sem;OSI_FromMeanCurve=value.metrics_from_mean_curve.osi;
T=table(ROI,DSI_MeanAcrossCycles,DSI_SEM,DSI_FromMeanCurve,OSI_MeanAcrossCycles,OSI_SEM,OSI_FromMeanCurve);file=string(fullfile(folder,name));writetable(T,file);
end

function population_curve(ax,angles,value,variant)
m=mean(value.stim_mean,1,'omitnan');sem=std(value.stim_mean,0,1,'omitnan')./sqrt(max(1,sum(isfinite(value.stim_mean),1)));errorbar(ax,angles,m,sem,'-o','LineWidth',1.5);hold(ax,'on');b=mean(value.nonstim_mean,1,'omitnan');bsem=std(value.nonstim_mean,0,1,'omitnan')./sqrt(max(1,sum(isfinite(value.nonstim_mean),1)));errorbar(ax,angles,b,bsem,'--o');xlabel(ax,'Direction (deg)');ylabel(ax,variant.units);title(ax,variant.label);legend(ax,{'Stim','Baseline'});
end
function metric_distribution(ax,value,name,label)
x=value.metrics_mean_across_cycles.(name).mean;y=value.metrics_from_mean_curve.(name);scatter(ax,x,y,38,'filled');hold(ax,'on');plot(ax,[0 1],[0 1],'k--');axis(ax,[0 1 0 1]);axis(ax,'square');xlabel(ax,'Mean cycle metric');ylabel(ax,'Metric from mean curve');title(ax,label);
end
function radar_plot(ax,angles,response,label)
values=max(double(response(:)'),0);mx=max(values,[],'omitnan');if ~isfinite(mx)||mx<=0,values(:)=0;else,values=values/mx;end;theta=deg2rad(mod(angles,360));[theta,order]=sort(theta);values=values(order);theta=[theta theta(1)];values=[values values(1)];cla(ax);hold(ax,'on');axis(ax,'equal');axis(ax,'off');for radius=.25:.25:1,plot(ax,radius*cos(linspace(0,2*pi,181)),radius*sin(linspace(0,2*pi,181)),'Color',[.85 .85 .85]);end;plot(ax,values.*cos(theta),values.*sin(theta),'-o','LineWidth',1.5);xlim(ax,[-1.2 1.2]);ylim(ax,[-1.2 1.2]);title(ax,label+" normalized tuning");
end
function tf=completed(value),tf=isstruct(value)&&isfield(value,'status')&&string(value.status)=="completed";end
function stem=trace_stem(mode,role,stage)
if mode=="dual"&&role=="voltage"&&stage=="sensitivity",stem="4_roi_cyc_hm_v_sens";elseif mode=="dual"&&role=="voltage"&&stage=="snr",stem="4_roi_cyc_hm_v_snr";else,stem="4_roi_cyc_"+extractBefore(role,2)+"_"+stage;end
end
function stem=grating_stem(mode,field)
if field=="voltage_peak_count",stem="8_grating";elseif field=="voltage_sensitivity",stem="8_grating_v_sens";elseif field=="calcium",stem="8_grating_calcium";else,stem="8_grating_"+field;end
if mode=="ap"&&field=="voltage_peak_count",stem="8_grating_peak_count";end
end
function stem=cycle_stem(mode,field)
if mode=="dual",stem="11_rg_v_ca_"+field;else,stem="11_rg_voltage_"+erase(field,"voltage_");end
end
function files=save_bundle(fig,folder,stem)
files=struct('fig_file',string(fullfile(folder,stem+".fig")),'png_file',string(fullfile(folder,stem+".png")));savefig(fig,files.fig_file);exportgraphics(fig,files.png_file,'Resolution',150);
end
function files=struct_files(value)
names=fieldnames(value);files=strings(0,1);for idx=1:numel(names),candidate=string(value.(names{idx}));if isscalar(candidate)&&strlength(candidate)>0,files(end+1,1)=candidate;end,end
end
function shade_ranges(ax,windows,role)
ranges=role_ranges(windows,role);if isempty(ranges),return;end;yl=ylim(ax);for idx=1:size(ranges,1),patch(ax,[ranges(idx,1) ranges(idx,2) ranges(idx,2) ranges(idx,1)],[yl(1) yl(1) yl(2) yl(2)],[.75 .75 .75],'FaceAlpha',.12,'EdgeColor','none');end;uistack(findobj(ax,'Type','line'),'top');
end
function shade_boundaries(ax,windows,role)
ranges=role_ranges(windows,role);for idx=1:size(ranges,1),xline(ax,ranges(idx,1),'w--');xline(ax,ranges(idx,2),'w--');end
end
function ranges=role_ranges(windows,role)
ranges=[];if ~isstruct(windows)||~isfield(windows,char(role)),return;end;value=windows.(char(role));if isfield(value,'stim_time_ranges'),ranges=value.stim_time_ranges;elseif isfield(value,'flash_time_ranges'),ranges=value.flash_time_ranges;end
end
function value=get_field(s,name,fallback),if isfield(s,name),value=s.(name);else,value=fallback;end,end
