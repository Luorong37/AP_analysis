function output_files = render_visualstim_results(stim_results, channels, output_path, spec)
%RENDER_VISUALSTIM_RESULTS Save AP/Dual stimulus figures and per-ROI tuning.
if nargin<4,spec=struct();end
output_path=string(output_path);if ~isfolder(output_path),mkdir(output_path);end
mode=lower(string(get_field(spec,'mode',"dual")));alpha=double(get_field(spec,'stim_shading_alpha',0.14));
output_files=struct();
if ~isfield(stim_results,'supported')||~stim_results.supported,return;end

n=numel(channels);fig=figure('Color','w','Visible','off');layout=tiledlayout(fig,n,2,'Padding','compact','TileSpacing','compact');
for idx=1:n
    role=lower(string(channels(idx).profile.role));
    for metric=["sensitivity","snr"]
        ax=nexttile(layout);[trace,stage]=safe_stage(channels(idx).results,metric,role);
        if isempty(trace),text(ax,.5,.5,"Unavailable",'HorizontalAlignment','center');axis(ax,'off');continue;end
        t=(1:size(trace,1))'/double(channels(idx).profile.frame_rate);plot(ax,t,offset(trace),'LineWidth',.65);hold(ax,'on');
        shade(ax,stim_results.windows.(char(role)).stim_time_ranges,alpha);title(ax,role+" "+stage,'Interpreter','none');xlabel(ax,'Time (s)');
    end
end
stem="4_"+mode+"_trace_summary_with_stim";output_files.trace=save_bundle(fig,output_path,stem);close(fig);

if isfield(stim_results,'tuning')
    role_names=string(fieldnames(stim_results.tuning));
    fig=figure('Color','w','Visible','off');layout=tiledlayout(fig,1,numel(role_names),'Padding','compact','TileSpacing','compact');
    for idx=1:numel(role_names),ax=nexttile(layout);plot_population(ax,stim_results.tuning.(role_names(idx)),role_names(idx));end
    output_files.tuning_summary=save_bundle(fig,output_path,'7_stim_tuning_summary');close(fig);
    roi_dir=fullfile(output_path,'7_stim_tuning_by_roi');if ~isfolder(roi_dir),mkdir(roi_dir);end
    nrois=max(arrayfun(@(r)size(stim_results.tuning.(r).response_by_condition,1),role_names));
    roi_files=strings(0,1);
    for roi=1:nrois
        fig=figure('Color','w','Visible','off');layout=tiledlayout(fig,1,numel(role_names),'Padding','compact','TileSpacing','compact');
        for idx=1:numel(role_names),role=role_names(idx);t=stim_results.tuning.(role);ax=nexttile(layout);plot_roi(ax,t,roi,role);end
        label=string(sprintf('roi_%03d_stim_tuning',roi));fig_file=string(fullfile(roi_dir,label+".fig"));png_file=string(fullfile(roi_dir,label+".png"));savefig(fig,fig_file);exportgraphics(fig,png_file,'Resolution',150);close(fig);roi_files=[roi_files;fig_file;png_file]; %#ok<AGROW>
    end
    output_files.tuning_by_roi=roi_files;
elseif isfield(stim_results,'metrics')&&isfield(stim_results.metrics,'snr')
    output_files.block_snr=condition_plot(stim_results.metrics.snr,stim_results.windows,output_path,'7_stim_delta_summary_snr','SNR Delta Response');
    if isfield(stim_results.metrics,'sensitivity')
        output_files.block_sensitivity=condition_plot(stim_results.metrics.sensitivity,stim_results.windows,output_path,'7_stim_delta_summary_sensitivity','Sensitivity Delta Response');
    end
end
end

function [trace,stage]=safe_stage(results,metric,role)
trace=[];stage="";if metric=="sensitivity"&&role=="calcium",names={'sensitivity_smoothed','sensitivity'};else,names={char(metric)};end
try,[trace,stage]=resolve_trace_stage(results,names);catch,end
end
function y=offset(x)
x=double(x);span=max(x,[],1,'omitnan')-min(x,[],1,'omitnan');gap=max(eps,median(span,'omitnan'));y=x+(0:size(x,2)-1)*gap;
end
function shade(ax,ranges,alpha)
if isempty(ranges),return;end;yl=ylim(ax);h=gobjects(size(ranges,1),1);for i=1:size(ranges,1),h(i)=patch(ax,[ranges(i,1),ranges(i,2),ranges(i,2),ranges(i,1)],[yl(1),yl(1),yl(2),yl(2)],[.35 .55 1],'FaceAlpha',alpha,'EdgeColor','none');end;uistack(h,'bottom');
end
function files=save_bundle(fig,folder,stem)
stem=string(stem);files=struct('fig_file',string(fullfile(folder,stem+".fig")),'png_file',string(fullfile(folder,stem+".png")));savefig(fig,files.fig_file);exportgraphics(fig,files.png_file,'Resolution',150);
end
function plot_population(ax,t,name)
angles=t.angles;m=mean(t.response_by_condition,1,'omitnan');sem=std(t.response_by_condition,0,1,'omitnan')./sqrt(max(1,sum(isfinite(t.response_by_condition),1)));errorbar(ax,angles,m,sem,'-o','LineWidth',1.2);hold(ax,'on');b=mean(t.nonstim_response_by_condition,1,'omitnan');plot(ax,angles,b,'--o');xlabel(ax,'Direction (deg)');ylabel(ax,'Response');title(ax,name);legend(ax,{'Stim','Baseline'});grid(ax,'on');
end
function plot_roi(ax,t,roi,name)
if roi>size(t.response_by_condition,1),axis(ax,'off');return;end;plot(ax,t.angles,t.response_by_condition(roi,:),'-o','LineWidth',1.2);hold(ax,'on');plot(ax,t.angles,t.nonstim_response_by_condition(roi,:),'--o');title(ax,name+" ROI "+roi);xlabel(ax,'Direction (deg)');grid(ax,'on');
end
function files=condition_plot(metrics,windows,folder,stem,title_text)
roles=string(fieldnames(metrics));fig=figure('Color','w','Visible','off');layout=tiledlayout(fig,1,numel(roles),'Padding','compact');for r=1:numel(roles),ax=nexttile(layout);m=metrics.(roles(r)).delta_mean;means=NaN(numel(windows.condition_labels),1);sem=means;for c=1:numel(means),x=m(:,windows.condition_index==c);means(c)=mean(x,'all','omitnan');sem(c)=std(x,0,'all','omitnan')/sqrt(max(1,sum(isfinite(x),'all')));end;errorbar(ax,1:numel(means),means,sem,'o-');xticks(ax,1:numel(means));xticklabels(ax,windows.condition_labels);title(ax,roles(r));ylabel(ax,'Stim - baseline');grid(ax,'on');end;sgtitle(layout,title_text);files=save_bundle(fig,folder,stem);close(fig);
end
function value=get_field(s,name,fallback)
if isstruct(s)&&isfield(s,name)&&~isempty(s.(name)),value=s.(name);else,value=fallback;end
end
