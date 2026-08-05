function output_files=render_dual_comparison(comparison,output_path,stim_windows)
%RENDER_DUAL_COMPARISON Save frozen Dual overlap population/per-ROI plots.
if nargin<3,stim_windows=struct();end;output_files=struct();
for metric=["sensitivity","snr"]
    if ~isfield(comparison,metric),continue;end;data=comparison.(metric);n=size(data.voltage_display,2);
    fig=figure('Color','w','Visible','off','Position',[100 100 1200 max(420,260*ceil(n/3))]);
    for roi=1:n,ax=subplot(ceil(n/3),min(3,n),roi);overlap(ax,data,roi,stim_windows);title(ax,"ROI "+roi);end;sgtitle("Dual "+upper(metric)+" overlap");
    stem="5_dual_"+metric+"_overlap";output_files.(metric)=save_bundle(fig,output_path,stem);close(fig);
    folder=fullfile(output_path,stem+"_rois");if ~isfolder(folder),mkdir(folder);end
    roi_files=strings(0,1);
    for roi=1:n,fig=figure('Color','w','Visible','off','Position',[100 100 1400 360]);ax=axes(fig);overlap(ax,data,roi,stim_windows);title(ax,"ROI "+roi+" "+metric);fig_file=string(fullfile(folder,sprintf('ROI_%03d.fig',roi)));png_file=string(fullfile(folder,sprintf('ROI_%03d.png',roi)));savefig(fig,fig_file);exportgraphics(fig,png_file,'Resolution',150);close(fig);roi_files=[roi_files;fig_file;png_file];end %#ok<AGROW>
    output_files.(metric).roi_files=roi_files;
end
end
function overlap(ax,data,roi,windows)
hold(ax,'on');shade(ax,windows);yyaxis(ax,'left');plot(ax,data.voltage_time,norm01(data.voltage_display(:,roi)),'r');ylim(ax,[0 1]);ylabel(ax,'Voltage (norm)');yyaxis(ax,'right');plot(ax,data.calcium_time,norm01(data.calcium_display(:,roi)),'g');ylim(ax,[0 1]);ylabel(ax,'Calcium (norm)');xlabel(ax,'Time (s)');grid(ax,'on');
end
function shade(ax,w)
if ~isstruct(w)||~isfield(w,'voltage')||~isfield(w.voltage,'stim_time_ranges'),return;end;r=w.voltage.stim_time_ranges;yl=[0 1];for i=1:size(r,1),patch(ax,[r(i,1) r(i,2) r(i,2) r(i,1)],[yl(1) yl(1) yl(2) yl(2)],[.5 .7 1],'FaceAlpha',.12,'EdgeColor','none');end
end
function y=norm01(x),x=double(x);lo=min(x,[],'omitnan');hi=max(x,[],'omitnan');if ~isfinite(lo)||~isfinite(hi)||hi<=lo,y=zeros(size(x));else,y=(x-lo)/(hi-lo);end,end
function f=save_bundle(fig,folder,stem),stem=string(stem);f=struct('fig_file',string(fullfile(folder,stem+".fig")),'png_file',string(fullfile(folder,stem+".png")));savefig(fig,f.fig_file);exportgraphics(fig,f.png_file,'Resolution',150);end
