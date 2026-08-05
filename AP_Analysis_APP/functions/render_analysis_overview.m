function output_files=render_analysis_overview(channels,output_path,spec,shared_results)
%RENDER_ANALYSIS_OVERVIEW Restore frozen Dual overview output bundle.
if nargin<4||isempty(shared_results),shared_results=struct();end
output_path=string(output_path);if ~isfolder(output_path),mkdir(output_path);end
mode=lower(string(get_field(spec,'mode',"dual")));output_files=struct();
if mode~="dual"
    fig=plot_trace_overview([channels.results],[channels.profile]);output_files.trace_overview=save_bundle(fig,output_path,'4_trace_overview');close(fig);return;
end
vi=find_role(channels,'voltage');ci=find_role(channels,'calcium');v=channels(vi);c=channels(ci);
[maps,rois]=resolve_map_roi(shared_results,v,c);output_files=render_maps(maps,rois,v,c,output_path,output_files);
output_files.raw_trace=raw_trace(v,c,output_path);
output_files=render_trace_qc(v,c,output_path,output_files,shared_results);
output_files.trace_summary=trace_summary(v,c,output_path,spec);
output_files.roi_heatmap=roi_heatmap(v,c,output_path,spec);
end

function [maps,rois]=resolve_map_roi(shared_results,v,c)
maps=struct();rois=struct();
if isfield(shared_results,'maps'),maps=shared_results.maps;end
for role=["voltage","calcium"]
    if isfield(maps,char(role))&&isstruct(maps.(char(role))) ...
            &&isfield(maps.(char(role)),'data')
        maps.(char(role))=maps.(char(role)).data;
    end
end
if ~isfield(maps,'voltage')&&isfield(v,'data')&&isfield(v.data,'movie_3d'),[maps.voltage,~]=create_map(v.data.movie_3d,v.profile);end
if ~isfield(maps,'calcium')&&isfield(c,'data')&&isfield(c.data,'movie_3d'),[maps.calcium,~]=create_map(c.data.movie_3d,c.profile);end
if isfield(shared_results,'roi')&&isfield(shared_results.roi,'rois'),rois=shared_results.roi.rois;end
end
function output_files=render_maps(maps,rois,v,c,folder,output_files)
if isfield(maps,'voltage')&&isfield(maps,'calcium')
    fig=figure('Color','w','Visible','off');subplot(1,2,1);imagesc(maps.voltage);axis image off;colorbar;title('Voltage sensitivity map');subplot(1,2,2);imagesc(maps.calcium);axis image off;colorbar;title('Calcium sensitivity map');output_files.sensitivity_map=save_bundle(fig,folder,'0_dual_sensitivity_map');close(fig);
    fig=figure('Color','w','Visible','off');subplot(1,2,1);imagesc(maps.voltage);axis image off;hold on;if isfield(rois,'bwmask'),contour(rois.bwmask,[.5 .5],'r');end;title('Voltage map + ROI');subplot(1,2,2);imagesc(maps.calcium);axis image off;hold on;if isfield(rois,'bwmask_ca'),contour(rois.bwmask_ca,[.5 .5],'r');end;title('Calcium map + ROI');output_files.map_with_roi=save_bundle(fig,folder,'0_dual_sensitivity_map_with_roi');close(fig);
end
if isfield(v.data,'movie_3d')&&isfield(c.data,'movie_3d')
    av=mean(v.data.movie_3d,3,'double');ac=mean(c.data.movie_3d,3,'double');av=scale01(av);ac=scale01(ac);rgb=cat(3,av,ac,zeros(size(av)));output_files.average_merge=write_image_pair(rgb,folder,'0_dual_average_color_merge');
    if isfield(rois,'bwmask')
        edge=bwperim(rois.bwmask>0);rgb_roi=rgb;rgb_roi(:,:,1)=max(rgb_roi(:,:,1),edge);rgb_roi(:,:,2)=rgb_roi(:,:,2).*(~edge);rgb_roi(:,:,3)=max(rgb_roi(:,:,3),edge);output_files.average_merge_roi=write_image_pair(rgb_roi,folder,'0_dual_average_color_merge_with_roi');
    end
end
end
function files=raw_trace(v,c,folder)
vr=fetch_trace_stage(v.results,'raw');cr=fetch_trace_stage(c.results,'raw');fig=figure('Color','w','Visible','off');subplot(2,1,1);plot_offset(gca,vr,v.profile.frame_rate);title('Voltage raw');subplot(2,1,2);plot_offset(gca,cr,c.profile.frame_rate);title('Calcium raw');files=save_bundle(fig,folder,'1_dual_raw_trace');close(fig);
end
function output_files=render_trace_qc(v,c,folder,output_files,shared_results)
if isfield(shared_results,'trace_details')&&isfield(shared_results.trace_details,'background')
    b=shared_results.trace_details.background;
    if isfield(b,'voltage')&&isfield(b,'calcium')
        fig=figure('Color','w','Visible','off');subplot(2,2,1);plot(b.voltage.background_raw);title('Voltage background');subplot(2,2,2);plot(b.voltage.background_fitted);title('Voltage fitted background');subplot(2,2,3);plot(b.calcium.background_raw);title('Calcium background');subplot(2,2,4);plot(b.calcium.background_fitted);title('Calcium fitted background');output_files.background=save_bundle(fig,folder,'1_dual_background_correction_summary');close(fig);
    end
end
try
    vb=fetch_trace_stage(v.results,'bleach_removed');cb=fetch_trace_stage(c.results,'bleach_removed');fig=figure('Color','w','Visible','off');subplot(2,1,1);plot_offset(gca,vb,v.profile.frame_rate);title('Voltage bleach removed');subplot(2,1,2);plot_offset(gca,cb,c.profile.frame_rate);title('Calcium bleach removed');output_files.bleach=save_bundle(fig,folder,'2_dual_bleach_correction_stacked');close(fig);
catch
end
try
    vb=fetch_trace_stage(v.results,'bleach_removed');vn=fetch_trace_stage(v.results,'noise_reference');cb=fetch_trace_stage(c.results,'bleach_removed');cn=fetch_trace_stage(c.results,'noise_reference');fig=figure('Color','w','Visible','off');subplot(2,1,1);plot_offset(gca,[vb vn],v.profile.frame_rate);title('Voltage bleach / noise reference');subplot(2,1,2);plot_offset(gca,[cb cn],c.profile.frame_rate);title('Calcium bleach / noise reference');output_files.noise_reference=save_bundle(fig,folder,'3_dual_noise_reference_comparison');close(fig);
catch
end
end
function files=trace_summary(v,c,folder,spec)
fig=figure('Color','w','Visible','off');names={{'raw'},{'sensitivity'},{'snr'},{'raw_smoothed','raw'},{'sensitivity_smoothed','sensitivity'},{'snr_smoothed','snr'}};channels={v,v,v,c,c,c};titles=["Voltage Raw","Voltage Sensitivity","Voltage SNR","Calcium Raw","Calcium Sensitivity","Calcium SNR"];
for i=1:6,ax=subplot(2,3,i);try,[x,~]=resolve_trace_stage(channels{i}.results,names{i});if i==2||i==3,x=double(get_field(spec,'voltage_polarity',1))*x;elseif i==5||i==6,x=double(get_field(spec,'calcium_polarity',1))*x;end;plot_offset(ax,x,channels{i}.profile.frame_rate);catch,text(ax,.5,.5,'Unavailable','HorizontalAlignment','center');end;title(ax,titles(i));end
files=save_bundle(fig,folder,'4_dual_trace_summary');close(fig);
end
function files=roi_heatmap(v,c,folder,spec)
files=struct();try,vs=fetch_trace_stage(v.results,'sensitivity');[cs,cstage]=resolve_trace_stage(c.results,{'sensitivity_smoothed','sensitivity'});catch,return;end
n=min(size(vs,2),size(cs,2));if n==0,return;end;vs=double(get_field(spec,'voltage_polarity',1))*vs(:,1:n);cs=double(get_field(spec,'calcium_polarity',1))*cs(:,1:n);base=prctile(cs,1,1);cz=max(0,cs-base);vr=max(vs,[],1,'omitnan')-min(vs,[],1,'omitnan');scale=max(vr);if ~isfinite(scale)||scale<=0,scale=1;end;mid=(max(vs,[],1,'omitnan')+min(vs,[],1,'omitnan'))/2;
fig=figure('Color','w','Visible','off');ax=axes(fig);tc=(1:size(cs,1))'/c.profile.frame_rate;tv=(1:size(vs,1))'/v.profile.frame_rate;imagesc(ax,tc,1:n,cz');set(ax,'YDir','normal');colormap(ax,ice_map(256,.6,.4));hold(ax,'on');for roi=1:n,plot(ax,tv,roi+(vs(:,roi)-mid(roi))/scale*1.5,'r','LineWidth',.45);end;xlabel(ax,'Time (s)');ylabel(ax,'ROI');title(ax,"Calcium "+cstage+" + voltage sensitivity",'Interpreter','none');colorbar(ax);files=save_bundle(fig,folder,'4_dual_roi_calcium_heatmap_voltage_trace');close(fig);
fig=figure('Color','w','Visible','off');x=linspace(0,1,256);plot(x,x.^.6,'--',x,interp1([0 .4 1],[0 .5 1],x.^.6),'b','LineWidth',1.2);grid on;xlabel('Normalized value');ylabel('Colormap index');title('Calcium heatmap colormap transfer');files.colormap_curve=save_bundle(fig,folder,'4_dual_roi_calcium_heatmap_colormap_curve');close(fig);
end
function cmap=ice_map(n,gamma,pivot)
base=zeros(n,3);base(:,3)=linspace(0,1,n);k=max(1,floor(n*.1));base(k:end,2)=linspace(0,1,n-k+1);k=max(1,floor(n*.9));base(k:end,1)=linspace(0,1,n-k+1);x=linspace(0,1,n);w=interp1([0 pivot 1],[0 .5 1],x.^gamma,'linear','extrap');cmap=interp1(x,base,min(max(w,0),1));
end
function plot_offset(ax,x,fps)
x=double(x);t=(1:size(x,1))'/double(fps);span=max(x,[],1,'omitnan')-min(x,[],1,'omitnan');gap=max(eps,median(span,'omitnan'));plot(ax,t,x+(0:size(x,2)-1)*gap);xlabel(ax,'Time (s)');
end
function y=scale01(x),lo=prctile(x(:),1);hi=prctile(x(:),99);y=min(1,max(0,(double(x)-lo)/max(eps,hi-lo)));end
function files=write_image_pair(rgb,folder,stem),stem=string(stem);files=struct('tif_file',string(fullfile(folder,stem+".tif")),'png_file',string(fullfile(folder,stem+".png")));imwrite(rgb,files.tif_file,'tif','Compression','none');imwrite(rgb,files.png_file);end
function files=save_bundle(fig,folder,stem),stem=string(stem);files=struct('fig_file',string(fullfile(folder,stem+".fig")),'png_file',string(fullfile(folder,stem+".png")));savefig(fig,files.fig_file);exportgraphics(fig,files.png_file,'Resolution',150);end
function idx=find_role(channels,role),idx=find(arrayfun(@(x)strcmpi(string(x.profile.role),role),channels),1);end
function value=get_field(s,name,fallback),if isstruct(s)&&isfield(s,name)&&~isempty(s.(name)),value=s.(name);else,value=fallback;end,end
