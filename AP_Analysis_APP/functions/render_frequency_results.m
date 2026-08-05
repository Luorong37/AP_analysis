function output_files=render_frequency_results(results,output_path,stim_windows)
%RENDER_FREQUENCY_RESULTS Save frozen FFT/wavelet summary and ROI pairs.
if nargin<3,stim_windows=struct();end;roles=["voltage","calcium"];output_files=struct();
fig=figure('Color','w','Visible','off');layout=tiledlayout(fig,2,2,'Padding','compact','TileSpacing','compact');
for idx=1:2,r=results.(roles(idx));ax=nexttile(layout);valid=r.fft.frequency<=r.parameters.max_freq_hz;plot(ax,r.fft.frequency(valid),mean(r.fft.amplitude(valid,:),2,'omitnan'),'LineWidth',1.2);xlabel(ax,'Hz');ylabel(ax,'Amplitude');title(ax,roles(idx)+" FFT");grid(ax,'on');ax=nexttile(layout);histogram(ax,r.roi.peak_frequency_hz);xlabel(ax,'Peak Hz');title(ax,roles(idx)+" ROI peak frequency");end
output_files.fourier_summary=save_bundle(fig,output_path,'8_fourier_summary');close(fig);
fig=figure('Color','w','Visible','off');layout=tiledlayout(fig,2,1,'Padding','compact','TileSpacing','compact');for idx=1:2,ax=nexttile(layout);plot_wavelet(ax,results.(roles(idx)).representative.wavelet,roles(idx)+" representative wavelet");end;output_files.wavelet_summary=save_bundle(fig,output_path,'8_wavelet_summary');close(fig);
folder=fullfile(output_path,'8_wavelet_roi_pairs');if ~isfolder(folder),mkdir(folder);end;n=min(size(results.voltage.trace_matrix,2),size(results.calcium.trace_matrix,2));
roi_files=strings(0,1);for roi=1:n,fig=figure('Color','w','Visible','off','Position',[100 100 1200 850]);layout=tiledlayout(fig,2,2,'Padding','compact','TileSpacing','compact');for idx=1:2,role=roles(idx);r=results.(role);x=r.trace_matrix(:,roi);ax=nexttile(layout);plot(ax,r.time,x);shade(ax,stim_windows,role);title(ax,role+" ROI "+roi);xlabel(ax,'Time (s)');ax=nexttile(layout);w=local_cwt(x,r.frame_rate,r.parameters);plot_wavelet(ax,w,role+" wavelet");end;fig_file=string(fullfile(folder,sprintf('ROI_%03d.fig',roi)));png_file=string(fullfile(folder,sprintf('ROI_%03d.png',roi)));savefig(fig,fig_file);exportgraphics(fig,png_file,'Resolution',150);close(fig);roi_files=[roi_files;fig_file;png_file];end %#ok<AGROW>
output_files.wavelet_summary.roi_pair_files=roi_files;
end
function w=local_cwt(x,fs,params)
w=struct('available',false,'time',(0:numel(x)-1)'/fs,'frequency',[],'coefficients',[]);if exist('cwt','file')~=2,return;end
x=fillmissing(double(x(:)),'linear','EndValues','nearest');x=x-mean(x,'omitnan');try,[wt,f]=cwt(x,char(string(params.wavelet_name)),fs,'VoicesPerOctave',double(params.wavelet_voices_per_octave));valid=f>=params.min_freq_hz&f<=params.max_freq_hz;w.available=true;w.frequency=f(valid);w.coefficients=wt(valid,:);catch,end
end
function plot_wavelet(ax,w,title_text)
if ~isfield(w,'available')||~w.available||isempty(w.coefficients),text(ax,.5,.5,'CWT unavailable','HorizontalAlignment','center');axis(ax,'off');return;end;imagesc(ax,w.time,w.frequency,abs(w.coefficients).^2);set(ax,'YDir','normal');xlabel(ax,'Time (s)');ylabel(ax,'Hz');title(ax,title_text);colorbar(ax);
end
function shade(ax,windows,role)
if ~isstruct(windows)||~isfield(windows,role)||~isfield(windows.(role),'stim_time_ranges'),return;end;hold(ax,'on');yl=ylim(ax);r=windows.(role).stim_time_ranges;for i=1:size(r,1),patch(ax,[r(i,1) r(i,2) r(i,2) r(i,1)],[yl(1) yl(1) yl(2) yl(2)],[.5 .7 1],'FaceAlpha',.12,'EdgeColor','none');end
end
function f=save_bundle(fig,folder,stem),stem=string(stem);f=struct('fig_file',string(fullfile(folder,stem+".fig")),'png_file',string(fullfile(folder,stem+".png")));savefig(fig,f.fig_file);exportgraphics(fig,f.png_file,'Resolution',150);end
