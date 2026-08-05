function [result,output_files] = analyze_record_frequency(record_average,spec,output_dir)
%ANALYZE_RECORD_FREQUENCY Compute and render Record-average FFT/CWT by role.

roles=string(record_average.info.roles(:));output_dir=string(output_dir);
result=struct('status',"completed",'channels',struct(),'parameters',spec);
output_files=struct('files',strings(0,1));
for role=roles'
    stage=record_average.channels.(char(role)).sensitivity;
    profile=struct('role',role,'frame_rate',stage.frame_rate);
    channel=analyze_frequency(stage.average,profile,spec);
    result.channels.(char(role))=channel;
    files=render_channel(channel,output_dir,"6_fft_"+role,"7_wav_"+role);
    output_files.(char(role))=files;
    output_files.files=[output_files.files;string(files.fft_fig);string(files.fft_png); ...
        string(files.wavelet_fig);string(files.wavelet_png)]; %#ok<AGROW>
end
output_files.files=unique(output_files.files,'stable');
end

function files=render_channel(value,folder,fft_stem,wavelet_stem)
fig=figure('Color','w','Visible','off','Position',[100 90 1200 720]);layout=tiledlayout(fig,1,2,'Padding','compact','TileSpacing','compact');
ax=nexttile(layout);plot(ax,value.fft.frequency,mean(value.fft.amplitude,2,'omitnan'),'LineWidth',1.5);xlim(ax,[value.parameters.min_freq_hz,min(value.parameters.max_freq_hz,value.frame_rate/2)]);xlabel(ax,'Frequency (Hz)');ylabel(ax,'Mean amplitude');title(ax,value.role+" Record FFT");
ax=nexttile(layout);histogram(ax,value.roi.peak_frequency_hz);xlabel(ax,'Peak frequency (Hz)');ylabel(ax,'ROI count');title(ax,'ROI peak frequency');
files.fft_fig=string(fullfile(folder,fft_stem+".fig"));files.fft_png=string(fullfile(folder,fft_stem+".png"));savefig(fig,files.fft_fig);exportgraphics(fig,files.fft_png,'Resolution',150);close(fig);

fig=figure('Color','w','Visible','off','Position',[100 90 1200 720]);wavelet=value.representative.wavelet;
if isfield(wavelet,'available')&&wavelet.available
    imagesc(wavelet.time,wavelet.frequency,abs(wavelet.coefficients));axis xy;colorbar;xlabel('Time (s)');ylabel('Frequency (Hz)');title(value.role+" representative ROI wavelet");
else
    axis off;text(.5,.5,'CWT unavailable','HorizontalAlignment','center');
end
files.wavelet_fig=string(fullfile(folder,wavelet_stem+".fig"));files.wavelet_png=string(fullfile(folder,wavelet_stem+".png"));savefig(fig,files.wavelet_fig);exportgraphics(fig,files.wavelet_png,'Resolution',150);close(fig);
end
