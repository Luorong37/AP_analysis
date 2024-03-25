function offset_plot(traces_corrected,t,rois,colors, ...
    peaks_index,peak_polarity,peak_threshold,peaks_sensitivity)
% Plot summary
fig = figure();
set(fig,'Position',get(0,'Screensize'));
offset_peak = 0;

if nargin > 4
for i = 1:length(rois)-1
    % plot trace
    trace = (traces_corrected(:,i)-1) * peak_polarity(i) + offset_peak;
    plot(t,  trace ,'Color',colors(i,:)); hold on;
    % plot threshold
    plot(t,ones(1,length(t)).*(1-peak_threshold(i)*peak_polarity(i)) + offset_peak,'Color',colors(i,:),'LineWidth',2); hold on;
    % plot peak
    dt = t(2)-t(1);
    plot(peaks_index{i}.*dt, (1-peaks_sensitivity{i}*peak_polarity(i))+ offset_peak,'v','Color',colors(i,:),'MarkerFaceColor',colors(i,:)); hold on;
    offset_peak_step = (max(trace(:))-min(trace(:)))*1; 
    offset_peak = offset_peak + offset_peak_step;
end

else
    for i = 1:length(rois)-1
    % plot trace
    trace = (traces_corrected(:,i)-1) + offset_peak;
    plot(t,  trace ,'Color',colors(i,:)); hold on;
    end
end
ylim tight
sgtitle('Peak finding');
hold on;
