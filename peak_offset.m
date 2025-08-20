function offset = peak_offset(AP_amp_cut, polarity,f)
if nargin <3
    f = false;
end
[~,peak_x,~] = findpeaks(AP_amp_cut*polarity,'Annotate','extents','MinPeakDistance',length(AP_amp_cut)*0.9,'WidthReference','halfheight');
AP_window_width = floor((numel(AP_amp_cut)-1)/2);
offset = peak_x - AP_window_width -1;
if f
     findpeaks(AP_amp_cut*polarity,'MinPeakDistance',length(AP_amp_cut)*0.9)
end
end
