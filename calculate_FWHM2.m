function FWHM = calculate_FWHM2(AP_amp, dt, polarity)
[~,peak_x,FWHM_index] = findpeaks(AP_amp*polarity,'Annotate','extents','MinPeakDistance',length(AP_amp)*0.9,'WidthReference','halfheight');
FWHM = FWHM_index*dt;
AP_window_width = floor((numel(AP_amp)-1)/2);

% false positive AP
if all((AP_window_width-1) <= peak_x) && all(peak_x <= (AP_window_width+3))

else
    FWHM = NaN;
end

if isempty(FWHM) || all(FWHM < 0.002)
        FWHM = NaN;
end


if isnan(FWHM)
     findpeaks(AP_amp*polarity,1/dt,'Annotate','extents','MinPeakDistance',length(AP_amp)*dt*0.9,'WidthReference','halfheight')

end


end
