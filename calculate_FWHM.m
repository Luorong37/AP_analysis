function FWHM = calculate_FWHM(AP_amp, dt, AP_window_width)
    % Calculate the baseline and half-maximum value
    AP_base = mean(AP_amp,'omitnan');
    HM = 0.5 * (AP_base + AP_amp(AP_window_width+1));

    % Calculate FWHM
    FWHM = NaN; % Initialize FWHM to NaN in case exact indices are not found

    % Find indices for half-max on both sides of the peak
    left_idx = find(AP_amp(1:AP_window_width+1) < HM, 1, 'last');
    right_idx = find(AP_amp(AP_window_width+1:end) < HM, 1, 'first') + AP_window_width - 1;

    % Linear interpolation to find exact points at half-maximum
    if ~isempty(left_idx) && left_idx < AP_window_width+1
        x1 = left_idx;
        x2 = left_idx + 1;
        y1 = AP_amp(x1);
        y2 = AP_amp(x2);
        left_HM_exact = ((HM - y1) * (x2 - x1) / (y2 - y1)) + x1;
    else
        left_HM_exact = NaN;
    end

    if ~isempty(right_idx) && right_idx > AP_window_width+1
        x1 = right_idx - 1;
        x2 = right_idx;
        y1 = AP_amp(x1);
        y2 = AP_amp(x2);
        right_HM_exact = ((HM - y1) * (x2 - x1) / (y2 - y1)) + x1;
    else
        right_HM_exact = NaN;
    end

    % Calculate FWHM if both left and right half-maximum points are found
    if ~isnan(left_HM_exact) && ~isnan(right_HM_exact)
        FWHM = (right_HM_exact - left_HM_exact) * dt;
    end
end
