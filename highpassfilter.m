function traces_filt = highpassfilter(traces, Hz, Fc_default)
% This function performs photobleaching correction on the intensity-time
% traces using a high-pass filter.
% Input arguments:
% traces: a matrix of size num_frames x num_rois containing the intensity-time traces for each ROI
% Hz: the frame rate in Hz
% Fc_default: the cutoff frequency for the high-pass filter in Hz. Default value is 1/3.

% Set default value for Fc
if nargin < 3
    Fc = 1/3;
else
    Fc = Fc_default;
end

% Calculate time axis
dt = 1 / Hz;

% Design high-pass filter
[b, a] = butter(3, Fc*dt, 'high'); % 3rd order Butterworth high-pass filter

% Apply high-pass filter to remove photobleaching from each ROI
traces_filt = zeros(size(traces));
for i = 1:size(traces, 2)
    traces_filt(:, i) = filtfilt(b, a, traces(:, i));
end

end

