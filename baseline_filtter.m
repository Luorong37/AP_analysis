!gitfunction baseline = baseline_filtter(t, trace)
    % extractBaseline - Extracts the baseline from a intensity-time series using low-pass filtering.
    %
    % Syntax:  baseline = baseline_filtter(t, trace)
    %
    % Inputs:
    %    t - A vector of time stamps
    %    trace- A vector of voltage measurements corresponding to the time stamps
    %
    % Outputs:
    %    baseline - A vector of the estimated baseline traces

    % Sampling rate calculation
    Fs = 1/mean(diff(t));% Average sampling frequency in Hz

    % Filter design
    cutoffFrequency = 20; % cutoff frequency in Hz 
    [b, a] = butter(3, cutoffFrequency/(Fs/2), 'low'); % 3rd order Butterworth filter

    % Filtering the data
    baseline = filtfilt(b, a, trace); % Zero-phase digital filtering

    % Plot the results
    figure;
    plot(t, trace, 'b', 'DisplayName', 'Original Signal');
    hold on;
    plot(t, baseline, 'r', 'LineWidth', 2, 'DisplayName', 'Extracted Baseline');
    legend('show');
    xlabel('Time (s)');
    ylabel('Intensity');
    title('Action Potential Signal and Extracted Baseline');
    grid on;
end
