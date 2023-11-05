function firing_rate_trace = calculate_firing_rate(peak_index, window, num_windows, window_length)
    % Calculate firing rate using a specified window (in seconds)
    % spikes: a vector containing the spike times in seconds
    % total_duration: the total duration of the recording in seconds
    % window: the window for calculating the firing rate in seconds
    
    % Initialize the firing rate vector with zeros
    firing_rate_trace = zeros(num_windows, 1);
    
    % Calculate the firing rate for each window
    for w = 1:num_windows
        % Define the window start and end times
        window_start = (w-1) *  window_length;
        window_end = w * window_length;
        
        % Count the number of spikes in the current window
        spikes_in_window = sum(peak_index >= window_start & peak_index < window_end);
        
        % Calculate the firing rate (spikes per second)
        firing_rate_trace(w) = spikes_in_window ./ window; % since window is 1 second, this is just spikes_in_window
    end
end
