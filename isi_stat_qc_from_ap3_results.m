function isi_stat_qc_summary = isi_stat_qc_from_ap3_results(results_path, params_override)
% isi_stat_qc_from_ap3_results
% Plot-only ISI statistic / firing phenotype QC for saved AP_analysis3 output.
%
% Usage:
%   isi_stat_qc_from_ap3_results('E:\path\to\AP3_result_folder')
%   isi_stat_qc_from_ap3_results('E:\path\to\AP3_result_folder', params_override)
%
% The function loads saved AP3 outputs from results_path, makes QC figures,
% and does not save a .mat result. "Saved parameters" here means the raw AP3
% context needed by this section, such as frame rate, frame count, traces,
% accepted peaks, and previously computed trend FR. ISI-stat thresholds are
% local debug parameters and can be changed below or through params_override.

if nargin < 1 || isempty(results_path)
    results_path = uigetdir(pwd, 'Select AP_analysis3 result folder');
    if isequal(results_path, 0)
        error('isi_stat_qc_from_ap3_results:NoPath', 'No result folder was selected.');
    end
end
if nargin < 2 || isempty(params_override)
    params_override = struct();
end

results_path = char(results_path);
if ~isfolder(results_path)
    error('isi_stat_qc_from_ap3_results:MissingFolder', ...
        'Result folder does not exist: %s', results_path);
end

%% Parameters
% Algorithm overview:
%   1. Use the final accepted spike indices as spike timestamps in seconds.
%   2. Compute basic firing and ISI metrics for each ROI.
%   3. Estimate firing-rate drift from binned FR and first/second-half FR.
%   4. Detect bursts as continuous runs of short ISIs.
%   5. Compute burst burden, burst periodicity, and non-burst local ISI regularity.
%   6. Assign a tentative phenotype label from the combined metrics.
% These labels are only QC-level rules; they are not fixed biological classes.
isi_stat_params = default_isi_stat_params();
isi_stat_params = merge_struct_fields(isi_stat_params, params_override);

%% Load AP3 context
explicit_results_file = fullfile(results_path, '-1_explicit_results.mat');
result_files = struct();
if isfile(explicit_results_file)
    explicit_data = load(explicit_results_file, 'results_summary');
    if isfield(explicit_data, 'results_summary') && ...
            isfield(explicit_data.results_summary, 'result_files')
        result_files = explicit_data.results_summary.result_files;
    end
end

movie_info_file = resolve_result_file(results_path, result_files, 'movie_info', 'movie_info.mat');
trace_results_file = resolve_result_file(results_path, result_files, 'trace_results', 'trace_results.mat');
peak_results_file = resolve_result_file(results_path, result_files, 'peak_results', 'peak_results.mat');

if ~isfile(movie_info_file)
    error('isi_stat_qc_from_ap3_results:MissingMovieInfo', ...
        'movie_info.mat was not found in %s.', results_path);
end
if ~isfile(trace_results_file)
    error('isi_stat_qc_from_ap3_results:MissingTraceResults', ...
        'trace_results.mat was not found in %s.', results_path);
end

s = load(movie_info_file, 'movie_info');
movie_info = s.movie_info;
freq = double(movie_info.frame_rate);
nframes = double(movie_info.frame_count);
recording_duration_s = nframes / freq;

s = load(trace_results_file, 'trace_results');
trace_results = s.trace_results;

peak_results = struct();
if isfile(peak_results_file)
    s = load(peak_results_file, 'peak_results');
    peak_results = s.peak_results;
end

[isi_stat_peak_indices, peak_source] = load_saved_peak_indices(results_path, peak_results);
[isi_stat_trace_for_qc, trace_source] = load_saved_trace_for_isi_qc(results_path, trace_results);

nrois = max([numel(isi_stat_peak_indices), size(isi_stat_trace_for_qc, 2)]);
if nrois == 0
    error('isi_stat_qc_from_ap3_results:NoRois', ...
        'No ROI traces or peak indices were found in %s.', results_path);
end

if numel(isi_stat_peak_indices) < nrois
    isi_stat_peak_indices(numel(isi_stat_peak_indices)+1:nrois) = {[]};
end

all_avgFR = [];
trendx = [];
trend_file = fullfile(results_path, 'Trend Analysis.mat');
if isfile(trend_file)
    trend_data = load(trend_file, 'all_avgFR', 'trendx');
    if isfield(trend_data, 'all_avgFR')
        all_avgFR = trend_data.all_avgFR;
    end
    if isfield(trend_data, 'trendx')
        trendx = trend_data.trendx;
    end
    if isempty(trendx) && ~isempty(all_avgFR) && isi_stat_params.reconstruct_missing_trendx
        trendbin_s = 5;
        if size(all_avgFR, 2) == floor(floor(recording_duration_s) / trendbin_s)
            trendx = (1:size(all_avgFR, 2)) * trendbin_s;
        else
            trendx = linspace(recording_duration_s / size(all_avgFR, 2), ...
                recording_duration_s, size(all_avgFR, 2));
        end
    elseif isempty(trendx) && ~isempty(all_avgFR) && isi_stat_params.reuse_existing_trend_FR
        error('isi_stat_qc_from_ap3_results:MissingTrendX', ...
            'Trend Analysis.mat contains all_avgFR but not trendx. Rerun AP3 after the updated Trend Analysis save, or enable reconstruct_missing_trendx.');
    end
end

%% Convert accepted peaks to spike timestamps
% Convert accepted spike frames into seconds. Non-finite and non-positive
% indices are removed, and duplicate frame indices are collapsed with unique.
spike_times_s = cell(1, nrois);
for roi_idx = 1:nrois
    idx = isi_stat_peak_indices{roi_idx};
    idx = idx(isfinite(idx) & idx > 0);
    spike_times_s{roi_idx} = unique(idx(:)') ./ freq;
end

isi_stat_qc_output_dir = fullfile(results_path, '10_isi_stat_phenotype_qc');
if ~isfolder(isi_stat_qc_output_dir)
    mkdir(isi_stat_qc_output_dir);
end

isi_stat_qc_summary = struct();
isi_stat_qc_summary.neuron_id = (1:nrois)';
isi_stat_qc_summary.spike_count = zeros(nrois, 1);
isi_stat_qc_summary.mean_firing_rate_hz = NaN(nrois, 1);
isi_stat_qc_summary.burst_spike_fraction = NaN(nrois, 1);
isi_stat_qc_summary.LvR_nonburst = NaN(nrois, 1);
isi_stat_qc_summary.CV_IBI = NaN(nrois, 1);
isi_stat_qc_summary.phenotype_label = strings(nrois, 1);
isi_stat_qc_summary.peak_source = string(peak_source);
isi_stat_qc_summary.trace_source = string(trace_source);

can_reuse_trend_FR = isi_stat_params.reuse_existing_trend_FR && ...
    ~isempty(all_avgFR) && ~isempty(trendx) && ...
    size(all_avgFR, 1) >= nrois && size(all_avgFR, 2) == numel(trendx);

%% Per-neuron ISI statistics and QC plots
for roi_idx = 1:nrois
    current_spike_times_s = sort(spike_times_s{roi_idx}(:)');
    spike_count = numel(current_spike_times_s);
    ISI_s = diff(current_spike_times_s);
    ISI_ms = ISI_s * 1000;

    % Basic metrics:
    % spike_count counts accepted spikes.
    % mean_firing_rate_hz is total spikes divided by full recording duration.
    % ISI_s is the interval between adjacent accepted spike timestamps.
    % median_ISI_ms summarizes the central adjacent-spike interval.
    isi_metrics = struct();
    isi_metrics.neuron_id = roi_idx;
    isi_metrics.spike_times_s = current_spike_times_s;
    isi_metrics.recording_duration_s = recording_duration_s;
    isi_metrics.spike_count = spike_count;
    isi_metrics.mean_firing_rate_hz = spike_count / recording_duration_s;
    isi_metrics.ISI_s = ISI_s;
    isi_metrics.median_ISI_ms = median(ISI_ms, 'omitnan');

    % Firing-rate drift:
    % Reuse the earlier Trend Analysis FR bins when dimensions match;
    % otherwise compute fixed-width FR bins here. FR_slope is the linear
    % slope of binned FR versus time in Hz/s. The first/second-half FR values
    % are independent coarse checks for slow state changes.
    if can_reuse_trend_FR
        isi_metrics.window_centers_s = trendx(:)';
        isi_metrics.window_FR_hz = all_avgFR(roi_idx, :);
        isi_metrics.FR_source = "reused_trend_analysis_all_avgFR";
    else
        window_starts = 0:isi_stat_params.window_step_s:(recording_duration_s - isi_stat_params.window_s);
        if isempty(window_starts)
            window_starts = 0;
        end
        isi_metrics.window_centers_s = window_starts + isi_stat_params.window_s / 2;
        isi_metrics.window_FR_hz = NaN(size(isi_metrics.window_centers_s));
        for window_idx = 1:numel(window_starts)
            in_window = current_spike_times_s >= window_starts(window_idx) & ...
                current_spike_times_s < window_starts(window_idx) + isi_stat_params.window_s;
            isi_metrics.window_FR_hz(window_idx) = sum(in_window) / isi_stat_params.window_s;
        end
        isi_metrics.FR_source = "computed_in_isi_stat_section";
    end
    if numel(isi_metrics.window_centers_s) >= 2 && sum(isfinite(isi_metrics.window_FR_hz)) >= 2
        fit_idx = isfinite(isi_metrics.window_FR_hz);
        p = polyfit(isi_metrics.window_centers_s(fit_idx), isi_metrics.window_FR_hz(fit_idx), 1);
        isi_metrics.FR_slope = p(1);
    else
        isi_metrics.FR_slope = NaN;
    end
    first_half = current_spike_times_s(current_spike_times_s <= recording_duration_s / 2);
    second_half = current_spike_times_s(current_spike_times_s > recording_duration_s / 2);
    isi_metrics.first_half_FR = numel(first_half) / max(recording_duration_s / 2, eps);
    isi_metrics.second_half_FR = numel(second_half) / max(recording_duration_s / 2, eps);
    isi_metrics.second_minus_first_FR = isi_metrics.second_half_FR - isi_metrics.first_half_FR;

    % Burst detection:
    % Mark each ISI shorter than burst_isi_threshold_ms. A continuous run of
    % short ISIs maps to one spike block: ISI k..m corresponds to spikes
    % k..m+1. Keeping each continuous run as one event prevents splitting one
    % burst into several smaller bursts.
    isi_bursts = struct('burst_onset_s', {}, 'burst_offset_s', {}, 'burst_duration_ms', {}, ...
        'spike_indices_in_burst', {}, 'spikes_per_burst', {}, 'intra_burst_ISI_ms', {});
    if spike_count >= isi_stat_params.min_spikes_per_burst
        short_isi = ISI_ms < isi_stat_params.burst_isi_threshold_ms;
        run_start = 1;
        burst_id = 0;
        while run_start <= numel(short_isi)
            if ~short_isi(run_start)
                run_start = run_start + 1;
                continue;
            end
            run_end = run_start;
            while run_end < numel(short_isi) && short_isi(run_end + 1)
                run_end = run_end + 1;
            end
            spike_indices = run_start:(run_end + 1);
            if numel(spike_indices) >= isi_stat_params.min_spikes_per_burst
                burst_id = burst_id + 1;
                burst_isi_ms = ISI_ms(run_start:run_end);
                isi_bursts(burst_id).burst_onset_s = current_spike_times_s(spike_indices(1)); %#ok<AGROW>
                isi_bursts(burst_id).burst_offset_s = current_spike_times_s(spike_indices(end));
                isi_bursts(burst_id).burst_duration_ms = ...
                    (isi_bursts(burst_id).burst_offset_s - isi_bursts(burst_id).burst_onset_s) * 1000;
                isi_bursts(burst_id).spike_indices_in_burst = spike_indices;
                isi_bursts(burst_id).spikes_per_burst = numel(spike_indices);
                isi_bursts(burst_id).intra_burst_ISI_ms = burst_isi_ms;
            end
            run_start = run_end + 1;
        end
    end

    % Burst burden and periodicity:
    % burst_spike_fraction is the fraction of all spikes assigned to bursts.
    % IBI_s is measured from burst onset to next burst onset. CV_IBI is only
    % meaningful when several bursts exist; sparse bursts remain descriptive.
    isi_metrics.burst_count = numel(isi_bursts);
    isi_metrics.burst_rate_per_min = isi_metrics.burst_count / recording_duration_s * 60;
    burst_spike_mask = false(1, spike_count);
    for burst_idx = 1:numel(isi_bursts)
        burst_spike_mask(isi_bursts(burst_idx).spike_indices_in_burst) = true;
    end
    isi_metrics.burst_spike_count = sum(burst_spike_mask);
    if spike_count > 0
        isi_metrics.burst_spike_fraction = isi_metrics.burst_spike_count / spike_count;
    else
        isi_metrics.burst_spike_fraction = NaN;
    end
    isi_metrics.median_spikes_per_burst = median([isi_bursts.spikes_per_burst], 'omitnan');
    isi_metrics.median_burst_duration_ms = median([isi_bursts.burst_duration_ms], 'omitnan');
    all_intra_burst_isi_ms = [isi_bursts.intra_burst_ISI_ms];
    isi_metrics.median_intra_burst_ISI_ms = median(all_intra_burst_isi_ms, 'omitnan');
    isi_metrics.median_intra_burst_frequency_hz = 1000 / isi_metrics.median_intra_burst_ISI_ms;
    isi_metrics.burst_onsets_s = [isi_bursts.burst_onset_s];
    isi_metrics.IBI_s = diff(isi_metrics.burst_onsets_s);
    isi_metrics.median_IBI_s = median(isi_metrics.IBI_s, 'omitnan');
    if numel(isi_metrics.IBI_s) >= 2 && mean(isi_metrics.IBI_s, 'omitnan') > 0
        isi_metrics.CV_IBI = std(isi_metrics.IBI_s, 'omitnan') / mean(isi_metrics.IBI_s, 'omitnan');
    else
        isi_metrics.CV_IBI = NaN;
    end

    % Non-burst local regularity:
    % Remove burst spikes first, then compute ISIs only within contiguous
    % non-burst spike-index segments. This avoids creating artificial ISIs
    % that jump across a removed burst. CV_ISI describes global non-burst ISI
    % spread, while LV and LvR describe local adjacent-ISI variability.
    nonburst_indices = find(~burst_spike_mask);
    nonburst_ISI_s = [];
    adjacent_isi_1 = [];
    adjacent_isi_2 = [];
    if numel(nonburst_indices) >= 3
        segment_breaks = [0 find(diff(nonburst_indices) > 1) numel(nonburst_indices)];
        for seg_idx = 1:(numel(segment_breaks) - 1)
            segment_indices = nonburst_indices(segment_breaks(seg_idx)+1:segment_breaks(seg_idx+1));
            segment_times = current_spike_times_s(segment_indices);
            segment_isi = diff(segment_times);
            nonburst_ISI_s = [nonburst_ISI_s segment_isi]; %#ok<AGROW>
            if numel(segment_isi) >= 2
                adjacent_isi_1 = [adjacent_isi_1 segment_isi(1:end-1)]; %#ok<AGROW>
                adjacent_isi_2 = [adjacent_isi_2 segment_isi(2:end)]; %#ok<AGROW>
            end
        end
    end
    isi_metrics.nonburst_ISI_s = nonburst_ISI_s;
    if numel(nonburst_ISI_s) >= 2 && mean(nonburst_ISI_s, 'omitnan') > 0
        isi_metrics.CV_ISI_nonburst = std(nonburst_ISI_s, 'omitnan') / mean(nonburst_ISI_s, 'omitnan');
    else
        isi_metrics.CV_ISI_nonburst = NaN;
    end
    isi_metrics.LV_nonburst = NaN;
    isi_metrics.LvR_nonburst = NaN;
    if ~isempty(adjacent_isi_1)
        isi_sum = adjacent_isi_1 + adjacent_isi_2;
        valid = isi_sum > 0;
        lv_terms = 3 * ((adjacent_isi_2(valid) - adjacent_isi_1(valid)) ./ isi_sum(valid)).^2;
        isi_metrics.LV_nonburst = mean(lv_terms, 'omitnan');
        R = isi_stat_params.LvR_refractory_R_s;
        lvr_terms = 3 * (1 - 4 * R ./ isi_sum(valid)) .* ...
            ((adjacent_isi_2(valid) - adjacent_isi_1(valid)) ./ isi_sum(valid)).^2;
        isi_metrics.LvR_nonburst = mean(lvr_terms, 'omitnan');
    end

    % Tentative phenotype label:
    % The label combines spike count, FR stationarity, burst burden, burst
    % periodicity, and non-burst LvR. It should be read together with the
    % continuous metrics and QC plots, especially before thresholds are tuned.
    if isi_metrics.spike_count < isi_stat_params.min_spike_count
        isi_metrics.phenotype_label = "insufficient_data";
    else
        slope_hz_per_min = isi_metrics.FR_slope * 60;
        half_diff_fraction = abs(isi_metrics.second_minus_first_FR) / max(isi_metrics.mean_firing_rate_hz, eps);
        is_rate_drifting = abs(slope_hz_per_min) >= isi_stat_params.rate_drift_slope_threshold_hz_per_min || ...
            half_diff_fraction >= isi_stat_params.rate_drift_half_diff_fraction;
        is_burst_like = isi_metrics.burst_spike_fraction >= isi_stat_params.high_burst_fraction_threshold;
        if is_rate_drifting && is_burst_like
            isi_metrics.phenotype_label = "state_changing_burst_like";
        elseif is_rate_drifting
            isi_metrics.phenotype_label = "rate_drifting_firing";
        elseif is_burst_like
            if isi_metrics.burst_count >= isi_stat_params.min_burst_count_for_periodicity && ...
                    isfinite(isi_metrics.CV_IBI) && isi_metrics.CV_IBI <= isi_stat_params.low_CV_IBI_threshold
                isi_metrics.phenotype_label = "periodic_bursting";
            elseif isi_metrics.burst_count >= isi_stat_params.min_burst_count_for_periodicity && ...
                    isfinite(isi_metrics.CV_IBI) && isi_metrics.CV_IBI >= isi_stat_params.high_CV_IBI_threshold
                isi_metrics.phenotype_label = "irregular_bursting";
            else
                isi_metrics.phenotype_label = "burst_like_unclassified_periodicity";
            end
        elseif isfinite(isi_metrics.LvR_nonburst) && isi_metrics.LvR_nonburst <= isi_stat_params.low_LvR_threshold
            isi_metrics.phenotype_label = "locally_regular_tonic";
        elseif isfinite(isi_metrics.LvR_nonburst) && isi_metrics.LvR_nonburst >= isi_stat_params.high_LvR_threshold
            isi_metrics.phenotype_label = "locally_irregular_tonic";
        else
            isi_metrics.phenotype_label = "tonic_unclassified_regularity";
        end
    end

    isi_stat_qc_summary.spike_count(roi_idx) = isi_metrics.spike_count;
    isi_stat_qc_summary.mean_firing_rate_hz(roi_idx) = isi_metrics.mean_firing_rate_hz;
    isi_stat_qc_summary.burst_spike_fraction(roi_idx) = isi_metrics.burst_spike_fraction;
    isi_stat_qc_summary.LvR_nonburst(roi_idx) = isi_metrics.LvR_nonburst;
    isi_stat_qc_summary.CV_IBI(roi_idx) = isi_metrics.CV_IBI;
    isi_stat_qc_summary.phenotype_label(roi_idx) = isi_metrics.phenotype_label;

    % QC figure:
    % 1. Trace with spike markers and burst shaded regions.
    % 2. Binned firing rate over time.
    % 3. log10(ISI_ms) histogram with the burst threshold marked.
    % 4. IBI histogram for burst-onset intervals.
    % 5. Spike timeline colored by burst membership.
    fig = figure('Visible', 'on', 'Position', [100, 100, 1100, 1200], ...
        'Color', 'w', 'Name', sprintf('Neuron %03d ISI statistic QC', roi_idx));

    subplot(5,1,1); hold on;
    if ~isempty(isi_stat_trace_for_qc) && size(isi_stat_trace_for_qc, 2) >= roi_idx
        trace_i = isi_stat_trace_for_qc(:, roi_idx);
        t_trace = (1:numel(trace_i)) ./ freq;
        plot(t_trace, trace_i, 'k');
        spike_frames = max(1, min(numel(trace_i), round(isi_metrics.spike_times_s * freq)));
        if ~isempty(spike_frames)
            plot(isi_metrics.spike_times_s, trace_i(spike_frames), 'rv', ...
                'MarkerFaceColor', 'r', 'MarkerSize', 4);
        end
        y_lim = ylim;
        for burst_idx = 1:numel(isi_bursts)
            patch([isi_bursts(burst_idx).burst_onset_s, isi_bursts(burst_idx).burst_offset_s, ...
                isi_bursts(burst_idx).burst_offset_s, isi_bursts(burst_idx).burst_onset_s], ...
                [y_lim(1), y_lim(1), y_lim(2), y_lim(2)], [1.0, 0.85, 0.55], ...
                'FaceAlpha', 0.35, 'EdgeColor', 'none');
        end
        uistack(findobj(gca, 'Type', 'line'), 'top');
    end
    xlim([0, recording_duration_s]);
    title(sprintf('Neuron %03d: %s', roi_idx, char(isi_metrics.phenotype_label)), 'Interpreter', 'none');
    xlabel('Time (s)');
    ylabel(trace_source, 'Interpreter', 'none');
    box off;

    subplot(5,1,2); hold on;
    plot(isi_metrics.window_centers_s, isi_metrics.window_FR_hz, '-o', ...
        'Color', [0.10, 0.35, 0.75], 'MarkerFaceColor', [0.10, 0.35, 0.75]);
    if numel(isi_metrics.window_centers_s) >= 2 && isfinite(isi_metrics.FR_slope)
        fit_line = polyval([isi_metrics.FR_slope, ...
            mean(isi_metrics.window_FR_hz, 'omitnan') - isi_metrics.FR_slope * mean(isi_metrics.window_centers_s, 'omitnan')], ...
            isi_metrics.window_centers_s);
        plot(isi_metrics.window_centers_s, fit_line, '--', 'Color', [0.85, 0.20, 0.20], 'LineWidth', 1.2);
    end
    xlim([0, recording_duration_s]);
    xlabel('Time (s)');
    ylabel('FR (Hz)');
    title(sprintf('Window FR (%s), slope %.4g Hz/s', char(isi_metrics.FR_source), isi_metrics.FR_slope), ...
        'Interpreter', 'none');
    box off;

    subplot(5,1,3); hold on;
    if ~isempty(ISI_ms)
        histogram(log10(ISI_ms), isi_stat_params.isi_histogram_num_bins, ...
            'FaceColor', [0.30, 0.30, 0.30], 'EdgeColor', 'none');
    end
    xline(log10(isi_stat_params.burst_isi_threshold_ms), 'r-', 'LineWidth', 1.5);
    if ~isempty(isi_stat_params.isi_histogram_xlim_log10_ms)
        xlim(isi_stat_params.isi_histogram_xlim_log10_ms);
    end
    xlabel('log10(ISI ms)');
    ylabel('Count');
    title(sprintf('ISI histogram, burst cutoff %.3g ms', isi_stat_params.burst_isi_threshold_ms));
    box off;

    subplot(5,1,4); hold on;
    if ~isempty(isi_metrics.IBI_s)
        histogram(isi_metrics.IBI_s, isi_stat_params.ibi_histogram_num_bins, ...
            'FaceColor', [0.35, 0.55, 0.80], 'EdgeColor', 'none');
    end
    if ~isempty(isi_stat_params.ibi_histogram_xlim_s)
        xlim(isi_stat_params.ibi_histogram_xlim_s);
    end
    xlabel('IBI (s)');
    ylabel('Count');
    title(sprintf('Burst IBI: count %d, CV %.3g', isi_metrics.burst_count, isi_metrics.CV_IBI));
    box off;

    subplot(5,1,5); hold on;
    if ~isempty(current_spike_times_s)
        nonburst_spikes = current_spike_times_s(~burst_spike_mask);
        burst_spikes = current_spike_times_s(burst_spike_mask);
        plot(nonburst_spikes, zeros(size(nonburst_spikes)), 'k|', 'MarkerSize', 10);
        plot(burst_spikes, zeros(size(burst_spikes)), 'r|', 'MarkerSize', 12);
    end
    xlim([0, recording_duration_s]);
    ylim([-1, 1]);
    xlabel('Time (s)');
    yticks([]);
    title(sprintf('Burst fraction %.3g, LvR non-burst %.3g', ...
        isi_metrics.burst_spike_fraction, isi_metrics.LvR_nonburst));
    box off;

    saveas(fig, fullfile(isi_stat_qc_output_dir, sprintf('neuron_%03d_isi_stat_qc.fig', roi_idx)), 'fig');
    saveas(fig, fullfile(isi_stat_qc_output_dir, sprintf('neuron_%03d_isi_stat_qc.png', roi_idx)), 'png');
end

%% Population summary figure
fig = figure('Visible', 'on', 'Position', [100, 100, 900, 800], ...
    'Color', 'w', 'Name', 'Population ISI statistic summary');

subplot(2,1,1); hold on;
scatter(isi_stat_qc_summary.burst_spike_fraction, isi_stat_qc_summary.LvR_nonburst, ...
    45, 'filled', 'MarkerFaceColor', [0.20, 0.45, 0.70], 'MarkerFaceAlpha', 0.75);
xline(isi_stat_params.high_burst_fraction_threshold, 'r--');
yline(isi_stat_params.low_LvR_threshold, 'b--');
yline(isi_stat_params.high_LvR_threshold, 'b--');
xlabel('Burst spike fraction');
ylabel('LvR non-burst');
title('Burst burden vs non-burst local regularity');
box off;

subplot(2,1,2); hold on;
burst_like = isi_stat_qc_summary.burst_spike_fraction >= isi_stat_params.high_burst_fraction_threshold;
scatter(isi_stat_qc_summary.burst_spike_fraction(burst_like), isi_stat_qc_summary.CV_IBI(burst_like), ...
    50, 'filled', 'MarkerFaceColor', [0.80, 0.25, 0.20], 'MarkerFaceAlpha', 0.75);
yline(isi_stat_params.low_CV_IBI_threshold, 'k--');
yline(isi_stat_params.high_CV_IBI_threshold, 'k--');
xlabel('Burst spike fraction');
ylabel('CV IBI');
title('Burst-like neurons: burst burden vs IBI variability');
box off;

saveas(fig, fullfile(isi_stat_qc_output_dir, 'population_isi_stat_summary.fig'), 'fig');
saveas(fig, fullfile(isi_stat_qc_output_dir, 'population_isi_stat_summary.png'), 'png');

fprintf('ISI statistic QC figures saved to %s\n', isi_stat_qc_output_dir);
fprintf('Loaded peaks from %s; loaded trace from %s.\n', peak_source, trace_source);

end

function params = default_isi_stat_params()
params = struct();
params.window_s = 10;                         % FR window length for time-varying rate check.
params.window_step_s = params.window_s;       % FR window step; equal to window_s gives non-overlap bins.
params.reuse_existing_trend_FR = true;        % Reuse Trend Analysis all_avgFR/trendx when available.
params.reconstruct_missing_trendx = true;     % Standalone fallback for old Trend Analysis.mat files.
params.min_spike_count = 30;                  % Below this count, label is insufficient_data.
% burst_isi_threshold_ms means the maximum allowed interval between two
% neighboring spikes inside one burst. With min_spikes_per_burst = 3, a
% burst needs at least two consecutive ISIs below this value.
params.burst_isi_threshold_ms = 1;            % ISI cutoff for burst detection; adjust for sensitivity checks.
params.min_spikes_per_burst = 3;              % Minimum consecutive spikes required for one burst event.
params.min_burst_count_for_periodicity = 5;   % Minimum burst count before CV_IBI is interpreted.
params.isi_histogram_num_bins = 100;          % Number of bins in log10(ISI_ms) histogram.
params.ibi_histogram_num_bins = 20;           % Number of bins in burst IBI_s histogram.
params.isi_histogram_xlim_log10_ms = [];      % Optional log10(ISI_ms) x limit, [] uses auto.
params.ibi_histogram_xlim_s = [];             % Optional IBI histogram x limit in seconds, [] uses auto.
params.rate_drift_slope_threshold_hz_per_min = 0.05; % Absolute FR slope threshold for drift label.
params.rate_drift_half_diff_fraction = 0.5;   % First/second half FR fractional difference threshold.
params.high_burst_fraction_threshold = 0.30;  % Burst-like burden threshold.
params.low_burst_fraction_threshold = 0.10;   % Low burst burden reference threshold.
params.low_CV_IBI_threshold = 0.50;           % Low IBI variability reference for periodic bursting.
params.high_CV_IBI_threshold = 1.00;          % High IBI variability reference for irregular bursting.
params.low_LvR_threshold = 0.80;              % Low non-burst local variability reference.
params.high_LvR_threshold = 1.20;             % High non-burst local variability reference.
params.LvR_refractory_R_s = 0.005;            % Refractory constant used in LvR.
end

function base = merge_struct_fields(base, override)
if ~isstruct(override)
    return;
end
names = fieldnames(override);
for i = 1:numel(names)
    base.(names{i}) = override.(names{i});
end
end

function file_path = resolve_result_file(results_path, result_files, field_name, default_name)
file_path = fullfile(results_path, default_name);
if isstruct(result_files) && isfield(result_files, field_name) && ...
        ~isempty(result_files.(field_name))
    candidate = char(result_files.(field_name));
    if isfile(candidate)
        file_path = candidate;
    elseif isfile(fullfile(results_path, candidate))
        file_path = fullfile(results_path, candidate);
    end
end
end

function [peaks_index, peak_source] = load_saved_peak_indices(results_path, peak_results)
peaks_index = {};
peak_source = 'none';

preferred_stages = {'manually_gated', 'accepted_for_events', ...
    'sensitivity_manually_edited', 'sensitivity_detected', 'volpy_detected'};
for i = 1:numel(preferred_stages)
    stage = preferred_stages{i};
    if isstruct(peak_results) && isfield(peak_results, stage) && ...
            isfield(peak_results.(stage), 'data') && isfield(peak_results.(stage).data, 'index')
        peaks_index = peak_results.(stage).data.index;
        peak_source = ['peak_results.' stage '.data.index'];
        return;
    end
end

candidate_files = { ...
    fullfile(results_path, '5_manual_peak_gate_results.mat'), 'peaks_index_gated'; ...
    fullfile(results_path, '5_accepted_peaks.mat'), 'peaks_index'; ...
    fullfile(results_path, '4_sensitivity_peak_edit_results.mat'), 'peaks_index'; ...
    fullfile(results_path, '4_sensitivity_peak_detect_results.mat'), 'peaks_index'};

for i = 1:size(candidate_files, 1)
    file_i = candidate_files{i, 1};
    var_i = candidate_files{i, 2};
    if isfile(file_i)
        s = load(file_i, var_i);
        if isfield(s, var_i)
            peaks_index = s.(var_i);
            [~, name_i, ext_i] = fileparts(file_i);
            peak_source = [name_i ext_i ':' var_i];
            return;
        end
    end
end

error('isi_stat_qc_from_ap3_results:MissingPeaks', ...
    'No saved peak indices were found in %s.', results_path);
end

function [trace_data, trace_source] = load_saved_trace_for_isi_qc(results_path, trace_results)
trace_data = get_trace_stage_data(trace_results, 'sensitivity', []);
trace_source = 'trace_results.sensitivity';
if isempty(trace_data)
    trace_data = get_trace_stage_data(trace_results, 'snr', []);
    trace_source = 'trace_results.snr';
end

if isempty(trace_data)
    metric_file = fullfile(results_path, '3_metric_results.mat');
    if isfile(metric_file)
        s = load(metric_file, 'traces_sensitivity', 'traces_SNR');
        if isfield(s, 'traces_sensitivity') && ~isempty(s.traces_sensitivity)
            trace_data = s.traces_sensitivity;
            trace_source = '3_metric_results.mat:traces_sensitivity';
        elseif isfield(s, 'traces_SNR') && ~isempty(s.traces_SNR)
            trace_data = s.traces_SNR;
            trace_source = '3_metric_results.mat:traces_SNR';
        end
    end
end

if isempty(trace_data)
    trace_data = [];
    trace_source = 'none';
end
end

function data = get_trace_stage_data(trace_results, stage_name, default_value)
data = default_value;
if isstruct(trace_results) && isfield(trace_results, 'trace_results') && ...
        isfield(trace_results.trace_results, stage_name) && ...
        isfield(trace_results.trace_results.(stage_name), 'data')
    data = trace_results.trace_results.(stage_name).data;
elseif isstruct(trace_results) && isfield(trace_results, stage_name) && ...
        isfield(trace_results.(stage_name), 'data')
    data = trace_results.(stage_name).data;
end
end
