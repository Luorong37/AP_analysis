function [ctx, outcome] = peak(ctx, action)
%PEAK Run or strictly reuse the shared voltage peak chain.
% Scientific detection/editing is mode-neutral. This section owns only
% result-stage composition and AP/Dual-compatible file names.

action = normalize_action(action);
outcome = struct( ...
    'section', "peak", ...
    'action', action, ...
    'state', "completed", ...
    'channels', "voltage", ...
    'message', "");
if action == "skip"
    outcome.state = "skipped";
    outcome.message = "Peak section skipped.";
    return;
end

% Load inputs: obtain voltage sensitivity or a compatible saved peak result.
[voltage_idx, sensitivity, output_path] = validate_context(ctx);
channel = ctx.channels(voltage_idx);
if action == "reuse"
    aaa.sections.emit_phase(ctx,"peak","reuse_load","running", ...
        "Loading and validating reusable accepted peaks.","voltage");
    [peak_results, source_file, compatibility_report] = load_reusable_peaks( ...
        ctx, size(sensitivity, 2), size(sensitivity, 1));
    peak_results.reuse_info = struct( ...
        'reused', true, ...
        'source_results_path', string(ctx.request.workflow.source_results_path), ...
        'destination_results_path', string(output_path), ...
        'validation', "accepted peak ROI count and frame indices match current sensitivity traces", ...
        'source_file', string(source_file), ...
        'reused_at', datetime("now"));
else
    % Run algorithm: detection and optional manual gating keep their names.
    peak_results = run_peak_chain(sensitivity, channel.profile,ctx);
    source_file = "";
    compatibility_report = struct();
end

% Return results: register the accepted peak chain and provenance.
channel.results.peak_results = peak_results;
aaa.sections.emit_phase(ctx,"peak","persist","running", ...
    "Preparing the unified peak result.","voltage");
channel.output_files.peak = strings(0,1);
ctx.channels(voltage_idx) = channel;

accepted = peak_results.accepted_for_events.data.index;
outcome.peak_count = sum(cellfun(@numel, accepted));
outcome.source_files = string(source_file);
outcome.compatibility_report = compatibility_report;
if action == "reuse"
    outcome.message = sprintf('Reused %d accepted voltage peaks.', ...
        outcome.peak_count);
else
    outcome.message = sprintf('Detected/edited %d accepted voltage peaks.', ...
        outcome.peak_count);
end
end

function action = normalize_action(action)
action = lower(strtrim(string(action)));
if ~isscalar(action) || ~ismember(action, ["run","reuse","skip"])
    error('AAA:Sections:InvalidPeakAction', ...
        'Peak action must be run, reuse, or skip.');
end
end

function [voltage_idx, sensitivity, output_path] = validate_context(ctx)
if ~isstruct(ctx) || ~isscalar(ctx) || ~isfield(ctx, 'channels') ...
        || isempty(ctx.channels) || ~isfield(ctx, 'mode') ...
        || ~ismember(string(ctx.mode), ["ap","dual"])
    error('AAA:Sections:InvalidPeakContext', ...
        'Peak context must be one AP/Dual section context.');
end
roles = strings(1, numel(ctx.channels));
for idx = 1:numel(ctx.channels)
    if ~isfield(ctx.channels(idx), 'profile') ...
            || ~isfield(ctx.channels(idx).profile, 'role')
        error('AAA:Sections:InvalidPeakContext', ...
            'Each channel must contain one bound profile.');
    end
    roles(idx) = lower(string(ctx.channels(idx).profile.role));
end
voltage_idx = find(roles == "voltage", 1);
if isempty(voltage_idx)
    error('AAA:Sections:MissingVoltageChannel', ...
        'Peak section requires one voltage channel.');
end
channel = ctx.channels(voltage_idx);
if ~isfield(channel, 'results') || ~isstruct(channel.results)
    error('AAA:Sections:PeakTraceMissing', ...
        'Voltage results are not available.');
end
sensitivity = fetch_trace_stage(channel.results, 'sensitivity');
output_path = string(ctx.output_path);
if strlength(output_path) == 0
    error('AAA:Sections:MissingOutputPath', ...
        'Peak section requires a nonempty output path.');
end
if ~isfolder(output_path)
    mkdir(output_path);
end
end

function peak_results = run_peak_chain(sensitivity, profile,ctx)
% Handler-local phase coordinator. Scientific atoms detect_peaks and
% edit_peaks receive only data plus the bound profile, never ctx/request.
aaa.sections.emit_phase(ctx,"peak","detect","running", ...
    "Denoising sensitivity traces and detecting peaks.","voltage");
detected = detect_peaks(sensitivity, profile);
peak_results = struct();
peak_results.sensitivity_detected.data = detected;
parameters = profile.peak;
parameters = rmfield(parameters, 'run_manual_edit');
peak_results.sensitivity_detected.info = struct( ...
    'result_name', 'sensitivity_detected', ...
    'parent_results', {{}}, ...
    'trace_result', 'sensitivity', ...
    'method', 'detect_sensitivity_peaks', ...
    'parameters', parameters, ...
    'created_at', datetime("now"));

if profile.peak.run_manual_edit
    aaa.sections.emit_phase(ctx,"peak","manual_gating","waiting", ...
        "Waiting for manual peak gating.","voltage");
else
    aaa.sections.emit_phase(ctx,"peak","manual_gating","skipped", ...
        "Manual peak gating is disabled.","voltage");
end
edited = edit_peaks(sensitivity, detected, profile);
manual_parameters = struct( ...
    'enabled', logical(profile.peak.run_manual_edit), ...
    'min_peak_distance_frames', profile.peak.min_peak_distance_frames, ...
    'frame_rate_hz', profile.frame_rate, ...
    'default_mode', "delete_box", ...
    'keys', "A add box, D delete box, C clear current ROI, F flip current ROI polarity, N/space next ROI, P/left previous ROI, R reset ROI, Q finish");
method = 'manual_add_delete_peak_edit';
if ~profile.peak.run_manual_edit
    method = 'manual_add_delete_peak_edit_skipped';
end
peak_results.sensitivity_manually_edited.data = edited;
peak_results.sensitivity_manually_edited.info = struct( ...
    'result_name', 'sensitivity_manually_edited', ...
    'parent_results', {{'sensitivity_detected'}}, ...
    'trace_result', 'sensitivity', ...
    'method', method, ...
    'parameters', manual_parameters, ...
    'created_at', datetime("now"));
peak_results.accepted_for_events = ...
    peak_results.sensitivity_manually_edited;
peak_results.accepted_for_events.info.result_name = 'accepted_for_events';
peak_results.accepted_for_events.info.parent_results = ...
    {'sensitivity_manually_edited'};
peak_results.accepted_for_events.info.method = ...
    'accepted_peaks_after_manual_edit';
peak_results.current_result = 'accepted_for_events';
end

function [peak_results, source_file, report] = load_reusable_peaks(ctx, nrois, nframes)
source_path = string(ctx.request.workflow.source_results_path);
if strlength(source_path) == 0 || ~isfolder(source_path)
    error('AAA:Sections:PeakReuseSourceMissing', ...
        'peak="reuse" requires workflow.source_results_path.');
end
[results,report]=load_analysis_results(source_path,struct( ...
    'mode',ctx.mode,'scope',ctx.scope,'runtime',ctx.runtime));
if ~isfield(results,'voltage') || ~isfield(results.voltage,'peak_results')
    error('AAA:Sections:PeakReuseResultMissing', ...
        'No reusable voltage peak result was found in %s.',source_path);
end
peak_results=results.voltage.peak_results;
[valid,reason]=validate_peak_result(peak_results,nrois,nframes);
if ~valid
    error('AAA:Sections:InvalidPeakReuse', ...
        'Saved voltage peaks cannot be reused: %s.',reason);
end
source_file="";
if ~isempty(report.absolute_source_files),source_file=report.absolute_source_files(1);end
end
