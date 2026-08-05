function entries = load_record_results(result_dirs, cycle_names, roles, progress)
%LOAD_RECORD_RESULTS Load AP/Dual Cycle outputs into one role-based schema.
% Saved result contents, rather than result-folder or MAT-file names alone,
% determine which channel and stimulus payloads are used.

if nargin<4,progress=[];end

result_dirs = string(result_dirs(:));
cycle_names = string(cycle_names(:));
roles = lower(string(roles(:)));
if numel(result_dirs) ~= numel(cycle_names) || isempty(result_dirs)
    error('AAA:Record:InvalidCycleSources', ...
        'result_dirs and cycle_names must be nonempty vectors of equal length.');
end
if isempty(roles) || any(~ismember(roles,["voltage","calcium"])) ...
        || numel(unique(roles)) ~= numel(roles)
    error('AAA:Record:InvalidRoles', ...
        'roles must contain unique voltage and/or calcium values.');
end

template = struct('cycle_name',"",'result_dir',"",'channels',struct(), ...
    'voltage',struct(),'calcium',struct(),'stim_results',struct(), ...
    'stim_windows',struct(),'grating_tuning',struct(), ...
    'stim_response',struct(),'stim_metrics',struct(), ...
    'stim_analysis_kind',"",'voltage_response',struct(), ...
    'alignment',struct(),'source_files',struct());
entries = repmat(template,numel(result_dirs),1);
for idx = 1:numel(result_dirs)
    folder = result_dirs(idx);
    report_progress(progress,sprintf('Loading Cycle result %d/%d: %s.', ...
        idx,numel(result_dirs),cycle_names(idx)));
    if ~isfolder(folder)
        error('AAA:Record:MissingCycleResult','Cycle result folder does not exist: %s',folder);
    end
    mode="ap";if numel(roles)>1,mode="dual";end
    [loaded,compatibility_report]=load_analysis_results(folder,struct( ...
        'mode',mode,'scope',"single",'roles',roles,'progress',progress));
    entry = template;
    entry.cycle_name = cycle_names(idx);
    entry.result_dir = folder;
    for role = roles'
        if ~isfield(loaded,char(role))
            error('AAA:Record:MissingChannelResult', ...
                'No semantic %s result was found under %s.',role,folder);
        end
        saved=loaded.(char(role));
        source_file=compatibility_report.absolute_source_files;
        channel = normalize_channel(saved,role,folder);
        entry.channels.(char(role)) = channel;
        entry.(char(role)) = channel; % frozen Dual_rec compatibility alias
        entry.source_files.(char(role)) = source_file;
    end
    if isfield(loaded,'stim')
        stim=loaded.stim;
        entry.stim_results = stim;
        entry.source_files.stim = compatibility_report.absolute_source_files;
        entry = attach_stimulus(entry,stim);
    end
    entry.alignment = build_alignment(entry,roles);
    entries(idx) = entry;
end
end

function report_progress(progress,message)
if ~isa(progress,'function_handle'),return;end
try,progress(string(message));catch,end
end

function channel = normalize_channel(results,role,folder)
[results,source_schema] = normalize_result_schema(results,role,folder);
channel = struct();
channel.raw = extract_stage(results,stage_candidates(role,"raw"));
channel.sensitivity = extract_stage(results,stage_candidates(role,"sensitivity"));
channel.snr = extract_stage(results,stage_candidates(role,"snr"));
if role == "voltage"
    channel.peaks = extract_peaks(results,folder);
end
channel.movie_info = results.movie_info;
channel.source_schema = source_schema;
end

function [results,source_schema] = normalize_result_schema(results,role,folder)
source_schema = "AAA channel results";
if ~isfield(results,'trace_results') || ~isstruct(results.trace_results)
    if has_stage_data(results,'raw'),results=struct('trace_results',results);
    else,error('AAA:Record:MissingTraceStages', ...
            'Saved result has no trace_results payload: %s.',folder);end
end
if isfield(results,'movie_info') && isstruct(results.movie_info) ...
        && isfield(results.movie_info,'frame_rate')
    validate_saved_frame_rate(double(results.movie_info.frame_rate));return;
end

raw_names = stage_candidates(role,"raw");
raw = struct();
for idx = 1:numel(raw_names)
    name = raw_names{idx};
    if isfield(results.trace_results,name) ...
            && isstruct(results.trace_results.(name)) ...
            && isfield(results.trace_results.(name),'data')
        raw = results.trace_results.(name);
        break;
    end
end
if isempty(fieldnames(raw)) || ~isfield(raw,'frame_rate')
    error('AAA:Record:MissingFrameRate', ...
        'Saved trace-stage result is missing raw.frame_rate: %s.',folder);
end
frame_rate = double(raw.frame_rate);
validate_saved_frame_rate(frame_rate);
movie_info = struct('frame_rate',frame_rate, ...
    'frame_count',size(raw.data,1),'role',string(role));
results.movie_info=movie_info;
source_schema = "AAA shared trace stages";
end

function tf=has_stage_data(value,name)
tf=isfield(value,name)&&isstruct(value.(name))&&isfield(value.(name),'data');
end

function validate_saved_frame_rate(frame_rate)
if ~isscalar(frame_rate) || ~isfinite(frame_rate) || frame_rate <= 0
    error('AAA:Record:InvalidFrameRate','Saved frame rate must be positive.');
end
end

function names = stage_candidates(role,metric)
if role == "calcium"
    switch metric
        case "raw", names = {'raw_smoothed','raw'};
        case "sensitivity", names = {'sensitivity_smoothed','sensitivity'};
        case "snr", names = {'snr_smoothed','snr'};
    end
else
    names = {char(metric)};
end
end

function stage = extract_stage(results,names)
if ~isfield(results,'movie_info') || ~isstruct(results.movie_info) ...
        || ~isfield(results.movie_info,'frame_rate')
    error('AAA:Record:MissingFrameRate', ...
        'Saved channel result is missing movie_info.frame_rate.');
end
for idx = 1:numel(names)
    name = names{idx};
    if isfield(results,'trace_results') ...
            && isfield(results.trace_results,name) ...
            && isfield(results.trace_results.(name),'data')
        data = double(results.trace_results.(name).data);
        frame_rate = double(results.movie_info.frame_rate);
        validate_saved_frame_rate(frame_rate);
        stage = struct('stage_name',string(name),'data',data, ...
            'frame_rate',frame_rate,'time',(1:size(data,1))'/frame_rate);
        return;
    end
end
error('AAA:Record:MissingTraceStage','Required stage is missing. Tried: %s.', ...
    strjoin(string(names),', '));
end

function peaks = extract_peaks(results,folder) %#ok<INUSD>
peaks = struct('status',"missing",'source_file',"",'index',{{}}, ...
    'amplitude',{{}},'polarity',{{}});
peak_results = struct();
if isfield(results,'peak_results') && isstruct(results.peak_results)
    peak_results = results.peak_results;
end
if ~isfield(peak_results,'accepted_for_events') ...
        || ~isfield(peak_results.accepted_for_events,'data'), return; end
accepted = peak_results.accepted_for_events.data;
if ~isfield(accepted,'index') || ~iscell(accepted.index), return; end
peaks.index = accepted.index;
if isfield(accepted,'amplitude') && iscell(accepted.amplitude)
    peaks.amplitude = accepted.amplitude;
else
    peaks.amplitude = cell(size(accepted.index));
end
if isfield(accepted,'polarity') && iscell(accepted.polarity)
    peaks.polarity = accepted.polarity;
else
    peaks.polarity = cell(size(accepted.index));
end
peaks.status = "loaded";
end

function entry = attach_stimulus(entry,stim)
if isfield(stim,'analysis_kind'), entry.stim_analysis_kind = string(stim.analysis_kind); end
if isfield(stim,'windows'), entry.stim_windows = stim.windows; end
if isfield(stim,'tuning'), entry.grating_tuning = stim.tuning; end
if isfield(stim,'response'), entry.stim_response = stim.response; end
if isfield(stim,'metrics'), entry.stim_metrics = stim.metrics; end
if isfield(stim,'voltage_response'), entry.voltage_response = stim.voltage_response; end

% Current AAA stores mean-sensitivity trials under metrics, while older Dual
% results may store them under response.sensitivity. Normalize both forms.
if isfield(entry.grating_tuning,'voltage') ...
        && ~isfield(entry.grating_tuning,'voltage_sensitivity') ...
        && isfield(entry.stim_windows,'orientations')
    metric = struct();
    if isfield(entry.stim_metrics,'sensitivity') ...
            && isfield(entry.stim_metrics.sensitivity,'voltage')
        metric = entry.stim_metrics.sensitivity.voltage;
    elseif isfield(entry.stim_response,'sensitivity') ...
            && isfield(entry.stim_response.sensitivity,'voltage')
        metric = entry.stim_response.sensitivity.voltage;
    end
    if isfield(metric,'stim_mean') && isfield(metric,'baseline_mean')
        entry.grating_tuning.voltage_sensitivity = analyze_grating_tuning( ...
            metric.stim_mean,metric.baseline_mean,entry.stim_windows.orientations);
    end
end
end

function alignment = build_alignment(entry,roles)
alignment = struct();
for role = roles'
    field = char(role + "_flash_onset_time");
    alignment.(field) = resolve_flash_onset(entry.stim_windows,role);
end
end

function onset = resolve_flash_onset(windows,role)
onset = NaN;
if ~isstruct(windows) || ~isfield(windows,char(role)), return; end
channel = windows.(char(role));
if isfield(channel,'flash_time_ranges') && ~isempty(channel.flash_time_ranges)
    onset = double(channel.flash_time_ranges(1,1));
elseif isfield(channel,'stim_time_ranges') && ~isempty(channel.stim_time_ranges)
    onset = double(channel.stim_time_ranges(1,1));
end
end
