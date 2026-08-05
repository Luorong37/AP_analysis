function record_average = build_record_average(entries,spec)
%BUILD_RECORD_AVERAGE Build AP/Dual Record traces with frozen Dual rules.
% Flash cycles align to the first flash onset. Grating cycles align to the
% recording start. Only the common overlapping interval is retained and
% each Cycle is linearly interpolated onto the first Cycle frame grid.

if nargin < 2 || ~isstruct(spec) || ~isscalar(spec), spec = struct(); end
roles = lower(string(get_field(spec,'roles',infer_roles(entries))));
roles = roles(:);
if isempty(entries) || isempty(roles)
    error('AAA:Record:MissingInput','Cycle entries and roles are required.');
end

stim_groups = strings(numel(entries),1);
for idx = 1:numel(entries), stim_groups(idx) = classify_stim(entries(idx)); end
if any(stim_groups == "unsupported")
    names = string({entries(stim_groups=="unsupported").cycle_name})';
    error('AAA:Record:UnsupportedStimulus', ...
        'Cannot identify flash/grating stimulus type for cycles: %s.',strjoin(names,', '));
end
types = unique(stim_groups);
if numel(types) ~= 1
    error('AAA:Record:MixedStimulus', ...
        'One Record must contain one stimulus type. Found: %s.',strjoin(types,', '));
end
stim_type = types(1);
alignment_mode = choose(stim_type=="flash","flash_onset","recording_start");

record_average = struct();
record_average.info = struct( ...
    'mode',lower(string(get_field(spec,'mode',"dual"))), ...
    'record_path',string(get_field(spec,'record_path',"")), ...
    'reference_roi_file',string(get_field(spec,'reference_roi_file',"")), ...
    'cycle_names',string({entries.cycle_name})', ...
    'result_dirs',string({entries.result_dir})', ...
    'cycle_count',numel(entries), ...
    'roles',roles, ...
    'stim_group_by_cycle',stim_groups, ...
    'stim_type',stim_type, ...
    'alignment_mode',alignment_mode, ...
    'created_at',datetime("now"), ...
    'alignment_rule',"flash onset for flash; recording start for grating", ...
    'time_range_rule',"common overlapping relative-time interval", ...
    'interpolation_rule',"linear; NaN outside source support", ...
    'frame_rate_rule',"saved movie_info.frame_rate; mismatches are errors", ...
    'roi_rule',"identical ROI count required across included cycles");
record_average.info.saved_frame_rates = struct();
record_average.info.source_movie_files = struct();
record_average.info.source_part_count = struct();
for role = roles'
    role_name = char(role);
    rates = NaN(numel(entries),1);
    source_files = cell(numel(entries),1);
    part_counts = NaN(numel(entries),1);
    for entry_idx = 1:numel(entries)
        info = entries(entry_idx).channels.(role_name).movie_info;
        rates(entry_idx) = double(info.frame_rate);
        if isfield(info,'source_files')
            source_files{entry_idx} = string(info.source_files(:));
            part_counts(entry_idx) = numel(source_files{entry_idx});
        elseif isfield(info,'source_path')
            source_files{entry_idx} = string(info.source_path);
            part_counts(entry_idx) = 1;
        else
            source_files{entry_idx} = strings(0,1);
            part_counts(entry_idx) = 0;
        end
    end
    record_average.info.saved_frame_rates.(role_name) = rates;
    record_average.info.source_movie_files.(role_name) = source_files;
    record_average.info.source_part_count.(role_name) = part_counts;
end
copy_names = {'voltage_polarity','calcium_polarity','calcium_smoothing_window'};
for idx = 1:numel(copy_names)
    name = copy_names{idx};
    if isfield(spec,name), record_average.info.(name) = spec.(name); end
end

record_average.channels = struct();
for role = roles'
    channel = struct();
    for stage = ["raw","sensitivity","snr"]
        channel.(char(stage)) = average_stage(entries,role,stage,alignment_mode);
    end
    record_average.channels.(char(role)) = channel;
    record_average.(char(role)) = channel; % frozen Dual_rec compatibility alias
end
record_average.stim_windows = resolve_average_windows(entries,record_average.channels,roles,alignment_mode);
record_average.visualizations = struct('status',"not_generated");
record_average.time_frequency = struct('status',"not_generated");
record_average.grating_tuning = struct('status',"not_generated");
end

function averaged = average_stage(entries,role,stage_name,alignment_mode)
sample = entries(1).channels.(char(role)).(char(stage_name));
nrois = size(sample.data,2);
frame_rate = sample.frame_rate;
stage_labels = strings(numel(entries),1);
relative_starts = NaN(numel(entries),1);
relative_ends = NaN(numel(entries),1);
for idx = 1:numel(entries)
    current = entries(idx).channels.(char(role)).(char(stage_name));
    if size(current.data,2) ~= nrois
        error('AAA:Record:ROICountMismatch', ...
            'ROI count mismatch while averaging %s %s across cycles.',role,stage_name);
    end
    if abs(current.frame_rate-frame_rate) > 1e-9
        error('AAA:Record:FrameRateMismatch', ...
            'Frame rate mismatch while averaging %s %s across cycles.',role,stage_name);
    end
    stage_labels(idx) = current.stage_name;
    anchor = alignment_time(entries(idx),role,alignment_mode);
    relative = current.time(:)-anchor;
    relative_starts(idx) = relative(1);
    relative_ends(idx) = relative(end);
end
common_start = max(relative_starts);
common_end = min(relative_ends);
if ~isfinite(common_start) || ~isfinite(common_end) || common_end <= common_start
    error('AAA:Record:NoCommonTime', ...
        'No overlapping time range remains after %s alignment for %s %s.', ...
        alignment_mode,role,stage_name);
end
dt = 1/frame_rate;
common_time = (common_start:dt:common_end)';
if numel(common_time)<2, common_time=[common_start;common_end]; end
per_cycle = NaN(numel(common_time),nrois,numel(entries));
anchors = NaN(numel(entries),1);
for idx = 1:numel(entries)
    current = entries(idx).channels.(char(role)).(char(stage_name));
    anchors(idx) = alignment_time(entries(idx),role,alignment_mode);
    relative = current.time(:)-anchors(idx);
    for roi = 1:nrois
        per_cycle(:,roi,idx) = interp1(relative,current.data(:,roi), ...
            common_time,'linear',NaN);
    end
end
cycle_roi_mean = squeeze(mean(per_cycle,2,'omitnan'));
cycle_roi_mean = reshape(cycle_roi_mean,numel(common_time),numel(entries));
averaged = struct('stage_name',stage_labels(1), ...
    'stage_name_by_cycle',stage_labels,'frame_rate',frame_rate, ...
    'time',common_time,'per_cycle',per_cycle, ...
    'cycle_roi_mean',cycle_roi_mean, ...
    'average',mean(per_cycle,3,'omitnan'), ...
    'cycle_names',string({entries.cycle_name})', ...
    'alignment_mode',string(alignment_mode), ...
    'alignment_time_by_cycle',anchors, ...
    'ncycles',numel(entries),'nrois',nrois);
end

function windows = resolve_average_windows(entries,channels,roles,alignment_mode)
windows = struct(); source_idx = [];
for idx = 1:numel(entries)
    candidate = entries(idx).stim_windows;
    if isstruct(candidate) && (~isfield(candidate,'supported') || candidate.supported)
        if ~isempty(fieldnames(candidate)), windows=candidate; source_idx=idx; break; end
    end
end
if isempty(source_idx), return; end
for role = roles'
    role_name = char(role);
    if ~isfield(windows,role_name), continue; end
    if alignment_mode == "recording_start"
        anchor = entries(source_idx).channels.(role_name).raw.time(1);
    else
        anchor = alignment_time(entries(source_idx),role,"flash_onset");
    end
    windows.(role_name) = shift_channel_windows( ...
        windows.(role_name),anchor,channels.(role_name).raw.time);
    if isfield(windows,'flash_windows')
        windows.(role_name).flash_windows = windows.flash_windows;
    end
end
end

function value = shift_channel_windows(value,anchor,common_time)
if ~isfinite(anchor)
    error('AAA:Record:MissingFlashOnset', ...
        'Flash onset is missing, so Record alignment cannot be built.');
end
range_fields = {'stim_time_ranges','flash_time_ranges','baseline_time_ranges', ...
    'response_time_ranges','nonstim_time_ranges'};
for idx = 1:numel(range_fields)
    name = range_fields{idx};
    if isfield(value,name) && isnumeric(value.(name))
        value.(name) = double(value.(name))-anchor;
    end
end
index_fields = {'stim_frame_ranges','flash_frame_ranges','baseline_frame_ranges', ...
    'response_frame_ranges','nonstim_frame_ranges'};
for idx = 1:numel(index_fields)
    name = index_fields{idx};
    if ~isfield(value,name), continue; end
    time_name = strrep(name,'frame','time');
    if ~isfield(value,time_name), continue; end
    ranges = value.(time_name);
    indices = NaN(size(ranges));
    for row = 1:size(ranges,1)
        [~,indices(row,1)] = min(abs(common_time-ranges(row,1)));
        [~,indices(row,2)] = min(abs(common_time-ranges(row,2)));
    end
    value.(name) = indices;
end
end

function anchor = alignment_time(entry,role,mode)
if mode == "recording_start"
    anchor = entry.channels.(char(role)).raw.time(1);
    return;
end
field = char(role+"_flash_onset_time");
anchor = entry.alignment.(field);
if ~isfinite(anchor)
    error('AAA:Record:MissingFlashOnset', ...
        'Flash onset is missing for %s in cycle %s.',role,entry.cycle_name);
end
end

function group = classify_stim(entry)
group = "unsupported";
if isstruct(entry.grating_tuning) && ~isempty(fieldnames(entry.grating_tuning))
    group = "grating"; return;
end
windows = entry.stim_windows;
if ~isstruct(windows), return; end
if isfield(windows,'is_grating') && logical(windows.is_grating)
    group = "grating";
elseif isfield(windows,'stim_type') ...
        && contains(lower(string(windows.stim_type)),"flash")
    group = "flash";
end
end

function roles = infer_roles(entries)
roles = string(fieldnames(entries(1).channels));
end
function value = get_field(s,name,fallback)
if isfield(s,name),value=s.(name);else,value=fallback;end
end
function value = choose(condition,a,b)
if condition,value=a;else,value=b;end
end
