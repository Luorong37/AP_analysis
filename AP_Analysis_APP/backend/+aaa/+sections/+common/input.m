function [ctx, outcome] = input(ctx, action)
%INPUT Load each channel through the shared movie-loading algorithm.
%   Channel profile is the only authority for role, fps, and transpose.
%   Dual channels are spatially matched by the established integer-ratio
%   average pooling rule; their frame counts are never aligned here.

action = normalize_action(action);
outcome = struct( ...
    'section', "input", ...
    'action', action, ...
    'state', "completed", ...
    'channels', strings(0, 1), ...
    'message', "");

if action == "skip"
    outcome.state = "skipped";
    outcome.message = "Input loading skipped.";
    return;
end
validate_context(ctx);
roles = channel_roles(ctx.channels);
outcome.channels = roles(:);

% Load inputs: resolve each logical channel source before reading pixels.
aaa.sections.emit_phase(ctx,"input","resolve_sources","running", ...
    "Resolving movie sources.");
[sources, input_layout] = resolve_sources(ctx, roles);

% Run algorithm: load each movie with its bound role profile.
for idx = 1:numel(ctx.channels)
    channel = ctx.channels(idx);
    role=string(channel.profile.role);
    aaa.sections.emit_phase(ctx,"input","load_movie","running", ...
        sprintf('Loading %s movie.',role),role);
    [movie_3d, movie_info] = load_movie( ...
        sources{idx}, channel.profile);
    aaa.sections.emit_phase(ctx,"input","load_movie","completed", ...
        sprintf('Loaded %s movie: %d frame(s).',role,movie_info.frame_count),role);
    channel.source = sources{idx};
    movie_info.camera_index = channel.camera_index;
    channel.data.movie_3d = movie_3d;
    channel.data.movie_info = movie_info;
    channel.results.movie_info = movie_info;
    if ~isfield(channel.results, 'trace_results')
        channel.results.trace_results = struct();
    end
    channel.profile.roi.frame_size = movie_info.analysis_frame_size;
    channel.data.time = (1:movie_info.frame_count)' / movie_info.frame_rate;
    channel.output_files.input = strings(0,1);
    ctx.channels(idx) = channel;
end

% Return results: register loaded arrays, source inventory, and geometry.
geometry = struct('action',"none",'common_size',[]);
if numel(ctx.channels) == 2
    aaa.sections.emit_phase(ctx,"input","match_geometry","running", ...
        "Matching Dual spatial geometry.");
    voltage_idx = find(roles == "voltage", 1, 'first');
    calcium_idx = find(roles == "calcium", 1, 'first');
    if isempty(voltage_idx) || isempty(calcium_idx)
        error('AAA:Sections:InvalidDualRoles', ...
            'A two-channel input must contain one voltage and one calcium profile.');
    end
    [ctx.channels(voltage_idx), ctx.channels(calcium_idx), geometry] = ...
        match_dual_geometry(ctx.channels(voltage_idx), ctx.channels(calcium_idx));
end
ctx.input_geometry = geometry;
ctx.input_layout = input_layout;
outcome.geometry = geometry;
outcome.input_layout = input_layout;
outcome.message = sprintf('Loaded %d channel movie(s).', numel(ctx.channels));
end

function action = normalize_action(action)
action = lower(strtrim(string(action)));
if ~isscalar(action) || ~ismember(action, ["run","skip"])
    error('AAA:Sections:InvalidInputAction', ...
        'Input action must be run or skip.');
end
end

function validate_context(ctx)
if ~isstruct(ctx) || ~isscalar(ctx) ...
        || ~isfield(ctx, 'channels') || isempty(ctx.channels) ...
        || ~isstruct(ctx.channels)
    error('AAA:Sections:InvalidInputContext', ...
        'Input context must contain a nonempty channels struct array.');
end
required = {'camera_index','source','profile','data','output_files'};
for channel_idx = 1:numel(ctx.channels)
    channel = ctx.channels(channel_idx);
    for field_idx = 1:numel(required)
        if ~isfield(channel, required{field_idx})
            error('AAA:Sections:IncompleteInputChannel', ...
                'Channel %d is missing %s.', channel_idx, required{field_idx});
        end
    end
    if ~isnumeric(channel.camera_index) || ~isscalar(channel.camera_index) ...
            || ~isfinite(channel.camera_index) ...
            || channel.camera_index < 1 ...
            || channel.camera_index ~= round(channel.camera_index)
        error('AAA:Sections:InvalidCameraIndex', ...
            'Channel %d camera_index must be a positive integer.', channel_idx);
    end
end
end

function roles = channel_roles(channels)
roles = strings(1, numel(channels));
for idx = 1:numel(channels)
    profile = channels(idx).profile;
    if ~isstruct(profile) || ~isscalar(profile) ...
            || ~isfield(profile, 'role') || strlength(string(profile.role)) == 0
        error('AAA:Sections:MissingChannelRole', ...
            'Channel %d profile must contain one role.', idx);
    end
    roles(idx) = lower(string(profile.role));
end
if numel(unique(roles)) ~= numel(roles)
    error('AAA:Sections:DuplicateChannelRole', ...
        'Each channel profile role must be unique.');
end
end

function [sources, layout] = resolve_sources(ctx, roles)
sources = cell(1, numel(ctx.channels));
layout = struct( ...
    'mode', "direct_channel_source", ...
    'description', "Each channel uses its bound source directly.");
if isscalar(ctx.channels)
    source_path = extract_source_path(ctx.channels(1).source);
    [sources,report] = aaa.helpers.common.resolve_movie_sources( ...
        source_path,channel_specs(ctx.channels));
    if ~report.valid
        error('AAA:Sections:UnresolvedCameraSource','%s', ...
            strjoin(report.errors,newline));
    end
    layout.mode = "resolved_single_channel_source";
    layout.description = "Resolved manifest/split parts before loading.";
    layout.sources = report.sources;
    return;
end

source_paths = strings(1, numel(ctx.channels));
for idx = 1:numel(ctx.channels)
    source_paths(idx) = extract_source_path(ctx.channels(idx).source);
end
if ~paths_equal(source_paths(1), source_paths(2))
    for idx = 1:numel(ctx.channels)
        [resolved,report] = aaa.helpers.common.resolve_movie_sources( ...
            source_paths(idx),channel_specs(ctx.channels(idx)));
        if ~report.valid
            error('AAA:Sections:UnresolvedCameraSource','%s', ...
                strjoin(report.errors,newline));
        end
        sources{idx} = resolved{1};
    end
    layout.mode = "resolved_explicit_channel_sources";
    layout.description = "Resolved each explicitly bound channel source.";
    layout.sources = vertcat(sources{:});
    return;
end

cycle_path = source_paths(1);
if ~isfolder(cycle_path)
    error('AAA:Sections:UnresolvedDualSource', ...
        ['Dual channels share one source that is not a Cycle/raw folder: %s. ' ...
         'Bind distinct explicit sources or provide a resolvable folder.'], ...
        cycle_path);
end

[sources, found, raw_layout] = resolve_raw_pair( ...
    cycle_path, ctx.channels, roles, ctx.request.params.dual);
if found
    layout = raw_layout;
    return;
end
[sources,report] = aaa.helpers.common.resolve_movie_sources( ...
    cycle_path,channel_specs(ctx.channels));
if ~report.valid
    error('AAA:Sections:UnresolvedCameraSource','%s', ...
        strjoin(report.errors,newline));
end
layout = struct('mode',"resolved_cycle_sources", ...
    'description',"Resolved manifest paths, split parts, then camera folders.", ...
    'sources',report.sources);
end

function specs=channel_specs(channels)
specs=repmat(struct('camera_index',NaN,'role',""),1,numel(channels));
for idx=1:numel(channels)
    specs(idx).camera_index=channels(idx).camera_index;
    specs(idx).role=string(channels(idx).profile.role);
end
end

function [sources, found, layout] = resolve_raw_pair( ...
        cycle_path, channels, roles, mode_params)
sources = cell(1, numel(channels));
found = false;
layout = struct( ...
    'mode', "paired_raw_channel_folders", ...
    'description', "Resolved one primary raw folder and its green sibling.");
required = {'raw_dual_green_suffix','raw_dual_primary_role', ...
    'raw_dual_green_role'};
if ~isstruct(mode_params) || ~all(isfield(mode_params, required))
    return;
end
suffix = string(mode_params.raw_dual_green_suffix);
primary_role = lower(string(mode_params.raw_dual_primary_role));
green_role = lower(string(mode_params.raw_dual_green_role));
if strlength(suffix) == 0 || primary_role == green_role ...
        || ~all(ismember([primary_role,green_role], roles))
    return;
end

[parent_path, folder_name, folder_ext] = fileparts(char(cycle_path));
folder_name = string([folder_name, folder_ext]);
primary_path = "";
green_path = "";
if endsWith(folder_name, suffix, 'IgnoreCase', true)
    green_path = cycle_path;
    primary_name = extractBefore( ...
        folder_name, strlength(folder_name) - strlength(suffix) + 1);
    candidate = string(fullfile(parent_path, primary_name));
    if isfolder(candidate)
        primary_path = candidate;
    end
else
    candidate = string(fullfile(parent_path, folder_name + suffix));
    if isfolder(candidate)
        primary_path = cycle_path;
        green_path = candidate;
    end
end
if strlength(primary_path) == 0 || strlength(green_path) == 0
    return;
end

primary_idx = find(roles == primary_role, 1, 'first');
green_idx = find(roles == green_role, 1, 'first');
sources{primary_idx} = struct( ...
    'path', primary_path, 'label', last_path_part(primary_path));
sources{green_idx} = struct( ...
    'path', green_path, 'label', last_path_part(green_path));
layout.primary_path = primary_path;
layout.green_path = green_path;
layout.green_suffix = suffix;
found = true;
end

function value = extract_source_path(source)
if isstruct(source) && isscalar(source) && isfield(source, 'path')
    value = string(source.path);
elseif (ischar(source) && isrow(source)) ...
        || (isstring(source) && isscalar(source))
    value = string(source);
else
    error('AAA:Sections:InvalidChannelSource', ...
        'Each channel source must be a path or a scalar struct with path.');
end
if strlength(value) == 0
    error('AAA:Sections:InvalidChannelSource', ...
        'Channel source path cannot be empty.');
end
end

function tf = paths_equal(left, right)
if ispc
    tf = strcmpi(char(left), char(right));
else
    tf = left == right;
end
end

function value = last_path_part(input_path)
[~, name, ext] = fileparts(char(input_path));
value = string([name, ext]);
end

function [voltage, calcium, geometry] = match_dual_geometry(voltage, calcium)
voltage_size = spatial_size(voltage.data.movie_3d);
calcium_size = spatial_size(calcium.data.movie_3d);
geometry = struct( ...
    'voltage_before', voltage_size, ...
    'calcium_before', calcium_size, ...
    'action', "none", ...
    'common_size', []);

if isequal(voltage_size, calcium_size)
    geometry.common_size = voltage_size;
else
    ratio_c_to_v = calcium_size ./ voltage_size;
    ratio_v_to_c = voltage_size ./ calcium_size;
    if all(mod(calcium_size, voltage_size) == 0) ...
            && ratio_c_to_v(1) == ratio_c_to_v(2)
        calcium.data.movie_3d = average_pool_movie( ...
            calcium.data.movie_3d, ratio_c_to_v(1));
        geometry.action = "downsample_calcium_to_voltage";
    elseif all(mod(voltage_size, calcium_size) == 0) ...
            && ratio_v_to_c(1) == ratio_v_to_c(2)
        voltage.data.movie_3d = average_pool_movie( ...
            voltage.data.movie_3d, ratio_v_to_c(1));
        geometry.action = "downsample_voltage_to_calcium";
    elseif all(mod(calcium_size, voltage_size) == 0) ...
            || all(mod(voltage_size, calcium_size) == 0)
        error('AAA:Sections:InconsistentGeometryRatio', ...
            'Dual resize ratios must be equal across dimensions.');
    else
        error('AAA:Sections:IncompatibleCameraGeometry', ...
            'Voltage/calcium sizes are incompatible: [%d %d] vs [%d %d].', ...
            voltage_size(1), voltage_size(2), ...
            calcium_size(1), calcium_size(2));
    end
    voltage = refresh_geometry(voltage);
    calcium = refresh_geometry(calcium);
    geometry.common_size = spatial_size(voltage.data.movie_3d);
end
end

function channel = refresh_geometry(channel)
[ncols, nrows, nframes] = size(channel.data.movie_3d);
channel.data.movie_info.analysis_frame_size = [ncols, nrows];
channel.data.movie_info.frame_count = nframes;
channel.data.movie_info.updated_at = datetime('now');
channel.results.movie_info = channel.data.movie_info;
channel.data.time = (1:nframes)' / channel.data.movie_info.frame_rate;
channel.profile.roi.frame_size = [ncols, nrows];
end

function movie_out = average_pool_movie(movie_in, ratio)
ratio = round(ratio);
[ncols, nrows, nframes] = size(movie_in);
movie_out = reshape(movie_in, ...
    ratio, ncols / ratio, ratio, nrows / ratio, nframes);
movie_out = squeeze(mean(mean(movie_out, 1), 3));
movie_out = reshape(movie_out, ncols / ratio, nrows / ratio, nframes);
end

function value = spatial_size(movie_3d)
value = [size(movie_3d, 1), size(movie_3d, 2)];
end
