function [sources, report] = resolve_movie_sources(input_path, channels)
%RESOLVE_MOVIE_SOURCES Resolve channel movies without reading movie arrays.
% Manifest paths are preferred. Missing logical TIFF paths may resolve to
% their NAME_stackNN parts. Camera-folder and unique-folder discovery are
% fallbacks. The returned source structs preserve every physical part.

input_path = strtrim(string(input_path));
channels = normalize_channels(channels);
sources = cell(1, numel(channels));
report = struct('valid',true,'input_path',input_path, ...
    'sources',repmat(empty_source(),0,1),'errors',strings(0,1), ...
    'warnings',strings(0,1));

if isfile(input_path)
    if numel(channels) ~= 1 || ~is_supported_movie_file(input_path)
        report.errors(end+1,1) = ...
            "A direct movie file can resolve only one supported channel.";
    else
        sources{1} = make_source(input_path, input_path, input_path, ...
            last_part(input_path), channels(1), "direct_movie_file");
    end
elseif ~isfolder(input_path)
    report.errors(end+1,1) = "Movie input path does not exist: " + input_path;
else
    manifest = read_cycle_manifest(input_path);
    for idx = 1:numel(channels)
        [source, message] = resolve_one(input_path, channels(idx), manifest, ...
            numel(channels));
        if strlength(message) > 0
            report.errors(end+1,1) = message; %#ok<AGROW>
        else
            sources{idx} = source;
        end
    end
end

report.errors = unique(report.errors,'stable');
report.valid = isempty(report.errors);
if report.valid
    report.sources = vertcat(sources{:});
end
end

function channels = normalize_channels(channels)
if ~isstruct(channels) || isempty(channels) ...
        || ~all(isfield(channels,{'camera_index','role'}))
    error('AAA:InputDiscovery:InvalidChannels', ...
        'channels must contain camera_index and role.');
end
channels = channels(:)';
for idx = 1:numel(channels)
    camera_index = double(channels(idx).camera_index);
    role = lower(strtrim(string(channels(idx).role)));
    if ~isscalar(camera_index) || ~isfinite(camera_index) ...
            || camera_index < 1 || camera_index ~= round(camera_index) ...
            || ~isscalar(role) || strlength(role) == 0
        error('AAA:InputDiscovery:InvalidChannels', ...
            'Every channel requires one positive camera index and role.');
    end
    channels(idx).camera_index = camera_index;
    channels(idx).role = role;
end
end

function [source, message] = resolve_one(cycle_path, channel, manifest, channel_count)
source = empty_source();
message = "";
[logical_path, label] = manifest_channel(manifest, channel.camera_index);
if strlength(logical_path) > 0
    [resolved_path, files, method] = resolve_logical_path(logical_path);
    if ~isempty(files)
        if strlength(label) == 0, label = last_part(resolved_path); end
        source = make_source(resolved_path, logical_path, files, label, ...
            channel, method);
        return;
    end
end

listing = dir(fullfile(cycle_path, sprintf('Cam%d_*',channel.camera_index)));
listing = listing(~startsWith({listing.name},'.'));
[source, found, ambiguous] = source_from_listing(listing, channel, label, ...
    "camera_folder_fallback");
if found, return; end
if ambiguous
    message = sprintf('Camera %d has multiple movie folders under %s.', ...
        channel.camera_index, cycle_path);
    return;
end

if channel_count == 1
    [files, resolved_path] = movie_files_in_path(cycle_path);
    if ~isempty(files)
        if strlength(label) == 0, label = last_part(resolved_path); end
        source = make_source(resolved_path, logical_path, files, label, ...
            channel, "cycle_folder_direct");
        return;
    end
    children = dir(cycle_path);
    children = children([children.isdir]);
    children = children(~ismember({children.name},{'.','..'}));
    [source, found, ambiguous] = source_from_listing(children, channel, ...
        label, "unique_movie_folder_fallback");
    if found, return; end
    if ambiguous
        message = "Multiple possible AP movie folders were found under " + cycle_path + ...
            "; camera identity is ambiguous.";
        return;
    end
end
message = sprintf('Cannot resolve camera %d (%s) movie input under %s.', ...
    channel.camera_index, channel.role, cycle_path);
end

function [source, found, ambiguous] = source_from_listing(listing, channel, label, method)
source = empty_source();
found = false;
ambiguous = false;
candidates = repmat(empty_source(),0,1);
for idx = 1:numel(listing)
    path_value = string(fullfile(listing(idx).folder,listing(idx).name));
    [files, resolved_path] = movie_files_in_path(path_value);
    if isempty(files), continue; end
    current_label = label;
    if strlength(current_label) == 0
        current_label = erase(string(listing(idx).name), ...
            sprintf('Cam%d_',channel.camera_index));
    end
    candidates(end+1,1) = make_source(resolved_path, "", files, ...
        current_label, channel, method); %#ok<AGROW>
end
if numel(candidates) == 1
    source = candidates(1);
    found = true;
elseif numel(candidates) > 1
    ambiguous = true;
end
end

function [resolved_path, files, method] = resolve_logical_path(logical_path)
resolved_path = "";
files = strings(0,1);
method = "manifest_exact";
if isfile(logical_path) || isfolder(logical_path)
    [files,resolved_path] = movie_files_in_path(logical_path);
    return;
end
[folder,stem,ext] = fileparts(logical_path);
if ismember(lower(string(ext)),[".tif",".tiff"]) && isfolder(folder)
    listing = dir(fullfile(folder, char(string(stem)+"_stack*"+string(ext))));
    listing = listing(~[listing.isdir]);
    if ~isempty(listing)
        files = sorted_paths(listing);
        resolved_path = string(folder);
        method = "manifest_split_tiff_parts";
    end
end
end

function [files, resolved_path] = movie_files_in_path(path_value)
path_value = string(path_value);
files = strings(0,1);
resolved_path = "";
if isfile(path_value)
    if is_supported_movie_file(path_value)
        files = path_value;
        resolved_path = path_value;
    end
    return;
end
if ~isfolder(path_value), return; end
listing = [dir(fullfile(path_value,'*.tif')); ...
    dir(fullfile(path_value,'*.tiff'))];
listing = listing(~[listing.isdir]);
if ~isempty(listing)
    files = sorted_paths(listing);
    resolved_path = path_value;
    return;
end
listing = dir(fullfile(path_value,'*.bin'));
listing = listing(~[listing.isdir]);
if numel(listing) == 1
    files = string(fullfile(listing.folder,listing.name));
    resolved_path = files;
    return;
end
listing = dir(fullfile(path_value,'*.mat'));
listing = listing(~[listing.isdir]);
for idx = 1:numel(listing)
    candidate = string(fullfile(listing(idx).folder,listing(idx).name));
    if is_supported_movie_file(candidate)
        files(end+1,1) = candidate; %#ok<AGROW>
    end
end
if numel(files) == 1
    resolved_path = files;
else
    files = strings(0,1);
end
end

function source = make_source(resolved_path, logical_path, files, label, channel, method)
files = string(files(:));
source = empty_source();
source.path = string(resolved_path);
source.label = string(label);
source.logical_path = string(logical_path);
source.files = files;
source.is_split = numel(files) > 1;
source.part_count = numel(files);
source.resolution_method = string(method);
source.camera_index = double(channel.camera_index);
source.role = string(channel.role);
end

function manifest = read_cycle_manifest(cycle_path)
manifest = [];
mat_file = fullfile(cycle_path,'cycle_manifest.mat');
if isfile(mat_file)
    try
        saved = load(mat_file,'manifest');
        if isfield(saved,'manifest'), manifest = saved.manifest; return; end
    catch
    end
end
json_file = fullfile(cycle_path,'cycle_manifest.json');
if isfile(json_file)
    try, manifest = jsondecode(fileread(json_file)); catch, manifest = []; end
end
end

function [path_value,label] = manifest_channel(manifest,camera_index)
path_value = ""; label = "";
if isempty(manifest) || ~isstruct(manifest), return; end
if isfield(manifest,'actual') && isfield(manifest.actual,'movie_paths') ...
        && numel(manifest.actual.movie_paths) >= camera_index
    path_value = indexed_text(manifest.actual.movie_paths,camera_index);
end
if isfield(manifest,'spec') && isfield(manifest.spec,'labels') ...
        && numel(manifest.spec.labels) >= camera_index
    label = indexed_text(manifest.spec.labels,camera_index);
end
end

function value = indexed_text(values,idx)
if iscell(values), value = string(values{idx}); else, value = string(values(idx)); end
end

function tf = is_supported_movie_file(path_value)
[~,~,ext] = fileparts(path_value);
ext = lower(string(ext));
tf = ismember(ext,[".tif",".tiff",".bin"]);
if ext == ".mat"
    try
        variables = whos('-file',char(path_value));
        tf = any(strcmp({variables.name},'movie'));
    catch
        tf = false;
    end
end
end

function paths = sorted_paths(listing)
paths = string(fullfile({listing.folder},{listing.name}))';
[~,order] = sort(lower(paths));
paths = paths(order);
end

function value = last_part(path_value)
[~,name,ext] = fileparts(char(path_value)); value = string([name,ext]);
end

function source = empty_source()
source = struct('path',"",'label',"",'logical_path',"", ...
    'files',strings(0,1),'is_split',false,'part_count',0, ...
    'resolution_method',"",'camera_index',NaN,'role',"");
end
