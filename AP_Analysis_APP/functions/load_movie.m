function [movie_3d, movie_info] = load_movie(source, profile)
%LOAD_MOVIE Load one movie using one bound channel profile.
%   This algorithm is mode-independent. PROFILE is the sole authority for
%   role, frame rate, and frame transpose. SOURCE is either a path or a
%   scalar struct containing a path field and an optional label.
% Inputs:
%   source  Path text, or scalar struct with path and optional label,
%           logical_path, files, and resolution_method fields.
%   profile Scalar bound profile requiring role, positive frame_rate, and
%           logical transpose_before_analysis. It contains no request.
% Outputs:
%   movie_3d Numeric [ncols x nrows x frames] movie after configured
%            transpose. Split parts are concatenated along frames.
%   movie_info Scalar source/frame inventory; frame_rate is in Hz.

[source_path, source_label, source_meta] = normalize_source(source);
profile = validate_profile(profile);

[movie_loaded, ncols_raw, nrows_raw, nframes, actual_files] = ...
    read_movie_source(source_path, source_meta.files);
if ismatrix(movie_loaded)
    movie_3d = reshape(movie_loaded, ncols_raw, nrows_raw, []);
else
    movie_3d = movie_loaded;
end
if profile.transpose_before_analysis
    movie_3d = pagetranspose(movie_3d);
end

[ncols, nrows, nframes_after] = size(movie_3d);
if nframes_after ~= nframes
    error('AAA:Algorithms:LoadMovieFrameCountChanged', ...
        'Loading/transposing changed the frame count from %d to %d.', ...
        nframes, nframes_after);
end

movie_info = struct( ...
    'role', string(profile.role), ...
    'label', source_label, ...
    'source_path', string(source_path), ...
    'logical_source_path', source_meta.logical_path, ...
    'source_files', string(actual_files(:)), ...
    'source_part_count', numel(actual_files), ...
    'source_is_split', numel(actual_files) > 1, ...
    'source_resolution_method', source_meta.resolution_method, ...
    'frame_rate', double(profile.frame_rate), ...
    'frame_count', nframes_after, ...
    'original_frame_size', [ncols_raw, nrows_raw], ...
    'analysis_frame_size', [ncols, nrows], ...
    'transpose_before_analysis', logical(profile.transpose_before_analysis), ...
    'motion', empty_motion_info(), ...
    'created_at', datetime('now'), ...
    'updated_at', datetime('now'));
end

function [source_path, source_label, meta] = normalize_source(source)
source_label = "";
meta = struct('logical_path',"",'files',strings(0,1), ...
    'resolution_method',"direct_source");
if isstruct(source) && isscalar(source)
    if ~isfield(source, 'path') || strlength(string(source.path)) == 0
        error('AAA:Algorithms:MissingMoviePath', ...
            'A struct movie source must contain a nonempty path field.');
    end
    source_path = char(string(source.path));
    if isfield(source, 'label') && strlength(string(source.label)) > 0
        source_label = string(source.label);
    end
    if isfield(source,'logical_path')
        meta.logical_path = string(source.logical_path);
    end
    if isfield(source,'files')
        meta.files = string(source.files(:));
    end
    if isfield(source,'resolution_method') ...
            && strlength(string(source.resolution_method)) > 0
        meta.resolution_method = string(source.resolution_method);
    end
elseif (ischar(source) && (isrow(source) || isempty(source))) ...
        || (isstring(source) && isscalar(source))
    source_path = char(string(source));
else
    error('AAA:Algorithms:InvalidMovieSource', ...
        'Movie source must be a path or a scalar struct with a path field.');
end

if isempty(source_path) || ~(isfile(source_path) || isfolder(source_path))
    error('AAA:Algorithms:MovieSourceNotFound', ...
        'Movie source does not exist: %s', source_path);
end
if strlength(source_label) == 0
    [~, name, ext] = fileparts(source_path);
    source_label = string([name, ext]);
end
if strlength(meta.logical_path) == 0
    meta.logical_path = string(source_path);
end
if ~isempty(meta.files) && any(~isfile(meta.files))
    error('AAA:Algorithms:MoviePartNotFound', ...
        'One or more resolved movie parts no longer exist.');
end
end

function profile = validate_profile(profile)
if ~isstruct(profile) || ~isscalar(profile)
    error('AAA:Algorithms:InvalidChannelProfile', ...
        'Channel profile must be a scalar struct.');
end
required = {'role','frame_rate','transpose_before_analysis'};
for idx = 1:numel(required)
    if ~isfield(profile, required{idx}) || isempty(profile.(required{idx}))
        error('AAA:Algorithms:IncompleteChannelProfile', ...
            'Channel profile is missing %s.', required{idx});
    end
end
if ~isscalar(string(profile.role)) || strlength(string(profile.role)) == 0
    error('AAA:Algorithms:InvalidChannelRole', ...
        'profile.role must be one nonempty string.');
end
if ~isnumeric(profile.frame_rate) || ~isscalar(profile.frame_rate) ...
        || ~isfinite(profile.frame_rate) || profile.frame_rate <= 0
    error('AAA:Algorithms:InvalidFrameRate', ...
        'profile.frame_rate must be one positive finite number.');
end
if ~(islogical(profile.transpose_before_analysis) ...
        && isscalar(profile.transpose_before_analysis))
    error('AAA:Algorithms:InvalidTransposeFlag', ...
        'profile.transpose_before_analysis must be one logical value.');
end
end

function [movie, ncols, nrows, nframes, actual_files] = ...
        read_movie_source(source_path, explicit_files)
if nargin < 2, explicit_files = strings(0,1); end
if ~isempty(explicit_files)
    [movie,ncols,nrows,nframes,actual_files] = ...
        read_explicit_files(explicit_files);
    return;
end
if isfolder(source_path)
    [movie, ncols, nrows, nframes, actual_files] = read_movie_folder(source_path);
    return;
end

[~, ~, extension] = fileparts(source_path);
switch lower(extension)
    case {'.tif','.tiff'}
        movie = read_tiff_stack(source_path);
        [ncols, nrows, nframes] = size(movie);
    case '.mat'
        [movie, ncols, nrows, nframes] = read_mat_movie(source_path);
    case '.bin'
        [movie, ncols, nrows, nframes] = read_binary_movie(source_path);
    otherwise
        error('AAA:Algorithms:UnsupportedMovieFormat', ...
            'Unsupported movie file extension: %s', extension);
end
actual_files = string(source_path);
end

function [movie, ncols, nrows, nframes, actual_files] = read_movie_folder(folder_path)
tiff_files = [dir(fullfile(folder_path, '*.tif')); ...
    dir(fullfile(folder_path, '*.tiff'))];
tiff_files = tiff_files(~[tiff_files.isdir]);
if ~isempty(tiff_files)
    tiff_files = sort_tiff_files(tiff_files);
    stacks = cell(numel(tiff_files), 1);
    frame_size = [];
    for idx = 1:numel(tiff_files)
        file_path = fullfile(tiff_files(idx).folder, tiff_files(idx).name);
        stacks{idx} = read_tiff_stack(file_path);
        current_size = [size(stacks{idx}, 1), size(stacks{idx}, 2)];
        if isempty(frame_size)
            frame_size = current_size;
        elseif ~isequal(current_size, frame_size)
            error('AAA:Algorithms:MovieFrameSizeMismatch', ...
                'TIFF files in one movie folder have different frame sizes.');
        end
    end
    movie = cat(3, stacks{:});
    [ncols, nrows, nframes] = size(movie);
    actual_files = string(fullfile({tiff_files.folder},{tiff_files.name}))';
    return;
end

mat_files = dir(fullfile(folder_path, '*.mat'));
mat_files = mat_files(~[mat_files.isdir]);
if isempty(mat_files)
    error('AAA:Algorithms:EmptyMovieFolder', ...
        'Movie folder contains no TIFF or MAT movie: %s', folder_path);
end
[movie, ncols, nrows, nframes] = read_mat_movie( ...
    fullfile(mat_files(1).folder, mat_files(1).name));
actual_files = string(fullfile(mat_files(1).folder,mat_files(1).name));
end

function [movie,ncols,nrows,nframes,actual_files] = read_explicit_files(files)
files = string(files(:));
if numel(files) == 1
    [movie,ncols,nrows,nframes,actual_files] = read_movie_source(char(files),strings(0,1));
    return;
end
extensions = strings(size(files));
for idx=1:numel(files), extensions(idx)=file_extension(files(idx)); end
extensions = lower(extensions);
if ~all(ismember(extensions,[".tif",".tiff"]))
    error('AAA:Algorithms:UnsupportedMovieParts', ...
        'A multi-part movie must contain only TIFF parts.');
end
parts = cell(numel(files),1); frame_size = [];
for idx=1:numel(files)
    parts{idx}=read_tiff_stack(char(files(idx)));
    current=[size(parts{idx},1),size(parts{idx},2)];
    if isempty(frame_size), frame_size=current;
    elseif ~isequal(frame_size,current)
        error('AAA:Algorithms:MovieFrameSizeMismatch', ...
            'TIFF parts in one movie have different frame sizes.');
    end
end
movie=cat(3,parts{:}); [ncols,nrows,nframes]=size(movie); actual_files=files;
end

function value=file_extension(path_value)
[~,~,value]=fileparts(char(path_value));
end

function files = sort_tiff_files(files)
names = string({files.name});
suffix_numbers = nan(size(names));
for idx = 1:numel(names)
    [~, stem] = fileparts(names(idx));
    token = regexp(stem, '([0-9]+)$', 'tokens', 'once');
    if ~isempty(token)
        suffix_numbers(idx) = str2double(token{1});
    end
end
if all(isfinite(suffix_numbers))
    [~, order] = sortrows([suffix_numbers(:),(1:numel(files))']);
else
    [~, order] = sort(lower(names));
end
files = files(order);
end

function movie = read_tiff_stack(file_path)
info = imfinfo(file_path);
if isempty(info)
    error('AAA:Algorithms:EmptyTiffMovie', ...
        'TIFF movie contains no frames: %s', file_path);
end
first_frame = imread(file_path, 1, 'Info', info);
if ~ismatrix(first_frame)
    error('AAA:Algorithms:NonGrayscaleMovie', ...
        'Movie frames must be two-dimensional grayscale images: %s', file_path);
end
movie = zeros(size(first_frame, 1), size(first_frame, 2), ...
    numel(info), 'like', first_frame);
movie(:, :, 1) = first_frame;
for frame_idx = 2:numel(info)
    frame = imread(file_path, frame_idx, 'Info', info);
    if ~isequal(size(frame), size(first_frame))
        error('AAA:Algorithms:MovieFrameSizeMismatch', ...
            'TIFF stack contains inconsistent frame sizes: %s', file_path);
    end
    movie(:, :, frame_idx) = frame;
end
end

function [movie, ncols, nrows, nframes] = read_mat_movie(file_path)
loaded = load(file_path);
if ~isfield(loaded, 'movie') || ...
        (~isnumeric(loaded.movie) && ~islogical(loaded.movie)) ...
        || isempty(loaded.movie) || ndims(loaded.movie) > 3
    error('AAA:Algorithms:InvalidMatMovie', ...
        'MAT movie must contain one nonempty numeric variable named movie: %s', ...
        file_path);
end
movie = loaded.movie;
[ncols, nrows, nframes] = size(movie);
end

function [movie, ncols, nrows, nframes] = read_binary_movie(file_path)
[nrows, ncols] = read_binary_dimensions(file_path);
file_id = fopen(file_path, 'r', 'ieee-le');
if file_id < 0
    error('AAA:Algorithms:BinaryMovieOpenFailed', ...
        'Cannot open binary movie: %s', file_path);
end
cleanup = onCleanup(@() fclose(file_id));
values = fread(file_id, Inf, '*uint16');
frame_pixels = nrows * ncols;
nframes = numel(values) / frame_pixels;
if ~isfinite(nframes) || nframes < 1 || nframes ~= round(nframes)
    error('AAA:Algorithms:InvalidBinaryMovieSize', ...
        'Binary movie size is not divisible by nrows*ncols.');
end
movie = reshape(values, [nrows, ncols, nframes]);
movie = permute(movie, [2, 1, 3]);
end

function [nrows, ncols] = read_binary_dimensions(file_path)
folder_path = fileparts(file_path);
metadata_files = {fullfile(folder_path, 'movie_info.txt'), ...
    fullfile(folder_path, 'movie.txt')};
nrows = NaN;
ncols = NaN;
for idx = 1:numel(metadata_files)
    if ~isfile(metadata_files{idx})
        continue;
    end
    text_data = fileread(metadata_files{idx});
    row_token = regexp(text_data, ...
        'nrow[^\r\n]*?([0-9]+(?:\.[0-9]+)?)', 'tokens', 'once');
    col_token = regexp(text_data, ...
        'ncol[^\r\n]*?([0-9]+(?:\.[0-9]+)?)', 'tokens', 'once');
    if ~isempty(row_token)
        nrows = str2double(row_token{1});
    end
    if ~isempty(col_token)
        ncols = str2double(col_token{1});
    end
    if isfinite(nrows) && isfinite(ncols)
        break;
    end
end

if ~(isfinite(nrows) && isfinite(ncols))
    fallback_file = fullfile(folder_path, 'output_data.mat');
    if isfile(fallback_file)
        loaded = load(fallback_file);
        if isfield(loaded, 'Device_Data') && numel(loaded.Device_Data) >= 3 ...
                && isfield(loaded.Device_Data{3}, 'ROI')
            roi_info = loaded.Device_Data{3}.ROI;
            nrows = roi_info(2);
            ncols = roi_info(4);
        end
    end
end
if ~isfinite(nrows) || ~isfinite(ncols) ...
        || nrows < 1 || ncols < 1 ...
        || nrows ~= round(nrows) || ncols ~= round(ncols)
    error('AAA:Algorithms:MissingBinaryMovieDimensions', ...
        'Cannot resolve integer nrows/ncols for binary movie: %s', file_path);
end
end

function info = empty_motion_info()
info = struct( ...
    'applied', false, ...
    'method', '', ...
    'shift_file', '', ...
    'source_shift_file', '', ...
    'parameter_file', '');
end
