function written_files = write_tiff_stack(movie, save_path, ncols, nrows)
%WRITE_TIFF_STACK Write an integer movie using the established stack split.
%   Each output stack stays below approximately 3.6 GiB and uses the
%   historical NAME_stackNN.tif naming convention.

if ~isnumeric(movie) && ~islogical(movie)
    error('AAA:IO:InvalidTiffMovie', ...
        'TIFF movie must be numeric or logical.');
end
if isempty(movie) || ndims(movie) > 3
    error('AAA:IO:InvalidTiffMovie', ...
        'TIFF movie must be a nonempty 2-D or 3-D array.');
end
if nargin >= 3
    if nargin < 4 || ~is_positive_integer(ncols) || ~is_positive_integer(nrows)
        error('AAA:IO:InvalidTiffFrameSize', ...
            'ncols and nrows must be positive integers when supplied.');
    end
    movie = reshape(movie, nrows, ncols, []);
end

[nrows, ncols, nframes] = size(movie);
[bit_depth, bytes_per_pixel] = integer_tiff_type(movie);
save_path = char(string(save_path));
if isempty(save_path)
    error('AAA:IO:MissingTiffPath', 'TIFF save path cannot be empty.');
end
[folder_path, file_name] = fileparts(save_path);
if isempty(folder_path)
    folder_path = pwd;
elseif ~isfolder(folder_path)
    mkdir(folder_path);
end

max_bytes_per_stack = 3.6 * 1024^3;
bytes_per_frame = double(nrows) * double(ncols) * bytes_per_pixel;
max_frames = max(1, floor(max_bytes_per_stack / bytes_per_frame));
stack_count = ceil(nframes / max_frames);
written_files = strings(stack_count, 1);

for stack_idx = 1:stack_count
    first_frame = (stack_idx - 1) * max_frames + 1;
    last_frame = min(stack_idx * max_frames, nframes);
    stack_path = fullfile(folder_path, ...
        sprintf('%s_stack%02d.tif', file_name, stack_idx));
    write_one_stack(movie(:, :, first_frame:last_frame), ...
        stack_path, nrows, ncols, bit_depth);
    written_files(stack_idx) = string(stack_path);
end
end

function write_one_stack(movie, stack_path, nrows, ncols, bit_depth)
writer = Tiff(stack_path, 'w');
cleanup = onCleanup(@() writer.close());
for frame_idx = 1:size(movie, 3)
    writer.setTag('ImageLength', nrows);
    writer.setTag('ImageWidth', ncols);
    writer.setTag('Photometric', Tiff.Photometric.MinIsBlack);
    writer.setTag('BitsPerSample', bit_depth);
    writer.setTag('SamplesPerPixel', 1);
    writer.setTag('RowsPerStrip', 32);
    writer.setTag('PlanarConfiguration', Tiff.PlanarConfiguration.Chunky);
    writer.setTag('Compression', Tiff.Compression.None);
    writer.write(movie(:, :, frame_idx));
    if frame_idx < size(movie, 3)
        writer.writeDirectory();
    end
end
end

function [bit_depth, bytes_per_pixel] = integer_tiff_type(movie)
switch class(movie)
    case {'logical','uint8'}
        bit_depth = 8;
        bytes_per_pixel = 1;
    case 'uint16'
        bit_depth = 16;
        bytes_per_pixel = 2;
    case 'uint32'
        bit_depth = 32;
        bytes_per_pixel = 4;
    otherwise
        error('AAA:IO:UnsupportedTiffClass', ...
            'TIFF writer supports logical, uint8, uint16, or uint32 movies.');
end
end

function tf = is_positive_integer(value)
tf = isnumeric(value) && isscalar(value) && isfinite(value) ...
    && value >= 1 && value == round(value);
end
