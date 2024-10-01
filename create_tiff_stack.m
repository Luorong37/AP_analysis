function create_tiff_stack(file_path, save_path, batch_size, max_file_size)
% CREATE_TIFF_STACK Loads single TIFF images and saves them as a stack of TIFF images.
% The resulting stack is split into multiple files if the size exceeds max_file_size.
%
% Parameters:
% file_path - String specifying the path to the directory containing TIFF files.
% batch_size - (Optional) Integer specifying the number of files to process in one batch. Default is 1000.
% max_file_size - (Optional) Maximum size of each TIFF file in bytes. Default is 4GB (4*1024^3 bytes).
%
% Example:
% create_tiff_stack('path/to/tiff/files/', 500, 4*1024^3);

[file_dir, file_base, ~] = fileparts(file_path);
stack_dir = file_dir;

if nargin < 2
    save_path = stack_dir;
end

if nargin < 3
    batch_size = 4000;
end
if nargin < 4
    max_file_size = 3.6 * 1024^3; % 4GB
end

% Get all TIFF file names in folder
file_list = dir(fullfile(file_path, '*.tif'));
file_names = {file_list.name};

% Extract numbers from file names using a regular expression
% Handle both cases: filenames with only numbers and filenames with _number
file_nums = cellfun(@(x) extractFileNumber(x), file_names,'UniformOutput', false);

try
    [~, idx] = sort(cell2mat(file_nums));
    sortable = true;
catch
    sortable = false;
end

if sortable
    file_names = file_names(idx);

    % Load first image to get dimensions
    im = imread(fullfile(file_path, file_names{1}));
    [nrows, ncols] = size(im);
    num_files = numel(file_list);

    % Initialize variables for stack creation
    batch_start = 1;
    file_counter = 1;
    stack_counter = 1;
    prev_percentage = -1;

    % Get directory path and base name for the stack files
    stack_base = [file_base '_stack'];

    tiff_file_name = fullfile(save_path, sprintf('%s%02d.tif', stack_base, file_counter));
    t = Tiff(tiff_file_name, 'w8'); % 'w' might can not save stack that more than 4GB.
    current_file_size = 0;

    % Loop through all TIF files and create TIFF stack
    tic;
    print_count = 0 ;
    while batch_start <= num_files
        % Load batch of TIF files
        batch_end = min(batch_start + batch_size - 1, num_files);
        batch_range = batch_start:batch_end;

        for i = batch_range
            try
                % Read the current image
                current_image = uint16(imread(fullfile(file_path, file_names{i})));% Convert to appropriate type

                % Get the size of the current image
                im_info = whos('current_image');
                im_size = im_info.bytes;

                % Check if adding this image exceeds max_file_size
                if current_file_size + im_size > max_file_size
                    t.close();
                    file_counter = file_counter + 1;
                    tiff_file_name = fullfile(save_path, sprintf('%s%02d.tif', stack_base, file_counter));
                    t = Tiff(tiff_file_name, 'w8');
                    current_file_size = 0;
                    stack_counter = 1;
                    fprintf('Creating new stack file: %s\n', tiff_file_name);
                    print_count = 0;
                end

                % Write the current image to the TIFF stack
                if stack_counter > 1
                    t.writeDirectory();
                end

                t.setTag('ImageLength', nrows);
                t.setTag('ImageWidth', ncols);
                t.setTag('Photometric', Tiff.Photometric.MinIsBlack);
                t.setTag('BitsPerSample', 16);
                t.setTag('SamplesPerPixel', 1);
                t.setTag('RowsPerStrip', 16);
                t.setTag('PlanarConfiguration', Tiff.PlanarConfiguration.Chunky);
                t.setTag('Compression', Tiff.Compression.None);
                t.setTag('Software', 'MATLAB');
                t.write(current_image);

                % Update file size and stack counter
                current_file_size = current_file_size + im_size;
                stack_counter = stack_counter + 1;

                % Calculate and display progress if percentage changes
                current_percentage = floor(((i - 1) / num_files) * 100);
                if current_percentage > prev_percentage
                    elapsed = toc;
                    remaining = elapsed / ((i-1) / num_files) - elapsed;
                    fprintf(repmat('\b',1,print_count))
                    print_count = fprintf('Processing %d/%d files (%d%% complete). Estimated time remaining: %.2f seconds\n', ...
                        i - 1, num_files, current_percentage, remaining);
                    prev_percentage = current_percentage;
                end
            catch ME
                t.close()
                fprintf('Error processing file %s: %s\n', file_names{i}, ME.message);
                continue;
            end
        end

        % Update batch start and end indices
        batch_start = batch_end + 1;
    end
    t.close();

fprintf('Stacked tif created.\n')

else
    fprintf('tifs are not sortable.\n')
end
end

function num = extractFileNumber(filename)
% Extract the numeric part of the filename
% Handle both cases: filenames with only numbers and filenames with _number
num = [];
[~, name, ~] = fileparts(filename);
if all(isstrprop(name, 'digit'))
    num = str2double(name);
else
    tokens = regexp(name, '_([0-9]+)$', 'tokens');
    if ~isempty(tokens)
        num = str2double(tokens{1}{1});
    end
end
end