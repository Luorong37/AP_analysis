function create_tiff_stack(file_path, save_path, batch_size, max_file_size)
% CREATE_TIFF_STACK Loads single TIFF images and saves them as a stack of TIFF images.
% The resulting stack is split into multiple files if the size exceeds max_file_size.
%
% Parameters:
% file_path - String specifying the path to the directory containing TIFF files.
% batch_size - (Optional) Integer specifying the number of files to process in one batch. 
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
    batch_size = 0;
end

if nargin < 4
    max_file_size = 4 * 1024^3; % 4GB
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

if isempty(idx) || length(cell2mat(file_nums)) < 400
    sortable = false;
end

if ~isempty(save_path)
    save_list = dir(fullfile(save_path, '*.tif'));
    save_names = {save_list.name};
    frames = 0;
    parfor i = 1: length(save_list)
        t = Tiff(fullfile(save_path,save_names{i}));
        tframes = numel(imfinfo(t.FileName));
        frames = frames + tframes;
        t.close;
    end
    if frames == length(file_nums)
        copyed = 1;
        fprintf('Stacked already copyed.\n')
    else
        copyed = 0;
    end
end
        

if sortable && ~copyed
    file_names = file_names(idx);

    % Load first image to get dimensions
    im = imread(fullfile(file_path, file_names{1}));
    [nrows, ncols] = size(im);
    num_files = numel(file_list);

    % Initialize variables for stack creation
    % batch_start = 1;
    % file_counter = 1;
    % stack_counter = 1;
    % prev_percentage = -1;

    % Get directory path and base name for the stack files
    stack_base = [file_base '_stack'];

    % tiff_file_name = fullfile(save_path, sprintf('%s%02d.tif', stack_base, file_counter));
    % t = Tiff(tiff_file_name, 'w'); % 'w' might can not save stack that more than 4GB.
    % current_file_size = 0;

    % Loop through all TIF files and create TIFF stack
    % tic;
    % print_count = 0;
    
    % 判断文件数量和大小（注：16位tiff一个像素点大小为2字节）
    single_tiff_size = nrows * ncols * 2;       % 计算每个单帧tiff文件的大小
    im_info = whos('im');
    im_size = im_info.bytes;
    if im_size ~= single_tiff_size              % 和读的第一帧实际大小check以下
        warning('Maybe not 16-bit tiff files');
    end

    if (nrows == ncols) == 512
        if num_files <= 32000
            stack_tiff_count_preference = 4000;
        end
        if 32000 < num_files && num_files <= 48000
            stack_tiff_count_preference = 6000;
        end
        if 48000 < num_files && num_files <= 64000
            stack_tiff_count_preference = 8000;
        end
        if 64000 < num_files && num_files <= 96000
            stack_tiff_count_preference = 6000;
        end
        if num_files > 96000
            stack_tiff_count_preference = 8000;
        end
    end
    
    if batch_size
        stack_tiff_count_preference = batch_size;
    end
    stack_tiff_max_count = min(stack_tiff_count_preference, floor(max_file_size / single_tiff_size));        % 选择生成的每个stack中tiff的最大页数（帧数为4000-8000最大化并行收益，且不超过4GB）（512*512的tiff，8000帧约3.9GB）
    
    parallel_create_tiff_stack(file_path, save_path, stack_base, num_files, stack_tiff_max_count, file_names);

%     while batch_start <= num_files
%         % Load batch of TIF files
%         batch_end = min(batch_start + batch_size - 1, num_files);
%         batch_range = batch_start:batch_end;
% 
%         for i = batch_range
%             try
%                 % Read the current image
%                 current_image = uint16(imread(fullfile(file_path, file_names{i})));% Convert to appropriate type
% 
%                 % Get the size of the current image
%                 im_info = whos('current_image');
%                 im_size = im_info.bytes;
% 
%                 % Check if adding this image exceeds max_file_size
%                 if current_file_size + im_size > max_file_size
%                     t.close();
%                     file_counter = file_counter + 1;
%                     tiff_file_name = fullfile(save_path, sprintf('%s%02d.tif', stack_base, file_counter));
%                     t = Tiff(tiff_file_name, 'w');
%                     current_file_size = 0;
%                     stack_counter = 1;
%                     fprintf('Creating new stack file: %s\n', tiff_file_name);
%                     print_count = 0;
%                 end
% 
%                 % Write the current image to the TIFF stack
%                 if stack_counter > 1
%                     t.writeDirectory();
%                 end
% 
%                 t.setTag('ImageLength', nrows);
%                 t.setTag('ImageWidth', ncols);
%                 t.setTag('Photometric', Tiff.Photometric.MinIsBlack);
%                 t.setTag('BitsPerSample', 16);
%                 t.setTag('SamplesPerPixel', 1);
%                 t.setTag('RowsPerStrip', 16);
%                 t.setTag('PlanarConfiguration', Tiff.PlanarConfiguration.Chunky);
%                 t.setTag('Compression', Tiff.Compression.None);
%                 t.setTag('Software', 'MATLAB');
%                 t.write(current_image);
% 
%                 % Update file size and stack counter
%                 current_file_size = current_file_size + im_size;
%                 stack_counter = stack_counter + 1;
% 
%                 % Calculate and display progress if percentage changes
%                 current_percentage = floor(((i - 1) / num_files) * 100);
%                 if current_percentage > prev_percentage
%                     elapsed = toc;
%                     remaining = elapsed / ((i-1) / num_files) - elapsed;
%                     fprintf(repmat('\b',1,print_count))
%                     print_count = fprintf('Processing %d/%d files (%d%% complete). Estimated time remaining: %.2f seconds\n', ...
%                         i - 1, num_files, current_percentage, remaining);
%                     prev_percentage = current_percentage;
%                 end
%             catch ME
%                 t.close()
%                 fprintf('Error processing file %s: %s\n', file_names{i}, ME.message);
%                 continue;
%             end
%         end
% 
%         % Update batch start and end indices
%         batch_start = batch_end + 1;
%     end
%     t.close();
% 
% fprintf('Stacked tif created.\n')

elseif ~sortable
    fprintf('tifs are not sortable. files will directly copyed.\n')
    for i = 1:length(file_names)
        source_file = fullfile(file_path, file_names{i});
        destination_file = fullfile(save_path, file_names{i});
        copyfile(source_file, destination_file);
    end
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

function parallel_create_tiff_stack(file_path, save_path, stack_base, num_files, stack_tiff_max_count, file_names)
% 用于在create_tiff_stack中分组并行读取多个tif文件并合并为tiff_stack
%

t0 = tic;

im = imread(fullfile(file_path, file_names{1}));
[nrows, ncols] = size(im);

stack_num = ceil(num_files / stack_tiff_max_count);         % 计算所需的stack数目
if stack_num > 0
    gcp;
end
stack_nframes(1:stack_num) = zeros(1, stack_num);           % 用stack_nframes(i)表示第i个stack包含的tiff页数
stack_nframes(1:(stack_num-1)) = stack_tiff_max_count;
stack_nframes(stack_num) = num_files - (stack_tiff_max_count * (stack_num-1));

% 进度跟踪：worker -> 客户端
q        = parallel.pool.DataQueue;
nDone    = 0;    % 累计已完成数
print_text = 0;
prevPcnt = 0;    % 上一次已打印的百分比
afterEach(q, @updateProgress);  % 嵌套函数，能直接用 nframes/t0/nDone/prevPcnt

parfor i = 1:stack_num
    
    tiff_file_name = fullfile(save_path, sprintf('%s%02d.tif', stack_base, i));
    fprintf('Creating new stack file: %s\n', tiff_file_name);
    t = Tiff(tiff_file_name, 'w');
    batch_start = 1 + (i - 1) * stack_tiff_max_count;
    batch_end = batch_start + stack_nframes(i) - 1;
    for j = batch_start:batch_end
        current_image = uint16(imread(fullfile(file_path, file_names{j})));         % 读取当前图片

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

        if j < batch_end
            t.writeDirectory();
        end

        send(q, 1);  % 每完成一个发一个通知（内容可忽略）

    end
    t.close();
end

% -------- 嵌套回调：只在百分比增加时打印 --------
function updateProgress(~)
    nDone  = nDone + 1;
    pcnt   = floor(nDone / num_files * 100);
    if pcnt > prevPcnt
        elapsed   = toc(t0);
        remaining = (elapsed / nDone) * (num_files - nDone);
        fprintf(repmat('\b',1,print_text));    % 擦除上次的进度信息
        print_text = fprintf('Processing %d/%d (%d%% complete). Estimated time remaining: %.1f s\n', ...    % 输出进度信息
                nDone, num_files, pcnt, remaining);
        prevPcnt = pcnt;
    end
end

t1 = toc(t0);
fprintf('Finished processing after %d s, ',round(t1));
fprintf('Stacked tif created in %s.\n', save_path);

end