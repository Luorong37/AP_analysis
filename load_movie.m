% Write by Liu-Yang Luorong and ChatGPT3.5
% 202303


function [intensity_time_series, num_cols,num_rows,num_frames] = load_movie(file_path,file_extension, batch_size)
%LOAD_TIF_MOVIE Load a TIFF movie or a directory containing TIFF files and
%return a 2D array containing intensity time series for each pixel.
%
%   intensity_time_series = LOAD_TIF_MOVIE(folder_path, file_extension)
%   returns a 2D array containing intensity time series for each pixel in a
%   TIFF movie. If the input folder_path is a path to a TIFF movie, the
%   function loads the movie and extracts intensity time series for each
%   pixel. If the input folder_path is a directory containing TIFF files,
%   the function loads all TIFF files in the directory, sorts them
%   numerically, and extracts intensity time series for each pixel from all
%   files.
%
%   folder_path: A string containing a path to a TIFF movie or a directory
%   containing TIFF files.
%
%   file_extension: A string containing a file extension of TIFF files to
%   load from the directory. This argument is ignored if folder_path is a
%   path to a TIFF movie. Default value is '*.tif'.
%
%   intensity_time_series: A 2D array containing intensity time series for
%   each pixel in the TIFF movie or directory. The rows of the array
%   represent pixels and the columns represent time points.

%   num_rows: An integer indicating the number of rows in the TIFF movie or
%   directory.
%
%   num_cols: An integer indicating the number of columns in the TIFF movie
%   or directory.

%   batch_size: An integer defining the number of TIFF files to load at a
%   time. Default value is 1000.

%   For cxd file, please download bio-formats in https://www.openmicroscopy.org/bio-formats/downloads/

%   Example:
%   intensity_time_series = load_tif_movie('D:\20230321-150800', '*.tif');

if nargin < 3
    batch_size = 1000;
end

num_rows = NaN;
num_cols = NaN;
num_frames = NaN;

% Check if folder_path is a folder or path to a tif movie
if isfolder(file_path)
    % Get all TIFF file names in folder
    file_list = dir(fullfile(file_path, ['*.' file_extension]));

    % Sort file names
    file_names = {file_list.name};
    file_nums = cellfun(@(x) sscanf(x,'%d.tif'), file_names);
    [~, idx] = sort(file_nums);
    file_names = file_names(idx);

    % Load first image to get dimensions
    im = imread(fullfile(file_path, file_names{1}));
    [num_rows, num_cols] = size(im);

    % Initialize parameter containing intensity time series
    %         if length(file_names) > 12000 % Memory save
    %             intensity_time_series = zeros(num_rows*num_cols, length(12000), 'double');
    %         else
    %
    %         end
    intensity_time_series = zeros(num_rows*num_cols, length(file_names), 'double');
    % Loop through all TIF files and populate intensity time series parameter
    num_files = numel(file_list);
    num_frames = num_files;
    batch_start = 1;
    while batch_start <= num_files
        %fprintf('Processing batch starting at file %d\n', batch_start);

        % Load batch of TIF files
        batch_end = min(batch_start + batch_size - 1, num_files);
        batch_files = file_names(batch_start:batch_end);
        num_batch_files = length(batch_files);
        batch_images = zeros(num_rows, num_cols, num_batch_files, 'uint16');
        for i = 1:num_batch_files
            fprintf('Processing %d/%d\n', i+batch_start-1, num_files);
            batch_images(:,:,i) = imread(fullfile(file_path, batch_files{i}));
        end

        % Reshape batch images into 2D array
        batch_images_2D = reshape(batch_images, num_rows*num_cols, []);

        % Append batch intensity time series to full intensity time series
        intensity_time_series(:,batch_start:batch_end) = double(batch_images_2D);

        % Update batch start and end indices
        batch_start = batch_end + 1;
        batch_end = min(batch_start + batch_size - 1, num_files);
    end


elseif isequal(file_extension,'tif')
    % Read TIFF movie
    info = imfinfo(file_path);
    [num_cols, num_rows] = size(imread(file_path, 1));
    num_frames = numel(info);
    % Initialize parameter containing intensity time series
    intensity_time_series = zeros(num_rows*num_cols, num_frames, 'double');
    % Loop through all frames and populate intensity time series parameter
    for i = 1:num_frames
        im = imread(file_path, i);
        intensity_time_series(:, i) = double(im(:));
        % Print progress
        fprintf('Processing frame %d for %d\n', i,num_frames)
    end

elseif isequal(file_extension,'tiff')
    % Read TIFF movie
    info = imfinfo(file_path);
    [num_rows, num_cols] = size(imread(file_path, 1));
    num_frames = numel(info);
    % Initialize parameter containing intensity time series
    intensity_time_series = zeros(num_rows*num_cols, num_frames, 'double');
    % Loop through all frames and populate intensity time series parameter
    for i = 1:num_frames
        im = imread(file_path, i);
        intensity_time_series(:, i) = double(im(:));
        % Print progress
        fprintf('Processing frame %d for %d\n', i,num_frames)
    end

elseif isequal(string(file_extension),'dcimg');
    image_data = hamamatsuOpen(file_path);
    open = dcimg_ndarray(image_data,12000);
    num_frames = numel(image_data);
    [num_rows, num_cols] = size(image_data);
    intensity_time_series = zeros(num_rows*num_cols, num_frames, 'double');
    for i = 1:num_frames
        im = dcimgread(folder_path, i);
        intensity_time_series(:, i) = double(im(:));
        % Print progress
        fprintf('Processing frame %d for %d\n', i,num_frames)
    end

elseif isequal(string(file_extension),'cxd');
    % Load cxd file using Bio-Formats toolbox
    reader = bfGetReader(file_path);
    num_frames = reader.getSizeT();
    num_rows = javaMethod('getSizeY', reader);
    num_cols = javaMethod('getSizeX', reader);

    % Initialize parameter containing intensity time series
    intensity_time_series = zeros(num_rows*num_cols, num_frames, 'double');
    % Loop through all frames and populate intensity time series parameter
    for i = 1:num_frames
        im = bfGetPlane(reader, i);
        intensity_time_series(:, i) = double(im(:));
        % Print progress
        fprintf('Processing frame %d for %d\n', i,num_frames);
    end
    %%

elseif isequal(string(file_extension),'bin');
    filename = [fileparts(file_path) '\movie.txt'];  % 指定文本文件的名称
    fileID = fopen(filename, 'r');  % 以读取模式打开文件
    % read ncol and nrow
    if fileID == -1
        error('Failed to open file.');
    end
    while ~feof(fileID)  % 继续读取，直到到达文件末尾
        line = fgetl(fileID);  % 读取文件的下一行
        if contains(line, 'nrow =')  % 检查该行是否包含 'x ='
            % 使用 str2double 和 strtrim 函数从字符串中提取数值
            num_rows = str2double(strtrim(extractAfter(line, 'nrow =')));
        elseif contains(line, 'ncol =')  % 检查该行是否包含 'x ='
            % 使用 str2double 和 strtrim 函数从字符串中提取数值
            num_cols = str2double(strtrim(extractAfter(line, 'ncol =')));
        end
        if ~isnan(num_rows) && ~isnan(num_cols)
            break;  % 找到数值后退出循环
        end
    end

    %open file readBinMov
    % read file into tmp vector
    Movid = fopen(file_path);                  % open file
    Mov = fread(Movid, '*uint16', 'l');       % uint16, little endian
    fclose(Movid);                            % close file

    % reshape vector into appropriately oriented, 3D array
    num_frames = length(Mov)/(num_rows*num_cols);
    intensity_time_series = reshape(Mov, [num_rows, num_cols, num_frames]);
    intensity_time_series = permute(intensity_time_series, [2 1 3]);
    intensity_time_series = reshape(intensity_time_series, [num_rows*num_cols, num_frames]);
    %        intensity_time_series = reshape(intensity_time_series, [num_cols,num_rows, num_frames]);


    fclose(fileID);


end
end