function [intensity_time_series,ncols,nrows,nframes] = ...
                    load_movie(file_path, file_extension, batch_size)

% ----------Write by Liu-Yang Luorong and ChatGPT----------
% ----------POWERED by Zoulab in Peking University----------
% Date: 23.11.16
% MATLAB Version: R2022b
% 
% LOAD_MOVIE Loads movie data from various file formats and computes the intensity time series.
%
%   This function is designed to load movie data from files in TIFF, binary, or MAT format. It can
%   handle data from a single file or a directory containing multiple files. The function reads the 
%   movie data, computes the intensity time series, and returns this data along with the dimensions 
%   and number of frames.
%
%   Syntax:
%   [intensity_time_series, ncols, nrows, nframes] = load_movie(file_path, file_extension)
%   [intensity_time_series, ncols, nrows, nframes] = load_movie(file_path, file_extension, batch_size)
%
%   Parameters:
%   file_path - String specifying the path to the movie file or directory containing movie files.
%   file_extension - String specifying the file format ('tif', 'tiff', 'bin', 'mat').
%   batch_size - (Optional) Integer specifying the number of files to process in one batch. Default is 1000.
%
%   Returns:
%   intensity_time_series - A 2D matrix where each column represents the intensity values of a frame.
%   ncols - Number of columns in each frame.
%   nrows - Number of rows in each frame.
%   nframes - Total number of frames in the movie.
%
%   Description:
%   The function processes TIFF files either as a batch from a directory or individually from a file. 
%   For binary files, it reads associated dimensions from a text file. For MAT files, it directly loads 
%   the stored variables. The function also includes basic error handling and progress updates.
%
%   Examples:
%   [its, nc, nr, nf] = load_movie('path/to/movie.tif', 'tif');
%   [its, nc, nr, nf] = load_movie('path/to/movies/', 'tif', 500);
%
%   Notes:
%   - The function is designed to manage memory efficiently by processing files in batches.
%   - Ensure that the file format assumptions match the actual structure of your data files.
%   - The function provides verbose output for progress tracking, especially useful for large datasets.
%
%   See also IMREAD, IMFINFO, FOPEN, FREAD, RESHAPE.


if nargin < 3
    batch_size = 1000;
end

nrows = NaN;
ncols = NaN;
nframes = NaN;

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
    [nrows, ncols] = size(im);

    % Initialize parameter containing intensity time series
    %         if length(file_names) > 12000 % Memory save
    %             intensity_time_series = zeros(num_rows*num_cols, length(12000), 'double');
    %         else
    %
    %         end
    intensity_time_series = zeros(nrows*ncols, length(file_names), 'double');
    % Loop through all TIF files and populate intensity time series parameter
    num_files = numel(file_list);
    nframes = num_files;
    batch_start = 1;
    while batch_start <= num_files
        %fprintf('Processing batch starting at file %d\n', batch_start);

        % Load batch of TIF files
        batch_end = min(batch_start + batch_size - 1, num_files);
        batch_files = file_names(batch_start:batch_end);
        num_batch_files = length(batch_files);
        batch_images = zeros(nrows, ncols, num_batch_files, 'uint16');
        for i = 1:num_batch_files
            fprintf('Processing %d/%d\n', i+batch_start-1, num_files);
            batch_images(:,:,i) = imread(fullfile(file_path, batch_files{i}));
        end

        % Reshape batch images into 2D array
        batch_images_2D = reshape(batch_images, nrows*ncols, []);

        % Append batch intensity time series to full intensity time series
        intensity_time_series(:,batch_start:batch_end) = double(batch_images_2D);

        % Update batch start and end indices
        batch_start = batch_end + 1;
        batch_end = min(batch_start + batch_size - 1, num_files);
    end


elseif isequal(file_extension,'tif') || isequal(file_extension,'tiff')
    % Read TIFF movie
    info = imfinfo(file_path);
    [ncols, nrows] = size(imread(file_path, 1));
    nframes = numel(info);
    % Initialize parameter containing intensity time series
    intensity_time_series = zeros(nrows*ncols, nframes, 'double');
    % Loop through all frames and populate intensity time series parameter
    for i = 1:nframes
        im = imread(file_path, i);
        intensity_time_series(:, i) = double(im(:));
        % Print progress
        fprintf('Processing frame %d for %d\n', i,nframes)
    end

elseif isequal(string(file_extension),'bin')
    filename = [fileparts(file_path) '\movie_info.txt'];  % 指定文本文件的名称
    fileID = fopen(filename, 'r');  % 以读取模式打开文件
    if fileID == -1
        filename = [fileparts(file_path) '\movie.txt'];  % change文本文件的名称
        fileID = fopen(filename, 'r');
        if fileID == -1
        error('Failed to open file. change the txt name manually');
        end
    end

    % read ncol and nrow
    while ~feof(fileID)  % 继续读取，直到到达文件末尾
        line = fgetl(fileID);  % 读取文件的下一行
        if contains(line, 'nrow =')  % 检查该行是否包含 'x ='
            % 使用 str2double 和 strtrim 函数从字符串中提取数值
            nrows = str2double(strtrim(extractAfter(line, 'nrow =')));
        elseif contains(line, 'ncol =')  % 检查该行是否包含 'y ='
            % 使用 str2double 和 strtrim 函数从字符串中提取数值
            ncols = str2double(strtrim(extractAfter(line, 'ncol =')));
        end
        if ~isnan(nrows) && ~isnan(ncols)
            break;  % 找到数值后退出循环
        end
    end

    %open file readBinMov
    % read file into tmp vector
    Movid = fopen(file_path);                  % open file
    Mov = fread(Movid, '*uint16', 'l');       % uint16, little endian
    fclose(Movid);                            % close file

    % reshape vector into appropriately oriented, 3D array
    nframes = length(Mov)/(nrows*ncols);
    intensity_time_series = reshape(Mov, [nrows, ncols, nframes]);
    intensity_time_series = permute(intensity_time_series, [2 1 3]);
    intensity_time_series = reshape(intensity_time_series, [nrows*ncols, nframes]);
    fclose(fileID);

elseif isequal(string(file_extension),'mat')
    
    fprintf('Loading saved data...\n')
    file = load(file_path);
    intensity_time_series = file.intensity_time_series;
    ncols = file.ncols;
    nrows = file.nrows;
    nframes = file.nframes;
end



end

% elseif isequal(string(file_extension),'cxd')
%     % Load cxd file using Bio-Formats toolbox
%     reader = bfGetReader(file_path);
%     nframes = reader.getSizeT();
%     nrows = javaMethod('getSizeY', reader);
%     ncols = javaMethod('getSizeX', reader);
% 
%     % Initialize parameter containing intensity time series
%     intensity_time_series = zeros(nrows*ncols, nframes, 'double');
%     % Loop through all frames and populate intensity time series parameter
%     for i = 1:nframes
%         im = bfGetPlane(reader, i);
%         intensity_time_series(:, i) = double(im(:));
%         % Print progress
%         fprintf('Processing frame %d for %d\n', i,nframes);
%     end