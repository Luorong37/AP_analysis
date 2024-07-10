function [movie,ncols,nrows,nframes] = ...
    load_movie(file_path, batch_size, integrate)

% ----------Write by Liu-Yang Luorong and ChatGPT----------
% ----------POWERED by Zoulab in Peking University----------
% Date: 24.07.10
% MATLAB Version: R2023a
%
% LOAD_MOVIE Loads movie data from various file formats and computes the intensity time series.
%
% This function is designed to load movie data from files in TIFF, binary, or MAT format. It can
% handle data from a single file or a directory containing multiple files. The function reads the
% movie data, computes the intensity time series, and returns this data along with the dimensions
% and number of frames.
%
% Syntax:
%   [movie, ncols, nrows, nframes] = load_movie(file_path)
%   [movie, ncols, nrows, nframes] = load_movie(file_path, batch_size)
%
% Parameters:
%   file_path - String specifying the path to the movie file or directory containing movie files.
%   batch_size - (Optional) Integer specifying the number of files to process in one batch. Default is 1000.
%
% Returns:
%   movie - A 2D matrix where each column represents the intensity values of a frame.
%   ncols - Number of columns in each frame.
%   nrows - Number of rows in each frame.
%   nframes - Total number of frames in the movie.
%
% Description:
%   The function processes TIFF files either as a batch from a directory or individually from a file.
%   For binary files, it reads associated dimensions from a text file. For MAT files, it directly loads
%   the stored variables. The function also includes basic error handling and progress updates.
%
% Examples:
%   [its, nc, nr, nf] = load_movie('path/to/movie.tif');
%   [its, nc, nr, nf] = load_movie('path/to/movies/', 500);
%
% Notes:
%   - The function is designed to manage memory efficiently by processing files in batches.
%   - Ensure that the file format assumptions match the actual structure of your data files.
%   - The function provides verbose output for progress tracking, especially useful for large datasets.
%
% See also IMREAD, IMFINFO, FOPEN, FREAD, RESHAPE.


if nargin < 2
    batch_size = 1000;
end

if nargin < 3
    integrate = false;
end

nrows = NaN;
ncols = NaN;
nframes = NaN;
max_file_size = 4 * 1024^3

% Check if folder_path is a folder or path to a tif movie
if isfolder(file_path)
    % Get all TIFF file names in folder
    file_list = dir(fullfile(file_path, '*.tif'));
    file_names = {file_list.name};
    

    % Extract numbers from file names using a regular expression
    % Handle both cases: filenames with only numbers and filenames with _number
    file_nums = cellfun(@(x) extractFileNumber(x), file_names);
    [~, idx] = sort(file_nums);
    file_names = file_names(idx);

    % Load first image to get dimensions
    im = imread(fullfile(file_path, file_names{1}));
    [nrows, ncols] = size(im);
    num_files = numel(file_list);
    movie = zeros(nrows*ncols, num_files, 'double');
    nframes = num_files;

    % Loop through all TIF files and populate intensity time series parameter
    batch_start = 1;
    prev_percentage = -1; % Initialize with -1 so the first update is always printed

    tic;
    while batch_start <= num_files
        % Load batch of TIF files
        batch_end = min(batch_start + batch_size - 1, num_files);
        batch_range = batch_start:batch_end;

        for i = batch_range
            % Read the current image, store the image directly in 'movie'
            current_image = imread(fullfile(file_path, file_names{i}));
            movie(:, i) = double(reshape(current_image, nrows*ncols, 1));

            % Calculate and display progress if percentage changes
            current_percentage = floor(((i - 1) / num_files) * 100);
            if current_percentage > prev_percentage
                elapsed = toc;
                remaining = elapsed / ((i-1) / num_files) - elapsed;
                fprintf('Processing %d/%d files (%d%% complete). Estimated time remaining: %.2f seconds\n', ...
                    i - 1, num_files, current_percentage, remaining);
                prev_percentage = current_percentage;
            end

        end

        % Update batch start and end indices
        batch_start = batch_end + 1;
    end

else
    [~, ~, file_extension] = fileparts(file_path);
    if isequal(file_extension,'.tif') || isequal(file_extension,'.tiff')
        % Read TIFF movie
        info = imfinfo(file_path);
        [ncols, nrows] = size(imread(file_path, 1));
        nframes = numel(info);
        % Initialize parameter containing intensity time series
        movie = zeros(nrows*ncols, nframes, 'double');
        % Loop through all frames and populate intensity time series parameter
        for i = 1:nframes
            im = imread(file_path, i);
            movie(:, i) = double(im(:));
            % Print progress
            fprintf('Processing frame %d for %d\n', i,nframes)
        end

    elseif isequal(file_extension,'.bin')
        filename = [fileparts(file_path) '\movie_info.txt'];  % 指定文本文件的名称
        fileID = fopen(filename, 'r');  % 以读取模式打开文件
        if fileID == -1
            filename = [fileparts(file_path) '\movie.txt'];  % 更改文本文件的名称
            fileID = fopen(filename, 'r');
        end
        if fileID == -1
            %error('Failed to open file. change the txt name manually');
            movie_info = load(fullfile(fileparts(file_path),'output_data.mat'));
            ROI_info = movie_info.Device_Data{3}.ROI;
            nrows = ROI_info(2);
            ncols = ROI_info(4);
        end

        % read ncol and nrow
        % while ~feof(fileID)  % 继续读取，直到到达文件末尾
        %     line = fgetl(fileID);  % 读取文件的下一行
        %     if contains(line, 'nrow =')  % 检查该行是否包含 'x ='
        %         % 使用 str2double 和 strtrim 函数从字符串中提取数值
        %         nrows = str2double(strtrim(extractAfter(line, 'nrow =')));
        %     elseif contains(line, 'ncol =')  % 检查该行是否包含 'y ='
        %         % 使用 str2double 和 strtrim 函数从字符串中提取数值
        %         ncols = str2double(strtrim(extractAfter(line, 'ncol =')));
        %     end
        %     if ~isnan(nrows) && ~isnan(ncols)
        %         break;  % 找到数值后退出循环
        %     end
        % end

        %open file readBinMov
        % read file into tmp vector
        Movid = fopen(file_path);                  % open file
        Mov = fread(Movid, '*uint16', 'l');       % uint16, little endian
        fclose(Movid);                            % close file

        % reshape vector into appropriately oriented, 3D array
        nframes = length(Mov)/(nrows*ncols);
        movie = reshape(Mov, [nrows, ncols, nframes]);
        movie = permute(movie, [2 1 3]);
        movie = reshape(movie, [nrows*ncols, nframes]);
        % fclose(fileID);

    elseif isequal(string(file_extension),'.mat')

        fprintf('Loading saved data...\n')
        file = load(file_path);
        movie = file.movie;
        ncols = file.ncols;
        nrows = file.nrows;
        nframes = file.nframes;
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