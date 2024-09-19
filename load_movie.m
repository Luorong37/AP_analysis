function [movie,ncols,nrows,nframes] = load_movie(file_path)

% ----------Write by Liu-Yang Luorong and ChatGPT----------
% ----------POWERED by Zoulab in Peking University----------
% Date: 24.07.14
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

nrows = NaN;
ncols = NaN;
nframes = NaN;

% Check if folder_path is a folder or path to a tif movie
if isfolder(file_path)
    % Get all TIFF file names in folder
    file_list = dir(fullfile(file_path, '*.tif'));
    file_names = {file_list.name};

    % Sort file names
    % other wise, will be like [4, 40, 400, 4000, 4001]
    file_nums = cellfun(@(x) extractFileNumber(x), file_names);
    [~, idx] = sort(file_nums);
    file_sortedaddress = fullfile(file_path,file_names(idx));

    % Load first image
    first_im = fullfile(file_path, file_names{1});
    t = Tiff(first_im, 'r');
    tifsize = gettifsize(t);
    nrows = tifsize(1);
    ncols = tifsize(2);

    % judge to merge single tifs or tif stacks
    if t.lastDirectory()
        [movie, nframes] = readsingletifs(file_sortedaddress,tifsize);
        fprintf('All single frame tifs movie loaded\n')
    else
        stacks_num = numel(file_sortedaddress);
        movie = [];
        nframes = 0;
        for i = 1:stacks_num
            fprintf('Loading %s\n',file_sortedaddress{i})
            [current_movie, current_nframes] = readstacktifs(file_sortedaddress{i},tifsize);
            movie = cat(2,movie,current_movie);
            nframes = nframes + current_nframes;
            fprintf('%s loaded\n',file_sortedaddress{i})
        end
        fprintf('All stacked frame tifs movie merged\n')
    end

else
    [~, ~, file_extension] = fileparts(file_path);
    switch file_extension
        case {'.tif','.tiff'}
            % Load first image
            t = Tiff(file_path, 'r');
            tifsize = gettifsize(t);
            nrows = tifsize(1);
            ncols = tifsize(2);
            t.close();
            [movie, nframes] = readstacktifs(file_path, tifsize);
            fprintf('Stacked frame tifs movie loaded\n')

        case '.bin'
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

        case '.mat'

            fprintf('Loading saved data...\n')
            file = load(file_path);
            movie = file.movie;
            ncols = file.ncols;
            nrows = file.nrows;
            nframes = file.nframes;
    end
end


end

function num = extractFileNumber(file_names)
    % Extract the numeric part of the filename
    % Handle both cases: filenames with only numbers and filenames with _number
    num = [];
    [~, name, ~] = fileparts(file_names);
    if all(isstrprop(name, 'digit'))
        num = str2double(name);
    else
        tokens = regexp(name, '([0-9]+)$', 'tokens');
        % tokens = regexp(name, '_([0-9]+)$', 'tokens');
        if ~isempty(tokens)
            num = str2double(tokens{1}{1});
        end
    end
end

function tifsize = gettifsize(t)
    tifsize = [t.getTag('ImageLength'),t.getTag('ImageWidth')];
end

function tifframe = gettifframe(t)
fprintf('Calculating frame...\n')
while true
    if t.lastDirectory()
        tifframe = t.currentDirectory();
        break
    else
        t.nextDirectory()
    end
end
end

function [movie, nframes]= readstacktifs(file_path,tifsize)
     % Initialize parameter
    t = Tiff(file_path, 'r');
    % 初始化帧数
    nframes = gettifframe(t);  % 初始为1，因为至少有一帧
    t.setDirectory(1);

    nrows = tifsize(1);
    ncols = tifsize(2);
    movie = zeros(nrows*ncols, nframes, 'uint16');
    prev_percentage = 0; % Initialize with -1 so the first update is always printed

    tic;
        for i = 1:nframes
            % Read TIFF movie
            % current_image = imread(file_path, i);
            % 跳到当前帧
            t.setDirectory(i);
            
            % 读取当前帧的图像
            current_image = t.read();
            movie(:,i) = uint16(reshape(current_image, nrows*ncols, 1));

            % Calculate and display progress if percentage changes
            current_percentage = floor(((i - 1) / nframes) * 100);
            if current_percentage > prev_percentage
                elapsed = toc;
                remaining = elapsed / ((i - 1) / nframes) - elapsed;
                fprintf('Processing %d/%d files (%d%% complete). Estimated time remaining: %.2f seconds\n', ...
                    i - 1, nframes, current_percentage, remaining);
                prev_percentage = current_percentage;
            end
        end
    % batchsize = 2000;
    % batchstart = 1;
    % while batchstart < nframes
    %     batchend = batchstart+batchsize;
    %     for i = batchstart:min(batchend,nframes)
    %         % Read TIFF movie
    %         % current_image = imread(file_path, i);
    %         % 跳到当前帧
    %         t.setDirectory(i);
    % 
    %         % 读取当前帧的图像
    %         current_image = t.read();
    %         movie(:,i) = uint16(reshape(current_image, nrows*ncols, 1));
    % 
    %         % Calculate and display progress if percentage changes
    %         current_percentage = floor((i / nframes) * 100);
    %         if current_percentage > prev_percentage
    %             elapsed = toc;
    %             remaining = elapsed / (i / nframes) - elapsed;
    %             fprintf('Processing %d/%d files (%d%% complete). Estimated time remaining: %.2f seconds\n', ...
    %                 i, nframes, current_percentage, remaining);
    %             prev_percentage = current_percentage;
    %         end
    %     end
    %     batchstart = batchend + 1;
    %     if batchstart > nframes
    %         break
    %     end
    % end
    t.close();
end

function [movie, nframes] = readsingletifs(file_sortedaddress, tifsize)

    % Loop through all TIF files and populate intensity time series parameter
    prev_percentage = 0; % Initialize with -1 so the first update is always printed
    nframes = numel(file_sortedaddress);
    nrows = tifsize(1);
    ncols = tifsize(2);
    movie = zeros(nrows*ncols,nframes, 'uint16');
    
    % Load batch of TIF files
    tic;
    for i = 1:nframes
        % Read the current image, store the image directly in 'movie'
        current_tif = file_sortedaddress{i};
        warning('off');
        t = Tiff(current_tif,'r');
        warning('on');
        current_image = t.read();
        movie(:,i) = uint16(reshape(current_image, nrows*ncols, 1));
        t.close();
        % Calculate and display progress if percentage changes
        current_percentage = floor((i / nframes) * 100);
        if current_percentage > prev_percentage
            elapsed = toc;
            remaining = (elapsed / i) * (nframes-i) ;
            fprintf('Processing %d/%d files (%d%% complete). Estimated time remaining: %.2f seconds\n', ...
                i, nframes, current_percentage, remaining);
            prev_percentage = current_percentage;
        end
    
    end
end






