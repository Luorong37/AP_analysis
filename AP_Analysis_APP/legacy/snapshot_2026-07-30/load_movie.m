function [movie,ncols,nrows,nframes,fext] = load_movie(file_path)

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
    [nrows, ncols, movie, nframes,fext] = readfoldertifs(file_path);

else
    [~, ~, file_extension] = fileparts(file_path);
    switch file_extension
        case {'.tif','.tiff'}
            % Load first image
            t = Tiff(file_path, 'r');
            tifsize = gettifsize(t);
            nrows = tifsize(2);
            ncols = tifsize(1);
            t.close();
            %[movie, nframes] = readstacktifs(file_path, tifsize);
            [movie, nframes] = readstacktifs(file_path);
            fprintf('Stacked frame tifs movie loaded\n')
            fext = '.tif';

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

            while ~feof(fileID)
                fileline =  fgetl(fileID);
                if contains(fileline,'nrow')
                    nrows = str2double(cell2mat(regexp(fileline,'\d*\.?\d*','match')));
                end
                if contains(fileline,'ncol')
                    ncols = str2double(cell2mat(regexp(fileline,'\d*\.?\d*','match')));
                end
            end
            fclose(fileID);
            %open file readBinMov
            % read file into tmp vector
            Movid = fopen(file_path);                  % open file
            Mov = fread(Movid, '*uint16', 'l');       % uint16, little endian
            fclose(Movid);                            % close file

            % reshape vector into appropriately oriented, 3D array
            nframes = length(Mov)/(nrows*ncols);
            movie = reshape(Mov, [nrows, ncols, nframes]);
            movie = permute(movie, [2 1 3]);
            % movie = reshape(movie, [nrows*ncols, nframes]);
            % fclose(fileID);
            fext = '.bin';

        case '.mat'

            fprintf('Loading saved data...\n')
            load(file_path);
            
            ncols = size(movie,1);
            nrows = size(movie,2);
            nframes = size(movie,3);
            fext = '.mat';
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

% function tifframe = gettifframe(t)
% count = 0;
% print_count = 0;
% try
% tifframe = numel(imfinfo(t));
% catch
%     ifcount = true;
% end
% 
% if isempty(tifframe)
%       ifcount = true;
% end
% 
% while ifcount
%     count = count+1;
%     if t.lastDirectory()
%         tifframe = t.currentDirectory();
%         fprintf(repmat('\b',1,print_count));   
%         fprintf('%d frames found\n', count);
%         break
%     else
%         fprintf(repmat('\b',1,print_count));   
%         print_count = fprintf('Calculating frame %d ...\n', count);
%         t.nextDirectory()
%     end
% end
% end
function tifframe = gettifframe(t)
try
    tifframe = numel(imfinfo(t.FileName));
catch
     tifframe = 1 ;
end

if tifframe == 1  % 如果 imfinfo 出错，使用循环方式估计帧数
    count = 1;
    print_count = 0;
    while ~t.lastDirectory()
        t.nextDirectory();
        count = count + 1;
        fprintf(repmat('\b',1,print_count));
        print_count = fprintf('Calculating frame %d ...\n', count);
    end
    tifframe = count;
    fprintf(repmat('\b',1,print_count));
    fprintf('%d frames found\n', count);
end
end

% % function [movie, nframes]= readstacktifs(file_path,tifsize)
%     t1 = tic;
%      % Initialize parameter
%     t = Tiff(file_path, 'r');
%     print_text = 0;
% 
%     % 初始化帧数
%     nframes = gettifframe(t);  % 初始为1，因为至少有一帧
%     t.setDirectory(1);
% 
%     nrows = tifsize(2);
%     ncols = tifsize(1);
%     movie = zeros(nrows*ncols, nframes, 'uint16');
%     prev_percentage = 0; % Initialize with -1 so the first update is always printed
% 
%     tic;
%         for i = 1:nframes
%             % Read TIFF movie
%             % current_image = imread(file_path, i);
%             % 跳到当前帧
%             t.setDirectory(i);
% 
%             % 读取当前帧的图像
%             current_image = t.read();
%             movie(:,i) = uint16(reshape(current_image, nrows*ncols, 1));
% 
%             % Calculate and display progress if percentage changes
%             current_percentage = floor(((i - 1) / nframes) * 100);
%             if current_percentage > prev_percentage
%                 elapsed = toc;
%                 remaining = elapsed / ((i - 1) / nframes) - elapsed;
%                 fprintf(repmat('\b',1,print_text));   
%                 print_text = fprintf('Processing %d/%d files (%d%% complete). Estimated time remaining: %.2f seconds\n', ...
%                     i - 1, nframes, current_percentage, remaining);
%                 prev_percentage = current_percentage;
%             end
%         end
%     % batchsize = 2000;
%     % batchstart = 1;
%     % while batchstart < nframes
%     %     batchend = batchstart+batchsize;
%     %     for i = batchstart:min(batchend,nframes)
%     %         % Read TIFF movie
%     %         % current_image = imread(file_path, i);
%     %         % 跳到当前帧
%     %         t.setDirectory(i);
%     % 
%     %         % 读取当前帧的图像
%     %         current_image = t.read();
%     %         movie(:,i) = uint16(reshape(current_image, nrows*ncols, 1));
%     % 
%     %         % Calculate and display progress if percentage changes
%     %         current_percentage = floor((i / nframes) * 100);
%     %         if current_percentage > prev_percentage
%     %             elapsed = toc;
%     %             remaining = elapsed / (i / nframes) - elapsed;
%     %             fprintf('Processing %d/%d files (%d%% complete). Estimated time remaining: %.2f seconds\n', ...
%     %                 i, nframes, current_percentage, remaining);
%     %             prev_percentage = current_percentage;
%     %         end
%     %     end
%     %     batchstart = batchend + 1;
%     %     if batchstart > nframes
%     %         break
%     %     end
%     % end
% %     t.close();
% %     t2 = toc(t1); % Get the elapsed time
% %     fprintf('Finished loading after %d s, ',round(t2))
% % end
function [movie, nframes] = readstacktifs(file_path, asVector)
% READSTACKTIFS - Robust reader for multi-frame TIFF stacks.
% movie: [nrows x ncols x nframes] by default, or [(nrows*ncols) x nframes] if asVector=true
% Also returns ncols, nrows, nframes for convenience.

    if nargin < 2
        asVector = false; % 默认返回 3D
    end

    info    = imfinfo(file_path);
    nframes = numel(info);
    nrows  = info(1).Height;
    ncols   = info(1).Width;
    % spp     = isfield(info(1),'SamplesPerPixel') * info(1).SamplesPerPixel;
    % if spp == 0, spp = 1; end

    % % 先读一帧以确定类和通道处理
    first = imread(file_path, 1, 'Info', info);
    % % 若多通道，默认转灰度（RGB 用 rgb2gray；>3 通道取第1通道，按需修改）
    % if ndims(first) == 3
    %     if size(first,3) == 3
    %         first = rgb2gray(first);
    %     else
    %         first = first(:,:,1);
    %     end
    % end
    className = class(first);

    % 预分配
    movie = zeros(nrows, ncols, nframes, className);
    pool = gcp('nocreate');
    if isempty(pool) || ~contains(class(pool), 'ProcessPool')
        fprintf('Reading multi-frame TIFF via serial read... %d frames\n', nframes);
        t0 = tic;
        progress_step = max(1, ceil(nframes / 100));
        for frame_idx = 1:nframes
            movie(:,:,frame_idx) = imread(file_path, frame_idx, 'Info', info);
            if mod(frame_idx, progress_step) == 0 || frame_idx == nframes
                elapsed_s = toc(t0);
                remaining_s = elapsed_s / frame_idx * (nframes - frame_idx);
                fprintf('Loaded %d/%d frames (%d%%) | elapsed %.1f s | ETA %.1f s\n', ...
                    frame_idx, nframes, floor(frame_idx / nframes * 100), elapsed_s, remaining_s);
            end
        end
        fprintf('Finished loading after %d s, ', round(toc(t0)));
        if asVector
            movie = reshape(movie, nrows*ncols, nframes);
        end
        return;
    end
    
    
    % 直接读单个tiff
    % print_text = fprintf('Reading multi-frame TIFF containing %d frames\n, please wait about %.0f seconds', nframes, nframes/1000 + 1);
    t0 = tic;
    % t  = Tiff(file_path, 'r');
    % for k = 1:nframes
    %     img = t.read();
    %     movie(:,:,k) = img;        % 预分配好类型/大小
    %     if k < nframes, t.nextDirectory(); end
    % end
    % t.close();
    
    % 并行分段读取
    % 分配工作负载：每个worker读取的帧的范围
    poolsize = pool.NumWorkers;
    frames_per_worker = ceil(nframes / poolsize);
    movie_cell = cell(1, poolsize);

    % 进度：worker -> 客户端
    q        = parallel.pool.DataQueue;
    nDone    = 0;    % 累计已完成数
    print_text = 0;
    prevPcnt = 0;    % 上一次已打印的百分比
    afterEach(q, @updateProgress);  % 嵌套函数，能直接用 nframes/t0/nDone/prevPcnt

    % --- 关闭警告（客户端 + 各 worker）并读取tif---
    
    % 1) 关闭
    pool = gcp;                        % 确保已开池
    warning('off','all');              % 客户端
    spmd, warning('off','all');  end    % worker

    parfor k = 1:poolsize
        % 每个worker的起始和结束帧
        start_frame = (k - 1) * frames_per_worker + 1;
        end_frame = min(k * frames_per_worker, nframes);
        % 创建每个 worker 的独立矩阵
        worker_movie = zeros(nrows, ncols, end_frame - start_frame + 1, className);
        % 读取该worker负责的帧
        t = Tiff(file_path, 'r');
        for frame = start_frame:end_frame
            t.setDirectory(frame); % 设置要读取的目录
            img = t.read();
            worker_movie(:,:,frame - start_frame + 1) = img;
            send(q, 1);  % 每完成一个发一个通知（内容可忽略）
        end
        t.close();
        % 将 worker 结果进行存储
        movie_cell{k} = worker_movie;
    end

    for i = 1:poolsize
        start_frame = (i - 1) * frames_per_worker + 1;
        end_frame = min(i * frames_per_worker, nframes);
        movie(:,:,start_frame:end_frame) = movie_cell{i};
    end

    % 2) 恢复原先警告状态
    spmd, warning('on','all');  end    % worker
    warning('on','all');               % 客户端

    % -------- 嵌套回调：只在百分比增加时打印 --------
    function updateProgress(~)
        nDone  = nDone + 1;
        pcnt   = floor(nDone / nframes * 100);
        if pcnt > prevPcnt
            elapsed   = toc(t0);
            remaining = (elapsed / nDone) * (nframes - nDone);
            fprintf(repmat('\b',1,print_text));    % 擦除上次的进度信息
            print_text = fprintf('Processing %d/%d (%d%% complete). Estimated time remaining: %.1f s\n', ...    % 输出进度信息
                    nDone, nframes, pcnt, remaining);
            prevPcnt = pcnt;
        end
    end

    % elapsed = toc(t0);
    % fprintf(repmat('\b',1,print_text));    % 擦除上次的进度信息
    % fprintf('Loaded %d frames in %.2f s\n', nframes, elapsed);
    t1 = toc(t0);
    fprintf('Finished loading after %d s, ',round(t1));

    

    % ##Old code
    % print_text = 0;
    % t1 = tic;
    % parfor i = 2:nframes
    %     f = imread(file_path, i, 'Info', info);
    % 
    %     % if ndims(f) == 3
    %     %     if size(f,3) == 3
    %     %         f = rgb2gray(f);
    %     %     else
    %     %         f = f(:,:,1);
    %     %     end
    %     % end % 强制转换到首帧类型，避免某些帧类名不同（极少见） if ~strcmp(class(f), className)
    %     %     f = cast(f, className);
    %     % end
    %     movie(:,:,i) = f;
    %     elapsed = toc;
    %     current_percentage = floor((i / nframes) * 100);
    %     remaining = (elapsed / i) * (nframes-i) ;
    %     print_text = fprintf('Processing %d/%d files (%d%% complete). Estimated time remaining: %.2f seconds\n', ...
    %             i, nframes, current_percentage, remaining);
    %     prev_percentage = current_percentage;
    % end

    % 可选：返回列向量堆叠
    if asVector
        movie = reshape(movie, nrows*ncols, nframes);
    end
end

function [movie, nframes] = readsingletifs(file_sortedaddress, tifsize)
    t0 = tic;

    nframes = numel(file_sortedaddress);
    nrows   = tifsize(2);
    ncols   = tifsize(1);
    movie   = zeros(nrows*ncols, nframes, 'uint16');
    pool = gcp('nocreate');
    if isempty(pool) || ~contains(class(pool), 'ProcessPool')
        fprintf('Reading single TIFFs via serial read... %d frames\n', nframes);
        t0 = tic;
        for frame_idx = 1:nframes
            current_image = imread(file_sortedaddress{frame_idx});
            movie(:,frame_idx) = uint16(reshape(current_image, nrows*ncols, 1));
        end
        fprintf('Loaded %d frames in %.2f s\n', nframes, toc(t0));
        return;
    end

    fprintf('Reading single TIFF via read (parallel)… %d frames\n', nframes);

    % 进度：worker -> 客户端
    q        = parallel.pool.DataQueue;
    nDone    = 0;    % 累计已完成数
    print_text = 0;
    prevPcnt = 0;    % 上一次已打印的百分比
    afterEach(q, @updateProgress);  % 嵌套函数，能直接用 nframes/t0/nDone/prevPcnt

    % --- 关闭警告（客户端 + 各 worker）并读取tif---
    
    % 1) 关闭
    pool = gcp;                        % 确保已开池
    warning('off','all');              % 客户端
    spmd, warning('off','all');  end    % worker

    % 并行读取
    parfor i = 1:nframes
        current_tif  = file_sortedaddress{i};
        t = Tiff(current_tif,'r');
        current_image = t.read(); 
        t.close();
        movie(:,i) = uint16(reshape(current_image, nrows*ncols, 1));

        send(q, 1);  % 每完成一个发一个通知（内容可忽略）
    end

    % 2) 恢复原先警告状态
    spmd, warning('on','all');  end    % worker
    warning('on','all');               % 客户端

    fprintf('Loaded %d frames in %.2f s\n', nframes, toc(t0));

    % -------- 嵌套回调：只在百分比增加时打印 --------
    function updateProgress(~)
        nDone  = nDone + 1;
        pcnt   = floor(nDone / nframes * 100);
        if pcnt > prevPcnt
            elapsed   = toc(t0);
            remaining = (elapsed / nDone) * (nframes - nDone);
            fprintf(repmat('\b',1,print_text));    % 擦除上次的进度信息
            print_text = fprintf('Processing %d/%d (%d%% complete). Estimated time remaining: %.1f s\n', ...    % 输出进度信息
                    nDone, nframes, pcnt, remaining);
            prevPcnt = pcnt;
        end
    end
end

% ##Old code
% function [movie, nframes] = readsingletifs(file_sortedaddress, tifsize)
%     t1 = tic;
%     % Loop through all TIF files and populate intensity time series parameter
%     nframes = numel(file_sortedaddress);
%     nrows = tifsize(2);
%     ncols = tifsize(1);
%     movie = zeros(nrows*ncols,nframes, 'uint16');
%     print_text = 0;
% 
%     % 预分配
%     fprintf('Reading single TIFF via read (parallel)… %d frames\n', nframes);
% 
%     % 进度队列（worker -> 客户端）
%     q = parallel.pool.DataQueue;
%     nDone = 0;
%     prev_percent = 0;
%     afterEach(q, @updateProgress);
% 
%     parfor i = 1:nframes
%        % Read the current image, store the image directly in 'movie'
%         current_tif = file_sortedaddress{i};
%         warning('off');
%         t = Tiff(current_tif,'r');
%         warning('on');
%         current_image = t.read();
%         movie(:,i) = uint16(reshape(current_image, nrows*ncols, 1));
%         t.close();
% 
%        send(q, 1); % 报告完成1个
%     end
% 
%     elapsed = toc(t1);
%     fprintf('Loaded %d frames in %.2f s\n', nframes-1, elapsed);
% 
%     % --- 客户端进度更新函数 ---
%     function updateProgress(i)
%         persistent prev_percent_inner
%         if isempty(prev_percent_inner)
%             prev_percent_inner = 0;
%         end
% 
%         nframes_local = evalin('base','nframes');
%         percent = floor(i / nframes_local * 100);
% 
%         if percent > prev_percent_inner
%             elapsed = toc(evalin('base','t0'));
%             remaining = (elapsed / i) * (nframes_local - i);
%             fprintf('Processing %d/%d (%d%% complete). ETA: %.1f s\n', ...
%                 i, nframes_local, percent, remaining);
%             prev_percent_inner = percent;
%         end
%     end
%     t2 = toc(t1);
%     fprintf('Finished loading after %d s, ',round(t2))
% end


function [nrows, ncols, movie, nframes,fext] = readfoldertifs(file_path)
% Get all TIFF file names in folder
file_list = dir(fullfile(file_path, '*.tif'));
file_names = {file_list.name};
if ~isempty(file_names)
% Sort file names
% other wise, will be like [4, 40, 400, 4000, 4001]
if length(file_names) > 400
try
file_nums = cellfun(@(x) extractFileNumber(x), file_names,'UniformOutput',false);
[~, idx] = sort(cell2mat(file_nums));
catch
warning('illegal stacked name.')
idx = 1:length(file_names);
end
else
    idx = 1:length(file_names);
end
file_sortedaddress = fullfile(file_path,file_names(idx));

% Load first image
first_im = fullfile(file_path, file_names{1});
t = Tiff(first_im, 'r');
tifsize = gettifsize(t);
nrows = tifsize(2);
ncols = tifsize(1);

% judge to merge single tifs or tif stacks
if t.lastDirectory()
    % thies means in a folder, the first tif is a single tif, means all tifs are single tif 
    [movie, nframes] = readsingletifs(file_sortedaddress,tifsize);
    fprintf('All single frame tifs movie loaded\n')
else
    stacks_num = numel(file_sortedaddress);
    movie = [];
    nframes = 0;
    for i = 1:stacks_num
        fprintf('Loading %s, ',file_sortedaddress{i})
        [current_movie, current_nframes] = readstacktifs(file_sortedaddress{i},tifsize);
        movie = cat(2,movie,current_movie);
        nframes = nframes + current_nframes;
        fprintf('%s loaded\n',file_sortedaddress{i})
    end
    fprintf('All stacked frame tifs movie merged\n')
end
fext = '.tif';
else
    file_list = dir(fullfile(file_path, '*.mat'));
    file_names = {file_list.name};
    load(fullfile(file_path, file_names{1}))
    nrows = size(movie,2);
    ncols = size(movie,1);
    nframes = size(movie,3);
    fext = '.mat';
end
end
% %function [nrows, ncols, movie, nframes] = readfoldertifs(imageDirectory2)
%     % 获取文件夹下所有 tif 文件
%     imageFiles = dir(fullfile(imageDirectory2, '*.tif'));
%     numImages = numel(imageFiles);
% 
%     if numImages == 0
%         error('No .tif files found in the specified directory.');
%     end
% 
%     % 读取第一个文件的信息以获取图像维度和数据类型
%     filename = fullfile(imageDirectory2, imageFiles(1).name);
%     imageInfo = imfinfo(filename);
%     imageData = imread(filename);
%     [nrows, ncols] = size(imageData);
%     dataSize = [imageInfo.Width, imageInfo.Height];
%     dataType = class(imageData);
% 
%     % 预分配 cell 数组用于 parfor 存储
%     Y = cell(1, numImages);
% 
%     % 并行读取每张图像数据
%     parfor i = 1:numImages
%         fname = fullfile(imageDirectory2, imageFiles(i).name);
%         info = imfinfo(fname);
%         try
%             offset = info.StripOffsets;
%             mm = memmapfile(fname, ...
%                 'Format', {dataType, dataSize, 'Data'}, ...
%                 'Offset', offset(1));
%             % 注意转置图像维度（必要时）
%             Y{i} = mm.Data.Data';
%         catch ME
%             warning('File %s cannot be read with memmapfile. Switching to imread. Reason: %s', fname, ME.message);
%             im = imread(fname);
%             Y{i} = im';
%         end
%     end
% 
%     % 将图像堆栈拼接成 3D 矩阵
%     movie3D = cat(3, Y{:});
% 
%     % 重塑为 (nrows*ncols, nframes)
%     movie = reshape(movie3D, [], numImages);
%     nframes = numImages;
% end
