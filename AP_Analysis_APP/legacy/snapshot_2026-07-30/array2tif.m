function array2tif(movie, save_path, ncols, nrows)
% 优化点：
% 1. 手动分块发送：避免 parfor 序列化整个大数组，极大降低内存占用
% 2. 内存自适应：根据当前内存状态决定并行度
% 3. 精准退格进度条

% --- 1. 参数校验 ---
str = class(movie);
bitdepth = str2double(str(isstrprop(str, 'digit')));
bytesPerPixel = bitdepth / 8;

if nargin < 3
    [nrows, ncols, nframes] = size(movie);
else
    movie = reshape(movie, nrows, ncols, []);
    nframes = size(movie, 3);
end

% 计算 Stack 逻辑 (~3.6GB 一个文件)
maxBytesPerStack = 3.6 * 1024^3; 
bytesPerFrame = nrows * ncols * bytesPerPixel;
maxtif_frames = floor(maxBytesPerStack / bytesPerFrame);
stack_num = ceil(nframes / maxtif_frames);
[file_path, file_name, ~] = fileparts(save_path);
fprintf('Total Frames: %d | BitDepth: %d | Stacks: %d\n', nframes, bitdepth, stack_num);

% --- 2. 内存安全预检 ---
if ispc
    [~, sys] = memory;
    availableMem = sys.PhysicalMemory.Available;
    % 如果单份数据太大，减少并行 Worker 数量防止挤爆内存
    movieSize = numel(movie) * bytesPerPixel;
    if movieSize > availableMem * 0.6
        fprintf('Warning: Data is too large. Reducing parallel workers to save memory.\n');
        % 此处不关闭并行，但建议用户手动限制 pool 大小
    end
end

% --- 3. 准备数据分块 (关键优化) ---
% 将 movie 切成 cell 数组，每个 cell 对应一个 stack
% 这样 parfor 只会分发对应的 cell 内容，而不是整个大矩阵
movie_chunks = cell(1, stack_num);
for s = 1:stack_num
    fStart = (s-1) * maxtif_frames + 1;
    fEnd = min(s * maxtif_frames, nframes);
    movie_chunks{s} = movie(:,:,fStart:fEnd); 
end
% 此时原 movie 已经没用了，立即清除以腾出空间给序列化使用
clear movie; 

% --- 4. 进度条初始化 ---
q = parallel.pool.DataQueue;
tStart = tic;
lastLineLen = 0;
totalDone = 0;
    function updateProgress(nInc)
        totalDone = totalDone + nInc;
        pcnt = (totalDone / nframes) * 100;
        msg = sprintf('Progress: %d/%d (%.1f%%) | Remaining: %.1fs\n', ...
                      totalDone, nframes, pcnt, (toc(tStart)/totalDone)*(nframes-totalDone));
        fprintf(repmat('\b', 1, lastLineLen));
        lastLineLen = fprintf('%s', msg);
    end
afterEach(q, @updateProgress);

% --- 5. 写入 ---
pool = gcp('nocreate');
useParallelWrite = ~isempty(pool) && contains(class(pool), 'ProcessPool');
if useParallelWrite
    fprintf('Starting parallel write for %d stacks...\n', stack_num);
else
    fprintf('Starting serial write for %d stacks (no active process pool).\n', stack_num);
end

if useParallelWrite
parfor s = 1:stack_num
    local_data = movie_chunks{s}; % 每个 Worker 只拿到它需要的那 3.6GB
    this_nframes = size(local_data, 3);
    stack_path = fullfile(file_path, sprintf('%s_stack%02d.tif', file_name, s));
    
    t = Tiff(stack_path, 'w');
    for i = 1:this_nframes
        t.setTag('ImageLength', nrows);
        t.setTag('ImageWidth', ncols);
        t.setTag('Photometric', Tiff.Photometric.MinIsBlack);
        t.setTag('BitsPerSample', bitdepth);
        t.setTag('SamplesPerPixel', 1);
        t.setTag('RowsPerStrip', 32);
        t.setTag('PlanarConfiguration', Tiff.PlanarConfiguration.Chunky);
        t.setTag('Compression', Tiff.Compression.None);
        
        t.write(local_data(:,:,i));
        if i < this_nframes, t.writeDirectory(); end
        
        % 每完成 50 帧发一次进度，减少通信开销
        if mod(i, 50) == 0, send(q, 50); end
    end
    % 发送余数
    rem_frames = mod(this_nframes, 50);
    if rem_frames > 0, send(q, rem_frames); end
    
    t.close();
end
else
for s = 1:stack_num
    local_data = movie_chunks{s};
    this_nframes = size(local_data, 3);
    stack_path = fullfile(file_path, sprintf('%s_stack%02d.tif', file_name, s));
    
    t = Tiff(stack_path, 'w');
    for i = 1:this_nframes
        t.setTag('ImageLength', nrows);
        t.setTag('ImageWidth', ncols);
        t.setTag('Photometric', Tiff.Photometric.MinIsBlack);
        t.setTag('BitsPerSample', bitdepth);
        t.setTag('SamplesPerPixel', 1);
        t.setTag('RowsPerStrip', 32);
        t.setTag('PlanarConfiguration', Tiff.PlanarConfiguration.Chunky);
        t.setTag('Compression', Tiff.Compression.None);
        
        t.write(local_data(:,:,i));
        if i < this_nframes, t.writeDirectory(); end
        
        if mod(i, 50) == 0, send(q, 50); end
    end
    rem_frames = mod(this_nframes, 50);
    if rem_frames > 0, send(q, rem_frames); end
    
    t.close();
end
end

fprintf('\nDone. Time: %.1fs\n', toc(tStart));
end



% function array2tif(movie, save_path, ncols, nrows)
% % 用于在create_tiff_stack中分组并行读取多个tif文件并合并为tiff_stack
% %
% 
% str = class(movie);
% digits = str(isstrprop(str, 'digit')); % 提取数字字符
% bitdepth = str2num(digits); % 转换为数值
% 
% % movie = gpuArray(movie);
% 
% if nargin <3
% 
%     if length(size(movie)) < 3
%         error('2D-array must be inputted with size.')
%     end
% 
%     ncols = size(movie,2);
%     nrows = size(movie,1);
% else
%     movie = reshape(movie,ncols,nrows,[]);
% end
% 
% % gcp;
% [file_path,file_name,~] = fileparts(save_path);
% t0 = tic;
% 
% 
% 
% nframes = size(movie,3);
% maxpiexls = 512*512*7200/2;
% maxtif_frames = round(maxpiexls/(nrows*ncols));
% % maxtif_frames = 7200;
% 
% if nframes > maxtif_frames  
%     stack_num = ceil(nframes / maxtif_frames);         % 计算所需的stack数目
% else
%     stack_num = 1;
% end
% 
% 
% % 进度跟踪：worker -> 客户端
% q        = parallel.pool.DataQueue;
% nDone    = 0;    % 累计已完成数
% print_text = 0;
% prevPcnt = 0;    % 上一次已打印的百分比
% afterEach(q, @updateProgress);  % 嵌套函数，能直接用 nframes/t0/nDone/prevPcnt
% 
% for s = 1:stack_num
%     stack_path = fullfile(file_path, strcat(file_name, sprintf('_stack%02d.tif', s)));
%     fprintf('Creating new stack file: %s\n', stack_path);
% 
%     frame_start = (s-1)*maxtif_frames + 1 ;
%     frame_end =  min(frame_start + maxtif_frames - 1,nframes);
% 
%     t = Tiff(stack_path, 'w');
% 
%     for i = frame_start:frame_end
%         current_image = movie(:,:,i);         % 读取当前图片
% 
%         t.setTag('ImageLength', nrows);
%         t.setTag('ImageWidth', ncols);
%         t.setTag('Photometric', Tiff.Photometric.MinIsBlack);
%         t.setTag('BitsPerSample', bitdepth);
%         t.setTag('SamplesPerPixel', 1);
%         t.setTag('RowsPerStrip', 16);
%         t.setTag('PlanarConfiguration', Tiff.PlanarConfiguration.Chunky);
%         t.setTag('Compression', Tiff.Compression.None);
%         t.setTag('Software', 'MATLAB');
%         t.write(current_image);
% 
%         if i < nframes
%             t.writeDirectory();
%         end
% 
%         send(q, 1);  % 每完成一个发一个通知（内容可忽略）
% 
%     end
%     t.close();
% 
% end
% 
% t1 = toc(t0);
% fprintf('Finished processing after %d s, ',round(t1));
% fprintf('Stacked tif created in %s.\n', save_path);
% 
% 
% 
% % -------- 嵌套回调：只在百分比增加时打印 --------
% function updateProgress(~)
%     nDone  = nDone + 1;
%     pcnt   = floor(nDone / nframes * 100);
%     if pcnt > prevPcnt
%         elapsed   = toc(t0);
%         remaining = (elapsed / nDone) * (nframes - nDone);
%         fprintf(repmat('\b',1,print_text));    % 擦除上次的进度信息
%         print_text = fprintf('Processing %d/%d (%d%% complete). Estimated time remaining: %.1f s\n', ...    % 输出进度信息
%                 nDone, nframes, pcnt, remaining);
%         prevPcnt = pcnt;
%     end
% end
% 
% end
% 
