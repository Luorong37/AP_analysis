function [quick_map] = create_map(movie, nrows, ncols, bin, mode)

% ----------Write by Liu-Yang Luorong and ChatGPT----------
% ----------POWERED by Zoulab in Peking University----------
% Date: 23.11.16
% MATLAB Version: R2022b
% Function to process movie data by binning, normalizing, and generating a quick_map.
%
% This function create a map by the extreme value for each pixel.
% It takes a 2D movie array and applies several processing steps to transform the movie data.
% It first adjusts the movie size to be divisible by a given bin size, then performs binning and normalization.
% The output is a quick_map representing normalized, binned data of the movie.
%
% Parameters:
% movie - A 2D array representing the movie, with dimensions [num_cols*num_rows, num_frames].
% bin - Integer representing the binning factor.
% num_rows - Number of rows in the movie.
% num_cols - Number of columns in the movie.
%
% Processing Steps:
% 1. Adjust movie size to make it divisible by the bin size.
% 2. Reshape the movie into a 5D array for binning.
% 3. Calculate the average intensity for each binned pixel across all frames.
% 4. Apply photobleaching correction using an exponential fitting on the average intensity trace.
% 5. Compute a quick_map based on normalized intensity extreme values.
% 6. Resize the quick_map to the original movie dimensions.
%
% Output:
% quick_map - A 2D array representing the processed and normalized data, resized to the original movie dimensions.
%
% Example:
% quick_map = create_map(movie, bin, num_rows, num_cols);
%
% Notes:
% - The function assumes the input movie is a 2D array with the third dimension representing time (frames).
% - Photobleaching correction is defined as single exponential function.
% - The map only find the extreme value in frames of binned pixels which may lost plateau.
%
% See also RESHAPE, MEAN, FIT, STD, IMRESIZE.

if nargin < 5
    mode = 'voltage'; % defined create sensitivity map
end

% cut the size
cut = false;
if mod(nrows,bin) ~= 0
    nrows = nrows - mod(nrows,bin);
    cut = true;
end
if mod(ncols,bin) ~= 0
    ncols = ncols - mod(ncols,bin);
    cut = true;
end
if cut == true
    temp = reshape(movie,ncols,nrows,[]);
    movie = temp(1:ncols,1:nrows,:);
end
movie_size = size(movie);
nframe = movie_size(end);

% Initialize matrix to store cell labels
quick_map = zeros(ncols * nrows,1);
% quick_map = zeros(ncols/bin * nrows/bin,1);
h = ones(3,3);
h(5) = 0;
movie_binned = reshape(movie,ncols,nrows,[]);
% movie_sum = zeros(size(movie_3D));
% Prefer the original parfor path when the caller has already opened a
% parallel pool. Pool startup is handled by the top-level analysis script.
useParallelMap = false;
if exist('gcp', 'file') == 2
    try
        useParallelMap = ~isempty(gcp('nocreate'));
    catch
        useParallelMap = false;
    end
end

if useParallelMap
    parfor i = 1:nframe
        current_frame = movie_binned(:,:,i)./8;
        sum_frame = imfilter(current_frame,h,'conv');
        movie_binned(:,:,i) = sum_frame;
    end
else
    for i = 1:nframe
        current_frame = movie_binned(:,:,i)./8;
        sum_frame = imfilter(current_frame,h,'conv');
        movie_binned(:,:,i) = sum_frame;
    end
end

if any(strcmpi(mode, {'voltage_chunked', 'paper'}))
    quick_map = create_voltage_chunked_map(movie_binned);
    return;
end

movie_binned = reshape(movie_binned,nrows*ncols,[]);
% reshape to high dimension for each bin
% movie_5D = reshape(movie,bin,ncols/bin,bin,nrows/bin,[]);
% % average across x y
% movie_ave = squeeze(mean(mean(movie_5D,1),3));
%
% % reshape to binned
% movie_binned = reshape(movie_ave,nrows/bin*ncols/bin,[]);
npixels = size(movie_binned,1);
%
% % Apply avg_intensity curve to correct photobleaching from each pixel (recommend)
%
%
% avg_trace = mean(movie_binned,1);
% x = (1:nframe)';
% f = fit(x,avg_trace','poly1');
% fitted_curves = x*f.p1+f.p2;
% f = fit(x,avg_trace','exp2');
% fitted_curves = f.a*exp(f.b*x)+f.c*exp(f.d*x);
% movie_binned_corrected = movie_binned ./ fitted_curves';
% movie_binned_corrected = wdenoise(movie_binned_corrected',DenoisingMethod="FDR")';
%
% print_count = 0;
% % % 定义滤波器参数
% fs = 400;                % 采样频率，单位为帧/秒（根据实际情况调整）
% cutoff_freq_time = 1/10; % 截止频率，单位为Hz（根据需要调整）
% order_time = 2;          % 滤波器阶数
%
% % 设计巴特沃斯滤波器
% [b_time, a_time] = butter(order_time, cutoff_freq_time/(fs/2), 'low');
%
% % 设置填充长度，将信号在头尾各填充5%的长度
% padlength = round(0.05 * size(movie_binned, 2));
%
% % 对信号进行填充，以减少边界效应
% padded_movie = [repmat(movie_binned(:, 1), 1,padlength), ...
%                  movie_binned, ...
%                  repmat(movie_binned(:,end), 1,padlength)];
%
% % 应用高通滤波器到所有像素的时间序列
% % filtfilt 会自动对每一列（对应一个像素的时间序列）进行滤波
%
% movie_binned_corrected = filtfilt(b_time, a_time, padded_movie');
% movie_binned_corrected = movie_binned_corrected';
% movie_binned_corrected = movie_binned_corrected(:,padlength+1:end-padlength);
% movie_binned_corrected = detrend(movie_binned',4)';

%  movie_binned_corrected_3D =  reshape(movie_binned_corrected,ncols/bin,nrows/bin,[]);
% select_ROI(movie_binned_corrected, ncols/bin, nrows/bin, [], [])

% 进度跟踪：worker -> 客户端
% 
% t0 = tic;
% q        = parallel.pool.DataQueue;
% nDone    = 0;    % 累计已完成数
% print_text = 0;
% prevPcnt = 0;    % 上一次已打印的百分比
% afterEach(q, @updateProgress);  % 嵌套函数，能直接用 nframes/t0/nDone/prevPcnt




if useParallelMap
parfor i = 1:npixels
    % fprintf('Processing %d / %d\n', i, npixels );
    % fprintf(repmat('\b',1,print_count));
    % % fprintf('Calculating %.2f %% \n', i/npixels*100);
    % print_count = fprintf('Calculating %.2f %% \n', i/npixels*100);
    % pixel_trace_corrected = movie_binned(i, :);
    % pixel_trace_corrected = movie_binned(i, :);
    % pixel_trace_corrected = wdenoise(movie_binned_corrected(i, :));
    % pixel_trace_corrected = detrend(double(movie_binned(i, :)),1);
    % pixel_trace_corrected =  movie_binned(i, :) - movie_binned_corrected (i, :);
    % baseline =mean(pixel_trace_corrected);
    % baseline = mean(pixel_trace_corrected); % filted
    if strcmpi(mode,'voltage')
        pixel_trace_corrected = detrend(double(movie_binned(i, :)),1);
        % abstrace = abs(pixel_trace_corrected - baseline);
        abstrace = abs(pixel_trace_corrected);
        maxpointabs =  max(abstrace);
        maxpointindex = abstrace == maxpointabs;
        pointdff  = pixel_trace_corrected(maxpointindex);
        % pointdff = ((maxpoint-baseline)/baseline);
        quick_map(i) = pointdff(1);
        % quick_map(i) = sign(pointdff(1)) .* (1 - exp(-k * abs(pointdff(1)*n).^n)) / (1 - exp(-k));

    elseif strcmpi(mode,'calcium')
        pixel_trace_corrected = detrend(double(movie_binned(i, :)),2);
        quick_map(i) = std(pixel_trace_corrected);
        % quick_map(i) = (max(pixel_trace_corrected)-baseline) /std(pixel_trace_corrected);
    end
    % send(q, 1);  % 每完成一个发一个通知（内容可忽略）
end
else
for i = 1:npixels
    % fprintf('Processing %d / %d\n', i, npixels );
    % fprintf(repmat('\b',1,print_count));
    % % fprintf('Calculating %.2f %% \n', i/npixels*100);
    % print_count = fprintf('Calculating %.2f %% \n', i/npixels*100);
    % pixel_trace_corrected = movie_binned(i, :);
    % pixel_trace_corrected = movie_binned(i, :);
    % pixel_trace_corrected = wdenoise(movie_binned_corrected(i, :));
    % pixel_trace_corrected = detrend(double(movie_binned(i, :)),1);
    % pixel_trace_corrected =  movie_binned(i, :) - movie_binned_corrected (i, :);
    % baseline =mean(pixel_trace_corrected);
    % baseline = mean(pixel_trace_corrected); % filted
    if strcmpi(mode,'voltage')
        pixel_trace_corrected = detrend(double(movie_binned(i, :)),1);
        % abstrace = abs(pixel_trace_corrected - baseline);
        abstrace = abs(pixel_trace_corrected);
        maxpointabs =  max(abstrace);
        maxpointindex = abstrace == maxpointabs;
        pointdff  = pixel_trace_corrected(maxpointindex);
        % pointdff = ((maxpoint-baseline)/baseline);
        quick_map(i) = pointdff(1);
        % quick_map(i) = sign(pointdff(1)) .* (1 - exp(-k * abs(pointdff(1)*n).^n)) / (1 - exp(-k));

    elseif strcmpi(mode,'calcium')
        pixel_trace_corrected = detrend(double(movie_binned(i, :)),2);
        quick_map(i) = std(pixel_trace_corrected);
        % quick_map(i) = (max(pixel_trace_corrected)-baseline) /std(pixel_trace_corrected);
    end
    % send(q, 1);  % 每完成一个发一个通知（内容可忽略）
end
end

% converse to Z score
% baseline = mean(quick_map);
% quick_map = (quick_map - baseline) ./ std(quick_map);
% quick_map = reshape(quick_map,ncols/bin,nrows/bin, 1);
quick_map = reshape(quick_map,ncols,nrows, []);
% quick_map = map;
quick_map(quick_map>0) = quick_map(quick_map>0) - mean(quick_map(quick_map>0));
quick_map(quick_map<0) = quick_map(quick_map<0) - mean(quick_map(quick_map<0));
% imagesc(quick_map)
% colorbar
% quick_map = imresize(quick_map,[ncols,nrows] );


% % -------- 嵌套回调：只在百分比增加时打印 --------
% function updateProgress(~)
% nDone  = nDone + 1;
% pcnt   = floor(nDone / npixels * 100);
% if pcnt > prevPcnt
%     elapsed   = toc(t0);
%     remaining = (elapsed / nDone) * (npixels - nDone);
%     fprintf(repmat('\b',1,print_text));    % 擦除上次的进度信息
%     print_text = fprintf('Processing %d/%d (%d%% complete). Estimated time remaining: %.1f s\n', ...    % 输出进度信息
%         nDone, npixels, pcnt, remaining);
%     prevPcnt = pcnt;
% end
% end

end

function quick_map = create_voltage_chunked_map(movie_3d)
% Paper-inspired summary map for voltage imaging.
% It aggregates short temporal chunks using average and max-minus-median images.

[ncols, nrows, nframe] = size(movie_3d);
chunk_size = 50;
min_chunk_frames = min(10, nframe);
spatial_sigma = 1.0;
avg_weight = 0.35;
active_weight = 0.65;

aggregate_map = zeros(ncols, nrows, 'single');
active_map = zeros(ncols, nrows, 'single');
chunk_counter = 0;

for start_idx = 1:chunk_size:nframe
    stop_idx = min(start_idx + chunk_size - 1, nframe);
    chunk = movie_3d(:,:,start_idx:stop_idx);
    if size(chunk, 3) < min_chunk_frames
        continue;
    end

    chunk = smooth_chunk_spatially(chunk, spatial_sigma);
    avg_image = normalize_robust(mean(chunk, 3));
    active_image = normalize_robust(max(chunk, [], 3) - median(chunk, 3));
    chunk_score = avg_weight .* avg_image + active_weight .* active_image;

    aggregate_map = max(aggregate_map, single(chunk_score));
    active_map = max(active_map, single(active_image));
    chunk_counter = chunk_counter + 1;
end

if chunk_counter == 0
    chunk = smooth_chunk_spatially(movie_3d, spatial_sigma);
    avg_image = normalize_robust(mean(chunk, 3));
    active_image = normalize_robust(max(chunk, [], 3) - median(chunk, 3));
    aggregate_map = avg_weight .* avg_image + active_weight .* active_image;
    active_map = active_image;
end

% Favor pixels that are strong in both structure and transient activity.
quick_map = normalize_robust(0.7 .* aggregate_map + 0.3 .* active_map);
end

function chunk_smoothed = smooth_chunk_spatially(chunk, sigma)
chunk_smoothed = zeros(size(chunk), 'single');
for frame_idx = 1:size(chunk, 3)
    chunk_smoothed(:,:,frame_idx) = imgaussfilt(single(chunk(:,:,frame_idx)), sigma);
end
end

function image_norm = normalize_robust(image_in)
image_in = single(image_in);
low_val = prctile(image_in(:), 1);
high_val = prctile(image_in(:), 99.5);

if ~isfinite(low_val) || ~isfinite(high_val) || high_val <= low_val
    image_norm = zeros(size(image_in), 'single');
    return;
end

image_norm = (image_in - low_val) ./ (high_val - low_val);
image_norm = min(max(image_norm, 0), 1);
end
