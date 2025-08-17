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
quick_map = zeros(ncols/bin * nrows/bin,1);

% reshape to high dimension for each bin
movie_5D = reshape(movie,bin,ncols/bin,bin,nrows/bin,[]);
% average across x y
movie_ave = squeeze(mean(mean(movie_5D,1),3));

% reshape to binned
movie_binned = reshape(movie_ave,nrows/bin*ncols/bin,[]);
npixels = size(movie_binned,1);
% 
% % Apply avg_intensity curve to correct photobleaching from each pixel (recommend)
% 
% 
avg_trace = mean(movie_binned,1);
x = (1:nframe)';
% f = fit(x,avg_trace','poly1');
% fitted_curves = x*f.p1+f.p2;
% f = fit(x,avg_trace','exp2');
% fitted_curves = f.a*exp(f.b*x)+f.c*exp(f.d*x);
% movie_binned_corrected = movie_binned ./ fitted_curves';
% movie_binned_corrected = wdenoise(movie_binned_corrected',DenoisingMethod="FDR")';
k = 5;
n = 3;

print_count = 0;
% % 定义滤波器参数
fs = 400;                % 采样频率，单位为帧/秒（根据实际情况调整）
cutoff_freq_time = 1/10; % 截止频率，单位为Hz（根据需要调整）
order_time = 2;          % 滤波器阶数

% 设计巴特沃斯滤波器
[b_time, a_time] = butter(order_time, cutoff_freq_time/(fs/2), 'low');

% 设置填充长度，将信号在头尾各填充5%的长度
padlength = round(0.05 * size(movie_binned, 2));

% 对信号进行填充，以减少边界效应
padded_movie = [repmat(movie_binned(:, 1), 1,padlength), ...
                 movie_binned, ...
                 repmat(movie_binned(:,end), 1,padlength)];

% 应用高通滤波器到所有像素的时间序列
% filtfilt 会自动对每一列（对应一个像素的时间序列）进行滤波

movie_binned_corrected = filtfilt(b_time, a_time, padded_movie');
movie_binned_corrected = movie_binned_corrected';
movie_binned_corrected = movie_binned_corrected(:,padlength+1:end-padlength);
% movie_binned_corrected = detrend(movie_binned',4)';

%  movie_binned_corrected_3D =  reshape(movie_binned_corrected,ncols/bin,nrows/bin,[]);
% select_ROI(movie_binned_corrected, ncols/bin, nrows/bin, [], [])

for i = 1:npixels
    % fprintf('Processing %d / %d\n', i, npixels );
    fprintf(repmat('\b',1,print_count));   
    % fprintf('Calculating %.2f %% \n', i/npixels*100);
    print_count = fprintf('Calculating %.2f %% \n', i/npixels*100);
    % pixel_trace_corrected = movie_binned(i, :);
    pixel_trace_corrected = movie_binned_corrected(i, :);
    % pixel_trace_corrected = wdenoise(movie_binned_corrected(i, :));
    % pixel_trace_corrected = detrend(movie_binned_corrected(i, :),4);
    % pixel_trace_corrected =  movie_binned(i, :) - movie_binned_corrected (i, :);
    % baseline =mean(pixel_trace_corrected);
    baseline =mean(movie_binned(i, :)); % filted
    if strcmp(mode,'voltage')
        abstrace = abs(pixel_trace_corrected - baseline');
        maxpoint =  max(abstrace);
        maxpointindex = abstrace == maxpoint;
        pointdff = ((pixel_trace_corrected(maxpointindex)-mean(baseline))/mean(baseline));
        % quick_map(i) = pointdff;
        quick_map(i) = sign(pointdff(1)) .* (1 - exp(-k * abs(pointdff(1)*n).^n)) / (1 - exp(-k));


    elseif strcmp(mode,'calcium')
        % quick_map(i) = std(pixel_trace_corrected);
        quick_map(i) = (max(pixel_trace_corrected)-baseline) /std(pixel_trace_corrected);
    end
end

% converse to Z score
% baseline = mean(quick_map);
% quick_map = (quick_map - baseline) ./ std(quick_map);
quick_map = reshape(quick_map,ncols/bin,nrows/bin, 1);
quick_map = imresize(quick_map,[ncols,nrows] );
end

