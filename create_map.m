function [quick_map] = create_map(movie, num_rows, num_cols, bin, mode)

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
nrow = num_rows;
ncol = num_cols;
if mod(nrow,bin) ~= 0
    nrow = nrow - mod(nrow,bin);
    cut = true;
end
if mod(ncol,bin) ~= 0
    ncol = ncol - mod(ncol,bin);
    cut = true;
end
if cut == true
    temp = reshape(movie,num_cols,num_rows,[]);
    movie = temp(1:ncol,1:nrow,:);
end
movie_size = size(movie);
nframe = movie_size(end);

% Initialize matrix to store cell labels
quick_map = zeros(ncol/bin * nrow/bin,1);

% reshape to high dimension for each bin
movie_5D = reshape(movie,bin,ncol/bin,bin,nrow/bin,[]);
% average across x y
movie_ave = squeeze(mean(mean(movie_5D,1),3));

% reshape to binned
movie_binned = reshape(movie_ave,nrow/bin*ncol/bin,nframe);
npixels = size(movie_binned,1);

% Apply avg_intensity curve to correct photobleaching from each pixel (recommend)
avg_trace = mean(movie_binned,1);
[fit_result, gof] = fit([1:nframe]' ,avg_trace', 'exp1' );
fit_trace = fit_result.a * exp(fit_result.b * [1: nframe]);

for i = 1:npixels
    fprintf('Processing %d / %d\n', i, npixels );
    pixel_trace_correct = movie_binned(i, :) ./ fit_trace;
    baseline = mean(pixel_trace_correct);
    if strcmp(mode,'voltage')
        quick_map(i) = (max(abs(pixel_trace_correct))- baseline) ./ baseline; % do not judge polarity
    elseif strcmp(mode,'calcium')
        quick_map(i) = std(pixel_trace_correct);
    end
end

% converse to Z score
baseline = mean(quick_map);
quick_map = (quick_map - baseline) ./ std(quick_map);
quick_map = reshape(quick_map,ncol/bin,nrow/bin, 1);
quick_map = imresize(quick_map,[ncol,nrow] );
end

