%% Find Cell Mask
% find cell mask according to plot the min SNR(minus) of each pixel.

function [quick_map] = quick_find_cell(intensity_time_series, num_rows, num_cols, bin, Hz)

% Design high-pass filter
% dt = 1 / Hz;
% Fc = 1/n;
% [b, a] = butter(3, Fc*dt, 'high'); % 3rd order Butterworth high-pass filter

% cut the size
cut = false;
nrow = num_rows;
ncol = num_cols;
if mod(nrow,bin) ~= 0
%     intensity_time_series = reshape(intensity_time_series,num_rows,num_cols,[]);
    nrow = nrow - mod(nrow,bin);
    cut = true;
end
if mod(ncol,bin) ~= 0
%     intensity_time_series = reshape(intensity_time_series,num_rows,num_cols,[]);
    ncol = ncol - mod(ncol,bin);
    cut = true;
end
if cut == true
    temp = reshape(intensity_time_series,num_cols,num_rows,[]);
    intensity_time_series = temp(1:ncol,1:nrow,:);
end
its_size = size(intensity_time_series);
nframe = its_size(end);
% Initialize matrix to store cell labels
quick_map = zeros(ncol/bin * nrow/bin,1);

% change the Image to bin16
% reshape to high dimension for each bin
its5 = reshape(intensity_time_series,bin,ncol/bin,bin,nrow/bin,[]);
% average across x y
its_ave = squeeze(mean(mean(its5,1),3));

% reshape to binned
its_binned = reshape(its_ave,nrow/bin*ncol/bin,nframe);
% figure()
% imagesc(mean(reshape(intensity_time_series_binned,ncol/bin,nrow/bin,[]),3))
% Apply avg_intensity curve to correct photobleaching from each pixel
npixels = size(its_binned,1);

% % collect extremum after highpass filter
% for i = 1:npixels
%     fprintf('Processing %d / %d\n', i, npixels );
%     filted_pixel_trace = highpassfilter(intensity_time_series_merged(i, :)', Hz);
%     % baseline = abs(mean(intensity_time_series_merged(i, :)))
%     extremum = max(abs(filted_pixel_trace));
%     quick_map(i) =  extremum ;
% end


% Apply avg_intensity curve to correct photobleaching from each pixel (recommend)
avg_trace = mean(its_binned,1);
[fit_result, gof] = fit([1:nframe]' ,avg_trace', 'exp1' );
fit_trace = fit_result.a * exp(fit_result.b * [1: nframe]);

% plot([1:nframe],fit_trace)
% hold on
for i = 1:npixels
    fprintf('Processing %d / %d\n', i, npixels );
    pixel_trace_correct = its_binned(i, :) ./ fit_trace;
%     plot([1:nframe],pixel_trace_correct)
%     hold on
%   pixel_trace_correct = intensity_time_series_merged(i, :) ./ avg_trace;
%     % High-pass filter
%     pixel_trace_correct = filtfilt(b, a, intensity_time_series_merged(i, :));
    baseline = mean(pixel_trace_correct);
    quick_map(i) = (min(pixel_trace_correct)- baseline) ./ baseline;
%     plot(i,quick_map(i),'o')
%     hold on
end

% converse to Z score
baseline = mean(quick_map);
quick_map = (quick_map - baseline) ./ std(quick_map);

quick_map = reshape(quick_map,ncol/bin,nrow/bin, 1);
% quick_map = permute(quick_map, [2, 1, 3]);
quick_map = imresize(quick_map,[ncol,nrow] );
end

% 
% figure()
% subplot(3,1,1)
% plot(pixel_trace_correct);
% title('pixel_trace_correct');
% subplot(3,1,2)
% plot(intensity_time_series(i, :)) ;
% title('intensity_time_series(i, :) ');
% subplot(3,1,3)
% plot(avg_trace) ;
% title('avg_trace');

