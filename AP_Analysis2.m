% AP_ANALYSIS_Auto - Quick Analysis of Voltage Signals
% This script is used for quickly analyzing voltage signals. It requires
% several functions and toolboxes to run properly.
%
% Requirements:
%   Functions: calculate_firing_rate, calculate_FWHM, create_map, calculate_SNR,
%              fit_exp1, highpassfilter, select_ROI
%   Toolboxes: Image Processing Toolbox, Curve Fitting Toolbox, Signal Processing Toolbox
%
% Usage:
%   1. Set the folder_path and file variables to point to your data.
%   2. Specify the frame rate (freq) of your data.
%   3. Run the script.
%
% Example:
%   folder_path = 'F:\20240321\';
%   file = '20240321.tif'; % Include the format in the file name.
%   freq = 400; % Frequency in Hz
%
% Author: Liu-Yang Luorong

% Version: 5.0
% Date: 2025.04.20
% GitHub: https://github.com/Luorong37/AP_analysis
%
% See also calculate_firing_rate, calculate_FWHM, create_map, calculate_SNR, fit_exp1, highpassfilter, select_ROI


clear; clc;
%## 写一个画图的段落
%% Loading raw data
nowtime = string(datetime( 'now'));
% Replace colons with hyphens to get the desired output format
nowtime = strrep(nowtime , ':', '-');
fprintf('Loading...\n')

% ↓↓↓↓↓-----------Prompt user for define path-----------↓↓↓↓↓
% support for folder, .tif, .tiff, .bin.
folder_path = 'E:\1_Data\Luorong\25.10.17_invivo\#77\';
file = '\R1_8h';  % must add format.do not add '\' at last
% ↓↓↓↓↓-----------Prompt user for frame rate------------↓↓↓↓↓
freq = 400; % Hz
gpu = true; % defined gpu open
% -----------------------------------------------------------

% Create an analysis folder
function [file_path, save_path] = create_folder(folder_path, file, nowtime)
% read path
file_path = fullfile(folder_path, file);
if isfolder(file_path)
    file_name = file;
    %file_dir = dir(file_path);
    %[~, ~, file_extension] = fileparts(file_dir(3).name);
else
    % [~, file_name, file_extension] = fileparts(file_path)
    [~, file_name, ~] = fileparts(file_path);
end

% create a folder for analysis
save_path = fullfile(folder_path, strcat(file_name, '_Analysis'), nowtime);
mkdir(save_path);
end

[file_path, save_path] = create_folder(folder_path, file, nowtime);

% Load image file
if gpu
    gcp;
end
[movie, ncols, nrows, nframes] = load_movie(file_path);

% Presetting
function [dt, colors, t, map, mask, avg_image, options] = presetting(freq, nframes, movie, ncols, nrows)

% Define parameters
dt = 1 / freq; % Calculate time axis
% colors = [lines(7);hsv(5);spring(3);winter(3);gray(3)];
colors = lines(100);
t = (1:nframes) * dt;
options.colors = colors;
map = [];
mask = [];
movie_vol_2D = reshape(mean(movie,2), ncols, nrows, []);
avg_image  = (movie_vol_2D - min(movie_vol_2D(:))) / (max(movie_vol_2D(:)) - min(movie_vol_2D(:)));

end

[dt, colors, t, map, mask, options]= presetting(freq, nframes, movie, ncols, nrows);
x = (1:nframes)' * dt;

% Save code
code_path = fullfile(save_path,'Code');
mkdir(code_path);
currentScript = which("AP_Analysis2.m");
% 获取当前脚本依赖的所有文件
[requiredFiles, ~] = matlab.codetools.requiredFilesAndProducts(currentScript);
% 复制当前脚本和所有依赖文件到目标文件夹
for k = 1:length(requiredFiles)
    [~, name, ext] = fileparts(requiredFiles{k});
    copyfile(requiredFiles{k}, fullfile(code_path, [name, ext]));
end
% 提示完成
fprintf('All codes have been copied to %s\n', code_path);
%% motion correction
% this section is drived from demo_1p_low_RAM of NoRMCorre
fprintf('Motion correcting...\n');
t1 = tic;
% [M1f,shifts1] = motion_correction(reshape(movie, ncols, nrows, []),'file_path',file_path);
d1 = ncols;
d2 = nrows;

h5_list = dir(fullfile(file_path, '*.h5'));
if isempty(h5_list)
    presave = true;
else
    presave = false;
end

bin_width = 200;
max_shift = [20 20];
grid_size = [128,128];
overlap_pre = [32,32];
mot_uf = 4;
iter = 1;
correct_bidir = false;
% appropriate name through options_r.tiff_filename or options_r.h5_filename


if isfolder(file_path)
    file_list = dir(fullfile(file_path, '*.tif'));
    file_nums = length(file_list);
    movie_mc = [];
    
    if file_nums ~= length(h5_list)
         presave = logical([zeros(length(h5_list),1),ones(file_nums-length(h5_list),1)]);
    end

    for i = 1:file_nums
        file = fullfile(file_path,file_list(i).name);
        if presave(i)
            [h5_name] = create_highpass_temple(file);
        else
            h5_name = fullfile(file_path,h5_list(i).name);
        end
        % [~, filename, ~] = fileparts(['\' file_list(i).name]);
        options_r = NoRMCorreSetParms('d1',d1,'d2',d2,'bin_width',bin_width,'max_shift',max_shift,'grid_size',grid_size,'overlap_pre', overlap_pre,'mot_uf' , mot_uf,'iter',iter,'correct_bidir',correct_bidir);

        %register using the high pass filtered data and apply shifts to original data
        tic; [~,shifts1,~] = normcorre_batch(h5_name,options_r); toc % register filtered data
        % exclude boundaries due to high pass filtering effects

        % if you save the file directly in memory make sure you save it with a
        % name that does not exist. Change options_r.tiff_filename
        % or options_r.h5_filename accordingly.
        %
        tic; Mr = apply_shifts(file,shifts1,options_r); toc % apply shifts to full dataset

        mc_filename = fullfile(save_path, ['0_motion_correction_',file_list(i).name]);
        save_tiffs(Mr,mc_filename)
        %save(mc_filename,'Mr','d1','d2','bin_width','max_shift','grid_size' ,'overlap_pre','mot_uf','iter','correct_bidir');

        movie_mc = cat(3,movie_mc,Mr);
    end
    movie_vol_2D = mean(movie_mc,3);
    avg_image  = (movie_vol_2D - min(movie_vol_2D(:))) / (max(movie_vol_2D(:)) - min(movie_vol_2D(:)));
    movie = reshape(movie_mc, ncols*nrows, []);
else
        if presave
            [h5_name] = create_highpass_temple(file_path);
        else
            h5_name = fullfile(file_path,h5_list.name);
        end
    [~, filename, ~] = fileparts(['\' file_path]);

    options_r = NoRMCorreSetParms('d1',d1,'d2',d2,'bin_width',bin_width,'max_shift',max_shift,'grid_size',grid_size, ...
        'overlap_pre', overlap_pre,'mot_uf' , mot_uf,'iter',iter,'correct_bidir',correct_bidir);

    % options_r = NoRMCorreSetParms('d1',d1,'d2',d2,'bin_width',200,'max_shift',[20 20],'grid_size' ,[128,128],'overlap_pre', [32,32],'mot_uf' , 4,'iter',1,'correct_bidir',false);

    %register using the high pass filtered data and apply shifts to original data
    tic; [~,shifts1,~] = normcorre_batch(h5_name,options_r); toc % register filtered data
    % exclude boundaries due to high pass filtering effects

    % if you save the file directly in memory make sure you save it with a
    % name that does not exist. Change options_r.tiff_filename
    % or options_r.h5_filename accordingly.
    %
    tic; Mr = apply_shifts(file_path,shifts1,options_r); toc % apply shifts to full dataset
    mc_filename = fullfile(save_path, ['0_motion_correction_',filename]);
    save_tiffs(Mr,mc_filename)
    %save(mc_filename,'Mr','d1','d2','bin_width','max_shift','grid_size' ,'overlap_pre','mot_uf','iter','correct_bidir');

    movie = reshape(Mr, ncols*nrows, []);
    movie_vol_2D = mean(Mr,3);
    avg_image  = (movie_vol_2D - min(movie_vol_2D(:))) / (max(movie_vol_2D(:)) - min(movie_vol_2D(:)));
end
%first try out rigid motion correction
%     % exclude boundaries due to high pass filtering effects
% options_r = NoRMCorreSetParms('d1',d1,'d2',d2,'bin_width',200,'max_shift',[20 20],'grid_size' ,[128,128],'overlap_pre', [32,32],'mot_uf' , 4,'iter',1,'correct_bidir',false);
%
% %register using the high pass filtered data and apply shifts to original data
% tic; [M1,shifts1,template1] = normcorre_batch(h5_name,options_r); toc % register filtered data
%     % exclude boundaries due to high pass filtering effects
%
% % if you save the file directly in memory make sure you save it with a
% % name that does not exist. Change options_r.tiff_filename
% % or options_r.h5_filename accordingly.
% %
% tic; Mr = apply_shifts(file_path,shifts1,options_r); toc % apply shifts to full dataset
%
% movie = reshape(Mr, ncols*nrows, []);
% movie_vol_2D = mean(Mr,3);
% avg_image  = (movie_vol_2D - min(movie_vol_2D(:))) / (max(movie_vol_2D(:)) - min(movie_vol_2D(:)));

t2 = toc(t1);
fprintf('Finished motion correction after %d s\n',round(t2))

function [h5_name] = create_highpass_temple(file_path)
gSig = 7; 
gSiz = 3*gSig; 
psf = fspecial('gaussian', round(2*gSiz), gSig);
ind_nonzero = (psf(:)>=max(psf(:,1)));
psf = psf-mean(psf(ind_nonzero));
psf(~ind_nonzero) = 0;   % only use pixels within the center disk
[filepath,file_name,~] = fileparts(file_path);
h5_name = fullfile(filepath,[file_name,'_filtered_data.h5']);
chunksize = 1000;    % read 500 frames at a time

cnt = 1;
while (1)  % read filter and save file in chunks
    Yf = single(read_file(file_path,cnt,chunksize));
    if isempty(Yf)
        break
    else
        Y = imfilter(Yf,psf,'symmetric');
        saveash5(Y,h5_name);
        cnt = cnt + size(Y,ndims(Y));
    end
    disp(cnt)
end
end

% save_tiffs(Mr,mc_filename)
function save_tiffs(movie,tiff_file_name)
fprintf('Saving results as tiff...\n')
tic;
t = Tiff(tiff_file_name, 'w');
batch_start = 1;
[ncols,nrows,batch_end] = size(movie);
    for j = batch_start:batch_end
        current_image = movie(:,:,j);         % 读取当前图片

        t.setTag('ImageLength', nrows);
        t.setTag('ImageWidth', ncols);
        t.setTag('Photometric', Tiff.Photometric.MinIsBlack);
        t.setTag('BitsPerSample', 16);
        t.setTag('SamplesPerPixel', 1);
        t.setTag('RowsPerStrip', 16);
        t.setTag('PlanarConfiguration', Tiff.PlanarConfiguration.Chunky);
        t.setTag('Compression', Tiff.Compression.None);
        t.setTag('Software', 'MATLAB');
        t.write(current_image);

        if j < batch_end
            t.writeDirectory();
        end

    end
    t.close();
 toc;
end

%% Create a map (optional)
t1 = tic; % Start a timer
fprintf('Creating a map...\n')
% if the map cannot figure out active cells, please large the bin.
bin = 8; % defined bin = 4

[quick_map] = create_map(movie, nrows, ncols, bin);
map = quick_map;

function plot_map(quick_map, save_path, map)
% Visualize correlation coefficients as heatmap
figure()
imagesc(quick_map);
title('Sensitivity Map');
axis image;
colorbar;
fig_filename = fullfile(save_path, '0_Sensitivity_Map.fig');
png_filename = fullfile(save_path, '0_Sensitivity_Map.png');
mat_filename = fullfile(save_path, '0_Sensitivity_Map.mat');

saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');
save(mat_filename, 'map');
end

plot_map(quick_map, save_path, map);

t2 = toc(t1); % Get the elapsed time
fprintf('Finished mask creating after %d s\n',round(t2))
%% Load Mask (optinal)
methods = questdlg('Load previous data?','load data','previous ROIs','cellpose','Cancel');
switch methods
    case 'previous ROIs'
        [roi_filename,roi_foldername] = uigetfile(save_path);
        rois_data = load(fullfile(roi_foldername,roi_filename));
    try
        rois = rois_data.rois;
        mask = rois.bwmask;
    catch ME
        rois.bwmask = bwmask;
        mask = rois.bwmask;
    end
    case 'cellpose'
    fprintf('cellpose running......\n');
    cp = cellpose(ExecutionEnvironment="gpu");
    avgdia = 15;
    gamma_image = imadjust(avg_image,[],[],1.2); % recommend raise the gamma factor from 1 to 4.
    mask = segmentCells2D(cp, gamma_image , ImageCellDiameter = avgdia,FlowErrorThreshold = 3,CellThreshold = -3);% CellThreshold = -2, ,  FlowErrorThreshold = 2
    fprintf('%d cells are found.\n',max(mask(:)));
    figure()
    % 使用labeloverlay函数显示图像
    overlayImage = labeloverlay(gamma_image, mask,'Transparency', 0.6);
    % 显示结果
    imshow(overlayImage);
    title('cellpose mask');
    fig_filename = fullfile(save_path, '0_cellpose_mask.fig');
    png_filename = fullfile(save_path, '0_cellpose_mask.png');
    mat_filename = fullfile(save_path, '0_cellpose_mask.mat');
    
    saveas(gcf, fig_filename, 'fig');
    saveas(gcf, png_filename, 'png');
    save(mat_filename, 'mask');
end
%% Select ROI
t1 = tic; % Start a timer
figure()
% with or without Mask and Map
if exist('methods','var')
    if any(strcmp(methods,{'No',''} ))
        [rois, traces] = select_ROI(movie, ncols, nrows, mask, map);
        nrois = max(rois.bwmask,[],'all');
    else
        [~, traces] = select_ROI(movie, ncols, nrows, mask, map);
        nois = max(mask(:));
        rois.bwmask = mask;
    end
else
    [rois, traces] = select_ROI(movie, ncols, nrows, mask, map);
    nrois = max(rois.bwmask,[],'all');
end

try
bwmask = rois.bwmask;
catch ME

end
traces_original = traces;

% if do not need a map, run the following code:
% map = [];
% [bwmask, traces] = select_ROI(movie, nrows, ncols, t, mask, map);

fig_filename = fullfile(save_path, '1_raw_trace.fig');
png_filename = fullfile(save_path, '1_raw_trace.png');
roi_filename = fullfile(save_path, '1_raw_ROI.mat');

saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');
save(roi_filename, 'rois','avg_image','traces');

h = msgbox('All rois saved.', 'Done', 'help');
uiwait(h);
%close(gcf);
%% Signal Process %%

% Background Correction
fprintf('Removing Background...\n')
[background, background_fitted, traces_bgcorr, traces_bgfitcorr, background_mask]...
    = remove_background(movie, ncols, nrows, rois);
roi_filename = fullfile(save_path, '1_background_ROI.mat');
save(roi_filename, 'background','background_fitted','background_mask','traces_bgcorr','traces_bgfitcorr');
fprintf('Finished\n')

%% Bleaching Correction

function plot_corrected(traces_corrected, traces, baseline, colors, save_path)
% plot
fig = figure();
set(fig,'Position',get(0,'Screensize'));
fit_axe = subplot(3,1,1);
fited_axe = subplot(3,1,2);
baseline_axe = subplot(3,1,3);

% plot
for i = 1: size(traces_corrected,2)
    plot(traces(:,i),'Color',colors(i,:),'Parent',fit_axe);
    hold(fit_axe, 'on');
    plot(traces_corrected(:,i),'Color',colors(i,:),'Parent',fited_axe);
    hold(fited_axe, 'on');
    plot( baseline(:,i),'Color',colors(i,:),'LineWidth',2,'Parent',baseline_axe);
    hold(baseline_axe, 'on');

end
hold off;

% note
title(fit_axe, 'Original and Fitted Curves');
title(fited_axe, 'Corrected Traces');
title(baseline_axe, 'Baseline');
legend(fit_axe, 'Original Trace', 'Fitted Curve');

fig_filename = fullfile(save_path, '2_fitted_trace.fig');
png_filename = fullfile(save_path, '2_fitted_trace.png');

saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');
end

fprintf('Correcting Bleaching...\n')

% highpass bleach
% fc = 1/200000;
fc = 1/t(end);
[traces_corrected,baseline] = highpass_bleach_remove(traces_bgfitcorr,freq,fc );
% 1/t for remove fake hyperpolarization

plot_corrected(traces_corrected, traces, baseline, colors, save_path);
fprintf('Finished\n');

%% Wavelet process

% wavelet降噪
Dnmethods = 'FDR';
Dnlevel = 8;
fprintf('Wavelet Denoising...\n')
traces_denoised = wdenoise(traces_corrected, Dnlevel ,DenoisingMethod=Dnmethods);

traces_filename = fullfile(save_path, '2_processed_traces.mat');
save(traces_filename,"traces_denoised", "traces_original", 'traces_corrected', 'baseline', 'Dnmethods','Dnlevel')
fprintf(' Finished\n');
%% Calculate sensitivity and SNR

traces_sensitivity = traces_corrected./baseline;

noise = traces_corrected-traces_denoised;
traces_SNR = traces_corrected./std(noise);

traces_filename = fullfile(save_path, '3_calculated_traces.mat');
save(traces_filename,"traces_sensitivity", 'traces_corrected', 'baseline', 'noise','traces_SNR')
figure()
%plot raw
subplot(1,3,1);
title('raw');
hold on;
[~] = offset_plot(traces,t);


% plot sensitivity
subplot(1,3,2);
title('Sensitivity');
hold on;
[~] = offset_plot(traces_sensitivity,t);

% plot SNR
subplot(1,3,3);
title('SNR');
hold on;
[~] = offset_plot(traces_SNR,t);

fig_filename = fullfile(save_path, '4_SNR.fig');
png_filename = fullfile(save_path, '4_SNR.png');
trace_filename = fullfile(save_path, '4_SNR.mat');

saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');
%% AP Processing %%

parts = 1;
MinPeakProminence_factor = 0.50;
MinPeakDistance_factor = 5;

findmode = 'dn';% find via denoised traces
%findmode = 'cr';% find via bleach corrected traces

% Peak finding
switch findmode
    case 'cr'
    [peaks_index, peaks_amplitude, peaks_polarity, parts_results] = peak_finding_auto(traces_corrected , save_path,'parts',parts, ...
    'MinPeakProminence_factor',MinPeakProminence_factor,'MinPeakHeight', 0, 'MinPeakDistance_factor' , MinPeakDistance_factor);
    case 'dn'
        [peaks_index, peaks_amplitude, peaks_polarity, parts_results] = peak_finding_auto(traces_denoised , save_path,'parts',parts, ...
    'MinPeakProminence_factor',MinPeakProminence_factor,'RawTraces',traces_corrected,'MinPeakHeight', 0, 'MinPeakDistance_factor' , MinPeakDistance_factor);
end
%% manually reselction
MinPeakProminence_factor = 0.3;
pr = 1;
rr = 19;
allr = true; % manually select all trace
if allr
    r0 = 1;
    p0 = 1;
else
    r0 = rr;
    p0 = pr;
end
for pr = p0:pr
    for rr = r0:rr
        switch findmode
             case 'cr'
                 current_traces = traces_corrected(parts_results.index{pr},rr);
                 case 'dn'
current_traces = traces_denoised(parts_results.index{pr},rr);
        end

[peaks_polarity_re, ~, peaks_index_re, peaks_amplitude_re, ~] = ...
                                            peak_finding(current_traces,MinPeakProminence_factor,save_path);
% saveas(gcf,fullfile(save_path,sprintf('repeakfinding of roi %d, p = %d.png',roi_re,part_re)));
% saveas(gcf,fullfile(save_path,sprintf('repeakfinding of roi %d, p = %d.fig',roi_re,part_re)));

parts_results.peaks_amplitude(pr,rr) = peaks_amplitude_re;
parts_results.peaks_index(pr,rr) = {peaks_index_re{1} + parts_results.index{pr}(1)-1};
parts_results.peaks_polarity(pr,rr) = peaks_polarity_re;

peaks_index= cell(1,nrois);
for i = 1:parts
    for j = 1:nrois
        peaks_index{j} = [peaks_index{j} ;parts_results.peaks_index{i,j}];
    end
end

peaks_amplitude = cell(1,nrois);
for i = 1:parts
    for j = 1:nrois
        peaks_amplitude{j} = [peaks_amplitude{j} ;parts_results.peaks_amplitude{i,j}];
    end
end
% 
% peaks_sensitivity= cell(1,nrois);
% for i = 1:parts
%     for j = 1:nrois
%         peaks_sensitivity{j} = [peaks_sensitivity{j} ;peaks_amplitude_part{i,j}];
%     end
% end

peaks_polarity = cell(1,nrois);
for i = 1:nrois
    polarity_index = find(abs(parts_results.peaks_polarity(:,i)) == max(abs(parts_results.peaks_polarity(:,i))));
    peaks_polarity{i} =parts_results.peaks_polarity(polarity_index(1),i);
end
    end
end


%% Statistic AP, FWHM gated
AP_window_width = 15; % number of frames to for AP window (defined = 40)
offset_width = 3;
% AP_list = AP_statistic(nrois, peaks_index, peaks_amplitude, traces_corrected, traces_sensitivity, traces_SNR, AP_window_width, nframes, dt, peaks_polarity, save_path);

% function [AP_list,peaks_index_corrected]  = AP_statistic(nrois, peaks_index, peaks_amplitude, traces_corrected, traces_sensitivity, traces_SNR, AP_window_width, nframes, dt, peaks_polarity, save_path)
AP_list = cell(1, nrois);

% each trace
for i = 1:nrois % i for trace
    peaks_num = length(peaks_index{i});
    each_trace_amp = traces_corrected(:,i);
    each_trace_sensitivity = traces_sensitivity(:,i);
    each_trace_SNR = traces_SNR(:,i);
    AP_list{i} = cell(1, length(peaks_index{i}));
    fprintf('ROI %d processing\n',i);
    j = 1;
    % each peak
    while j <= peaks_num % j for peak

        peak_index_ij = peaks_index{i}(j);
        peak_amp_ij = peaks_amplitude{i}(j);

        % keep in board
        AP_start_index = max(1, peak_index_ij - AP_window_width);
        AP_end_index = min(nframes, peak_index_ij + AP_window_width);
        AP_index = AP_start_index : AP_end_index;

        % search
        AP_amp = each_trace_amp(AP_start_index:AP_end_index)';
        AP_sensitivity = each_trace_sensitivity(AP_start_index:AP_end_index)';
        AP_SNR = each_trace_SNR(AP_start_index:AP_end_index)';

        % fill NaN
        if 1 > peak_index_ij - AP_window_width
            AP_amp = [NaN(1,0 - (peak_index_ij - AP_window_width)+1), AP_amp];
            AP_sensitivity = [NaN(1,0 - (peak_index_ij - AP_window_width)+1),AP_sensitivity];
            AP_SNR = [NaN(1,0 - (peak_index_ij - AP_window_width)+1),AP_SNR];
        elseif nframes < peak_index_ij + AP_window_width
            AP_amp = [AP_amp, NaN(1,peak_index_ij + AP_window_width - nframes)];
            AP_sensitivity = [AP_sensitivity, NaN(1,peak_index_ij + AP_window_width - nframes)];
            AP_SNR = [AP_SNR, NaN(1,peak_index_ij + AP_window_width - nframes)];
        end

        % Calculate;
        Amplitude = abs(peak_amp_ij);
        Sensitivity = AP_sensitivity(AP_window_width+1)*100 ;
        SNR = abs(AP_SNR(AP_window_width+1));
        f = false; % save each peak
        offset = peak_offset(AP_amp(AP_window_width-offset_width + 1:AP_window_width + offset_width+ 1), peaks_polarity{i});
        FWHM = calculate_FWHM(AP_amp, dt,  peaks_polarity{i},f);
        % FWHM = calculate_FWHM2(AP_amp, dt, peaks_polarity{i});
        if isempty(offset)
            offset = 99;
        end
        if all(abs(offset) > 2) ||FWHM <= 2.5
            FWHM = NaN;
        end

               
        if isnan(FWHM)
            FWHM = calculate_FWHM(AP_amp, dt,  peaks_polarity{i},true);
            title(sprintf('failed peak at noi %d, peaks %d',i,j))
            if ~isfolder(fullfile(save_path,'failed peaks',sprintf('roi %d',i)))
            mkdir(fullfile(save_path,'failed peaks',sprintf('roi %d',i)));
            end
            saveas(gcf,fullfile(save_path,'failed peaks',sprintf('roi %d',i),sprintf('false peak at noi %d, peaks %d.png',i,j)))
            close(gcf)
        elseif f
            title(sprintf('finded peak at noi %d, peaks %d',i,j))
            if ~isfolder(fullfile(save_path,'finded peaks',sprintf('roi %d',i)))
            mkdir(fullfile(save_path,'finded peaks',sprintf('roi %d',i)));
            end
            saveas(gcf,fullfile(save_path,'finded peaks',sprintf('roi %d',i),sprintf('finded  peak at noi %d, peaks %d.png',i,j)))
            close(gcf)
        end
        % sprintf('roi % d peaks %d FWHM:%d',i,j,FWHM);

        % save AP data
        each_AP = struct('Trace', i, 'AP_number', j, 'AP_index',AP_index, ...
            'AP_amp',AP_amp,'Amplitude', Amplitude,'FWHM',FWHM, ...
            'AP_sensitivity',AP_sensitivity,'Sensitivity',Sensitivity, ...
            'AP_SNR', AP_SNR, 'SNR', SNR);
        AP_list{i}{j} = each_AP;
        if offset ~= 0
            peaks_index{i}(j) = peaks_index{i}(j) + offset;
            fprintf('peaksindex %d in roi %d redirection\n',peaks_index{i}(j),i)
            j = j -1;
        end
        j = j + 1;
    end

    
end

save(fullfile(save_path,'FWHM gated peaks.mat'),'peaks_index','peaks_polarity','peaks_amplitude')
peaks_index_manually_gated = [];
% AP_window_width = 15; % number of frames to for AP window (defined = 40)
% offset_width = MinPeakDistance_factor;
% % AP_list = AP_statistic(nrois, peaks_index, peaks_amplitude, traces_corrected, traces_sensitivity, traces_SNR, AP_window_width, nframes, dt, peaks_polarity, save_path);
% 
% % function [AP_list,peaks_index_corrected]  = AP_statistic(nrois, peaks_index, peaks_amplitude, traces_corrected, traces_sensitivity, traces_SNR, AP_window_width, nframes, dt, peaks_polarity, save_path)
% AP_list = cell(1, nrois);
% 
% % each trace
% for i = 1:nrois % i for trace
%     peaks_num = length(peaks_index{i});
%     each_trace_amp = traces_corrected(:,i);
%     each_trace_sensitivity = traces_sensitivity(:,i);
%     each_trace_SNR = traces_SNR(:,i);
%     AP_list{i} = cell(1, length(peaks_index{i}));
% 
%     j = 1;
%     % each peak
%     while j <= peaks_num % j for peak
% 
%         peak_index_ij = peaks_index{i}(j);
%         peak_amp_ij = peaks_amplitude{i}(j);
% 
%         % keep in board
%         AP_start_index = max(1, peak_index_ij - AP_window_width);
%         AP_end_index = min(nframes, peak_index_ij + AP_window_width);
%         AP_index = AP_start_index : AP_end_index;
% 
%         % search
%         AP_amp = each_trace_amp(AP_start_index:AP_end_index)';
%         AP_sensitivity = each_trace_sensitivity(AP_start_index:AP_end_index)';
%         AP_SNR = each_trace_SNR(AP_start_index:AP_end_index)';
% 
%         % fill NaN
%         if 1 > peak_index_ij - AP_window_width
%             AP_amp = [NaN(1,0 - (peak_index_ij - AP_window_width)+1), AP_amp];
%             AP_sensitivity = [NaN(1,0 - (peak_index_ij - AP_window_width)+1),AP_sensitivity];
%             AP_SNR = [NaN(1,0 - (peak_index_ij - AP_window_width)+1),AP_SNR];
%         elseif nframes < peak_index_ij + AP_window_width
%             AP_amp = [AP_amp, NaN(1,peak_index_ij + AP_window_width - nframes)];
%             AP_sensitivity = [AP_sensitivity, NaN(1,peak_index_ij + AP_window_width - nframes)];
%             AP_SNR = [AP_SNR, NaN(1,peak_index_ij + AP_window_width - nframes)];
%         end
% 
%         % Calculate;
%         Amplitude = abs(peak_amp_ij);
%         Sensitivity = AP_sensitivity(AP_window_width+1)*100 ;
%         SNR = abs(AP_SNR(AP_window_width+1));
%         f = false; % save each peak
%         AP_amp_cut = AP_amp(AP_window_width-offset_width + 1:AP_window_width + offset_width+ 1);
%         if AP_amp(AP_window_width + 1)*peaks_polarity{i} < AP_amp(AP_window_width) *peaks_polarity{i}...
%             || AP_amp(AP_window_width + 1)*peaks_polarity{i} < AP_amp(AP_window_width+2)*peaks_polarity{i}
%             offset = peak_offset(AP_amp_cut, peaks_polarity{i});
%             FWHM = calculate_FWHM(AP_amp, dt,  peaks_polarity{i},f);
%         % FWHM = calculate_FWHM2(AP_amp, dt, peaks_polarity{i});
%         else
%             offset = 99;
%             FWHM = NaN;
%         end
% 
%         if isempty(offset) || abs(offset) > 2 || FWHM <= 2.5
%             FWHM = NaN;
%         end
% 
%         if isnan(FWHM)
%             try
%             FWHM = calculate_FWHM(AP_amp, dt,  peaks_polarity{i},true);
%             catch
%                 figure()
%                 plot(AP_amp); hold on;
%             end
%             title(sprintf('failed peak at noi %d, peaks %d',i,j))
%             if ~isfolder(fullfile(save_path,'failed peaks',sprintf('roi %d',i)))
%             mkdir(fullfile(save_path,'failed peaks',sprintf('roi %d',i)));
%             end
%             saveas(gcf,fullfile(save_path,'failed peaks',sprintf('roi %d',i),sprintf('false peak at noi %d, peaks %d.png',i,j)))
%             close(gcf)
%         elseif f
%             title(sprintf('finded peak at noi %d, peaks %d',i,j))
%             if ~isfolder(fullfile(save_path,'finded peaks',sprintf('roi %d',i)))
%             mkdir(fullfile(save_path,'finded peaks',sprintf('roi %d',i)));
%             end
%             saveas(gcf,fullfile(save_path,'finded peaks',sprintf('roi %d',i),sprintf('finded  peak at noi %d, peaks %d.png',i,j)))
%             close(gcf)
%         end
%         % sprintf('roi % d peaks %d FWHM:%d',i,j,FWHM);
% 
%         % save AP data
%         each_AP = struct('Trace', i, 'AP_number', j, 'AP_index',AP_index, ...
%             'AP_amp',AP_amp,'Amplitude', Amplitude,'FWHM',FWHM, ...
%             'AP_sensitivity',AP_sensitivity,'Sensitivity',Sensitivity, ...
%             'AP_SNR', AP_SNR, 'SNR', SNR);
%         AP_list{i}{j} = each_AP;
%         if ~isempty(offset)  && offset ~= 0  
%             peaks_index{i}(j) = peaks_index{i}(j) + offset;
%             fprintf('peaksindex %d in roi %d redirection\n',peaks_index{i}(j),i)
%             j = j -1;
%         end
%         j = j + 1;
%     end
% 
% 
% end
% % end
% save(fullfile(save_path,'FWHM gated peaks.mat'),'peaks_index','peaks_polarity','peaks_amplitude')
% peaks_index_manually_gated = [];
%% manually gate (optional)

peaks_index_manually_gated = peaks_index;

for i = 1:nrois
    figure()
    set(gcf,'Position',[0,0,2500,1800])
    
    hold on;
    plot(traces_SNR(:,i).*peaks_polarity{i});
    % peaks_x = peaks_index{i}(~isnan(AP_data.FWHM{i})); 
    peaks_x = peaks_index{i}; 
    peaks_y = traces_SNR(peaks_x,i).*peaks_polarity{i};
    plot(peaks_x, peaks_y,'v','MarkerFaceColor','r');
    title(sprintf('noi %d, %d peaks. DELETE peaks in rectangle',i,length(peaks_x)));
    % mkdir(fullfile(save_path,'FWHM gated peaks'));
    % saveas(gcf,fullfile(save_path,'FWHM gated peaks',sprintf('noi %d, peaks %d.png',i,length(peaks_x))))
    % saveas(gcf,fullfile(save_path,'FWHM gated peaks',sprintf('noi %d, peaks %d.fig',i,length(peaks_x))))

    rect = drawrectangle(gca);
    manual_gated_index = rect.Position(1) < peaks_x & peaks_x < rect.Position(1) + rect.Position(3) & ...
    rect.Position(2) < peaks_y & peaks_y < rect.Position(2) + rect.Position(4);
    
    fprintf(' %d peaks deleted.\n', sum(manual_gated_index));

    peaks_index_manually_gated{i}(manual_gated_index) = NaN;
    peaks_index_manually_gated{i}(~manual_gated_index) = 1;


    plot(peaks_x(manual_gated_index), peaks_y(manual_gated_index),'v','MarkerFaceColor','g');

    mkdir(fullfile(save_path,'Manually gated peaks'));
    saveas(gcf,fullfile(save_path,'Manually gated peaks',sprintf('noi %d, peaks %d.png',i,length(peaks_x))))
    saveas(gcf,fullfile(save_path,'Manually gated peaks',sprintf('noi %d, peaks %d.fig',i,length(peaks_x))))
    close(gcf)
end
save(fullfile(save_path,'Manually gated peaks.mat'),'peaks_index_manually_gated')
%% AP data statistic
% AP_data = AP_save(AP_list, peaks_index, nrois,save_path);


% function AP_data = AP_save(AP_list, peaks_index, nrois, save_path)
% write into excel
% 初始化平均值向量
% avg_amp = zeros(length(AP_list), 1);
avg_FWHM = zeros(length(AP_list), 1);
avg_sensitivity = zeros(length(AP_list), 1);
avg_SNR = zeros(length(AP_list), 1);
AP_number = zeros(length(AP_list), 1);
ROI_number = zeros(length(AP_list), 1);

AP_data.amp = {};
AP_data.FWHM = {};
AP_data.sensitivity = {};
AP_data.SNR = {};
AP_data.index = {};


% 存储在tables中
table_name = fullfile(save_path,'AP_data.xlsx');
for i = 1:length(AP_list)
    
    if cellfun('isempty',AP_list{i}) == 0
        AP_i = AP_list{i}; % 当前trace的所有APs

        % 初始化每个trace的数据向量
        number_i = zeros(length(AP_i), 1);
        amp_i = zeros(length(AP_i), 2*AP_window_width+1);
        FWHM_i = zeros(length(AP_i), 1);
        sensitivity_i = zeros(length(AP_i), 1);
        SNR_i = zeros(length(AP_i), 1);
        index_i = zeros(length(AP_i), 1);

        for j = 1:length(AP_i)
            if ~isempty(peaks_index_manually_gated)
                each_AP = AP_i{j};
                number_i(j) = each_AP.AP_number;
                amp_i(j,:)  = each_AP.AP_amp .* peaks_index_manually_gated{i}(j);
                FWHM_i(j)  = each_AP.FWHM*peaks_index_manually_gated{i}(j);
                sensitivity_i(j)  = each_AP.Sensitivity*peaks_index_manually_gated{i}(j);
                SNR_i(j)  = each_AP.SNR*peaks_index_manually_gated{i}(j);
                index_i(j) = peaks_index{i}(j)*peaks_index_manually_gated{i}(j);
            else
                each_AP = AP_i{j};
                number_i(j) = each_AP.AP_number;
                amp_i(j,:)  = each_AP.AP_amp;
                FWHM_i(j)  = each_AP.FWHM;
                sensitivity_i(j)  = each_AP.Sensitivity;
                SNR_i(j)  = each_AP.SNR;
                index_i(j) = peaks_index{i}(j);
            end
        end

        % 为当前trace创建一个表格
        T = table(number_i, amp_i(:,2*AP_window_width+1), FWHM_i, sensitivity_i, SNR_i, peaks_index{i}, ...
            'VariableNames', {'Number', 'Amplitude', 'FWHM (ms)', 'Sensitivity', 'SNR', 'Index'});

        % 将表格写入Excel的一个新工作表
        sheet_name = string(['ROI ' num2str(i)]);
        writetable(T,table_name, 'Sheet', sheet_name);

        % save average value
        % avg_amp(i) = mean(amp_i,'omitmissing');
        avg_FWHM(i) = mean(FWHM_i,'omitmissing');
        avg_sensitivity(i) = mean(sensitivity_i,'omitmissing');
        avg_SNR(i) = mean(SNR_i,'omitmissing');
        AP_number(i) = number_i(end) - sum(isnan(peaks_index_manually_gated{i}));
        ROI_number(i) = i;
        AP_data.amp{i} = amp_i;
        AP_data.FWHM{i} = FWHM_i;
        AP_data.sensitivity{i} = sensitivity_i;
        AP_data.SNR{i} = SNR_i;
        AP_data.index{i} = index_i;
    end
end
fprintf('Finished statistic AP\n')
%% Save results

figure()
FWHM_axe = subplot(1,3,3);hold on;xlim([0,nrois+1]);
sensitivity_axe = subplot(1,3,2);hold on;xlim([0,nrois+1]);
SNR_axe = subplot(1,3,1);hold on;xlim([0,nrois+1]);
sgtitle('AP statistic');

% 初始化数据向量和分组标签
allFWHM = [];
alldff = [];
allSNR = [];
Labels = [];

% 遍历每个 cell
for i = 1:nrois
    % 获取当前 cell 的数据
    currentFWHM = AP_data.FWHM{i};
    currentdff = AP_data.sensitivity{i};
    currentSNR = AP_data.SNR{i};
    % 合并数据
    allFWHM  = [allFWHM; currentFWHM(:)];
    alldff  = [alldff; currentdff(:)];
    allSNR  = [allSNR; currentSNR(:)];
    % 生成分组标签（例如：第1个cell标签为1，第2个为2，依此类推）
    Labels = [Labels; i * ones(length(currentFWHM), 1)];
end
boxchart(Labels, allFWHM, 'Parent',FWHM_axe,'MarkerStyle','x','JitterOutliers','on');

xlabel('ROI number','Parent',FWHM_axe);
ylabel('FWHM (ms)','Parent',FWHM_axe);

boxchart(Labels, alldff*-1, 'Parent',sensitivity_axe,'MarkerStyle','x','JitterOutliers','on');hold on;
xlabel('ROI number','Parent',sensitivity_axe);
ylabel('Sensitiviy','Parent',sensitivity_axe);

boxchart(Labels, allSNR,'Parent',SNR_axe,'MarkerStyle','x','JitterOutliers','on');hold on;
xlabel('ROI number','Parent',SNR_axe);
ylabel('SNR','Parent',SNR_axe);

fig_filename = fullfile(save_path, '5_AP statistic.fig');
png_filename = fullfile(save_path, '5_AP statistic.png');

saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');

T_ave = table(ROI_number, AP_number,  avg_FWHM, avg_sensitivity, avg_SNR, ...
    'VariableNames', {'ROI Number','AP Number',  'Average FWHM (ms)', 'Average Sensitivity', 'Average SNR'});
writetable(T_ave, table_name, 'Sheet', 'Average');
AP_data_filename = fullfile(save_path, 'AP_data.mat');
save(AP_data_filename, "AP_data",'AP_list')
fprintf('Finished statistic AP\n')
% end

fprintf('AP_data.xlsx saved.\n');
%% Trend
nframe = size(traces_SNR,1);
trendlength = floor(nframe/freq);
trendbin = 1;
trendpart = trendlength/trendbin;
trend_avgSNR = zeros(1,trendpart);
trend_stdSNR = zeros(1,trendpart);
all_avgSNR = zeros(nrois,trendpart);

trend_avgdff = zeros(1,trendpart);
trend_stddff = zeros(1,trendpart);
all_avgdff = zeros(nrois,trendpart);

trend_avgFWHM = zeros(1,trendpart);
trend_stdFWHM = zeros(1,trendpart);
all_avgFWHM = zeros(nrois,trendpart);

trend_avgFR = zeros(1,trendpart);
trend_stdFR = zeros(1,trendpart);
all_avgFR = zeros(nrois,trendpart);

trend_ISI = cell(nrois,trendpart);
all_ISI = cell(nrois,1);

for t = 1:trendpart
    startindex = (t-1) *freq + 1;
    endindex = t  *freq;
    current_avgSNR = zeros(1,trendpart);
    current_avgdff = zeros(1,trendpart);
    current_avgFWHM = zeros(1,trendpart);
    current_avgFR = zeros(1,trendpart);
    
    % current_stdSNR = zeros(1,24);
    for i = 1:nrois
        indice = find((AP_data.index{i} >= startindex) & (AP_data.index{i} <= endindex));
        current_avgSNR(i) = mean(AP_data.SNR{i}(indice),'omitmissing');
        % current_stdSNR(i) = std(AP_data.SNR{i}(indice),'omitmissing');
        all_avgSNR(i,t) = current_avgSNR(i);

        current_avgdff(i) = mean(AP_data.sensitivity{i}(indice),'omitmissing');
        all_avgdff(i,t) = current_avgdff(i);

        current_avgFWHM(i) = mean(AP_data.FWHM{i}(indice),'omitmissing');
        all_avgFWHM(i,t) = current_avgFWHM(i);

        current_avgFR(i) = sum(~isnan(indice));
        all_avgFR(i,t) = current_avgFR(i);
        

        all_ISI{i} = diff(AP_data.index{i})/freq;
        trend_ISI{i,t} = all_ISI{i}(indice(indice<=length(all_ISI{i})));
    end

    trend_avgSNR(t) = mean(current_avgSNR,'omitmissing');
    trend_stdSNR(t) = std(current_avgSNR,'omitmissing');

    trend_avgdff(t) = mean(current_avgdff,'omitmissing');
    trend_stddff(t) = std(current_avgdff,'omitmissing');
    
    trend_avgFWHM(t) = mean(current_avgFWHM,'omitmissing');
    trend_stdFWHM(t) = std(current_avgFWHM,'omitmissing');
    
    trend_avgFR(t) = mean(current_avgFR,'omitmissing');
    trend_stdFR(t) = std(current_avgFR,'omitmissing');


end

figure()
trendx = 1:trendbin:trendlength;
subplot(1,5,1)
title('SNR');hold on;
for i = 1:nrois
    plot(trendx, all_avgSNR(i,:),'Color',[0.8,0.8,0.8])
end
plot(trendx, trend_avgSNR,'k', 'LineWidth', 2);
errorbar(trendx,  trend_avgSNR, trend_stdSNR, 'k', 'LineStyle', 'none', 'LineWidth', 1,'CapSize',10); % 将误差转换为百分比，加粗误差线

xlabel('Time (s)')
ylabel('SNR')

subplot(1,5,2)
title('Sensitivity');hold on;
for i = 1:nrois
    plot(trendx, all_avgdff(i,:)*-1,'Color',[0.8,0.8,0.8])
end
plot(trendx, trend_avgdff*-1,'k', 'LineWidth', 2);
errorbar(trendx,  trend_avgdff*-1, trend_stddff*-1, 'k', 'LineStyle', 'none', 'LineWidth', 1,'CapSize',10); % 将误差转换为百分比，加粗误差线

xlabel('Time (s)')
ylabel('Sensitivity (-%)')


subplot(1,5,3)
title('FWHM');hold on;
for i = 1:nrois
    plot(trendx, all_avgFWHM(i,:)*1000,'Color',[0.8,0.8,0.8])
end
plot(trendx, trend_avgFWHM*1000,'k', 'LineWidth', 2);
errorbar(trendx,  trend_avgFWHM*1000, trend_stdFWHM*1000, 'k', 'LineStyle', 'none', 'LineWidth', 1,'CapSize',10); % 将误差转换为百分比，加粗误差线

xlabel('Time (s)')
ylabel('FWHM (ms)')

subplot(1,5,4)
title('Firing rate');hold on;
for i = 1:nrois
    plot(trendx, all_avgFR(i,:),'Color',[0.8,0.8,0.8])
end
plot(trendx, trend_avgFR,'k', 'LineWidth', 2);
errorbar(trendx,  trend_avgFR, trend_stdFR, 'k', 'LineStyle', 'none', 'LineWidth', 1,'CapSize',10); % 将误差转换为百分比，加粗误差线

xlabel('Time (s)')
ylabel('Firing rate (Hz)')

subplot(1,5,5)
title('Firing rate');hold on;
boxchart(all_avgFR','MarkerStyle','x');
xlabel('ROI number')
ylabel('Firing rate (Hz)')

saveas(gcf,fullfile(save_path,'Trend Analysis.png'))
saveas(gcf,fullfile(save_path,'Trend Analysis.fig'))

save(fullfile(save_path,'Trend Analysis.mat'),'all_avgdff', 'all_avgFR', 'all_avgFWHM', 'all_avgSNR', ...
    'trend_avgdff', 'trend_avgFR', 'trend_avgFWHM', 'trend_avgSNR', 'trend_avgSNR', 'trend_stddff', 'trend_stdFR', 'trend_stdFWHM', 'trend_stdSNR')

%% Plot average AP sensitivity
figure();
% 统计不为空的trace数目
plot_cols = sum(cellfun('isempty',AP_list)==0)+1;
plot_col = 0;
% Initialize arrays to store all traces for final average calculation
trace_AP_mean = [];

for i = 1:nrois % i for trace
    %判断是否为有AP的trace
    peaks_num = length(peaks_index{i});
    if cellfun(['isempt' ...
            'y'],AP_list{i}) == 0
        plot_col = plot_col + 1;
        subplot(2,ceil(plot_cols/2),plot_col);
        set(gca,'color','none');

        % set y axis direction
        if peaks_polarity{i} == -1
            set(gca,'YDir','reverse')
            hold on;
        end

        % get each AP
        AP_i = zeros(peaks_num, AP_window_width*2+1);
        for j = 1:peaks_num
            each_AP = AP_list{i}{j};
            AP_i(j,:) = each_AP.AP_sensitivity;
            % plot each AP
            % plot((1:AP_window_width*2+1)*dt, each_AP.AP_sensitivity','Color',[0.8 0.8 0.8]);
            hold on;
        end

        AP_mean = mean(AP_i, 1, 'omitnan');
        AP_sd = std(AP_i, 0, 1, 'omitnan');

        % plot average AP for each trace
        subplot(2,ceil(plot_cols/2),plot_col);
        plot((1:AP_window_width*2+1)*dt, AP_mean,'Color',colors(i,:),'LineWidth',2);
        hold on;
        title(sprintf('ROI %d\n',i));
        fill([(1:AP_window_width*2+1)*dt, fliplr((1:AP_window_width*2+1)*dt)], ...
            [AP_mean + AP_sd, fliplr(AP_mean - AP_sd)], ...
            colors(i,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        hold on;

        % plot average AP for average trace

        % plot((1:AP_window_width*2+1)*dt, mean(AP_i,1,'omitnan'),'Color',colors(i,:),'LineWidth',1);
        % hold on;
        trace_AP_mean = [trace_AP_mean; AP_mean*peaks_polarity{i}];
    end
end

if ~isempty(trace_AP_mean)
    overall_mean = mean(trace_AP_mean, 1, 'omitnan');
    overall_sem = std(trace_AP_mean, 0, 1, 'omitnan') / sum(~cellfun('isempty', AP_list));
    subplot(2, ceil(plot_cols / 2), plot_cols);set(gca,'color','none');hold on;
    fill([(1:AP_window_width*2+1)*dt, fliplr((1:AP_window_width*2+1)*dt)], ...
        [overall_mean + overall_sem, fliplr(overall_mean - overall_sem)], ...
        [0.8 0.8 0.8], 'EdgeColor', 'none');

    title('Averaged of All');
    plot((1:AP_window_width*2+1) * dt, overall_mean, 'Color', 'k', 'LineWidth', 1);
    hold on;
end

sgtitle('Averaged Sensitivity');
hold on;
fig_filename = fullfile(save_path, '6_average_AP_sensitivity.fig');
png_filename = fullfile(save_path, '6_average_AP_sensitivity.png');
mat_filename = fullfile(save_path, '6_average_AP_sensitivity.mat');
saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');
save(mat_filename, 'AP_list', 'peaks_index', 'nrois', 'AP_window_width', 'dt', ...
    'peaks_polarity', 'colors', 'trace_AP_mean', 'overall_mean', 'overall_sem');

%%
% Plot average AP SNR with SD
figure();
plot_cols = sum(cellfun('isempty', AP_list) == 0) + 1;
plot_col = 0;
% Initialize arrays to store all traces for final average calculation
trace_AP_SNR_mean = [];

for i = 1:nrois
    peaks_num = length(peaks_index{i});
    if ~isempty(AP_list{i})
        plot_col = plot_col + 1;
        subplot(2, ceil(plot_cols / 2), plot_col);
        set(gca,'color','none');

        % Set y axis direction if necessary
        if peaks_polarity{i} == -1
            set(gca, 'YDir', 'reverse');
            hold on;
        end

        % Get each AP and compute mean & SEM
        AP_i = zeros(peaks_num, AP_window_width*2 + 1);
        for j = 1:peaks_num
            each_AP = AP_list{i}{j};
            AP_i(j,:) = each_AP.AP_SNR;
            %plot((1:AP_window_width*2+1) * dt, each_AP.AP_SNR', 'Color', [0.8 0.8 0.8]);
            hold on;
        end

        AP_mean = mean(AP_i, 1, 'omitnan');
        AP_sd = std(AP_i, 0, 1, 'omitnan');

        % Plot mean with Sd
        plot((1:AP_window_width*2+1) * dt, AP_mean, 'Color', colors(i,:), 'LineWidth', 2);
        fill([(1:AP_window_width*2+1)*dt, fliplr((1:AP_window_width*2+1)*dt)], ...
            [AP_mean + AP_sd, fliplr(AP_mean - AP_sd)], ...
            colors(i,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        title(sprintf('ROI %d', i));
        hold on;

        % % Plot average AP for all traces
        % subplot(2, ceil(plot_cols / 2), plot_cols);
        % plot((1:AP_window_width*2+1) * dt, AP_mean*peaks_polarity(i), 'Color', [0.8 0.8 0.8]);
        % hold on;

        % Collect traces for final average calculation
        trace_AP_SNR_mean = [trace_AP_SNR_mean; AP_mean*peaks_polarity{i}];
    end

end

% Plot overall average and SEM in the last subplot
if ~isempty(trace_AP_SNR_mean)
    overall_mean = mean(trace_AP_SNR_mean, 1, 'omitnan');
    overall_sem = std(trace_AP_SNR_mean, 0, 1, 'omitnan') / sum(~cellfun('isempty', AP_list));
    subplot(2, ceil(plot_cols / 2), plot_cols);hold on;
    set(gca,'color','none');
    fill([(1:AP_window_width*2+1)*dt, fliplr((1:AP_window_width*2+1)*dt)], ...
        [overall_mean + overall_sem, fliplr(overall_mean - overall_sem)], ...
        [0.8 0.8 0.8], 'EdgeColor', 'none');

    title('Averaged of All');
    plot((1:AP_window_width*2+1) * dt, overall_mean, 'Color', 'k', 'LineWidth', 1);
    hold on;
end
sgtitle('Average SNR with SD');
fig_filename = fullfile(save_path, '6_average_AP_SNR_with_SD.fig');
png_filename = fullfile(save_path, '6_average_AP_SNR_with_SD.png');

saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');

%


%% Save parameter
% 定义保存路径和文件名
save_filename = fullfile(save_path, '-1_workspace_variables.mat');

% % 保存当前工作区中的所有变量到.mat文件
clear movie;
save(save_filename);



%%


