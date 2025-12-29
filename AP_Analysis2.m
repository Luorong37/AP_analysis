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
folder_path = 'V:\Luorong\Invivo\25.12.22 invivo dual-color\Cam2_Rec9_5%Red_dg_2025-12-22 20-01-22';
file = '\Cycle1';  % must add format.do not add '\' at las
bin = 1;
% ↓↓↓↓↓-----------Prompt user for frame rate------------↓↓↓↓↓
freq = 400; % Hz
gpu = true; % defined gpu open

if exist('movie','var')
    matim =true;
    [folder_path,file,exten] = fileparts(save_file);
    file = [file,exten];
else
    matim = false; % imaging via matlab
end
% loadtype = "mat"; % default "tif" = loading tif file, or "mat" = loading .mat file
% mat_name = 'LEDcyan_rec1_cycles1.mat';     % if loading .mat file

% -----------------------------------------------------------

% Create an analysis folder

file_path = fullfile(folder_path, file);
if isfolder(file_path)
    file_name = file;
    %file_dir = dir(file_path);
    %[~, ~, file_extension] = fileparts(file_dir(3).name);
else
    % [~, file_name, file_extension] = fileparts(file_path)
    % 定义位移参数保存路径
    [folder_path, file_name, fext] = fileparts(file_path);
    % [folder_path, file_name, ~] = fileparts(file_path);
end

% create a folder for analysis
save_path = fullfile(folder_path, strcat(file_name, '_Analysis'), nowtime);
mkdir(save_path);

% Load image file
if gpu
    gcp;
end

if ~matim
    fprintf("Start loading movie, please wait...\n");
    tload = tic;
    % if loadtype == "tif"
    [movie, ncols, nrows, nframes] = load_movie(file_path);
    % end
    % if loadtype == "mat"
    %     mat_path = fullfile(file_path, mat_name);
    %     load(mat_path);
    %     [ncols, nrows, nframes] = size(movie);
    % end
    fprintf("Finished loading movie after %d s\n", toc(tload));
else
    [ncols, nrows, nframes] = size(movie);
    movie = reshape(movie,ncols*nrows, []);
end

% Presetting
function [dt, colors, t, map, mask, options] = presetting(freq, nframes, movie, ncols, nrows)

% Define parameters
dt = 1 / freq; % Calculate time axis
% colors = [lines(7);hsv(5);spring(3);winter(3);gray(3)];
colors = lines(100);
t = (1:nframes) * dt;
options.colors = colors;
map = [];
mask = [];

end

try
    movie_vol_2D = reshape(mean(movie,2), ncols, nrows, []);
catch ME
    movie_vol_2D = squeeze(mean(movie,3));
end

avg_image  = (movie_vol_2D - min(movie_vol_2D(:))) / (max(movie_vol_2D(:)) - min(movie_vol_2D(:)));

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

%% Motion Correction (Support Save/Apply Shifts)
fprintf('Initializing Motion Correction...\n');
% --- 1. 配置参数 ---
apply_only = 0;      % 是否使用之前的运动校正shift参数
loadMC     = 0;      % 是否读取之前的运动校正结果
Norigid    = 0;      % 是否开启非刚性校正
hp         = 1;      % 是否开启高通滤波（用于辅助估算位移）默认开启
template   = [];     % 手动输入校正模板

% 确保 movie 是 3 维
if ismatrix(movie)
    movie = reshape(movie, ncols, nrows, []);
end

% NoRMCorre 基础配置
options_r = NoRMCorreSetParms('d1',ncols,'d2',nrows,'bin_width',200,'max_shift',10,'us_fac',30, ...
    'grid_size',[128,128],'overlap_pre',[32,32],'mot_uf',4,'max_dev', [5,5],'iter',1,'correct_bidir',false);

% 定义保存路径
[folder_path, file_name, fext] = fileparts(file_path);
shift_res_path = fullfile(save_path, 'motion_shifts_result.mat');
params_save_path = fullfile(save_path, 'motion_correction_para.mat');

% --- 2. 核心处理逻辑 ---
if loadMC
    [Mr, ncols, nrows, ~] = load_movie(mc_path);
else
    t1 = tic;
    
    % --- 内存优化：数据预处理 ---
    if ~isa(movie, 'single'), movie = single(movie); end
    movie = movie - min(movie(:)); % 原始数据保留在内存中

    if apply_only
        % --- 功能：直接应用位移 ---
        [shift_filename, shift_foldername] = uigetfile(save_path, '选择位移文件');
        if isequal(shift_filename,0), return; end
        S = load(fullfile(shift_foldername, shift_filename));
        
        fprintf(' -> Mode: Apply existing shifts...\n');
        Mr = apply_shifts(movie, S.shifts_r, options_r);
        if Norigid && isfield(S, 'shifts_nr'), Mr = apply_shifts(Mr, S.shifts_nr, S.options_nr); end
    else
        % --- 功能：重新估算并保存 ---
        if hp
            % 【高通滤波模式】
            fprintf(' -> High-pass mode: Filtering for better estimation...\n');
            Y = create_temp_highpass(movie); 
        else
            Y = movie;
        end

        fprintf(' -> Mode: Estimating shifts...\n');
        % 1. 刚体校正：用滤波后的图算位移(shifts_r)，但应用到原始 movie 上得到 Mr
        [~, shifts_r, template1] = normcorre_batch(Y, options_r);
        Mr = apply_shifts(movie, shifts_r, options_r);
        
        clear Y; % 估算完立即释放临时滤波数据

        % 2. 非刚体校正 (根据需要)
        shifts_nr = []; options_nr = [];
        if Norigid
            options_nr = NoRMCorreSetParms('d1',ncols,'d2',nrows,'bin_width',200,'grid_size',[128,128]);
            [~, shifts_nr, ~] = normcorre_batch(Y, options_nr, template1);
            Mr = apply_shifts(Mr, shifts_nr, options_nr);
      
        end
        
        % 保存
        fprintf(' -> Saving shifts and params...\n');
        save(shift_res_path, 'shifts_r', 'shifts_nr', 'template1', '-v7.3');
        save(params_save_path, 'options_r', 'options_nr', 'hp', 'Norigid');
    end

    % 保存校正后的 TIFF
    save_name = fullfile(folder_path, [file_name, '_motion_correction.tif']);
    array2tif(uint16(Mr), save_name);
    movie_vol_2D = mean(Mr, 3);
    fprintf('Finished in %d s\n', round(toc(t1)));
end

% --- 3. 后处理 ---
avg_image = (movie_vol_2D - min(movie_vol_2D(:))) ./ (max(movie_vol_2D(:)) - min(movie_vol_2D(:)));

% --- 辅助子函数 ---
function Y_hp = create_temp_highpass(movie)
    % 内存中快速创建高通滤波版本
    gSig = 7; gSiz = 17;
    psf = fspecial('gaussian', round(2*gSiz), gSig);
    ind = (psf >= max(psf(:,1)));
    psf = psf - mean(psf(ind)); psf(~ind) = 0;
    Y_hp = imfilter(movie, psf, 'symmetric');
end

% % this section is drived from demo_1p_low_RAM of NoRMCorre
% fprintf('Motion correcting...\n');
% loadMC = 0;
% Norigid = 0;
% hp = 0;
% template = [];
% % template = reshape(mean(movie(:,6300:6400),2),ncols,nrows);
% [folder_path, file_name, ~] = fileparts(file_path);
% if loadMC
%     mc_path = 'E:\1_Data\Luorong\25.12.11 invivo dual-color\Cam2_Rec7_25%Red_dg_2025-12-11 21-06-00\Cycle1_Analysis\2025-12-13 15-02-57\0_motion_correction_LEDred_rec7_cycles1_stack01_stack01.tif';
%     [movie, ncols, nrows, nframes] = load_movie(mc_path);
% 
% else
%     t1 = tic;
%     % [M1f,shifts1] = motion_correction(reshape(movie, ncols, nrows, []),'file_path',file_path);
%     d1 = ncols;
%     d2 = nrows;
% 
%     init_batch = 0;
%     bin_width = 100;
%     max_shift = [20 20];
%     grid_size = [128,128];
%     overlap_pre = [32,32];
%     mot_uf = 4;
%     iter = 1;
%     max_dev = [32,32];
%     us_fac = 50;
%     correct_bidir = false;
% 
%     border_nan = 'copy';
%     gSig_filt = [3,3];
%     shifts_method = 'cubic';
% 
%     % appropriate name through options_r.tiff_filename or options_r.h5_filename
%     options_r = NoRMCorreSetParms('d1',d1,'d2',d2,'bin_width',bin_width,'max_shift',max_shift,'grid_size',grid_size, ...
%         'overlap_pre', overlap_pre,'mot_uf' , mot_uf,'iter',iter,'max_dev',max_dev);
% 
%     options_nr = NoRMCorreSetParms('d1',d1,'d2',d2,'bin_width',20, ...
%         'grid_size',grid_size,'mot_uf',mot_uf,'correct_bidir',correct_bidir, ...
%         'overlap_pre',overlap_pre,'overlap_post',overlap_pre,'max_shift',max_shift);
% 
%     if isfolder(file_path)
%         file_list = dir(fullfile(file_path, ['*',fext]));
%         file_nums = length(file_list);
% 
%         movie_mc = [];
% 
%         % if file_nums ~= length(h5_list)
%         %     presave = logical([zeros(length(h5_list),1),ones(file_nums-length(h5_list),1)]);
%         % end
% 
%         for i = 1:file_nums
%             file = fullfile(file_path,file_list(i).name);
% 
%             Mr = normcorr(file,options_r,options_nr,Norigid,hp,template);
%             Mr = uint16(Mr);
%             % mc_filename = fullfile(save_path, ['0_motion_correction_',file_list(i).name]);
%             mc_filename = fullfile(folder_path, [file_list(i).name, '_motion_correction']);
%             % save_tiffs(Mr,mc_filename)
%             array2tif(Mr,mc_filename)
%             %save(mc_filename,'Mr','d1','d2','bin_width','max_shift','grid_size' ,'overlap_pre','mot_uf','iter','correct_bidir');
% 
%             movie_mc = cat(3,movie_mc,Mr);
% 
%             % movie = reshape(movie_mc, ncols*nrows, []);
%             % if presave(i)
%             %     [h5_name] = create_highpass_temple(file);
%             % else
%             %     h5_name = fullfile(file_path,h5_list(i).name);
%             % end
%             % % [~, filename, ~] = fileparts(['\' file_list(i).name]);
%             % options_r = NoRMCorreSetParms('d1',d1,'d2',d2,'bin_width',bin_width,'max_shift',max_shift,'grid_size',grid_size,'overlap_pre', overlap_pre,'mot_uf' , mot_uf,'iter',iter,'correct_bidir',correct_bidir);
%             %
%             % %register using the high pass filtered data and apply shifts to original data
%             % tic; [~,shifts1,template1] = normcorre_batch(h5_name,options_r); toc % register filtered data
%             % % exclude boundaries due to high pass filtering effects
%             %
%             % % if you save the file directly in memory make sure you save it with a
%             % % name that does not exist. Change options_r.tiff_filename
%             % % or options_r.h5_filename accordingly.
%             % %
%             % tic;
%             %
%             % if Norigid
%             %
%             %
%             %     tic; [~,shifts2,template2] = normcorre_batch(M1,options_nr,template1); toc % register filtered data
%             %     tic; Mr = apply_shifts(movie,shifts2,options_nr,0,0); toc % apply the shifts to the removed percentile
%             % else
%             %
%             %     Mr = apply_shifts(file,shifts1,options_r); toc % apply shifts to full dataset
%             % end
% 
% 
%         end
%         movie_vol_2D = mean(movie_mc,3);
%         Mr = movie_mc;
%     else
%         Mr = normcorr(file_path,options_r,options_nr,Norigid,hp,template);
% 
%         % if presave
%         %     [h5_name] = create_highpass_temple(file_path);
%         % else
%         %     h5_name = fullfile(folder_path,h5_list.name);
%         % end
%         % [~, filename, ~] = fileparts(['\' file_path]);
%         %
%         % options_r = NoRMCorreSetParms('d1',d1,'d2',d2,'bin_width',bin_width,'max_shift',max_shift,'grid_size',grid_size, ...
%         %     'overlap_pre', overlap_pre,'mot_uf' , mot_uf,'iter',iter,'correct_bidir',correct_bidir);
%         %
%         %
%         % %register using the high pass filtered data and apply shifts to original data
%         % tic; [M1,shifts1,~] = normcorre_batch(h5_name,options_r); toc % register filtered data
%         % % exclude boundaries due to high pass filtering effects
%         %
%         % % if you save the file directly in memory make sure you save it with a
%         % % name that does not exist. Change options_r.tiff_filename
%         % % or options_r.h5_filename accordingly.
%         % %
%         % tic;
%         %
%         % if Norigid
%         %     options_nr = NoRMCorreSetParms('d1',d1,'d2',d2,'bin_width',bin_width, ...
%         %         'grid_size',grid_size,'mot_uf',mot_uf,'correct_bidir',correct_bidir, ...
%         %         'overlap_pre',32,'overlap_post',32,'max_shift',max_shift);
%         %
%         %     tic; [~,shifts2,template2] = normcorre_batch(M1,options_nr,template1); toc % register filtered data
%         %     tic; Mr = apply_shifts(movie,shifts2,options_nr,0,0); toc % apply the shifts to the removed percentile
%         % else
%         %
%         %
%         %     Mr = apply_shifts(file_path,shifts1,options_r); toc % apply shifts to full dataset
%         % end
% 
%         movie_vol_2D = mean(Mr,3);
%         % mc_filename = fullfile(save_path, ['0_motion_correction_',file_name,'.tif']);
%         mc_filename = fullfile(folder_path, [file_name, '_motion_correction.tif']);
%         save_tiffs(Mr,mc_filename)
%         % movie = reshape(Mr, ncols*nrows, []);
%         %save(mc_filename,'Mr','d1','d2','bin_width','max_shift','grid_size' ,'overlap_pre','mot_uf','iter','correct_bidir');
%     end
%     %first try out rigid motion correction
%     %     % exclude boundaries due to high pass filtering effects
%     % options_r = NoRMCorreSetParms('d1',d1,'d2',d2,'bin_width',200,'max_shift',[20 20],'grid_size' ,[128,128],'overlap_pre', [32,32],'mot_uf' , 4,'iter',1,'correct_bidir',false);
%     %
%     % %register using the high pass filtered data and apply shifts to original data
%     % tic; [M1,shifts1,template1] = normcorre_batch(h5_name,options_r); toc % register filtered data
%     %     % exclude boundaries due to high pass filtering effects
%     %
%     % % if you save the file directly in memory make sure you save it with a
%     % % name that does not exist. Change options_r.tiff_filename
%     % % or options_r.h5_filename accordingly.
%     % %
%     % tic; Mr = apply_shifts(file_path,shifts1,options_r); toc % apply shifts to full dataset
%     %
% 
% 
%     % avg_image  = (movie_vol_2D - min(movie_vol_2D(:))) / (max(movie_vol_2D(:)) - min(movie_vol_2D(:)));
% 
%     t2 = toc(t1);
%     fprintf('Finished motion correction after %d s\n',round(t2))
% 
% 
% end
% 
% 
% 
% avg_image  = (movie_vol_2D - min(movie_vol_2D(:))) ./ (max(movie_vol_2D(:)) - min(movie_vol_2D(:)));
% 
% mcmat_filename = fullfile(save_path, '0_motion_correction.mat');
% 
% save(mcmat_filename, 'options_r', 'options_nr');
% 
% function [h5_name] = create_highpass_temple(file_path)
% 
% gSig = 7;
% gSiz = 3*gSig;
% psf = fspecial('gaussian', round(2*gSiz), gSig);
% ind_nonzero = (psf(:)>=max(psf(:,1)));
% psf = psf-mean(psf(ind_nonzero));
% psf(~ind_nonzero) = 0;   % only use pixels within the center disk
% [filepath,file_name,~] = fileparts(file_path);
% h5_name = fullfile(filepath,[file_name,'_filtered_data.h5']);
% chunksize = 1000;    % read 500 frames at a time
% 
% cnt = 1;
% while (1)  % read filter and save file in chunks
%     Yf = single(read_file(file_path,cnt,chunksize));
%     if isempty(Yf)
%         break
%     else
%         Y = imfilter(Yf,psf,'symmetric');
%         saveash5(Y,h5_name);
%         cnt = cnt + size(Y,ndims(Y));
%     end
%     disp(cnt)
% end
% end
% 
% % save_tiffs(Mr,mc_filename)
% % function save_tiffs(movie,tiff_file_name)
% % fprintf('Saving results as tiff...\n')
% % tic;
% % t = Tiff(tiff_file_name, 'w');
% % batch_start = 1;
% % [ncols,nrows,batch_end] = size(movie);
% % for j = batch_start:batch_end
% %     current_image = movie(:,:,j);         % 读取当前图片
% %
% %     t.setTag('ImageLength', nrows);
% %     t.setTag('ImageWidth', ncols);
% %     t.setTag('Photometric', Tiff.Photometric.MinIsBlack);
% %     t.setTag('BitsPerSample', 16);
% %     t.setTag('SamplesPerPixel', 1);
% %     t.setTag('RowsPerStrip', 16);
% %     t.setTag('PlanarConfiguration', Tiff.PlanarConfiguration.Chunky);
% %     t.setTag('Compression', Tiff.Compression.None);
% %     t.setTag('Software', 'MATLAB');
% %     t.write(current_image);
% %
% %     if j < batch_end
% %         t.writeDirectory();
% %     end
% %
% % end
% % t.close();
% % toc;
% % end
% 
% %
% 
% 
% function Mr = normcorr(file_path,options_r,options_nr,Norigid,hp,template)
% 
% [folder_path, file_name, ~] = fileparts(file_path);
% movie = load_movie(file_path);
% 
% if hp
%     h5_list = dir(fullfile(folder_path, [file_name,'_filtered_data.h5']));
%     if isempty(h5_list)
%         fprintf('creating high-pass file......\n')
%         presave = true;
%     else
%         fprintf('previous high-pass file found.\n')
%         presave = false;
%     end
% 
%     if presave
%         [h5_name] = create_highpass_temple(file_path);
%     else
%         h5_name = fullfile(folder_path,h5_list.name);
%     end
% end
% % [~, filename, ~] = fileparts(['\' file_path]);
% 
% %register using the high pass filtered data and apply shifts to original data
% 
% % f = figure;
% % plot(detrend(mean(reshape(movie,size(movie,1)*size(movie,2),[]),1)))
% % title('please select the range of temple\n')
% % % [x,~] = ginput(2);
% % % template = uint16(mean(movie(:,:,x(1):x(2)),3));
% % close(f);
% % imshow(imadjust(template))
% % title('this image will be the template of motion correction\n')
% % fprintf('%d frames will be the template of motion correction\n', x(2)-x(1));
% fprintf('Rigid correcting...\n')
% tic;
% if hp
%     fprintf('high-pass data was used for motion correction\n');
%     if ~isempty(template)
%         [Mr,shifts1,template1] = normcorre_batch(h5_name,options_r,template); toc % register filtered data
%     else
%         [Mr,shifts1,template1] = normcorre_batch(h5_name,options_r); toc % register filtered data
%     end
%     % Mr = apply_shifts(M1,shifts1,options_nr,0,0);
% else
%     fprintf('raw data was used for motion correction\n');
%     if ~isempty(template)
% 
%         [Mr,~,~] = normcorre_batch(single(movie)-min(single(movie(:))),options_r,template); toc % register filtered data
%     else
%         [Mr,~,template1] = normcorre_batch(single(movie)-min(single(movie(:))),options_r); toc % register filtered data
% 
%     end
% end
% 
%     if Norigid
%         fprintf('Non-rigid correcting...\n')
%         tic;
%         [Mr,~,~] = normcorre_batch(Mr,options_nr,template1); toc % register filtered data
%         tic;
%         % Mr = apply_shifts(M1,shifts2,options_nr,0,0);
%         % Mr = apply_shifts(movie,shifts2,options_nr,0,0); toc % apply the shifts to the removed percentile
%         % else
%         %     tic; Mr = apply_shifts(movie,shifts1,options_r); toc % apply shifts to full dataset
%     end
% end


movie = reshape(uint16(Mr), ncols*nrows, []);
%% Create a map (optional)
t1 = tic; % Start a timer
fprintf('Creating a map...\n')
% if the map cannot figure out active cells, please large the bin.
mapbin = 4; % defined bin = 4

[quick_map] = create_map(movie, ncols, nrows, mapbin);
map = quick_map;

% Visualize correlation coefficients as heatmap
figure()
imagesc(quick_map);
colormap turbo;
title('Sensitivity Map');
axis image;
colorbar;
fig_filename = fullfile(save_path, '0_Sensitivity_Map.fig');
png_filename = fullfile(save_path, '0_Sensitivity_Map.png');
mat_filename = fullfile(save_path, '0_Sensitivity_Map.mat');

saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');
save(mat_filename, 'map');


t2 = toc(t1); % Get the elapsed time
fprintf('Finished mask creating after %d s\n',round(t2))
%% Load Mask (optional)
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
        avgdia = 25;
        % gamma_image = imadjust(map,[],[],1.2); % recommend raise the gamma factor from 1 to 4.
        gamma_image = imadjust(avg_image,[],[],1.2); % recommend raise the gamma factor from 1 to 4.
        mask = segmentCells2D(cp, gamma_image , ImageCellDiameter = avgdia,FlowErrorThreshold = 3,CellThreshold = -6);% CellThreshold = -2, ,  FlowErrorThreshold = 2
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
% if exist('methods','var')
%     if any(strcmp(methods,{'No',''} ))
%         [rois, traces] = select_ROI(movie, ncols, nrows, mask, map);
%         nrois = max(rois.bwmask,[],'all');
%     else
%         [~, traces] = select_ROI(movie, ncols, nrows, mask, map);
%         nois = max(mask(:));
%         rois.bwmask = mask;
%     end
% else
%     [rois, traces] = select_ROI(movie, ncols, nrows, mask, map);
%     nrois = max(rois.bwmask,[],'all');
% end

[rois, traces] = select_ROI(movie, ncols, nrows, mask, map);
nrois = max(rois.bwmask,[],'all');

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


tif_filename = fullfile(save_path, '1_Averaged.tif');

tiffile = Tiff(tif_filename, 'w');

tiffile.setTag('ImageLength', nrows);
tiffile.setTag('ImageWidth', ncols);
tiffile.setTag('Photometric', Tiff.Photometric.MinIsBlack);
tiffile.setTag('BitsPerSample', 16);
tiffile.setTag('SamplesPerPixel', 1);
tiffile.setTag('RowsPerStrip', 16);
tiffile.setTag('PlanarConfiguration', Tiff.PlanarConfiguration.Chunky);
tiffile.setTag('Compression', Tiff.Compression.None);

tiffile.setTag('Software', 'MATLAB');
tiffile.write(uint16(movie_vol_2D));
tiffile.writeDirectory();
tiffile.close();


%close(gcf);
%% Signal Process %%
% Background Correction
fprintf('Removing Background...\n')
[background, background_fitted, traces_bgcorr, traces_bgfitcorr, background_mask]...
    = remove_background(movie, ncols, nrows, rois, freq, bin);
roi_filename = fullfile(save_path, '1_background_ROI.mat');
save(roi_filename, 'background','background_fitted','background_mask','traces_bgcorr','traces_bgfitcorr');
fprintf('Finished\n')

% --- 补全画图功能 ---
fprintf('Generating summary plots...\n')

% 1. 初始化设置
num_rois = max(rois.bwmask(:));
colors = lines(num_rois); % 生成颜色矩阵，确保每个ROI有唯一颜色

figure('Name', 'Background Correction Summary', 'Units', 'normalized', 'Position', [0.1 0.1 0.8 0.7]);

% 2. 绘制左图：原始图像 + ROI边界 + 背景掩码
subplot(1, 2, 1);
movie2D_mean = mean(reshape(movie, ncols, nrows, []), 3);
imshow(movie2D_mean, [], 'InitialMagnification', 'fit'); 
hold on;

for i = 1:num_rois
    % 获取当前 ROI 的颜色
    current_color = colors(i, :);
    
    % 绘制原始 ROI 边界 (使用实线)
    roi_bw = (rois.bwmask == i);
    roi_boundaries = bwboundaries(roi_bw);
    for k = 1:length(roi_boundaries)
        boundary = roi_boundaries{k};
        plot(boundary(:, 2), boundary(:, 1), 'Color', current_color, 'LineWidth', 1, 'DisplayName', ['ROI ', num2str(i)]);
    end
    
    % 绘制背景掩码区域 (使用点状或半透明填充)
    bg_bw = (background_mask == i);
    bg_boundaries = bwboundaries(bg_bw);
    for k = 1:length(bg_boundaries)
        boundary = bg_boundaries{k};
        plot(boundary(:, 2), boundary(:, 1), ':', 'Color', current_color, 'LineWidth', 1);
    end
    
    % 在 ROI 中心标序号
    stats = regionprops(roi_bw, 'Centroid');
    if ~isempty(stats)
        text(stats(1).Centroid(1), stats(1).Centroid(2), num2str(i), ...
            'Color', current_color, 'FontWeight', 'bold', 'FontSize', 12, 'HorizontalAlignment', 'center');
    end
end
title('Image: ROI (Solid) & Background (Dotted)');

% 3. 绘制右图：信号轨迹对比
subplot(1, 2, 2);
hold on;

% 为了避免图例太乱，我们先定义几个占位符用于显示图例
p1 = plot(nan, nan, 'Color', [0.7 0.7 0.7]); 
p2 = plot(nan, nan, 'k--', 'LineWidth', 1);
p3 = plot(nan, nan, 'k', 'LineWidth', 1);

for i = 1:num_rois
    current_color = colors(i, :);
    
    % 计算原始信号
    raw_signal = traces_bgfitcorr(:, i) + background_fitted(:, i);
    
    % 绘制原始信号 (浅色背景线)
    plot(t, raw_signal, 'Color', [current_color, 0.3], 'LineWidth', 0.5);
    
    % 绘制拟合背景 (虚线)
    plot(t, background_fitted(:, i), '--', 'Color', current_color, 'LineWidth', 1);
    
    % 绘制校正后信号 (深色主线)
    plot(t, traces_bgfitcorr(:, i), 'Color', current_color, 'LineWidth', 1);
end

xlabel('Times(s)');
ylabel('Intensity');
title('Traces: Raw (Faded), BG (Dash), Corrected (Solid)');
grid on;

% 添加统一图例
legend([p1, p2, p3], {'Raw Signal', 'Fitted Background', 'Corrected Signal'}, 'Location', 'northeast');

% 保存图表
saveas(gcf, fullfile(save_path, 'background_correction_summary.png'));
fprintf('Summary plot saved to: %s\n', save_path);
%% Bleaching Correction

bleachmode = 'linear';% 'linear' 'highpass' 'exp2' 

fprintf('Correcting Bleaching (Mode: %s)...\n', bleachmode);

switch bleachmode
    case 'highpass'
        fc = 0.5/t(end); 
        [traces_corrected, baseline] = highpass_bleach_remove(traces_bgfitcorr, freq, fc);
        
    case 'linear'
        % 使用 detrend 并计算基线
        traces_corrected = detrend(traces_bgfitcorr, 1);
        baseline = traces_bgfitcorr - traces_corrected;
        
    case 'exp2'
        % 双指数拟合提取基线
        [~, baseline] = fit_exp2(traces_bgfitcorr);
        traces_corrected = traces_bgfitcorr - baseline;
end

% 统一调用绘图函数
plot_corrected(traces_corrected, traces_bgfitcorr, baseline, t, save_path);

fprintf('Finished Bleaching Correction.\n');
function plot_corrected(traces_corrected, traces, baseline, t, save_path)
    % 获取 ROI 数量
    nrois = size(traces_corrected, 2);
    
    % 创建画布
    fig = figure('Name', 'Bleaching Correction Overview', 'Color', 'w');
    set(fig, 'Position', get(0, 'Screensize'));
    
    % --- 1. 左侧：原始数据 + 拟合基线 (Stacked) ---
    ax_fit = subplot(1, 2, 1); hold on;
    % 根据数据波动自动计算垂直间距
    spacing_raw = mean(std(traces, 0, 1)) * 5; 
    
    for r = 1:nrois
        offset = (nrois - r) * spacing_raw;
        % 原始数据 (彩色/薄线)
        plot(t, traces(:, r) + offset, 'LineWidth', 0.5, 'DisplayName', 'Original');
        % 拟合基线 (红色/粗线)
        plot(t, baseline(:, r) + offset, 'r', 'LineWidth', 1.2, 'DisplayName', 'Baseline');
        
        if mod(r, 5) == 0 || r == 1 || r == nrois
            text(t(1), offset, [' ROI ', num2str(r)], 'FontSize', 8, 'FontWeight', 'bold');
        end
    end
    title(['Original Traces & Fitted Baselines (N=', num2str(nrois), ')']);
    xlabel('Time (s)'); ylabel('Stacked Magnitude');
    grid on; axis tight;

    % --- 2. 右侧：校正后的信号 (Stacked) ---
    ax_corr = subplot(1, 2, 2); hold on;
    spacing_corr = mean(std(traces_corrected, 0, 1)) * 8;
    
    for r = 1:nrois
        offset = (nrois - r) * spacing_corr;
        % 校正后的信号 (黑色)
        plot(t, traces_corrected(:, r) + offset, 'k', 'LineWidth', 0.5);
        % 零位基准线 (浅蓝色)
        line([t(1) t(end)], [offset offset], 'Color', [0.3 0.7 1, 0.5], 'LineStyle', '--');
    end
    title('Corrected Traces (Residuals)');
    xlabel('Time (s)'); ylabel('Stacked Magnitude');
    grid on; axis tight;
    
    % 联动 X 轴（缩放左侧时右侧同步）
    linkaxes([ax_fit, ax_corr], 'x');

    % 保存图片
    fig_filename = fullfile(save_path, '2_bleach_correction_stacked.fig');
    png_filename = fullfile(save_path, '2_bleach_correction_stacked.png');
    saveas(gcf, fig_filename, 'fig');
    saveas(gcf, png_filename, 'png');
end
% %% Bleaching Correction
% 
% bleachmode = 'linear';% 'linear' 'highpass' 'exp2'
% 
%      plot_corrected1(traces_corrected, traces, baseline, save_path)
%         % plot
%         fig = figure();
%         set(fig,'Position',get(0,'Screensize'));
%         fit_axe = subplot(1,3,1);
%         fited_axe = subplot(1,3,2);
%         baseline_axe = subplot(1,3,3);
% 
%         % plot
%         for i = 1: size(traces_corrected,2)
%             plot(traces(:,i),'Parent',fit_axe);
%             hold(fit_axe, 'on');
%             plot(traces_corrected(:,i),'Parent',fited_axe);
%             hold(fited_axe, 'on');
% 
%             plot( baseline(:,i),'LineWidth',1,'Parent',baseline_axe);
%             hold(baseline_axe, 'on');
% 
% 
%         end
%         hold off;
% 
%         % note
%         title(fit_axe, 'Original and Fitted Curves');
%         title(fited_axe, 'Corrected Traces');
%         ylim(fited_axe,[-200,+200])
%         title(baseline_axe, 'Baseline');
%         legend(fit_axe, 'Original Trace', 'Fitted Curve');
% 
%         fig_filename = fullfile(save_path, '2_fitted_trace.fig');
%         png_filename = fullfile(save_path, '2_fitted_trace.png');
% 
%         saveas(gcf, fig_filename, 'fig');
%         saveas(gcf, png_filename, 'png');
%     end
% 
% fprintf('Correcting Bleaching...\n')
% 
% %highpass bleach
% switch bleachmode
%     case 'highpass'
%         fc = 15/t(end);%15
%         [traces_corrected,baseline] = highpass_bleach_remove(traces_bgfitcorr,freq,fc );
%     case 'linear'
%         traces_corrected = detrend(traces_bgfitcorr,1 );
%         baseline = traces_bgfitcorr-traces_corrected;
%     case 'exp2'
%         % [traces_corrected, fitted_curves, params, gof] = fit_exp2(traces_bgfitcorr);
%         [traces_corrected,fit_y_fit] = fit_bleach(traces_bgfitcorr,x);
%         baseline = traces_bgfitcorr-traces_corrected;
% end
% 
% 
% plot_corrected(traces_corrected, traces, baseline, save_path);
% fprintf('Finished\n');
%% Wavelet process


% wavelet降噪
Dnmethods = 'FDR';
Dnlevel = 8;
Wavename = 'bior6.8'; % 替代墨西哥帽的离散基
fprintf('Wavelet Denoising...\n')
traces_denoised = wdenoise(traces_corrected, Dnlevel ,DenoisingMethod=Dnmethods,Wavelet = Wavename);

traces_filename = fullfile(save_path, '2_processed_traces.mat');
save(traces_filename,"traces_denoised", "traces_original", 'traces_corrected', 'baseline', 'Dnmethods','Dnlevel')
fprintf('Finished\n');

% 全 ROI 降噪效果可视化
fprintf('Plotting all ROIs comparison...\n');

nrois = size(traces_denoised, 2);
t = (1:size(traces_denoised, 1))';


% --- 堆叠对比图 (Stacked Traces) ---

figure('Name', 'Stacked ROI Comparison', 'Color', 'w', 'Position', [150, 150, 1000, 900]);
max_display = min(10, nrois);
roi_to_show = round(linspace(1, nrois, max_display)); % 均匀选取 10 个 ROI

spacing = max(traces_corrected(:)) * 0.8; % 设置垂直间距

hold on;
for i = 1:length(roi_to_show)
    r_idx = roi_to_show(i);
    offset = (i-1) * spacing;

    % 原始信号（灰色）
    plot(t, traces_corrected(:, r_idx) + offset, 'Color', [0.7 0.7 0.7], 'HandleVisibility', 'off');
    % 降噪信号（彩色）
    plot(t, traces_denoised(:, r_idx) + offset, 'LineWidth', 1);
end

set(gca, 'YTick', (0:max_display-1)*spacing, 'YTickLabel', string(roi_to_show));
title('Comparison of Selected ROIs (Stacked View)');
xlabel('Time Points'); ylabel('ROI Index');
legend('Denoised Signal');
grid on; hold off;

fprintf('Finished plotting all requested ROI comparisons.\n');
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

parts = 2;
MinPeakProminence_factor = 0.4;
MinPeakDistance_factor = 4;

findmode = 'cr';% find via denoised traces
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
MinPeakProminence_factor = 0.25;
pr = 1;
rr = n_traces;
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



%% manually gate (optional) - 多边形选择，支持按R重新绘制（简洁版）

peaks_index_manually_gated = peaks_index;

for i = 1:nrois
    figure()
    set(gcf,'Position',[0,0,2500,1800])
    hold on;

    % 绘制信号和峰值
    plot(traces_SNR(:,i).*peaks_polarity{i});
    peaks_x = peaks_index{i};
    peaks_y = traces_SNR(peaks_x,i).*peaks_polarity{i};
    plot(peaks_x, peaks_y,'v','MarkerFaceColor','r');

    title(sprintf('ROI %d: 绘制多边形选择要删除的峰值 | 按R重新绘制', i));

    % 循环直到用户满意
    redraw = true;
    while redraw
        % 绘制多边形
        poly = drawpolygon(gca);
        wait(poly); % 等待多边形绘制完成（双击结束）

        % 获取多边形顶点
        xv = poly.Position(:,1);
        yv = poly.Position(:,2);

        % 使用inpolygon判断峰值是否在多边形内
        [in, ~] = inpolygon(peaks_x, peaks_y, xv, yv);

        % 可视化被选中的峰值
        selected_plot = plot(peaks_x(in), peaks_y(in), ...
            'v', 'MarkerFaceColor', 'g', 'MarkerSize', 10);

        % 询问用户是否确认
        choice = questdlg(sprintf('选择了 %d 个峰值。确认删除吗？', sum(in)), ...
            '确认选择', '确认', '重新绘制(R)', '取消', '重新绘制(R)');

        switch choice
            case '确认'
                manual_gated_index = in;
                redraw = false;
                fprintf('ROI %d: 确认删除 %d 个峰值\n', i, sum(in));

            case '重新绘制(R)'
                % 删除当前多边形和标记
                delete(poly);
                delete(selected_plot);
                fprintf('ROI %d: 重新绘制...\n', i);

            case '取消'
                % 取消整个选择
                delete(poly);
                delete(selected_plot);
                fprintf('ROI %d: 取消选择\n', i);
                manual_gated_index = false(size(peaks_x));
                redraw = false;
        end
    end

    % 应用删除标记
    peaks_index_manually_gated{i}(manual_gated_index) = NaN;
    peaks_index_manually_gated{i}(~manual_gated_index) = 1;

    % 保存结果
    mkdir(fullfile(save_path,'Manually gated peaks'));
    saveas(gcf, fullfile(save_path,'Manually gated peaks', ...
        sprintf('noi %d, peaks %d.png', i, length(peaks_x))))
    saveas(gcf, fullfile(save_path,'Manually gated peaks', ...
        sprintf('noi %d, peaks %d.fig', i, length(peaks_x))))
    close(gcf)
end

save(fullfile(save_path,'Manually gated peaks.mat'),'peaks_index_manually_gated')
% peaks_index_manually_gated = peaks_index;
%
% for i = 1:nrois
%     figure()
%     set(gcf,'Position',[0,0,2500,1800])
%
%     hold on;
%     plot(traces_SNR(:,i).*peaks_polarity{i});
%     % peaks_x = peaks_index{i}(~isnan(AP_data.FWHM{i}));
%     peaks_x = peaks_index{i};
%     peaks_y = traces_SNR(peaks_x,i).*peaks_polarity{i};
%     plot(peaks_x, peaks_y,'v','MarkerFaceColor','r');
%     title(sprintf('noi %d, %d peaks. DELETE peaks in rectangle',i,length(peaks_x)));
%     % mkdir(fullfile(save_path,'FWHM gated peaks'));
%     % saveas(gcf,fullfile(save_path,'FWHM gated peaks',sprintf('noi %d, peaks %d.png',i,length(peaks_x))))
%     % saveas(gcf,fullfile(save_path,'FWHM gated peaks',sprintf('noi %d, peaks %d.fig',i,length(peaks_x))))
%
%     rect = drawrectangle(gca);
%     manual_gated_index = rect.Position(1) < peaks_x & peaks_x < rect.Position(1) + rect.Position(3) & ...
%         rect.Position(2) < peaks_y & peaks_y < rect.Position(2) + rect.Position(4);
%
%     fprintf(' %d peaks deleted.\n', sum(manual_gated_index));
%
%     peaks_index_manually_gated{i}(manual_gated_index) = NaN;
%     peaks_index_manually_gated{i}(~manual_gated_index) = 1;
%
%
%     plot(peaks_x(manual_gated_index), peaks_y(manual_gated_index),'v','MarkerFaceColor','g');
%
%     mkdir(fullfile(save_path,'Manually gated peaks'));
%     saveas(gcf,fullfile(save_path,'Manually gated peaks',sprintf('noi %d, peaks %d.png',i,length(peaks_x))))
%     saveas(gcf,fullfile(save_path,'Manually gated peaks',sprintf('noi %d, peaks %d.fig',i,length(peaks_x))))
%     close(gcf)
% end
% save(fullfile(save_path,'Manually gated peaks.mat'),'peaks_index_manually_gated')
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
trendbin = 5;
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
%% 可配置时间窗口的ROI分析（SNR和灵敏度趋势）

% ==================== 参数设置 ====================
window_size_seconds =60;  % 时间窗口大小（秒）- 可以修改这个值
freq = 400;                % 采样频率（Hz）
nrois = length(AP_data.FWHM); % ROI数量

% ==================== 分析计算 ====================

% 计算总时间
nframe = size(traces_SNR,1);
total_time = nframe / freq; % 总时长（秒）

% 计算窗口数量
num_windows = floor(total_time / window_size_seconds);

% 初始化存储结构
roi_windows_data = struct();

% 为每个ROI创建分析
for roi_idx = 1:nrois
    fprintf('正在处理 ROI %d/%d...\n', roi_idx, nrois);

    % 获取当前ROI的数据
    ap_indices = AP_data.index{roi_idx}; % AP的时间索引
    ap_SNR = AP_data.SNR{roi_idx}; % AP的SNR值
    ap_sensitivity = AP_data.sensitivity{roi_idx}; % AP的灵敏度值

    % 转换AP索引为时间（秒）
    ap_times = ap_indices / freq;

    % 初始化该ROI的窗口数据
    roi_avg_SNR = zeros(1, num_windows);
    roi_avg_sensitivity = zeros(1, num_windows);
    window_centers = zeros(1, num_windows);

    % 分析每个时间窗口
    for win_idx = 1:num_windows
        % 计算窗口时间范围
        win_start = (win_idx - 1) * window_size_seconds;
        win_end = win_idx * window_size_seconds;
        window_centers(win_idx) = (win_start + win_end) / 2; % 窗口中心时间

        % 找到在该窗口内的AP
        in_window = ap_times >= win_start & ap_times < win_end;

        if any(in_window)
            % 计算窗口内的平均值
            roi_avg_SNR(win_idx) = mean(ap_SNR(in_window), 'omitnan');
            roi_avg_sensitivity(win_idx) = mean(ap_sensitivity(in_window), 'omitnan');
        else
            % 窗口内没有AP，设为NaN
            roi_avg_SNR(win_idx) = NaN;
            roi_avg_sensitivity(win_idx) = NaN;
        end
    end

    % 存储该ROI的数据
    roi_windows_data(roi_idx).ROI_number = roi_idx;
    roi_windows_data(roi_idx).window_centers = window_centers;
    roi_windows_data(roi_idx).avg_SNR = roi_avg_SNR;
    roi_windows_data(roi_idx).avg_sensitivity = roi_avg_sensitivity;
    roi_windows_data(roi_idx).num_APs_total = length(ap_indices);

    % 计算窗口内的AP数量
    ap_counts = zeros(1, num_windows);
    for win_idx = 1:num_windows
        win_start = (win_idx - 1) * window_size_seconds;
        win_end = win_idx * window_size_seconds;
        in_window = ap_times >= win_start & ap_times < win_end;
        ap_counts(win_idx) = sum(in_window);
    end
    roi_windows_data(roi_idx).AP_counts = ap_counts;
end

% ==================== 创建动态文件名 ====================
% 根据窗口大小生成文件名后缀
window_suffix = sprintf('%ds', window_size_seconds);

% 创建保存图形的文件夹（包含窗口大小信息）
trend_folder = fullfile(save_path, sprintf('%s_window_trends', window_suffix));
if ~exist(trend_folder, 'dir')
    mkdir(trend_folder);
end

% ==================== 绘制每个ROI的窗口趋势图 ====================

% 为每个ROI绘制单独的图
for roi_idx = 1:nrois
    figure('Position', [100, 100, 1400, 800], ...
        'Name', sprintf('ROI %d - %s窗口分析', roi_idx, window_suffix));

    % SNR趋势
    subplot(2, 2, 1);
    plot(roi_windows_data(roi_idx).window_centers, roi_windows_data(roi_idx).avg_SNR, ...
        'b-o', 'LineWidth', 2, 'MarkerFaceColor', 'b', 'MarkerSize', 6);
    xlabel('时间 (秒)');
    ylabel('平均 SNR');
    title(sprintf('ROI %d: SNR随时间变化 (%s窗口)', roi_idx, window_suffix));
    grid on;

    % 灵敏度趋势
    subplot(2, 2, 2);
    plot(roi_windows_data(roi_idx).window_centers, roi_windows_data(roi_idx).avg_sensitivity, ...
        'r-o', 'LineWidth', 2, 'MarkerFaceColor', 'r', 'MarkerSize', 6);
    xlabel('时间 (秒)');
    ylabel('平均灵敏度 (\DeltaF/F₀, %)');
    title(sprintf('ROI %d: 灵敏度随时间变化 (%s窗口)', roi_idx, window_suffix));
    grid on;

    % AP数量趋势
    subplot(2, 2, 3);
    bar(roi_windows_data(roi_idx).window_centers, roi_windows_data(roi_idx).AP_counts, ...
        'FaceColor', [0.2, 0.6, 0.2]);
    xlabel('时间 (秒)');
    ylabel('AP数量');
    title(sprintf('ROI %d: 每%s窗口内的AP数量', roi_idx, window_suffix));
    grid on;

    % 信息汇总
    subplot(2, 2, 4);
    text_str = {sprintf('ROI %d 信息汇总:', roi_idx), ...
        sprintf('时间窗口: %s', window_suffix), ...
        sprintf('总AP数量: %d', roi_windows_data(roi_idx).num_APs_total), ...
        sprintf('平均SNR: %.2f ± %.2f', ...
        mean(roi_windows_data(roi_idx).avg_SNR, 'omitnan'), ...
        std(roi_windows_data(roi_idx).avg_SNR, 'omitnan')), ...
        sprintf('平均灵敏度: %.2f%% ± %.2f%%', ...
        mean(roi_windows_data(roi_idx).avg_sensitivity, 'omitnan'), ...
        std(roi_windows_data(roi_idx).avg_sensitivity, 'omitnan'))};
    text(0.1, 0.5, text_str, 'FontSize', 11, 'VerticalAlignment', 'middle');
    axis off;

    % 保存图形（文件名包含窗口大小信息）
    fig_filename = fullfile(trend_folder, sprintf('ROI_%d_%s_window_trends.fig', roi_idx, window_suffix));
    png_filename = fullfile(trend_folder, sprintf('ROI_%d_%s_window_trends.png', roi_idx, window_suffix));

    saveas(gcf, fig_filename);
    saveas(gcf, png_filename);
    close(gcf);
end

% ==================== 绘制所有ROI的综合趋势图（平均） ====================

figure('Position', [100, 100, 1200, 800], ...
    'Name', sprintf('所有ROI - %s窗口综合趋势', window_suffix));

% 准备所有ROI的平均SNR和灵敏度数据
all_SNR_data = zeros(num_windows, nrois);
all_sensitivity_data = zeros(num_windows, nrois);

for roi_idx = 1:nrois
    all_SNR_data(:, roi_idx) = roi_windows_data(roi_idx).avg_SNR';
    all_sensitivity_data(:, roi_idx) = roi_windows_data(roi_idx).avg_sensitivity';
end

% 计算所有ROI的平均和标准差
mean_SNR = mean(all_SNR_data, 2, 'omitnan');
std_SNR = std(all_SNR_data, 0, 2, 'omitnan');

mean_sensitivity = mean(all_sensitivity_data, 2, 'omitnan');
std_sensitivity = std(all_sensitivity_data, 0, 2, 'omitnan');

window_centers = roi_windows_data(1).window_centers;

% 绘制平均SNR趋势
subplot(2, 1, 1);
hold on;
plot(window_centers, mean_SNR, 'b-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
fill([window_centers, fliplr(window_centers)], ...
    [mean_SNR' - std_SNR', fliplr(mean_SNR' + std_SNR')], ...
    'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
xlabel('时间 (秒)');
ylabel('平均 SNR');
title(sprintf('所有ROI的平均SNR随时间变化 (%s窗口)', window_suffix));
legend('平均值', '标准差范围', 'Location', 'best');
grid on;

% 绘制平均灵敏度趋势
subplot(2, 1, 2);
hold on;
plot(window_centers, mean_sensitivity, 'r-o', 'LineWidth', 2, 'MarkerFaceColor', 'r');
fill([window_centers, fliplr(window_centers)], ...
    [mean_sensitivity' - std_sensitivity', fliplr(mean_sensitivity' + std_sensitivity')], ...
    'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
xlabel('时间 (秒)');
ylabel('平均灵敏度 (\DeltaF/F₀, %)');
title(sprintf('所有ROI的平均灵敏度随时间变化 (%s窗口)', window_suffix));
legend('平均值', '标准差范围', 'Location', 'best');
grid on;

% 保存综合图形（文件名包含窗口大小信息）
fig_filename = fullfile(trend_folder, sprintf('All_ROIs_%s_window_trends.fig', window_suffix));
png_filename = fullfile(trend_folder, sprintf('All_ROIs_%s_window_trends.png', window_suffix));
saveas(gcf, fig_filename);
saveas(gcf, png_filename);

% ==================== 输出数据到Excel文件 ====================

% 创建Excel文件名（包含窗口大小信息）
excel_filename = fullfile(save_path, sprintf('%s_window_analysis.xlsx', window_suffix));

% 为每个ROI创建独立的sheet
for roi_idx = 1:nrois
    % 准备数据表
    time_points = roi_windows_data(roi_idx).window_centers';
    avg_SNR = roi_windows_data(roi_idx).avg_SNR';
    avg_sensitivity = roi_windows_data(roi_idx).avg_sensitivity';
    ap_counts = roi_windows_data(roi_idx).AP_counts';

    % 创建表格
    data_table = table(time_points, avg_SNR, avg_sensitivity, ap_counts, ...
        'VariableNames', {'Window_Center_Time_s', 'Avg_SNR', 'Avg_Sensitivity_percent', 'AP_Count'});

    % 写入Excel（每个ROI一个sheet）
    sheet_name = sprintf('ROI_%d', roi_idx);
    writetable(data_table, excel_filename, 'Sheet', sheet_name);

    % 添加汇总信息
    summary_data = {sprintf('ROI %d 汇总信息', roi_idx);
        sprintf('时间窗口大小: %s', window_suffix);
        sprintf('总AP数量: %d', roi_windows_data(roi_idx).num_APs_total);
        sprintf('平均SNR (全局): %.2f', mean(avg_SNR, 'omitnan'));
        sprintf('平均灵敏度 (全局): %.2f%%', mean(avg_sensitivity, 'omitnan'))};

    % 写入汇总信息到第二列
    writecell(summary_data, excel_filename, 'Sheet', sheet_name, 'Range', 'F1');
end

% 创建汇总sheet（所有ROI的数据）
summary_table = table();
for roi_idx = 1:nrois
    time_points = roi_windows_data(roi_idx).window_centers';
    avg_SNR = roi_windows_data(roi_idx).avg_SNR';
    avg_sensitivity = roi_windows_data(roi_idx).avg_sensitivity';

    % 添加ROI编号列
    roi_numbers = repmat(roi_idx, length(time_points), 1);

    % 创建该ROI的数据表
    roi_table = table(roi_numbers, time_points, avg_SNR, avg_sensitivity, ...
        'VariableNames', {'ROI_Number', 'Window_Center_Time_s', 'Avg_SNR', 'Avg_Sensitivity_percent'});

    % 合并到汇总表
    if roi_idx == 1
        summary_table = roi_table;
    else
        summary_table = [summary_table; roi_table];
    end
end

% 写入汇总sheet
writetable(summary_table, excel_filename, 'Sheet', 'All_ROIs_Summary');

% 创建统计摘要sheet
stats_summary = table();
roi_numbers = (1:nrois)';
total_APs = zeros(nrois, 1);
mean_SNR_all = zeros(nrois, 1);
mean_sensitivity_all = zeros(nrois, 1);

for roi_idx = 1:nrois
    total_APs(roi_idx) = roi_windows_data(roi_idx).num_APs_total;
    mean_SNR_all(roi_idx) = mean(roi_windows_data(roi_idx).avg_SNR, 'omitnan');
    mean_sensitivity_all(roi_idx) = mean(roi_windows_data(roi_idx).avg_sensitivity, 'omitnan');
end

stats_summary.ROI_Number = roi_numbers;
stats_summary.Total_AP_Count = total_APs;
stats_summary.Mean_SNR = mean_SNR_all;
stats_summary.Mean_Sensitivity_percent = mean_sensitivity_all;
stats_summary.Window_Size_s = repmat(window_size_seconds, nrois, 1);

writetable(stats_summary, excel_filename, 'Sheet', 'Statistics_Summary');

% ==================== 保存MAT文件 ====================
mat_filename = fullfile(save_path, sprintf('%s_window_analysis.mat', window_suffix));
save(mat_filename, 'roi_windows_data', 'window_size_seconds', 'num_windows', 'freq');

% ==================== 创建参数配置文件 ====================
% 保存分析参数，便于后续追溯
config_filename = fullfile(save_path, sprintf('%s_window_config.txt', window_suffix));
fid = fopen(config_filename, 'w');
fprintf(fid, '=== 时间窗口分析参数配置 ===\n');
fprintf(fid, '分析时间: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf(fid, '时间窗口大小: %d 秒\n', window_size_seconds);
fprintf(fid, '采样频率: %d Hz\n', freq);
fprintf(fid, 'ROI数量: %d\n', nrois);
fprintf(fid, '总记录时间: %.2f 秒\n', total_time);
fprintf(fid, '时间窗口数量: %d\n', num_windows);
fprintf(fid, '输出文件:\n');
fprintf(fid, '  1. Excel文件: %s\n', excel_filename);
fprintf(fid, '  2. MAT文件: %s\n', mat_filename);
fprintf(fid, '  3. 图形文件夹: %s\n', trend_folder);
fclose(fid);

% ==================== 显示完成信息 ====================
fprintf('\n=== 分析完成 ===\n');
fprintf('时间窗口大小: %d 秒\n', window_size_seconds);
fprintf('Excel文件已保存: %s\n', excel_filename);
fprintf('包含以下sheet:\n');
fprintf('  1. ROI_X (每个ROI的详细数据)\n');
fprintf('  2. All_ROIs_Summary (所有ROI的合并数据)\n');
fprintf('  3. Statistics_Summary (统计摘要)\n');
fprintf('图形已保存至: %s\n', trend_folder);
fprintf('MAT数据文件: %s\n', mat_filename);
fprintf('参数配置文件: %s\n\n', config_filename);

% ==================== 使用说明 ====================
fprintf('=== 使用说明 ===\n');
fprintf('如需更改时间窗口大小，请修改代码开头的参数:\n');
fprintf('  将 window_size_seconds = %d; 改为所需的值\n', window_size_seconds);
fprintf('然后重新运行此代码段即可。\n');
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
%% drafting (optional)
% plot SNR
figure;
title('SNR');
hold on;
[~] = offset_plot(traces_SNR,t);

cycle_gd = t(end)/32;

cycles = 32;
starts = 0;
for i = 1:cycles/2
    fill([starts,starts+cycle_gd,starts+cycle_gd,starts],[0,0,sum(max(traces_SNR))*3,sum(max(traces_SNR))*3],'k','FaceAlpha',0.2)
    starts = starts+cycle_gd*2 ;
end

fig_filename = fullfile(save_path, '4_grafting_SNR.fig');
png_filename = fullfile(save_path, '4_grafting_SNR.png');
trace_filename = fullfile(save_path, '4_grafting_SNR.mat');

saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');


%% Save parameter
% 定义保存路径和文件名
save_filename = fullfile(save_path, '-1_workspace_variables.mat');

% % 保存当前工作区中的所有变量到.mat文件
clear movie;
save(save_filename);



%%


