%% Nanobody Kinetics analysis (single-exponential, cleaned)
% Author: Luorong, 20250817 GPT5
% Notes:
% - Association (on):  y(t) = A*(1 - exp(-k*c*t)) + bg
% - Dissociation (off):y(t) = A*exp(-k*t) + bg
% All downstream steps (SBR gating, normalization, averages, export) use A/k/bg consistently.

clear; clc;
gcp;
%% ----------------------- Loading -----------------------
t1 = tic; % Start a timer
nowtime = string(datetime('now', 'Format', 'yyyy_MM_dd_HH_mm_ss.SSS'));
nowtime = strrep(nowtime , ':', '-');
fprintf('Loading...\n');

% ===== User inputs =====
folder_path = 'I:\1_Data\sx. Kinetics\off'; % ends in '\'
file_name   = 'C1-Koff-3.tif';                     % must add format.
mask_name   = ''; % no format

% concentration (for association "on")
c    = 20*1e-9;   % 200 nM
type = 'off';      % 'on' or 'off'
freq = 1/180;     % Hz (frame rate)

% -----------------------
colors = [lines(7);hsv(5);spring(3);winter(3);gray(3)]; %#ok<NASGU> % (not used but kept)

% Build paths
file_path  = fullfile(folder_path, file_name);
split_path = split(file_name, '.');
mask_path  = fullfile(folder_path, mask_name);

if numel(split_path) > 1
    file_extension = string(split_path(end)); %#ok<NASGU>
    baseName = split_path(1);
    save_path = fullfile(folder_path, sprintf('%s_Analysis', char(baseName)), char(nowtime));
else
    file_extension = 'tif'; %#ok<NASGU>
    save_path = fullfile(sprintf('%s_Analysis', folder_path) , char(nowtime));
end
if ~exist(save_path, 'dir'); mkdir(save_path); end

% Load image file  -> movie: [nrows x ncols x nframes] or reshaped by helper
[movie, ncols, nrows, nframes] = load_movie(file_path);

% Time axis (minutes)
dt = 1 / freq / 60;   % frame period in minutes
t  = (1:nframes)' * dt;

t2 = toc(t1);
fprintf('Finished loading movie after %d s\n', round(t2));
%% ----------------------- Load mask (optional) -----------------------
% This section supports .png / .txt / .mat mask files.
% If mask is already labeled (1..N), we keep labels. If it's binary, we label with bwlabel.
mask = [];
if exist([mask_path, '.png'], 'file')
    mask_data = imread([mask_path, '.png']);
    if size(mask_data, 3) == 3
        mask_data = rgb2gray(mask_data);
    end
    mask_data = double(mask_data); % ensure numeric
    mask = bwlabel(mask_data > 0);
elseif exist([mask_path, '.txt'], 'file')
    mask_data = load([mask_path, '.txt']);
    if ~islogical(mask_data)
        mask = bwlabel(mask_data > 0);
    else
        mask = bwlabel(mask_data);
    end
elseif exist([mask_path, '.mat'], 'file')
    S = load([mask_path, '.mat']);
    fns = fieldnames(S);
    mask_candidate = S.(fns{1});
    if ~islogical(mask_candidate)
        mask = bwlabel(mask_candidate > 0);
    else
        mask = bwlabel(mask_candidate);
    end
else
    fprintf('Mask file not found. Proceeding without pre-defined mask.\n');
    mask = []; % will allow manual ROI selection
end
%% ---------------------- Create mask （optional）---------------------
% create mask via cellpose
cp = cellpose(Model="cyto2",ExecutionEnvironment="gpu");
img = uint8(reshape(mean(movie,2),nrows ,ncols));
labels = segmentCells2D(cp, img, ImageCellDiameter= 15, CellThreshold= 3);

bw  = labels > 0;
mask = bwlabel(bw);
num_rois = max(mask(:));
fprintf('Cellpose segmentation done, %d ROIs found.\n', num_rois);

% 7) 可视化更直观：叠加伪彩
figure; imshow(labeloverlay(img, mask)); title(sprintf('Cellpose mask (N=%d)', num_rois));

% 如需保存
imwrite(uint16(mask), fullfile(save_path, 'cellpose_mask.png'));

%% ----------------------- Select ROI -----------------------
% with or without Mask and Map
[bwmask, traces] = select_ROI(movie, ncols, nrows, mask, []); %#ok<ASGLU>
% traces is [nframes x Nroi]
num_rois = size(traces, 2);
fprintf('Finished ROI tracing, ROI number = %d\n', num_rois);

% Save raw trace plot + ROI mask returned by UI
fig_filename = fullfile(save_path, '1_raw_trace.fig');
png_filename = fullfile(save_path, '1_raw_trace.png');
roi_filename = fullfile(save_path, '1_raw_ROI.mat');
if ~isempty(get(0,'CurrentFigure'))
    saveas(gcf, fig_filename, 'fig');
    saveas(gcf, png_filename, 'png');
end
save(roi_filename, 'bwmask');

%% ----------------------- FITTING (single exponential) -----------------------
% Models:
%   on : y = A*(1 - exp(-k*c*x)) + bg
%   off: y = A*exp(-k*x) + bg

A_all  = zeros(num_rois,1);
k_all  = zeros(num_rois,1);
bg_all = zeros(num_rois,1);
RMSDs  = zeros(num_rois,1);

fit_x_data = cell(num_rois,1);
fit_y_data = cell(num_rois,1);
fit_y_fit  = cell(num_rois,1);

switch string(type)
    case "on"
        kinetics = fittype('A*(1 - exp(-k*c*x)) + bg', ...
            'independent','x','coefficients',{'A','k','bg'},'problem','c');
    case "off"
        kinetics = fittype('A*exp(-k*x) + bg', ...
            'independent','x','coefficients',{'A','k','bg'});
    otherwise
        error('type must be "on" or "off".');
end

parfor i = 1:num_rois
    fprintf('Fitting %d of %d ROIs\n', i, num_rois);

    x = t;                 % column
    y = traces(:, i);      % column

    % Initial guesses (robust-ish)
    y0   = mean(y(1:min(5,numel(y))));
    yend = mean(y(max(1,end-4):end));
    amp_guess = c;
    switch string(type)
        case "on"
        k_guess   = 10^5; % rough time-constant
        case 'off'
        k_guess   = 10^-5; % rough time-constant
    end
    bg_guess  = min(y0, yend);
        
    opts = fitoptions(kinetics);
    opts.MaxIter = 1e5; 
    opts.StartPoint = [amp_guess, k_guess, bg_guess];

    opts.Lower = [0, 0, 0];   % A>=0, k>=0, bg无界
    opts.Upper = [Inf, Inf, mean(img(:))]; % k限制在一个合理范围，比如 <10 (单位与时间轴相关)
    try
        if string(type) == "on"
            [mdl, gof] = fit(x, y, kinetics, opts, 'problem', c);
        else
            [mdl, gof] = fit(x, y, kinetics, opts);
        end
    catch ME
        warning('Fit failed for ROI %d: %s', i, ME.message);
        continue;
    end

    A_all(i)  = mdl.A;
    k_all(i)  = mdl.k;
    bg_all(i) = mdl.bg;
    RMSDs(i)  = gof.rmse;

    fit_x_data{i} = x;
    fit_y_data{i} = y;
    fit_y_fit{i}  = feval(mdl, x);
end

% Plot all fits
figure();
for i = 1:num_rois
    plot(fit_x_data{i}, fit_y_data{i}, 'b'); hold on;
    plot(fit_x_data{i}, fit_y_fit{i}, 'r');
end
hold off
fig_filename = fullfile(save_path, '2_fit_trace.fig');
png_filename = fullfile(save_path, '2_fit_trace.png');
saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');
fprintf('Finished fitting trace\n');

%% ----------------------- Pick good ROIs and normalize -----------------------
% 依赖：A_all(num_rois,1), traces(nframes,num_rois), t, save_path, num_rois

% 1) 计算比例指标
ymax   = max(double(traces), [], 1)';     % 每列最大值 (num_rois x 1)
ratioA = double(A_all(:)) ./ (ymax + eps); % 比例 (A / ymax)

% 2) 阈值选择 —— 你可以改成固定阈值，比如 ratioA > 0.5
use_percentile = false;

% 3) 可视化：比例分布 + 阈值
% —— 选出“最接近 1 的 5%”并在 log 直方图上画边界线 ——
ratioA = max(ratioA, eps);      % 保护，确保>0
logr   = log10(ratioA);           % log(A/maxY)，中心值=0 对应 ratio=1
dist   = abs(logr);             % 与0的距离（越小越接近1）

N      = numel(ratioA);
if max(abs(logr(:))) < log10(2)
    perc = 1;
    hisnum = 5;
else
    perc = 0.8;
    hisnum = 80;
end
nkeep  = max(1, round(perc * N));
[~, order] = sort(dist, 'ascend');
idx_keep   = order(1:nkeep);
keep       = false(N,1); 
keep(idx_keep) = true;
idx_drop = find(~keep);

% 5%内最大偏差作为边界
dmax        = max(dist(idx_keep));
log_thr_low = -dmax;
log_thr_high=  dmax;
thr_low     = exp(log_thr_low);
thr_high    = exp(log_thr_high);

fprintf('Top-%d%% closest to 1 → boundaries: thr_low=%.6g, thr_high=%.6g (log bounds: [%.4g, %.4g])\n', ...
        perc*100, thr_low, thr_high, log_thr_low, log_thr_high);

% 可视化
figure; histogram(logr, hisnum); grid on; box on;
xlabel('log(A / max(trace))'); ylabel('Count');
title(sprintf('Distribution of log(A/max(trace)) with %d%% closeness boundaries', perc*100));
yl = ylim;
line([log_thr_low  log_thr_low ], yl, 'Color','r','LineStyle','--','LineWidth',1.2);
line([log_thr_high log_thr_high], yl, 'Color','r','LineStyle','--','LineWidth',1.2);
legend({'Histogram',sprintf('%d%% closeness boundaries', perc*100)}, 'Location','best');

saveas(gcf, fullfile(save_path,'ratioA_hist.fig'));
saveas(gcf, fullfile(save_path,'ratioA_hist.png'));

% 4) 原始曲线可视化（上：剔除(抽样, 红)，下：保留(全量, 蓝)）
max_show = 150;  % 抽样上限，避免曲线过多
figure;

% ---- 上图：Dropped ----
subplot(2,1,1); hold on;
if ~isempty(idx_drop)
    show = round(linspace(1, numel(idx_drop), min(max_show, numel(idx_drop))));
    plot(t, double(traces(:, idx_drop(show))), 'Color', [0.85 0.33 0.10], 'LineWidth', 0.5); % 橙红
    lg1 = {'Dropped (sampled)'}; %#ok<NASGU>
else
    % 占位以免空图显得奇怪
    plot(t, zeros(size(t)), 'w');
    lg1 = {};
end
grid on; box on;
xlabel('Time (min)'); ylabel('Intensity (a.u.)');
title(sprintf('Traces — Dropped (n=%d)', numel(idx_drop)));
if ~isempty(lg1), legend(lg1, 'Location','best'); end

% ---- 下图：Kept ----
subplot(2,1,2); hold on;
if ~isempty(idx_keep)
    plot(t, double(traces(:, idx_keep)), 'Color', [0 0.447 0.741], 'LineWidth', 0.8); % 蓝
    lg2 = {'Kept'};
else
    plot(t, zeros(size(t)), 'w');
    lg2 = {};
end
grid on; box on;
xlabel('Time (min)'); ylabel('Intensity (a.u.)');
title(sprintf('Traces — Kept (n=%d)', numel(idx_keep)));
if ~isempty(lg2), legend(lg2, 'Location','best'); end

saveas(gcf, fullfile(save_path, 'traces_ratioA_keep_drop.fig'));
saveas(gcf, fullfile(save_path, 'traces_ratioA_keep_drop.png'));

%5) 导出结果（含比值/对数/距离/排名与边界）
% 计算对数及到1的log距离（用于完整记录）
log_ratioA = log(max(ratioA, eps));
dist_log   = abs(log_ratioA);             % 与0的距离 (ratio=1 对应 log=0)

% 给所有 ROI 做“接近 1”的排名（1=最接近）
[~, order_all] = sort(dist_log, 'ascend');
rank_closeness = nan(num_rois,1);
rank_closeness(order_all) = 1:num_rois;

roi_id    = (1:num_rois).';
kept_flag = keep(:);

% 将阈值也写入表（每行相同，便于留痕）
thr_low_col  = repmat(thr_low,  num_rois, 1);
thr_high_col = repmat(thr_high, num_rois, 1);

T_ratioA = table( ...
    roi_id, kept_flag, ...
    A_all(:), ymax, ratioA(:), log_ratioA(:), dist_log(:), rank_closeness, ...
    thr_low_col, thr_high_col, ...
    'VariableNames', { ...
        'ROI','Kept', ...
        'A','Ymax','A_div_maxY','log_A_div_maxY','log_dist_to_1','rank_by_closeness', ...
        'thr_low','thr_high' ...
    });

writetable(T_ratioA, fullfile(save_path, 'ROI_keep_by_ratioA.csv'));

% —— 收集合格 ROI 的索引与参数（供后续归一化/导出/平均使用）——
% 若还未定义 idx_drop，一并给出
idx_drop = setdiff(1:num_rois, idx_keep);

% 被保留的 ROI 索引
picked_idx   = idx_keep(:);

% 各种“picked_*”变量（统一为 double，避免整型混算）
picked_traces = double(traces(:, picked_idx));     % [nframes x nkeep]
picked_A      = double(A_all(picked_idx));         % [nkeep x 1]
picked_k      = double(k_all(picked_idx));         % [nkeep x 1]
picked_bg     = double(bg_all(picked_idx));        % [nkeep x 1]  % 若你没有逐ROI拟合bg，这里会是0
picked_RMSDs  = double(RMSDs(picked_idx));         % [nkeep x 1]

% 记录被保留的数量
picked_rois   = size(picked_traces, 2);


%% ===== Normalize kept ROI: (trace - bg) / A =====
% 依赖：keep（逻辑向量）、traces(nframes x num_rois)、A_all(num_rois,1)
% 可选：bg_all(num_rois,1)；t、save_path

idx_keep = find(keep);
if isempty(idx_keep)
    warning('No ROI kept by ratioA rule. Skip normalization.');
else
    % 1) 准备数据（转 double，避免整数混算错误）
    traces_d = double(traces);
    A_used   = double(A_all(:));

    % 2) 背景：优先用拟合得到的 bg_all；没有就用前10%帧的中位数估计
    if exist('bg_all','var') && ~isempty(bg_all) && numel(bg_all)==size(traces,2)
        bg_used = double(bg_all(:));
    else
        n = size(traces,1);
        w = max(1, round(0.1*n));                   % 前10%帧窗口
        bg_used = median(traces_d(1:w, :), 1)';      % (num_rois x 1)
    end

    % 3) 取 kept 的子集
    picked_traces = traces_d(:, idx_keep);          % (nframes x nkeep)
    picked_A      = A_used(idx_keep);               % (nkeep x 1)
    picked_bg     = bg_used(idx_keep);              % (nkeep x 1)

    % 4) 标准化：(trace - bg) / A
    % R2016b+ 支持隐式扩展；老版用 bsxfun
    try
        normalized_traces = (picked_traces - picked_bg.') ./ (picked_A.' + eps);
    catch
        normalized_traces = bsxfun(@rdivide, bsxfun(@minus, picked_traces, picked_bg.'), picked_A.' + eps);
    end

    % 5) 可视化 & 保存
    figure(); hold on;
    plot(t, normalized_traces, 'LineWidth', 0.8);
    title('Normalized Traces (kept only)'); 
    xlabel('Time (min)'); ylabel('Normalized Intensity'); 
    grid on; box on;
    fig_filename = fullfile(save_path, '4_normalized_trace.fig');
    png_filename = fullfile(save_path, '4_normalized_trace.png');
    saveas(gcf, fig_filename, 'fig');
    saveas(gcf, png_filename, 'png');
end

% ----------------------- Average normalized trace -----------------------
trace_std     = std(normalized_traces,0,2);
average_trace = mean(normalized_traces,2);
figure
plot(t, average_trace + trace_std,'LineWidth',1 ,'Color','blue'); hold on
plot(t, average_trace - trace_std,'LineWidth',1 ,'Color','blue'); hold on
plot(t, average_trace,            'LineWidth',2 ,'Color','red');  hold off
title('Average Curve');
xlabel('Time (min)'); ylabel('Normalized Intensity'); grid on;
fig_filename = fullfile(save_path, '5_average_trace.fig');
png_filename = fullfile(save_path, '5_average_trace.png');
saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');

%% ----------------------- Average fit curve band (mean ± std of k) -----------------------
avg_k = mean(picked_k);
std_k = std(picked_k);

switch string(type)
    case "on"
        average_fit_curve        = -exp(-avg_k * c * t) + 1;
        average_fit_curve_upper  = -exp(-(avg_k+std_k) * c * t) + 1;
        average_fit_curve_lower  = -exp(-(avg_k-std_k) * c * t) + 1;
    case "off"
        % For dissociation, normalized ideal form is exp(-k*t).
        average_fit_curve        =  exp(-avg_k * t);
        average_fit_curve_upper  =  exp(-(avg_k-std_k) * t);
        average_fit_curve_lower  =  exp(-(avg_k+std_k) * t);
end

figure;
plot(t, average_fit_curve,       'LineWidth', 2, 'Color', 'red'); hold on;
plot(t, average_fit_curve_upper, 'LineWidth', 1, 'Color', 'b');   hold on;
plot(t, average_fit_curve_lower, 'LineWidth', 1, 'Color', 'b');   hold off;
title('Average Fit Curve');
xlabel('Time (min)'); ylabel('Normalized Intensity'); grid on;
fig_filename = fullfile(save_path, '6_average_fit_curve.fig');
png_filename = fullfile(save_path, '6_average_fit_curve.png');
saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');


%% ----------------------- Save data (aligned names) -----------------------
% Per-ROI table: unified column names
roi_ids_col = (1:picked_rois).';  % ROI index in the kept set (1..picked_rois)
Type_col    = repmat(string(type), picked_rois, 1);

T = table(roi_ids_col, Type_col, picked_A(:), picked_k(:), picked_bg(:), picked_RMSDs(:), ...
    'VariableNames', {'ROI','Type','A','k','bg','RMSD'});

% Average table: unified column names; c only meaningful for "on"
avg_picked_k = mean(picked_k);
std_picked_k = std(picked_k);

fprintf('averaged k = %.8f \n', avg_picked_k );

if strcmpi(string(type),'on')
    c_col = c;   % scalar concentration for this dataset
else
    c_col = NaN; % not applicable for dissociation
end
T_avg = table(string(type), avg_picked_k, std_picked_k, c_col, ...
    'VariableNames', {'Type','Average_k','Std_k','c'});

% Normalized traces (kept only)
ROI_column_names = arrayfun(@(x) ['ROI' num2str(x)], 1:picked_rois, 'UniformOutput', false);
T_normalized = array2table([t(:), normalized_traces], 'VariableNames', [{'Time'}, ROI_column_names]);

% Averaged normalized trace
T_avg_trace  = array2table([t(:), average_trace, trace_std], 'VariableNames', {'Time','Average','Std'});

% Protocol / model metadata
if strcmpi(string(type),'on')
    kinetics_formula = 'A*(1 - exp(-k*c*t)) + bg';
    problemParamsStr = 'c';
else
    kinetics_formula = 'A*exp(-k*t) + bg';
    problemParamsStr = '';
end
indepVarStr = 't';
depVarStr   = 'Intensity';
coeffsStr   = 'A, k, bg';

T_protocol = table({kinetics_formula}, {indepVarStr}, {depVarStr}, {coeffsStr}, {problemParamsStr}, ...
                   'VariableNames', {'Formula','Independent_Variables','Dependent_Variables','Coefficients','Constants'});

% Write all sheets
table_filename = fullfile(save_path, '7_Kinetics_params.xlsx');
writetable(T,            table_filename, 'Sheet','Fitting results');          % per-ROI with unified names
writetable(T_avg,        table_filename, 'Sheet','Average Fitting results');  % unified avg table
writetable(T_protocol,   table_filename, 'Sheet','Average Fitting results','Range','A6');
writetable(T_normalized, table_filename, 'Sheet','Normalized traces');
writetable(T_avg_trace,  table_filename, 'Sheet','Averaged traces');

fprintf('Saved kinetics table (aligned column names: A, k, bg, RMSD)\n');


%% ----------------------- Save code bundle -----------------------
code_path = fullfile(save_path,'Code');
if ~exist(code_path, 'dir'); mkdir(code_path); end
thisFile = mfilename('fullpath');
try
    [requiredFiles, ~] = matlab.codetools.requiredFilesAndProducts([thisFile '.m']);
    % Include this script itself if not captured
    requiredFiles = unique([requiredFiles; {[thisFile '.m']}]');
    for kf = 1:numel(requiredFiles)
        [~, name, ext] = fileparts(requiredFiles{kf});
        copyfile(requiredFiles{kf}, fullfile(code_path, [name, ext]));
    end
    fprintf('All files have been copied to %s\n', code_path);
catch ME
    warning('Could not enumerate/copy required files: %s', ME.message);
end
