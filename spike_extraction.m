%% 
load('D:\1_Data\2b. Dual-color imaging in SCN\2024.09.05_P2A-G8s\20240905-153142POA_Analysis\2024-10-05 12-01-08\4_Signal analysis.mat')% 假设信号数据已经存在 'trace' 变量中
trace = traces_corrected(:,8);
%% 
%% Photobleaching Correction
% Apply a bi-exponential fit to the raw trace to assess photobleaching
% This example uses a basic exponential fitting function, which may need adjustment based on your data.
time = t; % Use the previously defined time vector
[fitresult, ~] = fit(time', trace, 'exp2'); % Bi-exponential fitting
fit_trace = fitresult(time); % Evaluate the fit at each time point

% Normalize the trace by the fit function to correct for photobleaching
trace_corrected = trace ./ fit_trace;

% Apply a high-pass filter to remove low-frequency drift (0.5 Hz)
trace_corrected_hp = highpass(trace_corrected, 0.5, fs);

% Convert the trace to % ΔF/F0 with F0 as the mean signal
F0 = mean(trace_corrected_hp);
trace_dFF = (trace_corrected_hp - F0) / F0 * 100;

%% Spike Detection
% Initialize an empty array for detected spikes
spikes = zeros(size(trace_dFF));

% Apply the three metrics for spike detection

% Metric 1: High-pass filter (40 Hz)
trace_highpass = highpass(trace_dFF, 40, fs);

% Metric 2: Cumulative probability transform with erf
trace_transformed = zeros(size(trace_highpass));
for i = 1:length(trace_highpass)-N+1
    segment = trace_highpass(i:i+N-1);
    trace_transformed(i+N-1) = erf(mean(segment)); 
end

% Metric 3: Cumulative product of wavelet scales (2nd to 5th)
[~, L] = wavedec(trace_highpass, level, waveletType);
wavelet_product = ones(size(trace_highpass));
for i = 2:5
    detail = detcoef(C, L, i); 
    wavelet_product = wavelet_product .* upcoef('d', detail, waveletType, i, length(trace_highpass));
end

% Align wavelet product to original trace length
wavelet_product = wavelet_product(1:length(trace_highpass));

% Normalize and detect spikes with Z score > 3 in all three metrics
z_score_threshold = 3;
z_trace_highpass = zscore(trace_highpass);
z_trace_transformed = zscore(trace_transformed);
z_wavelet_product = zscore(wavelet_product);

spikes(z_trace_highpass > z_score_threshold & z_trace_transformed > z_score_threshold & z_wavelet_product > z_score_threshold) = 1;

%% Extract Spike Waveform Metrics
% Average waveform alignment on spike onset
spike_indices = find(spikes); % Indices of detected spikes
spike_waveforms = zeros(length(spike_indices), round(0.02 * fs) + 1); % 20 ms window around each spike

for j = 1:length(spike_indices)
    if spike_indices(j) + round(0.01 * fs) <= length(trace_dFF)
        spike_waveforms(j, :) = trace_dFF(spike_indices(j)-round(0.01*fs) : spike_indices(j)+round(0.01*fs));
    end
end
average_spike_waveform = mean(spike_waveforms, 1);

% Spike amplitude and FWHM
spike_amplitude = max(average_spike_waveform); % Peak amplitude from onset
half_max = spike_amplitude / 2;

% Find the indices where the average waveform is above half max
fwhm_idx = find(average_spike_waveform >= half_max);

% Check if fwhm_idx is not empty to calculate FWHM
if ~isempty(fwhm_idx)
    fwhm_time = (fwhm_idx(end) - fwhm_idx(1)) / fs; % FWHM in seconds
else
    fwhm_time = NaN; % Assign NaN if FWHM cannot be calculated
    warning('FWHM calculation failed: No values meet the half-maximum amplitude threshold.');
end

%% Subthreshold Fluctuation Analysis
% Bandpass filter for subthreshold fluctuations (0.1 - 30 Hz)
trace_subthreshold = bandpass(trace_dFF, [0.1 30], fs);

% Define cell up and down states
cell_up = mean(trace_subthreshold(spike_indices - round(0.005 * fs) : spike_indices + round(0.005 * fs)));
cell_down = prctile(trace_subthreshold, 1);

% Calculate subthreshold fluctuation
subthreshold_fluctuation = cell_up - cell_down;

%% Display Results
% Plot results
figure;
subplot(3, 1, 1);
plot(t, trace_dFF);
title('Photobleaching-Corrected Signal (% ΔF/F0)');
xlabel('Time (s)');
ylabel('% ΔF/F0');

subplot(3, 1, 2);
plot(t, trace_dFF, 'k'); hold on;
plot(t(spike_indices), trace_dFF(spike_indices), 'ro');
title('Detected Spikes');
xlabel('Time (s)');
ylabel('% ΔF/F0');

subplot(3, 1, 3);
plot(linspace(-10, 10, size(average_spike_waveform, 2)), average_spike_waveform);
title(['Average Spike Waveform (FWHM = ', num2str(fwhm_time*1000), ' ms)']);
xlabel('Time around Spike (ms)');
ylabel('% ΔF/F0');

disp(['Subthreshold Fluctuation: ', num2str(subthreshold_fluctuation)]);
disp(['Spike Amplitude: ', num2str(spike_amplitude)]);
disp(['FWHM of Spike: ', num2str(fwhm_time * 1000), ' ms']);



%% 
%% 
%% 
%% 
%% 
%% 
%% 
%% 
% 'trace' 是一个包含电压信号的数据向量
% 假设采样频率为 fs (单位：Hz)
fs = 400; % 示例采样频率，1kHz
t = (0:length(trace)-1) / fs; % 计算时间向量
trace_dn = lowpass_denoising(trace);


%% 第一步：高通滤波（40 Hz）
% 使用二阶Butterworth高通滤波器进行滤波
low_cutoff = 40; % 40 Hz 作为高通滤波的截止频率

% 使用 highpass 函数进行高通滤波
padlength = round(0.05 * length(trace));

% 对信号进行填充，以减少高通滤波的边界效应
padded_trace = [repmat(trace(1), padlength, 1); trace; repmat(trace(end), padlength, 1)];

% 使用highpass函数进行高通滤波
trace_padded_filtered = highpass(padded_trace, low_cutoff, fs);

% 去除填充部分，保留原始长度的滤波后信号
trace_filtered = trace_padded_filtered(padlength+1:end-padlength);

% 绘制滤波前后的信号
figure;
subplot(2, 1, 1);
plot(t, trace);
title('原始电压信号');
xlabel('时间（秒）');
ylabel('电压（V）');

subplot(2, 1, 2);
plot(t, trace_filtered);
title('高通滤波后的信号');
xlabel('时间（秒）');
ylabel('电压（V）');

%% 第二步：累计概率变换（使用误差函数）
% 使用标准的误差函数 (erf) 来进行变换
duration = 1.6e-3; % 持续时间 1.6 ms
N = round(duration * fs); % 计算需要的样本数
trace_transformed = zeros(size(trace_filtered));

% 对每个信号片段进行累计概率变换
for i = 1:length(trace_filtered)-N+1
    segment = trace_filtered(i:i+N-1);
    trace_transformed(i+N-1) = erf(mean(segment)); % 使用误差函数变换均值
end

% 绘制变换后的信号
figure;
plot(t, trace_transformed);
title('累计概率变换后的信号');
xlabel('时间（秒）');
ylabel('变换值');

%% 第三步：小波变换（Coiflet小波）
level = 5; % 分解层数
waveletType = 'coif1'; % 小波类型

% 对信号进行小波分解
[C, L] = wavedec(trace, level, waveletType);

% 提取每层的细节系数
cD = cell(1, level);
for i = 1:level
    % 确保L传递的是小波分解返回的长度数组
    cD{i} = detcoef(C, L, i); % 提取第i层的细节系数
end

% 绘制每层的细节系数
figure;
for i = 1:level
    subplot(level, 1, i);
    plot(cD{i});
    title(['第 ', num2str(i), ' 层的细节系数']);
    xlabel('样本点');
    ylabel('系数值');
end
%% 最后，合并结果并显示
% 你可以根据这些步骤调整和进一步优化检测算法，进一步提取峰值
