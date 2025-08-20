function [traces_corrected, baseline] = highpass_bleach_remove(traces, Fs, FcL)
% highPASS_DENOISING 对输入的电压信号矩阵进行高通滤波漂白
% 
% 输入参数：
%   traces - 一个矩阵，每一列是一个独立的电压信号序列
%   Fs - 采样频率 (Hz)，可选，默认值为500 Hz
%   FcL - 低通滤波的截止频率 (Hz)，可选，默认值为50 Hz
%
% 输出参数：
%   trace_dn - 去噪后的电压信号矩阵，每一列对应输入的一个信号

% 默认采样频率和截止频率
if nargin < 2 || isempty(Fs)
    Fs = 400; % 默认采样频率为400 Hz
end
if nargin < 3 || isempty(FcL)
    FcL = 0.1; % 默认0.1 Hz
end

% 设置填充长度，将信号在头尾各填充5%的长度
padlength = round(0.05 * size(traces, 1));

% 对信号进行填充，以减少边界效应
padded_traces = [repmat(traces(1, :), padlength, 1); ...
                 traces; ...
                 repmat(traces(end, :), padlength, 1)];

% 初始化存储滤波后信号的矩阵，大小与填充后的信号相同
filtered_traces = zeros(size(padded_traces));

% 对每一列（每一个独立信号）进行低通滤波
for col = 1:size(traces, 2)
    % disp(['Processing ROI',num2str(col)]);
    % 对填充后的信号应用高通滤波器
    filtered_traces(:, col) = highpass(padded_traces(:, col), FcL, Fs);
end

% 去掉填充部分，保留原始长度的滤波后信号
corrected_traces = filtered_traces(padlength+1:end-padlength, :);

% 输出去噪后的信号矩阵
traces_corrected = corrected_traces;
baseline = smoothdata(traces - traces_corrected);

end
