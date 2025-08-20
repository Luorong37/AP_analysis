function [peaks_index, peaks_amplitude, peaks_polarity, parts_results] = ...
    peak_finding_auto(traces, save_path, varargin)
%PEAK_FINDING_AUTO 自动检测信号中的峰值并保存可视化结果
%   本函数通过分块处理和多参数优化，自动识别神经钙信号中的动作电位峰值。
%   支持灵活的参数调节，并生成峰值检测的可视化报告。
%
%   Syntax语法:
%   [peaks_index, peaks_amplitude, peaks_polarity, peaks_sensitivity] = ...
%       peak_finding_auto(traces, save_path)
%   [peaks_index, peaks_amplitude, peaks_polarity, peaks_sensitivity] = ...
%       peak_finding_auto(traces, save_path, Name, Value)
%
% INPUT ARGUMENTS输入参数:
%   traces        - 待分析信号矩阵 [N x M]，N=时间点，M=ROI编号 (必需, numeric)
%   save_path     - 结果保存路径 (必需, string)
%
% NAME-VALUE PAIR ARGUMENTS名称-值对参数:
%   'parts'                   - 信号分段数量 (默认: 1, integer ≥1)
%   'MinPeakProminence_factor'- 最小峰值显著性因子 (默认: 0.3, scalar >0)
%   'MinPeakDistance_factor'  - 最小峰值间距因子 (默认: 5, scalar >0)
%   'MinPeakHeight'           - 绝对最小峰值高度 (默认: 0.03, scalar)
%
% OUTPUT ARGUMENTS输出参数:
%   peaks_index          - 峰值时间索引 [cell array]
%   peaks_amplitude      - 峰值幅度 [cell array]
%   peaks_polarity       - 峰值极性矩阵 [1:正峰, -1:负峰, 0:未检测]
%   peaks_sensitivity    - 峰值检测灵敏度指标 [cell array]
%
% EXAMPLES示例:
%   % 基础用法
%   [idx, amp] = peak_finding_auto(trace_data, 'D:\results');
%
%   % 高级参数调节
%   [idx, amp] = peak_finding_auto(trace_data, 'D:\results', ...
%       'parts', 3, ...
%       'MinPeakHeight', 0.05, ...
%       'MinPeakDistance_factor', 10);
%
% REFERENCES参考:
%   峰值检测算法基于MATLAB内置的 findpeaks 函数实现
%   Signal Processing Toolbox required
%
% See also findpeaks, inputParser, varargin

% 参数解析模块
% 1. 创建 inputParser 对象
p = inputParser;

% 2. 添加必需参数（traces, parts, data_path）
addRequired(p, 'traces', @isnumeric);     % 验证 traces 为数值
addRequired(p, 'data_path');     % 验证 data_path 为字符串

% 3. 添加可选参数（名称-值对，带默认值和类型验证）
addParameter(p, 'parts', 1, @isnumeric);      % 验证 parts 为数值
addParameter(p, 'MinPeakProminence_factor', 0.3, @isnumeric);
addParameter(p, 'MinPeakDistance_factor', 5, @isnumeric);
addParameter(p, 'MinPeakHeight', 0.03, @isnumeric);
addParameter(p, 'RawTraces', []);

% 4. 解析输入参数
parse(p, traces, save_path, varargin{:});

% 5. 提取参数值
MinPeakProminence_factor = p.Results.MinPeakProminence_factor;
MinPeakDistance_factor = p.Results.MinPeakDistance_factor;
MinPeakHeight = p.Results.MinPeakHeight;
parts = p.Results.parts;
if ~isempty(p.Results.RawTraces)
    raw = p.Results.RawTraces;
else
    raw = [];
end


nrois = size(traces,2);
frames = length(traces);
peaks_polarity_part = zeros(parts,nrois);
peaks_amplitude_part = cell(parts, nrois);
peaks_index_part = cell(parts, nrois);
% peaks_sensitivity_part = cell(parts, nrois);
% calculate part index
step_frame = round(frames/parts) ;
partindex = {};
end_frame = 0;
for i = 1:parts
    start_frame = 1 + end_frame;
    end_frame = round(min(start_frame + step_frame -1 ,frames));
    partindex{i} = start_frame:end_frame;
end

pffolder = sprintf('Auto denoised %.d p peakfinding ,Prominence = %.2f, Distance = %.2f, Height = %.2f' ,parts,MinPeakProminence_factor, MinPeakDistance_factor,MinPeakHeight);
mkdir (fullfile(save_path,pffolder))
%peakfinding
for p = 1:parts
    fprintf('Processing part %d',p)
    for i = 1:nrois


        current_trace = traces(partindex{p},i);

        % 确定峰值极性
        % abs(min(current_trace) - mean(current_trace)) < max(abs(current_trace) - mean(current_trace))
        if abs(min(current_trace) - mean(current_trace)) < abs(max(current_trace) - mean(current_trace))
            peaks_polarity_part(p,i) = 1;
            polarity = 'Positive';        
        else
            peaks_polarity_part(p,i) = -1;
            polarity = 'Negative';
        end
        plot_trace = current_trace * peaks_polarity_part(p,i) ;

        % 寻找峰值
        % MinPeakProminence = (max(current_trace)-min(current_trace)) * MinPeakProminence_factor;
        maxpeaksheight = (max(current_trace)-min(current_trace)) * MinPeakProminence_factor;
        MinPeakProminence = max(maxpeaksheight,MinPeakHeight);
        [peak_y, peak_x] = findpeaks(plot_trace, 'MinPeakProminence', ...
            MinPeakProminence,'MinPeakDistance',MinPeakDistance_factor,'MinPeakHeight',MinPeakHeight);
        
        peaks_amplitude_part{p,i}  = peak_y;
        peaks_index_part{p,i} = peak_x ;
        % peaks_sensitivity_part{p,i} = current_trace(peak_x)-mean(current_trace);
        peaks_amplitude_part{p,i} = peaks_amplitude_part{p,i} * peaks_polarity_part(p,i);


        figure()
        plot(plot_trace); hold on;
        if numel(peak_x) ~= 0
            plot(peak_x, (peaks_amplitude_part{p,i}* peaks_polarity_part(p,i)),'v','MarkerFaceColor','r');
        end
        title(sprintf('peakfinding of p = %d, roi %d.png',p,i));
        hold off;
        saveas(gcf,fullfile(save_path,pffolder,sprintf('peakfinding of roi %d, p = %d.png',i,p)));
        saveas(gcf,fullfile(save_path,pffolder,sprintf('peakfinding of roi %d, p = %d.fig',i,p)));
        close(gcf);

        if ~isempty(raw)
            current_raw = raw(partindex{p},i);
            plot_trace = current_raw * peaks_polarity_part(p,i) ;
            figure()
            plot(plot_trace); hold on;
            if numel(peak_x) ~= 0
                plot(peak_x, (current_raw(peak_x)* peaks_polarity_part(p,i)),'v','MarkerFaceColor','r');
            end
            title(sprintf('raw peakfinding of p = %d, roi %d.png',p,i));
            hold off;
            saveas(gcf,fullfile(save_path,pffolder,sprintf('raw peakfinding of roi %d, p = %d.png',i,p)));
            saveas(gcf,fullfile(save_path,pffolder,sprintf('raw peakfinding of roi %d, p = %d.fig',i,p)));
            close(gcf);

        end



        if numel(peak_x) == 0
            peaks_polarity_part(p,i) = 0;
        end
        fprintf('%d peaks found in ROI %d, Polarity : %s\n',numel(peak_x),i,polarity)
    end
    fprintf('All peaks found.\n')
end

peaks_index= cell(1,nrois);
for i = 1:parts
    for j = 1:nrois
        offset_index = peaks_index_part{i,j} + partindex{i}(1) -1;
        peaks_index{j} = [peaks_index{j} ;offset_index];
    end
end

peaks_amplitude = cell(1,nrois);
for i = 1:parts
    for j = 1:nrois
        peaks_amplitude{j} = [peaks_amplitude{j} ;peaks_amplitude_part{i,j}];
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
        count1 = sum(peaks_polarity_part(:,i) == 1);
        countNeg1 = sum(peaks_polarity_part(:,i) == -1);
    
    % 比较并返回结果
    if count1 > countNeg1
        peaks_polarity{i} = 1;
    elseif countNeg1 > count1
       peaks_polarity{i} = -1;
    else
        peaks_polarity{i} = -1; % 数量相等时默认返回1
    end
    % polarity_index = find(abs(peaks_polarity_part(:,i)) == max(abs(peaks_polarity_part(:,i))));
    % peaks_polarity{i} =peaks_polarity_part(polarity_index(1),i);
end

parts_results.index = partindex;
parts_results.peaks_polarity = peaks_polarity_part;
parts_results.peaks_amplitude = peaks_amplitude_part;
parts_results.peaks_index = peaks_index_part;



save(fullfile(save_path,pffolder,'peakfinding_denoised_version.mat'), 'MinPeakProminence_factor', 'MinPeakDistance_factor', 'parts', ...
    'peaks_polarity_part', 'peaks_index_part', 'peaks_amplitude_part', 'MinPeakHeight','peaks_index','peaks_polarity','peaks_amplitude','traces');
end