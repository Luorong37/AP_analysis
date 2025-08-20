%% subthreshold analysis

dire = -1;
pathname = 'E:\1_Data\JQH_wanglab\20250520\cell12';
% 第一次阈值选择：检测所有峰
figure();
findpeaks(dire*(trace), ...
    'MinPeakDistance', 10,'MinPeakProminence',0.02);
hold on;
title('左键单击选择所有峰的阈值');
[~, threshold_all_optical] = ginput(1);  % 注意变量名与电信号区分
plot([1 length(trace)], [threshold_all_optical threshold_all_optical], 'g', 'LineWidth', 2);
legend('原始峰', '所有峰阈值');

% 检测所有峰
[~, spikeT_all_optical] = findpeaks(...
    dire*trace, ...
    'MinPeakHeight', threshold_all_optical, ...
    'MinPeakDistance', 20);


% 第二次阈值选择：检测阈上峰
title('左键单击选择阈上峰的阈值');
[~, threshold_supra_optical] = ginput(1); 
plot([1 length(trace)], [threshold_supra_optical threshold_supra_optical], 'r', 'LineWidth', 2);
legend('原始峰', '所有峰阈值', '阈上峰阈值');

% 检测阈上峰
[~, spikeT_supra_optical] = findpeaks(...  % 注意变量名与电信号区分
    dire*(trace), ...
    'MinPeakHeight', threshold_supra_optical, ...
    'MinPeakDistance',20);


saveas(gca, [pathname '6 光信号峰检测.fig']);
saveas(gca, [pathname '6 光信号峰检测.png']);

% 计算阈下峰
spikeT_sub_optical = setdiff(spikeT_all_optical, spikeT_supra_optical);

%% 处理光信号三类峰

AP_width = 80;

% 所有峰
nspike_all_optical = length(spikeT_all_optical);
kernel_stack_all_optical = zeros(2*AP_width+1, nspike_all_optical);
amplitude_all_optical = zeros(1, nspike_all_optical);

for n = 1:nspike_all_optical
    kernel_stack_all_optical(:,n) = trace(spikeT_all_optical(n)-AP_width:spikeT_all_optical(n)+AP_width);
    amplitude_all_optical(n) = min(kernel_stack_all_optical(:,n)) - mean(kernel_stack_all_optical(1:10,n));
end

% 阈上峰
nspike_supra_optical = length(spikeT_supra_optical);  % 确保变量名唯一性
kernel_stack_supra_optical = zeros(2*AP_width+1, nspike_supra_optical);
amplitude_supra_optical = zeros(1, nspike_supra_optical);

for n = 1:nspike_supra_optical
    kernel_stack_supra_optical(:,n) = trace(spikeT_supra_optical(n)-AP_width:spikeT_supra_optical(n)+AP_width);
    amplitude_supra_optical(n) = min(kernel_stack_supra_optical(:,n)) - mean(kernel_stack_supra_optical(1:10,n));
end

% 阈下峰
nspike_sub_optical = length(spikeT_sub_optical);
kernel_stack_sub_optical = zeros(2*AP_width+1, nspike_sub_optical);
amplitude_sub_optical = zeros(1, nspike_sub_optical);

for n = 1:nspike_sub_optical
    kernel_stack_sub_optical(:,n) = trace(spikeT_sub_optical(n)-AP_width:spikeT_sub_optical(n)+AP_width);
    amplitude_sub_optical(n) = min(kernel_stack_sub_optical(:,n)) - mean(kernel_stack_sub_optical(1:10,n));
end

%%可视化光信号分类结果
figure();
subplot(3,1,1);
plot(dire*kernel_stack_all_optical, 'Color', [0.8 0.8 0.8]);
hold on;
plot(dire*mean(kernel_stack_all_optical,2), 'k', 'LineWidth', 2);
title(['所有光信号峰 (n=' num2str(nspike_all_optical) ')']);
box off; axis tight;

subplot(3,1,2);
plot(dire*kernel_stack_supra_optical, 'Color', [0.8 0.8 0.8]);
hold on;
plot(dire*mean(kernel_stack_supra_optical,2), 'r', 'LineWidth', 2);
title(['阈上光信号峰 (n=' num2str(nspike_supra_optical) ')']);
box off; axis tight;

subplot(3,1,3);
plot(dire*kernel_stack_sub_optical, 'Color', [0.8 0.8 0.8]);
hold on;
plot(dire*mean(kernel_stack_sub_optical,2), 'b', 'LineWidth', 2);
title(['阈下光信号峰 (n=' num2str(nspike_sub_optical) ')']);
box off; axis tight;

saveas(gca, [pathname '7 光信号峰分类对比.fig']);
saveas(gca, [pathname '7 光信号峰分类对比.png']);

%%保存光信号结果到Excel
% 定义列名
optical_headers = {'峰编号', '幅度', '基线', '峰值', '类型'};

% 合并数据
optical_data = [
    table((1:nspike_all_optical)', amplitude_all_optical', ...
    mean(kernel_stack_all_optical(1:10,:))', max(kernel_stack_all_optical)', ...
    repmat({'All'}, nspike_all_optical, 1), 'VariableNames', optical_headers);
    
    table((1:nspike_supra_optical)', amplitude_supra_optical', ...
    mean(kernel_stack_supra_optical(1:10,:))', max(kernel_stack_supra_optical)', ...
    repmat({'Supra'}, nspike_supra_optical, 1), 'VariableNames', optical_headers);
    
    table((1:nspike_sub_optical)', amplitude_sub_optical', ...
    mean(kernel_stack_sub_optical(1:10,:))', max(kernel_stack_sub_optical)', ...
    repmat({'Sub'}, nspike_sub_optical, 1), 'VariableNames', optical_headers)
];

writetable(optical_data, [pathname 'optical_analysis.xlsx'], 'Sheet', '原始数据');

% 保存统计结果
stats_table = table(...
    {'All'; 'Supra'; 'Sub'}, ...
    [nspike_all_optical; nspike_supra_optical; nspike_sub_optical], ...
    [mean(amplitude_all_optical); mean(amplitude_supra_optical); mean(amplitude_sub_optical)], ...
    [std(amplitude_all_optical)/sqrt(nspike_all_optical); ...
     std(amplitude_supra_optical)/sqrt(nspike_supra_optical); ...
     std(amplitude_sub_optical)/sqrt(nspike_sub_optical)], ...
    'VariableNames', {'类型', '数量', '平均幅度', '标准误'});

writetable(stats_table, [pathname 'optical_analysis.xlsx'], 'Sheet', '统计分析');