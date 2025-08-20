raw = traces_corrected(1:24000,1);
% raw = traces_corrected(:,27);
de_new = wdenoise(raw,10,DenoisingMethod="FDR");
de = wdenoise(raw,10);
subplot(3,1,1)
plot(raw);
title('raw')
subplot(3,1,2)
plot(de)
title('wavelet denoising SURE')
subplot(3,1,3)
plot(de_new)
title('wavelet denoising FDR')

%%
figure
cwt(raw,'amor',400);
figure
cwt(de_FDR,'amor',400);

%%
subplot(5,1,1)
plot(raw);
title('raw')
subplot(5,1,2)
plot(de)
title('wavelet denoising in bayes')
subplot(5,1,3)
plot(raw-de)
title('raw-denoised = noise')
subplot(5,1,4)
plot(calculate_SNR(raw))
title('SNR previous')
subplot(5,1,5)
plot(raw/std(raw-de))
title('SNR wavelet')

%%
subplot(4,1,1)
plot(raw);hold on;
title('denoising')
plot(de)
legend('raw')
legend('wavelet denoising in SURE')
subplot(4,1,2)
plot(raw);hold on;
plot(raw-de)
title('noise')
legend('raw')
legend('noise')
subplot(4,1,3)
plot(calculate_SNR(raw));
title('previous SNR calculation')
subplot(4,1,4)
plot(raw/std(raw-de));hold on;
title('new SNR calculation using wavelet')

%%
% 初始化变量
MinPeakProminence_factor = 0.5; % 定义最小峰值显著性因子
peaks_polarity = zeros(1,nrois);
peaks_threshold = zeros(1,nrois);
peaks_amplitude = cell(1, nrois);
peaks_index = cell(1, nrois);
peaks_sensitivity = cell(1, nrois);


fig = figure;
set(gcf,'Position',get(0,'Screensize'));
i = 0;


% plot trace;
subplot(2,1,1)
current_trace = raw*-1;
plot_trace = current_trace;
plot(plot_trace); hold on;
title('raw, without manual threshold setting')
% 寻找峰值
MinPeakProminence = (max(plot_trace)-mean(plot_trace)) * MinPeakProminence_factor;
[peak_y, peak_x] = findpeaks(plot_trace, 'MinPeakProminence', ...
    MinPeakProminence);
% plot peak
plot(peak_x, (peak_y),'v','MarkerFaceColor','r');
hold off;

% plot detrace;
subplot(2,1,2)
plot_trace = de*-1;
plot(plot_trace); hold on;
title('denoised == SURE, without manual threshold setting')
% 寻找峰值
MinPeakProminence = (max(plot_trace)-mean(plot_trace)) * MinPeakProminence_factor;
[peak_y, peak_x] = findpeaks(plot_trace, 'MinPeakProminence', ...
    MinPeakProminence);
% plot peak
plot(peak_x, (peak_y),'v','MarkerFaceColor','r');
hold off;

