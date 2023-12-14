% Calcium_ANALYSIS POWERED BY LIU-YANG LUORONG.
% version 1
% https://github.com/Luorong37/AP_analysis

% NEED FUNCTIONS:


% NEED TOOLBOX:
% Image Processing Toolbox, Curve Fitting Toolbox, Signal Processing Toolbox

% DO NOT CLEAN THE VAR FROM AP_ANALYSIS
%% Loading raw data
t1 = tic; % Start a timer
nowtime = string(datetime('now'));
% Replace colons with hyphens to get the desired output format
nowtime = strrep(nowtime , ':', '-');
fprintf('Loading...\n')

% ↓↓↓↓↓-----------Prompt user for define path-----------↓↓↓↓↓
% support for folder, .tif, .tiff, .bin.
folder_path = [folder_path '_Green']; % defined name. please change if not
file_name = '\';  % must add format.
% ↓↓↓↓↓-----------Prompt user for frame rate------------↓↓↓↓↓
freq = 10; % Hz (defined = 10)
% -----------------------------------------------------------

% read path
file_path = fullfile(folder_path, file_name);
split_path = split(file_name, '.');
if length(split_path)>1
    % when read a file
    file_extension = string(split_path(end));
    save_path =  fullfile(folder_path, [cell2mat(split_path(1)),'_Analysis'], nowtime);
else
    % when read a folder
    file_extension = 'tif';
    save_path = fullfile([folder_path, '_Analysis'] ,nowtime);
end
mkdir(save_path);

% Load image file
[movie, ncols, nrows, nframes] = load_movie(file_path,file_extension);

% Define parameters
dt = 1 / freq; % Calculate time axis
colors = [lines(7);hsv(5);spring(3);winter(3);gray(3)];
t = (1:nframes) * dt;
mask = [];

t2 = toc(t1); % Get the elapsed time
fprintf('Finished loading movie after %d s\n',round(t2))
%% Load saved ROI from AP_analysis

% ROI from AP_analysis. if the var has been clean, please
mask_filename = roi_filename;
mask = load(mask_filename);
mask = mask.rois;

fprintf('Finished mask loading\n')
%% Photobleaching correction
traces_input = zeros(size(traces,1),size(traces,2)-1);
background = traces(:,end);

% remove background
for i = 1: size(traces,2)-1
    traces_input(:,i) = traces(:,i) - background;
end

% fit
[traces_corrected, fitted_curves] = fit_exp1(traces_input, freq);

% plot
fig = figure();
set(fig,'Position',get(0,'Screensize'));
fit_axe = subplot(2,1,1);
fited_axe = subplot(2,1,2);

for i = 1: size(traces_corrected,2)
    plot(t,traces_input(:,i),'Color',colors(i,:),'Parent',fit_axe);
    hold(fit_axe, 'on');
    plot(t,fitted_curves(:,i),'Color',colors(i,:),'LineWidth',2,'Parent',fit_axe);
    hold(fit_axe, 'on');
    plot(t,traces_corrected(:,i),'Color',colors(i,:),'Parent',fited_axe);
    hold(fited_axe, 'on');
end

% note
title(fit_axe, 'Original and Fitted Curves');
title(fited_axe, 'Corrected Traces');
legend(fit_axe, 'Original Trace', 'Fitted Curve');

fig_filename = fullfile(save_path, '2_fitted_trace.fig');
png_filename = fullfile(save_path, '2_fitted_trace.png');

saveas(gcf, fig_filename, 'fig');
saveas(gcf, png_filename, 'png');

fprintf('Finished exp1 fit\n')
