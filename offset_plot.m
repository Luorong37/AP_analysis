function [offset_array] = offset_plot(traces_corrected,t)

colors = lines(100);
offset_array = zeros(size(traces_corrected,2));
legendlist = {};
for i = 1:size(traces_corrected,2)
    
    % plot trace
    % if polarity(i) == -1
    trace = traces_corrected(:,i);
    if i == 1
        % offset_current = - min(trace);
        offset_current = 0;
        offset = offset_current;
        offset_array(i) = offset;
    else
        offset_current = mean(traces_corrected(:,i)) - min(traces_corrected(:,i)) + max(traces_corrected(:,i-1)) - mean(traces_corrected(:,i-1));
        % offset_current = max(traces(:,i-1)) - min(traces(:,i-1)) + mean(trace)- 2 * min(trace) ;
        % offset = offset_array(i-1) + min(traces(:,i-1)) + offset_current;
        offset = offset_current;
        offset_array(i) = offset + offset_array(i-1);
    end
    
    plot(t,  traces_corrected(:,i) + offset_array(i),'Color',colors(i,:)); hold on;
    legendlist{end+1} = sprintf('ROI %d',i);
    % elseif polarity(i) == 1
    %     trace = (traces(:,i) - mean(traces(:,i))) + offset;
    %     plot(t,  trace ,'Color',colors(i,:)); hold on;
    %     offset_current = (max(trace(:))-min(trace(:)));
    %     offset = offset + offset_current;
end
legend(legendlist);hold on;
end
% 创建stackedplot


% traces_flipped =  traces .* polarity + 2 - polarity;
% sp = stackedplot(t, traces_flipped);
% 
% % 统一y轴刻度范围
% minY = min(traces_flipped,[],'all'); % 统一的最小值
% maxY = max(traces_flipped,[],'all'); % 统一的最大值
% 
% % 获取所有的AxesProperties
% allAxesProperties = sp.AxesProperties;
% 
% % 设置每个轴的y轴范围
% for k = 1:length(allAxesProperties)
%     allAxesProperties(k).YLimits = [minY maxY];
% end
