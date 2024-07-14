function [offset_array] = offset_plot(traces,t)



colors = lines(100);
offset_array = zeros(size(traces,2));

for i = 1:size(traces,2)
    
    % plot trace
    % if polarity(i) == -1
    trace = traces(:,i);
    if i == 1
        offset_current = - min(trace);
        offset = offset_current;
    else
        offset_current = max(traces(:,i-1)) - min(traces(:,i-1)) + mean(trace)- 2 * min(trace) ;
        offset = offset_array(i-1) + min(traces(:,i-1)) + offset_current;
    end
    offset_array(i) = offset;
    

    plot(t,  trace + offset,'Color',colors(i,:)); hold on;
    % elseif polarity(i) == 1
    %     trace = (traces(:,i) - mean(traces(:,i))) + offset;
    %     plot(t,  trace ,'Color',colors(i,:)); hold on;
    %     offset_current = (max(trace(:))-min(trace(:)));
    %     offset = offset + offset_current;
end
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
