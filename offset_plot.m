function [offset_array] = offset_plot(traces,t,options)

if nargin > 2
    if isfield(options, 'colors')
        colors= options.colors;
    else
        colors = lines; % 默认值
    end

end
offset_array = zeros(size(traces,2));
offset = 0;
for i = 1:size(traces,2)
    offset_array(i) = offset;
    % plot trace
    trace = (traces(:,i)-1)*-1 + offset;
    plot(t,  trace ,'Color',colors(i,:)); hold on;
    offset_current = (max(trace(:))-min(trace(:)))*1;
    offset = offset + offset_current;
end
ylim tight
xlim tight
sgtitle('Peak finding');
hold on;
end